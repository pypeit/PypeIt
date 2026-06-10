"""
Refactored B-spline classes with C-order memory layout and scipy-based linear
algebra.

This module provides :class:`BSpline` (1D) and :class:`BSpline2D` (quasi-2D with
a polynomial second-variable dependence) as self-contained classes with no
dependency on PypeIt infrastructure.  FITS serialisation will be provided by
thin wrapper classes in a separate follow-on module.

Key design choices vs the original :class:`~pypeit.bspline.bspline.bspline`:

- All arrays use native C (row-major) memory order.  No ``order='F'``
  allocations, ``flatten('F')``, or ``reshape(..., order='F')`` appear anywhere.
- ``coeff`` for 2D fits has shape ``(nc, npoly)`` (knot index first) rather than
  the original ``(npoly, nc)`` (Fortran-column order).
- Matrix multiplications use the ``@`` operator throughout.
- The banded Cholesky solve is delegated to
  :func:`scipy.linalg.cholesky_banded` / :func:`scipy.linalg.cho_solve_banded`
  (LAPACK ``dpbtrf`` / ``dpbtrs``).
- The design matrix is cached between :meth:`BSpline.fit` calls so that repeated
  calls with the same ``x`` (e.g. sigma-clipping loops) do not recompute it.
- The ``fit`` and ``workit`` methods of the original are merged into a single
  :meth:`BSpline.fit`.
- The quasi-2D behaviour (polynomial basis in a second variable) is isolated in
  :class:`BSpline2D`, which subclasses :class:`BSpline`.
- ``flatten()`` is not used anywhere; ``reshape()`` is used only where it
  produces a view (i.e. the array is contiguous).

References
----------
Original implementation ported from IDL/PYDL:
  https://doi.org/10.5281/zenodo.1095150
"""

import warnings

import numpy as np
from scipy.linalg import cholesky_banded, cho_solve_banded, LinAlgError

from pypeit.core import basis


# ---------------------------------------------------------------------------
# Knots — knot-vector construction strategy
# ---------------------------------------------------------------------------

class Knots:
    """
    Construct and store breakpoints for B-spline fitting.

    The breakpoint knot vector is built lazily when :meth:`build` is called and
    then cached in :attr:`breakpoints`.

    When building the knots for any given fit, the parameters of this class are
    with respect to to the independent coordinate (``x``) and the order
    (``nord``).  Note that ``x`` and ``nord`` to not need to be provided at
    instantiation; the breakpoints can be (re)constructed later using
    :meth:`build`

    Parameters
    ----------
    interior : :class:`numpy.ndarray`, optional
        Interior knots supplied directly.  Any points outside the range of ``x``
        are omitted.
    spread : float, optional
        Scale factor used for spacing the phantom knots placed at the beginning
        and end of the fitting range.
    spacing : float, optional
        Fixed spacing in the ``x`` coordinate value between interior knots.
    count : int, optional
        Number of interior knots spanning the range of ``x``.
    stride : int or float, optional
        Place knots separated by this *number* of values in ``x``.
    full : :class:`numpy.ndarray`, optional
        Pre-built full padded knot vector.  Sorted and cast to float on input.
        If its length is less than twice the order of the fit (``nord``), it is
        padded by :meth:`_pad`.  When provided, all other strategy parameters
        are ignored.
    x : :class:`numpy.ndarray`, optional
        Independent variable used to determine interior breakpoints.  Ignored if
        :attr:`full` is not None, otherwise this must be provided and cannot be
        None.
    nord : int, optional
        B-spline order.  If provided, the breakpoints are built immediately
        using :meth:`build`.  If that method fails (e.g., because both ``x`` and
        ``full`` are not provided), the object is instantiated as if ``nord``
        was not provided.

    Notes
    -----
    The arguments ``spacing``, ``count``, and ``stride`` are mutually exclusive
    and given priority in that order.  That is, if ``spacing`` is provided, the
    other two are ignored, etc.

    The difference between ``spacing`` and ``stride`` is that ``spacing``
    requests a specific change in the ``x`` value, whereas ``stride`` simply
    requests a spacing between number of ``x`` values, regardless of whether or
    not ``x`` is a regular grid.
    """

    def __init__(
        self, interior=None, spread=1.0, spacing=None, count=None, stride=4, full=None, x=None,
        nord=None,
    ):
        self.interior = interior
        self.spread = spread
        self.spacing = spacing
        self.count = count
        self.stride = stride
        self.full = full
        self._breakpoints = None

        if nord is not None:
            try:
                self.build(x, nord)
            except ValueError as e:
                warnings.warn(
                    'Unable to build breakpoints at instantiation.  Error from build function: '
                    f'{e}'
                )

    @property
    def breakpoints(self):
        """
        Full padded knot vector.

        Returns
        -------
        :class:`numpy.ndarray`
            The 1D array of breakpoints; if they are unavailable
            (i.e.,:meth:`build` has not yet been called successfully), None is
            returned.
        """
        return self._breakpoints

    def build(self, x, nord):
        """
        Construct and store the breakpoint vector.

        .. important::

            This *always* rebuilds the set of breakpoints.  If you want to
            access a pre-existing set of breakpoints, use :attr:`breakpoints`.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable used to determine interior breakpoints.
            Ignored if :attr:`full` is not None, otherwise this must be provided
            and cannot be None.
        nord : int
            B-spline order.

        Raises
        ------
        ValueError
            If :attr:`full` and ``x`` ar both ``None``, the breakpoints cannot
            be defined.
        """
        if self.full is not None:
            _full = np.sort(self.full, kind='heapsort').astype(float)
            if _full.size < 2 * nord:
                _full = Knots._pad(_full, nord, self.spread)
            self._breakpoints = _full
            return

        if x is None:
            raise ValueError(
                'A prebuilt full knot vector was not provided during instantiation.  The build '
                'function requires the x argument.'
            )

        sx = np.amin(x)
        ex = np.amax(x)
        if self.interior is None:
            if self.spacing is not None:
                if self.spacing >= ex - sx:
                    _interior = np.array([sx, ex])
                else:
                    _count = int((ex - sx) / self.spacing) + 1
                    _interior = np.linspace(sx, ex, _count)
            elif self.count is not None:
                _interior = np.linspace(sx, ex, max(self.count, 2))
            elif self.stride is not None:
                if self.stride < x.size:
                    _count = max(x.size / self.stride, 2.)
                    indx = (x.size / _count) * np.arange(_count)
                    _interior = np.interp(indx, np.arange(x.size, dtype=float), x)
                else:
                    _interior = np.array([sx, ex])
            else:
                raise ValueError('Insufficient information to set breakpoints.')
        else:
            _interior = np.sort(self.interior, kind='heapsort')
            w = (_interior >= sx) & (_interior <= ex)
            _interior = np.array([sx, ex]) if np.sum(w) < 2 else _interior[w]

        if _interior.size < 2:
            _interior = np.array([sx, ex])
        if _interior[0] > sx:
            _interior[0] = sx
        if _interior[-1] < ex:
            _interior[-1] = ex

        self._breakpoints = Knots._pad(_interior, nord, self.spread).astype(float)

    @staticmethod
    def _pad(knots, nord, spread):
        """
        Pad a knot vector with ``nord - 1`` phantom knots at each end.

        The phantom spacing is ``(knots[1] - knots[0]) * spread``.

        Parameters
        ----------
        knots : :class:`numpy.ndarray`
            Interior knots (at least 2 elements).
        nord : int
            B-spline order.
        spread : float
            Scale factor for phantom-knot spacing.

        Returns
        -------
        :class:`numpy.ndarray`
            Full padded knot vector of length ``len(knots) + 2 * (nord - 1)``.
        """
        spacing = (knots[1] - knots[0]) * spread
        indx = np.arange(1, nord)
        return np.concatenate((
            knots[0] - spacing * indx[::-1], knots, knots[-1] + spacing * indx
        ))

    def copy(self):
        """
        Return a copy of this instance.

        Strategy parameters are copied by value.  The built breakpoints
        array is deep-copied via :func:`numpy.copy`.

        Returns
        -------
        Knots
            A new instance with the same strategy parameters and a copy
            of the built breakpoints.
        """
        new = Knots.__new__(Knots)
        new.interior = self.interior
        new.spread = self.spread
        new.spacing = self.spacing
        new.count = self.count
        new.stride = self.stride
        new.full = self.full
        new._breakpoints = (
            np.copy(self._breakpoints) if self._breakpoints is not None else None
        )
        return new


# ---------------------------------------------------------------------------
# BSpline — 1D weighted least-squares B-spline
# ---------------------------------------------------------------------------

class BSpline:
    r"""
    1D weighted least-squares B-spline fitting.

    Fits the model

    .. math::

        \hat{y}_i = \sum_k c_k B_k(x_i)

    in the weighted least-squares sense

    .. math::

        \min_c \sum_i w_i (y_i - \hat{y}_i)^2

    where :math:`w_i = \varepsilon^{-2}_i`, :math:`\varepsilon_i` are the
    1-sigma uncertainties, and :math:`B_k` are B-spline basis functions of order
    ``nord``.

    The normal equations :math:`(A^\top W A) c = A^\top W y` are assembled in
    banded storage and solved via LAPACK's banded Cholesky routines
    (:func:`scipy.linalg.cholesky_banded` /
    :func:`scipy.linalg.cho_solve_banded`).

    All arrays use native C (row-major) memory order.

    Parameters
    ----------
    x : :class:`numpy.ndarray`, optional
        Independent variable used to set breakpoints.  Passed to
        :meth:`Knots.build` when ``knots`` does not already hold a built
        knot vector.

    knots : :class:`Knots`, :class:`numpy.ndarray`, optional
        Fit breakpoints, which are defined by this parameter in one of the
        following ways:

            - If ``knots=None``, a default :class:`Knots` instance is created.
              Note that if ``x`` is provided, the code will attempt to construct
              the knots during instantiation.  See :class:`Knots` for the
              default construction parameters.

            - If a :class:`numpy.ndarray` is provided, this is treated as the
              pre-built, fully padded knot vector and used to instantiate
              :attr:`knots` as a :class:`Knots` instance (``knots`` is passed to
              :class:`Knots` using its ``full`` argument).

            - If the object is a :class:`Knots` instance, it is used as is.

    nord : int, optional
        B-spline order (default 4 = cubic).
    """

    def __init__(self, x=None, knots=None, nord=4):
        self.nord = nord

        if isinstance(knots, np.ndarray):
            self.knots = Knots(full=knots)
        elif knots is None:
            self.knots = Knots() if x is None else Knots(x=x, nord=nord)
        else:
            self.knots = knots
            
        if not isinstance(self.knots, Knots):
            raise TypeError(
                'The knots provided to BSpline must be None, a numpy.ndarray, or a Knots '
                f'instance, not a {type(knots).__name__}.'
            )

        if self.knots.breakpoints is None:
            try:
                self.knots.build(x, self.nord)
            except ValueError:
                # Assume the user will provide the information necessary to
                # instantiate the breakpoints when calling the fit function.
                pass

        self._init_result_storage()

    def _init_result_storage(self):
        """
        Initialize the arrays used to store the results of the fit
        """
        self._cached_design = None
        self._cached_x_shape = None

        if self.breakpoints is None:
            # Empty instance (e.g. for copy / deserialization)
            self.bkpt_gpm = None
            self.coeff = None
            self.icoeff = None
            return

        self.bkpt_gpm = np.ones(self.breakpoints.size, dtype=bool)

        self.reinit_coeff()

    def reinit_coeff(self):
        """
        Reset coefficient arrays to zero.

        Does *not* reset the breakpoints or mask.  Does *not* invalidate the
        design matrix cache (breakpoints have not changed).
        """
        nc = self.breakpoints.size - self.nord
        self.coeff = np.zeros(nc, dtype=float)
        self.icoeff = np.zeros(nc, dtype=float)

    def copy(self):
        """
        Return a deep copy of this instance.

        The design matrix cache is *not* copied (the copy starts with a
        cold cache).

        Returns
        -------
        BSpline
            A new instance with copies of all stored arrays.
        """
        new = BSpline.__new__(BSpline)
        new.nord = self.nord
        new.knots = self.knots.copy()
        new.bkpt_gpm = np.copy(self.bkpt_gpm)
        new.coeff = np.copy(self.coeff)
        new.icoeff = np.copy(self.icoeff)
        new._cached_design = None
        new._cached_x_shape = None
        return new

    # ------------------------------------------------------------------
    # Breakpoints property
    # ------------------------------------------------------------------

    @property
    def breakpoints(self):
        """
        Full padded knot vector.

        Returns
        -------
        :class:`numpy.ndarray
            Delegates to :attr:`knots.breakpoints <Knots.breakpoints>`.
        """
        return self.knots.breakpoints

    # ------------------------------------------------------------------
    # Static helper — span index lookup
    # ------------------------------------------------------------------

    @staticmethod
    def _find_spans(x, breakpoints, nord):
        """
        Find the B-spline interval index for each value in ``x``.

        For each ``x[i]``, returns the index ``j`` such that ``breakpoints[j] <=
        x[i] < breakpoints[j + 1]``, clipped to ``[nord - 1, n - 1]`` where ``n
        = breakpoints.size - nord``.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Data values.  **Must be monotonically non-decreasing.**
        breakpoints : :class:`numpy.ndarray`
            Active (unmasked) knot vector.
        nord : int
            B-spline order.

        Returns
        -------
        :class:`numpy.ndarray` of int, shape (N,)
            Span index for each element of ``x``.
        """
        n = breakpoints.size - nord
        return np.clip(np.searchsorted(breakpoints, x, side='right') - 1, nord - 1, n - 1)
        # indx = np.zeros(x.size, dtype=int)
        # ileft = nord - 1
        # for i in range(x.size):
        #     while x[i] > breakpoints[ileft + 1] and ileft < n - 1:
        #         ileft += 1
        #     indx[i] = ileft
        # return indx

    # ------------------------------------------------------------------
    # Static helper — unique run-end indices
    # ------------------------------------------------------------------

    @staticmethod
    def _uniq(x):
        """
        Return the index of the last occurrence of each unique value
        in a sorted array.

        Replicates the IDL ``UNIQ()`` behaviour used internally for
        building the design-matrix span boundaries.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            A sorted array (ascending).  Must be non-empty.

        Returns
        -------
        :class:`numpy.ndarray`
            Indices of the last occurrence of each unique value.

        Raises
        ------
        ValueError
            If ``x`` is empty.
        """
        if len(x) == 0:
            raise ValueError('No unique elements in an empty array.')
        return np.flatnonzero(
            np.concatenate(([True], x[1:] != x[:-1], [True]))
        )[1:] - 1

    # ------------------------------------------------------------------
    # Private algorithmic methods
    # ------------------------------------------------------------------

    def _poly_scale(self, n):
        """
        Convert a knot-span count to a total unknown count.

        For the 1D case this is the identity.  :class:`BSpline2D`
        overrides this to return ``n * self.npoly``.

        Parameters
        ----------
        n : int
            Number of active knot spans (or a span index).

        Returns
        -------
        int
            ``n`` (base class).
        """
        return n

    def _dedup_bad_cols(self, bad_cols):
        """
        Map bad Cholesky column indices to unique knot-span indices.

        For the 1D case (``npoly = 1``) column indices and span indices
        are identical.  :class:`BSpline2D` overrides this to fold
        multiple columns that belong to the same span via integer
        division by ``self.npoly``.

        Parameters
        ----------
        bad_cols : :class:`numpy.ndarray`
            Bad column indices returned by :meth:`_solve_banded`.

        Returns
        -------
        :class:`numpy.ndarray`
            Unique span indices corresponding to the bad columns.
        """
        return bad_cols[self._uniq(bad_cols)]

    def _bspline_basis(self, x, ileft):
        r"""
        Evaluate B-spline basis functions via the de Boor triangular recursion.

        For each data point ``x[i]``, evaluates the ``nord`` non-zero B-spline
        basis functions active at that point.  The recurrence relation is:

        .. math::

            B_{i,1}(x) &= 1 \\
            B_{i,k}(x) &= \frac{x - t_i}{t_{i+k-1} - t_i} B_{i,k-1}(x)
                        + \frac{t_{i+k} - x}{t_{i+k} - t_{i+1}} B_{i+1,k-1}(x)

        Parameters
        ----------
        x : :class:`numpy.ndarray`, shape (N,)
            Data values (sorted).
        ileft : :class:`numpy.ndarray` of int, shape (N,)
            Span index for each ``x[i]``, as returned by :meth:`_find_spans`.

        Returns
        -------
        :class:`numpy.ndarray`, shape (N, nord)
            B-spline basis values in C order.  ``vnikx[i, k]`` is the ``k``-th
            active basis function evaluated at ``x[i]``.
        """
        t = self.breakpoints[self.bkpt_gpm]
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltap = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltam = np.zeros((x.size, self.nord), dtype=x.dtype)
        j = 0
        vnikx[:, 0] = 1.0
        while j < self.nord - 1:
            deltap[:, j] = t[ileft + j + 1] - x
            deltam[:, j] = x - t[ileft - j]
            vmprev = 0.0
            for l in range(j + 1):
                vm = vnikx[:, l] / (deltap[:, l] + deltam[:, j - l])
                vnikx[:, l] = vm * deltap[:, l] + vmprev
                vmprev = vm * deltam[:, j - l]
            j += 1
            vnikx[:, j] = vmprev
        return vnikx

    def _build_design_matrix(self, x):
        """
        Construct the 1D B-spline design matrix.

        Returns the banded design matrix ``A`` of shape ``(N, nord)`` together with
        span boundary arrays ``lower`` and ``upper`` such that data points in the
        closed interval ``[lower[k], upper[k]]`` fall within B-spline span ``k``.

        Parameters
        ----------
        x : :class:`numpy.ndarray`, shape (N,)
            Independent variable values.  **Must be sorted in ascending order.**

        Returns
        -------
        A : :class:`numpy.ndarray`, shape (N, nord)
            Design matrix.  ``A[i, :]`` are the ``nord`` non-zero B-spline basis values
            at ``x[i]``.
        lower : :class:`numpy.ndarray` of int, shape (n - nord + 1,)
            First data index (inclusive) falling in each B-spline span.
        upper : :class:`numpy.ndarray` of int, shape (n - nord + 1,)
            Last data index (inclusive) falling in each B-spline span.
            ``upper[k] < lower[k]`` signals an empty span.

        Warns
        -----
        UserWarning
            If the number of active breakpoints is less than ``2 * nord``.
        """
        nbkpt = self.bkpt_gpm.sum()
        if nbkpt < 2 * self.nord:
            warnings.warn(f'Order ({self.nord}) too low for {nbkpt} breakpoints.')
            return None, None, None

        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros(n - self.nord + 1, dtype=int)
        upper = np.zeros(n - self.nord + 1, dtype=int) - 1

        indx = BSpline._find_spans(x, self.breakpoints[self.bkpt_gpm], self.nord)
        A = self._bspline_basis(x, indx)

        aa = self._uniq(indx)
        upper[indx[aa] - self.nord + 1] = aa
        rindx = indx[::-1]
        bb = self._uniq(rindx)
        lower[rindx[bb] - self.nord + 1] = nx - bb - 1

        return A, lower, upper

    def _assemble_normal_equations(self, A, y, w, lower, upper):
        r"""
        Assemble the banded normal equations :math:`A^\top W A\,c = A^\top W y`.

        Constructs the symmetric positive-definite banded matrix ``alpha`` and
        right-hand-side vector ``beta`` for the weighted least-squares problem.
        ``alpha`` uses upper-triangular banded storage compatible with
        :func:`scipy.linalg.cholesky_banded`:

        - ``alpha[0, j]`` = diagonal element :math:`(A^\top W A)_{jj}`

        - ``alpha[k, j]`` = ``k``-th superdiagonal element :math:`(A^\top W A)_{j-k,j}`

        The matrix is assembled span-by-span.  For each span ``k``, the
        contribution is the ``(bw, bw)`` block :math:`A_k^\top W_k A_k` (and
        similarly for ``beta``), which is accumulated into the banded storage
        via pre-computed flat index arrays.

        Parameters
        ----------
        A : :class:`numpy.ndarray`, shape (N, bw)
            Design matrix with bandwidth ``bw = nord * npoly``.
        y : :class:`numpy.ndarray`, shape (N,)
            Dependent variable.
        w : :class:`numpy.ndarray`, shape (N,)
            Inverse-variance weights.
        lower : :class:`numpy.ndarray` of int
            First data index (inclusive) in each B-spline span.
        upper : :class:`numpy.ndarray` of int
            Last data index (inclusive) in each B-spline span.

        Returns
        -------
        alpha : :class:`numpy.ndarray`, shape (bw, nfull + bw)
            Banded normal-equations matrix (padded by ``bw`` zero columns for
            compatibility with the original banded Cholesky routine).
        beta : :class:`numpy.ndarray`, shape (nfull + bw,)
            Right-hand side vector (padded by ``bw`` zeros).
        """
        goodbk = self.bkpt_gpm[self.nord:]
        nn = goodbk.sum()
        bw = A.shape[1]
        nfull = self._poly_scale(nn)

        sqrt_w = np.sqrt(w)
        a2 = A * sqrt_w[:, np.newaxis]  # whitened design matrix

        alpha = np.zeros((bw, nfull + bw), dtype=float)
        beta = np.zeros(nfull + bw, dtype=float)

        # bi: flat indices of the upper triangle (incl. diagonal) of a (bw, bw) block
        # bo: corresponding flat indices in alpha.T (offset to the current span)
        # Together they map each upper-triangle entry of (A_k^T W_k A_k) into alpha.
        bi = np.concatenate([np.arange(i) + (bw - i) * (bw + 1) for i in range(bw, 0, -1)])
        bo = np.concatenate([np.arange(i) + (bw - i) * bw for i in range(bw, 0, -1)])

        for k in range(nn - self.nord + 1):
            if upper[k] < lower[k]:
                continue
            sl = slice(lower[k], upper[k] + 1)
            itop = self._poly_scale(k)
            ibottom = min(itop, nfull) + bw - 1
            work = a2[sl, :].T @ a2[sl, :]          # (bw, bw) Gram block
            alpha.T.flat[bo + itop * bw] += work.flat[bi]
            beta[itop:ibottom + 1] += (y[sl] * sqrt_w[sl]) @ a2[sl, :]

        return alpha, beta

    def _solve_banded(self, alpha, beta, mininf):
        r"""
        Solve the banded normal equations via LAPACK banded Cholesky.

        Uses :func:`scipy.linalg.cholesky_banded` (LAPACK ``dpbtrf``) and
        :func:`scipy.linalg.cho_solve_banded` (LAPACK ``dpbtrs``).

        Before calling LAPACK, a pre-check is made for diagonal entries of
        ``alpha`` that are :math:`\leq \mathrm{mininf}` or non-finite.  If any
        are found, the bad column indices are returned immediately without
        attempting a decomposition.  This replicates the pre-check in the
        original ``cholesky_band`` implementation.

        Parameters
        ----------
        alpha : :class:`numpy.ndarray`, shape (bw, nfull + bw)
            Banded normal-equations matrix from
            :meth:`_assemble_normal_equations`.
        beta : :class:`numpy.ndarray`, shape (nfull + bw,)
            Right-hand side vector.
        mininf : float
            Threshold below which a diagonal entry is treated as degenerate.

        Returns
        -------
        sol : :class:`numpy.ndarray` or None
            Solution vector of length ``nfull``, or ``None`` on failure.
        chol : :class:`numpy.ndarray` or None
            Cholesky factor (same banded format as ``alpha[:, :nfull]``), or
            ``None`` on failure.  ``chol[0, :]`` is the diagonal of the factor
            (stored as :attr:`icoeff` via :meth:`_update_coefficients`).
        bad_cols : :class:`numpy.ndarray`
            ``np.array([-1])`` signals success; otherwise the indices of
            degenerate diagonal entries.
        """
        bw = alpha.shape[0]
        nfull = alpha.shape[1] - bw

        # Pre-check for degenerate diagonal entries
        neg_mask = (alpha[0, :nfull] <= mininf) | ~np.isfinite(alpha[0, :nfull])
        if neg_mask.any():
            return None, None, neg_mask.nonzero()[0]

        try:
            chol = cholesky_banded(alpha[:, :nfull], lower=True)
            sol = cho_solve_banded((chol, True), beta[:nfull])
        except LinAlgError:
            # Fallback: identify the first zero or negative diagonal
            bad = (alpha[0, :nfull] <= 0).nonzero()[0]
            if bad.size == 0:
                bad = np.array([0])
            return None, None, bad

        return sol, chol, np.array([-1])

    def _update_coefficients(self, sol, chol, goodbk_idx):
        """
        Store solution coefficients and Cholesky diagonal in instance attributes.

        Parameters
        ----------
        sol : :class:`numpy.ndarray`, shape (nfull,) = (nn,) for 1D
            Solution vector from :meth:`_solve_banded`.
        chol : :class:`numpy.ndarray`, shape (bw, nfull)
            Cholesky factor; ``chol[0, :]`` is the diagonal (stored in
            :attr:`icoeff`).
        goodbk_idx : :class:`numpy.ndarray` of int
            Indices within :attr:`coeff` of the currently active breakpoints.
        """
        nn = len(goodbk_idx)
        self.coeff[goodbk_idx] = sol[:nn].astype(self.coeff.dtype)
        self.icoeff[goodbk_idx] = chol[0, :nn].astype(self.icoeff.dtype)

    def _evaluate_model(self, A, lower, upper):
        r"""
        Evaluate the fitted B-spline model at the points encoded in ``A``.

        For each active span ``k``, computes:

        .. math::

            \hat{y}[\mathtt{lower}[k]:\mathtt{upper}[k]+1]
            = A[\mathtt{lower}[k]:\mathtt{upper}[k]+1, :]\; c_k

        where :math:`c_k = \mathtt{coeff}[k:k+\mathtt{nord}]` is the local
        coefficient block for span ``k``.

        Parameters
        ----------
        A : :class:`numpy.ndarray`, shape (N, nord)
            Design matrix from :meth:`_build_design_matrix`.
        lower : :class:`numpy.ndarray` of int
            First data index (inclusive) in each span.
        upper : :class:`numpy.ndarray` of int
            Last data index (inclusive) in each span.

        Returns
        -------
        :class:`numpy.ndarray`, shape (N,)
            Fitted model values.
        """
        n = self.bkpt_gpm.sum() - self.nord
        coeffbk = self.bkpt_gpm[self.nord:].nonzero()[0]
        goodcoeff = self.coeff[coeffbk]  # (nn,)

        yfit = np.zeros(A.shape[0], dtype=float)
        for k in range(n - self.nord + 1):
            if lower[k] <= upper[k]:
                sl = slice(lower[k], upper[k] + 1)
                yfit[sl] = A[sl, :] @ goodcoeff[k:k + self.nord]
        return yfit

    def _mask_breakpoints(self, bad_cols):
        """
        Mask breakpoints corresponding to degenerate Cholesky columns.

        Maps the bad column indices returned by :meth:`_solve_banded` back to
        breakpoint indices and masks a neighbourhood of size ``nord`` centred on
        each problem location.  After masking, the design matrix cache is
        invalidated so that the next :meth:`fit` call rebuilds it.

        Parameters
        ----------
        bad_cols : :class:`numpy.ndarray`
            Column indices flagged by :meth:`_solve_banded`.

        Returns
        -------
        int
            Value is -1 if masking succeeded (caller should retry :meth:`fit`),
            or -2 if too few breakpoints remain (caller should abort).
        """
        if not isinstance(bad_cols, np.ndarray):
            bad_cols = np.array([bad_cols])
        goodbkpt = np.where(self.bkpt_gpm)[0]
        nbkpt = len(goodbkpt)
        if nbkpt <= 2 * self.nord:
            warnings.warn('Fewer good break points than order of b-spline. Returning...')
            return -2

        hmm = self._dedup_bad_cols(bad_cols)

        n = nbkpt - self.nord
        if np.any(hmm >= n):
            warnings.warn('Not enough unique points in Cholesky decomposition. Returning...')
            return -2
        test = np.zeros(nbkpt, dtype=bool)
        for jj in range(-int(np.ceil(self.nord / 2)), int(self.nord / 2.)):
            foo = np.where((hmm + jj) > 0, hmm + jj, np.zeros(hmm.shape, dtype=hmm.dtype))
            inside = np.where(
                foo + self.nord < n - 1, foo + self.nord,
                np.zeros(hmm.shape, dtype=hmm.dtype) + n - 1,
            )
            if len(inside) > 0:
                test[inside] = True
        if test.any():
            reality = goodbkpt[test]
            if self.bkpt_gpm[reality].any():
                self.bkpt_gpm[reality] = False
                self._cached_design = None  # Mask changed — invalidate design cache
                return -1
            return -2
        return -2

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(self, x, y, ivar=None, reset_knots=False):
        """
        Fit a weighted least-squares B-spline to ``(x, y)``.

        The design matrix is cached: if ``x`` has the same shape as the
        previous call *and* the breakpoint mask has not changed, the cached
        matrix is reused.  This makes repeated calls (e.g. sigma-clipping loops)
        efficient.

        .. note::

            ``x`` must be **sorted in ascending order** before being passed
            to this method.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable (sorted ascending).
        y : :class:`numpy.ndarray`
            Dependent variable.
        ivar : :class:`numpy.ndarray`, optional
            Inverse variance of ``y``.  Zero entries are effectively masked.
            If not provided, uniform unit weights are used.
        reset_knots : bool, optional
            Regardless of any existing breakpoints, reset the breakpoints using
            :attr:`knots` and ``x``.

        Returns
        -------
        err : int
            Set to 0 on success; -1 if breakpoints were masked (Cholesky was
            degenerate — caller should retry); and -2 on failure (too few active
            breakpoints).
        yfit : :class:`numpy.ndarray`
            Fitted B-spline evaluated at ``x``.
        """
        if ivar is None:
            ivar = np.ones(y.size, dtype=float)

        if reset_knots:
            self.knots.build(x, self.nord)
            self._init_result_storage()

        goodbk = self.bkpt_gpm[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            return -2, np.zeros(y.shape, dtype=float)

        nfull = self._poly_scale(nn)
        mininf = 1.0e-10 * ivar.sum() / nfull

        # Build / retrieve cached design matrix.  This if statement should
        # always be true if reset_knots is true because of the call to
        # _init_result_storage.
        if self._cached_design is None or x.shape != self._cached_x_shape:
            self._cached_design = self._build_design_matrix(x)
            self._cached_x_shape = x.shape
        A, lower, upper = self._cached_design
        if A is None:
            return -2, np.zeros(y.shape, dtype=float)

        alpha, beta = self._assemble_normal_equations(A, y, ivar, lower, upper)
        sol, chol, bad_cols = self._solve_banded(alpha, beta, mininf)

        if bad_cols[0] != -1:
            return self._mask_breakpoints(bad_cols), self._evaluate_model(A, lower, upper)

        self._update_coefficients(sol, chol, goodbk.nonzero()[0])
        return 0, self._evaluate_model(A, lower, upper)

    def value(self, x):
        """
        Evaluate the fitted B-spline at arbitrary ``x`` positions.

        Sorts ``x`` internally so that :meth:`_find_spans` receives a monotone
        input, then un-sorts the result before returning.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable at which to evaluate the fit.

        Returns
        -------
        yfit : :class:`numpy.ndarray`
            Fitted model values at ``x``.
        mask : :class:`numpy.ndarray` of bool
            ``True`` where the evaluation is within the valid fitting range
            (between the outermost good breakpoints, excluding any gaps created
            by masked breakpoints).
        """
        xsort = x.argsort(kind='stable')
        A, lower, upper = self._build_design_matrix(x[xsort])

        n = self.bkpt_gpm.sum() - self.nord
        yfit = self._evaluate_model(A, lower, upper)

        mask = np.ones(x.shape, dtype=bool)
        goodbk = self.bkpt_gpm.nonzero()[0]
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord - 1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size > 0:
            for jj in range(hmm.size):
                mask[
                    (x >= self.breakpoints[goodbk[hmm[jj]]])
                    & (x <= self.breakpoints[goodbk[hmm[jj] + 1] - 1])
                ] = False

        return yfit[np.argsort(xsort, kind='stable')], mask


# ---------------------------------------------------------------------------
# BSpline2D — B-spline with polynomial second-variable dependence
# ---------------------------------------------------------------------------

class BSpline2D(BSpline):
    r"""
    B-spline with polynomial dependence on a second variable.

    Extends :class:`BSpline` to handle quasi-2D fitting where the model is:

    .. math::

        \hat{y}_i = \sum_{j,k} c_{jk}\, B_j(x_i)\, P_k(\xi_i)

    where :math:`B_j` are B-spline basis functions, :math:`P_k` is the ``k``-th
    polynomial basis function (Legendre, Chebyshev, or power), and :math:`\xi_i
    = 2(x_{2,i} - x_{\min}) / (x_{\max} - x_{\min}) - 1` is the normalised
    second variable.

    ``x2`` is a **required positional argument** in :meth:`fit` and :meth:`value`.  A
    :class:`BSpline2D` object always operates on two-variable data.

    The 2D coefficient array :attr:`coeff` has shape ``(nc, npoly)`` (knot index
    first), stored in C order.

    .. warning::

        ``x2`` must be statistically independent of ``x``.  If ``x2`` is a smooth,
        monotonic function of ``x``, the polynomial-B-spline product columns become
        nearly linearly dependent within each span, making the normal equations
        ill-conditioned.

    Parameters
    ----------
    x : :class:`numpy.ndarray`, optional
        Independent variable used to set breakpoints.
    npoly : int
        Polynomial order in the second variable.
    xmin : float, optional
        Minimum value of ``x2`` for normalisation (default 0.0).
    xmax : float, optional
        Maximum value of ``x2`` for normalisation (default 1.0).
    funcname : str, optional
        Polynomial basis type.  One of ``'legendre'`` (default),
        ``'chebyshev'``, ``'poly'`` (monomial:
        :math:`1, x_2, x_2^2, \ldots`), or ``'poly1'`` (monomial without
        constant term; ill-conditioned — use with caution).
    knots : :class:`Knots` or :class:`numpy.ndarray` or None
        Breakpoint specification forwarded to :class:`BSpline.__init__`.
    nord : int, optional
        B-spline order (default 4 = cubic).
    """

    def __init__(
        self, x=None, npoly=1, xmin=0.0, xmax=1.0, funcname='legendre', knots=None, nord=4
    ):
        if funcname == 'poly1':
            warnings.warn(
                "BSpline2D funcname='poly1' produces a basis with no constant x2 term "
                "and may be ill-conditioned.  Consider using 'legendre' instead."
            )

        # WARNING: self.npoly must be defined before running super().__init__()
        # so that it can be used by the call to reinit_coeffs() during
        # instantion of the base class.

        self.npoly = npoly
        self.xmin = xmin
        self.xmax = xmax
        self.funcname = funcname

        super().__init__(x=x, knots=knots, nord=nord)

    def _init_result_storage(self):
        """
        Initialize the arrays used to store the results of the fit
        """
        super()._init_result_storage()
        self._cached_x2_shape = None

    def reinit_coeff(self):
        """
        Reset 2D coefficient arrays to zero.
        """
        nc = self.breakpoints.size - self.nord
        self.coeff = np.zeros((nc, self.npoly), dtype=float)
        self.icoeff = np.zeros((nc, self.npoly), dtype=float)

    def copy(self):
        """
        Return a deep copy of this :class:`BSpline2D` instance.

        Returns
        -------
        BSpline2D
            A new instance with copies of all stored arrays.
        """
        new = BSpline2D.__new__(BSpline2D)
        new.nord = self.nord
        new.npoly = self.npoly
        new.knots = self.knots.copy()
        new.bkpt_gpm = np.copy(self.bkpt_gpm)
        new.coeff = np.copy(self.coeff)
        new.icoeff = np.copy(self.icoeff)
        new.xmin = self.xmin
        new.xmax = self.xmax
        new.funcname = self.funcname
        new._cached_design = None
        new._cached_x_shape = None
        new._cached_x2_shape = None
        return new

    # ------------------------------------------------------------------
    # Private methods — 2D-specific
    # ------------------------------------------------------------------

    def _poly_scale(self, n):
        """
        Convert a knot-span count to a total unknown count.

        Overrides the base-class identity with ``n * self.npoly``.

        Parameters
        ----------
        n : int
            Number of active knot spans (or a span index).

        Returns
        -------
        int
            ``n * self.npoly``.
        """
        return n * self.npoly

    def _dedup_bad_cols(self, bad_cols):
        """
        Map bad Cholesky column indices to unique knot-span indices.

        Overrides the base-class implementation to account for the
        ``npoly`` polynomial terms packed into each knot span.

        Parameters
        ----------
        bad_cols : :class:`numpy.ndarray`
            Bad column indices returned by :meth:`_solve_banded`.

        Returns
        -------
        :class:`numpy.ndarray`
            Unique span indices corresponding to the bad columns.
        """
        return bad_cols[self._uniq(bad_cols // self.npoly)] // self.npoly

    def _normalize_x2(self, x2):
        """
        Map ``x2`` linearly to the interval ``[-1, +1]``.

        Uses the stored :attr:`xmin` and :attr:`xmax` as the normalisation
        range.

        Parameters
        ----------
        x2 : :class:`numpy.ndarray`
            Second variable values.

        Returns
        -------
        :class:`numpy.ndarray`
            Normalised second variable in ``[-1, +1]``.
        """
        return 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0

    def _poly_basis(self, x2norm):
        """
        Evaluate the polynomial basis in the normalised second variable.

        Dispatches on :attr:`funcname` to build the ``(N, npoly)`` polynomial
        basis matrix.

        Parameters
        ----------
        x2norm : :class:`numpy.ndarray`, shape (N,)
            Second variable values normalised to ``[-1, +1]``.

        Returns
        -------
        :class:`numpy.ndarray`, shape (N, npoly)
            Polynomial basis values.  ``P[i, k]`` is the ``k``-th basis function
            evaluated at ``x2norm[i]``.

        Raises
        ------
        ValueError
            If :attr:`funcname` is not one of the recognised types.
        """
        nx = x2norm.size
        match self.funcname:
            case 'poly':
                P = np.ones((nx, self.npoly), dtype=float)
                for i in range(1, self.npoly):
                    P[:, i] = P[:, i - 1] * x2norm
            case 'poly1':
                P = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
                for i in range(1, self.npoly):
                    P[:, i] = P[:, i - 1] * x2norm
            case 'chebyshev':
                P = basis.fchebyshev(x2norm, self.npoly)
            case 'legendre':
                P = basis.flegendre(x2norm, self.npoly)
            case _:
                raise ValueError(
                    f"Unknown funcname '{self.funcname}'.  "
                    "Use 'legendre', 'chebyshev', 'poly', or 'poly1'."
                )
        return P

    def _build_design_matrix(self, x, x2):
        r"""
        Construct the 2D B-spline design matrix via a vectorised outer product.

        The design matrix has shape ``(N, bw)`` where ``bw = nord * npoly`` and:

        .. math::

            A_{i,\, ii \cdot \mathtt{npoly} + jj}
            = B_{i, ii} \cdot P_{i, jj}

        This is computed without any nested loops by broadcasting:

        .. code-block:: python

            A_3d = B[:, :, np.newaxis] * P[:, np.newaxis, :]  # (N, nord, npoly)
            A = A_3d.reshape(N, nord * npoly)                  # view, no copy

        Parameters
        ----------
        x : :class:`numpy.ndarray`, shape (N,)
            Independent variable (sorted ascending).
        x2 : :class:`numpy.ndarray`, shape (N,)
            Second variable.

        Returns
        -------
        A : :class:`numpy.ndarray`, shape (N, bw)
            Design matrix.
        lower : :class:`numpy.ndarray` of int
            First data index (inclusive) in each B-spline span.
        upper : :class:`numpy.ndarray` of int
            Last data index (inclusive) in each B-spline span.

        Raises
        ------
        ValueError
            If ``x2.size != x.size``.
        """
        if x2.size != x.size:
            raise ValueError('Dimensions of x and x2 do not match.')

        nbkpt = self.bkpt_gpm.sum()
        if nbkpt < 2 * self.nord:
            warnings.warn(f'Order ({self.nord}) too low for {nbkpt} breakpoints.')
            return None, None, None

        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros(n - self.nord + 1, dtype=int)
        upper = np.zeros(n - self.nord + 1, dtype=int) - 1

        indx = BSpline._find_spans(x, self.breakpoints[self.bkpt_gpm], self.nord)
        B = self._bspline_basis(x, indx)  # (N, nord)

        aa = self._uniq(indx)
        upper[indx[aa] - self.nord + 1] = aa
        rindx = indx[::-1]
        bb = self._uniq(rindx)
        lower[rindx[bb] - self.nord + 1] = nx - bb - 1

        x2norm = self._normalize_x2(x2)
        P = self._poly_basis(x2norm)  # (N, npoly)

        # Vectorised outer product; reshape returns a view (no copy)
        A_3d = B[:, :, np.newaxis] * P[:, np.newaxis, :]   # (N, nord, npoly)
        A = A_3d.reshape(nx, self.nord * self.npoly)         # (N, bw)

        return A, lower, upper

    def _evaluate_model(self, A, lower, upper):
        r"""
        Evaluate the fitted 2D B-spline model.

        For each active span ``k``, computes:

        .. math::

            \hat{y}[\mathtt{sl}]
            = \sum_{ii, jj} A[\mathtt{sl}, ii \cdot \mathtt{npoly}+jj]
              \cdot \mathtt{coeff}[k+ii, jj]

        implemented via ``einsum`` (no ``flatten``):

        .. code-block:: python

            yfit[sl] = np.einsum('nij,ij->n',
                                 A[sl, :].reshape(-1, nord, npoly),
                                 coeff[k:k+nord, :])

        The ``reshape`` call produces a *view* (no memory copy) because ``A`` is a
        contiguous C-order array.

        Parameters
        ----------
        A : :class:`numpy.ndarray`, shape (N, bw)
            Design matrix with ``bw = nord * npoly``.
        lower : :class:`numpy.ndarray` of int
            First data index (inclusive) in each span.
        upper : :class:`numpy.ndarray` of int
            Last data index (inclusive) in each span.

        Returns
        -------
        :class:`numpy.ndarray`, shape (N,)
            Fitted model values.
        """
        n = self.bkpt_gpm.sum() - self.nord
        coeffbk = self.bkpt_gpm[self.nord:].nonzero()[0]
        goodcoeff = self.coeff[coeffbk, :]  # (nn, npoly)

        yfit = np.zeros(A.shape[0], dtype=float)
        for k in range(n - self.nord + 1):
            if lower[k] <= upper[k]:
                sl = slice(lower[k], upper[k] + 1)
                # reshape is a view (A is C-contiguous); coeff slice is (nord, npoly)
                yfit[sl] = np.einsum(
                    'nij,ij->n',
                    A[sl, :].reshape(-1, self.nord, self.npoly),
                    goodcoeff[k:k + self.nord, :],
                )
        return yfit

    def _update_coefficients(self, sol, chol, goodbk_idx):
        """
        Store 2D solution coefficients and Cholesky diagonal.

        The solution vector ``sol`` is ordered with polynomial index fastest:
        ``sol[k * npoly + jj] = coeff[k, jj]``.  The C-order reshape
        ``sol.reshape(nn, npoly)`` recovers this without any Fortran-order
        tricks.

        Parameters
        ----------
        sol : :class:`numpy.ndarray`, shape (nfull,) = (nn * npoly,)
            Solution vector from :meth:`_solve_banded`.
        chol : :class:`numpy.ndarray`, shape (bw, nfull)
            Cholesky factor; ``chol[0, :]`` is the diagonal stored as :attr:`icoeff`.
        goodbk_idx : :class:`numpy.ndarray` of int
            Active breakpoint indices within the coefficient array.
        """
        nn = len(goodbk_idx)
        nfull = nn * self.npoly
        # C-order reshape: sol[k*npoly + jj] -> coeff_block[k, jj]
        self.coeff[goodbk_idx, :] = sol[:nfull].reshape(nn, self.npoly).astype(self.coeff.dtype)
        self.icoeff[goodbk_idx, :] = chol[0, :nfull].reshape(nn, self.npoly).astype(
            self.icoeff.dtype
        )

    # ------------------------------------------------------------------
    # Public API — overrides with required x2
    # ------------------------------------------------------------------

    def fit(self, x, y, x2, ivar=None, reset_knots=False):
        """
        Fit a weighted least-squares 2D B-spline.

        .. note::

            Both ``x`` and ``x2`` must be **sorted in ascending order by the
            same permutation** before being passed to this method.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable (sorted ascending).
        y : :class:`numpy.ndarray`
            Dependent variable.
        x2 : :class:`numpy.ndarray`
            Second variable (required).  Must be statistically independent of
            ``x``; see class-level warning.
        ivar : :class:`numpy.ndarray`, optional
            Inverse variance of ``y``.  Zero entries are effectively masked.
            If not provided, uniform unit weights are used.
        reset_knots : bool, optional
            Regardless of any existing breakpoints, reset the breakpoints using
            :attr:`knots` and ``x``.

        Returns
        -------
        err : int
            0 on success; -1 if breakpoints were masked; -2 on failure.
        yfit : :class:`numpy.ndarray`
            Fitted model at ``x``.
        """
        if ivar is None:
            ivar = np.ones(y.size, dtype=float)

        if reset_knots:
            self.knots.build(x, self.nord)
            self._init_result_storage()

        goodbk = self.bkpt_gpm[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            return -2, np.zeros(y.shape, dtype=float)

        nfull = self._poly_scale(nn)
        mininf = 1.0e-10 * ivar.sum() / nfull

        # Design matrix depends on both x and x2; invalidate cache if shapes change
        if (
            self._cached_design is None
            or x.shape != self._cached_x_shape
            or x2.shape != self._cached_x2_shape
        ):
            self._cached_design = self._build_design_matrix(x, x2)
            self._cached_x_shape = x.shape
            self._cached_x2_shape = x2.shape
        A, lower, upper = self._cached_design
        if A is None:
            return -2, np.zeros(y.shape, dtype=float)

        alpha, beta = self._assemble_normal_equations(A, y, ivar, lower, upper)
        sol, chol, bad_cols = self._solve_banded(alpha, beta, mininf)

        if bad_cols[0] != -1:
            return self._mask_breakpoints(bad_cols), self._evaluate_model(A, lower, upper)

        goodbk_idx = goodbk.nonzero()[0]
        self._update_coefficients(sol, chol, goodbk_idx)
        return 0, self._evaluate_model(A, lower, upper)

    def value(self, x, x2):
        """
        Evaluate the fitted 2D B-spline at arbitrary ``(x, x2)`` positions.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`
            Second variable.

        Returns
        -------
        yfit : :class:`numpy.ndarray`
            Fitted model values.
        mask : :class:`numpy.ndarray` of bool
            ``True`` where the evaluation is within the valid fitting range.
        """
        xsort = x.argsort(kind='stable')
        A, lower, upper = self._build_design_matrix(x[xsort], x2[xsort])

        n = self.bkpt_gpm.sum() - self.nord
        yfit = self._evaluate_model(A, lower, upper)

        mask = np.ones(x.shape, dtype=bool)
        goodbk = self.bkpt_gpm.nonzero()[0]
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord - 1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size > 0:
            for jj in range(hmm.size):
                mask[
                    (x >= self.breakpoints[goodbk[hmm[jj]]])
                    & (x <= self.breakpoints[goodbk[hmm[jj] + 1] - 1])
                ] = False

        return yfit[np.argsort(xsort, kind='stable')], mask
