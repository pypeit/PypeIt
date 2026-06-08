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
# Module-level helper
# ---------------------------------------------------------------------------

def _uniq(x):
    """
    Return the indices of the *last* occurrence of each unique value in a sorted
    array.

    Replicates the IDL ``UNIQ()`` behaviour used internally for building the
    design-matrix span boundaries.

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
    return np.flatnonzero(np.concatenate(([True], x[1:] != x[:-1], [True])))[1:] - 1


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
        Independent variable used to set breakpoints.  Either ``x`` or
        ``fullbkpt`` must be provided.
    fullbkpt : :class:`numpy.ndarray`, optional
        Full pre-built padded knot vector.  Sorted and cast to float on input.
        If its length is less than ``2 * nord`` it is padded by
        :meth:`_fill_bkpt`.
    nord : int, optional
        B-spline order (default 4 = cubic).
    bkpt : :class:`numpy.ndarray`, optional
        Interior breakpoints supplied directly.  Points outside the range of
        ``x`` are omitted.  If fewer than 2 remain, endpoints of ``x`` are used.
    bkspread : float, optional
        Scale factor for the phantom-knot spacing.
    bkspace : float, optional
        Fixed spacing between interior breakpoints.
    nbkpts : int, optional
        Number of breakpoints spanning the range of ``x``.
    everyn : int or float, optional
        Place a breakpoint at every ``everyn``-th value of ``x``.
    """

    def __init__(
        self, x=None, fullbkpt=None, nord=4, bkpt=None, bkspread=1.0, bkspace=None, nbkpts=None,
        everyn=None
    ):
        self.nord = nord
        self.npoly = 1  # Always 1 for the base class

        if x is None and fullbkpt is None:
            # Empty instance (e.g. for copy / deserialization)
            self.breakpoints = None
            self.mask = None
            self.coeff = None
            self.icoeff = None
            self._cached_design = None
            self._cached_x_shape = None
            return

        self.breakpoints = BSpline._build_breakpoints(
            x=x, fullbkpt=fullbkpt, nord=nord, bkpt=bkpt, bkspread=bkspread,
            bkspace=bkspace, nbkpts=nbkpts, everyn=everyn,
        )
        nc = self.breakpoints.size - nord
        self.mask = np.ones(self.breakpoints.size, dtype=bool)
        self.coeff = np.zeros(nc, dtype=float)
        self.icoeff = np.zeros(nc, dtype=float)
        self._cached_design = None
        self._cached_x_shape = None

    # ------------------------------------------------------------------
    # Static helpers — breakpoint construction
    # ------------------------------------------------------------------

    @staticmethod
    def _fill_bkpt(bkpt, nord, bkspread):
        """
        Pad a breakpoint vector with ``nord - 1`` phantom knots at each end.

        The phantom spacing is ``(bkpt[1] - bkpt[0]) * bkspread``.

        Parameters
        ----------
        bkpt : :class:`numpy.ndarray`
            Interior breakpoints (at least 2 elements).
        nord : int
            B-spline order.
        bkspread : float
            Scale factor for phantom-knot spacing.

        Returns
        -------
        :class:`numpy.ndarray`
            Full padded knot vector of length ``len(bkpt) + 2 * (nord - 1)``.
        """
        bkspace = (bkpt[1] - bkpt[0]) * bkspread
        indx = np.arange(1, nord)
        return np.concatenate([bkpt[0] - bkspace * indx[::-1], bkpt,
                                bkpt[-1] + bkspace * indx])

    @staticmethod
    def _build_breakpoints(
        x=None, fullbkpt=None, nord=4, bkpt=None, bkspread=1.0, bkspace=None, nbkpts=None,
        everyn=None
    ):
        """
        Construct the full padded knot vector.

        One of ``fullbkpt`` or ``x`` must be provided.  When ``fullbkpt`` is given it is
        sorted, cast to float, and padded if its length is less than ``2 * nord``.
        Otherwise the interior breakpoints are derived from ``x`` using one of the
        ``bkspace`` / ``nbkpts`` / ``everyn`` / ``bkpt`` strategies and then passed through
        :meth:`_fill_bkpt`.

        Parameters
        ----------
        x : :class:`numpy.ndarray`, optional
            Independent variable.
        fullbkpt : :class:`numpy.ndarray`, optional
            Pre-built full breakpoint vector.
        nord : int, optional
            B-spline order.
        bkpt : :class:`numpy.ndarray`, optional
            Interior breakpoints supplied directly.
        bkspread : float, optional
            Phantom-knot spacing scale factor.
        bkspace : float, optional
            Fixed spacing between interior breakpoints.
        nbkpts : int, optional
            Number of interior breakpoints.
        everyn : int or float, optional
            Spacing between breakpoints in terms of ``x`` array indices.

        Returns
        -------
        :class:`numpy.ndarray`
            Full padded knot vector.

        Raises
        ------
        ValueError
            If neither ``x`` nor ``fullbkpt`` is provided, or if no breakpoint
            strategy can be determined.
        """
        if fullbkpt is not None:
            _fullbkpt = np.sort(fullbkpt, kind='heapsort').astype(float)
            if _fullbkpt.size < 2 * nord:
                _fullbkpt = BSpline._fill_bkpt(_fullbkpt, nord, bkspread)
            return _fullbkpt

        if x is None:
            raise ValueError('Must provide x to determine breakpoints.')

        sx = np.amin(x)
        ex = np.amax(x)
        if bkpt is None:
            if bkspace is not None:
                if bkspace >= ex - sx:
                    _bkpt = np.array([sx, ex])
                else:
                    _nbkpts = int((ex - sx) / bkspace) + 1
                    _bkpt = np.linspace(sx, ex, _nbkpts)
            elif nbkpts is not None:
                _bkpt = np.linspace(sx, ex, max(nbkpts, 2))
            elif everyn is not None:
                if everyn < x.size:
                    _nbkpts = max(x.size / everyn, 2.)
                    indx = (x.size / _nbkpts) * np.arange(_nbkpts)
                    _bkpt = np.interp(indx, np.arange(x.size, dtype=float), x)
                else:
                    _bkpt = np.array([sx, ex])
            else:
                raise ValueError('Insufficient information to set breakpoints.')
        else:
            _bkpt = np.sort(bkpt, kind='heapsort')
            w = (_bkpt >= sx) & (_bkpt <= ex)
            _bkpt = np.array([sx, ex]) if np.sum(w) < 2 else _bkpt[w]

        if _bkpt.size < 2:
            _bkpt = np.array([sx, ex])
        if _bkpt[0] > sx:
            _bkpt[0] = sx
        if _bkpt[-1] < ex:
            _bkpt[-1] = ex

        return BSpline._fill_bkpt(_bkpt, nord, bkspread).astype(float)

    # ------------------------------------------------------------------
    # Static helper — span index lookup
    # ------------------------------------------------------------------

    @staticmethod
    def _find_spans(x, breakpoints, nord):
        """
        Find the B-spline interval index for each value in ``x``.

        For each ``x[i]``, returns the index ``ileft`` such that
        ``breakpoints[ileft] <= x[i] < breakpoints[ileft + 1]``, clamped to
        ``[nord - 1, n - 1]`` where ``n = breakpoints.size - nord``.

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
        indx = np.zeros(x.size, dtype=int)
        ileft = nord - 1
        for i in range(x.size):
            while x[i] > breakpoints[ileft + 1] and ileft < n - 1:
                ileft += 1
            indx[i] = ileft
        return indx

    # ------------------------------------------------------------------
    # Private algorithmic methods
    # ------------------------------------------------------------------

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
        bkpt = self.breakpoints[self.mask]
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltap = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltam = np.zeros((x.size, self.nord), dtype=x.dtype)
        j = 0
        vnikx[:, 0] = 1.0
        while j < self.nord - 1:
            ipj = ileft + j + 1
            deltap[:, j] = bkpt[ipj] - x
            imj = ileft - j
            deltam[:, j] = x - bkpt[imj]
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
        nbkpt = self.mask.sum()
        if nbkpt < 2 * self.nord:
            warnings.warn(f'Order ({self.nord}) too low for {nbkpt} breakpoints.')
            return None, None, None

        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros(n - self.nord + 1, dtype=int)
        upper = np.zeros(n - self.nord + 1, dtype=int) - 1

        indx = BSpline._find_spans(x, self.breakpoints[self.mask], self.nord)
        A = self._bspline_basis(x, indx)

        aa = _uniq(indx)
        upper[indx[aa] - self.nord + 1] = aa
        rindx = indx[::-1]
        bb = _uniq(rindx)
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
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        bw = A.shape[1]
        nfull = nn * self.npoly

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
            itop = k * self.npoly
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
        n = self.mask.sum() - self.nord
        coeffbk = self.mask[self.nord:].nonzero()[0]
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
        goodbkpt = np.where(self.mask)[0]
        nbkpt = len(goodbkpt)
        if nbkpt <= 2 * self.nord:
            warnings.warn('Fewer good break points than order of b-spline. Returning...')
            return -2

        hmm = bad_cols[_uniq(bad_cols // self.npoly)] // self.npoly

        n = nbkpt - self.nord
        if np.any(hmm >= n):
            warnings.warn('Not enough unique points in Cholesky decomposition. Returning...')
            return -2
        test = np.zeros(nbkpt, dtype=bool)
        for jj in range(-int(np.ceil(self.nord / 2)), int(self.nord / 2.)):
            foo = np.where((hmm + jj) > 0, hmm + jj,
                           np.zeros(hmm.shape, dtype=hmm.dtype))
            inside = np.where(
                (foo + self.nord) < n - 1,
                foo + self.nord,
                np.zeros(hmm.shape, dtype=hmm.dtype) + n - 1,
            )
            if len(inside) > 0:
                test[inside] = True
        if test.any():
            reality = goodbkpt[test]
            if self.mask[reality].any():
                self.mask[reality] = False
                self._cached_design = None  # Mask changed — invalidate design cache
                return -1
            return -2
        return -2

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def fit(self, xdata, ydata, invvar):
        """
        Fit a weighted least-squares B-spline to ``(xdata, ydata)``.

        The design matrix is cached: if ``xdata`` has the same shape as the
        previous call *and* the breakpoint mask has not changed, the cached
        matrix is reused.  This makes repeated calls (e.g. sigma-clipping loops)
        efficient.

        .. note::

            ``xdata`` must be **sorted in ascending order** before being passed
            to this method.

        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable (sorted ascending).
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of ``ydata``.  Zero entries are effectively masked.

        Returns
        -------
        err : int
            Set to 0 on success; -1 if breakpoints were masked (Cholesky was
            degenerate — caller should retry); and -2 on failure (too few active
            breakpoints).
        yfit : :class:`numpy.ndarray`
            Fitted B-spline evaluated at ``xdata``.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            return -2, np.zeros(ydata.shape, dtype=float)

        nfull = nn * self.npoly
        mininf = 1.0e-10 * invvar.sum() / nfull

        # Build / retrieve cached design matrix
        if self._cached_design is None or xdata.shape != self._cached_x_shape:
            self._cached_design = self._build_design_matrix(xdata)
            self._cached_x_shape = xdata.shape
        A, lower, upper = self._cached_design
        if A is None:
            return -2, np.zeros(ydata.shape, dtype=float)

        alpha, beta = self._assemble_normal_equations(A, ydata, invvar, lower, upper)
        sol, chol, bad_cols = self._solve_banded(alpha, beta, mininf)

        if bad_cols[0] != -1:
            yfit = self._evaluate_model(A, lower, upper)
            return self._mask_breakpoints(bad_cols), yfit

        goodbk_idx = goodbk.nonzero()[0]
        self._update_coefficients(sol, chol, goodbk_idx)
        yfit = self._evaluate_model(A, lower, upper)
        return 0, yfit

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

        n = self.mask.sum() - self.nord
        yfit = self._evaluate_model(A, lower, upper)

        mask = np.ones(x.shape, dtype=bool)
        goodbk = self.mask.nonzero()[0]
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord - 1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size > 0:
            for jj in range(hmm.size):
                mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
                     & (x <= self.breakpoints[goodbk[hmm[jj] + 1] - 1])] = False

        return yfit[np.argsort(xsort, kind='stable')], mask

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

        The design matrix cache is *not* copied (the copy starts with a cold cache).

        Returns
        -------
        BSpline
            A new instance with copies of all stored arrays.
        """
        new = BSpline.__new__(BSpline)
        new.nord = self.nord
        new.npoly = self.npoly
        new.breakpoints = np.copy(self.breakpoints)
        new.mask = np.copy(self.mask)
        new.coeff = np.copy(self.coeff)
        new.icoeff = np.copy(self.icoeff)
        new._cached_design = None
        new._cached_x_shape = None
        return new


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
        ``'chebyshev'``, ``'poly'`` (monomial: :math:`1, x_2, x_2^2, \ldots`),
        or ``'poly1'`` (monomial without constant term; ill-conditioned — use
        with caution).
    **kwargs
        Remaining keyword arguments forwarded to :class:`BSpline.__init__`.
    """

    def __init__(self, x=None, npoly=1, xmin=0.0, xmax=1.0, funcname='legendre',
                 fullbkpt=None, nord=4, bkpt=None, bkspread=1.0,
                 bkspace=None, nbkpts=None, everyn=None):
        super().__init__(x=x, fullbkpt=fullbkpt, nord=nord, bkpt=bkpt, bkspread=bkspread,
                         bkspace=bkspace, nbkpts=nbkpts, everyn=everyn)

        self.npoly = npoly
        self.xmin = xmin
        self.xmax = xmax
        self.funcname = funcname
        self._cached_x2_shape = None

        if funcname == 'poly1':
            warnings.warn(
                "BSpline2D funcname='poly1' produces a basis with no constant x2 term "
                "and may be ill-conditioned.  Consider using 'legendre' instead."
            )

        # Override coeff / icoeff to 2D shape (nc, npoly)
        if self.breakpoints is not None:
            nc = self.breakpoints.size - nord
            self.coeff = np.zeros((nc, npoly), dtype=float)
            self.icoeff = np.zeros((nc, npoly), dtype=float)

    # ------------------------------------------------------------------
    # Private methods — 2D-specific
    # ------------------------------------------------------------------

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
        if self.funcname == 'poly':
            P = np.ones((nx, self.npoly), dtype=float)
            for i in range(1, self.npoly):
                P[:, i] = P[:, i - 1] * x2norm
        elif self.funcname == 'poly1':
            P = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
            for i in range(1, self.npoly):
                P[:, i] = P[:, i - 1] * x2norm
        elif self.funcname == 'chebyshev':
            P = basis.fchebyshev(x2norm, self.npoly)
        elif self.funcname == 'legendre':
            P = basis.flegendre(x2norm, self.npoly)
        else:
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

        nbkpt = self.mask.sum()
        if nbkpt < 2 * self.nord:
            warnings.warn(f'Order ({self.nord}) too low for {nbkpt} breakpoints.')
            return None, None, None

        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros(n - self.nord + 1, dtype=int)
        upper = np.zeros(n - self.nord + 1, dtype=int) - 1

        indx = BSpline._find_spans(x, self.breakpoints[self.mask], self.nord)
        B = self._bspline_basis(x, indx)  # (N, nord)

        aa = _uniq(indx)
        upper[indx[aa] - self.nord + 1] = aa
        rindx = indx[::-1]
        bb = _uniq(rindx)
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
        n = self.mask.sum() - self.nord
        coeffbk = self.mask[self.nord:].nonzero()[0]
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
        self.coeff[goodbk_idx, :] = sol[:nfull].reshape(nn, self.npoly).astype(
            self.coeff.dtype)
        self.icoeff[goodbk_idx, :] = chol[0, :nfull].reshape(nn, self.npoly).astype(
            self.icoeff.dtype)

    # ------------------------------------------------------------------
    # Public API — overrides with required x2
    # ------------------------------------------------------------------

    def fit(self, xdata, ydata, invvar, x2):
        """Fit a weighted least-squares 2D B-spline.

        .. note::

            Both ``xdata`` and ``x2`` must be **sorted in ascending order by the
            same permutation** before being passed to this method.

        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable (sorted ascending).
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of ``ydata``.
        x2 : :class:`numpy.ndarray`
            Second variable (required).  Must be statistically independent of
            ``xdata``; see class-level warning.

        Returns
        -------
        err : int
            0 on success; -1 if breakpoints were masked; -2 on failure.
        yfit : :class:`numpy.ndarray`
            Fitted model at ``xdata``.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            return -2, np.zeros(ydata.shape, dtype=float)

        nfull = nn * self.npoly
        mininf = 1.0e-10 * invvar.sum() / nfull

        # Design matrix depends on both x and x2; invalidate cache if shapes change
        if (self._cached_design is None
                or xdata.shape != self._cached_x_shape
                or x2.shape != self._cached_x2_shape):
            self._cached_design = self._build_design_matrix(xdata, x2)
            self._cached_x_shape = xdata.shape
            self._cached_x2_shape = x2.shape
        A, lower, upper = self._cached_design
        if A is None:
            return -2, np.zeros(ydata.shape, dtype=float)

        alpha, beta = self._assemble_normal_equations(A, ydata, invvar, lower, upper)
        sol, chol, bad_cols = self._solve_banded(alpha, beta, mininf)

        if bad_cols[0] != -1:
            yfit = self._evaluate_model(A, lower, upper)
            return self._mask_breakpoints(bad_cols), yfit

        goodbk_idx = goodbk.nonzero()[0]
        self._update_coefficients(sol, chol, goodbk_idx)
        yfit = self._evaluate_model(A, lower, upper)
        return 0, yfit

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

        n = self.mask.sum() - self.nord
        yfit = self._evaluate_model(A, lower, upper)

        mask = np.ones(x.shape, dtype=bool)
        goodbk = self.mask.nonzero()[0]
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord - 1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size > 0:
            for jj in range(hmm.size):
                mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
                     & (x <= self.breakpoints[goodbk[hmm[jj] + 1] - 1])] = False

        return yfit[np.argsort(xsort, kind='stable')], mask

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
        new.breakpoints = np.copy(self.breakpoints)
        new.mask = np.copy(self.mask)
        new.coeff = np.copy(self.coeff)
        new.icoeff = np.copy(self.icoeff)
        new.xmin = self.xmin
        new.xmax = self.xmax
        new.funcname = self.funcname
        new._cached_design = None
        new._cached_x_shape = None
        new._cached_x2_shape = None
        return new
