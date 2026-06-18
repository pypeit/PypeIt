"""
FITS-serialisable wrappers for the :class:`~pypeit.bspline.refactor.BSpline` and
:class:`~pypeit.bspline.refactor.BSpline2D` classes.

Each container inherits from both :class:`~pypeit.datamodel.DataContainer` (for FITS
I/O) and its respective algorithmic base class (for all fitting and evaluation
methods).  The dual inheritance requires a specific construction order:
:class:`~pypeit.datamodel.DataContainer.__init__` must run first to establish the
attribute-enforcement mechanism, after which the base-class constructor runs against
the already-prepared ``__dict__``.

.. note::

    The :attr:`~pypeit.bspline.refactor.BSpline.breakpoints` property defined on
    :class:`~pypeit.bspline.refactor.BSpline` is a read-only data descriptor that
    delegates to ``self.knots.breakpoints``.  Using ``breakpoints`` as a datamodel
    key would create an irresolvable naming conflict (the property intercepts all
    attribute-level reads, making the stored value invisible).  The breakpoint
    array is therefore stored under the distinct key ``bkpt_full``, and the
    property continues to serve all algorithmic code via ``self.knots``.

.. include:: ../include/links.rst
"""

import numpy as np

from pypeit import datamodel
from pypeit.core.bspline.refactor import BSpline, BSpline2D, Knots


class BSplineContainer(datamodel.DataContainer, BSpline):
    """
    FITS-serialisable wrapper for :class:`~pypeit.bspline.refactor.BSpline`.

    Inherits all fitting and evaluation methods from
    :class:`~pypeit.bspline.refactor.BSpline` and adds FITS I/O via
    :class:`~pypeit.datamodel.DataContainer`.  Instances can be used wherever a
    plain :class:`~pypeit.bspline.refactor.BSpline` is expected.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_bsplinecontainer.rst

    Parameters
    ----------
    x : :class:`numpy.ndarray`, optional
        Independent variable used to build the initial breakpoint grid.  Passed
        directly to :class:`~pypeit.bspline.refactor.BSpline`.
    knots : :class:`~pypeit.bspline.refactor.Knots` or :class:`numpy.ndarray` or None, optional
        Knot specification.  Passed directly to
        :class:`~pypeit.bspline.refactor.BSpline`.
    nord : int, optional
        Order of the b-spline (default 4).

    See Also
    --------
    pypeit.bspline.refactor.BSpline
    """

    version = '2.0.0'

    datamodel = {
        'nord': dict(
            otype=int,
            descr='Order of the b-spline.'
        ),
        'bkpt_full': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Full padded knot vector stored for serialisation.'
        ),
        'bkpt_gpm': dict(
            otype=np.ndarray,
            atype=np.bool_,
            descr='Active-breakpoint mask; True where a breakpoint participates in the fit.'
        ),
        'coeff': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Fitted b-spline coefficients, shape (nc,).'
        ),
        'icoeff': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Inverse-covariance diagonal, shape (nc,).'
        ),
    }

    internals = ['knots', 'x', 'yfit', '_cached_design']

    def __init__(self, x=None, knots=None, nord=4):
        datamodel.DataContainer.__init__(self)
        BSpline.__init__(self, x=x, knots=knots, nord=nord)
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints

    @classmethod
    def from_bspline(cls, spl):
        """
        Construct a :class:`BSplineContainer` from an existing
        :class:`~pypeit.bspline.refactor.BSpline`.

        Parameters
        ----------
        spl : :class:`~pypeit.bspline.refactor.BSpline`
            Source B-spline.  Its
            :attr:`~pypeit.bspline.refactor.BSpline.breakpoints` must not be
            ``None``.

        Returns
        -------
        :class:`BSplineContainer`
            A new container wrapping a copy of the spline state.  The transient
            internals ``x`` and ``yfit`` are not copied; call
            :meth:`~pypeit.bspline.refactor.BSpline.fit` or use
            :meth:`copy` if those are needed.
        """
        return cls.from_dict(d={
            'nord': spl.nord,
            'bkpt_full': np.copy(spl.knots.breakpoints),
            'bkpt_gpm': np.copy(spl.bkpt_gpm),
            'coeff': None if spl.coeff is None else np.copy(spl.coeff),
            'icoeff': None if spl.icoeff is None else np.copy(spl.icoeff),
        })

    def copy(self):
        """
        Return a deep copy as a new :class:`BSplineContainer`.

        The design-matrix cache is not copied (the copy starts with a cold
        cache), matching the behaviour of
        :meth:`pypeit.bspline.refactor.BSpline.copy`.

        Returns
        -------
        :class:`BSplineContainer`
        """
        new = type(self).from_bspline(self)
        new.x = self.x
        new.yfit = None if self.yfit is None else np.copy(self.yfit)
        return new

    def fit(self, x, y, ivar=None, reset_knots=False):
        """
        Fit a weighted least-squares B-spline.

        This delegates entirely to :meth:`~pypeit.bspline.refactor.BSpline.fit`,
        but updates :attr:`bkpt_full`, as necessary.  See documentation for
        :meth:`~pypeit.bspline.refactor.BSpline.fit`.
        """
        result = super().fit(x, y, ivar=ivar, reset_knots=reset_knots)
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints
        return result

    def _validate(self):
        if self.bkpt_full is not None and self.knots is None:
            self.knots = Knots(full=self.bkpt_full, nord=self.nord)
            if self.bkpt_gpm is None:
                self.bkpt_gpm = np.ones(self.knots.breakpoints.size, dtype=bool)

    def _bundle(self):
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints
        return super()._bundle(ext='BSPLINE')


class BSpline2DContainer(datamodel.DataContainer, BSpline2D):
    """
    FITS-serialisable wrapper for :class:`~pypeit.bspline.refactor.BSpline2D`.

    Inherits all fitting and evaluation methods from
    :class:`~pypeit.bspline.refactor.BSpline2D` and adds FITS I/O via
    :class:`~pypeit.datamodel.DataContainer`.  Instances can be used wherever a
    plain :class:`~pypeit.bspline.refactor.BSpline2D` is expected.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_bspline2dcontainer.rst

    The ``basis`` datamodel field is populated only when
    :meth:`~pypeit.bspline.refactor.BSpline2D.fit` was called with a
    :class:`numpy.ndarray` ``basis`` argument (i.e., when ``funcname`` is
    ``None`` after fitting).  Storing the basis enables evaluation via
    :meth:`~pypeit.bspline.refactor.BSpline2D.value` without re-supplying it
    after loading from a file.  When the fit was performed with a string basis
    (e.g. ``'legendre'``), ``basis`` is ``None`` and the caller must supply
    ``basis_x`` to :meth:`~pypeit.bspline.refactor.BSpline2D.value`.

    Parameters
    ----------
    x : :class:`numpy.ndarray`, optional
        Independent variable used to build the initial breakpoint grid.
    knots : :class:`~pypeit.bspline.refactor.Knots` or :class:`numpy.ndarray` or None, optional
        Knot specification.
    nord : int, optional
        Order of the b-spline (default 4).

    See Also
    --------
    pypeit.bspline.refactor.BSpline2D
    """

    version = '2.0.0'

    datamodel = {
        'nord': dict(
            otype=int,
            descr='Order of the b-spline.'
        ),
        'bkpt_full': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Full padded knot vector stored for serialisation.'
        ),
        'bkpt_gpm': dict(
            otype=np.ndarray,
            atype=np.bool_,
            descr='Active-breakpoint mask; True where a breakpoint participates in the fit.'
        ),
        'coeff': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Fitted 2D b-spline coefficients, shape (nc, npoly).'
        ),
        'icoeff': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Inverse-covariance diagonal, shape (nc, npoly).'
        ),
        'npoly': dict(
            otype=int,
            descr='Number of polynomial basis functions in the second dimension.'
        ),
        'xmin': dict(
            otype=float,
            descr='Lower normalisation bound for the polynomial basis coordinate.'
        ),
        'xmax': dict(
            otype=float,
            descr='Upper normalisation bound for the polynomial basis coordinate.'
        ),
        'funcname': dict(
            otype=str,
            descr='Polynomial family name (e.g. ``legendre``); None when fit was '
                  'called with an array basis.'
        ),
        'basis': dict(
            otype=np.ndarray,
            atype=np.floating,
            descr='Polynomial basis matrix, shape (N, npoly); stored only when fit '
                  'was called with a numpy array basis (funcname is None).'
        ),
    }

    internals = ['knots', 'x', 'yfit', '_cached_design', 'P', 'basis_x']

    def __init__(self, x=None, knots=None, nord=4):
        datamodel.DataContainer.__init__(self)
        BSpline2D.__init__(self, x=x, knots=knots, nord=nord)
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints

    @classmethod
    def from_bspline2d(cls, spl):
        """
        Construct a :class:`BSpline2DContainer` from an existing
        :class:`~pypeit.bspline.refactor.BSpline2D`.

        Parameters
        ----------
        spl : :class:`~pypeit.bspline.refactor.BSpline2D`
            Source B-spline.  Its
            :attr:`~pypeit.bspline.refactor.BSpline.breakpoints` must not be
            ``None``.

        Returns
        -------
        :class:`BSpline2DContainer`
            A new container wrapping a copy of the spline state.  The ``basis``
            datamodel field is populated only when ``spl.funcname is None`` and
            ``spl.P is not None`` (array-basis fit); otherwise it is ``None``.
        """
        return cls.from_dict(d={
            'nord': spl.nord,
            'bkpt_full': np.copy(spl.knots.breakpoints),
            'bkpt_gpm': np.copy(spl.bkpt_gpm),
            'coeff': None if spl.coeff is None else np.copy(spl.coeff),
            'icoeff': None if spl.icoeff is None else np.copy(spl.icoeff),
            'npoly': spl.npoly,
            'xmin': None if spl.xmin is None else float(spl.xmin),
            'xmax': None if spl.xmax is None else float(spl.xmax),
            'funcname': spl.funcname,
            'basis': (None if spl.funcname is not None or spl.P is None
                      else np.copy(spl.P)),
        })

    def copy(self):
        """
        Return a deep copy as a new :class:`BSpline2DContainer`.

        The design-matrix cache is not copied (the copy starts with a cold
        cache).  The polynomial basis ``P`` is copied when it exists: for
        array-basis fits (``funcname is None``) it is already restored from
        :attr:`basis` by :meth:`_validate`; for string-basis fits it is copied
        explicitly.

        Returns
        -------
        :class:`BSpline2DContainer`
        """
        new = type(self).from_bspline2d(self)
        new.x = self.x
        new.basis_x = self.basis_x
        new.yfit = None if self.yfit is None else np.copy(self.yfit)
        # _validate() restores P from basis when funcname is None;
        # copy P explicitly for the string-basis case.
        if self.funcname is not None and self.P is not None:
            new.P = np.copy(self.P)
        return new

    def fit(
        self, x, y, ivar=None, basis='legendre', npoly=1, basis_x=None,
        xmin=None, xmax=None, reset_knots=False
    ):
        """
        Fit a weighted least-squares 2D B-spline.

        This delegates entirely to
        :meth:`~pypeit.bspline.refactor.BSpline2D.fit`, but updates
        :attr:`bkpt_full`, as necessary.  See documentation for
        :meth:`~pypeit.bspline.refactor.BSpline2D.fit`.
        """
        result = super().fit(
            x, y, ivar=ivar, basis=basis, npoly=npoly, basis_x=basis_x,
            xmin=xmin, xmax=xmax, reset_knots=reset_knots
        )
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints
        return result

    def _validate(self):
        if self.bkpt_full is not None and self.knots is None:
            self.knots = Knots(full=self.bkpt_full, nord=self.nord)
            if self.bkpt_gpm is None:
                self.bkpt_gpm = np.ones(self.knots.breakpoints.size, dtype=bool)
        if self.basis is not None and self.P is None:
            self.P = self.basis

    def _bundle(self):
        if self.knots is not None and self.knots.breakpoints is not None:
            self.bkpt_full = self.knots.breakpoints
        if self.funcname is None and self.P is not None:
            self.basis = self.P
        return super()._bundle(ext='BSPLINE2D')
