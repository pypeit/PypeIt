""" Module for fitting codes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
# TODO -- Consider moving the Object out of core

import numpy as np
import inspect
from matplotlib import pyplot as plt


from scipy.optimize import curve_fit

from pypeit import log
from pypeit import PypeItError
from pypeit.core import bspline
from pypeit.core import pydl
from pypeit.datamodel import DataContainer

from IPython import embed


class PypeItFitCollection:
    """
    A collection of 1D fits to a set of data.

    The provided data to be fit (``xpos``, ``ypos``) can be provided as 1D or 2D
    arrays, but their shape must match.

        - If 1D, only one fit is performed, and this is effectively identical to
          a single instance of :class:`~pypeit.core.fitting.PypeItFit`.

        - If 2D, fits are performed along the 2nd axis; i.e., a model is fit to
          ``(xpos[0],ypos[0])`` vectors, then to the ``(xpos[1],ypos[1])``
          vectors, etc.

    This class uses :func:`~pypeit.core.fitting.robust_fit` to perform all the
    fits.

    Parameters
    ----------
    xpos : :class:`numpy.ndarray`
        The x positions of the data to be fit.  Can be 1D or 2D; if the latter,
        fits are performed along the 2nd axis (see class description).
    ypos : :class:`numpy.ndarray`
        The y positions of the data to be fit.  Must have the same shape as
        ``xpos``.
    ivar : :class:`numpy.ndarray`, optional
        The inverse variance in the ``ypos`` data.  Must have the same shape as
        ``ypos``.
    gpm : :class:`numpy.ndarray`, optional
        Good pixel mask.  Must have the same shape as ``ypos``.  If None and
        ``ivar`` is None, all pixels are considered good.  If None and ``ivar``
        is provided, the data with ``ivar > 0`` are considered good.
    func : str, optional
        The functional form to use for the fit.  Must be one of
        'polynomial', 'legendre', or 'chebyshev'.
    order : int, optional
        The order of the polynomial to be fit.
    xmin : float, optional
        The minimum x value to be used for the fit.  If None, this is set to the
        minimum of ``xpos`` (i.e., the *entire* array).  Generally, you should
        *not* provide this, and just let the code determine it from ``xpos``.
    xmax : float, optional
        The maximum x value to be used for the fit.  If None, this is set to the
        maximum of ``xpos`` (i.e., the *entire* array).  Generally, you should
        *not* provide this, and just let the code determine it from ``xpos``.
    maxiter : :obj:`int`, optional
        Maximum number of rejection iterations; see
        :func:`~pypeit.core.fitting.robust_fit`.
    maxdev : :obj:`int`, :obj:`float`, optional
        An absolute-difference threshold for rejecting outliers; see
        :func:`~pypeit.core.fitting.robust_fit`.
    lower : :obj:`int`, :obj:`float`, optional
        A sigma-rejection threshold for data with values less than the model;
        see :func:`~pypeit.core.fitting.robust_fit`.
    upper : :obj:`int`, :obj:`float`, optional
        A sigma-rejection threshold for data with values greater than the model;
        see :func:`~pypeit.core.fitting.robust_fit`.
    """

    allowed_functions = ['polynomial', 'legendre', 'chebyshev']
    """
    Allowed functional forms for fitting.
    """

    def __init__(
        self, xpos, ypos, ivar=None, gpm=None, func='legendre', order=3, xmin=None, xmax=None,
        maxiter=10, maxdev=None, lower=None, upper=None
    ):
        
        self.xpos = xpos
        self.nfit = xpos.shape[0]

        if ypos.shape != self.xpos.shape:
            raise PypeItError(
                'Shape of the ypos array must match the xpos array in PypeItFitCollection.'
            )
        self.ypos = ypos

        if ivar is not None and ivar.shape != self.xpos.shape:
            raise PypeItError(
                'Shape of the ivar array must match the xpos array in PypeItFitCollection.'
            )
        self.ivar = ivar 

        if gpm is None:
            self.gpm = (
                np.ones(self.ypos.shape, dtype=bool)
                if self.ivar is None
                else self.ivar > 0.0
            )
        else:
            self.gpm = gpm
        if self.gpm.shape != self.xpos.shape:
            raise PypeItError(
                'Shape of the gpm array must match the xpos array in PypeItFitCollection.'
            )

        self.func = func
        self.order = order
        self.xmin = np.min(self.xpos) if xmin is None else xmin
        self.xmax = np.max(self.xpos) if xmax is None else xmax

        self.maxdev = maxdev
        self.maxiter = maxiter
        self.lower = lower
        self.upper = upper

        self.coeff = np.zeros((self.nfit, self.order+1), dtype=float)
        self.out_gpm = self.gpm.copy()
        self.xnorm = scale_minmax(self.xpos, minx=self.xmin, maxx=self.xmax)[0]
        self.yfit = np.zeros(self.ypos.shape, dtype=self.ypos.dtype)
        self.pypeitFits = [None]*self.nfit
        for i in range(self.nfit):

            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # TODO: The use of xnorm below IS A BUG!!  However, this reproduces
            # the behavior of the old TraceSet class.  We need to fix this, but
            # it will likely cause havoc with our tests, and the default order
            # we use for edge tracing.
            # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            self.pypeitFits[i] = robust_fit(
                self.xnorm[i], self.ypos[i], self.order,
                function=self.func, maxiter=self.maxiter,
                in_gpm=self.gpm[i], invvar=None if self.ivar is None else self.ivar[i],
                lower=self.lower, upper=self.upper, minx=self.xmin, maxx=self.xmax,
                maxdev=self.maxdev, grow=0, use_mad=False, sticky=False
            )

            self.yfit[i] = self.pypeitFits[i].eval(self.xnorm[i])
            self.coeff[i] = self.pypeitFits[i].fitc
            self.out_gpm[i] = self.pypeitFits[i].gpm

    def eval(self, xpos=None, copy=True):
        """
        Evaluate the fits at the provided coordinates.

        Parameters
        ----------
        xpos : array-like, optional
            Positions at which to evaluate the fits.  This can be 1D or 2D.  If
            1D, the same positions will be used for all fits.  If 2D, the length
            of the first axis *must* be the same as :attr:`nfit`.  If None, the
            x positions are the same as :attr:`xpos` and the fit values are
            identically :attr:`yfit`.
        copy : bool, optional
            Only relevant if ``xpos`` is None.  If True, the returned x and y
            positions are copies of the internal attributes; otherwise, the
            attributes themselves are returned.

        Returns
        -------
        x : :class:`numpy.ndarray`
            Sampled x positions
        y : :class:`numpy.ndarray`
            Evaluated fits at the sampled x positions
        """
        if xpos is None:
            return (
                self.xpos.copy() if copy else self.xpos,
                self.yfit.copy() if copy else self.yfit
            )
        if xpos.ndim == 1:
            _xpos = np.tile(xpos, (self.nfit, 1))
        else:
            if xpos.shape[0] != self.nfit:
                raise PypeItError(
                    f'First axis of a 2D xpos array must be {self.nfit}, not {xpos.shape[0]}.'
                )
            _xpos = xpos

        # TODO: When we fix the use of xnorm in the fit call above, we need to
        # fix it here, as well. 
        _xnorm = scale_minmax(_xpos, minx=self.xmin, maxx=self.xmax)[0]
        return _xpos, np.vstack([self.pypeitFits[i].eval(_xnorm[i]) for i in range(self.nfit)])


class PypeItFit(DataContainer):
    """
    General fitting class used by PypeIt.

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_pypeitfit.rst

    When written to an output-file HDU, all `numpy.ndarray`_ elements are
    bundled into an `astropy.io.fits.BinTableHDU`_, and the other elements are
    written as header keywords.  Any datamodel elements that are None are *not*
    included in the output.

    """

    # Set the version of this class
    version = '1.0.0'

    datamodel = {'xval': dict(otype=np.ndarray, atype=np.floating, descr='x inputs'),
                 'yval': dict(otype=np.ndarray, atype=np.floating, descr='y inputs'),
                 'order': dict(otype=np.ndarray, atype=np.integer,
                               descr='The order of the polynomial to be used in the fitting. '
                                     'This is a 2d array for 2d fits'),
                 'x2': dict(otype=np.ndarray, atype=np.floating,
                            descr='x2 inputs, second independent variable'),
                 'weights': dict(otype=np.ndarray, atype=np.floating, descr='Weights.  Often the same as invvar'),
                 'fitc': dict(otype=np.ndarray, atype=np.floating, descr='Fit coefficients'),
                 'fitcov': dict(otype=np.ndarray, atype=np.floating,
                                descr='Covariance of the coefficients'),
                 # TODO: Can we make this boolean?
                 'gpm': dict(otype=np.ndarray, atype=np.integer, descr='Mask (1=good)'),
                 'success': dict(otype=int,
                                 descr='Flag indicating whether fit was successful (success=1) '
                                       'or if it failed (success=0)'),
                 'func': dict(otype=str,
                              descr='Fit function (polynomial, legendre, chebyshev, polynomial2d,'
                                    ' legendre2d)'),
                 'minx': dict(otype=float,
                              descr='minimum value in the array (or the left limit for a '
                                    'legendre / chebyshev polynomial)'),
                 'maxx': dict(otype=float,
                              descr='maximum value in the array (or the right limit for a '
                                    'legendre / chebyshev polynomial)'),
                 'minx2': dict(otype=float,
                               descr='Same as minx for the second independent variable x2'),
                 'maxx2': dict(otype=float,
                               descr='Same as maxx for the second independent variable x2')}

    # This needs to contain all datamodel items
    # TODO: It depends on how you use it, but the above statement isn't
    # strictly true; see, e.g., TracePCA as one example.
    def __init__(self, xval=None, yval=None, order=None, x2=None, weights=None, fitc=None,
                 fitcov=None, func=None, minx=None, maxx=None, minx2=None,
                 maxx2=None, gpm=None, success=None):
        # Setup the DataContainer
        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = {k: values[k] for k in args[1:]}
        # Init
        super(PypeItFit, self).__init__(d=_d)

    def _bundle(self, ext='PYPEITFIT'):
        """
        Bundle the data in preparation for writing to a fits file.

        See :func:`pypeit.datamodel.DataContainer._bundle`. Data is
        always written to a 'PYPEITFIT' extension.
        """
        return super(PypeItFit, self)._bundle(ext=ext)

    def to_hdu(self, **kwargs):
        """
        Over-ride :func:`pypeit.datamodel.DataContainer.to_hdu` to force to
        a BinTableHDU

        See that func for Args and Returns
        """
        if 'force_to_bintbl' in kwargs and not kwargs['force_to_bintbl']:
            log.warning('PypeItFits objects must always be forced to a BinaryTableHDU for writing.')
        kwargs['force_to_bintbl'] = True
        return super(PypeItFit, self).to_hdu(**kwargs)

    @property
    def bool_gpm(self):
        """
        Generate a bool version of gpm which is int
        for I/O

        Returns:
            `numpy.ndarray`_ or None: bool version of self.gpm or None

        """
        return self.gpm.astype(bool) if self.gpm is not None else None

    def fit(self):
        """
        Perform the fit, either in 1D or 2D depending on the
        data and model.

        Returns:
            int: Flag indicating whether fit was successful (1) or if it failed (0)
        """

        # Init
        self.fitcov = None

        # If the user provided an gpm apply it. The logic below of evaluating the fit only at the non-masked
        # pixels is preferable to the other approach of simply setting the weights to zero. The reason for that is that
        # the fits use a least-square optimization approach using matrix algebra, and lots of zero weights are
        # 1) more costly, and 2) will not produce exactly the same result (due to roundoff error) as actually
        # removing the locations you want to mask.

        # This block ensures sensible zero coefficient outputs are returned if the fits was successful
        if self.bool_gpm is not None and not np.any(self.bool_gpm):
            if self.func == "gaussian":
                self.fitc = np.zeros(self.order[0]).astype(float)
            elif '2d' in self.func:
                self.fitc = np.zeros(self.order[0] + 1, self.order[1] + 1).astype(float)
            else:
                self.fitc = np.zeros(self.order[0] + 1).astype(float)
            log.warning('Input gpm is masked everywhere. Fit is probably probelmatic')
            self.success = 0
            return self.success

        if self.bool_gpm is not None:
            x_out = self.xval[self.bool_gpm]
            y_out = self.yval[self.bool_gpm]
            if self.x2 is not None:
                x2_out = self.x2[self.bool_gpm]
            else:
                x2_out = None
            if self.weights is not None:
                w_out = self.weights[self.bool_gpm]
            else:
                w_out = None
        else:
            x_out = self.xval
            y_out = self.yval
            if self.x2 is not None:
                x2_out = self.x2
            else:
                x2_out = None
            if self.weights is not None:
                w_out = self.weights
            else:
                w_out = None

        # For two-d fits x = x, y = x2, y = z
        if ('2d' in self.func) and (x2_out is not None):
            # Is this a 2d fit?
            self.fitc, self.minx, self.maxx, self.minx2, self.maxx2 = polyfit2d_general(
                x_out, x2_out, y_out, self.order, w=w_out, function=self.func[:-2],
                minx=self.minx, maxx=self.maxx,
                miny=self.minx2, maxy=self.maxx2)
        elif self.func == "polynomial":
            self.fitc = np.polynomial.polynomial.polyfit(
                x_out, y_out, self.order[0], 
                w=np.sqrt(w_out) if w_out is not None else None) # numpy convention
        elif self.func == "legendre" or self.func == "chebyshev":
            xv, self.minx, self.maxx = scale_minmax(x_out, minx=self.minx, maxx=self.maxx)
            self.fitc = np.polynomial.legendre.legfit(xv, y_out, self.order[0], 
                                                      w=np.sqrt(w_out) if w_out is not None else None) \
                if self.func == "legendre" else np.polynomial.chebyshev.chebfit(
                    xv, y_out, self.order[0], 
                    w=np.sqrt(w_out) if w_out is not None else None) # numpy convention
        else:
            raise PypeItError(
                f"Fitting function '{self.func}' is not implemented yet\nPlease choose from "
                "'polynomial', 'legendre', 'chebyshev', 'polynomial2d', 'legendre2d', 'chebyshev2d'"
            )

        self.success = 1
        return self.success

    def eval(self, x, x2=None):
        """
        Return the evaluated fit at locations x
        (and x2, if 2D)

        Args:
            x (`numpy.ndarray`_):
            x2 (`numpy.ndarray`_, optional):
                For 2D fits

        Returns:
            `numpy.ndarray`_:

        """
        return evaluate_fit(self.fitc, self.func, x, x2=x2, minx=self.minx,
                            maxx=self.maxx, minx2=self.minx2, maxx2=self.maxx2)

    def calc_fit_rms(self, apply_mask=True, x2=None):
        """ Simple RMS calculation for the fit on the data.

        Args:
            apply_mask (bool, optional):
                If true, apply mask to data before calculating RMS.
            x2 (`numpy.ndarray`_, optional):
                x locations for 2D fits

        Returns:
            float: Root mean square

        """
        msk = self.bool_gpm
        
        if self.weights is None:
            weights = np.ones(self.xval.size)
        else:
            weights = self.weights
        if apply_mask:
            xval = self.xval[msk]
            yval = self.yval[msk]
            weights = weights[msk]
            x2_val = x2[msk] if x2 is not None else None
        else:
            xval = self.xval.copy()
            yval = self.yval.copy()
            x2_val = x2
        # Normalise
        weights /= np.sum(weights)
        values = self.eval(xval, x2=x2_val)
        # RMS
        return np.sqrt(np.sum(weights * (yval - values) ** 2))


def evaluate_fit(fitc, func, x, x2=None, minx=None,
                 maxx=None, minx2=None, maxx2=None):
    """
    Return the evaluated fit at the x locations

    Args:
        fitc (`numpy.ndarray`_):
            Fit coefficients
        func (str):
            Name of the functional form to fit
        x (`numpy.ndarray`_):
            x locations for the evaluation
        x2 (`numpy.ndarray`_, optional):
            x2 locations for 2D fits
        minx (float, optional):
            Minimum x value for the fit used to normalise the x values
        maxx (float, optional):
            Maximum x value for the fit used to normalise the x values
        minx2 (float, optional):
            Minimum x value for the fit used to normalise the x2 values
        maxx2 (float, optional):
            Maximum x value for the fit used to normalise the x2 values

    Returns:
        `numpy.ndarray`_:  Evaluated fit at the x (and x2) locations

    """
    if func is None:
        return None
    # For two-d fits x = x, y = x2, y = z
    if ('2d' in func) and (x2 is not None):
        # Is this a 2d fit?
        if func[:-2] == "polynomial":
            return np.polynomial.polynomial.polyval2d(x, x2, fitc)
        elif func[:-2] in ["legendre", "chebyshev"]:
            # Scale x-direction
            xv, _, _ = scale_minmax(x, minx=minx, maxx=maxx)
            # Scale x2-direction
            x2v, _, _ = scale_minmax(x2, minx=minx2, maxx=maxx2)
            return (np.polynomial.legendre.legval2d(xv, x2v, fitc) if func[:-2] == "legendre"
                    else np.polynomial.chebyshev.chebval2d(xv, x2v, fitc))
        else:
            raise PypeItError("Function {0:s} has not yet been implemented for 2d fits".format(func))
        # TODO: Why is this return here?  The code will never reach this point
        # because of the if/elif/else above.  What should the behavior be, raise
        # an exception or return None?
        return None
    elif func == "polynomial":
        return np.polynomial.polynomial.polyval(x, fitc)
    elif func == "legendre" or func == "chebyshev":
        xv, _, _ = scale_minmax(x, minx=minx, maxx=maxx)
        return (np.polynomial.legendre.legval(xv, fitc) if func == "legendre"
                else np.polynomial.chebyshev.chebval(xv, fitc))
    else:
        raise PypeItError(
            f"Fitting function '{func}' is not implemented yet\nPlease choose from "
            "'polynomial', 'legendre', 'chebyshev', 'polynomial2d', 'legendre2d', 'chebyshev2d'"
        )


def robust_fit(xarray, yarray, order, x2=None, function='polynomial',
               minx=None, maxx=None, minx2=None, maxx2=None,
               maxiter=10, in_gpm=None, weights=None, invvar=None,
               lower=None, upper=None, maxdev=None, maxrej=None, groupdim=None,
               groupsize=None, groupbadpix=False, grow=0, sticky=True, use_mad=True,
               verbose=True, show_fit=False):
    """
    A robust fit is performed to the xarray, yarray pairs ``mask[i] = 1`` are
    good values, if provided.

    The underlying method(s) are the numpy fitting routines,
    e.g. polyfit, legfit.

    Args:
        xarray (`numpy.ndarray`_):
            independent variable values
        yarray (`numpy.ndarray`_):
            dependent variable values
        order (:obj:`int` or numpy.ndarray`_):
            the order of the polynomial to be used in the fitting. This is an integer for 1d fits and must
            be a tuple or 2d array for 2d fits (i.e. using x2 as the second independent variable).
        x2  (`numpy.ndarray`_, optional):
            Do a 2d fit? This is the second independent variable for 2d fits.
        function (str):
            which function should be used in the fitting.
            (valid inputs are:
            'polynomial', 'legendre', 'chebyshev', 'polynomial2d', 'legendre2d')
        minx (float, optional):
            minimum value in the array (or the left limit for a
            legendre/chebyshev polynomial)
        maxx (float, optional):
            maximum value in the array (or the right limit for a
            legendre/chebyshev polynomial)
        minx2 (float, optional):
            Same as minx for second independent variable x2.
        maxx2 (float, optional):
            Same as maxx for second independent variable x2.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations, default 10.  Set
            this to zero to disable rejection and simply do a fit.
        in_gpm (`numpy.ndarray`_, optional):
            Input mask.  Bad points are marked with a value that
            evaluates to ``False``.  Must have the same number of
            dimensions as ``data``. Points masked as bad "False" in the
            inmask will also always evaluate to "False" in the outmask.
        invvar (:obj:`float`, `numpy.ndarray`_, optional):
            Inverse variance of the data, used to reject points based on
            the values of `upper` and `lower`.  This can either be a
            single float for the entire yarray or a ndarray with the
            same shape as the yarray.
        weights (`numpy.ndarray`_, optional): shape same as xarray and yarray
            If input the code will do a weighted fit. If not input, the
            code will use invvar as the weights. If both invvar and
            weights are input. The fit will be done with weights, but
            the rejection will be based on::

                chi = (data-model) * np.sqrt(invvar)

        lower (:obj:`int` or :obj:`float`, optional):
            If set, reject points with ``data < model - lower * sigma``,
            where ``sigma = 1.0/sqrt(invvar)``.
        upper (:obj:`int` or :obj:`float`, optional):
            If set, reject points with ``data > model + upper * sigma``,
            where ``sigma = 1.0/sqrt(invvar)``.
        maxdev (:obj:`int` or :obj:`float`, optional):
            If set, reject points with ``abs(data-model) > maxdev``.  It is
            permitted to set all three of `lower`, `upper` and `maxdev`.
        maxrej (:obj:`int`, :obj:`numpy.ndarray`, optional):
            Maximum number of points to reject in this iteration.  If
            `groupsize` or `groupdim` are set to arrays, this should be
            an array as well.
        groupdim (:obj:`int`, optional):
            Dimension along which to group the data; set to 1 to group
            along the 1st dimension, 2 for the 2nd dimension, etc.  If
            data has shape ``[100,200]``, then setting ``GROUPDIM=2`` is
            equivalent to grouping the data with ``groupsize=100``.  In
            either case, there are 200 groups, specified by ``[*,i]``. NOT
            WELL TESTED IN PYTHON!
        groupsize (:obj:`int`, optional):
            If this and maxrej are set, then reject a maximum of maxrej
            points per group of groupsize points.  If groupdim is also
            set, then this specifies sub-groups within that. NOT WELL
            TESTED IN PYTHON!!
        groupbadpix (:obj:`bool`, optional):
            If set to ``True``, consecutive sets of bad pixels are
            considered groups, overriding the values of `groupsize`.
        grow (:obj:`int`, optional, default = 0):
            If set to a non-zero integer, N, the N nearest neighbors of
            rejected pixels will also be rejected.
        sticky (:obj:`bool`, optional, default is True):
            If set to ``True``, pixels rejected in one iteration remain
            rejected in subsequent iterations, even if the model
            changes. If
        use_mad (:obj:`bool`, optional, default = False):
            It set to ``True``, compute the median of the maximum
            absolute deviation between the data and use this for the
            rejection instead of the default which is to compute the
            standard deviation of the yarray - modelfit. Note that it is
            not possible to specify use_mad=True and also pass in values
            invvar, and the code will return an error if this is done.
        verbose (:obj:`bool`, optional, default = True):
            If set to ``True``, increase the verbosity to print additional messages.
        show_fit (:obj:`bool`, optional, default = False):
            If set to ``True``, show the final fit and the data in a plot.
            This is only for debugging purposes and should not be used regularly.

    Returns:
        PypeItFit or None:
            Object containing the inputs to the fit and the
            fit itself
    """

    # Setup the initial mask
    if in_gpm is None:
        in_gpm = np.ones(xarray.size, dtype=bool)

    if weights is None:
        if invvar is not None:
            weights = np.copy(invvar)
        else:
            weights = np.ones(xarray.size, dtype=float)

    # Iterate, and mask out new values on each iteration
    iIter = 0
    qdone = False
    this_gpm = np.copy(in_gpm)
    #mskcnt = np.sum(this_gpm)
    #pypeitFit = None
    while (not qdone) and (iIter < maxiter):
        if np.sum(this_gpm) <= np.sum(order) + 1:
            log.warning("More parameters than data points - fit might be undesirable")
        if not np.any(this_gpm):
            log.warning("All points were masked. Returning current fit and masking all points. Fit is likely undesirable")
        pypeitFit = PypeItFit(xval=xarray.astype(float), yval=yarray.astype(float),
                              func=function, order=np.atleast_1d(order),
                              x2=x2.astype(float) if x2 is not None else x2,
                              weights=weights.astype(float), gpm=this_gpm.astype(int),
                              minx=float(minx) if minx is not None else minx,
                              maxx=float(maxx) if maxx is not None else maxx,
                              minx2=float(minx2) if minx2 is not None else minx2,
                              maxx2=float(maxx2) if maxx2 is not None else maxx2)
        pypeitFit.fit()
        ymodel = pypeitFit.eval(xarray, x2=x2)
        # TODO Add nrej and nrej_tot as in robust_optimize below?
        this_gpm, qdone = pydl.djs_reject(yarray, ymodel, outmask=this_gpm, inmask=in_gpm, invvar=invvar,
                                          lower=lower, upper=upper, maxdev=maxdev, maxrej=maxrej,
                                          groupdim=groupdim, groupsize=groupsize, groupbadpix=groupbadpix, grow=grow,
                                          use_mad=use_mad, sticky=sticky)
        # Update the iteration
        iIter += 1
    if (iIter == maxiter) & (maxiter != 0) & verbose:
        log.warning(f'Maximum number of iterations maxiter={maxiter} reached in robust_polyfit_djs')

    # Do the final fit
    pypeitFit = PypeItFit(xval=xarray.astype(float), yval=yarray.astype(float),
                          func=function, order=np.atleast_1d(order),
                          x2=x2.astype(float) if x2 is not None else x2,
                          weights=weights.astype(float), gpm=this_gpm.astype(int),
                          minx=float(minx) if minx is not None else minx,
                          maxx=float(maxx) if maxx is not None else maxx,
                          minx2=float(minx2) if minx2 is not None else minx2,
                          maxx2=float(maxx2) if maxx2 is not None else maxx2)
    pypeitFit.fit()

    if show_fit and x2 is None:
        plt.figure(figsize=(10, 6))
        plt.plot(xarray, yarray, 'ko', label='Data')
        plt.plot(xarray[this_gpm], yarray[this_gpm], 'go', label='Good data')
        plt.plot(xarray[~this_gpm], yarray[~this_gpm], 'ro', label='Rejected data')
        x_fit = np.linspace(np.min(xarray), np.max(xarray), 1000)
        y_fit = pypeitFit.eval(x_fit)
        plt.plot(x_fit, y_fit, 'b-', label='Fit')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.legend()
        plt.show()
    # Return
    return pypeitFit


def robust_optimize(ydata, fitfunc, arg_dict, maxiter=10, inmask=None, invvar=None,
                    lower=None, upper=None, maxdev=None, maxrej=None, groupdim=None,
                    groupsize=None, groupbadpix=False, grow=0, sticky=True, use_mad=False,
                    verbose=False,
                    **kwargs_optimizer):
    """
    A routine to perform robust optimization. It is completely analogous
    to :func:`robust_fit`, but is more general in that it allows
    one to fit a more general model using the optimizer of the users
    choice. If you are fitting simple functions like Chebyshev or
    Legednre polynomials using a linear least-squares algorithm, you
    should use :func:`robust_fit` instead of this function.

    Args:
        ydata (`numpy.ndarray`_):
            Data to fit.
        fitfunc (callable):
            The callable object used to perform the fitting.  The
            calling sequence must be::

                ret_tuple = fitfunc(ydata, inmask, arg_dict, **kwargs_optimizer)

            See the descriptions of `ydata`, `inmask`, `arg_dict`, and
            `kwargs_optimizer`.  The returned object ret_tuple that can
            have two or three elements.  If it has two elements (result,
            ymodel):

                - `result`: Object returned by the specific
                  scipy.optimize method used to perform the fit.
                - `ymodel`: A `numpy.ndarray` with the model fit to
                  `ydata` and with the same shape.

            If it has three elements (result, ymodel, newivar):

                - `newivar`: new inverse variance for the ydata ymodel
                  comparison, in other words chi = (ydata -
                  ymodel)*np.sqrt(newivar). This functionality allows
                  one to deal with cases where the noise of the
                  data-model comaprison is model dependent.

        arg_dict (:obj:`dict`):
            Dictionary containing the other variables needed to evaluate
            the model fit.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations.  Set this to zero to
            disable rejection and simply do a fit.
        inmask (`numpy.ndarray`_, optional):
            Input mask.  Bad points are marked with a value that
            evaluates to `False`.  Must have the same number of
            dimensions as `ydata`.  Points masked as `False` in `inmask`
            will also always evaluate to `False` in the output mask.
        invvar (:obj:`float`, `numpy.ndarray`_, optional):
            Inverse variance of the data, used to reject points based on
            the values of `upper` and `lower`.  This can either be a
            single float for the entire yarray or a ndarray with the
            same shape as the yarray.
        lower (:obj:`int`, :obj:`float`, optional):
            If set, reject points with ``data < model - lower * sigma``, where
            ``sigma = 1/sqrt(invvar)``
        upper (:obj:`int`, :obj:`float`, optional):
            If set, reject points with ``data > model + upper * sigma``, where
            ``sigma = 1/sqrt(invvar)``.
        maxdev (:obj:`int` or :obj:`float`, optional):
            If set, reject points with ``abs(data-model) > maxdev``.  It
            is permitted to set all three of `lower`, `upper` and
            `maxdev`.
        maxrej (:obj:`int`, `numpy.ndarray`_, optional):
            Maximum number of points to reject in this iteration.  If
            `groupsize` or `groupdim` are set to arrays, this should be
            an array, as well.
        groupdim (:obj:`int`, optional):
            Dimension along which to group the data. Set to 1 to group
            along the 1st dimension, 2 for the 2nd dimension, etc.  For
            example, if data has shape [100,200], then setting
            `groupdim=2` is equivalent to grouping the data with
            `groupsize=100`.  In either case, there are 200 groups,
            specified by `[*,i]`.  This functionality is **not well
            tested in python**!
        groupsize (:obj:`int`, optional):
            If this and `maxrej` are set, then reject a maximum of
            `maxrej` points per group of `groupsize` points.  If
            `groupdim` is also set, then this specifies sub-groups
            within that.  This functionality is **not well tested in
            python**!
        groupbadpix (:obj:`bool`, optional):
            If `True`, consecutive sets of bad pixels are considered
            groups, overriding the values of `groupsize`.
        grow (:obj:`int`, optional):
            If set to a non-zero integer, N, the N nearest neighbors of
            rejected pixels will also be rejected.
        sticky (:obj:`bool`, optional):
            If `True`, pixels rejected in one iteration remain rejected
            in subsequent iterations, even if the model changes.
        use_mad (:obj:`bool`, optional):
            It `True`, compute the median of the maximum absolute
            deviation between the data and use this for the rejection
            instead of the default, which is to compute the standard
            deviation of `ydata - modelfit`. Note that it is not
            possible to specify `use_mad=True` and also pass in a value for
            `invvar`, and the code will return an error if this is done.
        **kwargs_optimizer:
            Optional parameters passed to the optimizer.

    Returns:
        tuple:
            - The object returned by the `scipy.optimize` function used
              by the fitter.  See `fitfunc`.
            - A `numpy.ndarray`_ with the model value fit to `ydata` and
              has its same shape.
            - Boolean `numpy.ndarray`_ with the same shape as data
              indicating which pixels were masked in the final fit.
              Convention is that `True` are good values where `False`
              indicates bad values.

    """
    # Setup the initial mask
    if inmask is None:
        inmask = np.ones(ydata.size, dtype=bool)

    nin_good = np.sum(inmask)
    iIter = 0
    qdone = False
    thismask = np.copy(inmask)

    # If init_from_last is not None, the fitfunc will initialize from the previous iteration's fit, which
    # results in signficant speedup for e.g. differential_evolution optimization. Thus
    # init_from_last is None on the first iteration and then is updated in the iteration loop.
    init_from_last = None
    while (not qdone) and (iIter < maxiter):
        ret_tuple = fitfunc(ydata, thismask, arg_dict, init_from_last=init_from_last, **kwargs_optimizer)
        if (len(ret_tuple) == 2):
            result, ymodel = ret_tuple
            invvar_use = invvar
        elif (len(ret_tuple) == 3):
            result, ymodel, invvar_use = ret_tuple
        else:
            raise PypeItError('Invalid return value from fitfunc')
        # Update the
        init_from_last = result
        thismask_iter = thismask.copy()
        thismask, qdone = pydl.djs_reject(ydata, ymodel, outmask=thismask, inmask=inmask, invvar=invvar_use,
                                          lower=lower, upper=upper, maxdev=maxdev, maxrej=maxrej,
                                          groupdim=groupdim, groupsize=groupsize, groupbadpix=groupbadpix, grow=grow,
                                          use_mad=use_mad, sticky=sticky)
        nrej = np.sum(thismask_iter & np.logical_not(thismask))
        nrej_tot = np.sum(inmask & np.logical_not(thismask))
        if verbose:
            log.info(
                'Iteration #{:d}: nrej={:d} new rejections, nrej_tot={:d} total rejections out of ntot={:d} '
                'total pixels'.format(iter, nrej, nrej_tot, nin_good))
        iIter += 1

    if (iIter == maxiter) & (maxiter != 0):
        log.warning('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in robust_optimize')
    outmask = np.copy(thismask)
    if np.sum(outmask) == 0:
        log.warning('All points were rejected!!! The fits will be zero everywhere.')

    # Perform a final fit using the final outmask if new pixels were rejected on the last iteration
    if qdone is False:
        ret_tuple = fitfunc(ydata, outmask, arg_dict, init_from_last=init_from_last, **kwargs_optimizer)

    return ret_tuple + (outmask,)

    #return result, ymodel, outmask


def scale_minmax(x, minx=None, maxx=None):
    """
    Scale in the input array

    Args:
        x (`numpy.ndarray`_): x values to be scaled
        minx (float, optional): Minimum value for scaling
        maxx (float, optional): Maximum value for scaling

    Returns:
        tuple:
            - the scaled x-values in a `numpy.ndarray`_
            - xmin as a float
            - xmax as a float

    """
    xmin = (-1.0 if np.size(x)==1 else np.min(x)) if minx is None else minx
    xmax = ( 1.0 if np.size(x)==1 else np.max(x)) if maxx is None else maxx

    xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
    return xv, xmin, xmax


def moffat(x,p0,p1,p2):
    """
    Moffat profile
    This 3 parameter formulation assumes the trace is known

    Args:
        x (float or `numpy.ndarray`_): x values
        p0 (float): Amplitude
        p1 (float):
          Width scaling
        p2 : float

    Returns:
        float or `numpy.ndarray`_: Evaluated Moffat
    """
    return p0 / (1+(x/p1)**2)**p2


def fit_gauss(x_out, y_out, guesses=None, w_out=None, nparam=3, maxfev=0):
    """
    Fit a 3 or 4 parameter gaussian

    Args:
        x_out (`numpy.ndarray`_):
            x values to be fit
        y_out (`numpy.ndarray`_):
            y values to be fit
        guesses (tuple, optional):
            ampl, cent, sigma, [floor] guesses for the Gaussian; each as floats
        w_out (`numpy.ndarray_`):
            Weights.  1./sqrt(ivar) is expected
        nparam (int, optional):
            Number of parameters in the Gaussian
            Only options are 3 or 4 where the latter includes
            a floor in the fit.
        maxfev (:obj:`int`, optional):
            Maximum number of function evaluations.  Passed directly to
            `scipy.optimize.curve_fit`_.  Note that setting ``maxfev`` to 0 uses
            the default value set by `scipy.optimize.leastsq`_.

    Returns:
        tuple: Fit coefficients, fit covariance from numpy's curve_fit

    """
    if guesses is None:
        ampl, cent, sigma, floor = guess_gauss(x_out, y_out)
    else:
        ampl, cent, sigma, floor = guesses
    # Error
    if w_out is not None:
        sig_y = 1. / w_out
    else:
        sig_y = None
    
    # 3 param values
    p0=[ampl, cent, sigma] 
    func = gauss_3deg

    if nparam == 4:
        p0 = [floor] + p0
        func = gauss_4deg

    return curve_fit(func, x_out, y_out, p0=p0, sigma=sig_y, maxfev=maxfev)


def gauss_3deg(x,ampl,cent,sigm):
    """  Generate a simple 3-parameter Gaussian

    Args:
        x (float or `numpy.ndarray`_): x values
        ampl (float): Amplitude
        cent (float): Centroid
        sigm (float): sigma

    Returns:
        float or `numpy.ndarray`_: Evaluated Gausssian
    """
    return ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)


def gauss_4deg(x,b, ampl,cent,sigm):
    """  Generate a simple 4-parameter Gaussian

    Args:
        x (float or `numpy.ndarray`_): x values
        b (float): Floor
        ampl (float): Amplitude
        cent (float): Centroid
        sigm (float): sigma

    Returns:
        float or `numpy.ndarray`_: Evaluated Gausssian
    """
    return b + ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)


def guess_gauss(x,y):
    """
    Guesses Gaussian parameters with basic stats

    Args:
        x (`numpy.ndarray`_): x-values
        y (`numpy.ndarray`_): y-values

    Returns:
        tuple:  Amplitude, centroid, sigma, floor all as :obj:`float`

    """
    ypos = y - y.min()
    cent = np.sum(ypos*x)/np.sum(ypos)
    sigma = np.sqrt(np.abs(np.sum((x-cent)**2*ypos)/np.sum(ypos))) # From scipy doc
    # Calculate ampl from pixels within +/- sigma/2
    cen_pix= np.abs(x-cent)<sigma/2
    if np.any(cen_pix):
        ampl = np.median(y[cen_pix])
    else:
        ampl = y.max()
    # Floor
    floor = np.median(np.percentile(y,50))
    # Return
    return ampl, cent, sigma, floor


def polyfit2d_general(x, y, z, deg, w=None, function='polynomial',
                      minx=None, maxx=None, miny=None, maxy=None):
    """
    2D Polynomimal fit

    Args:
        x (`numpy.ndarray`_): x-values
        y (`numpy.ndarray`_): y-values
        z (`numpy.ndarray`_): value of data at each (x,y) coordinate
        deg (tuple): degree of polynomial fit in the form [nx,ny]
        w (`numpy.ndarray`_, optional):
            weights.  Often invvar
        function (str, optional):
            2D function to fit.  Options are 'polynomial', 'chebyshev' or 'legendre'
        minx (float, optional):
            Minimum x value for the fit used to normalise the x values
        maxx (float, optional):
            Maximum x value for the fit used to normalise the x values
        miny (float, optional):
            Minimum value for the fit used to normalise the y values
        maxy (float, optional):
            Maximum value for the fit used to normalise the y values

    Returns:
        tuple:
            - The coefficients of the polynomial fit as a `numpy.ndarray`_
            - minx, maxx, miny, maxy: min and max values for the fit as :obj:`float`

    """
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    deg = np.asarray(deg)
    # Vander
    if function == 'polynomial':
        vander = np.polynomial.polynomial.polyvander2d(x, y, deg)
    elif function == 'legendre' or function == 'chebyshev':
        xv, minx, maxx = scale_minmax(x, minx=minx, maxx=maxx)
        yv, miny, maxy = scale_minmax(y, minx=miny, maxx=maxy)
        vander = np.polynomial.legendre.legvander2d(xv, yv, deg) if function == 'legendre' \
            else np.polynomial.chebyshev.chebvander2d(xv, yv, deg)
    else:
        raise PypeItError("Not ready for this type of {:s}".format(function))
    # Weights
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            log.debug("fitting.polyfit2d - Expected 1D vector for weights")
        if len(x) != len(w) or len(y) != len(w) or len(x) != len(y):
            log.debug("fitting.polyfit2d - Expected x, y and weights to have same length")
        z = z * w
        vander = vander * w[:,np.newaxis]
    # Reshape
    vander = vander.reshape((-1,vander.shape[-1]))
    z = z.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, z, rcond=None)[0]
    return c.reshape(deg+1), minx, maxx, miny, maxy


def twoD_Gaussian(tup, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    A 2D Gaussian to be used to fit the cross-correlation

    Args:
        tup (tuple):
            A two element tuple containing the (x,y) coordinates where the 2D Gaussian will be evaluated
        amplitude (float):
            The amplitude of the 2D Gaussian
        xo (float):
            The centre of the Gaussian in the x direction
        yo (float):
            The centre of the Gaussian in the y direction
        sigma_x (float):
            The dispersion of the Gaussian in the x direction
        sigma_y (float):
            The dispersion of the Gaussian in the y direction
        theta (float):
            The angle of the major axis relative to the horizontal
        offset (float):
            Constant additive term

    Returns:
        `numpy.ndarray`_: The value of the 2D Gaussian at the given coordinates
    """
    (x, y) = tup
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()


def bspline_qa(xdata, ydata, sset, gpm, yfit, xlabel=None, ylabel=None, title=None, show=True):
    """
    Construct a QA plot of the B-spline fit.

    .. warning::

        This function has not been tested with
        :class:`~pypeit.core.bspline.BSpline2D` fits.  The call to
        :meth:`~pypeit.core.bspline.BSpline.value` at the breakpoints
        may fail or produce unexpected results for 2D fits.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.  Regardless of shape, data is treated as
        one-dimensional.
    ydata : :class:`numpy.ndarray`
        Dependent variable.  Regardless of shape, data is treated as
        one-dimensional.
    sset : :class:`~pypeit.core.bspline.BSpline`
        Fitted B-spline object, as returned by
        :func:`~pypeit.core.fitting.iterative_bspline_fit`.
    gpm : :class:`numpy.ndarray`
        Boolean array with the same size as ``xdata``.  Points
        rejected during the fit have ``gpm=False``.
    yfit : :class:`numpy.ndarray`
        Best-fitting model sampled at ``xdata``, as returned by
        :func:`~pypeit.core.fitting.iterative_bspline_fit`.
    xlabel : str, optional
        Label for the abscissa.  If None, no label is added.
    ylabel : str, optional
        Label for the ordinate.  If None, no label is added.
    title : str, optional
        Title for the plot.  If None, no title is added.
    show : bool, optional
        If True, display the plot with the legend, axis labels, and
        title applied.  If False, return the
        :class:`~matplotlib.axes.Axes` instance immediately, before
        the legend, axis labels, or title are added.

    Returns
    -------
    ax : :class:`matplotlib.axes.Axes`
        Axes instance with the data, model, and breakpoints plotted.
        Only returned when ``show`` is False.
    """
    goodbk = sset.bkpt_gpm
    bkpt, _ = sset.value(sset.breakpoints[goodbk])
    was_fit_and_masked = np.logical_not(gpm)

    plt.clf()
    ax = plt.gca()
    ax.plot(xdata, ydata, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
            linestyle='None', label='data')
    ax.plot(xdata[was_fit_and_masked], ydata[was_fit_and_masked], color='red', marker='+',
            markersize=1.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
    ax.plot(xdata, yfit, color='cornflowerblue', label='fit')
    ax.plot(sset.breakpoints[goodbk], bkpt, color='lawngreen', marker='o', markersize=2.0,
            mfc='lawngreen', fillstyle='full', linestyle='None', label='bspline breakpoints')
    ax.set_ylim(0.99 * np.amin(yfit), 1.01 * np.amax(yfit))
    if not show:
        return ax

    plt.legend()
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.show()


def iterative_bspline_fit(
    x, y, ivar=None, gpm=None, nord=4, basis=None, npoly=1, basis_x=None, xmin=None, xmax=None,
    kwargs_knots={}, relative=None, upper=5, lower=5, maxiter=25, kwargs_reject={}
):
    """
    Fit a B-spline with iterative sigma-clipping rejection.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        Independent variable.  Should be sorted in ascending order;
        :meth:`~pypeit.core.bspline.BSpline.fit` sorts internally, but
        the default ``stride`` knot strategy places breakpoints at
        evenly-spaced indices of ``x``, which is only meaningful when
        ``x`` is sorted.
    y : :class:`numpy.ndarray`
        Dependent variable.
    ivar : :class:`numpy.ndarray`, optional
        Inverse variance of ``y``.  Points with ``ivar <= 0`` are
        excluded from the fit regardless of ``gpm``.
    gpm : :class:`numpy.ndarray` of bool, optional
        Additional input mask; True for pixels to include in the fit.
        When None (default), all pixels with positive ``ivar`` are
        used; if ``ivar`` is also None, all pixels are used.
    nord : int, optional
        B-spline order.
    basis : str or :class:`numpy.ndarray` or None, optional
        Polynomial basis specification for the second variable.  When a
        string, it must be one of ``'legendre'``, ``'chebyshev'``,
        ``'poly'``, or ``'poly1'``; ``basis_x`` must also be provided.
        When a :class:`numpy.ndarray`, it is used directly as the
        pre-built polynomial basis matrix; a 1D array of size
        ``x.size * npoly`` is reshaped to ``(x.size, npoly)``
        automatically.  When ``None`` (default), a 1D
        :class:`~pypeit.core.bspline.BSpline` fit is performed.
    npoly : int, optional
        Number of polynomial terms; forwarded to
        :meth:`~pypeit.core.bspline.BSpline2D.fit` and ignored when
        ``basis`` is an array or ``None``.
    basis_x : :class:`numpy.ndarray` or None, optional
        Second variable.  Required when ``basis`` is a string; ignored
        when ``basis`` is ``None`` or a pre-built array.
    xmin : float, optional
        Minimum value of ``basis_x`` for normalisation; forwarded to
        :meth:`~pypeit.core.bspline.BSpline2D.fit` and ignored when
        ``basis`` is an array or ``None``.
    xmax : float, optional
        Maximum value of ``basis_x`` for normalisation; forwarded to
        :meth:`~pypeit.core.bspline.BSpline2D.fit` and ignored when
        ``basis`` is an array or ``None``.
    kwargs_knots : dict, optional
        Keyword arguments forwarded to :class:`~pypeit.core.bspline.Knots`
        to control knot placement.  See :class:`~pypeit.core.bspline.Knots`
        for the available strategies and their default values.
    relative : array-like or None, optional
        Index array for computing the relative chi-squared used to
        rescale the rejection thresholds.  When ``None`` (default) the
        thresholds are not rescaled.
    upper : float, optional
        Upper sigma-clipping threshold.
    lower : float, optional
        Lower sigma-clipping threshold.
    maxiter : int, optional
        Maximum number of fit-reject iterations.
    kwargs_reject : dict, optional
        Additional keyword arguments forwarded to
        :func:`~pypeit.core.pydl.djs_reject`.

    Returns
    -------
    bspl : :class:`~pypeit.core.bspline.BSpline`, :class:`~pypeit.core.bspline.BSpline2D`
        Fitted spline object.  A
        :class:`~pypeit.core.bspline.BSpline2D` is returned when
        ``basis`` is provided; a
        :class:`~pypeit.core.bspline.BSpline` is returned otherwise.
        ``None`` is returned on catastrophic failure (exit status 4
        with no valid input points).
    outmask : :class:`numpy.ndarray` of bool
        Final good-pixel mask after all rejection iterations.
    yfit : :class:`numpy.ndarray`
        Best-fit model values at ``x``.
    reduced_chi : float
        Reduced chi-squared of the final fit.
    exit_status : int
        Convergence status code that will be one of the following:

            - 0 — converged normally
            - 1 — maximum iterations reached
            - 2 — too few good points during iteration
            - 3 — singular or degenerate fit (all breakpoints lost)
            - 4 — fewer good points than ``nord`` on entry

    """
    _basis = basis
    if basis is None:
        _npoly = 1
        bspl_cls = bspline.BSpline
    elif isinstance(basis, np.ndarray):
        _basis = np.asarray(basis)
        if basis.ndim == 1:
            nx = x.size
            _npoly = basis.size // nx
            if basis.size != nx * _npoly:
                raise PypeItError('basis array size is not a multiple of x.size.')
            _basis = _basis.reshape(nx, _npoly).copy()
        else:
            _npoly = _basis.shape[1]
        bspl_cls = bspline.BSpline2D
    else:
        _npoly = npoly
        bspl_cls = bspline.BSpline2D

    outmask = np.ones(y.shape, dtype=bool)
    maskwork = outmask.copy()
    if gpm is not None:
        maskwork &= gpm
    if ivar is not None:
        maskwork &= ivar > 0
    if not maskwork.any():
        return None, outmask, np.zeros(y.shape), 0., 4

    bspl = bspl_cls(x=x[maskwork], knots=bspline.Knots(**kwargs_knots), nord=nord)
    if maskwork.sum() < bspl.nord:
        log.warning(
            f'B-spline fit failed: Number of good data points ({maskwork.sum()}) fewer than '
            f'order of the fit ({bspl.nord}).'
        )
        # TODO: Should this return outmask or maskwork?  Seems odd to return a
        # fully True gpm for a failed fit.
        return bspl, outmask, np.zeros(y.shape), 0., 4

    err = -1
    qdone = False
    iiter = 0
    exit_status = 0
    reduced_chi = 0.
    relative_factor = 1.0
    nrel = 0 if relative is None else len(relative)
    _gpm = None if gpm is None else gpm.copy()

    termwidth = 80 - 13
    log.info('B-spline fit:')
    log.info(f'    npoly = {_npoly} profile basis functions')
    log.info(f'    ngood = {maskwork.sum()}/{maskwork.size} measurements')
    log.info(
        f' {"Iter":>4}  {"Chi^2":>8}  {"N Rej":>7}  {"R. Fac":>6} '.center(termwidth, '*')
    )
    log.info(f' {"-"*4}  {"-"*8}  {"-"*7}  {"-"*6} '.center(termwidth))
    nullval = f'  {"--":>8}  {"--":>7}  {"--":>6} '

    while (err != 0 or not qdone) and iiter <= maxiter and exit_status == 0:
        if maskwork.sum() <= 1:
            exit_status = 2
            break

        _ivar = maskwork.astype(float)
        if ivar is not None:
            _ivar *= ivar

        if _basis is None:
            err, yfit = bspl.fit(x=x, y=y, ivar=_ivar)
        else:
            err, yfit = bspl.fit(
                x=x, y=y, ivar=_ivar, basis=_basis, basis_x=basis_x,
                npoly=_npoly, xmin=xmin, xmax=xmax,
            )

        # WARNING: Because we count the iteration here, before assessing err,
        # loop iterations include fits that are redone because of breakpoints
        # being masked.
        # TODO: Consider altering `BSpline(2D).fit` such that it accepts a
        # number of iterations it is allowed to perform while masking bad
        # breakpoints.
        iiter += 1

        if err == -2:
            log.warning('B-spline fit failed: All break points lost!')
            return bspl, np.zeros(x.shape, dtype=bool), np.zeros(x.shape), reduced_chi, 3
        if err != 0:
            log.info(f' {iiter:4d}{nullval}'.center(termwidth))
            continue

        ngood = maskwork.sum()
        goodbk_count = bspl.bkpt_gpm[bspl.nord:].sum()
        chi_array = (y - yfit) * np.sqrt(_ivar)
        reduced_chi = np.sum(chi_array**2) / (ngood - _npoly * (goodbk_count + bspl.nord) - 1)
        if relative is not None:
            if nrel == 1:
                this_chi2 = reduced_chi
            else:
                this_chi2 = np.sum(chi_array[relative]**2) / (
                    nrel - (goodbk_count + bspl.nord) - 1
                )
            relative_factor = max(np.sqrt(this_chi2), 1.0)
        maskwork, qdone = pydl.djs_reject(
            y, yfit, invvar=_ivar, inmask=_gpm, outmask=maskwork,
            upper=upper * relative_factor, lower=lower * relative_factor,
            **kwargs_reject
        )
        _gpm = np.copy(maskwork)
        log.info((
            f' {iiter:4d}  {reduced_chi:8.3f}  {maskwork.sum():7d}  {relative_factor:6.2f} '
        ).center(termwidth))

    log.info((
        f' {"DONE":>4}  {reduced_chi:8.3f}  {maskwork.sum():7d}  {relative_factor:6.2f} '
    ).center(termwidth))

    if iiter == maxiter + 1:
        exit_status = 1
    return bspl, maskwork.copy(), yfit, reduced_chi, exit_status
