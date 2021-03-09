""" Module for fitting codes

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import numpy as np
import inspect
from matplotlib import pyplot as plt


from scipy.optimize import curve_fit

from pypeit.core import pydl
from pypeit import bspline
from pypeit import msgs
from pypeit.datamodel import DataContainer

from IPython import embed

class PypeItFit(DataContainer):
    # Set the version of this class
    version = '1.0.0'

    datamodel = {'xval': dict(otype=np.ndarray, atype=np.floating, descr='x inputs'),
                 'yval': dict(otype=np.ndarray, atype=np.floating, descr='y inputs'),
                 'order': dict(otype=np.ndarray, atype=np.integer,
                               descr='The order of the polynomial to be used in the fitting. '
                                     'This is a 2d array for 2d fits'),
                 'x2': dict(otype=np.ndarray, atype=np.floating,
                            descr='x2 inputs, second independent variable'),
                 'weights': dict(otype=np.ndarray, atype=np.floating, descr='Weights'),
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
            msgs.warn('PypeItFits objects must always be forced to a BinaryTableHDU for writing.')
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
        Perform the fit
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
            msgs.warn('Input gpm is masked everywhere. Fit is probably probelmatic')
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
            self.fitc = np.polynomial.polynomial.polyfit(x_out, y_out, self.order[0], w=w_out)
        elif self.func == "legendre" or self.func == "chebyshev":
            xv, self.minx, self.maxx = scale_minmax(x_out, minx=self.minx, maxx=self.maxx)
            self.fitc = np.polynomial.legendre.legfit(xv, y_out, self.order[0], w=w_out)  \
                if self.func == "legendre" else np.polynomial.chebyshev.chebfit(xv, y_out, self.order[0], w=w_out)
        else:
            msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                       "Please choose from 'polynomial', 'legendre', 'chebyshev','polynomial2d', 'legendre2d'")

        self.success = 1
        return self.success

    def eval(self, x, x2=None):
        """
        Return the evaluated fit

        Args:
            x (`numpy.ndarray`_, optional):
            x2 (`numpy.ndarray`_, optional):
                For 2D fits

        Returns:
            `numpy.ndarray`_:

        """
        return evaluate_fit(self.fitc, self.func, x, x2=x2, minx=self.minx,
                            maxx=self.maxx, minx2=self.minx2, maxx2=self.maxx2)

    def calc_fit_rms(self, apply_mask=True, x2=None):
        """ Simple RMS calculation

        Args:
            apply_mask (bool, optional):
                Apply mask?
            x2 (`numpy.ndarray`_, optional):
                For 2D fits

        Returns:
            float: RMS

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


def evaluate_fit(fitc, func, x, x2=None, minx=None, maxx=None, minx2=None, maxx2=None):
    """
    Return the evaluated fit

    Args:
        x (`numpy.ndarray`_, optional):
        x2 (`numpy.ndarray`_, optional):
            For 2D fits

    Returns:
        `numpy.ndarray`_:

    """
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
            msgs.error("Function {0:s} has not yet been implemented for 2d fits".format(func))
        return None
    elif func == "polynomial":
        return np.polynomial.polynomial.polyval(x, fitc)
    elif func == "legendre" or func == "chebyshev":
        xv, _, _ = scale_minmax(x, minx=minx, maxx=maxx)
        return (np.polynomial.legendre.legval(xv, fitc) if func == "legendre"
                else np.polynomial.chebyshev.chebval(xv, fitc))
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev', 'polynomial2d', 'legendre2d', 'chebyshev2d'")


def robust_fit(xarray, yarray, order, x2=None, function='polynomial',
               minx=None, maxx=None, minx2=None, maxx2=None,
               maxiter=10, in_gpm=None, weights=None, invvar=None,
               lower=None, upper=None, maxdev=None, maxrej=None, groupdim=None,
               groupsize=None, groupbadpix=False, grow=0, sticky=True, use_mad=True,
               debug=False):
    """
    A robust fit is performed to the xarray, yarray pairs
    ``mask[i] = 1`` are good values.

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
        function:
            which function should be used in the fitting (valid inputs:
            'polynomial', 'legendre', 'chebyshev', 'polynomial2d', 'legendre2d')
        minx:
            minimum value in the array (or the left limit for a
            legendre/chebyshev polynomial)
        maxx:
            maximum value in the array (or the right limit for a
            legendre/chebyshev polynomial)
        minx2:
            Same as minx for second independent variable x2.
        maxx2:
            Same as maxx for second independent variable x2.
        maxiter (:class:`int`, optional):
            Maximum number of rejection iterations, default 10.  Set
            this to zero to disable rejection and simply do a fit.
        in_gpm (:class:`numpy.ndarray`, optional):
            Input mask.  Bad points are marked with a value that
            evaluates to ``False``.  Must have the same number of
            dimensions as `data`. Points masked as bad "False" in the
            inmask will also always evaluate to "False" in the outmask
        invvar (:class:`float`, `numpy.ndarray`, optional):
            Inverse variance of the data, used to reject points based on
            the values of `upper` and `lower`.  This can either be a
            single float for the entire yarray or a ndarray with the
            same shape as the yarray.
        weights (np.ndarray): shape same as xarray and yarray
            If input the code will do a weighted fit. If not input, the
            code will use invvar as the weights. If both invvar and
            weights are input. The fit will be done with weights, but
            the rejection will be based on::

                chi = (data-model) * np.sqrt(invvar)

        lower (:class:`int`, :class:`float`, optional):
            If set, reject points with ``data < model - lower * sigma``,
            where ``sigma = 1.0/sqrt(invvar)``.
        upper (:class:`int`, :class:`float`, optional):
            If set, reject points with ``data > model + upper * sigma``,
            where ``sigma = 1.0/sqrt(invvar)``.
        maxdev (:class:`int`, :class:`float`, optional):
            If set, reject points with ``abs(data-model) > maxdev``.  It is
            permitted to set all three of `lower`, `upper` and `maxdev`.
        maxrej (:class:`int`, :class:`numpy.ndarray`, optional):
            Maximum number of points to reject in this iteration.  If
            `groupsize` or `groupdim` are set to arrays, this should be
            an array as well.
        groupdim (:class:`int`):
            Dimension along which to group the data; set to 1 to group
            along the 1st dimension, 2 for the 2nd dimension, etc.  If
            data has shape ``[100,200]``, then setting ``GROUPDIM=2`` is
            equivalent to grouping the data with ``groupsize=100``.  In
            either case, there are 200 groups, specified by ``[*,i]``. NOT
            WELL TESTED IN PYTHON!
        groupsize (:class:`int`):
            If this and maxrej are set, then reject a maximum of maxrej
            points per group of groupsize points.  If groupdim is also
            set, then this specifies sub-groups within that. NOT WELL
            TESTED IN PYTHON!!
        groupbadpix (:class:`bool`, optional):
            If set to ``True``, consecutive sets of bad pixels are
            considered groups, overriding the values of `groupsize`.
        grow (:class:`int`, optional, default = 0):
            If set to a non-zero integer, N, the N nearest neighbors of
            rejected pixels will also be rejected.
        sticky (:class:`bool`, optional, default is True):
            If set to ``True``, pixels rejected in one iteration remain
            rejected in subsequent iterations, even if the model
            changes. If
        use_mad (:class:`bool`, optional, default = False):
            It set to ``True``, compute the median of the maximum
            absolute deviation between the data and use this for the
            rejection instead of the default which is to compute the
            standard deviation of the yarray - modelfit. Note that it is
            not possible to specify use_mad=True and also pass in values
            invvar, and the code will return an error if this is done.

    Returns:
        PypeItFit or None:
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
            msgs.warn("More parameters than data points - fit might be undesirable")
        if not np.any(this_gpm):
            msgs.warn("All points were masked. Returning current fit and masking all points. Fit is likely undesirable")
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
    if (iIter == maxiter) & (maxiter != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in robust_polyfit_djs')

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

    # Return
    return pypeitFit


def scale_minmax(x, minx=None, maxx=None):
    """
    Scale in the input array

    Args:
        x (`numpy.ndarray`_): x-values
        minx (float, optional): Minimum value for scaling
        maxx (float, optional): Maximum value for scaling

    Returns:
        `numpy.ndarray`_: Scaled x values

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
        x (float or ndarray): x values
        p0 (float): Amplitude
        p1 (float):
          Width scaling
        p2 : float

    Returns:
        float or ndarray: Evaluated Moffat
    """
    return p0 / (1+(x/p1)**2)**p2


def fit_gauss(x_out, y_out, guesses=None, w_out=None):
    """
    Fit a 3 parameter gaussian

    NOTE:  THIS WILL BE REPLACED SOON BY SOMETHING FAST

    Args:
        x_out (`numpy.ndarray_`):
        y_out (`numpy.ndarray_`):
        guesses (tuple, optional):
            ampl, cent, sigma guesses for the Gaussian; each as floats
        w_out (`numpy.ndarray_`):
            Weights

    Returns:
        tuple: Fit coefficients, fit covariance

    """
    if guesses is None:
        ampl, cent, sigma = guess_gauss(x_out, y_out)
    else:
        ampl, cent, sigma = guesses
    # Error
    if w_out is not None:
        sig_y = 1. / w_out
    else:
        sig_y = None
    return curve_fit(gauss_3deg, x_out, y_out, p0=[ampl, cent, sigma], sigma=sig_y)


def gauss_2deg(x,ampl,sigm):
    """
    Simple 2 parameter Gaussian (amplitude, sigma)

    Args:
        x
        ampl
        sigm

    Returns:
        float or ndarray: Evaluated Gausssian
    """
    return ampl*np.exp(-1.*x**2/2./sigm**2)


def gauss_3deg(x,ampl,cent,sigm):
    """  Simple 3 parameter Gaussian

    Args:
        x (float or ndarray): x-valus
        ampl (float): Amplitude
        cent (float): Centroid
        sigm (float): sigma

    Returns:
        float or ndarray: Evaluated Gausssian
    """
    return ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)


def gauss_4deg(x,b, ampl,cent,sigm):
    """  Simple 3 parameter Gaussian

    Args:
        x
        b (float): Floor
        ampl (float): Amplitude
        cent (float): Centroid
        sigm (float): sigma

    Returns:
        float or ndarray: Evaluated Gausssian
    """
    return b + ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)


def gauss_5deg(x,m, b, ampl,cent,sigm):
    """  Simple 3 parameter Gaussian

    Args:
        x
        m (float): Slope of floor
        b (float): Floor
        ampl (float): Amplitude
        cent (float): Centroid
        sigm (float): sigma

    Returns:
        float or ndarray: Evaluated Gausssian
    """
    return b + m*x + ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)


def guess_gauss(x,y):
    """
    Guesses Gaussian parameters with basic stats

    Args:
        x (ndarray): x-values
        y (ndarray): y-values

    Returns:
        float, float, float:  Amplitude, centroid, sigma

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
    # Return
    return ampl, cent, sigma



def polyfit2d_general(x, y, z, deg, w=None, function='polynomial',
                      minx=None, maxx=None, miny=None, maxy=None):
    """
    2D Polynomimal fit

    Args:
        x (`numpy.ndarray`_):
        y (`numpy.ndarray`_):
        z (`numpy.ndarray`_): value of data at each (x,y) coordinate
        deg (tuple): degree of polynomial fit in the form [nx,ny]
        w (`numpy.ndarray`_):
        function:
        minx:
        maxx:
        miny:
        maxy:

    Returns:
        `numpy.ndarray`_:

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
        msgs.error("Not ready for this type of {:s}".format(function))
    # Weights
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            msgs.bug("fitting.polyfit2d - Expected 1D vector for weights")
        if len(x) != len(w) or len(y) != len(w) or len(x) != len(y):
            msgs.bug("fitting.polyfit2d - Expected x, y and weights to have same length")
        z = z * w
        vander = vander * w[:,np.newaxis]
    # Reshape
    vander = vander.reshape((-1,vander.shape[-1]))
    z = z.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, z, rcond=None)[0]
    return c.reshape(deg+1), minx, maxx, miny, maxy



def twoD_Gaussian(tup, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """ A 2D Gaussian to be used to fit the cross-correlation

    Args:
        tup (tuple):
            A two element tuple containing the (x,y) coordinates where the 2D Gaussian will be evaluated
        amplitude (float):
            The amplitude of the 2D Gaussian
        xo (float):
            The centre of the Gaussian in the x direction
        yo (float:
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
        model (`numpy.ndarray`_)
    """
    (x, y) = tup
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

# Below here are codes related to b-spline fitting
def iterfit(xdata, ydata, invvar=None, inmask=None, upper=5, lower=5, x2=None,
            maxiter=10, nord=4, bkpt=None, fullbkpt=None, kwargs_bspline={}, kwargs_reject={}):
    """Iteratively fit a b-spline set to data, with rejection. This is a utility function that allows
    the bspline to use via a direct function call.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.
    ydata : :class:`numpy.ndarray`
        Dependent variable.
    invvar : :class:`numpy.ndarray`
        Inverse variance of `ydata`.  If not set, it will be calculated based
        on the standard deviation.
    upper : :class:`int` or :class:`float`
        Upper rejection threshold in units of sigma, defaults to 5 sigma.
    lower : :class:`int` or :class:`float`
        Lower rejection threshold in units of sigma, defaults to 5 sigma.
    x2 : :class:`numpy.ndarray`, optional
        Orthogonal dependent variable for 2d fits.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations, default 10.  Set this to
        zero to disable rejection.

    Returns
    -------
    :func:`tuple`
        A tuple containing the fitted bspline object and an output mask.
    """
    # from .math import djs_reject
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if invvar is not None:
        if invvar.size != nx:
            raise ValueError('Dimensions of xdata and invvar do not agree.')
    else:
        #
        # This correction to the variance makes it the same
        # as IDL's variance()
        #
        var = ydata.var() * (float(nx) / float(nx - 1))
        if var == 0:
            var = 1.0
        invvar = np.ones(ydata.shape, dtype=ydata.dtype) / var

    if inmask is None:
        inmask = invvar > 0.0

    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & inmask & (invvar > 0.0))[xsort]
    if 'oldset' in kwargs_bspline:
        sset = kwargs_bspline['oldset']
        sset.mask[:] = True
        sset.coeff[:] = 0.
    else:
        if not maskwork.any():
            raise ValueError('No valid data points.')
            # return (None,None)
        # JFH comment this out for now
        #        if 'fullbkpt' in kwargs:
        #            fullbkpt = kwargs['fullbkpt']
        else:
            sset = bspline.bspline(xdata[xsort[maskwork]], nord=nord, bkpt=bkpt, fullbkpt=fullbkpt, **kwargs_bspline)
            if maskwork.sum() < sset.nord:
                print('Number of good data points fewer than nord.')
                return (sset, outmask)
            if x2 is not None:
                if 'xmin' in kwargs_bspline:
                    xmin = kwargs_bspline['xmin']
                else:
                    xmin = x2.min()
                if 'xmax' in kwargs_bspline:
                    xmax = kwargs_bspline['xmax']
                else:
                    xmax = x2.max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs_bspline:
                    sset.funcname = kwargs_bspline['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if x2 is not None:
        x2work = x2[xsort]
    else:
        x2work = None
    iiter = 0
    error = -1
    qdone = False
    while (error != 0 or qdone is False) and iiter <= maxiter:
        goodbk = sset.mask.nonzero()[0]
        if maskwork.sum() <= 1 or not sset.mask.any():
            sset.coeff[:] = 0.
            iiter = maxiter + 1  # End iterations
        else:
            if 'requiren' in kwargs_bspline:
                i = 0
                while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx - 1:
                    i += 1
                ct = 0
                for ileft in range(sset.nord, sset.mask.sum() - sset.nord + 1):
                    while (xwork[i] >= sset.breakpoints[goodbk[ileft]] and
                           xwork[i] < sset.breakpoints[goodbk[ileft + 1]] and
                           i < nx - 1):
                        ct += invwork[i] * maskwork[i] > 0
                        i += 1
                    if ct >= kwargs_bspline['requiren']:
                        ct = 0
                    else:
                        sset.mask[goodbk[ileft]] = False
            error, yfit = sset.fit(xwork, ywork, invwork * maskwork,
                                   x2=x2work)
        iiter += 1
        inmask_rej = maskwork
        if error == -2:

            return (sset, outmask)
        elif error == 0:
            # ToDO JFH by setting inmask to be tempin which is maskwork, we are basically implicitly enforcing sticky rejection
            # here. See djs_reject.py. I'm leaving this as is for consistency with the IDL version, but this may require
            # further consideration. I think requiring stick to be set is the more transparent behavior.
            maskwork, qdone = pydl.djs_reject(ywork, yfit, invvar=invwork, inmask=inmask_rej, outmask=maskwork,
                                              upper=upper, lower=lower, **kwargs_reject)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset, outmask)


def bspline_profile(xdata, ydata, invvar, profile_basis, ingpm=None, upper=5, lower=5, maxiter=25,
                    nord=4, bkpt=None, fullbkpt=None, relative=None, kwargs_bspline={},
                    kwargs_reject={}, quiet=False):
    """
    Fit a B-spline in the least squares sense with rejection to the
    provided data and model profiles.

    .. todo::
        Fully describe procedure.

    Parameters
    ----------
    xdata : `numpy.ndarray`_
        Independent variable.
    ydata : `numpy.ndarray`_
        Dependent variable.
    invvar : `numpy.ndarray`_
        Inverse variance of `ydata`.
    profile_basis : `numpy.ndarray`_
        Model profiles.
    ingpm : `numpy.ndarray`_, optional
        Input good-pixel mask. Values to fit in ``ydata`` should be
        True.
    upper : :obj:`int`, :obj:`float`, optional
        Upper rejection threshold in units of sigma, defaults to 5
        sigma.
    lower : :obj:`int`, :obj:`float`, optional
        Lower rejection threshold in units of sigma, defaults to 5
        sigma.
    maxiter : :obj:`int`, optional
        Maximum number of rejection iterations, default 10. Set this
        to zero to disable rejection.
    nord : :obj:`int`, optional
        Order of B-spline fit
    bkpt : `numpy.ndarray`_, optional
        Array of breakpoints to be used for the b-spline
    fullbkpt : `numpy.ndarray`_, optional
        Full array of breakpoints to be used for the b-spline,
        without letting the b-spline class append on any extra bkpts
    relative : `numpy.ndarray`_, optional
        Array of integer indices to be used for computing the reduced
        chi^2 of the fits, which then is used as a scale factor for
        the upper,lower rejection thresholds
    kwargs_bspline : :obj:`dict`, optional
        Keyword arguments used to instantiate
        :class:`pypeit.bspline.bspline`
    kwargs_reject : :obj:`dict`, optional
        Keyword arguments passed to :func:`pypeit.core.pydl.djs_reject`
    quiet : :obj:`bool`, optional
        Suppress output to the screen

    Returns
    -------
    sset : :class:`pypeit.bspline.bspline`
        Result of the fit.
    gpm : `numpy.ndarray`_
        Output good-pixel mask which the same size as ``xdata``. The
        values in this array for the corresponding data are not used in
        the fit, either because the input data was masked or the data
        were rejected during the fit, if they are False. Data
        rejected during the fit (if rejection is performed) are::

            rejected = ingpm & np.logical_not(gpm)

    yfit : `numpy.ndarray`_
        The best-fitting model; shape is the same as ``xdata``.
    reduced_chi : :obj:`float`
        Reduced chi-square of the best-fitting model.
    exit_status : :obj:`int`
        Indication of the success/failure of the fit.  Values are:

            - 0 = fit exited cleanly
            - 1 = maximum iterations were reached
            - 2 = all points were masked
            - 3 = all break points were dropped
            - 4 = Number of good data points fewer than nord

    """
    # Checks
    nx = xdata.size
    if ydata.size != nx:
        msgs.error('Dimensions of xdata and ydata do not agree.')

    # TODO: invvar and profile_basis should be optional

    # ToDO at the moment invvar is a required variable input
    #    if invvar is not None:
    #        if invvar.size != nx:
    #            raise ValueError('Dimensions of xdata and invvar do not agree.')
    #        else:
    #            #
    #            # This correction to the variance makes it the same
    #            # as IDL's variance()
    #            #
    #            var = ydata.var()*(float(nx)/float(nx-1))
    #            if var == 0:
    #                var = 1.0
    #            invvar = np.ones(ydata.shape, dtype=ydata.dtype)/var

    npoly = int(profile_basis.size / nx)
    if profile_basis.size != nx * npoly:
        msgs.error('Profile basis is not a multiple of the number of data points.')

    # Init
    yfit = np.zeros(ydata.shape)
    reduced_chi = 0.

    # TODO: Instanting these place-holder arrays can be expensive.  Can we avoid doing this?
    outmask = True if invvar.size == 1 else np.ones(invvar.shape, dtype=bool)

    if ingpm is None:
        ingpm = invvar > 0

    if not quiet:
        termwidth = 80 - 13
        msgs.info('B-spline fit:')
        msgs.info('    npoly = {0} profile basis functions'.format(npoly))
        msgs.info('    ngood = {0}/{1} measurements'.format(np.sum(ingpm), ingpm.size))
        msgs.info(' {0:>4}  {1:>8}  {2:>7}  {3:>6} '.format(
            'Iter', 'Chi^2', 'N Rej', 'R. Fac').center(termwidth, '*'))
        hlinestr = ' {0}  {1}  {2}  {3} '.format('-' * 4, '-' * 8, '-' * 7, '-' * 6)
        nullval = '  {0:>8}  {1:>7}  {2:>6} '.format('-' * 2, '-' * 2, '-' * 2)
        msgs.info(hlinestr.center(termwidth))

    maskwork = outmask & ingpm & (invvar > 0)
    if not maskwork.any():
        msgs.error('No valid data points in bspline_profile!.')

    # Init bspline class
    sset = bspline.bspline(xdata[maskwork], nord=nord, npoly=npoly, bkpt=bkpt, fullbkpt=fullbkpt,
                   funcname='Bspline longslit special', **kwargs_bspline)
    if maskwork.sum() < sset.nord:
        if not quiet:
            msgs.warn('Number of good data points fewer than nord.')
        # TODO: Why isn't maskwork returned?
        return sset, outmask, yfit, reduced_chi, 4

    # This was checked in detail against IDL for identical inputs
    # KBW: Tried a few things and this was about as fast as you can get.
    outer = np.outer(np.ones(nord, dtype=float), profile_basis.flatten('F')).T
    action_multiple = outer.reshape((nx, npoly * nord), order='F')
    # --------------------
    # Iterate spline fit
    iiter = 0
    error = -1  # Indicates that the fit should be done
    qdone = False  # True if rejection iterations are done
    exit_status = 0
    relative_factor = 1.0
    nrel = 0 if relative is None else len(relative)
    # TODO: Why do we need both maskwork and tempin?
    tempin = np.copy(ingpm)
    while (error != 0 or qdone is False) and iiter <= maxiter and exit_status == 0:
        ngood = maskwork.sum()
        goodbk = sset.mask.nonzero()[0]
        if ngood <= 1 or not sset.mask.any():
            sset.coeff[:] = 0.
            exit_status = 2  # This will end iterations
        else:
            # Do the fit. Return values from workit for error are as follows:
            #    0 if fit is good
            #   -1 if some breakpoints are masked, so try the fit again
            #   -2 if everything is screwed

            # we'll do the fit right here..............
            if error != 0:
                bf1, laction, uaction = sset.action(xdata)
                if np.any(bf1 == -2) or bf1.size != nx * nord:
                    msgs.error("BSPLINE_ACTION failed!")
                action = np.copy(action_multiple)
                for ipoly in range(npoly):
                    action[:, np.arange(nord) * npoly + ipoly] *= bf1
                del bf1  # Clear the memory

            if np.any(np.logical_not(np.isfinite(action))):
                msgs.error('Infinities in action matrix.  B-spline fit faults.')

            error, yfit = sset.workit(xdata, ydata, invvar * maskwork, action, laction, uaction)

        iiter += 1

        if error == -2:
            if not quiet:
                msgs.warn('All break points lost!!  Bspline fit failed.')
            exit_status = 3
            return sset, np.zeros(xdata.shape, dtype=bool), np.zeros(xdata.shape), reduced_chi, \
                   exit_status

        if error != 0:
            if not quiet:
                msgs.info((' {0:4d}'.format(iiter) + nullval).center(termwidth))
            continue

        # Iterate the fit -- next rejection iteration
        chi_array = (ydata - yfit) * np.sqrt(invvar * maskwork)
        reduced_chi = np.sum(np.square(chi_array)) / (ngood - npoly * (len(goodbk) + nord) - 1)

        relative_factor = 1.0
        if relative is not None:
            this_chi2 = reduced_chi if nrel == 1 \
                else np.sum(np.square(chi_array[relative])) \
                     / (nrel - (len(goodbk) + nord) - 1)
            relative_factor = max(np.sqrt(this_chi2), 1.0)

        # Rejection

        # TODO: JFH by setting ingpm to be tempin which is maskwork, we
        #  are basically implicitly enforcing sticky rejection here. See
        #  djs_reject.py. I'm leaving this as is for consistency with
        #  the IDL version, but this may require further consideration.
        #  I think requiring sticky to be set is the more transparent
        #  behavior.
        maskwork, qdone = pydl.djs_reject(ydata, yfit, invvar=invvar, inmask=tempin,
                                          outmask=maskwork, upper=upper * relative_factor,
                                          lower=lower * relative_factor, **kwargs_reject)
        tempin = np.copy(maskwork)
        if not quiet:
            msgs.info(' {0:4d}  {1:8.3f}  {2:7d}  {3:6.2f} '.format(iiter,
                                                                    reduced_chi, np.sum(maskwork == 0),
                                                                    relative_factor).center(termwidth))

    if iiter == (maxiter + 1):
        exit_status = 1

    # Exit status:
    #    0 = fit exited cleanly
    #    1 = maximum iterations were reached
    #    2 = all points were masked
    #    3 = all break points were dropped
    #    4 = Number of good data points fewer than nord

    if not quiet:
        msgs.info(' {0:>4}  {1:8.3f}  {2:7d}  {3:6.2f} '.format('DONE',
                                                                reduced_chi, np.sum(maskwork == 0),
                                                                relative_factor).center(termwidth))
        msgs.info('*' * termwidth)

    # Finish
    # TODO: Why not return maskwork directly
    outmask = np.copy(maskwork)
    # Return
    return sset, outmask, yfit, reduced_chi, exit_status


def bspline_qa(xdata, ydata, sset, gpm, yfit, xlabel=None, ylabel=None, title=None, show=True):
    """
    Construct a QA plot of the bspline fit.

    Args:
        xdata (`numpy.ndarray`_):
            Array with the independent variable. Regardless of shape,
            data is treated as one-dimensional.
        ydata (`numpy.ndarray`_):
            Array with the dependent variable. Regardless of shape,
            data is treated as one-dimensional.
        sset (:class:`pypeit.bspline.bspline`):
            Object with the results of the fit. (First object
            returned by :func:`bspline_profile`).
        gpm (`numpy.ndarray`_):
            Boolean array with the same size as ``xdata``.
            Measurements rejected during the fit have ``gpm=False``.
            (Second object returned by :func:`bspline_profile`).
        yfit (`numpy.ndarray`_):
            Best-fitting model sampled at ``xdata``. (Third object
            returned by :func:`bspline_profile`).
        xlabel (:obj:`str`, optional):
            Label for the ordinate.  If None, none given.
        ylabel (:obj:`str`, optional):
            Label for the abcissa.  If None, none given.
        title (:obj:`str`, optional):
            Label for the plot.  If None, none given.
        show (:obj:`bool`, optional):
            Plot the result. If False, the axis instance is returned.
            This is done before any labels or legends are added to
            the plot.

    Returns:
        `matplotlib.axes.Axes`_: Axes instance with the data, model,
        and breakpoints.  Only returned if ``show`` is False.
    """
    goodbk = sset.mask
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

