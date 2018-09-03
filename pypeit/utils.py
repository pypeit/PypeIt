from __future__ import absolute_import, division, print_function, unicode_literals

import os
import warnings
import itertools

import numpy as np

from scipy.optimize import curve_fit
from scipy import interpolate

from astropy import units
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel

#from pydl.pydlutils import math
#from pydl.pydlutils import bspline

from pypeit.core import pydl

from pypeit import msgs
from pypeit import debugger

def quicksave(data,fname):
    """
    Save a fits file (quickly) -- overwrite is forced, and no quality control checks
    """
    hdu = fits.PrimaryHDU(data)
    hdulist = fits.HDUList([hdu])
    if os.path.exists(fname):
        os.remove(fname)
    hdulist.writeto(fname)
    return

def bspline_inner_knots(all_knots):
    """Trim to the inner knots.  Used in bspline_magfit
    Might be useful elsewhere
    """
    diff = all_knots - np.roll(all_knots,1)
    pos = np.where(diff>0.)[0]
    i0=pos[0]
    i1=pos[-1]
    return all_knots[i0:i1]

def fast_running_median(seq, window_size):
    """
      Compute the median of sequence of numbers with a running window. The boundary conditions are identical to the
      scipy 'reflect' boundary codition:

         'reflect' (`d c b a | a b c d | d c b a`)
         The input is extended by reflecting about the edge of the last pixel.

      This code has been confirmed to produce identical results to scipy.ndimage.filters.median_filter with the reflect
      boundary condition, but is ~ 100 times faster.

      Parameters
      ----------
      seq : list or 1-d numpy array of numbers.

      window_size = size of running window.

      Returns
      -------
      ndarray of median filtered values

      Code contributed by Peter Otten, made to be consistent with scipy.ndimage.filters.median_filter by Joe Hennawi.

      See discussion at:
      http://groups.google.com/group/comp.lang.python/browse_thread/thread/d0e011c87174c2d0
      """



    # pad the array for the reflection
    seq_pad = np.concatenate((seq[0:window_size][::-1],seq,seq[-1:(-1-window_size):-1]))

    window_size= int(window_size)
    seq_pad = iter(seq_pad)
    d = deque()
    s = []
    result = []
    for item in islice(seq_pad, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
    m = window_size // 2
    for item in seq_pad:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])

    # This takes care of the offset produced by the original code deducec by trial and error comparison with
    # scipy.ndimage.filters.medfilt

    result = np.roll(result, -window_size//2 + 1)
    return result[window_size:-window_size]


# TODO JFH: This is the old bspline_fit which shoul be deprecated. I think some codes still use it though. We should transtion to pydl everywhere
def bspline_fit(x,y,order=3,knots=None,everyn=20,xmin=None,xmax=None,w=None,bkspace=None):
    ''' bspline fit to x,y
    Should probably only be called from func_fit

    Parameters:
    ---------
    x: ndarray
    y: ndarray
    func: str
      Name of the fitting function:  polynomial, legendre, chebyshev, bspline
    deg: int 
      deg of the spline.  Default=3 (cubic)
    xmin: float, optional
      Minimum value in the array  [both must be set to normalize]
    xmax: float, optional
      Maximum value in the array  [both must be set to normalize]
    w: ndarray, optional
      weights to be used in the fitting (weights = 1/sigma)
    knots: ndarray, optional
      Internal knots only.  External ones are added by scipy
    everyn: int 
      Knot everyn good pixels, if used
    bkspace: float 
      Spacing of breakpoints in units of x

    Returns:
    ---------
    tck : tuple
      describes the bspline
    '''
    task = 0  # Default of splrep
    if w is None:
        ngd = x.size
        gd = np.arange(ngd)
        weights = None
    else:
        gd = np.where(w > 0.)[0]
        weights = w[gd]
        ngd = len(gd)
    # Make the knots
    if knots is None:
        if bkspace is not None: 
            xrnge = (np.max(x[gd]) - np.min(x[gd]))
            startx = np.min(x[gd])
            nbkpts = max(int(xrnge/bkspace) + 1,2)
            tempbkspace = xrnge/(nbkpts-1)
            knots = np.arange(1, nbkpts-1)*tempbkspace + startx
            # Remove cases where two knots have no data between them
            keep_knots = np.array([True]*len(knots))
            for ii in range(1,len(knots)): # Ugly for loop..
                if not np.any((x[gd] > knots[ii-1]) & (x[gd] < knots[ii])):
                    keep_knots[ii] = False
            knots = knots[keep_knots]
        elif everyn is not None:
            # A knot every good N pixels
            idx_knots = np.arange(everyn//2, ngd-everyn//2, everyn)
            knots = x[gd[idx_knots]]
        else:
            msgs.error("No method specified to generate knots")
    else:
        task = -1
    # Generate spline
    try:
        tck = interpolate.splrep(x[gd], y[gd], w=weights, k=order, xb=xmin, xe=xmax, t=knots, task=task)
    except ValueError:
        # Knot problem (usually)
        msgs.warn("Problem in the bspline knot")
        raise ValueError("Crashing out of bspline fitting")
    return tck

#ToDo consider adding an inmask here as a keyword argument. Also I would prefer to remove the kwargs_bspline and
# and make them explicit
def bspline_profile(xdata, ydata, invvar, profile_basis, inmask = None, upper=5, lower=5,
                    maxiter=25, nord = 4, bkpt=None, fullbkpt=None,
                    relative=None, kwargs_bspline={}, kwargs_reject={}):
    """
    Create a B-spline in the least squares sense with rejection, using a model profile

     Parameters
     ----------
     xdata : :class:`numpy.ndarray`
         Independent variable.
     ydata : :class:`numpy.ndarray`
         Dependent variable.
     invvar : :class:`numpy.ndarray`
         Inverse variance of `ydata`.
     profile_basis : :class:`numpy.ndarray`
         model profiles
     upper : :class:`int` or :class:`float`, optional
         Upper rejection threshold in units of sigma, defaults to 5 sigma.
     lower : :class:`int` or :class:`float`, optional
         Lower rejection threshold in units of sigma, defaults to 5 sigma.
     maxiter : :class:`int`, optional
         Maximum number of rejection iterations, default 10.  Set this to
         zero to disable rejection.
     nord : :class:`int`, optional
         Order of B-spline fit
     bkpt : :class:`numpy.ndarray`
         Array of breakpoints to be used for the b-spline
     fullbkpt : :class:`numpy.ndarray`
         Full array of breakpoints to be used for the b-spline, without letting the b-spline class append on any extra bkpts
     relative : class:`numpy.ndarray`
        Array of integer indices to be used for computing the reduced chi^2 of the fits, which then is used as a scale factor for
         the upper,lower rejection thresholds
     kwargs_bspline : dict
       Passed to bspline
     kwargs_reject : dict
       Passed to djs_reject

     Returns
     -------
     :func:`tuple`
         A tuple containing the (sset, outmask, yfit, reduced_chi), where

            sset: object
               bspline object
            outmask: : :class:`numpy.ndarray`
               output mask which the same size as xdata, such that rejected points have outmask set to False
            yfit  : :class:`numpy.ndarray`
               result of the bspline fit (same size as xdata)
            reduced_chi: float
               value of the reduced chi^2
     """
    # Checks
    nx = xdata.size
    if ydata.size != nx:
        msgs.error('Dimensions of xdata and ydata do not agree.')

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

    npoly = int(profile_basis.size/nx)
    if profile_basis.size != nx*npoly:
        msgs.error('Profile basis is not a multiple of the number of data points.')

    msgs.info("Fitting npoly =" + "{:3d}".format(npoly) + " profile basis functions, nx=" + "{:3d}".format(nx) + " pixels")
    msgs.info("****************************  Iter  Chi^2  # rejected  Rel. fact   ****************************")
    msgs.info("                              ----  -----  ----------  --------- ")

    # Init
    yfit = np.zeros(ydata.shape)
    reduced_chi = 0.

    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')

    if inmask is None:
        inmask = (invvar > 0)
    maskwork = outmask & inmask & (invvar > 0)
    if not maskwork.any():
        msgs.error('No valid data points in bspline_profile!.')
    else:
        #from IPython import embed
        #embed()
        # Init bspline class
        sset = pydl.bspline(xdata[maskwork], nord=nord, npoly=npoly, bkpt=bkpt, fullbkpt=fullbkpt,
                       funcname='Bspline longslit special', **kwargs_bspline)
        if maskwork.sum() < sset.nord:
            msgs.warn('Number of good data points fewer than nord.')
            return (sset, outmask, yfit, reduced_chi)

    # This was checked in detail against IDL for identical inputs
    outer = (np.outer(np.ones(nord, dtype=float), profile_basis.flatten('F'))).T
    action_multiple = outer.reshape((nx, npoly * nord), order='F')
    #--------------------
    # Iterate spline fit
    iiter = 0
    error = -1
    qdone = False

    relative_factor = 1.0
    tempin = np.copy(inmask)
    while (error != 0 or qdone is False) and iiter <= maxiter:
        ngood = maskwork.sum()
        goodbk = sset.mask.nonzero()[0]
        if ngood <= 1 or not sset.mask.any():
            sset.coeff = 0
            iiter = maxiter + 1 # End iterations
        else:
            # Do the fit. Return values from workit for error are as follows:
            #    0 if fit is good
            #   -1 if some breakpoints are masked, so try the fit again
            #   -2 if everything is screwed

            # we'll do the fit right here..............
            if error != 0:
                bf1, laction, uaction = sset.action(xdata)
                if(np.any(bf1 == -2) or (bf1.size !=nx*nord)):
                    msgs.error("BSPLINE_ACTION failed!")
                action = np.copy(action_multiple)
                for ipoly in range(npoly):
                    action[:, np.arange(nord)*npoly + ipoly] *= bf1
                del bf1 # Clear the memory
            if np.sum(np.isfinite(action) is False) > 0:
                msgs.error("Infinities in action matrix, wavelengths may be very messed up!!!")
            error, yfit = sset.workit(xdata, ydata, invvar*maskwork,action, laction, uaction)
        iiter += 1
        if error == -2:
            msgs.warn(" All break points have been dropped!!")
            return (sset, outmask, yfit, reduced_chi)
        elif error == 0:
            # Iterate the fit -- next rejection iteration
            chi_array = (ydata - yfit)*np.sqrt(invvar * maskwork)
            reduced_chi = np.sum(chi_array**2)/(ngood - npoly*(len(goodbk) + nord)-1)
            relative_factor = 1.0
            # JFH -- What is
            if relative is not None:
                nrel = len(relative)
                if nrel == 1:
                    relative_factor = np.sqrt(reduced_chi)
                else:
                    this_chi2 = np.sum(chi_array[relative]**2)/(nrel - (len(goodbk) + nord) - 1)
                    relative_factor = np.sqrt(this_chi2)
                relative_factor = max(relative_factor,1.0)
            # Rejection
            maskwork, qdone = pydl.djs_reject(ydata, yfit, invvar=invvar,
                                         inmask=tempin, outmask=maskwork,
                                         upper=upper*relative_factor,
                                         lower=lower*relative_factor, **kwargs_reject)
            tempin = np.copy(maskwork)
            msgs.info("                             {:4d}".format(iiter) + "{:8.3f}".format(reduced_chi) +
                      "  {:7d}".format((maskwork == 0).sum()) + "      {:6.2f}".format(relative_factor))

        else:
            msgs.info("                             {:4d}".format(iiter) + "    ---    ---    ---    ---")


    msgs.info("***********************************************************************************************")
    msgs.info(
        "Final fit after " + "{:2d}".format(iiter) + " iterations: reduced_chi = " + "{:8.3f}".format(reduced_chi) +
        ", rejected = " + "{:7d}".format((maskwork == 0).sum()) + ", relative_factor = {:6.2f}".format(relative_factor))
    # Finish
    outmask = np.copy(maskwork)
    # Return
    return sset, outmask, yfit, reduced_chi



def calc_ivar(varframe):
    """ Calculate the inverse variance based on the input array
    """
    ivar = (varframe > 0.) / (np.abs(varframe) + (varframe == 0))
    return ivar


def calc_offset(raA, decA, raB, decB, distance=False):
    """
    Calculate the offset in arcseconds between two sky coordinates
    All coordinates must be in decimal degrees.
    """
    delRA  = 3600.0*(raA-raB)*np.cos(decA*np.pi/180.0)
    delDEC = 3600.0*(decA-decB)
    if distance==True:
        return np.sqrt(delRA**2 + delDEC**2)
    else:
        return delRA, delDEC

def dummy_fitsdict(nfile=10, spectrograph='shane_kast_blue', directory='./'):
    pass



def func_der(coeffs, func, nderive=1):
    if func == "polynomial":
        return np.polynomial.polynomial.polyder(coeffs, m=nderive)
    elif func == "legendre":
        return np.polynomial.legendre.legder(coeffs, m=nderive)
    elif func == "chebyshev":
        return np.polynomial.chebyshev.chebder(coeffs, m=nderive)
    else:
        msgs.error("Functional derivative '{0:s}' is not implemented yet"+msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev'")


def func_fit(x, y, func, deg, minv=None, maxv=None, w=None, guesses=None,
             bspline_par=None):
    """ General routine to fit a function to a given set of x,y points

    Parameters
    ----------
    x : ndarray
    y : ndarray
    func : str
      polynomial, legendre, chebyshev, bspline, gaussian
    deg : int
      degree of the fit
    minv : float, optional
    maxv
    w
    guesses : tuple
    bspline_par : dict
      Passed to bspline_fit()

    Returns
    -------
    coeff : ndarray or tuple
      ndarray for standard function fits
      tuple for bspline

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyfit(x, y, deg, w=w)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legfit(xv, y, deg, w=w)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebfit(xv, y, deg, w=w)
    elif func == "bspline":
        if bspline_par is None:
            bspline_par = {}
        # TODO -- Deal with this kwargs-like kludge
        return bspline_fit(x, y, order=deg, w=w, **bspline_par)
    elif func in ["gaussian"]:
        # Guesses
        if guesses is None:
            mx, cent, sigma = guess_gauss(x, y)
        else:
            if deg == 2:
                mx, sigma = guesses
            elif deg == 3:
                mx, cent, sigma = guesses
        # Error
        if w is not None:
            sig_y = 1./w
        else:
            sig_y = None
        if deg == 2:  # 2 parameter fit
            popt, pcov = curve_fit(gauss_2deg, x, y, p0=[mx, sigma], sigma=sig_y)
        elif deg == 3:  # Standard 3 parameters
            popt, pcov = curve_fit(gauss_3deg, x, y, p0=[mx, cent, sigma],
                                   sigma=sig_y)
        else:
            msgs.error("Not prepared for deg={:d} for Gaussian fit".format(deg))
        # Return
        return popt
    elif func in ["moffat"]:
        # Guesses
        if guesses is None:
            mx, cent, sigma = guess_gauss(x, y)
            p0 = mx
            p2 = 3. # Standard guess
            p1 = (2.355*sigma)/(2*np.sqrt(2**(1./p2)-1))
        else:
            p0,p1,p2 = guesses
        # Error
        if w is not None:
            sig_y = 1./w
        else:
            sig_y = None
        if deg == 3:  # Standard 3 parameters
            popt, pcov = curve_fit(moffat, x, y, p0=[p0,p1,p2], sigma=sig_y)
        else:
            msgs.error("Not prepared for deg={:d} for Moffat fit".format(deg))
        # Return
        return popt
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev','bspline'")


def func_val(c, x, func, minv=None, maxv=None):
    """ Generic routine to return an evaluated function
    Functional forms include:
      polynomial, legendre, chebyshev, bspline, gauss

    Parameters
    ----------
    c : ndarray
      coefficients
    x
    func
    minv
    maxv

    Returns
    -------
    values : ndarray

    """
    if func == "polynomial":
        return np.polynomial.polynomial.polyval(x, c)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legval(xv, c)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebval(xv, c)
    elif func == "bspline":
        return interpolate.splev(x, c, ext=1)
    elif func == "gaussian":
        if len(c) == 2:
            return gauss_2deg(x, c[0], c[1])
        elif len(c) == 3:
            return gauss_3deg(x, c[0], c[1], c[2])
        else:
            msgs.error("Not ready for this type of gaussian")
    elif func == "moffat":
        if len(c) == 3:
            return moffat(x, c[0], c[1], c[2])
        else:
            msgs.error("Not ready for this type of Moffat")
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev', 'bspline'")


def calc_fit_rms(xfit, yfit, fit, func, minv=None, maxv=None):
    """ Simple RMS calculation

    Parameters
    ----------
    xfit : ndarray
    yfit : ndarray
    fit : coefficients
    func : str
    minv : float, optional
    maxv : float, optional

    Returns
    -------
    rms : float

    """
    values = func_val(fit, xfit, func, minv=minv, maxv=maxv)
    rms = np.std(yfit-values)
    # Return
    return rms


def func_vander(x, func, deg, minv=None, maxv=None):
    if func == "polynomial":
        return np.polynomial.polynomial.polyvander(x, deg)
    elif func == "legendre":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legvander(xv, deg)
    elif func == "chebyshev":
        if minv is None or maxv is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minv, maxv
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebvander(xv, deg)
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev'")


def mask_polyfit(xarray,yarray,order,maxone=True,sigma=3.0):
    mask = np.zeros(xarray.size,dtype=np.int)
    mskcnt=0
    while True:
        w = np.where(mask==0)
        xfit = xarray[w]
        yfit = yarray[w]
        ct = np.polyfit(xfit,yfit,order)
        yrng = np.polyval(ct,xarray)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if xarray.size-np.sum(mask) <= order+2:
            msgs.warn("More parameters than data points - fit might be undesirable")
            break # More data was masked than allowed by order
        if maxone: # Only remove the most deviant point
            tst = np.abs(yarray[w]-yrng[w])
            m = np.argmax(tst)
            if tst[m] > sigma*sigmed:
                mask[w[0][m]] = 1
        else:
            w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break # No new values have been included in the mask
        mskcnt = np.sum(mask)
    return mask


def robust_meanstd(array):
    """
    Determine a robust measure of the mean and dispersion of array
    :param array: an array of values
    :return: median of the array and a robust estimate of the standand deviation (assuming a symmetric distribution)
    """
    med = np.median(array)
    mad = np.median(np.abs(array-med))
    return med, 1.4826*mad


def perturb(covar, bparams, nsim=1000):
    cvsize = covar.shape[0]
    # Generate a new set of starting parameters from the covariance matrix
    X_covar_fit = np.matrix(np.random.standard_normal((cvsize, nsim)))
    C_covar_fit = np.matrix(covar)
    U_covar_fit = np.linalg.cholesky(C_covar_fit)
    Y_covar_fit = U_covar_fit * X_covar_fit
    newpar = bparams.reshape(cvsize, 1).repeat(nsim, axis=1)
    # Find the new best-fitting model parameters
    newpar += Y_covar_fit
    return newpar


def polyfitter2d(data, mask=None, order=2):
    x, y = np.meshgrid(np.linspace(0.0, 1.0, data.shape[1]), np.linspace(0.0, 1.0, data.shape[0]))
    if isinstance(mask, (float, int)):
        # mask is the value that should be masked in data
        w = np.where(data != mask)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    elif mask is None or mask.size == 0:
        # There are no masks
        xf = x.flatten()
        yf = y.flatten()
        m = polyfit2d(xf, yf, data.T.flatten(), order)
    elif len(mask.shape) == 1:
        # mask is applied along one axis
        mskar = np.ones((data.shape[0], data.shape[1]))
        mskar[mask, :] = 0
        w = np.where(mskar == 1)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    elif mask.shape[0] == data.shape[0] and mask.shape[1] == data.shape[1]:
        # mask is an array that indicates the masked data
        w = np.where(mask == 0)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    # Return the best model
    return m, polyval2d(x, y, m).T
    # print m
    # model = polyval2d(x,y,m)
    # # Plot it
    # zmin, zmax = arplot.zscale(data)
    # plt.subplot(131)
    # implt = plt.imshow(data,aspect='auto',interpolation='none')
    # plt.colorbar()
    # implt.set_clim(zmin, zmax)
    # plt.subplot(132)
    # implt = plt.imshow(model,aspect='auto',interpolation='none')
    # plt.colorbar()
    # implt.set_clim(zmin, zmax)
    # plt.subplot(133)
    # zmin, zmax = arplot.zscale(data-model)
    # implt = plt.imshow(data-model,aspect='auto',interpolation='none')
    # plt.colorbar()
    # implt.set_clim(zmin, zmax)
    # plt.show()
    # return m, polyval2d(x,y,m).T


def polyfit2d(x, y, z, order=3):
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, null, null, null = np.linalg.lstsq(G, z)
    return m


def polyval2d(x, y, m):
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def moffat(x,p0,p1,p2):
    """  Moffat profile
    This 3 parameter formulation assumes the trace is known
    Parameters
    ----------
    x
    p0 : float
      Amplitude
    p1 : float
      Width scaling
    p2 : float

    Returns
    -------
    Evaluated Moffat
    """
    return p0 / (1+(x/p1)**2)**p2

def gauss_2deg(x,ampl,sigm):
    """  Simple 2 parameter Gaussian (amplitude, sigma)
    Parameters
    ----------
    x
    ampl
    sigm

    Returns
    -------
    Evaluated Gausssian
    """
    return ampl*np.exp(-1.*x**2/2./sigm**2)


def gauss_3deg(x,ampl,cent,sigm):
    """  Simple 3 parameter Gaussian
    Parameters
    ----------
    x
    ampl
    cent
    sigm

    Returns
    -------
    Evaluated Gausssian
    """
    return ampl*np.exp(-1.*(cent-x)**2/2/sigm**2)

def guess_gauss(x,y):
    """ Guesses Gaussian parameters with basic stats

    Parameters
    ----------
    x
    y

    Returns
    -------

    """
    cent = np.sum(y*x)/np.sum(y)
    sigma = np.sqrt(np.abs(np.sum((x-cent)**2*y)/np.sum(y))) # From scipy doc
    # Calculate mx from pixels within +/- sigma/2
    cen_pix = np.where(np.abs(x-cent)<sigma/2)
    mx = np.median(y[cen_pix])
    # Return
    return mx, cent, sigma


def poly_to_gauss(coeffs):
    try:
        sigm = np.sqrt(-0.5/coeffs[2])
        cent = -0.5*coeffs[1]/coeffs[2]
        ampl = np.exp( coeffs[0] + 0.5*(cent/sigm)**2 )
    except:
        return [0.0, 0.0, 0.0], True
    return [ampl, cent, sigm], False


def polyfit2d_general(x, y, z, deg, w=None, function='polynomial',
                      minx=None, maxx=None, miny=None, maxy=None):
    """
    :param x: array of x values
    :param y: array of y values
    :param z: value of data at each (x,y) coordinate
    :param deg: degree of polynomial fit in the form [nx,ny]
    :param w: weights
    :return: coefficients
    """
    msgs.work("Generalize to different polynomial types")
    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)
    deg = np.asarray(deg)
    # Vander
    if function == 'polynomial':
        vander = np.polynomial.polynomial.polyvander2d(x, y, deg)
    elif function == 'legendre':
        xv = scale_minmax(x, minx=minx, maxx=maxx)
        yv = scale_minmax(y, minx=miny, maxx=maxy)
        vander = np.polynomial.legendre.legvander2d(xv, yv, deg)
    else:
        msgs.error("Not read for this type of {:s}".format(function))
    # Weights
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            msgs.bug("arutils.polyfit2d - Expected 1D vector for weights")
        if len(x) != len(w) or len(y) != len(w) or len(x) != len(y):
            msgs.bug("arutils.polyfit2d - Expected x, y and weights to have same length")
        z = z * w
        vander = vander * w[:,np.newaxis]
    # Reshape
    vander = vander.reshape((-1,vander.shape[-1]))
    z = z.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, z)[0]
    return c.reshape(deg+1)

def scale_minmax(x, minx=None, maxx=None):
    if minx is None or maxx is None:
        if np.size(x) == 1:
            xmin, xmax = -1.0, 1.0
        else:
            xmin, xmax = np.min(x), np.max(x)
    else:
        xmin, xmax = minx, maxx
    xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
    return xv

def polyval2d_general(c, x, y, function="polynomial", minx=None, maxx=None, miny=None, maxy=None):
    if function == "polynomial":
        xx, yy = np.meshgrid(x, y)
        return np.polynomial.polynomial.polyval2d(xx, yy, c)
    elif function in ["legendre", "chebyshev"]:
        # Scale x-direction
        xv = scale_minmax(x, minx=minx, maxx=maxx)
        # Scale y-direction
        yv = scale_minmax(y, minx=miny, maxx=maxy)
        #if miny is None or maxy is None:
        #    if np.size(y) == 1:
        #        ymin, ymax = -1.0, 1.0
        #    else:
        #        ymin, ymax = np.min(y), np.max(y)
        #else:
        #    ymin, ymax = miny, maxy
        #yv = 2.0 * (y-ymin)/(ymax-ymin) - 1.0
        xx, yy = np.meshgrid(xv, yv)
        if function == "legendre":
            return np.polynomial.legendre.legval2d(xx, yy, c)
        elif function == "chebyshev":
            return np.polynomial.chebyshev.chebval2d(xx, yy, c)
    else:
        msgs.error("Function {0:s} has not yet been implemented".format(function))
    return None


def polyfit_integral(x, y, dx, deg, rcond=None, full=False, w=None):
    order = int(deg) + 1
    x = np.asarray(x)
    y = np.asarray(y)

    # check arguments.
    if deg < 0:
        msgs.bug("Expected deg >= 0")
        msgs.error("Input of function arutils.polyfit_integral is incorrect")
    if x.ndim != 1:
        msgs.bug("Expected 1D vector for x")
        msgs.error("Input of function arutils.polyfit_integral is incorrect")
    if x.size == 0:
        msgs.bug("Expected non-empty vector for x")
        msgs.error("Input of function arutils.polyfit_integral is incorrect")
    if y.ndim < 1 or y.ndim > 2:
        msgs.bug("Expected 1D or 2D array for y")
        msgs.error("Input of function arutils.polyfit_integral is incorrect")
    if len(x) != len(y):
        msgs.bug("Expected x and y to have same length")
        msgs.error("Input of function arutils.polyfit_integral is incorrect")

    # set up the least squares matrices in transposed form
    lhst = np.polynomial.polynomial.polyvander(x+dx/2.0, deg+1) - np.polynomial.polynomial.polyvander(x-dx/2.0, deg+1)
    div = np.arange(1., deg+2.).reshape(1, deg+1).repeat(x.size, axis=0)
    lhs = (lhst[:, 1:]/(dx.reshape(dx.size, 1).repeat(deg+1, axis=1)*div)).T
    rhs = y.T
    if w is not None:
        w = np.asarray(w)
        if w.ndim != 1:
            msgs.bug("Expected 1D vector for weights in arutils.polyfit2d")
        if len(x) != len(w):
            msgs.bug("Expected x and weights to have same length in arutils.polyfit2d")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None:
        rcond = len(x)*np.finfo(x.dtype).eps

    # Determine the norms of the design matrix columns.
    if issubclass(lhs.dtype.type, np.complexfloating):
        scl = np.sqrt((np.square(lhs.real) + np.square(lhs.imag)).sum(1))
    else:
        scl = np.sqrt(np.square(lhs).sum(1))
    scl[scl == 0] = 1

    # Solve the least squares problem.
    c, resids, rank, s = np.linalg.lstsq(lhs.T/scl, rhs.T, rcond)
    c = (c.T/scl).T

    # warn on rank reduction
    if rank != order and not full:
        msgs.warn("The fit result of the function arutils.polyfit_integral may be poorly conditioned")

    if full:
        return c, [resids, rank, s, rcond]
    else:
        return c


def poly_iterfit(x,y,ordr,maxrej=5):
    xfit = x.copy()
    yfit = y.copy()
    r = 0
    chisqn = None
    while r < maxrej:
        chisqp = chisqn
        wrng = np.arange(xfit.size)
        chisq = np.zeros(xfit.size)
        for i in range(xfit.size):
            sel = np.delete(wrng,i)
            c=np.polyfit(xfit[sel],yfit[sel],ordr)
            m=np.polyval(c,xfit[sel])
            chisq[i] = np.sum(((yfit[sel]-m)/m)**2)
        csa = np.argmin(chisq)
        chisqn = chisq[csa]
        if chisqp is not None:
            if chisqp-0.001 < chisqn:
                break
        xfit = np.delete(xfit,csa)
        yfit = np.delete(yfit,csa)
        r += 1
    msgs.info("Robust regression identified {0:d} outliers".format(r))
    return c


def rebin(frame, newshape):
    shape = frame.shape
    lenShape = len(shape)
    factor = np.asarray(shape)/np.asarray(newshape)
    evList = ['frame.reshape('] + \
             ['int(newshape[%d]),int(factor[%d]),'% (i, i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)' % (i+1) for i in range(lenShape)] + \
             ['/factor[%d]' % i for i in range(lenShape)]
    return eval(''.join(evList))


def robust_polyfit(xarray, yarray, order, weights=None, maxone=True, sigma=3.0,
                   function="polynomial", initialmask=None, forceimask=False,
                   minv=None, maxv=None, guesses=None, bspline_par=None, verbose=True):
    """
    A robust (equally weighted) polynomial fit is performed to the xarray, yarray pairs
    mask[i] = 1 are masked values

    :param xarray: independent variable values
    :param yarray: dependent variable values
    :param order: the order of the polynomial to be used in the fitting
    :param weights: weights to be used in the fitting (weights = 1/sigma)
    :param maxone: If True, only the most deviant point in a given iteration will be removed
    :param sigma: confidence interval for rejection
    :param function: which function should be used in the fitting (valid inputs: 'polynomial', 'legendre', 'chebyshev', 'bspline')
    :param initialmask: a mask can be supplied as input, these values will be masked for the first iteration. 1 = value masked
    :param forceimask: if True, the initialmask will be forced for all iterations
    :param minv: minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    :param maxv: maximum value in the array (or the right limit for a legendre/chebyshev polynomial)
    :return: mask, ct -- mask is an array of the masked values, ct is the coefficients of the robust polyfit.
    """
    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size, dtype=np.int)
        if forceimask:
            msgs.warn("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    mskcnt = np.sum(mask)
    # Iterate, and mask out new values on each iteration
    ct = guesses
    while True:
        w = np.where(mask == 0)
        xfit = xarray[w]
        yfit = yarray[w]
        if weights is not None:
            wfit = weights[w]
        else:
            wfit = None
        ct = func_fit(xfit, yfit, function, order, w=wfit,
                      guesses=ct, minv=minv, maxv=maxv, bspline_par=bspline_par)
        yrng = func_val(ct, xarray, function, minv=minv, maxv=maxv)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if xarray.size-np.sum(mask) <= order+2:
            if verbose:
                msgs.warn("More parameters than data points - fit might be undesirable")
            break  # More data was masked than allowed by order
        if maxone:  # Only remove the most deviant point
            tst = np.abs(yarray[w]-yrng[w])
            m = np.argmax(tst)
            if tst[m] > sigma*sigmed:
                mask[w[0][m]] = 1
        else:
            if forceimask:
                w = np.where((np.abs(yarray-yrng) > sigma*sigmed) | (initialmask==1))
            else:
                w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break  # No new values have been included in the mask
        mskcnt = np.sum(mask)
    # Final fit
    w = np.where(mask == 0)
    xfit = xarray[w]
    yfit = yarray[w]
    if weights is not None:
        wfit = weights[w]
    else:
        wfit = None
    ct = func_fit(xfit, yfit, function, order, w=wfit, minv=minv, maxv=maxv, bspline_par=bspline_par)
    return mask, ct


def robust_regression(x, y, ordr, outfrac, maxiter=100, function='polynomial', min=None, max=None):
    """
    Deprecated
    """
    msgs.bug("PYPIT using deprecated function")
    msgs.error("Please contact the authors")
    xsize=x.size
    infrac = 1.0-outfrac
    if infrac < 0.5: infrac = 0.5
    slct = int(xsize*infrac)
    if slct == xsize: slct = xsize-1
    if ordr+1 >= slct:
        if xsize <= 1:
            msgs.error("At least 2 points are required for a statistical fit")
        elif xsize == 2:
            msgs.warn("Only a constant can be fit to 2 points")
            msgs.info("Fitting a constant instead")
            return func_fit(x,y,function,0)
        elif  ordr+1 >= xsize:
            msgs.warn("Not enough points ({0:d}) for a {1:d}th order fit".format(xsize,ordr))
            ordr = xsize-3
            slct = ordr+2
            msgs.info("Changing order to a {0:d} order {1:s} fucntion".format(ordr,function))
        else:
            slct = ordr
    indx = np.arange(xsize)
    np.random.shuffle(indx)
    ind = indx[:slct]
    i=0
    while True:
        tc = func_fit(x[ind],y[ind],function,ordr)
        diff = np.abs(y[ind]-func_val(tc,x[ind],function))
        mad = np.median(diff)
        w=np.argsort(diff)
        inds=-1
        for j in range(0,xsize-slct):
            temp = ind[w[-1]]
            ind[w[-1]] = indx[slct+j]
            indx[slct+j] = temp
            diff = np.abs(y[ind]-func_val(tc,x[ind],function))
            if np.median(diff) < mad:
                inds = j
                mad = np.median(diff)
            # Switch it back
            temp = ind[w[-1]]
            ind[w[-1]] = indx[slct+j]
            indx[slct+j] = temp
        if inds == -1 or i>maxiter: break
        temp = ind[w[-1]]
        ind[w[-1]] = indx[slct+inds]
        indx[slct+inds] = temp
        i += 1
    return tc


def subsample(frame):
    newshape = (2*frame.shape[0], 2*frame.shape[1])
    slices = [slice(0, old, float(old)/new) for old, new in zip(frame.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')
    return frame[tuple(indices)]



def yamlify(obj, debug=False):
    """Recursively process an object so it can be serialised for yaml.
    Based on jsonify in `linetools <https://pypi.python.org/pypi/linetools>`_.

    Note: All string-like keys in :class:`dict` s are converted to
    :class:`str`.

    Also found in desiutils

    Parameters
    ----------
    obj : :class:`object`
        Any object.
    debug : :class:`bool`, optional
        Print extra information if requested.

    Returns
    -------
    :class:`object`
       An object suitable for yaml serialization.  For example
       :class:`numpy.ndarray` is converted to :class:`list`,
       :class:`numpy.int64` is converted to :class:`int`, etc.
    """
    if isinstance(obj, (np.float64, np.float32)):
        obj = float(obj)
    elif isinstance(obj, (np.int32, np.int64, np.int16)):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
#    elif isinstance(obj, bytes):
#        obj = obj.decode('utf-8')
    elif isinstance(obj, (np.string_, str)):
        obj = str(obj)
    elif isinstance(obj, units.Quantity):
        try:
            obj = obj.value.tolist()
        except AttributeError:
            obj = obj.value
    elif isinstance(obj, np.ndarray):  # Must come after Quantity
        obj = obj.tolist()
    elif isinstance(obj, dict):
        # First convert keys
        nobj = {}
        for key, value in obj.items():
            if isinstance(key, str):
                nobj[str(key)] = value
            else:
                nobj[key] = value
        # Now recursive
        obj = nobj
        for key, value in obj.items():
            obj[key] = yamlify(value, debug=debug)
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            obj[i] = yamlify(item, debug=debug)
    elif isinstance(obj, tuple):
        obj = list(obj)
        for i, item in enumerate(obj):
            obj[i] = yamlify(item, debug=debug)
        obj = tuple(obj)
    # elif isinstance(obj, Unit):
    #     obj = obj.name
    # elif obj is units.dimensionless_unscaled:
    #     obj = 'dimensionless_unit'
    if debug:
        print(type(obj))
    return obj

###########
def fit_min(xarr, yarr, xguess, width=None):

    errcode = 0
    # Edges
    if width is None:
        xleft, xright = np.min(xarr), np.max(xarr)
    else:
        xleft = xguess - width
        xright = xguess + width
    idx = np.where((xarr >= xleft) & (xarr <= xright))[0]

    # Setup
    thisx = xarr[idx]
    thisy = yarr[idx]

    # Guess for Gaussian
    guess = np.max(thisy), 0., width/2.

    # Fit with Gaussian
    try:
        coeff = func_fit(thisx-xguess, thisy, 'gaussian', 3, guesses=guess)
    except RuntimeError:  # Bad fit
        errcode = -1
        return xguess, 0., errcode
    sigma = coeff[2]
    xbest = xguess + coeff[1]

    # Could/should add a bunch of sanity checks
    # Insist on it being a minimum
    if coeff[0] > 0.:
        errcode = -4
    if (xbest < xleft) or (xbest > xright):
        errcode = -6
    # Return
    return xbest, sigma, errcode


def find_nminima(yflux, xvec=None, nfind=10, nsmooth=None, minsep=5, width=5):
    """ Find minima in an input 1D array
    Parameters
    ----------
    yflux : ndarray
    xvec : ndarray, optional
      Assumed to be ascending
    nfind : int, optional
      Number of peaks to find in the input array
    nsmooth : int, optional
      Smooth by a Gaussian with kenrel of nsmooth
    minsep : int, optional
      Minimum separation between peaks
    width : int, optional
      Width around a putative peak to fit a Gaussian

    Returns
    -------
    peaks: ndarray
      x values of the peaks 
    sigmas: ndarray
      sigma widths of the Gaussian fits to each peak
    ledges: ndarray
      left edges of each peak;  defined to be at least minsep away
      from the peak and where the slope of the data switches 
    redges: ndarray
      right edges of each peak;  defined to be at least minsep away
      from the peak and where the slope of the data switches 
    """
    # Init
    if xvec is None:
        xvec = np.arange(len(yflux))
    # Gaussian smooth
    if nsmooth is not None:
        yflux = convolve(yflux, Gaussian1DKernel(nsmooth))#, **kwargs)

    # ycopy, yderiv, ydone
    ycopy = yflux.copy()
    yderiv = np.roll(ycopy,1)-ycopy
    yderiv[0] = 0.
    yderiv[-1] = 0.
    ydone = np.max(ycopy)

    # Find first one
    peaks, sigmas, ledges, redges = [], [], [], []
    npeak = 0
    for kk in range(nfind):
        imin = np.argmin(ycopy)
        xbest, sigma, errcode = fit_min(xvec, ycopy, xvec[imin], width=width)
        #
        noldpeak = npeak
        npeak = len(peaks)
        # Find edges and
        # Block out pixels within minsep and 2*minsep
        x1 = (xvec < xvec[imin]-minsep) & (np.roll(yderiv,1) < 0.)
        if np.any(x1):
            ix1 = np.where(x1)[0][-1]
        else:
            ix1 = 0
        x2 = (xvec > xvec[imin]+minsep) & (yderiv > 0.)  # Scans until increasing
        if np.any(x2):
            ix2 = np.where(x2)[0][0]
        else:
            ix2 = len(xvec)
        ycopy[ix1:ix2] = ydone
        # Save
        if npeak == 0:  # Always grab at least one
            peaks.append(xbest)
            sigmas.append(sigma)
            ledges.append(ix1)
            redges.append(ix2-1)
        else:  # Check it is minsep away (seems like it will always be)
            xmin = np.min(np.abs(np.array(peaks-xbest)))
            if (xmin > minsep) & (errcode >= 0):
                peaks.append(xbest)
                sigmas.append(sigma)
                ledges.append(ix1)
                redges.append(ix2-1)
        # Any more to look for?
        if not np.any(ycopy < ydone):
            npeak = nfind
    return np.array(peaks), np.array(sigmas), np.array(ledges), np.array(redges)


def unravel_specobjs(specobjs):
    """
    Method to unwrap nested specobjs objects into a single list

    Parameters
    ----------
    specobjs : list of lists or list of SpecObj

    Returns
    -------
    all_specobj : list of SpecObj

    """
    # Wrapped is all None and lists
    ans = [isinstance(ispec, (list, type(None))) for ispec in specobjs]
    if np.all(ans):
        all_specobj = []
        for det in range(len(specobjs)):           # detector loop
            if specobjs[det] is None:
                continue
            for sl in range(len(specobjs[det])):   # slit loop
                for spobj in specobjs[det][sl]:    # object loop
                    all_specobj.append(spobj)
    else:
        all_specobj = specobjs
    # Return
    return all_specobj
