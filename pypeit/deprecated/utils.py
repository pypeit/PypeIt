
import itertools
import matplotlib

import numpy as np

from scipy.optimize import curve_fit
from scipy import interpolate

from astropy import units
from astropy.io import fits
from astropy.convolution import convolve, Gaussian1DKernel
from matplotlib import pyplot as plt

# Imports for fast_running_median
from collections import deque
from itertools import islice
from bisect import insort, bisect_left

#from pydl.pydlutils import math
#from pydl.pydlutils import bspline



from pypeit.core import pydl

from pypeit import msgs


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


def func_fit(x, y, func, deg, x2=None, minx=None, maxx=None, minx2=None, maxx2=None,
             w=None, inmask=None, guesses=None, return_errors=False):
    """

    Args:
        x (`numpy.ndarray`_):
        y (`numpy.ndarray`_):
        func (:obj:`str`):
            polynomial, legendre, chebyshev, gaussian
        deg (:obj:`int`):
            degree of the fit
        x2 (`numpy.ndarray`_, optional):
        minx:
        maxx:
        minx2:
        maxx2:
        w (`numpy.ndarray`_, optional):
        inmask (`numpy.ndarray`_, optional):
        guesses:
        return_errors:

    Returns:
        PypeItFit:

    """
    if return_errors:
        msgs.error("Need to deal with this")
    # Init
    pcov = None

    # If the user provided an inmask apply it. The logic below of evaluating the fit only at the non-masked
    # pixels is preferable to the other approach of simply setting the weights to zero. The reason for that is that
    # the fits use a least-square optimization approach using matrix algebra, and lots of zero weights are
    # 1) more costly, and 2) will not produce exactly the same result (due to roundoff error) as actually
    # removing the locations you want to mask.
    if inmask is not None:
        x_out = x[inmask]
        y_out = y[inmask]
        if x2 is not None:
            x2_out = x2[inmask]
        else:
            x2_out = None
        if w is not None:
            w_out = w[inmask]
        else:
            w_out = None
    else:
        x_out = x
        y_out = y
        if x2 is not None:
            x2_out = x2
        else:
            x2_out = None
        if w is not None:
            w_out = w
        else:
            w_out = None

    # For two-d fits x = x, y = x2, y = z
    if ('2d' in func) and (x2_out is not None):
        # Is this a 2d fit?
        fitc = polyfit2d_general(x_out, x2_out, y_out, deg, w=w_out, function=func[:-2],minx=minx, maxx=maxx, miny=minx2, maxy=maxx2)
    elif func == "polynomial":
        fitc = np.polynomial.polynomial.polyfit(x_out, y_out, deg, w=w_out)
    elif func == "legendre":
        if minx is None or maxx is None:
            if np.size(x_out) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x_out), np.max(x_out)
        else:
            xmin, xmax = minx, maxx
        xv = 2.0 * (x_out-xmin)/(xmax-xmin) - 1.0
        fitc = np.polynomial.legendre.legfit(xv, y_out, deg, w=w_out)
    elif func == "chebyshev":
        if minx is None or maxx is None:
            if np.size(x_out) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x_out), np.max(x_out)
        else:
            xmin, xmax = minx, maxx
        xv = 2.0 * (x_out-xmin)/(xmax-xmin) - 1.0
        fitc = np.polynomial.chebyshev.chebfit(xv, y_out, deg, w=w_out)
    elif func == "bspline":
        msgs.error("Need to update whatever is calling this to call bspline instead..")
        #if bspline_par is None:
        #    bspline_par = {}
        ## TODO -- Deal with this kwargs-like kludge
        #fitc = bspline_fit(x_out, y_out, order=deg, w=w_out, **bspline_par)
    elif func == "gaussian":
        # Guesses
        if guesses is None:
            ampl, cent, sigma = guess_gauss(x_out, y_out)
            # As first guess choose slope and intercept to be zero
            b = 0
            m = 0
        else:
            if deg == 2:
                ampl, sigma = guesses
            elif deg == 3:
                ampl, cent, sigma = guesses
            elif deg == 4:
                b, ampl, cent, sigma = guesses
            elif deg == 5:
                m, b, ampl, cent, sigma = guesses
        # Error
        if w_out is not None:
            sig_y = 1./w_out
        else:
            sig_y = None
        if deg == 2:  # 2 parameter fit
            fitc, pcov = curve_fit(gauss_2deg, x_out, y_out, p0=[ampl, sigma], sigma=sig_y)
        elif deg == 3:  # Standard 3 parameters
            fitc, pcov = curve_fit(gauss_3deg, x_out, y_out, p0=[ampl, cent, sigma],
                                   sigma=sig_y)
        elif deg == 4:  # 4 parameters
            fitc, pcov = curve_fit(gauss_4deg, x_out, y_out, p0=[b, ampl, cent, sigma],sigma=sig_y)
        elif deg == 5:  # 5 parameters
            fitc, pcov = curve_fit(gauss_5deg, x_out, y_out, p0=[m, b, ampl, cent, sigma],sigma=sig_y)
        else:
            msgs.error("Not prepared for deg={:d} for Gaussian fit".format(deg))
    elif func == "moffat":
        # Guesses
        if guesses is None:
            ampl, cent, sigma = guess_gauss(x_out, y_out)
            p0 = ampl
            p2 = 3. # Standard guess
            p1 = (2.355*sigma)/(2*np.sqrt(2**(1./p2)-1))
        else:
            p0,p1,p2 = guesses
        # Error
        if w_out is not None:
            sig_y = 1./w_out
        else:
            sig_y = None
        if deg == 3:  # Standard 3 parameters
            fitc, pcov = curve_fit(moffat, x_out, y_out, p0=[p0,p1,p2], sigma=sig_y)
        else:
            msgs.error("Not prepared for deg={:d} for Moffat fit".format(deg))
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet" + msgs.newline() +
                   "Please choose from 'polynomial', 'legendre', 'chebyshev','bspline'")
    # DataContainer
    pypeitFit = PypeItFit(xval=x, yval=y, weights=w, fitc=fitc, fitcov=pcov, func=func,
                          minx=minx, maxx=maxx, minx2=minx2, maxx2=maxx2)
    if inmask is not None:
        pypeitFit.gpm = inmask.astype(int)
    return pypeitFit


def bspline_fit(x,y,order=3,knots=None,everyn=20,xmin=None,xmax=None,w=None,bkspace=None):
    """ bspline fit to x,y
    Should probably only be called from func_fit

    Parameters
    ----------
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

    Returns
    -------
    tck : tuple
        describes the bspline
    """
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


# JFH Testing
def zerocross1d(x, y, getIndices=False):
    """
      Find the zero crossing points in 1d data.

      Find the zero crossing events in a discrete data set.
      Linear interpolation is used to determine the actual
      locations of the zero crossing between two data points
      showing a change in sign. Data point which are zero
      are counted in as zero crossings if a sign change occurs
      across them. Note that the first and last data point will
      not be considered whether or not they are zero.

      Parameters
      ----------
      x, y : arrays
          Ordinate and abscissa data values.
      getIndices : boolean, optional
          If True, also the indicies of the points preceding
          the zero crossing event will be returned. Defeualt is
          False.

      Returns
      -------
      xvals : array
          The locations of the zero crossing events determined
          by linear interpolation on the data.
      indices : array, optional
          The indices of the points preceding the zero crossing
          events. Only returned if `getIndices` is set True.
    """

    # Check sorting of x-values
    if np.any((x[1:] - x[0:-1]) <= 0.0):
        msgs.error("The x-values must be sorted in ascending order!. Sort the data prior to calling zerocross1d.")

    # Indices of points *before* zero-crossing
    indi = np.where(y[1:] * y[0:-1] < 0.0)[0]

    # Find the zero crossing by linear interpolation
    dx = x[indi + 1] - x[indi]
    dy = y[indi + 1] - y[indi]
    zc = -y[indi] * (dx / dy) + x[indi]

    # What about the points, which are actually zero
    zi = np.where(y == 0.0)[0]
    # Do nothing about the first and last point should they
    # be zero
    zi = zi[np.where((zi > 0) & (zi < x.size - 1))]
    # Select those point, where zero is crossed (sign change
    # across the point)
    zi = zi[np.where(y[zi - 1] * y[zi + 1] < 0.0)]

    # Concatenate indices
    zzindi = np.concatenate((indi, zi))
    # Concatenate zc and locations corresponding to zi
    zz = np.concatenate((zc, x[zi]))

    # Sort by x-value
    sind = np.argsort(zz)
    zz, zzindi = zz[sind], zzindi[sind]

    if not getIndices:
        return zz
    else:
        return zz, zzindi


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


def poly_to_gauss(coeffs):
    try:
        sigm = np.sqrt(-0.5/coeffs[2])
        cent = -0.5*coeffs[1]/coeffs[2]
        ampl = np.exp( coeffs[0] + 0.5*(cent/sigm)**2 )
    except:
        return [0.0, 0.0, 0.0], True
    return [ampl, cent, sigm], False

def robust_regression(x, y, ordr, outfrac, maxiter=100, function='polynomial', min=None, max=None):
    """
    Deprecated
    """
    msgs.bug("PypeIt using deprecated function")
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
