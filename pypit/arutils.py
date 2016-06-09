import os
import astropy.io.fits as pyfits
from astropy.stats import sigma_clip as sigma_clip
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy import interpolate
import itertools
import numpy as np
import armsgs
try:
    import arcyutils
except ImportError:  # Travis testing only (Cython)
    pass
else:
    import arcyarc
import warnings

#from xastropy.xutils import xdebug as xdb

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

try:
    import ds9
except ImportError:
    warnings.warn("ds9 module not installed")
else:
    def ds9plot(array):
        # Set up a ds9 instance
        d = ds9.ds9()
        # Load the image
        d.set_np2arr(array)
        # Zoom to fit
        d.set('zoom to fit')
        # Change the colormap and scaling
        d.set('cmap gray')
        d.set('scale log')
        null = raw_input("TEST PLOT GENERATED: Press enter to continue...")
        return

def quicksave(data,fname):
    """
    Save a fits file (quickly) -- overwrite is forced, and no quality control checks
    """
    hdu = pyfits.PrimaryHDU(data)
    hdulist = pyfits.HDUList([hdu])
    if os.path.exists(fname):
        os.remove(fname)
    hdulist.writeto(fname)
    return

def bspline_inner_knots(all_knots):
    '''Trim to the inner knots.  Used in bspline_magfit
    Might be useful elsewhere
    '''
    diff = all_knots - np.roll(all_knots,1)
    pos = np.where(diff>0.)[0]
    i0=pos[0]
    i1=pos[-1]
    return all_knots[i0:i1]

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
    everyn: int 
      Knot everyn good pixels, if used
    bkspace: float 
      Spacing of breakpoints in units of x

    Returns:
    ---------
    fit_dict: dict  
      dict describing the bspline fit 
    ''' 
    #
    if w is None:
        ngd = x.size
        gd = np.arange(ngd)
        weights = None
    else:
        gd = np.where(w > 0.)[0]
        weights = w[gd]
    # Make the knots
    if knots is None:
        if bkspace is not None: 
            xrnge = (np.max(x[gd]) - np.min(x[gd]))
            startx = np.min(x[gd])
            nbkpts = max(int(xrnge/bkspace) + 1,2)
            tempbkspace = xrnge/(nbkpts-1)
            knots = np.arange(1,nbkpts-1)*tempbkspace + startx
        elif everyn is not None:
            idx_knots = np.arange(10, ngd-10, everyn) # A knot every good N pixels
            knots = x[gd[idx_knots]]
        else:
            msgs.error("No method specified to generate knots")
    # Generate spline
    try:
        tck = interpolate.splrep(x[gd], y[gd], w=weights, k=order, t=knots)
    except ValueError: # Knot problem
        debugger.set_trace()
    return tck


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


def dummy_self(pypitdir=None):
    """
    Generate a dummy self class for testing
    Parameters:
    -----------
    pypitdir : str, optional
      Path to the PYPIT main directory
    Returns:
    --------
    slf
    """
    # Dummy Class
    slf = type('Dummy', (object,), {"_argflag": {}, "_spect": {}})
    slf._argflag['run'] = {}
    if pypitdir is not None:
        slf._argflag['run']['pypitdir'] = pypitdir
    #
    slf._spect['mosaic'] = {}
    #
    return slf


def erf_func(x):
    """
    This is a very good approximation to the erf function.
    A probability is calculated as such:
    prob = erf( N / np.sqrt(2.0) )
    where N is the level of significance.
    e.g. erf( 1.0 / np.sqrt(2.0) ) = 0.68268947
    """
    # Constants
    a1 =  0.254829592
    a2 = -0.284496736
    a3 =  1.421413741
    a4 = -1.453152027
    a5 =  1.061405429
    p  =  0.3275911
    # Save the sign of x
    sign = np.ones(x.size)
    sign[np.where(x<0.0)] *= -1
    x = np.abs(x)
    # A&S formula 7.1.26
    t = 1.0/(1.0 + p*x)
    y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*np.exp(-x*x)
    return sign*y


def find_peaks(yval, siglev=10.):
    """
    Find peaks in a input array
    yerror is calculated from the array itself (RMS of clipped data)

    Parameters:
    -------------
    yval: ndarray
      Vector of y-values
    siglev: float, optional
      Sigma level for detection
    """
    # Calculate RMS
    yclipped = sig_clip(yval)
    rms = np.std(yclipped)
    #
    tpixt, num = arcyarc.detections_sigma(yval,np.array([rms]*mid_row.shape[0]),np.zeros(yval.shape[0],dtype=np.int),siglev/2.0,siglev) 
    # Remove similar
    pixt = arcyarc.remove_similar(tpixt, num)
    pixt = pixt[np.where(pixt!=-1)].astype(np.int)
    # Return
    return pixt


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
             **kwargs):
    """ General routine to fit a function to a given set of x,y points

    Parameters
    ----------
    x
    y
    func : str
      polynomial, legendre, chebyshev, bspline, gauss
    deg
    minv
    maxv
    w
    guesses
    kwargs

    Returns
    -------

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
        return bspline_fit(x, y, order=deg, w=w, **kwargs)
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
    c
    x
    func
    minv
    maxv

    Returns
    -------

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


def get_splknots(xarr, yarr, num, minv=None, maxv=None, maxknots=None):
    """
    Determine the best location for the knots of the scipy function LSQUnivariateSpline
    :param xarr, yarr: Input (x,y) arrays which are used to highlight were the strongest gradients are in the fitted
                        function. yarr should be a reduced (in size) and approximate representation of the data being
                        fitted, and xarr should be the corresponding x values.
    :param yarr: Input array which is used to highlight were the strongest gradients are in the fitted function.
                    This array should be a reduced (in size) and approximate representation of the data being fitted
    :param num: Number of knot locations to use
    :param minv: The minimum x-value of the knots
    :param maxv: The maximum x-value of the knots
    :return: knots
    """
    # First determine the derivative
    if minv is None: minv = np.min(xarr)
    if maxv is None: maxv = np.max(xarr)
    tck = interpolate.splrep(xarr, np.sqrt(np.abs(yarr)), xb=minv, xe=maxv, s=0)
    deriv = np.abs(interpolate.splev(xarr, tck, der=2))
    drvsum = np.cumsum(np.abs(deriv))
    drvsum *= num/drvsum[-1] # drvsum represents the cumulative number of knots to be placed at a given coordinate.
    # Now undo the cumulative sum
    drv = np.append(drvsum[0], drvsum[1:]-drvsum[:-1])
    drv = np.ceil(drv).astype(np.int)
    drv[np.where(drv<2)] = 2
    if maxknots is not None: drv[np.where(drv>maxknots)] = maxknots
    print drv
    knots = arcyutils.get_splknots(xarr, drv, minv, maxv, np.sum(drv))
    msgs.info("Generated {0:d} knots".format(knots.size))
    return knots


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
    if mask is None or mask.size == 0:
        xf = x.flatten()
        yf = y.flatten()
        m = polyfit2d(xf, yf, data.T.flatten(), order)
    else:
        mskar = np.ones((data.shape[0], data.shape[1]))
        mskar[mask,:] = 0
        w = np.where(mskar == 1)
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


def gauss_lsqfit(x,y,pcen):
    """

    :param x:
    :param y:
    :param pcen: An estimate of the Gaussian mean
    :return:
    """
    def gfunc(x,ampl,cent,sigm,cons,tilt):
        df = (x[1:]-x[:-1])/2.0
        df = np.append(df,df[-1])
        dff = (x[1:]**2 - x[:-1]**2)/2.0
        dff = np.append(dff,dff[-1])
        sqt = sigm*np.sqrt(2.0)
        return cons*df*2.0 + tilt*dff + ampl*0.5*np.sqrt(np.pi)*sqt*(erf((x+df-cent)/sqt) - erf((x-df-cent)/sqt))
        #return cons + ampl*np.exp(-0.5*((x-cent)/sigm)**2)

    if np.any(y<0.0):
        return [0.0, 0.0, 0.0], True
    # Obtain a quick first guess at the parameters
    ampl, cent, sigm, good = arcyarc.fit_gauss(x, y, np.zeros(3,dtype=np.float), 0, x.size, float(pcen))
    if good == 0:
        return [0.0, 0.0, 0.0], True
    elif np.any(np.isnan([ampl, cent, sigm])):
        return [0.0, 0.0, 0.0], True
    else:
        # Perform a least squares fit
        try:
            popt, pcov = curve_fit(gfunc, x, y, p0=[ampl, cent, sigm, 0.0, 0.0], maxfev=100)
            #popt, pcov = curve_fit(gfunc, x, y, p0=[0.0,ampl, cent, sigm], maxfev=100)
        except:
            return [0.0, 0.0, 0.0], True
        return [popt[0], popt[1], popt[2]], False


def gauss_fit(x, y, pcen):
    # dx = np.ones(x.size)*np.mean(x[1:]-x[:-1])
    # coeffs = polyfit_integral(x, y, dx, 2)
    # return poly_to_gauss(coeffs)
    try:
        if np.any(y<0.0):
            return [0.0, 0.0, 0.0], True
        ampl, cent, sigm, good = arcyarc.fit_gauss(x, y, np.zeros(3,dtype=np.float), 0, x.size, float(pcen))
        if good == 0:
            return [0.0, 0.0, 0.0], True
        elif np.any(np.isnan([ampl, cent, sigm])):
            return [0.0, 0.0, 0.0], True
        else:
            return [ampl, cent, sigm], False
    except:
        return [0.0, 0.0, 0.0], True


def poly_to_gauss(coeffs):
    try:
        sigm = np.sqrt(-0.5/coeffs[2])
        cent = -0.5*coeffs[1]/coeffs[2]
        ampl = np.exp( coeffs[0] + 0.5*(cent/sigm)**2 )
    except:
        return [0.0, 0.0, 0.0], True
    return [ampl, cent, sigm], False


def polyfit2d_general(x, y, z, deg, w=None):
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
    vander = np.polynomial.polynomial.polyvander2d(x, y, deg)
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            msgs.bug("arutils.polyfit2d - Expected 1D vector for weights")
        if len(x) != len(w) or len(y) != len(w) or len(x) != len(y):
            msgs.bug("arutils.polyfit2d - Expected x, y and weights to have same length")
        z = z * w
        vander = vander * w[:,np.newaxis]

    vander = vander.reshape((-1,vander.shape[-1]))
    z = z.reshape((vander.shape[0],))
    c = np.linalg.lstsq(vander, z)[0]
    return c.reshape(deg+1)


def polyval2d_general(c, x, y, function="polynomial", minx=None, maxx=None, miny=None, maxy=None):
    if function == "polynomial":
        xx, yy = np.meshgrid(x, y)
        return np.polynomial.polynomial.polyval2d(xx, yy, c)
    elif function in ["legendre", "chebyshev"]:
        # Scale x-direction
        if minx is None or maxx is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = minx, maxx
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        # Scale y-direction
        if miny is None or maxy is None:
            if np.size(y) == 1:
                ymin, ymax = -1.0, 1.0
            else:
                ymin, ymax = np.min(y), np.max(y)
        else:
            ymin, ymax = miny, maxy
        yv = 2.0 * (y-ymin)/(ymax-ymin) - 1.0
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
    x = np.asarray(x) + 0.0
    y = np.asarray(y) + 0.0

    # check arguments.
    if deg < 0 :
        raise ValueError("expected deg >= 0")
    if x.ndim != 1:
        raise TypeError("expected 1D vector for x")
    if x.size == 0:
        raise TypeError("expected non-empty vector for x")
    if y.ndim < 1 or y.ndim > 2 :
        raise TypeError("expected 1D or 2D array for y")
    if len(x) != len(y):
        raise TypeError("expected x and y to have same length")

    # set up the least squares matrices in transposed form
    lhst = np.polynomial.polynomial.polyvander(x+dx/2.0, deg+1) - np.polynomial.polynomial.polyvander(x-dx/2.0, deg+1)
    div = np.arange(1.,deg+2.).reshape(1,deg+1).repeat(x.size,axis=0)
    lhs = (lhst[:,1:]/(dx.reshape(dx.size,1).repeat(deg+1,axis=1)*div)).T
    rhs = y.T
    if w is not None:
        w = np.asarray(w) + 0.0
        if w.ndim != 1:
            msgs.bug("Expected 1D vector for weights in arutils.polyfit2d")
        if len(x) != len(w):
            msgs.bug("Expected x and weights to have same length in arutils.polyfit2d")
        # apply weights. Don't use inplace operations as they
        # can cause problems with NA.
        lhs = lhs * w
        rhs = rhs * w

    # set rcond
    if rcond is None :
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
        msg = "The fit may be poorly conditioned"
        warnings.warn(msg, pu.RankWarning)

    if full :
        return c, [resids, rank, s, rcond]
    else :
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
        for i in xrange(xfit.size):
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
             ['newshape[%d],factor[%d],'%(i,i) for i in xrange(lenShape)] + \
             [')'] + ['.sum(%d)'%(i+1) for i in xrange(lenShape)] + \
             ['/factor[%d]'%i for i in xrange(lenShape)]
    return eval(''.join(evList))


def robust_polyfit(xarray, yarray, order, weights=None, maxone=True, sigma=3.0,
                   function="polynomial", initialmask=None, forceimask=False,
                   minv=None, maxv=None, guesses=None, **kwargs):
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
                      guesses=ct, minv=minv, maxv=maxv, **kwargs)
        yrng = func_val(ct, xarray, function, minv=minv, maxv=maxv)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        if xarray.size-np.sum(mask) <= order+2:
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
    ct = func_fit(xfit, yfit, function, order, w=wfit, minv=minv, maxv=maxv, **kwargs)
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
        for j in xrange(0,xsize-slct):
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


def spline_coeffs(a, b, y, alpha=0.0, beta=0.0):
    """
    Return spline coefficients
    a : float
        lower bound of the grid.
    b : float
        upper bound of the grid.
    y : ndarray
        actual function value at grid points.
    c : (y.shape[0] + 2, ) ndarray, optional
        ndarry to be written
    alpha : float
        Second-order derivative at a. Default is 0.
    beta : float
        Second-order derivative at b. Default is 0.
    """
    n = y.shape[0] - 1
    h = (b - a)/n

    c = np.zeros(n + 3)
    c[1] = 1.0/6.0 * (y[0] - (alpha * h**2)/6.0)
    c[n + 1] = 1.0/6.0 * (y[n] - (beta * h**2)/6.0)

    # ab matrix here is just compressed banded matrix
    ab = np.ones((3, n - 1))
    ab[0, 0] = 0
    ab[1, :] = 4
    ab[-1, -1] = 0

    B = y[1:-1].copy()
    B[0] -= c[1]
    B[-1] -=  c[n + 1]

    c[2:-2] = np.linalg.solve(arcyutils.unband(ab, 1, 1), B)
    c[0] = alpha * h**2/6 + 2 * c[1] - c[2]
    c[-1] = beta * h**2/6 + 2 * c[-2] - c[-3]
    return c


def spline_interp(xnew,xold,yold):
    # Calculate the coefficients
    c = spline_coeffs(xold[0], xold[-1], yold)
    ynew = arcyutils.spline_interpolate(xnew, c, xold[0], xold[-1])
    return ynew


def subsample(frame):
    newshape = (2*frame.shape[0], 2*frame.shape[1])
    slices = [slice(0, old, float(old)/new) for old, new in zip(frame.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')
    return frame[tuple(indices)]
