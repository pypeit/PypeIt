from __future__ import absolute_import, division, print_function, unicode_literals

import os
import astropy.io.fits as pyfits
from astropy.stats import sigma_clip as sigma_clip
from scipy.optimize import curve_fit
from scipy.special import erf
from scipy import interpolate
import itertools
import numpy as np
from pypit import armsgs
#from pypit import arcyutils
#from pypit import arcyarc
import warnings

# Logging
msgs = armsgs.get_logger()

from pypit import ardebug as debugger

try:
    basestring
except NameError:  # For Python 3
    basestring = str

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
    """Trim to the inner knots.  Used in bspline_magfit
    Might be useful elsewhere
    """
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
    #
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
        debugger.set_trace()
    return tck


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
    """
    Parameters
    ----------
    nfile : int, optional
      Number of files to mimic
    spectrograph : str, optional
      Name of spectrograph to mimic

    Returns
    -------

    """
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    fitsdict['utc'] = ['2015-01-23']*nfile
    fitsdict['directory'] = [directory]*nfile
    fitsdict['filename'] = ['b{:03d}.fits'.format(i) for i in range(nfile)]
    fitsdict['date'] = ['2015-01-23T00:{:02d}:11.04'.format(i) for i in range(nfile)]  # Will fail at 60
    fitsdict['time'] = [(1432085758+i*60)/3600. for i in range(nfile)]
    fitsdict['target'] = ['Dummy']*nfile
    fitsdict['ra'] = ['00:00:00']*nfile
    fitsdict['dec'] = ['+00:00:00']*nfile
    fitsdict['exptime'] = [300.] * nfile
    fitsdict['naxis0'] = [2048] * nfile
    fitsdict['naxis1'] = [2048] * nfile
    fitsdict['dispname'] = ['600/4310'] * nfile
    fitsdict['dichroic'] = ['560'] * nfile
    fitsdict['dispangle'] = ['none'] * nfile
    fitsdict["binning"] = ['1x1']*nfile
    fitsdict["airmass"] = [1.0]*nfile
    #
    if spectrograph == 'shane_kast_blue':
        fitsdict['numamplifiers'] = [1] * nfile
        fitsdict['naxis0'] = [2112] * nfile
        fitsdict['naxis1'] = [2048] * nfile
        fitsdict['slitwid'] = [1.] * nfile
        fitsdict['slitlen'] = ['none'] * nfile
        # Lamps
        for i in range(1,17):
            fitsdict['lampstat{:02d}'.format(i)] = ['off'] * nfile
        fitsdict['exptime'][0] = 0        # Bias
        fitsdict['lampstat06'][1] = 'on'  # Arc
        fitsdict['exptime'][1] = 30       # Arc
        fitsdict['lampstat01'][2] = 'on'  # Trace, pixel, slit flat
        fitsdict['lampstat01'][3] = 'on'  # Trace, pixel, slit flat
        fitsdict['exptime'][2] = 30     # flat
        fitsdict['exptime'][3] = 30     # flat
        fitsdict['ra'][4] = '05:06:36.6'  # Standard
        fitsdict['dec'][4] = '52:52:01.0'
        fitsdict['airmass'][4] = 1.2
        fitsdict['ra'][5] = '07:06:23.45' # Random object
        fitsdict['dec'][5] = '+30:20:50.5'
        fitsdict['decker'] = ['0.5 arcsec'] * nfile
    elif spectrograph == 'none':
        pass
    # arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])
    # Return
    return fitsdict


def dummy_self(inum=0, fitsdict=None, nfile=10):
    """ Generate a dummy self class for testing
    Parameters:
    -----------
    inum : int, optional
      Index in sciexp
    Returns:
    --------
    slf
    """
    from pypit import arsciexp
    # Dummy fitsdict
    if fitsdict is None:
        fitsdict = dummy_fitsdict(nfile=nfile)
    # Dummy Class
    slf = arsciexp.ScienceExposure(inum, fitsdict, do_qa=False)
    return slf


def dummy_settings(pypitdir=None, nfile=10, spectrograph='shane_kast_blue',
                   set_idx=True):
    """ Generate settings classes
    Parameters
    ----------
    pypitdir
    nfile
    spectrograph
    set_idx : bool, optional
      Set dummy index values for science and calibs

    Returns
    -------

    """
    from pypit import arparse
    # Dummy argflag
    if spectrograph != 'shane_kast_blue':
        msgs.error("Only setup for Kast Blue")  # You will need to fuss with scidx
    argf = arparse.get_argflag_class(("ARMLSD", spectrograph))
    argf.init_param()
    if pypitdir is None:
        pypitdir = __file__[0:__file__.rfind('/')]
    # Run specific
    argf.set_param('run pypitdir {0:s}'.format(pypitdir))
    argf.set_param('run spectrograph {:s}'.format(spectrograph))
    argf.set_param('run directory science ./')
    # Dummy spect
    spect = arparse.get_spect_class(("ARMLSD", spectrograph, "dummy"))
    lines = spect.load_file(base=True)  # Base spectrograph settings
    spect.set_paramlist(lines)
    lines = spect.load_file()
    spect.set_paramlist(lines)
    if set_idx:
        for jj, key in enumerate(spect._spect.keys()):
            if key in ['det']:
                continue
            if 'index' in spect._spect[key].keys():
                if spectrograph == 'shane_kast_blue':  # Science frames from idx = 5 to 9
                    assert nfile == 10
                for kk in [5,6,7,8,9]:
                    if key == 'science':
                        spect._spect[key]['index'] += [np.array([kk])]
                    elif key == 'arc':
                        spect._spect[key]['index'] += [np.array([1])]
                    elif key == 'standard':
                        spect._spect[key]['index'] += [np.array([4])]
                    elif key == 'bias':
                        spect._spect[key]['index'] += [np.array([0])]
                    elif key == 'trace':
                        spect._spect[key]['index'] += [np.array([2,3])]
                    elif key == 'pixelflat':
                        spect._spect[key]['index'] += [np.array([2,3])]
    arparse.init(argf, spect)
    return


def dummy_specobj(fitsdict, det=1, extraction=True):
    """ Generate dummy specobj classes
    Parameters
    ----------
    fitsdict : dict
      Expecting the fitsdict from dummy_fitsdict
    Returns
    -------

    """
    from astropy import units as u
    from pypit import arspecobj
    shape = fitsdict['naxis1'][0], fitsdict['naxis0'][0]
    config = 'AA'
    scidx = 5 # Could be wrong
    xslit = (0.3,0.7) # Center of the detector
    ypos = 0.5
    xobjs = [0.4, 0.6]
    specobjs = []
    for xobj in xobjs:
        specobj = arspecobj.SpecObjExp(shape, config, scidx, det, xslit, ypos, xobj)
        # Dummy extraction?
        if extraction:
            npix = 2001
            specobj.boxcar['wave'] = np.linspace(4000., 6000., npix)*u.AA
            specobj.boxcar['counts'] = 50.*(specobj.boxcar['wave'].value/5000.)**-1.
            specobj.boxcar['var']  = specobj.boxcar['counts'].copy()
        # Append
        specobjs.append(specobj)
    # Return
    return specobjs

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
    kwargs

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


def gauss_lsqfit(x,y,pcen):
    """

    :param x:
    :param y:
    :param pcen: An estimate of the Gaussian mean
    :return:
    """
    from pypit import arcyarc
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
    from pypit import arcyarc
    try:
        if np.any(y<0.0):
            return [0.0, 0.0, 0.0], True
        if x.size <= 3:
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
            try:
                msgs.warn("More parameters than data points - fit might be undesirable")
            except AttributeError:
                print("More parameters than data points - fit might be undesirable")
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


def trace_gweight(fimage, xcen, ycen, sigma, invvar=None, maskval=-999999.9):
    """ Determines the trace centroid by weighting the flux by the integral
    of a Gaussian over a pixel
    Port of SDSS trace_gweight algorithm

    Parameters
    ----------
    fimage : ndarray
      image to centroid on
    xcen : ndarray
      guess of centroids in x (column) dimension
    ycen : ndarray (usually int)
      guess of centroids in y (rows) dimension
    sigma : float
      Width of gaussian
    invvar : ndarray, optional
    maskval : float, optional
      Value for masking

    Returns
    -------
    xnew : ndarray
      New estimate for trace in x-dimension
    xerr : ndarray
      Error estimate for trace.  Rejected points have maskval

    """
    # Setup
    nx = fimage.shape[1]
    xnew = np.zeros_like(xcen)
    xerr = maskval*np.ones_like(xnew)

    if invvar is None:
        invvar = np.ones_like(fimage)

    # More setting up
    x_int = np.round(xcen).astype(int)
    nstep = 2*int(3.0*sigma) - 1

    weight = np.zeros_like(xcen)
    numer  = np.zeros_like(xcen)
    meanvar = np.zeros_like(xcen)
    bad = np.zeros_like(xcen).astype(bool)

    for i in range(nstep):
        xh = x_int - nstep//2 + i
        xtemp = (xh - xcen - 0.5)/sigma/np.sqrt(2.0)
        g_int = (erf(xtemp+1./sigma/np.sqrt(2.0)) - erf(xtemp))/2.
        xs = np.minimum(np.maximum(xh,0),(nx-1))
        cur_weight = fimage[ycen, xs] * (invvar[ycen, xs] > 0) * g_int * ((xh >= 0) & (xh < nx))
        weight += cur_weight
        numer += cur_weight * xh
        meanvar += cur_weight * cur_weight * (xcen-xh)**2 / (
                            invvar[ycen, xs] + (invvar[ycen, xs] == 0))
        bad = np.any([bad, xh < 0, xh >= nx], axis=0)

    # Masking
    good = (~bad) & (weight > 0)
    if np.sum(good) > 0:
        xnew[good] = numer[good]/weight[good]
        xerr[good] = np.sqrt(meanvar[good])/weight[good]
    # Return
    return xnew, xerr

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
    import numpy as np
    from astropy.units import Quantity

    if isinstance(obj, (np.float64, np.float32)):
        obj = float(obj)
    elif isinstance(obj, (np.int32, np.int64, np.int16)):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
    elif isinstance(obj, (np.string_, basestring)):
        obj = str(obj)
    elif isinstance(obj, Quantity):
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
            if isinstance(key, basestring):
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
    # elif obj is u.dimensionless_unscaled:
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
    # Imports
    from astropy.convolution import convolve, Gaussian1DKernel
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

def get_atm_template(theta, trans):
    """ Scales the original transmission spectrum
    Parameters
    ----------
    theta: float
        Scale factor
    trans: ndarray
        Atm. transmission spectrum

    Returns
    -------
    template: ndarray
    """
    return 1 + theta*(trans - 1)

def lnlike_tf(theta, data, trans):
    """ Liklihood function to maximize the fit between
    the data and the template
    Parameters
    ----------
    theta: float
        Scale factor
    data: 3D ndarray
    tran: 2D ndaarray
    Returns
    -------
    liklihood: float
    """
    i = ((data[:,0] >= 9320) & (data[:,0] <= 9380))
    tmp = get_atm_template(theta, trans[:,1])

    return (-0.5 * np.sum(np.log(2 * np.pi * data[:,2][i]**2)
            + (data[:,1][i] - tmp[i])**2 / data[:,2][i]**2))

def opposite_lnlike_tf(theta, data, trans):
    """ Opposite of the liklihood function, to be used with
    minimization algorithms
    Parameters
    ----------
    theta: float
        Scale factor
    data: 3D ndarray
    tran: 2D ndaarray
    Returns
    -------
    opposite liklihood: float
    """
    i = ((data[:,0] >= 9320) & (data[:,0] <= 9380))
    return -1. * lnlike_tf(theta, data, trans)

def get_fwhm_pix(sp, inval=10000., outval=150):
    """ Find FWHM (pix) to convolve a spec1D by
    Parameters
    ----------
    sp:     XSpectrum1D object
        Spectrum to be convolved
    inval:  int
        Resolution of sp (R)
    outval: float
        Resolution to convolve sp to (km/s)
    Returns
    -------
    fwhm_dff: ndarray
        Array of FWHM (pix) values to convolve
        a spectrum by
    """
    clight = 299792.458

    # Get FWHM of spectrum being smoothed
    fwhm_in_ = sp.wavelength.value/inval
    psize    = np.diff(sp.wavelength.value)
    psize    = np.append(psize,
                    sp.wavelength[-1].value - sp.wavelength[-2].value)

    # Get FWHM of resolution thaat spectrum is to be smoothed to
    fwhm_out = outval*2.*np.sqrt(2*np.log(2))
    fwhm_out = (fwhm_out / clight) * sp.wavelength.value

    fwhm_diff = (fwhm_out - fwhm_in_)/psize

    return fwhm_diff

