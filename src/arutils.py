import os
import astropy.io.fits as pyfits
import itertools
import numpy as np
import armsgs as msgs
import arcyutils
import arcyarc
#import matplotlib.pyplot as plt
#import arplot

try:
    import ds9
except ImportError:
    msgs.warn('ds9 module not installed')
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

def erf(x):
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

def func_der(coeffs,func,nderive=1):
    if func == "polynomial":
        return np.polynomial.polynomial.polyder(coeffs,m=nderive)
    elif func == "legendre":
        return np.polynomial.legendre.legder(coeffs,m=nderive)
    elif func == "chebyshev":
        return np.polynomial.chebyshev.chebder(xv,y,deg)
    else:
        msgs.error("Functional derivative '{0:s}' is not implemented yet"+msgs.newline()+"Please choose from 'polynomial', 'legendre', 'chebyshev'")

def func_fit(x,y,func,deg,min=None,max=None):
    if func == "polynomial":
        return np.polynomial.polynomial.polyfit(x,y,deg)
    elif func == "legendre":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legfit(xv,y,deg)
    elif func == "chebyshev":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebfit(xv,y,deg)
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet"+msgs.newline()+"Please choose from 'polynomial', 'legendre', 'chebyshev'")

def func_val(c,x,func,min=None,max=None):
    if func == "polynomial":
        return np.polynomial.polynomial.polyval(x,c)
    elif func == "legendre":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legval(xv,c)
    elif func == "chebyshev":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebval(xv,c)
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet"+msgs.newline()+"Please choose from 'polynomial', 'legendre', 'chebyshev'")

def func_vander(x,func,deg,min=None,max=None):
    if func == "polynomial":
        return np.polynomial.polynomial.polyvander(x,deg)
    elif func == "legendre":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.legendre.legvander(xv,deg)
    elif func == "chebyshev":
        if min is None or max is None:
            if np.size(x) == 1:
                xmin, xmax = -1.0, 1.0
            else:
                xmin, xmax = np.min(x), np.max(x)
        else:
            xmin, xmax = min, max
        xv = 2.0 * (x-xmin)/(xmax-xmin) - 1.0
        return np.polynomial.chebyshev.chebvander(xv,deg)
    else:
        msgs.error("Fitting function '{0:s}' is not implemented yet"+msgs.newline()+"Please choose from 'polynomial', 'legendre', 'chebyshev'")


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


def perturb(covar, bparams, nsim=1000):
    cvsize = covar.shape[0]
    # Generate a new set of starting parameters from the covariance matrix
    X_covar_fit=np.matrix(np.random.standard_normal((cvsize,nsim)))
    C_covar_fit=np.matrix(covar)
    U_covar_fit=np.linalg.cholesky(C_covar_fit)
    Y_covar_fit=U_covar_fit * X_covar_fit
    newpar = bparams.reshape(cvsize,1).repeat(nsim,axis=1)
    # Find the new best-fitting model parameters
    newpar += Y_covar_fit
    return newpar


def polyfitter2d(data,mask=None,order=2):
    x, y = np.meshgrid(np.linspace(0.0,1.0,data.shape[1]),np.linspace(0.0,1.0,data.shape[0]))
    if mask is None or mask.size==0:
        xf = x.flatten()
        yf = y.flatten()
        m = polyfit2d(xf,yf,data.T.flatten(),order)
    else:
        mskar = np.ones((data.shape[0],data.shape[1]))
        mskar[mask,:] = 0
        w = np.where(mskar==1)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf,yf,data[w].T.flatten(),order)
    # Return the best model
    return m, polyval2d(x,y,m).T
    #print m
    #model = polyval2d(x,y,m)
    ## Plot it
    #zmin, zmax = arplot.zscale(data)
    #plt.subplot(131)
    #implt = plt.imshow(data,aspect='auto',interpolation='none')
    #plt.colorbar()
    #implt.set_clim(zmin, zmax)
    #plt.subplot(132)
    #implt = plt.imshow(model,aspect='auto',interpolation='none')
    #plt.colorbar()
    #implt.set_clim(zmin, zmax)
    #plt.subplot(133)
    #zmin, zmax = arplot.zscale(data-model)
    #implt = plt.imshow(data-model,aspect='auto',interpolation='none')
    #plt.colorbar()
    #implt.set_clim(zmin, zmax)
    #plt.show()
    #return m, polyval2d(x,y,m).T

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
    for a, (i,j) in zip(m, ij):
        z += a * x**i * y**j
    return z

def gauss_fit(x,y,pcen):
#	dx = np.ones(x.size)*np.mean(x[1:]-x[:-1])
#	coeffs = polyfit_integral(x, y, dx, 2)
#	return poly_to_gauss(coeffs)
    if np.any(y<0.0):
        return [0.0, 0.0, 0.0], True
    ampl, cent, sigm, good = arcyarc.fit_gauss(x, y, np.zeros(3,dtype=np.float), 0, x.size, float(pcen))
    if good == 0:
        return [0.0, 0.0, 0.0], True
    else:
        return [ampl, cent, sigm], False


def poly_to_gauss(coeffs):
    sigm = np.sqrt(-0.5/coeffs[2])
    cent = -0.5*coeffs[1]/coeffs[2]
    ampl = np.exp( coeffs[0] + 0.5*cent/sigm**2 )
    return [ampl, cent, sigm]

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
            raise TypeError("expected 1D vector for w")
        if len(x) != len(w):
            raise TypeError("expected x and w to have same length")
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


def robust_polyfit(xarray, yarray, order, maxone=True, sigma=3.0, function="polynomial", initialmask=None, forceimask=False, min=None, max=None):
    """
    A robust (equally weighted) polynomial fit is performed to the xarray, yarray pairs
    mask[i] = 1 are masked values
    """
    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size,dtype=np.int)
        if forceimask:
            msgs.warn("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    mskcnt=np.sum(mask)
    # Iterate, and mask out new values on each iteration
    while True:
        w = np.where(mask==0)
        xfit = xarray[w]
        yfit = yarray[w]
        ct = func_fit(xfit,yfit,function,order,min=min,max=max)
        yrng = func_val(ct,xarray,function,min=min,max=max)
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
            if forceimask:
                w = np.where((np.abs(yarray-yrng) > sigma*sigmed) | (initialmask==1))
            else:
                w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break # No new values have been included in the mask
        mskcnt = np.sum(mask)
        w = np.where(mask==0)
    xfit = xarray[w]
    yfit = yarray[w]
    ct = func_fit(xfit,yfit,function,order,min=min,max=max)
    return mask, ct


def robust_regression(x,y,ordr,outfrac,maxiter=100,function='polynomial',min=None,max=None):
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
