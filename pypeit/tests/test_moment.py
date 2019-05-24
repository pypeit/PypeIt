import pytest

import numpy as np

from scipy.special import erf

from pypeit import moment

def test_center_uniform():
    """
    Test the centering algorithm with uniform weighting
    """

    ngau = 10
    npix = 200

    # Make a Gaussian comb with some random subpixel offsets
    xi = np.arange(ngau,npix,npix//ngau, dtype=float)
    xt = xi + np.round(np.random.normal(scale=1., size=xi.size), decimals=1)
    img = np.zeros((1,npix), dtype=float)
    x = np.arange(npix) 
    for c in xt:
        img[0,:] += (erf((x-c+0.5)/np.sqrt(2)/2) - erf((x-c-0.5)/np.sqrt(2)/2))/2.

    #xr, xe, bad = new_trace.recenter_moment(img, xi, ycen=np.zeros(xi.size), width=10.)
    xr, xe, bad = moment.moment1d(img, xi, 10., row=0, order=1)

    assert np.mean(np.absolute(xi-xt)) - np.mean(np.absolute(xr-xt)) > 0, \
                'Recentering did not improve the mean difference'
    assert np.std(xi-xt) - np.std(xr-xt) > 0, \
                'Recentering did not improve the standard deviation in the difference'

def test_center_gaussian():
    """
    Test the centering algorithm with Gaussian weighting
    """
    ngau = 10
    npix = 200

    # Make a Gaussian comb with some random subpixel offsets
    xi = np.arange(ngau,npix,npix//ngau, dtype=float)
    xt = xi + np.round(np.random.normal(scale=1., size=xi.size), decimals=1)
    img = np.zeros((1,npix), dtype=float)
    x = np.arange(npix) 
    for c in xt:
        img[0,:] += (erf((x-c+0.5)/np.sqrt(2)/2) - erf((x-c-0.5)/np.sqrt(2)/2))/2.

    #xr, xe, bad = new_trace.recenter_moment(img, xi, ycen=np.zeros(xi.size), width=4.,
    #                                        weighting='gaussian')
    xr, xe, bad = moment.moment1d(img, xi, 4., row=0, weighting='gaussian', order=1)

    assert np.mean(np.absolute(xi-xt)) - np.mean(np.absolute(xr-xt)) > 0, \
                'Recentering did not improve the mean difference'
    assert np.std(xi-xt) - np.std(xr-xt) > 0, \
                'Recentering did not improve the standard deviation in the difference'

def test_extract():
    """
    Test a simple zeroth moment flux extraction
    """
    ngau = 10
    npix = 2000
    sig = 10.
    delt = sig

    # Make a Gaussian comb with some random subpixel offsets
    xi = np.arange(npix//ngau//2,npix,npix//ngau, dtype=float)
    xt = xi + np.round(np.random.normal(scale=1., size=xi.size), decimals=1)
    img = np.zeros((1,npix), dtype=float)
    x = np.arange(npix) 
    for c in xt:
        img[0,:] += (erf((x-c+0.5)/np.sqrt(2)/sig) - erf((x-c-0.5)/np.sqrt(2)/sig))/2.

    #xr, xe, bad = new_trace.extract_aperture(img, xt-delt, xt+delt, ycen=np.zeros(xi.size))
    xr, xe, bad = moment.moment1d(img, xt, delt*2, row=0, order=0)

    truth = (erf(delt/np.sqrt(2)/sig) - erf(-delt/np.sqrt(2)/sig))/2.
    assert np.mean(np.absolute(xr-truth)) < 1e-3, 'Extraction inaccurate'


def test_width():
    """
    Test the measurement of the second moment
    """
    ngau = 10
    npix = 2000
    sig = 10.
    delt = 3*sig

    # Make a Gaussian comb with some random subpixel offsets
    xi = np.arange(npix//ngau//2,npix,npix//ngau, dtype=float)
    xt = xi + np.round(np.random.normal(scale=1., size=xi.size), decimals=1)
    img = np.zeros((1,npix), dtype=float)
    x = np.arange(npix) 
    for c in xt:
#        img[0,:] += np.exp(-np.square((x-c)/sig)/2)
        img[0,:] += (erf((x-c+0.5)/np.sqrt(2)/sig) - erf((x-c-0.5)/np.sqrt(2)/sig))/2.

    xr, xe, bad = moment.moment1d(img, xt, delt*2, row=0, order=2)

    assert np.absolute(np.mean(xr/sig)-1) < 0.02, 'Second moment should be good to better than 2%'


