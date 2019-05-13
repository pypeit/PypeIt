"""
Module to run tests on ampsec definition
"""
import pytest

import numpy as np

from scipy.special import erf

from pypeit import new_trace

def test_recenter_moment_uniform():
    """
    Test the recentering algorithm
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

    xr, xe, bad = new_trace.recenter_moment(img, xi, ycen=np.zeros(xi.size), width=10.)

    assert np.mean(np.absolute(xi-xt)) - np.mean(np.absolute(xr-xt)) > 0, \
                'Recentering did not improve the mean difference'
    assert np.std(xi-xt) - np.std(xr-xt) > 0, \
                'Recentering did not improve the standard deviation in the difference'

def test_recenter_moment_gaussian():
    """
    Test the recentering algorithm
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

    xr, xe, bad = new_trace.recenter_moment(img, xi, ycen=np.zeros(xi.size), width=4.,
                                            weighting='gaussian')

    assert np.mean(np.absolute(xi-xt)) - np.mean(np.absolute(xr-xt)) > 0, \
                'Recentering did not improve the mean difference'
    assert np.std(xi-xt) - np.std(xr-xt) > 0, \
                'Recentering did not improve the standard deviation in the difference'

def test_extract():
    """
    Test the flux extraction
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

    print(xt)
    xr, xe, bad = new_trace.extract(img, xt-delt, xt+delt, ycen=np.zeros(xi.size))

    test = (erf(delt/np.sqrt(2)/sig) - erf(-delt/np.sqrt(2)/sig))/2.
    import pdb; pdb.set_trace()
    
    recenter_moment(img, xi, ycen=np.zeros(xi.size), width=10.)

    assert np.mean(np.absolute(xi-xt)) - np.mean(np.absolute(xr-xt)) > 0, \
                'Recentering did not improve the mean difference'
    assert np.std(xi-xt) - np.std(xr-xt) > 0, \
                'Recentering did not improve the standard deviation in the difference'

if __name__ == '__main__':
    test_extract()

