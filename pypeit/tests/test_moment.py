import pytest

import numpy as np

from scipy.special import erf

from pypeit.core import moment

def test_basics():
    c = [45,50,55]
    img = np.zeros((len(c),100), dtype=float)
    x = np.arange(100)
    sig = 5.
    for i,_c in enumerate(c):
        img[i,:] = (erf((x-c[i]+0.5)/np.sqrt(2)/sig) 
                    - erf((x-c[i]-0.5)/np.sqrt(2)/sig))/2.

    # Calculate all moments at one column and row
    mu, mue, flag = moment.moment1d(img, 50, 40., row=0, order=[0,1,2])
    assert np.allclose(mu, np.array([0.99858297, 45.02314924,  4.97367636]))
    
    # Calculate all moments at one column for all rows
    assert np.allclose(moment.moment1d(img, 50, 40., order=[0,1,2])[0],
                       np.array([[ 0.99858297,  0.99993125,  0.99858297],
                                 [45.02314924, 50.        , 54.97685076],
                                 [ 4.97367636,  5.00545947,  4.97367636]]))
    
    # Calculate zeroth moments in all rows centered at column 50
    assert np.allclose(moment.moment1d(img, 50, 40., order=0)[0],
                       np.array([0.99858297, 0.99993125, 0.99858297]))

    # Calculate zeroth moments in all rows for three column positions
    assert np.allclose(moment.moment1d(img, [45,50,55], 40., order=0, mesh=True)[0],
                       np.array([[0.99993125, 0.99858297, 0.97670951],
                                 [0.99858297, 0.99993125, 0.99858297],
                                 [0.97670951, 0.99858297, 0.99993125]]))

    # Calculate zeroth moments in each row with one column center per row
    assert np.allclose(moment.moment1d(img, [45,50,55], 40., row=[0,1,2], order=0)[0],
                       np.array([0.99993125, 0.99993125, 0.99993125]))

    # Calculate the first moment in a column for all rows
    assert np.allclose(moment.moment1d(img, 50, 40., row=[0,1,2], order=1)[0],
                       np.array([45.02314924, 50.        , 54.97685076]))

    # Or pick a column unique to each row
    assert np.allclose(moment.moment1d(img, [43,52,57], 40., row=[0,1,2], order=1)[0],
                       np.array([44.99688181, 50.00311819, 55.00311819]))


def test_bounded():
    c = [45,50,55]
    img = np.zeros((len(c),100), dtype=float)
    x = np.arange(100)
    sig = 5.
    for i,_c in enumerate(c):
        img[i,:] = (erf((x-c[i]+0.5)/np.sqrt(2)/sig) 
                    - erf((x-c[i]-0.5)/np.sqrt(2)/sig))/2.

    # Should find first moment is less than 49, meaning that the moment
    # is set to the input and flagged
    mu, mue, flag = moment.moment1d(img, 50, 40., row=0, order=[0,1,2],
                                    bounds=([0.9,-1.,4],[1.1,1.,6]))
    assert mu[1] == 50. and flag[1]


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


