"""
Module to run tests on ararclines
"""
import os

import numpy as np
import pytest

from pypeit import utils
from pypeit import msgs

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)




def test_calc_ivar():
    """ Run the parameter setup script
    """
    x = np.array([-1.0, -0.1, 0.0, 0.1, 1.0])
    res = utils.inverse(x)
    assert np.array_equal(res, np.array([0.0, 0.0, 0.0, 10.0, 1.0]))
    assert np.array_equal(utils.calc_ivar(res), np.array([0.0, 0.0, 0.0, 0.1, 1.0]))


def test_nearest_unmasked():
    arr = np.ma.MaskedArray(np.arange(10))
    arr[3] = np.ma.masked
    arr[8] = np.ma.masked
    nearest = utils.nearest_unmasked(arr)
    assert np.array_equal(nearest, np.array([1, 0, 1, 2, 5, 4, 5, 6, 7, 7])), \
            'Closest indices did not match expected result'
    assert np.array_equal(nearest, utils.nearest_unmasked(arr, use_indices=True)), \
            'Result should be independent of use_indices for this array' 


def test_boxcar_smooth_rows():
    # Build a test image ...
    nrows = 31
    ncols = 11
    nave = 11
    img = np.zeros((nrows,ncols), dtype=float)
    img[nrows//2-3:nrows//2+4,:] = 1.
    img[0,:] = 1.
    img[-1,:] = 1.
    # ... and a good pixel mask
    gpm = np.ones(img.shape, dtype=float)
    gpm[nrows//2,:] = 0.

    # Use the function both without ...
    smimg = utils.boxcar_smooth_rows(img, nave)
    # ... and with the mask
    smmimg = utils.boxcar_smooth_rows(img, nave, wgt=gpm)

    # Setup for a brute-force calculation
    #   - Image with repeated rows
    _img = np.zeros((nrows+2*nave,ncols), dtype=float)
    _img[nave:nrows+nave,:] = img
    _img[:nave,:] = img[0,None,:]
    _img[nrows+nave:,:] = img[-1,None,:]
    #   - good pixel mask
    _gpm = np.zeros((nrows+2*nave,ncols), dtype=float)
    _gpm[nave:nrows+nave,:] = gpm
    _gpm[:nave,:] = gpm[0,None,:]
    _gpm[nrows+nave:,:] = gpm[-1,None,:]
    #   - weighted image
    _wimg = _gpm * _img
    #   - image used for averaging
    left = np.arange(nrows+nave)
    right = np.arange(nrows+nave)+nave
    pix = np.arange(nrows+2*nave)
    avg = (pix[:,None] >= left[None,:]) & (pix[:,None] < right[None,:])

    # Perform a brute-force calculation w/ and w/o the gpm
    _smimg = np.zeros(img.shape, dtype=float)
    _smmimg = np.zeros(img.shape, dtype=float)
    for j in range(ncols):
        m = np.sum(avg * _img[:,None,j], axis=0)/nave
        _smimg[:,j] = m[nave//2+1:-nave//2+1]

        m = np.sum(avg * _wimg[:,None,j], axis=0)/np.sum(avg * _gpm[:,None,j], axis=0)
        _smmimg[:,j] = m[nave//2+1:-nave//2+1]

    # Should be the same within the numerical precision.  Test here is
    # much larger than that.
    assert np.allclose(smimg, _smimg), 'Difference with brute-force approach unmasked.'
    assert np.allclose(smmimg, _smmimg), 'Difference with brute-force approach masked.'



