"""
Module to run tests on methods in utils.py
"""
import os

from IPython import embed

import yaml

import numpy as np

from pypeit import utils
from pypeit import msgs
from pypeit.tests.tstutils import data_path
from pypeit import io


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



def test_yamlify():
    """ This tests the yamlify method and also the approach to 
    writing and reading the Setup block of PypeIt"""

    obj = dict(a=1., b='acb', datasec='[2:23,:2048]', d=dict(c=3))

    new_obj = utils.yamlify(obj)

    # Write
    tst_file = data_path('tst.yaml')
    with open(tst_file, 'w') as f:
        setup_lines = io.dict_to_lines(new_obj, level=1)
        f.write('\n'.join(setup_lines)+'\n')

    # Read
    with open(tst_file, 'r') as f:
        lines = f.readlines()

    # Strip white space
    lines = [line.strip() for line in lines]
    # Add back in \n
    ystr = '\n'.join(lines)
    sdict = yaml.safe_load(ystr)

    # Clean up
    os.remove(tst_file)


def test_add_sub_dict():
    d = {}
    utils.add_sub_dict(d, 'test')
    assert d == {'test': {}}, 'add_sub_dict failure'
    d['test'] = 'this'
    utils.add_sub_dict(d, 'test')
    assert d == {'test': 'this'}, 'add_sub_dict failure'
    utils.add_sub_dict(d, 'and')
    d['and'] = 'that'
    assert d == {'test': 'this', 'and': 'that'}, 'add_sub_dict failure'


def test_recursive_update():
    d = {}
    d['rdx'] = dict(spectrograph='shane_kast_blue')
    u = {}
    u['rdx'] = dict(detnum=3)

    d = utils.recursive_update(d, u)
    assert sorted(list(d['rdx'].keys())) == ['detnum', 'spectrograph'], 'Missing merged keys'

