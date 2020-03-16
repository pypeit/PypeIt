"""
Module to run tests on SpecObj
"""
import numpy as np
import sys
import os
import pytest


#import pypeit

from astropy.table import Table
from astropy.io import fits

from pypeit import spec2dobj

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def init_dict():
    sciimg = np.ones((500,500)).astype(float)
    sdict = dict(sciimg = sciimg,
                 ivarraw = 0.1 * np.ones_like(sciimg),
                 skymodel = 0.95 * np.ones_like(sciimg),
                 objmodel = np.ones_like(sciimg),
                 ivarmodel = 0.05 * np.ones_like(sciimg),
                 mask = np.ones_like(sciimg).astype(int),
                 det = 1,
                 detector = None,
        )
    return sdict

def test_init(init_dict):
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    # Check
    assert spec2DObj.hdu_prefix == 'DET01-'

# Testing of AllSpec2DObj
def test_all2dobj_hdr(init_dict):
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    #


