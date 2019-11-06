"""
Module to run tests loading
"""
import os
import pytest
from pypeit import specobjs
from pypeit import specobj

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)



def test_load_specobjs():
    spec_file = data_path('spec1d_r153-J0025-0312_KASTr_2015Jan23T025323.850.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(spec_file)

    # Test
    assert isinstance(sobjs, specobjs.SpecObjs)
    assert len(sobjs[0].BOX_COUNTS) == 1200

    assert isinstance(sobjs[0], specobj.SpecObj)

