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


# IF THIS TEST IS FAILING BECAUSE THE SPEC1D DATAMODEL UPDATED,
#   RUN python copy_spec1d.py in files/
#   *After* you have run the Dev Suite tests
def test_load_specobjs():
    spec_file = data_path('spec1d_r153-J0025-0312_KASTr_2015Jan23T025323.850.fits')
    sobjs = specobjs.SpecObjs.from_fitsfile(spec_file)

    # Test
    assert isinstance(sobjs, specobjs.SpecObjs)
    assert len(sobjs[0].BOX_COUNTS) == 1200

    assert isinstance(sobjs[0], specobj.SpecObj)

