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
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests import tstutils

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

####################################################3
# Testing of Spec2DObj


def test_init(init_dict):
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    # Check
    assert spec2DObj.hdu_prefix == 'DET01-'


def test_spec2dobj_io(init_dict):
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    spec2DObj.detector = tstutils.get_kastb_detector()
    # Write
    ofile = data_path('tst_spec2d.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    spec2DObj.to_file(ofile)
    # Read
    _spec2DObj = spec2dobj.Spec2DObj.from_file(ofile, init_dict['det'])
    os.remove(ofile)


####################################################3
# Testing of AllSpec2DObj

def test_all2dobj_hdr(init_dict):
    # Build one
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['ir_redux'] = False
    allspec2D[1] = spec2DObj
    #
    kast_file = data_path('b1.fits.gz')
    header = fits.getheader(kast_file)
    spectrograph = load_spectrograph('shane_kast_blue')
    # Do it
    hdr = allspec2D.build_primary_hdr(header, spectrograph, master_dir=data_path(''))
    # Test it
    assert hdr['SKYSUB'] == 'MODEL'


def test_all2dobj_write(init_dict):
    # Build one
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['ir_redux'] = False
    allspec2D[1] = spec2DObj
    allspec2D[1].detector = tstutils.get_kastb_detector()
    # Write
    ofile = data_path('tst_allspec2d.fits')
    allspec2D.write_to_fits(ofile)
    # Read
    _allspec2D = spec2dobj.AllSpec2DObj.from_fits(ofile)

