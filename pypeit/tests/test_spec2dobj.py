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
from pypeit.tests import test_wavetilts
from pypeit import wavetilts
from pypeit import slittrace
from pypeit import pypmsgs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def init_dict():
    sciimg = np.ones((1000,1000)).astype(float)
    # Slits
    left = np.full((1000, 3), 2, dtype=float)
    right = np.full((1000, 3), 8, dtype=float)
    left[:,1] = 15.
    right[:,1] = 21.
    left[:,2] = 25.
    right[:,2] = 31.
    slits = slittrace.SlitTraceSet(left, right, 'MultiSlit',
                                   nspat=1000, PYP_SPEC='dummy')
    # Construct table of spectral flexure
    spec_flex_table = Table()
    spec_flex_table['spat_id'] = slits.spat_id
    spec_flex_table['sci_spec_flexure'] = np.zeros(left.shape[1])
    #
    sdict = dict(sciimg = sciimg,
                 ivarraw = 0.1 * np.ones_like(sciimg),
                 skymodel = 0.95 * np.ones_like(sciimg),
                 objmodel = np.ones_like(sciimg),
                 ivarmodel = 0.05 * np.ones_like(sciimg),
                 scaleimg = np.ones_like(sciimg),
                 waveimg = 1000 * np.ones_like(sciimg),
                 bpmmask=np.ones_like(sciimg).astype(int),
                 det=1,
                 detector=None,
                 slits=slits,
                 tilts=np.ones_like(sciimg).astype(float),
                 #tilts=wavetilts.WaveTilts(**test_wavetilts.instant_dict),
                 sci_spat_flexure=3.5,
                 sci_spec_flexure=spec_flex_table,
                 vel_type='HELIOCENTRIC',
                 vel_corr=1.0+1.0e-5
                 )
    return sdict

'''
from IPython import embed

dpath = '/home/xavier/Projects/PypeIt-development-suite/REDUX_OUT/keck_lris_blue/multi_300_5000_d680'
new_spec2dfile = os.path.join(dpath, 'Science', 'spec2d_b170816_0076-E570_LRISb_2017Aug16T071652.378.fits')
orig_spec2dfile = os.path.join(dpath, 'Science', 'Orig', 'spec2d_b170816_0076-E570_LRISb_2017Aug16T071652.378.fits')
new_spec2DObj = spec2dobj.Spec2DObj.from_file(new_spec2dfile, 1)
orig_spec2DObj = spec2dobj.Spec2DObj.from_file(orig_spec2dfile, 1)

orig_spec2DObj.update_slits(new_spec2DObj)
'''


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

def test_spec2dobj_update_slit(init_dict):
    # Build two
    spec2DObj1 = spec2dobj.Spec2DObj(**init_dict)
    spec2DObj2 = spec2dobj.Spec2DObj(**init_dict)

    # Checks
    spec2DObj2.det = 2
    with pytest.raises(pypmsgs.PypeItError):
        spec2DObj1.update_slits(spec2DObj2)

    # Update
    spec2DObj2.det = 1
    spec2DObj2.sciimg = spec2DObj1.sciimg.copy()*2.
    spec2DObj2.slits.mask[1:] = 1

    spec2DObj1.update_slits(spec2DObj2)

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
    if os.path.isfile(ofile):
        os.remove(ofile)
    allspec2D.write_to_fits(ofile)
    # Read
    _allspec2D = spec2dobj.AllSpec2DObj.from_fits(ofile)
    # Write again
    os.remove(ofile)
    _allspec2D.write_to_fits(ofile)

    os.remove(ofile)

def test_all2dobj_update_image(init_dict):
    # Build two
    spec2DObj1 = spec2dobj.Spec2DObj(**init_dict)
    spec2DObj2 = spec2dobj.Spec2DObj(**init_dict)
    spec2DObj2.det = 2
    #
    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['ir_redux'] = False
    allspec2D[1] = spec2DObj1
    allspec2D[2] = spec2DObj2

    # Write
    ofile = data_path('tst_allspec2d.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    allspec2D.write_to_fits(ofile)

    # Update
    _allspec2D = spec2dobj.AllSpec2DObj()
    spec2DObj1.sciimg = spec2DObj1.sciimg.copy()*2.
    _allspec2D['meta']['ir_redux'] = False
    _allspec2D[1] = spec2DObj1
    _allspec2D.write_to_fits(ofile, update_det=1, overwrite=True)

    # Check
    allspec2D_2 = spec2dobj.AllSpec2DObj.from_fits(ofile)
    assert np.array_equal(allspec2D_2[1].sciimg, spec2DObj1.sciimg)

    os.remove(ofile)
