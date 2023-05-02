"""
Module to run tests on SpecObj
"""
import numpy as np
import sys
import os
from copy import deepcopy
import pytest

from IPython import embed

from astropy.table import Table
from astropy.io import fits

from pypeit import spec2dobj
from pypeit.spectrographs.util import load_spectrograph
from pypeit.tests import tstutils
from pypeit import slittrace
from pypeit import pypmsgs
from pypeit.images import imagebitmask


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
    return dict(sciimg = sciimg,
                ivarraw = 0.1 * np.ones_like(sciimg),
                skymodel = 0.95 * np.ones_like(sciimg),
                objmodel = np.ones_like(sciimg),
                ivarmodel = 0.05 * np.ones_like(sciimg),
                scaleimg = np.ones_like(sciimg),
                waveimg = 1000 * np.ones_like(sciimg),
                bpmmask=imagebitmask.ImageBitMaskArray(sciimg.shape),
                slits=slits,
                wavesol=None,
                maskdef_designtab=None,
                tilts=np.ones_like(sciimg).astype(float),
                #tilts=wavetilts.WaveTilts(**test_wavetilts.instant_dict),
                sci_spat_flexure=3.5,
                sci_spec_flexure=spec_flex_table,
                vel_type='HELIOCENTRIC',
                vel_corr=1.0+1.0e-5)

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
    init_dict['detector'] = tstutils.get_kastb_detector()
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    # Check
    assert spec2DObj.hdu_prefix == 'DET01-'


def test_spec2dobj_io(init_dict):
    init_dict['detector'] = tstutils.get_kastb_detector()
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    # Write
    ofile = tstutils.data_path('tst_spec2d.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    spec2DObj.to_file(ofile)
    # Read
    _spec2DObj = spec2dobj.Spec2DObj.from_file(ofile, spec2DObj.detname)
    os.remove(ofile)


def test_spec2dobj_update_slit(init_dict):
    # Build two
    spec2DObj1 = spec2dobj.Spec2DObj(**init_dict,
                                     detector=load_spectrograph('keck_deimos').get_detector_par(1))
    spec2DObj2 = spec2dobj.Spec2DObj(**init_dict,
                                     detector=load_spectrograph('keck_deimos').get_detector_par(2))

    # WARNING: The instantiation of the two objects above using the same
    # dictionary means their components point to *the same objects*.  I.e.,
    # ``spec2DObj1.sciimg is spec2DObj2.sciimg`` is True!  That means any
    # alterations to complex attributes of spec2DObj1 are also made to
    # spec2DObj2.  I made some changes below to ensure this isn't true for this
    # test, but we need to beware of any instantiations like the above in the
    # code itself!
    
    # Checks
    with pytest.raises(pypmsgs.PypeItError):
        spec2DObj1.update_slits(spec2DObj2)

    # Update
    spec2DObj2.detector = load_spectrograph('keck_deimos').get_detector_par(1)
    spec2DObj2.sciimg = spec2DObj1.sciimg.copy()*2.
    spec2DObj2.slits = deepcopy(init_dict['slits'])
    spec2DObj2.slits.mask[1:] = 1
    # WARNING: This barfs!!
    # spec2DObj2.bpmmask = deepcopy(init_dict['bpmmask'])
    spec2DObj2.bpmmask = imagebitmask.ImageBitMaskArray(init_dict['sciimg'].shape)
    spec2DObj2.bpmmask[...] = 1

    spec2DObj1.update_slits(spec2DObj2)


####################################################3
# Testing of AllSpec2DObj

def test_all2dobj_hdr(init_dict):
    # Build one
    init_dict['detector'] = tstutils.get_kastb_detector()
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['bkg_redux'] = False
    allspec2D['meta']['find_negative'] = False
    allspec2D[spec2DObj.detname] = spec2DObj
    #
    kast_file = tstutils.data_path('b1.fits.gz')
    header = fits.getheader(kast_file)
    spectrograph = load_spectrograph('shane_kast_blue')
    # Do it
    hdr = allspec2D.build_primary_hdr(header, spectrograph, calib_dir=tstutils.data_path(''))
    # Test it
    assert hdr['SKYSUB'] == 'MODEL'


def test_all2dobj_write(init_dict):
    # Build one
    init_dict['detector'] = tstutils.get_kastb_detector()
    spec2DObj = spec2dobj.Spec2DObj(**init_dict)
    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['bkg_redux'] = False
    allspec2D['meta']['find_negative'] = False
    detname = spec2DObj.detname
    allspec2D[detname] = spec2DObj
    # Write
    ofile = tstutils.data_path('tst_allspec2d.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    allspec2D.write_to_fits(ofile)
    # Read
    _allspec2D = spec2dobj.AllSpec2DObj.from_fits(ofile)
    # Check
    assert allspec2D.detectors == _allspec2D.detectors, 'Bad read: detector mismatch'
    assert allspec2D['meta'] == _allspec2D['meta'], 'Bad read: meta mismatch'
    # Try to update it
    _allspec2D['meta']['bkg_redux'] = True
    _allspec2D[detname].vel_corr = 2.
    _allspec2D.write_to_fits(ofile, update_det='DET01')

    __allspec2D = spec2dobj.AllSpec2DObj.from_fits(ofile)
    assert __allspec2D['meta'] == _allspec2D['meta'], 'Bad read: meta mismatch'
    assert __allspec2D['meta'] != allspec2D['meta'], 'Bad read: meta mismatch'
    assert __allspec2D[detname].vel_corr == 2., 'Bad update'
    os.remove(ofile)


def test_all2dobj_update_image(init_dict):

    allspec2D = spec2dobj.AllSpec2DObj()
    allspec2D['meta']['bkg_redux'] = False
    allspec2D['meta']['find_negative'] = False
    for i in range(2):
        d = load_spectrograph('keck_deimos').get_detector_par(i+1)
        allspec2D[d.name] = spec2dobj.Spec2DObj(detector=d, **init_dict)

    # Write
    ofile = tstutils.data_path('tst_allspec2d.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    allspec2D.write_to_fits(ofile)

    _allspec2D = spec2dobj.AllSpec2DObj()
    _allspec2D['meta']['bkg_redux'] = False
    _allspec2D['meta']['find_negative'] = False
    d = load_spectrograph('keck_deimos').get_detector_par(2)
    detname = d.name
    init_dict['sciimg'] = allspec2D[detname].sciimg.copy() * 2
    _allspec2D[detname] = spec2dobj.Spec2DObj(detector=d, **init_dict)

    _allspec2D.write_to_fits(ofile, update_det=_allspec2D.detectors, overwrite=True)

    # Check
    allspec2D_2 = spec2dobj.AllSpec2DObj.from_fits(ofile)
    assert np.array_equal(allspec2D_2[detname].sciimg, _allspec2D[detname].sciimg), 'Bad update'
    assert np.array_equal(allspec2D_2[detname].sciimg, allspec2D[detname].sciimg*2), 'Bad update'

    os.remove(ofile)

