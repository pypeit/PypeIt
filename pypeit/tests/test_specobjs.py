"""
Module to run tests on SpecObjs
"""
import os

from IPython import embed

import pytest

import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit import specobjs
from pypeit import specobj
from pypeit import io
from pypeit.tests import tstutils


@pytest.fixture
def sobj1():
    return specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
@pytest.fixture
def sobj2():
    return specobj.SpecObj('MultiSlit', 'DET02', SLITID=1)
@pytest.fixture
def sobj3():
    return specobj.SpecObj('MultiSlit', 'DET03', SLITID=0)
@pytest.fixture
def sobj4():
    return specobj.SpecObj('MultiSlit', 'DET01', SLITID=10)


def test_init(sobj1, sobj2):
    """ Run the parameter setup script
    """
    # Null
    sobjs1 = specobjs.SpecObjs()

    # With a few objs
    sobjs2 = specobjs.SpecObjs([sobj1,sobj2])
    assert sobjs2.nobj == 2


def test_access(sobj1, sobj2):
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    #
    assert sobjs[0]['PYPELINE'] == 'MultiSlit'
    assert len(sobjs['PYPELINE']) == 2

def test_add_rm(sobj1, sobj2, sobj3):
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    sobjs.add_sobj(sobj3)
    assert sobjs.nobj == 3
    # Remove
    sobjs.remove_sobj(2)
    assert len(sobjs.specobjs) == 2

    # Numpy 18
    sobjs1 = specobjs.SpecObjs()
    sobjs2 = specobjs.SpecObjs()
    sobjs2.add_sobj(sobjs1)


def test_set(sobj1, sobj2, sobj3):
    sobjs = specobjs.SpecObjs([sobj1,sobj2,sobj3])
    # All
    sobjs.DET = 'DET03'
    assert np.all(sobjs[:].DET == np.array(['DET03','DET03','DET03']))
    sobjs[:].DET = 'DET04'
    assert np.all(sobjs[:].DET == np.array(['DET04','DET04','DET04']))
    # Slice
    sobjs[1:2].DET = 'DET02'
    assert sobjs.DET[1] == 'DET02'
    # With logic
    det2 = sobjs.DET == 'DET02'
    sobjs[det2].PYPELINE = 'BLAH'
    assert sobjs.PYPELINE[1] == 'BLAH'
    assert sobjs.PYPELINE[0] == 'MultiSlit'


def test_io(sobj1, sobj2, sobj3, sobj4):
    sobjs = specobjs.SpecObjs([sobj1,sobj2,sobj3,sobj4])
    sobjs[0]['BOX_WAVE'] = np.arange(1000).astype(float)
    sobjs[1]['BOX_WAVE'] = np.arange(1000).astype(float)
    sobjs[2]['BOX_WAVE'] = np.arange(1000).astype(float)
    #sobjs[0]['BOX_COUNTS'] = np.ones_like(sobjs[0].BOX_WAVE)  # This tests single array
    sobjs[1]['BOX_COUNTS'] = np.ones_like(sobjs[0].BOX_WAVE)
    sobjs[2]['BOX_COUNTS'] = np.ones_like(sobjs[0].BOX_WAVE)
    # Detector
    sobjs[0]['DETECTOR'] = tstutils.get_kastb_detector()
    tmp = tstutils.get_kastb_detector()

    tmp['det'] = 2
    sobjs[1]['DETECTOR'] = tmp
    # Write
    header = fits.PrimaryHDU().header
    header['TST'] = 'TEST'
    ofile = tstutils.data_output_path('tst_specobjs.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)
    sobjs.write_to_fits(header, ofile, overwrite=False)
    # Read
    hdul = io.fits_open(ofile)
    assert len(hdul) == 7  # Primary + 4 Obj + 2 Detectors
    assert hdul[0].header['NSPEC'] == 4
    hdul.close()
    #
    _sobjs = specobjs.SpecObjs.from_fitsfile(ofile)
    assert _sobjs.nobj == 4
    assert np.array_equal(sobjs[0].BOX_WAVE, _sobjs[0].BOX_WAVE)
    assert np.array_equal(sobjs[1].BOX_WAVE, _sobjs[1].BOX_WAVE)
    _sobjs.write_to_fits(header, ofile, overwrite=True)

    # Detector
    assert _sobjs[0].DETECTOR is not None, '1st object started with Detector'
    assert _sobjs[1].DETECTOR is not None, '2nd object has DET=1 so should get decorated'
    assert _sobjs[2].DETECTOR is None

    # Now try updates!
    sobjs1 = specobjs.SpecObjs([sobj1])
    sobjs[0]['BOX_WAVE'] = np.arange(2000).astype(float)
    sobjs1[0]['DETECTOR'] = tstutils.get_kastb_detector()
    header1 = fits.PrimaryHDU().header
    sobjs1.write_to_fits(header1, ofile, overwrite=True, update_det='DET01')

    # Test
    _sobjs1 = specobjs.SpecObjs.from_fitsfile(ofile)
    assert _sobjs1.nobj == 3
    assert _sobjs1[2].BOX_WAVE.size == 2000
    os.remove(ofile)


def _sobjs_with_offsets(dets, first=0):
    """
    SpecObjs with one named object per detector and a table of fiber offsets.
    """
    sobjs = specobjs.SpecObjs()
    for i, det in enumerate(dets):
        sobj = specobj.SpecObj('MultiSlit', det, SLITID=first+i)
        sobj.BOX_WAVE = np.arange(10).astype(float)
        sobj.NAME = f'SPAT{first+i:04d}-{det}'
        sobjs.add_sobj(sobj)
    sobjs.fiber_offsets = Table(dict(
        NAME=[sobj.NAME for sobj in sobjs], TARGETX_ASEC=np.arange(len(dets)) + 0.1 + first
    ))
    return sobjs


def test_fiber_offsets_io(tmp_path):
    ofile = tmp_path / 'spec1d.fits'
    hdr = fits.Header()
    sobjs = _sobjs_with_offsets(['DET01', 'DET02'])
    sobjs.write_to_fits(hdr, str(ofile))
    # Read and rewrite, as done by, e.g., the flux calibration
    for _ in range(2):
        with io.fits_open(ofile) as hdul:
            assert 'FIBER_OFFSETS' in hdul, 'the offsets must be in their own extension'
            assert hdul[0].header['NSPEC'] == 2, 'the table must not be counted as a spectrum'
        _sobjs = specobjs.SpecObjs.from_fitsfile(ofile)
        assert _sobjs.nobj == 2, 'the table must not be read as a spectrum'
        assert list(_sobjs.fiber_offsets['NAME']) == list(sobjs.NAME), \
            'names must match the spectra after reading'
        np.testing.assert_allclose(
            _sobjs.fiber_offsets['TARGETX_ASEC'], sobjs.fiber_offsets['TARGETX_ASEC'],
            err_msg='offsets must be unchanged after reading'
        )
        _sobjs.write_to_fits(hdr, str(ofile))


def test_fiber_offsets_kept_by_slice_and_copy():
    sobjs = _sobjs_with_offsets(['DET01', 'DET02'])
    for new in [sobjs[sobjs.DET == 'DET01'], sobjs.copy()]:
        assert new.fiber_offsets is not None, 'slices and copies must keep the fiber offsets'


def test_fiber_offsets_only_written_for_spectra(tmp_path):
    ofile = tmp_path / 'spec1d.fits'
    sobjs = _sobjs_with_offsets(['DET01', 'DET02'])
    sobjs[sobjs.DET == 'DET01'].write_to_fits(fits.Header(), str(ofile))
    tab = specobjs.SpecObjs.from_fitsfile(ofile).fiber_offsets
    assert list(tab['NAME']) == ['SPAT0000-DET01'], \
        'offsets must only be written for spectra in the file'


def test_fiber_offsets_update_det(tmp_path):
    ofile = tmp_path / 'spec1d.fits'
    hdr = fits.Header()
    _sobjs_with_offsets(['DET01', 'DET02']).write_to_fits(hdr, str(ofile))
    _sobjs_with_offsets(['DET02'], first=5).write_to_fits(hdr, str(ofile), update_det='DET02')
    tab = specobjs.SpecObjs.from_fitsfile(ofile).fiber_offsets
    assert list(tab['NAME']) == ['SPAT0000-DET01', 'SPAT0005-DET02'], \
        'updating a detector must replace its offsets and keep the others'
    np.testing.assert_allclose(
        tab['TARGETX_ASEC'], [0.1, 5.1], err_msg='offsets must stay with their fibers'
    )


def test_no_fiber_offsets_extension_by_default(tmp_path, sobj1):
    sobj1.BOX_WAVE = np.arange(10).astype(float)
    ofile = tmp_path / 'spec1d.fits'
    specobjs.SpecObjs([sobj1]).write_to_fits(fits.Header(), str(ofile))
    with io.fits_open(ofile) as hdul:
        assert 'FIBER_OFFSETS' not in hdul, 'files without fiber offsets must be unchanged'
    assert specobjs.SpecObjs.from_fitsfile(ofile).fiber_offsets is None, \
        'fiber_offsets must be None for files without the table'
