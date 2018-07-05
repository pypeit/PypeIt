# Module to run tests on FlatField class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from astropy.table import Table

from pypit.core import arsort

from pypit.spectrographs import spectro_utils
from pypit import calibrations

from pypit.tests import tstutils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
    fitstbl = arsort.dummy_fitstbl(directory=data_path(''))
else:
    skip_test=False
    # MultiSlit
    fitstbl = arsort.dummy_fitstbl(
        directory=os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA/Shane_Kast_blue/600_4310_d55/'))
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{:d}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'


def chk_for_files(root):
    files = glob.glob(root+'*')
    if len(files) == 0:
        return False
    else:
        return True


@pytest.fixture
def multi_caliBrate():
    setup = 'A_01_aa'
    det = 1
    sci_ID = 1
    settings = tstutils.load_kast_blue_masters(get_settings=True)[0]
    multi_caliBrate= calibrations.MultiSlitCalibrations(fitstbl, save_masters=False, write_qa=False)
    #
    multi_caliBrate.reset(setup, det, sci_ID, settings)
    # Extra settings
    multi_caliBrate.settings['masters']['reuse'] = False
    multi_caliBrate.settings['bias'] = {}
    multi_caliBrate.settings['bias']['useframe'] = 'overscan'
    multi_caliBrate.settings['reduce'] = {}
    multi_caliBrate.settings['reduce']['badpix'] = False
    multi_caliBrate.settings['reduce']['trim'] = True
    multi_caliBrate.settings['reduce']['calibrate'] = {}
    multi_caliBrate.settings['reduce']['calibrate']['wavelength'] = True
    # Kludge me
    from pypit import traceslits
    multi_caliBrate.settings['trace'] = traceslits.default_settings()['trace']
    # Kludge me
    from pypit import wavecalib
    multi_caliBrate.settings['arc'] = {}
    multi_caliBrate.settings['arc']['calibrate'] = wavecalib.default_settings()['calibrate']
    # Yet another
    from pypit import wavetilts
    multi_caliBrate.settings['trace']['slits']['tilts'] = wavetilts.default_settings()['tilts'].copy()
    # Yup
    from pypit import flatfield
    from pypit import processimages
    multi_caliBrate.settings['reduce']['flatfield'] = flatfield.default_settings()['flatfield']
    multi_caliBrate.settings['reduce']['flatfield']['perform'] = True
    multi_caliBrate.settings['reduce']['flatfield']['useframe'] = 'pixelflat'
    multi_caliBrate.settings['reduce']['slitprofile'] = flatfield.default_settings()['slitprofile']
    multi_caliBrate.settings['pixelflat'] = {}
    multi_caliBrate.settings['pixelflat']['combine'] = processimages.default_settings()['combine']
    # Return
    return multi_caliBrate


def test_instantiate():
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)
    print(caliBrate)


def test_datasec(multi_caliBrate):
    if skip_test:
        assert True
        return
    datasec_img, naxis0, naxis1 = multi_caliBrate.get_datasec_img()
    # Test
    assert naxis0 == 2112
    assert naxis1 == 350

def test_pixlocn(multi_caliBrate):
    multi_caliBrate.shape = (2048,350)
    pixlocn = multi_caliBrate.get_pixlocn()
    # Test
    assert pixlocn.shape == (2048,350,4)

def test_bias(multi_caliBrate):
    #
    # Build
    multi_caliBrate.get_bias()


def test_arc(multi_caliBrate):
    multi_caliBrate.msbias = 'overscan'
    # Build
    arc = multi_caliBrate.get_arc()
    assert arc.shape == (2048,350)


def test_bpm(multi_caliBrate):
    # Prep
    multi_caliBrate.shape = (2048,350)
    # Build
    bpm = multi_caliBrate.get_bpm()
    assert bpm.shape == (2048,350)
    assert np.sum(bpm) == 0.


def test_slits(multi_caliBrate):
    if skip_test:
        assert True
        return
    # Setup
    multi_caliBrate.shape = (2048,350)
    _ = multi_caliBrate.get_pixlocn()
    _ = multi_caliBrate.get_datasec_img()
    _ = multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    # Run
    tslits_dict, maskslits = multi_caliBrate.get_slits()
    # Test
    assert isinstance(tslits_dict, dict)
    assert isinstance(maskslits, np.ndarray)


def test_wv_calib(multi_caliBrate):
    if skip_test:
        assert True
        return
    # Setup
    multi_caliBrate.shape = (2048,350)
    _ = multi_caliBrate.get_pixlocn()
    _ = multi_caliBrate.get_datasec_img()
    _ = multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    _ = multi_caliBrate.get_slits()
    _ = multi_caliBrate.get_arc()
    # Run
    wv_calib, maskslits = multi_caliBrate.get_wv_calib()
    assert isinstance(wv_calib, dict)
    assert isinstance(maskslits, np.ndarray)


def test_tilts(multi_caliBrate):
    if skip_test:
        assert True
        return
    # Setup
    multi_caliBrate.shape = (2048,350)
    _ = multi_caliBrate.get_pixlocn()
    _ = multi_caliBrate.get_datasec_img()
    _ = multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    _ = multi_caliBrate.get_slits()
    _ = multi_caliBrate.get_arc()
    _ = multi_caliBrate.get_wv_calib()
    # Run
    mstilts, maskslits = multi_caliBrate.get_tilts()
    assert mstilts.shape == (2048,350)
    assert isinstance(maskslits, np.ndarray)


def test_flat(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    _ = multi_caliBrate.get_pixlocn()
    _ = multi_caliBrate.get_datasec_img()
    _ = multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    _ = multi_caliBrate.get_slits()
    _ = multi_caliBrate.get_arc()
    _ = multi_caliBrate.get_wv_calib()
    _ = multi_caliBrate.get_tilts()
    # Run
    mspixflatnrm, slitprof = multi_caliBrate.get_pixflatnrm()
    assert mspixflatnrm.shape == (2048,350)
    assert slitprof.shape == (2048,350)


def test_waveimg(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    _ = multi_caliBrate.get_pixlocn()
    _ = multi_caliBrate.get_datasec_img()
    _ = multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    _ = multi_caliBrate.get_slits()
    _ = multi_caliBrate.get_arc()
    _ = multi_caliBrate.get_wv_calib()
    _ = multi_caliBrate.get_tilts()
    _ = multi_caliBrate.get_pixflatnrm()
    # Run
    mswave = multi_caliBrate.get_wave()
    assert mswave.shape == (2048,350)
