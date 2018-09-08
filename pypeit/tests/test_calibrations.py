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

from pypeit.core import fsort

from pypeit.spectrographs.util import load_spectrograph
from pypeit import calibrations
from pypeit.par import pypeitpar

from pypeit.tests import tstutils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
    fitstbl = fsort.dummy_fitstbl(directory=data_path(''))
    fitstbl['filename'][1] = 'b1.fits.gz'
else:
    skip_test=False
    # MultiSlit
    fitstbl = fsort.dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                          'Shane_Kast_blue', '600_4310_d55'))
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
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
    spectrograph = tstutils.load_kast_blue_masters(get_spectrograph=True)[0]

    # Only changing the defaults
    # TODO: (KBW) I agree FrameGroupPar is a bit onerous...
    calib_par = pypeitpar.CalibrationsPar(badpix=False,
                                         biasframe=pypeitpar.FrameGroupPar('bias',
                                                                          useframe='overscan'))

    redux_path = data_path('') if os.getenv('PYPEIT_DEV') is None \
        else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked')
    #master_root = data_path('MF') if os.getenv('PYPEIT_DEV') is None \
    #                else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF')
    
    multi_caliBrate= calibrations.MultiSlitCalibrations(fitstbl, spectrograph=spectrograph,
                                                        par=calib_par, redux_path=redux_path,
                                                        save_masters=False, write_qa=False)
    multi_caliBrate.reset(setup, det, sci_ID, calib_par)
    return multi_caliBrate


def test_instantiate():
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)
    print(caliBrate)


def test_datasec(multi_caliBrate):
    if skip_test:
        assert True
        return
    multi_caliBrate.par['trim'] = False
    datasec_img = multi_caliBrate.get_datasec_img()
    naxis0, naxis1 = datasec_img.shape
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
    if skip_test:
        assert True
        return
    # Setup
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
    _ = multi_caliBrate.get_tilts()
    _ = multi_caliBrate.get_pixflatnrm()
    # Run
    mswave = multi_caliBrate.get_wave()
    assert mswave.shape == (2048,350)

