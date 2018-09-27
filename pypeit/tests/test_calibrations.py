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

from pypeit import metadata
from pypeit.core import fsort

from pypeit.spectrographs.util import load_spectrograph
from pypeit import calibrations
from pypeit.par import pypeitpar

from pypeit.tests import tstutils
from pypeit.tests.tstutils import dev_suite_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def fitstbl():
    if os.getenv('PYPEIT_DEV') is None:
        fitstbl = metadata.dummy_fitstbl(directory=data_path(''))
        fitstbl['filename'][1] = 'b1.fits.gz'
        return fitstbl

    fitstbl = metadata.dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                           'Shane_Kast_blue', '600_4310_d55'))
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'

    return fitstbl

#def chk_for_files(root):
#    return len(glob.glob(root+'*')) != 0

@pytest.fixture
def multi_caliBrate(fitstbl):
    setup = 'A_01_aa'
    det = 1
    sci_ID = 1
    spectrograph = tstutils.load_kast_blue_masters(get_spectrograph=True)[0]

    # Only changing the defaults
    calib_par = pypeitpar.CalibrationsPar(badpix=False,
                                biasframe=pypeitpar.FrameGroupPar('bias', useframe='overscan'))

    redux_path = data_path('') if os.getenv('PYPEIT_DEV') is None \
                                else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked')

    multi_caliBrate= calibrations.MultiSlitCalibrations(fitstbl, spectrograph=spectrograph,
                                                        par=calib_par, redux_path=redux_path,
                                                        save_masters=False, write_qa=False)
    multi_caliBrate.reset(setup, det, sci_ID, calib_par)
    return multi_caliBrate


def test_instantiate(fitstbl):
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)
    print(caliBrate)

#@dev_suite_required
#def test_datasec(multi_caliBrate):
#    multi_caliBrate.par['trim'] = False
#    datasec_img = multi_caliBrate.get_datasec_img()
#    naxis0, naxis1 = datasec_img.shape
#    # Test
#    assert naxis0 == 2112
#    assert naxis1 == 350

def test_pixlocn(multi_caliBrate):
    multi_caliBrate.shape = (2048,350)
    pixlocn = multi_caliBrate.get_pixlocn()
    # Test
    assert pixlocn.shape == (2048,350,4)


def test_bias(multi_caliBrate):
    multi_caliBrate.get_bias()


def test_arc(multi_caliBrate):
    multi_caliBrate.msbias = 'overscan'
    arc = multi_caliBrate.get_arc()
    assert arc.shape == (2048,350)


def test_bpm(multi_caliBrate):
    # Prep
    multi_caliBrate.shape = (2048,350)
    # Build
    bpm = multi_caliBrate.get_bpm()
    assert bpm.shape == (2048,350)
    assert np.sum(bpm) == 0.


@dev_suite_required
def test_slits(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_pixlocn()
#    multi_caliBrate.get_datasec_img()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    # Run
    tslits_dict, maskslits = multi_caliBrate.get_slits()
    # Test
    assert isinstance(tslits_dict, dict)
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_wv_calib(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_pixlocn()
#    multi_caliBrate.get_datasec_img()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_slits()
    multi_caliBrate.get_arc()
    # Run
    wv_calib, maskslits = multi_caliBrate.get_wv_calib()
    assert isinstance(wv_calib, dict)
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_tilts(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_pixlocn()
#    multi_caliBrate.get_datasec_img()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_slits()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_wv_calib()
    # Run
    tilts_dict, maskslits = multi_caliBrate.get_tilts()
    # TODO: Include a test of tilts_dict!
#    assert mstilts.shape == (2048,350)
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_flat(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_pixlocn()
#    multi_caliBrate.get_datasec_img()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_slits()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    # Run
    mspixflatnrm, msillumflat = multi_caliBrate.get_flats()
    assert mspixflatnrm.shape == (2048,350)
    assert msillumflat.shape == (2048,350)


@dev_suite_required
def test_waveimg(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    multi_caliBrate.get_pixlocn()
#    multi_caliBrate.get_datasec_img()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_slits()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    multi_caliBrate.get_flats()
    # Run
    mswave = multi_caliBrate.get_wave()
    assert mswave.shape == (2048,350)


