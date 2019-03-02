"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit import calibrations
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

from pypeit.tests.tstutils import dev_suite_required, dummy_fitstbl

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
def fitstbl():
    if os.getenv('PYPEIT_DEV') is None:
        fitstbl = dummy_fitstbl(directory=data_path(''))
        fitstbl['filename'][1] = 'b1.fits.gz'
        fitstbl['filename'][5] = 'b27.fits.gz'
        return fitstbl

    fitstbl = dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                           'Shane_Kast_blue', '600_4310_d55'))
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'

    return fitstbl


@pytest.fixture
def multi_caliBrate(fitstbl):
    det = 1
    spectrograph = load_spectrograph('shane_kast_blue')
    # Grab a science file for configuration specific parameters
    for idx, row in enumerate(fitstbl):
        if 'science' in row['frametype']:
            sci_file = os.path.join(row['directory'], row['filename'])
            break
    # Par
    par = spectrograph.config_specific_par(sci_file)
    #
    calib_par = par['calibrations']
    calib_par['badpix'] = False
    calib_par['biasframe']['useframe'] = 'overscan'

    multi_caliBrate= calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph)
    # Find the first science row
    frame = fitstbl.find_frames('science', index=True)[0]
    # Set
    multi_caliBrate.set_config(frame, det, par=calib_par)
    return multi_caliBrate


def test_instantiate(fitstbl):
    par = pypeitpar.PypeItPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl, par['calibrations'], spectrograph)


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
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    # Run
    tslits_dict, maskslits = multi_caliBrate.get_slits(write_qa=False)
    # Test
    assert isinstance(tslits_dict, dict)
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_wv_calib(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_arc()
    multi_caliBrate.get_slits(write_qa=False)
    # Run
    wv_calib, maskslits = multi_caliBrate.get_wv_calib()
    assert isinstance(wv_calib, dict)
    assert wv_calib['0'] is not None
    assert wv_calib['0']['rms'] < 0.2
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_tilts(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_arc()
    multi_caliBrate.get_slits(write_qa=False)
    multi_caliBrate.get_wv_calib()
    # Run
    tilts_dict, maskslits = multi_caliBrate.get_tilts()
    assert tilts_dict['tilts'].shape == (2048,350)
    assert isinstance(maskslits, np.ndarray)


@dev_suite_required
def test_flat(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_arc()
    multi_caliBrate.get_slits(write_qa=False)
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
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.get_arc()
    multi_caliBrate.get_slits(write_qa=False)
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    multi_caliBrate.get_flats()
    # Run
    mswave = multi_caliBrate.get_wave()
    assert mswave.shape == (2048,350)


