"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import shutil

import numpy as np

from pypeit import calibrations
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from IPython import embed

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
    calib_par['biasframe']['useframe'] = 'none' # Only use overscan
    calib_par['slitedges']['sync_predict'] = 'nearest'

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph)
    return reset_calib(multi_caliBrate)


def reset_calib(calib):
    # Find the first science row
    frame = calib.fitstbl.find_frames('science', index=True)[0]
    # Set
    det = 1
    calib.set_config(frame, det)
    return calib


@pytest.fixture
def multi_caliBrate_reuse(multi_caliBrate):
    multi_caliBrate.reuse_masters=True
    multi_caliBrate.master_dir = data_path('Masters')
    multi_caliBrate.save_masters = True
    return multi_caliBrate

###################################################
# TESTS BEGIN HERE


def test_instantiate(fitstbl):
    par = pypeitpar.PypeItPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl, par['calibrations'], spectrograph)


def test_bias(multi_caliBrate):
    """
    This should only overscan subtract

    Returns:

    """
    multi_caliBrate.get_bias()


def test_arc(multi_caliBrate):
    arc = multi_caliBrate.get_arc()
    assert arc.image.shape == (2048,350)


def test_tiltimg(multi_caliBrate):
    tilt = multi_caliBrate.get_tiltimg()
    assert tilt.image.shape == (2048,350)

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
    # Run
    tslits_dict = multi_caliBrate.get_slits(write_qa=False)
    # Test
    assert isinstance(tslits_dict, dict)
    assert isinstance(tslits_dict['maskslits'], np.ndarray)


@dev_suite_required
def test_wv_calib(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_slits(write_qa=False)
    # Run
    wv_calib = multi_caliBrate.get_wv_calib()
    assert isinstance(wv_calib, dict)
    assert wv_calib['0'] is not None
    assert wv_calib['0']['rms'] < 0.2
    assert isinstance(multi_caliBrate.tslits_dict['maskslits'], np.ndarray)


@dev_suite_required
def test_tilts(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_tiltimg()
    multi_caliBrate.get_slits(write_qa=False)
    multi_caliBrate.get_wv_calib()
    # Run
    tilts_dict = multi_caliBrate.get_tilts()
    assert tilts_dict['tilts'].shape == (2048,350)
    assert isinstance(multi_caliBrate.tslits_dict['maskslits'], np.ndarray)


@dev_suite_required
def test_flat(multi_caliBrate):
    # Setup
    multi_caliBrate.shape = (2048,350)
    #multi_caliBrate.get_pixlocn()
    multi_caliBrate.get_bpm()
    multi_caliBrate.get_arc()
    multi_caliBrate.get_tiltimg()
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
    multi_caliBrate.get_arc()
    multi_caliBrate.get_tiltimg()
    multi_caliBrate.get_slits(write_qa=False)
    multi_caliBrate.get_wv_calib()
    multi_caliBrate.get_tilts()
    multi_caliBrate.get_flats()
    # Run
    mswave = multi_caliBrate.get_wave()
    assert mswave.shape == (2048,350)


@dev_suite_required
def test_reuse(multi_caliBrate_reuse):
    """
    Test that Calibrations appropriately reuses existing calibrations frames.
    """
    # In case of previous data or failures
    if os.path.isdir(multi_caliBrate_reuse.master_dir):
        shutil.rmtree(multi_caliBrate_reuse.master_dir)

    os.makedirs(multi_caliBrate_reuse.master_dir)

    # Perform the calibrations and check that the data are correctly
    # stored in memory
    multi_caliBrate_reuse.shape = (2048,350)
    multi_caliBrate_reuse.get_bpm()
    assert list(multi_caliBrate_reuse.calib_dict.keys()) == ['A_1_01'], 'Incorrect master key'
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) == ['bpm'], \
                'Incorrect list of master types in memory'
    msarc = multi_caliBrate_reuse.get_arc()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) == ['bpm', 'arc'], \
                'Incorrect list of master types in memory'
    msarc = multi_caliBrate_reuse.get_tiltimg()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) == ['bpm', 'arc', 'tiltimg'], \
                'Incorrect list of master types in memory'
    multi_caliBrate_reuse.get_slits(write_qa=False)
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) == ['bpm', 'arc', 'tiltimg', 'trace'], \
                'Incorrect list of master types in memory'
    multi_caliBrate_reuse.get_wv_calib()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) \
                == ['bpm', 'arc', 'tiltimg','trace', 'wavecalib', 'wvmask'], \
                'Incorrect list of master types in memory'
    multi_caliBrate_reuse.get_tilts()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) \
                == ['bpm', 'arc', 'tiltimg', 'trace', 'wavecalib', 'wvmask', 'tilts_dict', 'wtmask'], \
                'Incorrect list of master types in memory'
    multi_caliBrate_reuse.get_flats()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) \
                == ['bpm', 'arc', 'tiltimg', 'trace', 'wavecalib', 'wvmask', 'tilts_dict', 'wtmask',
                    'pixelflat', 'illumflat'], \
                'Incorrect list of master types in memory'
    mswave = multi_caliBrate_reuse.get_wave()
    assert list(multi_caliBrate_reuse.calib_dict['A_1_01'].keys()) \
                == ['bpm', 'arc', 'tiltimg', 'trace', 'wavecalib', 'wvmask', 'tilts_dict', 'wtmask',
                    'pixelflat', 'illumflat', 'wave'], \
                'Incorrect list of master types in memory'
    assert mswave.shape == (2048,350)

    # Reset
    reset_calib(multi_caliBrate_reuse)
    multi_caliBrate_reuse.save_masters = False  # Don't overwrite them

    # Redo the calibrations
    #   - These don't source a master file
    multi_caliBrate_reuse.shape = (2048,350)
    multi_caliBrate_reuse.get_bpm()
    #   - The arc is the first sourced master
    assert 'arc' not in multi_caliBrate_reuse.master_key_dict.keys(), \
            'arc master key should not be defined yet'
    _msarc = multi_caliBrate_reuse.get_arc()
    assert multi_caliBrate_reuse._cached('arc',
                    multi_caliBrate_reuse.master_key_dict['arc']), 'Should find cached data.'
    _msarc = multi_caliBrate_reuse.get_tiltimg()
    assert multi_caliBrate_reuse._cached('tiltimg',
                    multi_caliBrate_reuse.master_key_dict['arc']), 'Should find cached data.'
    assert os.path.isfile(multi_caliBrate_reuse.arcImage.master_file_path), \
            'Should find master file.'
    assert multi_caliBrate_reuse.arcImage.load() is not None, \
            'Load should not return None'
    # TODO: Not a great test because this should be true regardless of
    # whether or not the master was actually reused...
    assert np.array_equal(msarc, _msarc), 'Arrays not equal!'
    #   - Make sure the rest of the steps complete
    multi_caliBrate_reuse.get_slits(write_qa=False)
    multi_caliBrate_reuse.get_wv_calib()
    multi_caliBrate_reuse.get_tilts()
    multi_caliBrate_reuse.get_flats()
    mswave = multi_caliBrate_reuse.get_wave()

    # Clean-up
    shutil.rmtree(multi_caliBrate_reuse.master_dir)

