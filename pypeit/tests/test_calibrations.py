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

from pypeit.tests.tstutils import dummy_fitstbl

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


#@pytest.fixture
def fitstbl():
    if os.getenv('PYPEIT_DEV') is None:
        fitstbl = dummy_fitstbl(directory=data_path(''))
        fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
        fitstbl['filename'][1] = 'b1.fits.gz'
        fitstbl['filename'][5] = 'b27.fits.gz'
        return fitstbl

    fitstbl = dummy_fitstbl(directory=os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                   'shane_kast_blue', '600_4310_d55'))
    # Set the Bias to known
    fitstbl['framebit'][0] = fitstbl.type_bitmask.turn_off(fitstbl['framebit'][0], flag='bias')
    fitstbl['filename'][1] = 'b1.fits.gz'
    for ii in range(2,5):
        fitstbl['filename'][ii] = 'b{0}.fits.gz'.format(ii)
    fitstbl['filename'][5] = 'b27.fits.gz'

    return fitstbl


#@pytest.fixture
def multi_caliBrate(fitstbl):
    # Grab a science file for configuration specific parameters
    for idx, row in enumerate(fitstbl):
        if 'science' in row['frametype']:
            sci_file = os.path.join(row['directory'], row['filename'])
            break
    # Par
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.config_specific_par(sci_file)
    turn_off = dict(use_biasimage=False)
    par.reset_all_processimages_par(**turn_off)
    #
    calib_par = par['calibrations']
    calib_par['bpm_usebias'] = False
    #calib_par['biasframe']['useframe'] = 'none' # Only use overscan
    calib_par['slitedges']['sync_predict'] = 'nearest'

    caldir = data_path('Calibrations')

    multi_caliBrate = calibrations.MultiSlitCalibrations(fitstbl, calib_par, spectrograph, caldir)
    return reset_calib(multi_caliBrate)


def reset_calib(calib):
    # Find the first science row
    frame = calib.fitstbl.find_frames('science', index=True)[0]
    # Set
    det = 1
    calib.set_config(frame, det)
    return calib


###################################################
# TESTS BEGIN HERE

def test_abstract_init(fitstbl):
    par = pypeitpar.CalibrationsPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caldir = data_path('Calibrations')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.MultiSlitCalibrations), 'Wrong calibration object type'
    spectrograph = load_spectrograph('keck_nires')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.MultiSlitCalibrations), 'Wrong calibration object type'
    spectrograph = load_spectrograph('keck_kcwi')
    calib = calibrations.Calibrations.get_instance(fitstbl, par, spectrograph, caldir)
    assert isinstance(calib, calibrations.IFUCalibrations), 'Wrong calibration object type'


def test_instantiate(fitstbl):
    par = pypeitpar.CalibrationsPar()
    spectrograph = load_spectrograph('shane_kast_blue')
    caldir = data_path('Calibrations')
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl, par, spectrograph, caldir)


def test_bias(multi_caliBrate):
    """
    This should produce nothing as we have no bias frames

    Returns:

    """
    multi_caliBrate.get_bias()

test_bias(multi_caliBrate(fitstbl()))


def test_arc(multi_caliBrate):
    arc = multi_caliBrate.get_arc()
    assert arc.image.shape == (2048,350)

    # Cleanup
    shutil.rmtree(multi_caliBrate.master_dir)


def test_tiltimg(multi_caliBrate):
    tilt = multi_caliBrate.get_tiltimg()
    assert tilt.image.shape == (2048,350)

    # Cleanup
    shutil.rmtree(multi_caliBrate.master_dir)

def test_bpm(multi_caliBrate):
    # Prep
    multi_caliBrate.shape = (2048,350)
    # Build
    bpm = multi_caliBrate.get_bpm()
    assert bpm.shape == (2048,350)
    assert np.sum(bpm) == 0.
