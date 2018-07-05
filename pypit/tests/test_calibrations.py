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
    multi_caliBrate= calibrations.MultiSlitCalibrations(fitstbl)
    #
    multi_caliBrate.reset(setup, det, sci_ID, settings)
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
    multi_caliBrate.settings['bias'] = {}
    multi_caliBrate.settings['bias']['useframe'] = 'overscan'
    # Build
    multi_caliBrate.get_bias()


def test_arc(multi_caliBrate):
    multi_caliBrate.msbias = 'overscan'
    multi_caliBrate.settings['masters']['reuse'] = False
    # Build
    arc = multi_caliBrate.get_arc()
    assert arc.shape == (2048,350)


def test_bpm(multi_caliBrate):
    # Prep
    multi_caliBrate.settings['reduce'] = {}
    multi_caliBrate.settings['reduce']['badpix'] = False
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
    pixlocn = multi_caliBrate.get_pixlocn()
    datasec_img, naxis0, naxis1 = multi_caliBrate.get_datasec_img()
    multi_caliBrate.settings['reduce'] = {}
    multi_caliBrate.settings['reduce']['badpix'] = False
    _ = multi_caliBrate.get_bpm()
    # Settings -- To be replaced
    from pypit import traceslits
    multi_caliBrate.settings['trace'] = traceslits.default_settings()['trace']
    # Run
    tslits_dict, maskslits = multi_caliBrate.get_slits()
    # Test
    assert isinstance(tslits_dict, dict)
    assert isinstance(maskslits, np.ndarray)
