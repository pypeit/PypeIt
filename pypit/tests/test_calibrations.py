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

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def chk_for_files(root):
    files = glob.glob(root+'*')
    if len(files) == 0:
        return False
    else:
        return True

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

# MultiSlit
fitstbl = arsort.dummy_fitstbl(directory=data_path(''))
fitstbl['filename'][1] = 'b1.fits.gz'

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
