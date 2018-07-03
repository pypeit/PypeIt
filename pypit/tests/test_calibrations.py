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
fitstbl = arsort.dummy_fitstbl()
setup = 'A_01_aa'
det = 1
sci_ID = 1
settings = tstutils.load_kast_blue_masters(get_settings=True)[0]

def test_instantiate():
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)

def test_bias():
    caliBrate = calibrations.MultiSlitCalibrations(fitstbl)
    #
    settings['bias'] = {}
    settings['bias']['useframe'] = 'overscan'
    caliBrate.reset(setup, det, sci_ID, settings)
    # Build
    caliBrate.get_bias()

