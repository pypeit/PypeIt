# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import sys
import os

import pdb
import pytest

import numpy as np

from linetools.spectra.io import readspec

import pypeit
#from pypeit import msgs
from pypeit.core import wave


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_flex_shift():
    # Dummy slf
    settings.dummy_settings()
    settings.argflag['reduce']['flexure']['maxshift'] = 50
    # Read spectra
    obj_spec = readspec(data_path('obj_lrisb_600_sky.fits'))
    arx_file = pypeit.__path__[0]+'/data/sky_spec/sky_LRISb_600.fits'
    arx_spec = readspec(arx_file)
    # Call
    #msgs._debug['flexure'] = True
    flex_dict = wave.flex_shift(obj_spec, arx_spec)
    assert np.abs(flex_dict['shift'] - 15.7) < 0.1


