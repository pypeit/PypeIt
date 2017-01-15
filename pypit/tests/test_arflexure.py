# Module to run tests on simple fitting routines for arrays

### TEST_UNICODE_LITERALS

import numpy as np
import sys
import os, pdb
import pytest

from astropy import units as u

from linetools.spectra import io as lsio

from pypit import arutils as arut
from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arparse as settings  # Has to come after the logger
import pypit
from pypit import arwave

msgs._debug['testing'] = True

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_flex_shift():
    # Dummy slf
    arut.dummy_settings()
    settings.argflag['reduce']['flexure']['maxshift'] = 50
    slf = arut.dummy_self()
    # Read spectra
    obj_spec = lsio.readspec(data_path('obj_lrisb_600_sky.fits'))
    arx_file = pypit.__path__[0]+'/data/sky_spec/sky_LRISb_600.fits'
    arx_spec = lsio.readspec(arx_file)
    # Call
    #msgs._debug['flexure'] = True
    flex_dict = arwave.flex_shift(slf, 1, obj_spec, arx_spec)
    assert np.abs(flex_dict['shift'] - 43.7) < 0.1


