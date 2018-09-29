# Module to run tests on arcoadd
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import os

import pytest
import numpy as np

from astropy import units
from linetools.spectra.utils import collate
from linetools.spectra.xspectrum1d import XSpectrum1D

from pypeit.core import coadd
from pypeit import msgs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def test_cooked_version():
    if skip_test:
        assert True
        return
    # Load up the version
    v_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'version')
    with open(v_file) as f:
        tmp = f.readlines()
    value = float(tmp[-1].strip())
    # Test
    assert value >= 0.82

