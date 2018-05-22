# Module to run tests on TraceImage class
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

from pypit import traceimage

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def deimos_flat_files():
    if not skip_test:
        deimos_flat_files = [os.getenv('PYPIT_DEV') + '/RAW_DATA/Keck_DEIMOS/830G_L/' + ifile for ifile in [  # Longslit in dets 3,7
            'd0914_0014.fits', 'd0914_0015.fits']]
    else:
        deimos_flat_files = None
    return deimos_flat_files


def test_instantiate(deimos_flat_files):
    if skip_test:
        assert True
        return
    # Empty
    Timage = traceimage.TraceImage(deimos_flat_files)
    assert Timage.nfiles == 2


def test_process(deimos_flat_files):
    if skip_test:
        assert True
        return
    # Instantiate
    Timage = traceimage.TraceImage(deimos_flat_files, spectrograph='keck_deimos')
    # Run
    mstrace = Timage.process(bias_subtract='overscan', trim=True)
    assert isinstance(mstrace, np.ndarray)
    assert Timage.steps[-1] == 'combine'
    assert Timage.steps[-2] == 'bias_subtract'

