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

from pypeit import traceimage

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def deimos_flat_files():
    if not skip_test:
        deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS',
                                          '830G_L_8400', ifile)
                                for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    else:
        deimos_flat_files = None
    return deimos_flat_files


def test_instantiate(deimos_flat_files):
    if skip_test:
        assert True
        return
    # Empty
    traceImage = traceimage.TraceImage('keck_deimos')


def test_process(deimos_flat_files):
    if skip_test:
        assert True
        return
    # Instantiate
    traceImage = traceimage.TraceImage('keck_deimos', file_list=deimos_flat_files)
    # Run
    assert traceImage.nfiles == 2
    mstrace = traceImage.process(bias_subtract='overscan', trim=True)
    assert isinstance(mstrace, np.ndarray)
    assert traceImage.steps[-1] == 'combine'
    assert traceImage.steps[-2] == 'bias_subtract'

