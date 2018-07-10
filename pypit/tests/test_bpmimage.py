# Module to run tests on BPMImage class
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

from pypit import bpmimage
from pypit.armsgs import PypitError

from pypit.spectrographs.util import load_spectrograph

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_instantiate():
    # Empty
    # Can no longer be empty!
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage()
    bpm = bpmimage.BPMImage(spectrograph='keck_lris_red')
    assert bpm.spectrograph.spectrograph == 'keck_lris_red'

    # These will no longer error!
#    with pytest.raises(PypitError):
#        _ = bpmimage.BPMImage('keck_lris_red', det=1)
#    with pytest.raises(PypitError):
#        _ = bpmimage.BPMImage('keck_lris_red')
#    with pytest.raises(PypitError):
#        _ = bpmimage.BPMImage('keck_lris_red', binning='1,1')
#    with pytest.raises(PypitError):
#        _ = bpmimage.BPMImage('keck_deimos')
#   TODO: reduce_badpix is no longer an argument.  This functionality is
#   performed when calling `build` when msbias is not None
#    with pytest.raises(PypitError):
#        _ = bpmimage.BPMImage(reduce_badpix='bias')

def test_dummy_image():
    # Simple
    shape=(2048,2048)
    bpmImage = bpmimage.BPMImage(spectrograph='shane_kast_blue', shape=shape)
    bpm = bpmImage.build()
    assert isinstance(bpm, np.ndarray)
    assert bpm.shape == shape
    assert np.sum(bpm) == 0

def test_keck_lris_red():
    # TODO: Now requires a file:
    if skip_test:
        return
    example_file = os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits')
    # Simple
    bpmImage = bpmimage.BPMImage(spectrograph='keck_lris_red', filename=example_file, det=2)
    bpm = bpmImage.build()
    assert np.sum(bpm) > 0

def test_keck_deimos():
    if skip_test:
        return
    example_file = os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L',
                                'd0914_0002.fits')
    # Simple
    bpmImage = bpmimage.BPMImage(spectrograph='keck_deimos', filename=example_file, det=4)
    bpm = bpmImage.build()
    assert bpm[0,0] == 1

def test_bpm_from_bias():
    bias = np.full((1024,1024), 1000, dtype=float)
    bias[512,512] += 50.
    bpmImage = bpmimage.BPMImage(msbias=bias)
    bpm = bpmImage.build()
    # Test
    assert np.isclose(bpm[512,512],1)


