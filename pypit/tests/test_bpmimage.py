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

from pypit.spectrographs import bpmimage
from pypit.armsgs import PypitError

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
    bpmImage = bpmimage.BPMImage()
    # Errors
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage(spectrograph='keck_lris_red', det=1)
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage(spectrograph='keck_lris_red')
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage(spectrograph='keck_lris_red', binning='1,1')
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage(spectrograph='keck_deimos')
    with pytest.raises(PypitError):
        _ = bpmimage.BPMImage(reduce_badpix='bias')


def test_dummy_image():
    # Simple
    shape=(2048,2048)
    bpmImage = bpmimage.BPMImage(shape=shape)
    bpm = bpmImage.build()
    assert isinstance(bpm, np.ndarray)
    assert bpm.shape == shape
    assert np.sum(bpm) == 0.


def test_keck_lris_red():
    # Simple
    bpmImage = bpmimage.BPMImage(spectrograph='keck_lris_red', binning='1,1', det=2)
    bpm = bpmImage.build()
    assert np.sum(bpm) > 0.

def test_keck_deimos():
    # Simple
    bpmImage = bpmimage.BPMImage(spectrograph='keck_deimos', det=4)
    bpm = bpmImage.build()
    assert bpm[0,0] == 1.

def test_bpm_from_bias():
    bias = np.ones((1024,1024))*1000
    bias[512,512] += 50.
    settings = dict(detector={})
    settings['detector']['numamplifiers'] = 1
    settings['detector']['datasec01'] = [[0,1024], [0,1024]]
    #
    bpmImage = bpmimage.BPMImage(reduce_badpix='bias', msbias=bias, settings=settings)
    bpm = bpmImage.build()
    # Test
    assert np.isclose(bpm[512,512],1.)
