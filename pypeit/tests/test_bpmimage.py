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

from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs import util
from pypeit.core import procimg


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)



def test_dummy_image():
    # Simple
    shape=(2048,2048)
    spectrograph = util.load_spectrograph('shane_kast_blue')
    bpm = spectrograph.bpm(shape=shape)#, trim=False)
    assert isinstance(bpm, np.ndarray)
    assert bpm.shape == shape
    assert np.sum(bpm) == 0


@dev_suite_required
def test_keck_lris_red():
    # Spectrograph
    spectrograph = util.load_spectrograph('keck_lris_red')
    #
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    # Get the shape
    dsec_img = spectrograph.get_datasec_img(example_file, det=2)
    shape = procimg.trim_frame(dsec_img, dsec_img < 1).shape
    # Simple
    bpm = spectrograph.bpm(shape=shape, filename=example_file, det=2)
    assert np.sum(bpm) > 0


@dev_suite_required
def test_keck_deimos():
    spectrograph = util.load_spectrograph('keck_deimos')
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L_8400',
                                'd0914_0002.fits.gz')
    # Get the shape
    dsec_img = spectrograph.get_datasec_img(example_file, det=2)
    shape = procimg.trim_frame(dsec_img, dsec_img < 1).shape
    # Simple
    bpm = spectrograph.bpm(shape=shape,det=4)
    assert bpm[0,0] == 1


# This is too experimental
'''
def test_bpm_from_bias():
    bias = np.full((1024,1024), 1000, dtype=float)
    bias[512,512] += 50.
    bpmImage = bpmimage.BPMImage(msbias=bias)
    bpm = bpmImage.build()
    # Test
    assert np.isclose(bpm[512,512],1)
'''


