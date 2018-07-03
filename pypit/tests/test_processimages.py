# Module to run tests on ProcessImages class
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

from pypit import processimages

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

@pytest.fixture
def deimos_flat_files():
    if not skip_test:
        # Longslit in dets 3,7
        deimos_flat_files = [os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA/Keck_DEIMOS/830G_L',
                                        ifile) for ifile in ['d0914_0014.fits', 'd0914_0015.fits']]
        assert len(deimos_flat_files) == 2
    else:
        deimos_flat_files = None
    return deimos_flat_files

@pytest.fixture
def kast_blue_bias_files():
    if not skip_test:
        kast_blue_bias_files = glob.glob(os.path.join(os.getenv('PYPIT_DEV'),
                                            'RAW_DATA/Shane_Kast_blue/600_4310_d55/b1?.fits*'))
    else:
        kast_blue_bias_files = None
    return kast_blue_bias_files


@pytest.fixture
def kast_settings():
    # Instantiate
    kast_settings = processimages.default_settings()
    kast_settings['detector']['dataext'] = 0
    kast_settings['detector']['dataext01'] = 0
    kast_settings['detector']['datasec01'] = [[0, 1024], [0, 0]]
    kast_settings['detector']['datasec02'] = [[1024, 2048], [0, 0]]
    kast_settings['detector']['oscansec01'] = [[2049, 2080], [0, 0]]
    kast_settings['detector']['oscansec02'] = [[2080, 2111], [0, 0]]
    return kast_settings


def test_instantiate():
    proc_img = processimages.ProcessImages([])
    assert proc_img.nfiles == 0


def test_load(deimos_flat_files, kast_blue_bias_files, kast_settings):
    if skip_test:
        assert True
        return
    # DEIMOS
    deimos_flats = processimages.ProcessImages(deimos_flat_files, spectrograph='keck_deimos')
    # Load
    deimos_flats.load_images()
    # Test
    assert deimos_flats.nloaded == 2
    assert deimos_flats.steps == ['load_images']

    # Kast blue
    kastb_bias = processimages.ProcessImages(kast_blue_bias_files, spectrograph='shane_kast_blue', settings=kast_settings)
    # Load
    kastb_bias.load_images()
    # Check datasec
    assert kastb_bias.datasec[0][0] == [0,1024]


def test_bias_subtract(deimos_flat_files):
    if skip_test:
        assert True
        return
    # DEIMOS
    deimos_flats = processimages.ProcessImages(deimos_flat_files, spectrograph='keck_deimos')
    # Load
    deimos_flats.load_images()
    # Bias subtract (and trim)
    deimos_flats.bias_subtract('overscan')
    # Test
    assert isinstance(deimos_flats.proc_images, np.ndarray)
    assert deimos_flats.proc_images.shape == (4096,2048,2)


def test_combine(deimos_flat_files):
    if skip_test:
        assert True
        return
    # DEIMOS
    deimos_flats = processimages.ProcessImages(deimos_flat_files, spectrograph='keck_deimos')
    # Load
    deimos_flats.load_images()
    # Bias subtracgt
    deimos_flats.bias_subtract('overscan')
    # Combine
    stack = deimos_flats.combine()
    # Test
    assert isinstance(deimos_flats.stack, np.ndarray)
    assert deimos_flats.stack.shape == (4096,2048)
