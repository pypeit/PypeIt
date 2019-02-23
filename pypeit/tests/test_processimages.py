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

from pypeit import processimages
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar

par = pypeitpar.ProcessImagesPar()

@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS',
                                      '830G_L_8400', ifile) 
                            for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Shane_Kast_blue',
                                  '600_4310_d55', 'b1?.fits*'))


def test_instantiate():
    proc_img = processimages.ProcessImages('shane_kast_blue', par)
    assert proc_img.nfiles == 0


@dev_suite_required
def test_load(deimos_flat_files, kast_blue_bias_files):
    # DEIMOS
    deimos_flats = processimages.ProcessImages('keck_deimos', par, files=deimos_flat_files)
    # Load
    deimos_flats.load_images()
    # Test
    assert deimos_flats.nloaded == 2
    assert deimos_flats.steps == ['load_images']

    # Kast blue
    kastb_bias = processimages.ProcessImages('shane_kast_blue', par, files=kast_blue_bias_files)
    # Load
    kastb_bias.load_images()
    # Check datasec
    assert kastb_bias.datasec[0][0][0] == slice(0,1024,None)


@dev_suite_required
def test_bias_subtract(deimos_flat_files):
    # DEIMOS
    deimos_flats = processimages.ProcessImages('keck_deimos', par, files=deimos_flat_files)
    # Load
    deimos_flats.load_images()
    # Bias subtract (and trim)
    deimos_flats.bias_subtract('overscan')
    # Test
    assert isinstance(deimos_flats.proc_images, np.ndarray)
    assert deimos_flats.proc_images.shape == (4096,2048,2)


@dev_suite_required
def test_combine(deimos_flat_files):
    # DEIMOS
    deimos_flats = processimages.ProcessImages('keck_deimos', par, files=deimos_flat_files)
    # Load
    deimos_flats.load_images()
    # Bias subtracgt
    deimos_flats.bias_subtract('overscan')
    # Combine
    stack = deimos_flats.combine()
    # Test
    assert isinstance(deimos_flats.stack, np.ndarray)
    assert deimos_flats.stack.shape == (4096,2048)


