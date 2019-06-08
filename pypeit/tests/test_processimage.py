"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import processimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

par = pypeitpar.ProcessImagesPar()

kast_blue = load_spectrograph('shane_kast_blue')

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
    processImage = processimage.ProcessImage(None, kast_blue, 1, par)
    for key in processImage.steps.keys():
        assert processImage.steps[key] is False


@dev_suite_required
def test_load(deimos_flat_files, kast_blue_bias_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flat = processimage.ProcessImage(one_file, spectograph, 3, par, frametype='pixelflat')
    # Load
    img, head0 = deimos_flat.load_rawimage()
    # Test
    assert isinstance(img, np.ndarray)
    assert deimos_flat.rawdatasec_img.shape == (4096, 2128)

    # Kast blue
    one_file = kast_blue_bias_files[0]
    spectograph2 = load_spectrograph('shane_kast_blue')
    kastb_bias = processimage.ProcessImage(one_file, spectograph2, 1, par)
    # Load
    kastb_bias.load_rawimage()
    # Check datasec
    assert kastb_bias.datasec_img.shape == (350,2112)  # Not oriented yet


@dev_suite_required
def test_overscan_subtract(deimos_flat_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flat = processimage.ProcessImage(one_file, spectograph, 3, par, frametype='pixelflat')
    # Load
    deimos_flat.load_rawimage()
    # Bias subtract
    pre_sub = deimos_flat.image.copy()
    _ = deimos_flat.subtract_overscan()
    oscan = np.median(pre_sub-deimos_flat.image)
    assert np.isclose(oscan, 1001.2, rtol=0.01)
    # Trim
    _ = deimos_flat.trim()
    # Test
    assert deimos_flat.steps['subtract_overscan']
    assert deimos_flat.steps['trim']
    assert deimos_flat.image.shape == (4096,2048)


