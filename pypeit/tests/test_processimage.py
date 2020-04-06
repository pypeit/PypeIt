"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import rawimage
from pypeit.images import pypeitimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

par = pypeitpar.ProcessImagesPar()

kast_blue = load_spectrograph('shane_kast_blue')

@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                      '830G_L_8400', ifile) 
                            for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue',
                                  '600_4310_d55', 'b1?.fits*'))


@dev_suite_required
def test_instantiate(deimos_flat_files, kast_blue_bias_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    det = 3
    rawImage = rawimage.RawImage(one_file, spectograph, det)
    # Test
    assert isinstance(rawImage.image, np.ndarray)
    assert rawImage.datasec_img.shape == (4096, 2128)

    # Kast blue
    det2 = 1
    one_file = kast_blue_bias_files[0]
    spectograph2 = load_spectrograph('shane_kast_blue')
    rawImage2 = rawimage.RawImage(one_file, spectograph2, det2)
    assert isinstance(rawImage2.image, np.ndarray)


@dev_suite_required
def test_overscan_subtract(deimos_flat_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    det = 3
    rawImage = rawimage.RawImage(one_file, spectograph, det)
    rawImage.par = spectograph.default_pypeit_par()['scienceframe']['process']
    # Bias subtract
    pre_sub = rawImage.image.copy()
    rawImage.subtract_overscan()
    oscan = np.median(pre_sub-rawImage.image)
    assert np.isclose(oscan, 1001.2, rtol=0.01)
    # Trim
    rawImage.trim()
    # Test
    assert rawImage.steps['subtract_overscan']
    assert rawImage.steps['trim']
    assert rawImage.image.shape == (4096,2048)


