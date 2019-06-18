"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import processrawimage
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


@dev_suite_required
def test_instantiate(deimos_flat_files, kast_blue_bias_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flat = processrawimage.ProcessRawImage(one_file, spectograph, 3, par)
    # Test
    assert isinstance(deimos_flat.image, np.ndarray)
    assert deimos_flat.rawdatasec_img.shape == (4096, 2128)

    # Kast blue
    one_file = kast_blue_bias_files[0]
    spectograph2 = load_spectrograph('shane_kast_blue')
    kastb_bias = processrawimage.ProcessRawImage(one_file, spectograph2, 1, par)


@dev_suite_required
def test_overscan_subtract(deimos_flat_files):
    one_file = deimos_flat_files[0]
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flat = processrawimage.ProcessRawImage(one_file, spectograph, 3, par)
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


