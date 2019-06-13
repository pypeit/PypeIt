"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import calibrationimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import procimg

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

'''

def test_instantiate():
    combinedImage = combinedimage.CombinedImage(kast_blue, 1, par)
    assert combinedImage.nfiles == 0


@dev_suite_required
def test_load(deimos_flat_files, kast_blue_bias_files):
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flats = combinedimage.CombinedImage(spectograph, 3, par, files=deimos_flat_files)
    assert deimos_flats.nfiles == 2
    assert deimos_flats.nimages == 0
    # Load
    deimos_flats.load_images()
    # Test
    assert deimos_flats.nimages == 2

    # Kast blue
    kastb_bias = combinedimage.CombinedImage(kast_blue, 1, par, files=kast_blue_bias_files)
    # Load
    kastb_bias.load_images()
    # Check datasec
    assert kastb_bias.nimages == 10


@dev_suite_required
def test_process(deimos_flat_files):
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flats = combinedimage.CombinedImage(spectograph, 3, par, files=deimos_flat_files)
    # Load
    deimos_flats.load_images()
    # Process steps
    psteps = procimg.init_process_steps(None, deimos_flats.proc_par)
    assert 'subtract_overscan' in psteps
    psteps += ['trim', 'orient']
    # Bias subtract (and trim)
    deimos_flats.process_images(psteps)
    # Test
    assert deimos_flats.pimages[0].images.shape == (4096,2048)


'''

@dev_suite_required
def test_combine(deimos_flat_files):
    spectograph = load_spectrograph('keck_deimos')
    # DEIMOS
    deimos_flats = calibrationimage.CalibrationImage(spectograph, 3, par, files=deimos_flat_files)
    # Process steps
    psteps = procimg.init_process_steps(None, deimos_flats.proc_par)
    assert 'subtract_overscan' in psteps
    psteps += ['trim', 'orient', 'apply_gain']
    deimos_flats.process_steps = psteps
    # Bias subtract (and trim)
    deimos_flats.build_image()
    # Test
    assert deimos_flats.image.shape == (4096,2048)


