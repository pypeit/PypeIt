"""
Module to run tests on SciImgStack class
Requires files in Development suite
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import dev_suite_required, cooked_required
from pypeit.tests.tstutils import load_kast_blue_masters
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import scienceimage
from pypeit.images import rawimage
from pypeit.images import processrawimage
from pypeit.core import procimg


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

kast_blue = load_spectrograph('shane_kast_blue')
kast_par = kast_blue.default_pypeit_par()
keck_nires = load_spectrograph('keck_nires')
nires_par = keck_nires.default_pypeit_par()

@pytest.fixture
@dev_suite_required
def shane_kast_blue_sci_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Shane_Kast_blue', '600_4310_d55', ifile)
            for ifile in ['b27.fits.gz', 'b28.fits.gz']]

@pytest.fixture
@dev_suite_required
def nires_sci_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRES', 'NIRES', ifile)
            for ifile in ['s180604_0089.fits.gz', 's180604_0092.fits.gz']]

@pytest.fixture
@dev_suite_required
def nires_bg_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRES', 'NIRES', ifile)
            for ifile in ['s180604_0090.fits.gz', 's180604_0091.fits.gz']]


def test_standard_instantiate():
    """
    Simple instantiate
    """
    # Empty
    bpm = np.zeros((100,100))
    det = 1
    sciImg = scienceimage.ScienceImage(kast_blue, det, kast_par['scienceframe']['process'],
                                       np.zeros_like(bpm),
                                       np.zeros_like(bpm),
                                       bpm)


@cooked_required
def test_instantiate_from_one(shane_kast_blue_sci_files):
    """
    Run on a single science frame
    """
    #
    det = 1
    # Load calibrations
    pixelflat = load_kast_blue_masters(pixflat=True)[0]
    bpm = kast_blue.empty_bpm(shane_kast_blue_sci_files[0], det)
    # Process steps
    bias = None
    par = kast_par['scienceframe']['process']
    process_steps = procimg.init_process_steps(bias, par)
    process_steps += ['trim', 'apply_gain', 'orient']
    process_steps += ['flatten']
    process_steps += ['extras']
    process_steps += ['crmask']
    # Load
    rawImage = rawimage.RawImage(shane_kast_blue_sci_files[0], kast_blue, det)
    processRawImage = processrawimage.ProcessRawImage(rawImage,
                                                      kast_par['scienceframe']['process'])
    pypeItImage = processRawImage.process(process_steps, pixel_flat=pixelflat)
    # Do it
    sciImg = scienceimage.ScienceImage(kast_blue, det, kast_par['scienceframe']['process'],
                                       pypeItImage.image, pypeItImage.ivar, bpm)


@cooked_required
def test_from_list(shane_kast_blue_sci_files):
    """
    Run on two frames
    """
    det = 1
    # Load calibrations
    pixelflat = load_kast_blue_masters(pixflat=True)[0]
    bpm = kast_blue.empty_bpm(shane_kast_blue_sci_files[0], det)
    # Do it
    sciImg = scienceimage.build_from_file_list(kast_blue, det,
                                               kast_par['scienceframe']['process'],
                                               bpm, shane_kast_blue_sci_files,
                                               None, pixelflat)
    # Test
    assert isinstance(sciImg, scienceimage.ScienceImage)


@dev_suite_required
def test_proc_diff(nires_sci_files, nires_bg_files):
    """
    Run on near-IR frames
    """
    # Setup
    det = 1
    bpm = np.zeros((2048,1024))
    pixelflat = np.ones_like(bpm)

    # Sci image
    sciImg = scienceimage.build_from_file_list(keck_nires, det, nires_par['scienceframe']['process'], bpm,
                                       nires_sci_files, None, pixelflat)
    # Bg image
    bgImg = scienceimage.build_from_file_list(keck_nires, det, nires_par['scienceframe']['process'], bpm,
                                                      nires_bg_files, None, pixelflat)
    # Difference
    sciImg = sciImg - bgImg
    # Test
    assert isinstance(sciImg, scienceimage.ScienceImage)
