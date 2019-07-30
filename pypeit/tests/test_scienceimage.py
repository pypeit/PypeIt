"""
Module to run tests on SciImgStack class
Requires files in Development suite
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import dev_suite_required
from pypeit.tests.tstutils import load_kast_blue_masters
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import scienceimage


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
    sciImg = scienceimage.ScienceImage(kast_blue, [], kast_par['scienceframe']['process'], bpm)
    assert sciImg.nfiles == 0


@dev_suite_required
def test_instantiate_from_one(shane_kast_blue_sci_files):
    """
    Run on a single science frame
    """
    det = 1
    # Load calibrations
    tslits_dict, mstrace, tilts_dict, pixelflat = load_kast_blue_masters(
        tslits=True, tilts=True, pixflat=True)
    bpm = kast_blue.empty_bpm(shane_kast_blue_sci_files[0], det)
    # Do it
    sciImg = scienceimage.ScienceImage.from_single_file(kast_blue, det, kast_par['scienceframe']['process'],
                                            bpm, shane_kast_blue_sci_files[0], None, pixelflat)
    # Test
    assert sciImg.nfiles == 1


@dev_suite_required
def test_from_list(shane_kast_blue_sci_files):
    """
    Run on two frames
    """
    det = 1
    # Load calibrations
    tslits_dict, mstrace, tilts_dict, pixelflat = load_kast_blue_masters(
        tslits=True, tilts=True, pixflat=True)
    bpm = kast_blue.empty_bpm(shane_kast_blue_sci_files[0], det)
    # Do it
    sciImg = scienceimage.ScienceImage.from_file_list(kast_blue, det, kast_par['scienceframe']['process'],
                                                        bpm, shane_kast_blue_sci_files, None, pixelflat)
    # Test
    assert sciImg.nfiles == 2


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
    sciImg = scienceimage.ScienceImage.from_file_list(keck_nires, det, nires_par['scienceframe']['process'], bpm,
                                       nires_sci_files, None, pixelflat)
    # Bg image
    bgImg = scienceimage.ScienceImage.from_file_list(keck_nires, det, nires_par['scienceframe']['process'], bpm,
                                                      nires_bg_files, None, pixelflat)
    # Difference
    sciImg = sciImg - bgImg
