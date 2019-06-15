"""
Module to run tests on SciImgStack class
Requires files in Development suite
"""
import os

import pytest
import glob
import numpy as np

from pypeit import sciimgstack
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
    # Empty
    sciImgStack = sciimgstack.SciImgStack(kast_blue, [], kast_par['scienceframe'])
    assert sciImgStack.nfiles == 0
    assert 'nfiles' in sciImgStack.__repr__()

@dev_suite_required
def test_proc_one(shane_kast_blue_sci_files):
    """
    Instantiate and run on a single sience frame
    """
    # Load calibrations
    tslits_dict, mstrace, tilts_dict, pixelflat = load_kast_blue_masters(
        tslits=True, tilts=True, pixflat=True)
    bpm = kast_blue.empty_bpm(shape=pixelflat.shape)
    # Instantiate
    sciImgStack = sciimgstack.SciImgStack(kast_blue, [shane_kast_blue_sci_files[0]], kast_par['scienceframe'])
    # Do it
    sciImg = sciImgStack.proc(None, pixelflat, bpm)
    # Test
    assert isinstance(sciImg, scienceimage.ScienceImage)
    assert len(sciImg.files) == 1

@dev_suite_required
def test_proc_list(shane_kast_blue_sci_files):
    """
    Instantiate and run on two frames
    """
    # Load calibrations
    tslits_dict, mstrace, tilts_dict, pixelflat = load_kast_blue_masters(
        tslits=True, tilts=True, pixflat=True)
    bpm = kast_blue.empty_bpm(shape=pixelflat.shape)
    # Instantiate
    sciImgStack = sciimgstack.SciImgStack(kast_blue, shane_kast_blue_sci_files, kast_par['scienceframe'])
    # Do it
    sciImg = sciImgStack.proc(None, pixelflat, bpm)
    # Test
    assert len(sciImg.files) == 2

@dev_suite_required
def test_proc_diff(nires_sci_files, nires_bg_files):
    """
    Run on near-IR frames
    """
    # Instantiate
    sciImgStack = sciimgstack.SciImgStack(keck_nires, nires_sci_files, nires_par['scienceframe'],
                                          bg_file_list=nires_bg_files, ir_redux=True)
    # Dummy calibs
    bpm = np.zeros((2048,1024))
    pixelflat = np.ones_like(bpm)
    # Do It!
    sciImg = sciImgStack.proc(None, pixelflat, bpm)
    assert isinstance(sciImg, scienceimage.ScienceImage)
    assert len(sciImg.files) == 4
