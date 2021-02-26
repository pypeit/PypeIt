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
from pypeit.images import buildimage
from pypeit.images import rawimage
from pypeit.images import pypeitimage
from pypeit.core import procimg
from pypeit import flatfield


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
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue', '600_4310_d55',
                         ifile) for ifile in ['b27.fits.gz', 'b28.fits.gz']]

@pytest.fixture
@dev_suite_required
def nires_sci_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_nires', 'NIRES', ifile)
            for ifile in ['s180604_0089.fits.gz', 's180604_0092.fits.gz']]

@pytest.fixture
@dev_suite_required
def nires_bg_files():
    return [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_nires', 'NIRES', ifile)
            for ifile in ['s180604_0090.fits.gz', 's180604_0091.fits.gz']]


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
    # Process steps -- Set in PypeItPar
    frame_par = kast_par['scienceframe']
    frame_par['process']['use_illumflat'] = False
    frame_par['process']['use_biasimage'] = False
    # Load
    rawImage = rawimage.RawImage(shane_kast_blue_sci_files[0], kast_blue, det)
    flatImages = flatfield.FlatImages(pixelflat_norm=pixelflat)
    pypeItImage = rawImage.process(frame_par['process'], flatimages=flatImages)




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
    flatImages = flatfield.FlatImages(pixelflat_norm=pixelflat)
    kast_par['scienceframe']['process']['use_illumflat'] = False
    kast_par['scienceframe']['process']['use_biasimage'] = False
    sciImg = buildimage.buildimage_fromlist(kast_blue, det, kast_par['scienceframe'],
                                            shane_kast_blue_sci_files, bpm=bpm,
                                            bias=None, flatimages=flatImages)
    # Test
    assert isinstance(sciImg, pypeitimage.PypeItImage)


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
    flatImages = flatfield.FlatImages(pixelflat_norm=pixelflat)
    sciImg = buildimage.buildimage_fromlist(keck_nires, det, nires_par['scienceframe'],
                                            nires_sci_files, bias=None, bpm=bpm,
                                            flatimages=flatImages)
    # Bg image
    bgImg = buildimage.buildimage_fromlist(keck_nires, det, nires_par['scienceframe'],
                                           nires_bg_files, bias=None, bpm=bpm,
                                           flatimages=flatImages)
    # Difference
    sciImg = sciImg.sub(bgImg, nires_par['scienceframe']['process'])
    # Test
    assert isinstance(sciImg, pypeitimage.PypeItImage)
