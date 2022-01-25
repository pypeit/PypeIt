"""
Module to run tests on BiasFrame class
Requires files in Development suite and an Environmental variable
"""
import os

from IPython import embed

import pytest
import glob
import numpy as np

from pypeit.images import buildimage
from pypeit.images.mosaic import Mosaic
from pypeit.tests.tstutils import dev_suite_required, data_path
from pypeit.spectrographs.util import load_spectrograph
from pypeit import masterframe

# Init a few things
shane_kast_blue = load_spectrograph('shane_kast_blue')
frame_par = shane_kast_blue.default_pypeit_par()['calibrations']['biasframe']
master_key = 'A_1_DET01'
master_dir = data_path('')

@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'),
                                      '/RAW_DATA/keck_deimos/830G_L/', ifile) \
                            for ifile in ['d0914_0014.fits', 'd0914_0015.fits'] ]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    kast_blue_files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                  'shane_kast_blue', '600_4310_d55', 'b1?.fits*'))
    kast_blue_files.sort()
    # Trim to bias
    kast_blue_bias_files = kast_blue_files[5:]
    return kast_blue_bias_files

#@dev_suite_required
#def test_instantiate(kast_blue_bias_files):
#    # Empty
#    bias_frame0 = biasframe.BiasFrame(shane_kast_blue)
#    assert bias_frame0.nfiles == 0
#    #
#    bias_frame1 = biasframe.BiasFrame(shane_kast_blue, files=kast_blue_bias_files)
#    assert bias_frame1.nfiles == 5

@dev_suite_required
def test_process(kast_blue_bias_files):
    # Instantiate
    bias_img = buildimage.buildimage_fromlist(shane_kast_blue, 1, frame_par,
                                              kast_blue_bias_files)
    assert isinstance(bias_img.image, np.ndarray)


@dev_suite_required
def test_io(kast_blue_bias_files):
    #
    outfile = masterframe.construct_file_name(buildimage.BiasImage, master_key,
                                              master_dir=master_dir)
    if os.path.isfile(outfile):
        os.remove(outfile)
    # Build
    msbias = buildimage.buildimage_fromlist(shane_kast_blue, 1, frame_par,
                                            kast_blue_bias_files)
    # Save as a master frame
    master_filename = masterframe.construct_file_name(msbias, master_key, master_dir=master_dir)
    msbias.to_master_file(master_filename)

    assert os.path.isfile(outfile), 'Error writing MasterBias'
    # Load master frame
    biasImage = buildimage.BiasImage.from_file(outfile)
    assert np.array_equal(biasImage.image, msbias.image), 'Image changed'
    assert np.array_equal(biasImage.ivar, msbias.ivar), 'Inverse-variance changed'
    # Clean up
    os.remove(outfile)


@dev_suite_required
def test_process_multidet():
    files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                                   'GN_HAM_R400_885', 'N20190205S024*.fits'))
    files.sort()
    spec = load_spectrograph('gemini_gmos_north_ham')
    frame_par = spec.default_pypeit_par()['calibrations']['biasframe']

    det = 1
    bias_img_det1 = buildimage.buildimage_fromlist(spec, det, frame_par, files)

    det = (1,2,3)
    bias_img = buildimage.buildimage_fromlist(spec, det, frame_par, files)

    assert np.array_equal(bias_img_det1.image, bias_img.image[0]) \
            and not np.array_equal(bias_img_det1.image, bias_img.image[1]) \
            and not np.array_equal(bias_img_det1.image, bias_img.image[2]), \
                'Bad multi-detector processing'

@dev_suite_required
def test_mosaic_io():
    outfile = data_path('test_bias_mosaic.fits')
    if os.path.isfile(outfile):
        os.remove(outfile)
    files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'gemini_gmos',
                                   'GN_HAM_R400_885', 'N20190205S024*.fits'))
    files.sort()
    spec = load_spectrograph('gemini_gmos_north_ham')
    frame_par = spec.default_pypeit_par()['calibrations']['biasframe']

    det = (1,2,3)
    bias = buildimage.buildimage_fromlist(spec, det, frame_par, files)
    bias.to_file(outfile)

    _bias = buildimage.BiasImage.from_file(outfile)
    assert isinstance(_bias.detector, Mosaic), 'Bad detector type'
    assert _bias.detector.ndet == bias.detector.ndet, 'Bad number of detectors'
    assert np.array_equal(_bias.detector.tform, bias.detector.tform), 'Bad number of detectors'

    # Clean up
    os.remove(outfile)


