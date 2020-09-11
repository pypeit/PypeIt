"""
Module to run tests on BiasFrame class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph
from pypeit import masterframe

def data_root():
    return os.path.join(os.path.dirname(__file__), 'files')

# Init a few things
shane_kast_blue = load_spectrograph('shane_kast_blue')
frame_par = shane_kast_blue.default_pypeit_par()['calibrations']['biasframe']
master_key = 'A_1_01'
master_dir = data_root()

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
    msbias.to_master_file(master_filename)#master_dir, master_key,  # Naming
                               #shane_kast_blue.spectrograph,  # Header
                               #steps=msbias.process_steps,
                               #raw_files=kast_blue_bias_files)

    assert os.path.isfile(outfile), 'Error writing MasterBias'
    # Load master frame
    biasImage = buildimage.BiasImage.from_file(outfile)
    assert np.array_equal(biasImage.image, msbias.image)
    # Instantiate from master frame
    #bias_frame2 = biasframe.BiasFrame.from_master_file(bias_frame.master_file_path)
    #assert np.array_equal(pypeitImage.image, bias_frame2.pypeitImage.image)
    # Clean up
    os.remove(outfile)

