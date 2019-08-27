"""
Module to run tests on BiasFrame class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit import biasframe
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph

shane_kast_blue = load_spectrograph('shane_kast_blue')

def data_root():
    return os.path.join(os.path.dirname(__file__), 'files')

@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'),
                                      '/RAW_DATA/Keck_DEIMOS/830G_L/', ifile) \
                            for ifile in ['d0914_0014.fits', 'd0914_0015.fits'] ]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    kast_blue_bias_files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                  'Shane_Kast_blue', '600_4310_d55', 'b1?.fits*'))
    return kast_blue_bias_files

@dev_suite_required
def test_instantiate(kast_blue_bias_files):
    # Empty
    bias_frame0 = biasframe.BiasFrame(shane_kast_blue)
    assert bias_frame0.nfiles == 0
    #
    bias_frame1 = biasframe.BiasFrame(shane_kast_blue, files=kast_blue_bias_files)
    assert bias_frame1.nfiles == 10

@dev_suite_required
def test_process(kast_blue_bias_files):
    # Instantiate
    bias_frame = biasframe.BiasFrame(shane_kast_blue, files=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.build_image()
    assert isinstance(bias_img.image, np.ndarray)
    assert isinstance(bias_frame.pypeitImage.image, np.ndarray)
    #assert bias_frame.steps[-1] == 'combine'


@dev_suite_required
def test_io(kast_blue_bias_files):
    # Instantiate
    bias_frame = biasframe.BiasFrame(shane_kast_blue, files=kast_blue_bias_files,
                                     master_dir=data_root(), master_key='A_01_1', 
                                     reuse_masters=True)
    # In case of previous test failure
    if os.path.isfile(bias_frame.master_file_path):
        os.remove(bias_frame.master_file_path)
    # Run
    bias_frame.build_image()
    # Save as a master frame
    bias_frame.save()
    assert os.path.isfile(bias_frame.master_file_path), 'Error writing MasterBias'
    # Load master frame
    pypeitImage = bias_frame.load()
    assert np.array_equal(pypeitImage.image, bias_frame.pypeitImage.image)
    # Instantiate from master frame
    bias_frame2 = biasframe.BiasFrame.from_master_file(bias_frame.master_file_path)
    assert np.array_equal(pypeitImage.image, bias_frame2.pypeitImage.image)
    # Clean up
    os.remove(bias_frame.master_file_path)


@dev_suite_required
def test_run_and_master(kast_blue_bias_files):
    # Instantiate
    bias_frame = biasframe.BiasFrame(shane_kast_blue, files=kast_blue_bias_files,
                                     master_key='A_1_01', master_dir=data_root())
    assert bias_frame.frametype == 'bias'
    # In case of previous test failure
    if os.path.isfile(bias_frame.master_file_path):
        os.remove(bias_frame.master_file_path)

    # Run
    msbias = bias_frame.build_image()
    #assert bias_frame.steps[-1] == 'combine'
    # Save it
    bias_frame.save()

    # Run with reuse (should simply load the file)
    bias_frame2 = biasframe.BiasFrame(shane_kast_blue, master_key='A_1_01',
                                      master_dir=data_root(), reuse_masters=True)
    bias2 = bias_frame2.load()
    assert isinstance(bias2.image, np.ndarray)
    assert len(bias_frame2.process_steps) == 3
    assert np.array_equal(bias2.image, bias_frame.pypeitImage.image)

    # Clean up
    os.remove(bias_frame.master_file_path)

# Should probably test overscan

