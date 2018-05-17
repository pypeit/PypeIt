# Module to run tests on ProcessImages class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from pypit import processimages
from pypit import biasframe

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def deimos_flat_files():
    deimos_flat_files = [os.getenv('PYPIT_DEV') + '/RAW_DATA/Keck_DEIMOS/830G_L/' + ifile for ifile in [  # Longslit in dets 3,7
        'd0914_0014.fits', 'd0914_0015.fits']]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
def kast_blue_bias_files():
    kast_blue_bias_files = glob.glob(os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b1?.fits*')
    return kast_blue_bias_files


@pytest.fixture
def kast_settings():
    kast_settings = processimages.default_settings.copy()
    kast_settings['detector']['dataext'] = 0
    kast_settings['bias'] = {}  # This is a kludge
    kast_settings['bias']['combine'] = kast_settings['combine']  # This is a kludge
    kast_settings['bias']['useframe'] = 'bias'  # For the run() method only
    return kast_settings


def test_instantiate(kast_blue_bias_files, kast_settings):
    # Empty
    bias_frame0 = biasframe.BiasFrame()
    assert bias_frame0.nfiles == 0
    #
    bias_frame1 = biasframe.BiasFrame(settings=kast_settings, file_list=kast_blue_bias_files)
    assert bias_frame1.nfiles == 10

def test_process(kast_blue_bias_files, kast_settings):
    if skip_test:
        assert True
        return
    # Instantiate
    bias_frame = biasframe.BiasFrame(spectrograph='shane_kast_blue', settings=kast_settings, file_list=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.process_bias()
    assert isinstance(bias_img, np.ndarray)
    assert isinstance(bias_frame.stack, np.ndarray)
    assert bias_frame.steps[-1] == 'combine'


def test_read_write(kast_blue_bias_files, kast_settings):
    if skip_test:
        assert True
        return
    # Instantiate
    bias_frame = biasframe.BiasFrame(spectrograph='shane_kast_blue', settings=kast_settings, file_list=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.process_bias()
    # Write
    bias_frame.write_stack_to_fits(data_path('tmp.fits'))
    # Read
    bias_frame2 = biasframe.BiasFrame.from_fits(data_path('tmp.fits'))
    assert np.array_equal(bias_frame2.stack, bias_img)


def test_master(kast_blue_bias_files, kast_settings):
    if skip_test:
        assert True
        return
    # Additional settings
    kast_settings['reduce'] = {}
    kast_settings['reduce']['masters'] = {}
    kast_settings['reduce']['masters']['reuse'] = False
    kast_settings['reduce']['masters']['force'] = False
    kast_settings['reduce']['masters']['loaded'] = []
    #
    kast_settings['run'] = {}
    kast_settings['run']['spectrograph'] = 'shane_kast_blue'
    kast_settings['run']['directory'] = {}
    kast_settings['run']['directory']['master'] = data_path('MF')
    setup = 'A_01_aa'
    # Instantiate
    kast_settings['reduce']['masters']['reuse'] = False
    bias_frame = biasframe.BiasFrame(settings=kast_settings, file_list=kast_blue_bias_files, setup=setup)
    # Run
    _ = bias_frame.build_master()
    assert bias_frame.steps[-1] == 'combine'
    # Run with reuse (should simply load the file)
    kast_settings['reduce']['masters']['reuse'] = True
    bias_frame2 = biasframe.BiasFrame(settings=kast_settings, setup=setup)
    bias2 = bias_frame2.build_master()
    assert isinstance(bias_frame2.stack, np.ndarray)
    assert len(bias_frame2.steps) == 0
    # Load (not kept in the Object!)
    bias_frame3 = biasframe.BiasFrame(settings=kast_settings, setup=setup)
    bias3, _, _ = bias_frame3.load_master_frame()
    assert bias_frame3.stack is None
    assert np.array_equal(bias2, bias3)

# Should probaby test overscan
