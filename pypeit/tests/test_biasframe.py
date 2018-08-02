# Module to run tests on BiasFrame class
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

from pypeit import processimages
from pypeit import biasframe

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@pytest.fixture
def deimos_flat_files():
    if not skip_test:
        # Longslit in dets 3,7
        deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'),
                                          '/RAW_DATA/Keck_DEIMOS/830G_L/', ifile) \
                                    for ifile in ['d0914_0014.fits', 'd0914_0015.fits'] ]
        assert len(deimos_flat_files) == 2
    else:
        deimos_flat_files = None
    return deimos_flat_files

@pytest.fixture
def kast_blue_bias_files():
    if not skip_test:
        kast_blue_bias_files = glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA',
                                                      'Shane_Kast_blue', '600_4310_d55',
                                                      'b1?.fits*'))
    else:
        kast_blue_bias_files = None
    return kast_blue_bias_files

def test_instantiate(kast_blue_bias_files):
    if skip_test:
        assert True
        return
    # Empty
    bias_frame0 = biasframe.BiasFrame('shane_kast_blue')
    assert bias_frame0.nfiles == 0
    #
    bias_frame1 = biasframe.BiasFrame('shane_kast_blue', file_list=kast_blue_bias_files)
    assert bias_frame1.nfiles == 10

def test_process(kast_blue_bias_files):
    if skip_test:
        assert True
        return
    # Instantiate
    bias_frame = biasframe.BiasFrame('shane_kast_blue', file_list=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.build_image()
    assert isinstance(bias_img, np.ndarray)
    assert isinstance(bias_frame.stack, np.ndarray)
    assert bias_frame.steps[-1] == 'combine'


def test_read_write(kast_blue_bias_files):
    if skip_test:
        assert True
        return
    # Instantiate
    bias_frame = biasframe.BiasFrame('shane_kast_blue', file_list=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.build_image()
    # Write
    bias_frame.write_stack_to_fits(data_path('tmp.fits'))
    # TODO: This is now more complicated.  Will need a whole PR to allow
    # for the `from_fits` function
#    # Read
#    bias_frame2 = biasframe.BiasFrame.from_fits(data_path('tmp.fits'))
#    assert np.array_equal(bias_frame2.stack, bias_img)


def test_run_and_master(kast_blue_bias_files):
    if skip_test:
        assert True
        return

    root_path = data_path('MF')
    setup = 'A_01_aa'
    # Instantiate
    bias_frame = biasframe.BiasFrame('shane_kast_blue', file_list=kast_blue_bias_files,
                                     setup=setup, root_path=root_path)
    assert bias_frame.frametype == 'bias'

    # Run
    msbias = bias_frame.build_image()
    bias_frame.save_master(msbias)
    assert bias_frame.steps[-1] == 'combine'

    # Run with reuse (should simply load the file)
    bias_frame2 = biasframe.BiasFrame('shane_kast_blue', setup=setup, root_path=root_path,
                                      mode='reuse')
    bias2 = bias_frame2.master()
    assert isinstance(bias_frame2.stack, np.ndarray)
    assert len(bias_frame2.steps) == 0

    # Load (not kept in the Object!)
    bias_frame3 = biasframe.BiasFrame('shane_kast_blue', setup=setup, root_path=root_path,
                                      mode='reuse')
    bias3, _, _ = bias_frame3.load_master_frame()
    assert bias_frame3.stack is None
    assert np.array_equal(bias2, bias3)

# Should probably test overscan
