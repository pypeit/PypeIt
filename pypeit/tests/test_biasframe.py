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
from pypeit.tests.tstutils import dev_suite_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

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
    bias_frame0 = biasframe.BiasFrame('shane_kast_blue')
    assert bias_frame0.nfiles == 0
    #
    bias_frame1 = biasframe.BiasFrame('shane_kast_blue', files=kast_blue_bias_files)
    assert bias_frame1.nfiles == 10

@dev_suite_required
def test_process(kast_blue_bias_files):
    # Instantiate
    bias_frame = biasframe.BiasFrame('shane_kast_blue', files=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.build_image()
    assert isinstance(bias_img, np.ndarray)
    assert isinstance(bias_frame.stack, np.ndarray)
    assert bias_frame.steps[-1] == 'combine'

#test_traceimage.py

@dev_suite_required
def test_read_write(kast_blue_bias_files):
    # Instantiate
    bias_frame = biasframe.BiasFrame('shane_kast_blue', files=kast_blue_bias_files)
    # Run
    bias_img = bias_frame.build_image()
    # Write
    bias_frame.write_stack_to_fits(data_path('tmp.fits'))
    # TODO: This is now more complicated.  Will need a whole PR to allow
    # for the `from_fits` function
#    # Read
#    bias_frame2 = biasframe.BiasFrame.from_fits(data_path('tmp.fits'))
#    assert np.array_equal(bias_frame2.stack, bias_img)


@dev_suite_required
def test_run_and_master(kast_blue_bias_files):

    root_path = data_path('MF')
    master_key = 'A_1_01'
    # Instantiate
    master_dir = root_path+'_shane_kast_blue'
    bias_frame = biasframe.BiasFrame('shane_kast_blue', files=kast_blue_bias_files,
                                     master_key=master_key, master_dir=master_dir)
    assert bias_frame.frametype == 'bias'

    # Run
    msbias = bias_frame.build_image()
    bias_frame.save_master(msbias)
    assert bias_frame.steps[-1] == 'combine'

    # Run with reuse (should simply load the file)
    bias_frame2 = biasframe.BiasFrame('shane_kast_blue', master_key=master_key, master_dir=master_dir,reuse_masters=True)
    bias2 = bias_frame2.master()
    assert isinstance(bias_frame2.msframe, np.ndarray)
    assert len(bias_frame2.steps) == 0

    # Load (not kept in the Object!)
    bias_frame3 = biasframe.BiasFrame('shane_kast_blue', master_key=master_key, master_dir=master_dir,reuse_masters=True)
    bias3, head = bias_frame3.load_master(bias_frame3.ms_name)
    assert bias_frame3.msframe is None
    assert np.array_equal(bias2, bias3)

# Should probably test overscan

