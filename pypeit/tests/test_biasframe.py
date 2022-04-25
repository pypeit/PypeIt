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
