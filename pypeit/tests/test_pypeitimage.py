"""
Module to run tests on PypeItImage class
"""
import os

import pytest
import glob
import numpy as np

from astropy.io import fits

from pypeit.images import pypeitimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs.util import load_spectrograph

kast_blue = load_spectrograph('shane_kast_blue')

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Shane_Kast_blue',
                                  '600_4310_d55', 'b1?.fits*'))

def test_instantiate():
    pypeitImage = pypeitimage.PypeItImage(kast_blue, 1)
    assert pypeitImage.head0 is None

@dev_suite_required
def test_load(kast_blue_bias_files):
    one_file = kast_blue_bias_files[0]
    #
    pypeitImage = pypeitimage.PypeItImage(kast_blue, 1)


