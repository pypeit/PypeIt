"""
Module to run tests on PypeItImage class
"""
import os

import pytest
import glob
import numpy as np

from astropy.io import fits

from pypeit.images import pypeitimage
from pypeit import arcimage
from pypeit.tests.tstutils import dev_suite_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.mask.fullmask = np.zeros((1000, 1000), dtype=np.int64)

    # Now the arcimage
    arcImage = arcimage.ArcImage.from_pypeitimage(pypeitImage)
    pytest.set_trace()
