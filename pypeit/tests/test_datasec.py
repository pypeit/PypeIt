# Module to run tests on ampsec definition
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import os
import numpy as np

from pypeit import arpixels
from pypeit.core import arsort
from pypeit.core import arprocimg
from pypeit.spectrographs.util import load_spectrograph

@pytest.fixture
def spectrograph():
    return load_spectrograph(spectrograph='shane_kast_blue')

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_ampsec(spectrograph):
    """ Test sort_data
    """
    datasec_img = spectrograph.get_datasec_img(data_path('b1.fits.gz'), det=1)
    datasec_img = arprocimg.trim_frame(datasec_img, datasec_img < 1)
    # Test
    assert datasec_img.shape == (2048, 350)
    #assert np.sum(np.isclose(datasec_img, 1)) == 2162688  # Data region
    #assert np.sum(np.isclose(datasec_img, 2)) == 2162688  # second amp
    assert np.sum(np.isclose(datasec_img, 1)) == 358400  # Data region
    assert np.sum(np.isclose(datasec_img, 2)) == 358400  # second amp
    #assert settings.spect[dnum]['oscansec01'] == [[0, 0], [2049, 2080]]
    #assert settings.spect[dnum]['datasec01'] == [[0, 0], [0, 1024]]

