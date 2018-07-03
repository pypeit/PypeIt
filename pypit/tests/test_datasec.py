# Module to run tests on ampsec definition
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import os
import numpy as np

from pypit import arparse as settings
from pypit import arpixels
from pypit.spectrographs import spectro_utils
from pypit.core import arsort


@pytest.fixture
def fitstbl():
    return arsort.dummy_fitstbl()

@pytest.fixture
def spec_class():
    return spectro_utils.load_spec_class(spectrograph='generic')

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_ampsec(spec_class):
    """ Test sort_data
    """
    spectrograph='shane_kast_blue'
    settings.dummy_settings(spectrograph=spectrograph)
    # Run
    dnum = 'det01'
    settings_det = settings.spect[dnum].copy()
    settings_det['datasec01'] = [[0, 1024], [0, 0]]
    settings_det['datasec02'] = [[1024, 2048], [0, 0]]
    settings_det['oscansec01'] = [[2049, 2080], [0, 0]]
    settings_det['oscansec02'] = [[2080, 2111], [0, 0]]
    settings_det['dispaxis'] = 1
    #datasec_img, naxis0, naxis1 = arprocimg.get_datasec_trimmed(
    #    settings.argflag['run']['spectrograph'], None, det, settings_det,
    #    naxis0=fitstbl['naxis0'][scidx],
    #    naxis1=fitstbl['naxis1'][scidx])
    datasec, _, naxis0, naxis1 = spec_class.get_datasec(data_path('b1.fits.gz'),
                                                1, settings_det)
    datasec_img = arpixels.pix_to_amp(naxis0, naxis1,
                                      datasec, settings_det['numamplifiers'])
    # Test
    assert datasec_img.shape == (2048, 350)
    #assert np.sum(np.isclose(datasec_img, 1)) == 2162688  # Data region
    #assert np.sum(np.isclose(datasec_img, 2)) == 2162688  # second amp
    assert np.sum(np.isclose(datasec_img, 1)) == 358400  # Data region
    assert np.sum(np.isclose(datasec_img, 2)) == 358400  # second amp
    #assert settings.spect[dnum]['oscansec01'] == [[0, 0], [2049, 2080]]
    #assert settings.spect[dnum]['datasec01'] == [[0, 0], [0, 1024]]

