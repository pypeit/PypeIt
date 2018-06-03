# Module to run tests on ampsec definition
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy as np

from pypit import arparse as settings
from pypit.core import arprocimg
from pypit.core import arsort


@pytest.fixture
def fitstbl():
    return arsort.dummy_fitstbl()


def test_ampsec(fitstbl):
    """ Test sort_data
    """
    settings.dummy_settings(spectrograph='shane_kast_blue')
    # Run
    namp, det, scidx = 2, 1, 5
    dnum = 'det01'
    settings_det = settings.spect[dnum].copy()  # Should include naxis0, naxis1 in this
    datasec_img, naxis0, naxis1 = arprocimg.get_datasec_trimmed(
        settings.argflag['run']['spectrograph'], None, namp, det, settings_det,
        naxis0=fitstbl['naxis0'][scidx],
        naxis1=fitstbl['naxis1'][scidx])
    # Test
    assert datasec_img.shape == (2112, 2048)
    assert np.sum(np.isclose(datasec_img, 1)) == 2162688  # Data region
    assert np.sum(np.isclose(datasec_img, 2)) == 2162688  # second amp
    assert settings.spect[dnum]['oscansec01'] == [[0, 0], [2049, 2080]]
    assert settings.spect[dnum]['datasec01'] == [[0, 0], [0, 1024]]

