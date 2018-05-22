# Module to run tests on ampsec definition
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy as np

from pypit import arparse as settings
from pypit import arproc
from pypit.core import arsort
from pypit import arsciexp


@pytest.fixture
def fitstbl():
    return arsort.dummy_fitstbl()


def test_ampsec(fitstbl):
    """ Test sort_data
    """
    settings.dummy_settings(spectrograph='shane_kast_blue')
    slf = arsciexp.dummy_self()
    # Run
    det, scidx = 1, 5
    datasec_img = arproc.get_datasec_trimmed(slf, fitstbl, det, scidx)
    # Test
    assert datasec_img.shape == (2112, 2048)
    assert np.sum(np.isclose(datasec_img, 1)) == 2162688  # Data region
    assert np.sum(np.isclose(datasec_img, 2)) == 2162688  # second amp
    assert settings.spect['det01']['oscansec01'] == [[0, 0], [2049, 2080]]
    assert settings.spect['det01']['datasec01'] == [[0, 0], [0, 1024]]

