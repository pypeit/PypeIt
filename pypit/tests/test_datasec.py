# Module to run tests on ampsec definition
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy as np

from pypit import arparse as settings
from pypit import arpixels
from pypit.spectrographs import io
from pypit.core import arsort


@pytest.fixture
def fitstbl():
    return arsort.dummy_fitstbl()


def test_ampsec(fitstbl):
    """ Test sort_data
    """
    spectrograph='shane_kast_blue'
    settings.dummy_settings(spectrograph=spectrograph)
    # Run
    dnum = 'det01'
    settings_det = settings.spect[dnum].copy()
    settings_det['naxis0'] = 2112  # Raw frame, with overscan
    settings_det['naxis1'] = 350
    #datasec_img, naxis0, naxis1 = arprocimg.get_datasec_trimmed(
    #    settings.argflag['run']['spectrograph'], None, det, settings_det,
    #    naxis0=fitstbl['naxis0'][scidx],
    #    naxis1=fitstbl['naxis1'][scidx])
    datasec, _ = io.get_datasec(spectrograph, filename=None, det_settings=settings_det,
                                numamplifiers=settings_det['numamplifiers'], det=1)
    pytest.set_trace()
    datasec_img = arpixels.pix_to_amp(settings_det['naxis0'],
                                      settings_det['naxis1'],
                                      datasec, settings_det['numamplifiers'])
    # Test
    assert datasec_img.shape == (2112, 2048)
    assert np.sum(np.isclose(datasec_img, 1)) == 2162688  # Data region
    assert np.sum(np.isclose(datasec_img, 2)) == 2162688  # second amp
    assert settings.spect[dnum]['oscansec01'] == [[0, 0], [2049, 2080]]
    assert settings.spect[dnum]['datasec01'] == [[0, 0], [0, 1024]]

