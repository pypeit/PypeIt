# Module to run tests on ampsec definition

import pytest

import numpy as np
from pypit import pyputils
msgs = pyputils.get_dummy_logger()#develop=True)
from pypit import arproc
from pypit import arutils
from pypit import arparse as settings

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


@pytest.fixture
def fitsdict():
    return arutils.dummy_fitsdict()

def test_ampsec(fitsdict):
    """ Test sort_data
    """
    arutils.dummy_settings(spectrograph='kast_blue')
    slf = arutils.dummy_self()
    # Run
    det, scidx = 1, 5
    arproc.get_ampsec_trimmed(slf, fitsdict, det, scidx)
    # Test
    assert slf._datasec[det-1].shape == (2112, 2048)
    assert np.sum(np.isclose(slf._datasec[0],1)) == 2162688  # Data region
    assert np.sum(np.isclose(slf._datasec[0],2)) == 2162688  #   second amp
    assert settings.spect['det01']['oscansec01'] == [[0, 0], [2049, 2080]]
    assert settings.spect['det01']['datasec01'] == [[0, 0], [0, 1024]]

