# Module to run tests on ararclines


import numpy as np
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arparse as settings
from pypit import ararc as pyarc
from pypit import arutils as arut


#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_setup_param():
    """ Run the parameter setup script
    Returns
    -------

    """
    # Initialize some settings
    arut.dummy_settings()
    # Load Dummy self
    slf = arut.dummy_self()
    settings.argflag['run']['spectrograph'] = 'kast_blue'
    settings.spect['arc'] = {}
    settings.spect['arc']['index'] = [[0]]
    fitsdict = arut.dummy_fitsdict()
    # Run
    arcparm = pyarc.setup_param(slf, 0, 1, fitsdict)
    for key in ['llist','disp','wvmnx']:
        assert key in arcparm

def test_detect_lines():
