# Module to run tests on armsld module

# TEST_UNICODE_LITERALS

import numpy as np
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_instr_config():
    from pypit import arutils as arut
    from pypit import armlsd
    # Dummy self
    slf = arut.dummy_self()
    # Grab
    fits_dict = {'slitwid': [0.5], 'dichroic': ['d55'],
                 'dispname': ['B600/400'], 'dispangle': [11000.]}
    det, scidx = 1, 0
    #
    config = armlsd.instconfig(slf, det, scidx, fits_dict)
    # Test
    assert config == 'S05-D55-G600400-T110000-B11'

