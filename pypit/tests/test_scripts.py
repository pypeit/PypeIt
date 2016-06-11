# Module to run tests on scripts

import matplotlib
matplotlib.use('Agg')  # For Travis

# TEST_UNICODE_LITERALS

import os
import pytest

from pypit import pyputils

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_arcid_plot():
    # Dummy self
    ver,upd = pyputils.get_version()
    assert isinstance(ver,basestring)



