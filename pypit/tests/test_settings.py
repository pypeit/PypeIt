# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest
import glob
import filecmp

import pypit
from pypit import pyputils

msgs = pyputils.get_dummy_logger()

sys.path.append(os.path.abspath(pypit.__path__[0]+"/data/settings"))
import settings as py_sett

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_baseargflag():
    """ Test that the current settings.baseargflag matches the
    most recent archived one.  This avoids our changing the former
    without careful consideration (I hope)
    """
    # Base
    settings_path = pypit.__path__[0]+'/data/settings/'
    baseargf_file = settings_path+'settings.baseargflag'
    # Archive
    archive_path = pypit.__path__[0]+'/data/settings/archive/'
    arch_file = py_sett.current_arch_file(archive_path)
    # Test
    assert filecmp.cmp(baseargf_file, arch_file)

