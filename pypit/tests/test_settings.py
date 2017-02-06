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


def test_base_settings():
    """ Test that the current settings.base files match the
    most recent archived one.  This avoids our changing the former
    without careful consideration (I hope)
    """
    # Base
    settings_path = pypit.__path__[0]+'/data/settings/'
    archive_path = pypit.__path__[0]+'/data/settings/archive/'
    for ftype in ['argflag', 'spect']:
        base_file = settings_path+'settings.base{:s}'.format(ftype)
        # Archive
        arch_file = py_sett.current_sett_file(archive_path, ftype)
        # Test
        assert filecmp.cmp(base_file, arch_file)
