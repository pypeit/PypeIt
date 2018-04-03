# Module to run tests on scripts
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest
from glob import glob
import filecmp

import pypit
from pypit import pyputils

msgs = pyputils.get_dummy_logger()

sys.path.append(os.path.abspath(pypit.__path__[0]+"/data/settings"))
import settings as py_sett


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_settings_vs_archive():
    """ Test that the current settings.base files match the
    most recent archived one.  This avoids our changing the former
    without careful consideration (I hope)

    If this test fails, you likely need to run data/settings/settings.archive()
    """
    # Base
    settings_path = pypit.__path__[0]+'/data/settings/'
    archive_path = pypit.__path__[0]+'/data/settings/archive/'
    sett_files = glob(settings_path+'settings.*')
    for sfile in sett_files:
        # Extension
        ext = sfile.split('.')[-1]
        if ext in ['py', 'pyc']:
            continue
        # Archive
        arch_file = py_sett.current_sett_file(archive_path, sfile)
        # Test
        assert filecmp.cmp(sfile, arch_file)
