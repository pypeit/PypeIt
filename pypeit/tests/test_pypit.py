# Module to run tests on arsave
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import pytest

from pypit import pypit
from pypit import msgs
from pypit.par.util import make_pypit_file
from pypit import pypitsetup

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_initialization():
    """ Load input PYPIT file
    """
    # Generate a PYPIT file
    pypit_file = data_path('test.pypit')
    make_pypit_file(pypit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)
    # Perform the setup
    setup = pypitsetup.PypitSetup.from_pypit_file(pypit_file)
    par, spectrograph, fitstbl, setup_dict = setup.run()
    # Test
    assert spectrograph.spectrograph == 'shane_kast_blue'
    assert len(fitstbl) == 2
    assert par['calibrations']['arcframe']['number'] == 1

