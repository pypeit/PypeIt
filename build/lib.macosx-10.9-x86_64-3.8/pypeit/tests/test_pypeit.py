"""
Module to run tests on arsave
"""
import os

import numpy as np

import pytest

from pypeit import msgs
from pypeit.par.util import make_pypeit_file
from pypeit import pypeitsetup
from pypeit.pypeit import PypeIt

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_initialization():
    """ Load input PYPIT file
    """
    # Generate a PYPIT file
    pypit_file = data_path('test.pypeit')
    make_pypeit_file(pypit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypit_file)
    par, spectrograph, fitstbl = setup.run(sort_dir=data_path(''))

    # Test
    assert spectrograph.name == 'shane_kast_blue'
    assert len(fitstbl) == 2

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))

def test_select_detectors():
    # Generate a PYPIT file
    pypit_file = data_path('test.pypeit')
    make_pypeit_file(pypit_file, 'shane_kast_blue', [data_path('b*fits.gz')], setup_mode=True)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypit_file)
    par, spectrograph, fitstbl = setup.run(sort_dir=data_path(''))

    assert PypeIt.select_detectors(detnum=par['rdx']['detnum'], ndet=spectrograph.ndet) == [1], \
            'Incorrect detectors selected.'

    # Clean-up
    os.remove(data_path('test.calib'))
    os.remove(data_path('test.pypeit'))

    assert np.array_equal(PypeIt.select_detectors(), [1]), 'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(detnum=3, ndet=5), [3]), \
            'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(ndet=5), [1,2,3,4,5]), \
            'Incorrect detectors selected.'
    assert np.array_equal(PypeIt.select_detectors(detnum=[1,3]), [1,3]), \
            'Incorrect detectors selected.'

