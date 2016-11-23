# Module to run tests on arsave

import os
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_load_input():
    """ Load input PYPIT file
    """
    from pypit import pypit
    # Generate a PYPIT file
    pyp_file = data_path('test.pypit')
    pyputils.make_pypit_file(pyp_file, 'kast_blue', data_path('b'), 'fits')
    parlines, datlines, spclines, dfnames = pypit.load_input(pyp_file, msgs)
    # Test
    assert len(parlines) == 4
    assert 'b1.fits' in datlines[0]
    assert 'arc number 1' in spclines[1]

