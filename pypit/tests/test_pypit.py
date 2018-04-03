# Module to run tests on arsave
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

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
    pyputils.make_pypit_file(pyp_file, 'shane_kast_blue', [data_path('b*fits.gz')])
    pyp_dict = pypit.load_input(pyp_file, msgs)
    parlines, datlines, spclines, dfnames = [pyp_dict[ii] for ii in ['par','dat','spc','dfn']]
    # Test
    assert len(parlines) == 3
    assert len(datlines) == 2
    assert 'arc number 1' in spclines[0]

