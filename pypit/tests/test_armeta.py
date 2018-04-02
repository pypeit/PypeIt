# Module to run tests on ararclines


import os
import numpy as np
import pytest

from pypit import armeta


#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_instr_list():
    """ Instrument list
    """
    instr_list = armeta.instr_list()
    assert isinstance(instr_list, list)
    assert 'shane_kast_red' in instr_list

