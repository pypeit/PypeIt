# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import numpy as np
import sys
import os, pdb
import pytest

#import pypit

from pypit import arspecobj
from pypit import pyputils

msgs = pyputils.get_dummy_logger()

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

# def test_load_specobj -- See test_arload.py

def test_objnm_to_dict():
    idict = arspecobj.objnm_to_dict('O968-S5387-D01-I0026')
    assert 'O' in idict.keys()
    assert idict['O'] == 968
    assert 'S' in idict.keys()
    assert idict['S'] == 5387
    assert 'D' in idict.keys()
    assert idict['D'] == 1
    assert 'I' in idict.keys()
    assert idict['I'] == 26
    # List
    idict2 = arspecobj.objnm_to_dict(['O968-S5387-D01-I0026', 'O967-S5397-D01-I0026'])
    assert len(idict2['O']) == 2
    assert idict2['O'] == [968, 967]


def test_findobj():
    objects = ['O968-S5387-D01-I0026', 'O967-S5397-D01-I0027']
    mtch_obj, indices = arspecobj.mtch_obj_to_objects('O965-S5390-D01-I0028', objects)
    assert mtch_obj == objects
    assert indices == [0,1]
    # Now hit only 1
    mtch_obj2, _ = arspecobj.mtch_obj_to_objects('O965-S5338-D01-I0028', objects)
    assert mtch_obj2[0] == 'O968-S5387-D01-I0026'


def test_instr_config():
    # Make dummy fitsdict
    fitsdict = {'slitwid': [0.5], 'dichroic': ['d55'],
                 'dispname': ['B600/400'], 'dispangle': [11000.]}
    det, scidx = 1, 0
    #
    config = arspecobj.instconfig(det, scidx, fitsdict)
    # Test
    assert config == 'S05-D55-G600400-T110000-B11'

