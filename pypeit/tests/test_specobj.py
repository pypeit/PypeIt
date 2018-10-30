# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

### TEST_UNICODE_LITERALS

import numpy as np
import sys
import os
import pytest

#import pypeit

from astropy.table import Table

from pypeit import specobjs
from pypeit import msgs

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

# def test_load_specobj -- See test_arload.py

objnm1 = 'SPAT0968-SLIT0001-DET01-SCI023'
objnm2 = 'SPAT0967-SLIT0001-DET01-SCI023'
objnm3 = 'SPAT0965-SLIT0001-DET01-SCI028'

def test_objnm_to_dict():
    idict = specobjs.objnm_to_dict(objnm1)
    assert 'SPAT' in idict.keys()
    assert idict['SPAT'] == 968
    assert 'SLIT' in idict.keys()
    assert idict['SLIT'] == 1
    assert 'DET' in idict.keys()
    assert idict['DET'] == 1
    assert 'SCI' in idict.keys()
    assert idict['SCI'] == 23
    # List
    idict2 = specobjs.objnm_to_dict([objnm1,objnm2])
    assert len(idict2['SPAT']) == 2
    assert idict2['SPAT'] == [968, 967]


def test_findobj():
    objects = [objnm1, objnm2]
    mtch_obj, indices = specobjs.mtch_obj_to_objects(objnm3, objects)
    assert mtch_obj == objects
    assert indices == [0,1]
    # Now hit only 1
    #mtch_obj2, _ = specobjs.mtch_obj_to_objects('O965-S5338-D01-I0028', objects)
    #assert mtch_obj2[0] == 'O968-S5387-D01-I0026'


def test_instr_config():
    # Make dummy fitsdict
    fitsdict = {'slitwid': [0.5], 'dichroic': ['d55'],
                 'dispname': ['B600/400'], 'dispangle': [11000.]}
    fitstbl = Table(fitsdict)
    det, scidx = 1, 0
    #
    config = specobjs.instconfig(fitstbl[scidx])
    # Test
    assert config == 'S05-D55-G600400-T110000-B11'

