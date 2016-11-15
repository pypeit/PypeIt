# Module to run tests on simple fitting routines for arrays

### TEST_UNICODE_LITERALS

import numpy as np
import sys
import os, pdb
import pytest

from pypit import pyputils
import pypit
msgs = pyputils.get_dummy_logger()

#from xastropy.xutils import afits as xafits
#from xastropy.xutils import xdebug as xdb

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_objnm_to_dict():
    from pypit import arspecobj as aspobj
    idict = aspobj.objnm_to_dict('O968-S5387-D01-I0026')
    assert 'O' in idict.keys()
    assert idict['O'] == 968
    # List
    idict2 = aspobj.objnm_to_dict(['O968-S5387-D01-I0026', 'O967-S5397-D01-I0026'])
    assert len(idict2['O']) == 2
    assert idict2['O'] == [968, 967]

def test_findobj():
    from pypit import arspecobj as aspobj
    objects = ['O968-S5387-D01-I0026', 'O967-S5397-D01-I0027']
    mtch_obj, indices = aspobj.mtch_obj_to_objects('O965-S5390-D01-I0028', objects)
    assert mtch_obj == objects
    assert indices == [0,1]
    # Now hit only 1
    mtch_obj2, _ = aspobj.mtch_obj_to_objects('O965-S5338-D01-I0028', objects)
    assert mtch_obj2[0] == 'O968-S5387-D01-I0026'

