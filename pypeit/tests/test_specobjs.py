"""
Module to run tests on SpecObjs
"""
import os

import numpy as np
import pytest

from pypeit import msgs
from pypeit import specobjs
from pypeit import specobj
msgs.reset(verbosity=2)

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

sobj1 = specobj.SpecObj('MultiSlit', 1, slitid=0)
sobj2 = specobj.SpecObj('MultiSlit', 1, slitid=1)
sobj3 = specobj.SpecObj('MultiSlit', 1, slitid=2)


def test_init():
    """ Run the parameter setup script
    """
    # Null
    sobjs1 = specobjs.SpecObjs()

    # With a few objs
    sobjs2 = specobjs.SpecObjs([sobj1,sobj2])
    assert sobjs2.nobj == 2


def test_access():
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    #
    assert sobjs[0]['PYPELINE'] == 'MultiSlit'
    assert len(sobjs['PYPELINE']) == 2

def test_add_rm():
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    sobjs.add_sobj(sobj3)
    assert sobjs.nobj == 3
    # Remove
    sobjs.remove_sobj(2)
    assert len(sobjs.specobjs) == 2

def test_set():
    sobjs = specobjs.SpecObjs([sobj1,sobj2,sobj3])
    # All
    sobjs.DET = 3
    assert np.all(sobjs[:].DET == np.array([3,3,3]))
    sobjs[:].DET = 4
    assert np.all(sobjs[:].DET == np.array([4,4,4]))
    # Slice
    sobjs[1:2].DET = 2
    assert sobjs.DET[1] == 2
    # With logic
    det2 = sobjs.DET == 2
    sobjs[det2].PYPELINE = 'BLAH'
    assert sobjs.PYPELINE[1] == 'BLAH'
    assert sobjs.PYPELINE[0] == 'MultiSlit'

