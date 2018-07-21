# Module to run tests on SpecObjs
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import numpy as np
import pytest

from pypeit import msgs
from pypeit import specobjs
msgs.reset(verbosity=2)

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

shape = (1024,2048)
sobj1 = specobjs.SpecObj(shape, (800., 810.), 1240.)
sobj2 = specobjs.SpecObj(shape, (400., 430.), 1240.)
sobj3 = specobjs.SpecObj(shape, (200., 220.), 1240.)


def test_init():
    """ Run the parameter setup script
    """
    # Null
    sobjs1 = specobjs.SpecObjs()

    # With a few objs
    sobjs2 = specobjs.SpecObjs([sobj1,sobj2])


def test_access():
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    #
    assert sobjs[0]['shape'] == shape
    assert len(sobjs['shape']) == 2

def test_add_rm():
    sobjs = specobjs.SpecObjs([sobj1,sobj2])
    sobjs.add_sobj(sobj3)
    assert len(sobjs.specobjs) == 3
    # Remove
    sobjs.remove_sobj(2)
    assert len(sobjs.specobjs) == 2
    assert len(sobjs.summary) == 2
