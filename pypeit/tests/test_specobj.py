"""
Module to run tests on SpecObj
"""
import numpy as np
import sys
import os
import pytest

#import pypeit

from astropy.table import Table

from pypeit import specobj
from pypeit import msgs

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_init():
    sobj = specobj.SpecObj('MultiSlit', 1, slitid=0)
    # Test
    assert sobj.PYPELINE == 'MultiSlit'
    assert sobj['PYPELINE'] == 'MultiSlit'
    assert sobj.name == 'SPAT-----SLIT0000-DET01'
    assert len(sobj._data.keys()) == 0


def test_assignment():
    sobj = specobj.SpecObj('MultiSlit', 1, slitid=0)
    #
    sobj.PYPELINE = 'Blah'
    #
    with pytest.raises(OSError):
        sobj.PYPELINE = 2
    #
    sobj.SPAT_PIXPOS = 523.0
    sobj.PYPELINE = 'MultiSlit'
    sobj.set_name()
    assert sobj.name == 'SPAT0523-SLIT0000-DET01'


def test_data():
    sobj = specobj.SpecObj('MultiSlit', 1, slitid=0)
    #
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    sobj['BOX_COUNTS'] = np.ones_like(sobj.BOX_WAVE)
    sobj['TRACE_SPAT'] = np.arange(100) * 2.
    # Test
    assert isinstance(sobj._data, Table)
    assert 'TRACE_SPAT' in sobj._data.keys()
    assert 'SLITID' in sobj._data.meta.keys()

def test_io():
    sobj = specobj.SpecObj('MultiSlit', 1, slitid=0)
    #
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    sobj['BOX_COUNTS'] = np.ones_like(sobj.BOX_WAVE)
    sobj['TRACE_SPAT'] = np.arange(100) * 2.

    # Write table
    sobj._data.write(data_path('tmp.fits'), overwrite=True)
    tbl = Table.read(data_path('tmp.fits'))
    sobj2 = specobj.SpecObj.from_table(tbl)
    #
    assert isinstance(sobj2, specobj.SpecObj)
