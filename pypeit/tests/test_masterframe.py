"""
Module to run tests on armasters
"""
import os
import numpy as np
import pytest

from pypeit import masterframe

def data_root():
    return os.path.join(os.path.dirname(__file__), 'files')

# Could break this into a bunch of separate tests...
def test_master_io():
    """
    Test the methods in MasterFrame
    """
    # Init
    master_dir = data_root()
    mf = masterframe.MasterFrame('Test', master_dir=master_dir, master_key='A_1_01',
                                 reuse_masters=True)
    # Check the filename
    assert mf.file_name == 'MasterTest_A_1_01.fits', 'Incorrect master frame file name'

    # Generate a header
    hdr = mf.build_master_header(steps=['one', 'two', 'three'])
    assert hdr['STEPS'] == 'one,two,three'
    for key in ['MSTRTYP', 'MSTRDIR', 'MSTRKEY']:
        assert key in hdr.keys()

