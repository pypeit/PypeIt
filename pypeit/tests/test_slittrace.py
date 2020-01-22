"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.slittrace import SlitTraceSet

def test_init():

    slits = SlitTraceSet(left=np.full((1000,3), 2, dtype=float),
                         right=np.full((1000,3), 8, dtype=float), nspat=10, spectrograph='dummy',
                         master_key='dummy', master_dir=os.getcwd())

    assert np.all(slits.center == 5), 'Bad center'
    assert slits.file_name == 'MasterSlits_dummy.fits.gz', 'Bad master name'


def test_io():

    slits = SlitTraceSet(left=np.full((1000,3), 2, dtype=float),
                         right=np.full((1000,3), 8, dtype=float), nspat=10, spectrograph='dummy',
                         master_key='dummy', master_dir=os.getcwd())
    
    if slits.exists:
        os.remove(slits.master_file_path)

    slits.save()
    assert slits.exists, 'File not written'

    slits.load()
    assert np.array_equal(slits.left, np.full((1000,3), 2, dtype=float)), 'Bad left read'

    other_ofile = 'test.fits.gz'
    if os.path.isfile(other_ofile):
        os.remove(other_ofile)

    # Test write
    slits.save(ofile=other_ofile)

    # Test overwrite
    slits.save(ofile=other_ofile, overwrite=True)

    # Test from_file
    _slits = SlitTraceSet.from_file(other_ofile)
    assert np.array_equal(slits.right, _slits.right), 'Bad read from_file'
    
    # Clean up
    os.remove(slits.master_file_path)
    os.remove(other_ofile)
    

