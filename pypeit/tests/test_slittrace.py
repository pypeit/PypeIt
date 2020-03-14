"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.slittrace import SlitTraceSet

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():

    slits = SlitTraceSet(left=np.full((1000,3), 2, dtype=float),
                         right=np.full((1000,3), 8, dtype=float), nspat=10, spectrograph='dummy',
                         master_key='dummy', master_dir=os.getcwd())

    assert np.all(slits.center == 5), 'Bad center'
    assert slits.file_name == 'MasterSlits_dummy.fits.gz', 'Bad master name'


def test_io():

    slits = SlitTraceSet(np.full((1000,3), 2, dtype=float), np.full((1000,3), 8, dtype=float),
                         nspat=10, spectrograph='dummy', master_key='dummy',
                         master_dir=os.getcwd())

    # Remove any existing file from previous runs that were interrupted
    if slits.exists:
        os.remove(slits.master_file_path)

    # Try to save it
    slits.to_master()
    assert slits.exists, 'File not written'

    # Instantiation default is to not reuse existing master frames
    assert slits.chk_load_master(None) is None, 'Should not try to load the master'

    # Instantiate an empty SlitTraceSet with the same master file, and
    # indicate it should be reused
    assert SlitTraceSet.from_master('dummy', os.getcwd(), reuse=False) is None, \
            'Should not load master'
    _slits = SlitTraceSet.from_master('dummy', os.getcwd())
    assert np.array_equal(_slits.left, np.full((1000,3), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left, slits.left), 'Bad left read'

    # Try to read/write to a custom file name
    # Remove existing file from previous runs that were interrupted
    other_ofile = 'test.fits.gz'
    if os.path.isfile(other_ofile):
        os.remove(other_ofile)

    # Test write
    slits.to_file(other_ofile)

    # Test overwrite
    slits.to_file(other_ofile, overwrite=True)

    # Test from_file
    _slits = SlitTraceSet.from_file(other_ofile)
    assert np.array_equal(slits.right, _slits.right), 'Bad read from_file'

    # Clean up
    os.remove(slits.master_file_path)
    os.remove(other_ofile)


def test_io_single():
    slits = SlitTraceSet(np.full((1000, 1), 2, dtype=float), np.full((1000, 1), 8, dtype=float),
                         nspat=10, spectrograph='dummy', master_key='dummy',
                         master_dir=os.getcwd())

    # Remove any existing file from previous runs that were interrupted

    # Try to save it
    slits.to_file(data_path('tst_slittrace.fits'))

    _slits = SlitTraceSet.from_file(data_path('tst_slittrace.fits'))

    assert np.array_equal(_slits.left, np.full((1000, 3), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left, slits.left), 'Bad left read'

    # Try to read/write to a custom file name
    # Remove existing file from previous runs that were interrupted
    other_ofile = 'test.fits.gz'
    if os.path.isfile(other_ofile):
        os.remove(other_ofile)

    # Test write
    slits.to_file(other_ofile)

    # Test overwrite
    slits.to_file(other_ofile, overwrite=True)

    # Test from_file
    _slits = SlitTraceSet.from_file(other_ofile)
    assert np.array_equal(slits.right, _slits.right), 'Bad read from_file'

    # Clean up
    os.remove(slits.master_file_path)
    os.remove(other_ofile)


