"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.slittrace import SlitTraceSet, SlitTraceBitMask
from pypeit import masterframe

master_key = 'dummy'
master_dir = os.getcwd()

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_bits():
    # Make sure bits are correct
    bm = SlitTraceBitMask()
    assert bm.bits['USERIGNORE'] == 2, 'Bits changed'
    assert bm.bits['BADFLATCALIB'] == 6, 'Bits changed'

def test_init():

    slits = SlitTraceSet(left_init=np.full((1000,3), 2, dtype=float),
                         right_init=np.full((1000,3), 8, dtype=float),
                         pypeline='MultiSlit',
                         nspat=10, PYP_SPEC='dummy')

    left, right, _ = slits.select_edges()
    center = (left+right)/2
    assert np.all(center == 5), 'Bad center'

def test_io():

    slits = SlitTraceSet(np.full((1000,3), 2, dtype=float), np.full((1000,3), 8, dtype=float),
                         'MultiSlit',
                         nspat=10, PYP_SPEC='dummy')
    master_file = masterframe.construct_file_name(slits, master_key, master_dir=master_dir)

    # Remove any existing file from previous runs that were interrupted
    if os.path.isfile(master_file):
        os.remove(master_file)

    # Try to save it
    slits.to_master_file(master_file) #master_dir, master_key,  'dummy_spectrograph_name')
    assert os.path.isfile(master_file), 'File not written'

    # Instantiate an empty SlitTraceSet with the same master file, and
    # indicate it should be reused
    #_slits = SlitTraceSet.from_master('dummy', os.getcwd())
    _slits = SlitTraceSet.from_file(master_file)
    assert np.array_equal(_slits.left_init, np.full((1000,3), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

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
    assert np.array_equal(slits.right_init, _slits.right_init), 'Bad read from_file'

    # Clean up
    os.remove(master_file)
    os.remove(other_ofile)


def test_io_single():
    # NOTE this is just a test string. The file itself is not actually read, so it is not required to run the test.
    file = '/home/xavier/Projects/PypeIt-development-suite/RAW_DATA/keck_deimos/830G_M_8500/DE.20100913.57006.fits.gz'
    slits = SlitTraceSet(np.full((1000, 1), 2, dtype=float), np.full((1000, 1), 8, dtype=float),
                         'MultiSlit',
                         nspat=10, PYP_SPEC='dummy',
                         maskfile=file)

    # Remove any existing file from previous runs that were interrupted
    tst_file = data_path('tst_slittrace.fits')
    if os.path.isfile(tst_file):
        os.remove(tst_file)

    # Try to save it
    slits.to_file(tst_file)

    _slits = SlitTraceSet.from_file(tst_file)

    assert np.array_equal(_slits.left_init, np.full((1000, 1), 2, dtype=float)), 'Bad left read'
    # And that it's the same as the existing one
    assert np.array_equal(_slits.left_init, slits.left_init), 'Bad left read'

    # Try to read/write to a custom file name
    # Remove existing file from previous runs that were interrupted
    if os.path.isfile(tst_file):
        os.remove(tst_file)

    # Test write
    slits.to_file(tst_file)

    # Test overwrite
    slits.to_file(tst_file, overwrite=True)

    # Test from_file
    _slits = SlitTraceSet.from_file(tst_file)
    assert np.array_equal(slits.right_init, _slits.right_init), 'Bad read from_file'

    # Clean up
    os.remove(tst_file)


