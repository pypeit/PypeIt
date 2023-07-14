"""
Module to run tests on arparse
"""
import numpy as np
import pytest

from pypeit.core import parse

def test_parse_binning():
    """ Test parse binning algorithm
    """
    bin1, bin2 = parse.parse_binning('2,2')
    assert bin1 == 2
    assert bin2 == 2
    # Other input
    bin1, bin2 = parse.parse_binning((2,2))
    assert bin1 == 2
    assert bin2 == 2


def test_sec2slice():
    sub = ':10,10:'
    subslice = parse.sec2slice(sub, require_dim=2)
    assert subslice[0].start is None

    subslice = parse.sec2slice(sub, include_end=True, require_dim=2)
    assert subslice[0].stop == 11

    subslice = parse.sec2slice(sub, one_indexed=True, require_dim=2)
    assert subslice[0].stop == 9


def test_str2list():
    """
    Test the conversion of the string to a list of integers
    """
    assert np.array_equal(parse.str2list('all', length=10), [0,1,2,3,4,5,6,7,8,9])
    assert np.array_equal(parse.str2list(':4', length=10), [0,1,2,3])
    assert np.array_equal(parse.str2list('3:5,8:', length=10), [3,4,8,9])
    assert np.array_equal(parse.str2list('3,1:5,6', length=10), [1,2,3,4,6])
    assert np.array_equal(parse.str2list('3,1:5,8:', length=10), [1,2,3,4,8,9])

    
