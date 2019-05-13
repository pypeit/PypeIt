"""
Module to run tests on core.procimg functions.
"""
import pytest
import numpy as np

from pypeit.core import procimg

def test_replace_columns():
    y = np.zeros((10,3), dtype=float)
    y[:,2] = 2
    bad_col = np.array([False, True, False])
    assert np.array_equal(procimg.replace_columns(y, bad_col, copy=True, replace_with='mean'),
                          procimg.replace_columns(y, bad_col, copy=True, replace_with='linear')), \
                'Interpolation and mean should provide the same result.'

    bad_col = np.array([False, True, True])
    assert np.array_equal(procimg.replace_columns(y, bad_col, copy=True, replace_with='mean'),
                          np.zeros_like(y)), 'Should set everything to 0.'

    bad_col = np.array([True, True, False])
    assert np.array_equal(procimg.replace_columns(y, bad_col, copy=True, replace_with='mean'),
                          np.full_like(y, 2)), 'Should set everything to 2.'

    y = np.zeros((10,4), dtype=float)
    y[:,3] = 3
    bad_col = np.array([False, True, True, False])
    assert np.array_equal(procimg.replace_columns(y, bad_col, copy=True, replace_with='linear'),
                          np.repeat(np.arange(4),10).reshape(4,10).T), \
                'Interpolation failed.'

