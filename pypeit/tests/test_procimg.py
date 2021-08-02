"""
Module to run tests on core.procimg functions.
"""
from IPython import embed
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

def test_rn2_frame():
    # Bogus image
    datasec = np.ones((10,10), dtype=int)
    datasec[5:] = 2

    rn = np.array([2.5, 3.5])
    gain = np.array([1.2, 1.5])

    rnvar = procimg.rn2_frame(datasec, gain, rn, digitization=False)
    assert rnvar.shape == datasec.shape, 'Shape mismatch'
    assert np.array_equal(np.unique(rnvar), rn**2), 'Bad RN variance calculation'

    rnvar = procimg.rn2_frame(datasec, gain, rn, units='ADU', digitization=False)
    assert np.allclose(np.unique(rnvar), (rn/gain)**2), 'Bad RN variance calculation'


def test_sub_overscan():
    datasec = np.zeros((10,10), dtype=int)
    datasec[:5,:-3] = 1 
    datasec[5:,:-3] = 2 

    oscan = np.zeros((10,10), dtype=int)
    oscan[:5,-3:] = 1 
    oscan[5:,-3:] = 2 

    raw = np.zeros((10,10), dtype=float)
    raw[datasec == 1] = 10.
    raw[datasec == 2] = 20.
    raw[oscan == 1] = 9.
    raw[oscan == 2] = 19.

    raw_sub = procimg.subtract_overscan(raw, datasec, oscan, method='median')
    assert np.array_equal(raw_sub[datasec > 0], np.ones(np.sum(datasec > 0), dtype=float)), \
            'Bad overscan subtraction'

    var = np.ones((10,10), dtype=float)
    raw_sub, var_sub = procimg.subtract_overscan(raw, datasec, oscan, method='median', var=var)
    assert np.array_equal(var_sub[datasec > 0],
                          np.ones(np.sum(datasec > 0), dtype=float) + np.pi/2), \
            'Bad variance calculation'

def test_trim():
    datasec = np.zeros((10,10), dtype=int)
    datasec[:5,:-3] = 1 
    datasec[5:,:-3] = 2 

    _datasec = procimg.trim_frame(datasec, datasec < 1)
    assert _datasec.shape == (10,7), 'Trimming error'
    assert np.array_equal(datasec[datasec > 0], _datasec.flat), 'Values changed'




