"""
Module to run tests on core.procimg functions.
"""
import os

from IPython import embed

import numpy as np
import pytest
import scipy.signal

from astropy.convolution import convolve

from pypeit.core import procimg
from pypeit.pkg.exceptions import PypeItError
from pypeit import utils

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

    rnvar = procimg.rn2_frame(datasec, rn, digitization=False)
    assert rnvar.shape == datasec.shape, 'Shape mismatch'
    assert np.array_equal(np.unique(rnvar), rn**2), 'Bad RN variance calculation'

    rnvar = procimg.rn2_frame(datasec, rn, units='ADU', gain=gain, digitization=False)
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

    raw_sub, _ = procimg.subtract_overscan(raw, datasec, oscan, method='median')
    assert np.array_equal(raw_sub[datasec > 0], np.ones(np.sum(datasec > 0), dtype=float)), \
            'Bad overscan subtraction'

    var = np.ones((10,10), dtype=float)
    raw_sub, var_sub = procimg.subtract_overscan(raw, datasec, oscan, method='median', var=var)
    assert np.array_equal(var_sub[datasec > 0],
                          np.full(np.sum(datasec > 0), np.pi/2/15, dtype=float)), \
            'Bad variance calculation'


def _overscan_test_frame():
    """Build a two-amplifier test frame with row-dependent overscan levels."""
    datasec = np.zeros((10,10), dtype=int)
    datasec[:5,:-3] = 1
    datasec[5:,:-3] = 2

    oscan = np.zeros((10,10), dtype=int)
    oscan[:5,-3:] = 1
    oscan[5:,-3:] = 2

    level = np.arange(10, dtype=float)[:,None]
    raw = np.zeros((10,10), dtype=float)
    raw[datasec > 0] = (level + 1. + np.zeros((10,10)))[datasec > 0]
    raw[oscan > 0] = (level + np.zeros((10,10)))[oscan > 0]
    return raw, datasec, oscan


def test_overscan_stats_scalar():
    overscan = np.arange(12, dtype=float).reshape(4,3)
    osvar = np.full(overscan.shape, 2., dtype=float)
    os_bpm = np.zeros(overscan.shape, dtype=bool)

    # No masked pixels: identical to the unmasked calculation
    osfit, var = procimg.overscan_stats(overscan, os_bpm, osvar=osvar)
    assert osfit == np.median(overscan), 'Bad median with no masked pixels'
    assert var == np.pi/2*(np.sum(osvar)/osvar.size**2), 'Bad variance with no masked pixels'

    # Masked pixels are excluded from both the median and the variance
    os_bpm[3,:] = True
    osfit, var = procimg.overscan_stats(overscan, os_bpm, osvar=osvar)
    assert osfit == np.median(overscan[:3]), 'Masked pixels should be excluded from the median'
    assert np.isclose(var, np.pi/2*(2.*9/9**2)), \
            'Masked pixels should be excluded from the variance'

    # Variance is only computed if requested
    _, var = procimg.overscan_stats(overscan, os_bpm)
    assert var is None, 'Variance should be None if osvar is not provided'


def test_overscan_stats_vector():
    nrow, ncol = 6, 4
    overscan = np.repeat(np.arange(nrow, dtype=float)[:,None], ncol, axis=1)
    osvar = np.ones(overscan.shape, dtype=float)
    os_bpm = np.zeros(overscan.shape, dtype=bool)
    # One masked pixel in row 1, and all pixels masked in row 3
    overscan[1,0] = 1000.
    os_bpm[1,0] = True
    os_bpm[3,:] = True
    overscan[3,:] = -1000.

    osfit, var = procimg.overscan_stats(overscan, os_bpm, osvar=osvar, axis=1)
    # Row 3 is interpolated from rows 2 and 4
    assert np.array_equal(osfit, np.arange(nrow, dtype=float)), \
            'Masked pixels should be excluded and fully masked rows interpolated'
    expected_var = np.full(nrow, np.pi/2*(ncol/ncol**2))
    expected_var[1] = np.pi/2*((ncol-1)/(ncol-1)**2)
    assert np.allclose(var, expected_var), \
            'Variance should be computed per row from its good pixels'

    # Collapsing along the other axis gives one value per column
    osfit, _ = procimg.overscan_stats(overscan, os_bpm, axis=0)
    assert osfit.shape == (ncol,), 'Collapsing along axis 0 should give one value per column'


def test_overscan_stats_all_masked():
    overscan = np.ones((3,3), dtype=float)
    with pytest.raises(PypeItError):
        procimg.overscan_stats(overscan, np.ones(overscan.shape, dtype=bool))


def test_sub_overscan_undefined_median():
    raw, datasec, oscan = _overscan_test_frame()
    raw[datasec == 1] = 10.
    raw[datasec == 2] = 20.
    raw[oscan == 1] = 9.
    raw[oscan == 2] = 19.
    # A NaN in the overscan (e.g., a FITS BLANK pixel) should be ignored
    raw[2,-1] = np.nan
    # A pixel flagged by the bpm should be ignored, even if it is finite
    bpm = np.zeros(raw.shape, dtype=bool)
    raw[7,-1] = 1000.
    bpm[7,-1] = True

    var = np.ones((10,10), dtype=float)
    raw_sub, var_sub = procimg.subtract_overscan(raw, datasec, oscan, method='median', var=var,
                                                 bpm=bpm)
    assert np.array_equal(raw_sub[datasec > 0], np.ones(np.sum(datasec > 0), dtype=float)), \
            'Masked overscan pixels should not affect the overscan subtraction'
    # Each amplifier has 14 good overscan pixels
    assert np.allclose(var_sub[datasec > 0], np.pi/2/14), 'Bad variance calculation'


def test_sub_overscan_undefined_rows():
    raw, datasec, oscan = _overscan_test_frame()
    # Overscan-subtracted data should be all ones
    expected = np.ones(np.sum(datasec > 0), dtype=float)
    _raw = raw.copy()
    # One undefined overscan pixel
    _raw[1,-2] = np.nan
    # An entire overscan row that is undefined; its level is interpolated from
    # the adjacent rows
    _raw[6,-3:] = np.nan
    for method, params in zip(['savgol', 'chebyshev'], [[1, 3], [1]]):
        raw_sub, var_sub = procimg.subtract_overscan(_raw, datasec, oscan, method=method,
                                                     params=params,
                                                     var=np.ones(raw.shape, dtype=float))
        assert np.all(np.isfinite(raw_sub[datasec > 0])) \
                and np.all(np.isfinite(var_sub[datasec > 0])), \
                f'Undefined overscan pixels should not propagate ({method})'
        assert np.allclose(raw_sub[datasec > 0], expected), \
                f'Bad overscan subtraction with undefined pixels ({method})'


def test_sub_overscan_undefined_odd_even():
    raw, datasec, oscan = _overscan_test_frame()
    _raw = raw.copy()
    _raw[3,-1] = np.nan
    raw_sub, _ = procimg.subtract_overscan(_raw, datasec, oscan, method='odd_even')
    assert np.all(np.isfinite(raw_sub[datasec > 0])), \
            'Undefined overscan pixels should not propagate'
    assert np.allclose(raw_sub[datasec > 0], 1.), 'Bad overscan subtraction'


def test_sub_overscan_odd_even_variance():
    raw, datasec, oscan = _overscan_test_frame()
    var = np.ones(raw.shape, dtype=float)
    # Flag one pixel in the even stripe of the second row of amplifier 1
    bpm = np.zeros(raw.shape, dtype=bool)
    bpm[1,-3] = True
    _, var_sub = procimg.subtract_overscan(raw, datasec, oscan, method='odd_even', var=var,
                                           bpm=bpm)
    # Each row of amplifier 1 has 1 odd and 2 even overscan pixels (columns 7-9
    # of the frame are overscan stripe columns 0-2).  The variance in the median
    # of each row is pi/2 times the variance in the mean of the good pixels in
    # that row.
    odd_var = np.pi/2*(1./1**2)
    even_var = np.full(5, np.pi/2*(2./2**2))
    even_var[1] = np.pi/2*(1./1**2)
    amp_var = var_sub[:5,:-3]
    assert np.allclose(amp_var[:,1::2], odd_var), 'Bad odd-column variance'
    assert np.allclose(amp_var[:,0::2], even_var[:,None]), \
            'Bad even-column variance; must be computed per row from the good pixels'


def test_sub_overscan_unmasked_unchanged():
    # Results must be identical to the unmasked calculation if no overscan
    # pixels are flagged, even if pixels in the data region are flagged
    raw, datasec, oscan = _overscan_test_frame()
    rng = np.random.default_rng(seed=8)
    raw += rng.normal(scale=0.1, size=raw.shape)
    var = rng.uniform(low=0.5, high=1.5, size=raw.shape)
    bpm = np.zeros(raw.shape, dtype=bool)
    bpm[0,0] = True

    # Unmasked reference calculations for amplifier 1
    data = raw[:5,:-3]
    overscan = raw[:5,-3:]
    osvar = var[:5,-3:]
    ref = {'median': data - np.median(overscan),
           'savgol': data - scipy.signal.savgol_filter(np.median(overscan, axis=1), 3, 1)[:,None]}
    ref_var = {'median': np.full(data.shape, np.pi/2*(np.sum(osvar)/osvar.size**2)),
               'savgol': np.repeat(np.pi/2*(np.sum(osvar, axis=1)/osvar.shape[1]**2)[:,None],
                                   data.shape[1], axis=1)}
    for method in ['median', 'savgol']:
        test, test_var = procimg.subtract_overscan(raw, datasec, oscan, method=method,
                                                   params=[1, 3], var=var, bpm=bpm)
        assert np.array_equal(test[:5,:-3], ref[method]), \
                f'Unmasked overscan subtraction changed ({method})'
        assert np.array_equal(test_var[:5,:-3], ref_var[method]), \
                f'Unmasked overscan variance changed ({method})'

    # Odd/even stripes
    test, _ = procimg.subtract_overscan(raw, datasec, oscan, method='odd_even', bpm=bpm)
    expected = data.copy()
    expected[:,1::2] -= np.median(overscan[:,1::2], axis=1)[:,None]
    expected[:,0::2] -= np.median(overscan[:,0::2], axis=1)[:,None]
    assert np.array_equal(test[:5,:-3], expected), 'Unmasked odd/even subtraction changed'


def test_sub_overscan_all_undefined():
    raw, datasec, oscan = _overscan_test_frame()
    raw[oscan == 1] = np.nan
    with pytest.raises(PypeItError):
        procimg.subtract_overscan(raw, datasec, oscan, method='median')


def test_trim():
    datasec = np.zeros((10,10), dtype=int)
    datasec[:5,:-3] = 1 
    datasec[5:,:-3] = 2 

    _datasec = procimg.trim_frame(datasec, datasec < 1)
    assert _datasec.shape == (10,7), 'Trimming error'
    assert np.array_equal(datasec[datasec > 0], _datasec.flat), 'Values changed'


def test_var_model():
    # Bogus image
    datasec = np.ones((10,10), dtype=int)
    datasec[5:] = 2

    rn = np.array([2.5, 3.5])

    rnvar = procimg.rn2_frame(datasec, rn)

    assert np.array_equal(rnvar, procimg.base_variance(rnvar)), \
        'Variance model with only rnvar is just rnvar'

    counts = np.full(rnvar.shape, 10., dtype=float)

    assert np.array_equal(rnvar, procimg.variance_model(rnvar)), \
        'Variance model should just return input if no optional parameters are provided.'

    base = procimg.base_variance(rnvar, darkcurr=10.)
    base_t = procimg.base_variance(rnvar, darkcurr=5., exptime=2.*3600)

    assert np.all(procimg.variance_model(rnvar, counts=counts) > rnvar), \
        'Shot noise should increase the variance'
    assert np.all(procimg.variance_model(base, counts=counts) > base), \
        'Shot noise should increase the variance'
    assert np.array_equal(
                procimg.variance_model(base, counts=counts),
                procimg.variance_model(base_t, counts=counts)), \
        'Dark current should be equivalent'
    assert np.all(procimg.base_variance(rnvar, proc_var=10.) > rnvar), \
        'Processing variance should increase the total variance'

    assert np.all(procimg.variance_model(rnvar, counts=counts, count_scale=0.5) <
                  procimg.variance_model(rnvar, counts=counts)), \
        'Scaling should have decreased the noise.'

    assert np.all(procimg.variance_model(rnvar, counts=counts, noise_floor=0.1) > rnvar), \
        'Noise floor should have increased the variance.'


def test_grow_mask():
    mask = np.zeros((9,9), dtype=bool)
    mask[4,4] = True
    grw_msk = procimg.grow_mask(mask, 2.)
    _grw_msk = np.zeros((9,9), dtype=bool)
    _grw_msk[3:-3,3] = True
    _grw_msk[2:-2,4] = True
    _grw_msk[3:-3,5] = True
    _grw_msk[3,3:-3] = True
    _grw_msk[4,2:-2] = True
    _grw_msk[5,3:-3] = True
    assert np.array_equal(grw_msk, _grw_msk), 'Bad mask growth'


def test_boxcar():
    a = np.arange(100).reshape(10,10).astype(float)
    arep = procimg.boxcar_replicate(a, 2)
    assert np.array_equal(procimg.boxcar_average(arep, 2), a), 'Bad replicate/average'
    assert np.array_equal(utils.rebinND(arep, a.shape), a), 'Bad replicate/rebin'


def test_lacosmic_spatial_median_residual():
    """Diagonal CR trail must be flagged without harming horizontal arc lines."""
    rng = np.random.default_rng(7)
    image = 10.0 + rng.normal(scale=1.0, size=(120, 160))
    image[35, :] += 500.0
    image[85, :] += 300.0
    xs = np.arange(20, 135)
    ys = np.rint(15 + 0.48 * xs).astype(int)
    image[ys, xs] += 1000.0

    crmask = procimg.lacosmic_spatial_median_residual(
        image, median_width=51,
        sigclip=10.0, sigfrac=0.3, objlim=0.0,
        remove_compact_obj=False, maxiter=2, grow=2.0,
    )

    # Most of the diagonal trail must be flagged.
    assert crmask[ys, xs].mean() > 0.5, 'Diagonal CR trail not detected.'

    # The bulk of each bright arc row must remain unflagged -- the point of
    # subtracting the row-local median is to leave arc lines intact.
    for row in (35, 85):
        unflagged_frac = float((~crmask[row]).mean())
        assert unflagged_frac > 0.9, \
            f'Arc row {row} over-flagged: only {unflagged_frac:.1%} unflagged'


def test_lacosmic_spatial_median_residual_disabled():
    """median_width <= 1 returns an empty mask."""
    rng = np.random.default_rng(0)
    image = rng.normal(size=(20, 30))
    crmask = procimg.lacosmic_spatial_median_residual(image, median_width=1)
    assert not np.any(crmask)
    assert crmask.shape == image.shape

