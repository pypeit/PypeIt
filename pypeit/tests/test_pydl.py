from IPython import embed

import numpy as np
import pytest

from pypeit.core import pydl


# NOTE: Most of the tests in this module were written with the help of AI.

def make_test_arrays(n=1000, seed=99):
    """
    Create test arrays for pydl tests using fixed size and seed.

    Parameters
    ----------
    n : int, optional
        Length of the data array to create.
    seed : int, `np.random.Generator`, optional
        RNG seed for reproducibility.

    Returns
    -------
    data : np.ndarray
        Normal-distributed data with unity standard deviation and length 1000.
    invvar : np.ndarray
        Inverse variance array of ones.
    mask : np.ndarray
        Boolean mask array set to True everywhere.
    """
    rng = np.random.default_rng(seed)
    data = rng.normal(loc=0.0, scale=1.0, size=1000)
    invvar = np.ones_like(data)
    mask = np.ones_like(data, dtype=bool)
    return data, invvar, mask


def add_outliers(data, n_outliers=100, sigma=10.0, seed=99):
    """Return a copy of `data` with `n_outliers` randomly replaced by
    samples drawn from a normal distribution with standard deviation
    `sigma`.

    Parameters
    ----------
    data : `numpy.ndarray`
        Input 1-D data array.
    n_outliers : int, optional
        Number of values to replace.
    sigma : float, optional
        Standard deviation of the outlier normal distribution.
    seed : int, `np.random.Generator`, optional
        RNG seed for reproducibility.

    Returns
    -------
    new_data : `numpy.ndarray`
        Copy of input data with outliers inserted.
    bpm : `numpy.ndarray`
        Bad-pixel mask used to select pixels with outliers in ``new_data``.
    """
    rng = np.random.default_rng(seed)
    n = data.size
    if n_outliers >= n:
        raise ValueError('n_outliers must be less than data length')
    indx = rng.choice(n, size=n_outliers, replace=False)
    new_data = data.copy()
    new_data[indx] = rng.normal(loc=0.0, scale=sigma, size=n_outliers)
    bpm = np.zeros(n, dtype=bool)
    bpm[indx] = True
    return new_data, bpm


def test_no_rejection_when_equal():
    data, invvar, mask = make_test_arrays()
    model = data.copy()
    outmask_in = np.ones_like(data, dtype=bool)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=outmask_in.copy(), invvar=invvar, lower=3., upper=3.)
    assert qdone, 'Expected qdone True when data equals model'
    assert np.all(outmask), 'Expected all points to be kept when data equals model'


def test_detect_injected_outliers():
    rng = np.random.default_rng(99)
    data, invvar, mask = make_test_arrays(seed=rng)
    data, injected = add_outliers(data, n_outliers=100, sigma=10.0, seed=rng)
    model = np.zeros_like(data)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar, lower=5., upper=5.
    )
    # At least one injected outlier should be detected as rejected
    detected = np.logical_not(outmask) & injected
    assert np.sum(detected) == 63, 'Should have rejected 63 injected outliers'
    assert np.sum(np.logical_not(outmask) & np.logical_not(injected)) == 0, \
        'Should not reject any values that were not injected'


def test_lower_upper_rejection():
    data, invvar, mask = make_test_arrays()
    model = np.zeros_like(data)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar, lower=1.
    )
    assert np.all(data[np.logical_not(outmask)] < 0.0), 'Should only reject lower outliers'
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar, upper=1.
    )
    assert np.all(data[np.logical_not(outmask)] > 0.0), 'Should only reject upper outliers'


def test_use_mad():
    rng = np.random.default_rng(99)
    data, invvar, mask = make_test_arrays(seed=rng)
    data, injected = add_outliers(data, n_outliers=100, sigma=10.0, seed=rng)
    model = np.zeros_like(data)

    with pytest.raises(ValueError):
        # Cannot use MAD when providing the inverse variance
        pydl.djs_reject(data, model, invvar=invvar, use_mad=True)

    outmask, qdone = pydl.djs_reject(data, model, use_mad=True, upper=5.0, lower=5.0)
    assert np.sum(np.logical_not(outmask)) == 60, 'Should have rejected 60 outliers'


def test_grow_and_sticky():
    rng = np.random.default_rng(99)
    data, invvar, mask = make_test_arrays(seed=rng)
    # inject an outlier in the middle
    mid = data.size // 2
    data[mid] = 10.0
    model = np.zeros_like(data)
    outmask_in = np.ones_like(data, dtype=bool)

    # Reject the outlier and grow by 1 should also reject its neighbors
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=outmask_in.copy(), invvar=invvar, upper=3, grow=1
    )
    assert not np.any(outmask[mid-1:mid+2]), 'Should reject outlier and its immediate neighbors'

    # Now test sticky
    # Reject some additional points at random
    newrej = rng.choice(data.size, size=10, replace=False)
    outmask[newrej] = False
    _outmask, qdone = pydl.djs_reject(
        data, model, outmask=outmask, invvar=invvar, upper=3, sticky=True
    )
    # The sticky False at index 0 must remain
    assert not np.any(outmask[newrej]), \
        'Expected previously rejected pixels to remain rejected when sticky=True'


def test_outmask_none_initialization():
    rng = np.random.default_rng(99)
    data, invvar, mask = make_test_arrays(seed=rng)
    data, injected = add_outliers(data, n_outliers=100, sigma=10.0, seed=rng)
    model = np.zeros_like(data)
    outmask, qdone = pydl.djs_reject(data, model, outmask=None, invvar=invvar, upper=5)
    assert outmask.size == data.size, 'Outmask has the wrong size'
    assert np.any(np.logical_not(outmask)), 'Expected some points to be rejected'


def test_percentile_flag_rejection():
    rng = np.random.default_rng(99)
    data, invvar, mask = make_test_arrays(seed=rng)
    data, injected = add_outliers(data, n_outliers=100, sigma=10.0, seed=rng)
    model = np.zeros_like(data)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar, percentile=True,
        upper=95
    )
    detected = np.logical_not(outmask) & injected
    assert np.sum(detected) > 0, \
        'Expected percentile=True to reject some injected high-value outliers'


def test_maxdev_rejection():
    data, invvar, mask = make_test_arrays()
    data[5] = 100.0
    model = np.zeros_like(data)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar, maxdev=10.0
    )
    assert not outmask[5], 'Expected large absolute deviation to be rejected by maxdev'


def test_inmask_propagation():
    data, invvar, mask = make_test_arrays()
    data[0] = 50.0
    model = np.zeros_like(data)
    inmask = np.ones_like(data, dtype=bool)
    inmask[0] = False
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), inmask=inmask, invvar=invvar
    )
    assert not outmask[0], \
        'Expected inmask=False to propagate into outmask and mark pixel rejected'


# NOTE: Used to confirm consistency between old and new versions.  Test is now
# obsolete, but I've kept it here as a historical record.
#def test_djs_laxisgen_numpy_matches_original():
#    # Compare 1D, 2D, 3D outputs between original and numpy implementations
#    cases = {
#        1: (5,),
#        2: (4, 3),
#        3: (3, 4, 2),
#    }
#    for nd, dims in cases.items():
#        for iaxis in range(nd):
#            a1 = pydl.djs_laxisgen(list(dims), iaxis=iaxis)
#            a2 = pydl.arange_ndim(list(dims), axis=iaxis)
#            assert a1.shape == a2.shape, (
#                f'Shape mismatch for dims={dims} iaxis={iaxis}'
#            )
#            assert np.array_equal(a1, a2), (
#                f'djs_laxisgen and arange_ndim differ for dims={dims} '
#                f'iaxis={iaxis}'
#            )


def test_arange_ndim():
    n = 4
    arr = pydl.arange_ndim(n)
    assert np.array_equal(arr, np.arange(n, dtype=int)), '1D output should be the same as arange'
    dims = (4,)
    _arr = pydl.arange_ndim(dims)
    assert np.array_equal(arr, _arr), 'Output for an int and a 1-tuple should be the same'

    dims = (4,3,2)
    for axis in range(len(dims)):
        arr = pydl.arange_ndim(dims, axis=axis)
        for i in np.arange(dims[axis]):
            slc = [slice(None)] * len(dims)
            slc[axis] = i
            assert np.all(arr[tuple(slc)] == i), (
                 f'All elements at axis={axis} index={i} should be equal to {i}'
            )


def test_maxdev_lower_upper_interplay():
    # Use a small array for clarity
    data, invvar, mask = make_test_arrays(n=20, seed=99)
    # Create a large positive and a large negative deviation
    data[:] = 0.0
    data[0] = 20.0
    data[1] = -20.0
    model = np.zeros_like(data)

    # maxdev alone should reject the large positive deviation
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar,
        maxdev=10.0
    )
    assert not outmask[0], 'maxdev should reject a value with abs(diff)>maxdev'

    # maxdev + high upper: upper does not reject, maxdev still does
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar,
        upper=30.0, maxdev=10.0
    )
    assert not outmask[0], (
        'Expected rejection when maxdev triggers even if upper threshold is loose'
    )

    # upper low enough to reject even if maxdev would not (reverse interplay)
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar,
        upper=5.0, maxdev=25.0
    )
    assert not outmask[0], 'upper threshold should reject the positive outlier here'

    # lower alone should reject the large negative deviation
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar,
        lower=15.0
    )
    assert not outmask[1], 'lower should reject a value with chi < -lower'

    # lower too large to trigger but maxdev small triggers instead
    outmask, qdone = pydl.djs_reject(
        data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar,
        lower=25.0, maxdev=10.0
    )
    assert not outmask[1], (
        'Expected maxdev to reject the negative outlier when lower threshold is not met'
    )


def test_model_shape_mismatch_raises():
    data, invvar, mask = make_test_arrays()
    # model of wrong size should raise
    model = np.zeros(data.size - 10)
    with pytest.raises(ValueError):
        pydl.djs_reject(data, model, outmask=np.ones_like(data, dtype=bool), invvar=invvar)


def test_outmask_shape_mismatch_raises():
    data, invvar, mask = make_test_arrays()
    model = np.zeros_like(data)
    bad_outmask = np.ones(data.size - 5, dtype=bool)
    with pytest.raises(ValueError):
        pydl.djs_reject(data, model, outmask=bad_outmask, invvar=invvar)


def test_inmask_shape_mismatch_raises():
    data, invvar, mask = make_test_arrays()
    model = np.zeros_like(data)
    bad_inmask = np.ones(data.size - 2, dtype=bool)
    with pytest.raises(ValueError):
        pydl.djs_reject(
            data, model, outmask=np.ones_like(data, dtype=bool), inmask=bad_inmask, invvar=invvar
        )
