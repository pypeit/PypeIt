
import time
import pytest

from IPython import embed

import numpy as np

from pypeit.tests.tstutils import bspline_ext_required, data_path

@bspline_ext_required
def test_cholesky_band_versions():
    # Import only when the test is performed
    from pypeit.bspline.utilpy import cholesky_band as cholesky_band_py
    from pypeit.bspline.utilc import cholesky_band as cholesky_band_c

    # Read data
    d = np.load(data_path('cholesky_band_l.npz'))

    # Run python version
    pytime = time.perf_counter()
    e, l = cholesky_band_py(d['l'], mininf=d['mininf'])
    pytime = time.perf_counter() - pytime

    # Run C version
    ctime = time.perf_counter()
    e, _l = cholesky_band_c(d['l'], mininf=d['mininf'])
    ctime = time.perf_counter() - ctime

    # Check time and output arrays
    assert ctime < pytime, 'C is less efficient!'
    assert np.allclose(l, _l), 'Differences in cholesky_band'


@bspline_ext_required
def test_cholesky_solve_versions():
    # Import only when the test is performed
    from pypeit.bspline.utilpy import cholesky_solve as cholesky_solve_py
    from pypeit.bspline.utilc import cholesky_solve as cholesky_solve_c

    # Read data
    d = np.load(data_path('cholesky_solve_abb.npz'))

    # Run python version
    pytime = time.perf_counter()
    e, b = cholesky_solve_py(d['a'], d['bb'])
    pytime = time.perf_counter() - pytime

    # Run C version
    ctime = time.perf_counter()
    t = time.perf_counter()
    e, _b = cholesky_solve_c(d['a'], d['bb'])
    ctime = time.perf_counter() - ctime

    # Check time and output arrays
    assert ctime < pytime, 'C is less efficient!'
    assert np.allclose(b, _b), 'Differences in cholesky_solve'


