"""
Module to run tests on fitting
"""
import os

from IPython import embed

import pytest

import numpy as np

from pypeit import dataPaths
from pypeit import PypeItError
from pypeit.core import fitting
from pypeit.tests.tstutils import data_output_path


def test_pypeitfit():
    out_file = data_output_path('test_fit.fits')
    if os.path.isfile(out_file):
        os.remove(out_file)
    pypeitFit = fitting.PypeItFit(fitc=np.arange(100).astype(float))
    # Write
    pypeitFit.to_file(out_file)
    # Read
    pypeitFit2 = fitting.PypeItFit.from_file(out_file)
    assert np.array_equal(pypeitFit.fitc, pypeitFit2.fitc)
    pypeitFit2.to_file(out_file, overwrite=True)
    # Finish
    os.remove(out_file)


def test_polynomial():
    """ Run the parameter setup script
    """
    x = np.pi*np.linspace(0, 1., 100)
    y = np.sin(x)
    # Polynomial
    pypeitFit = fitting.PypeItFit(xval=x, yval=y, func='polynomial', order=np.array([3]))
    pypeitFit.fit()
    #pypeitFit = fitting.func_fit(x, y, 'polynomial', 3)
    np.testing.assert_allclose(pypeitFit.fitc, np.array([ -4.74660344e-02,   1.30745471e+00,
                                                  -4.16175760e-01, 3.08557167e-18]), atol=1e-9)
    # Evaluate
    val = pypeitFit.eval(x)
    assert np.isclose(val[0], -0.04746603), 'Bad value'


def test_legendre():
    x = np.pi*np.linspace(0, 1., 100)
    y = np.sin(x)
    # Legendre
    pypeitFit = fitting.PypeItFit(xval=x, yval=y, func='legendre', order=np.array([3]))
    pypeitFit.fit() # = fitting.func_fit(x, y, 'legendre', 3)
    np.testing.assert_allclose(pypeitFit.fitc, np.array([  6.37115652e-01,   6.83317251e-17,
                                                   -6.84581686e-01, -7.59352737e-17]), atol=1e-9)


@pytest.mark.remote_data
def test_pypeitfitcollection():

    # These are results generated from the old TraceSet fitter.  This just
    # checks that the results are the same.
    test_data_file = dataPaths.tests.get_file_path('trace_fit_data.npz')
    db = np.load(test_data_file)

    trace_coo = db['trace_coo']
    cen = db['cen']
    trace_fit_ivar = db['trace_fit_ivar']
    xmin = float(db['xmin'])
    xmax = float(db['xmax'])
    function = str(db['function'])
    order = int(db['order'])
    maxdev = float(db['maxdev'])
    maxiter = int(db['maxiter'])
    trace_fit = db['trace_fit']

    fit = fitting.PypeItFitCollection(
        trace_coo, cen, ivar=trace_fit_ivar, func=function, order=order, xmin=xmin, xmax=xmax,
        maxdev=maxdev, maxiter=maxiter
    )

    # NOTE: The array_equal test passed on my machine, but I changed it to
    # allclose, just to be on the safe side.
#    assert np.array_equal(fit.yfit, trace_fit), 'Fit changed!'
    assert np.allclose(fit.yfit, trace_fit), 'Fit changed!'


def _polyfit2d_both_paths(x, y, z, deg, **kwargs):
    # Solve once with the full design matrix and once in chunks
    full_chunk = fitting.POLYFIT2D_CHUNK
    try:
        fitting.set_polyfit2d_chunk(x.size)
        c_full = fitting.polyfit2d_general(x, y, z, deg, **kwargs)[0]
        fitting.set_polyfit2d_chunk(1000)
        c_chunk = fitting.polyfit2d_general(x, y, z, deg, **kwargs)[0]
    finally:
        fitting.set_polyfit2d_chunk(full_chunk)
    return c_full, c_chunk


def test_set_polyfit2d_chunk():
    full_chunk = fitting.POLYFIT2D_CHUNK
    try:
        fitting.set_polyfit2d_chunk(1000)
        assert fitting.POLYFIT2D_CHUNK == 1000, 'Chunk size not set'
        with pytest.raises(PypeItError):
            fitting.set_polyfit2d_chunk(0)
        assert fitting.POLYFIT2D_CHUNK == 1000, 'Invalid chunk size changed the setting'
    finally:
        fitting.set_polyfit2d_chunk(full_chunk)


@pytest.mark.parametrize('function, deg', [('legendre', (5, 6)), ('chebyshev', (3, 4)),
                                           ('polynomial', (3, 3))])
def test_polyfit2d_chunked(function, deg):
    rng = np.random.default_rng(42)
    npts = 10500     # Deliberately not a multiple of the chunk size
    x = rng.uniform(0, 4096, npts)
    y = rng.uniform(0, 4112, npts)
    z = np.sin(x/700) + 0.3*np.cos(y/900) + 0.01*rng.standard_normal(npts)
    w = rng.uniform(0.5, 2.0, npts)
    c_full, c_chunk = _polyfit2d_both_paths(x, y, z, deg, w=w, function=function)
    assert c_chunk.shape == c_full.shape, 'Chunked fit has the wrong shape'
    if function == 'polynomial':
        # The unscaled power basis is so poorly conditioned (cond ~ 1e23, with
        # half the singular values below the lstsq cutoff) that the
        # coefficients are only determined up to BLAS-dependent round-off, so
        # compare the fitted models instead.
        m_full = np.polynomial.polynomial.polyval2d(x, y, c_full)
        m_chunk = np.polynomial.polynomial.polyval2d(x, y, c_chunk)
        assert np.allclose(m_chunk, m_full, rtol=0, atol=1e-8), \
            'Chunked fit differs from the full fit'
    else:
        assert np.allclose(c_chunk, c_full, rtol=0, atol=1e-10*np.max(np.abs(c_full))), \
            'Chunked fit differs from the full fit'


def test_polyfit2d_chunked_rank_deficient():
    # With y constant, only the x terms are constrained; both paths must return
    # the same minimum-norm solution.
    rng = np.random.default_rng(42)
    npts = 10500
    x = rng.uniform(0, 4096, npts)
    y = np.full(npts, 7.0)
    z = np.sin(x/700) + 0.01*rng.standard_normal(npts)
    c_full, c_chunk = _polyfit2d_both_paths(x, y, z, (3, 2), function='legendre', miny=0., maxy=10.)
    assert np.allclose(c_chunk, c_full, rtol=0, atol=1e-10), \
        'Chunked fit differs from the full fit for a rank-deficient system'


@pytest.mark.parametrize('function', ['legendre2d', 'chebyshev2d', 'polynomial2d'])
def test_evaluate_fit_2d_chunked(function):
    rng = np.random.default_rng(42)
    c = rng.standard_normal((4, 5))
    val2d = {'legendre2d': np.polynomial.legendre.legval2d,
             'chebyshev2d': np.polynomial.chebyshev.chebval2d,
             'polynomial2d': np.polynomial.polynomial.polyval2d}[function]
    # 2D input, so the output shape must be preserved, and not a multiple of
    # the chunk size
    x, y = np.meshgrid(np.linspace(-1, 1, 101), np.linspace(-1, 1, 99))
    full = val2d(x, y, c)
    for chunk in (1000, 7, x.size, 10*x.size):
        chunked = fitting._chunked_val2d(val2d, x, y, c, chunk=chunk)
        assert chunked.shape == full.shape, 'Chunked evaluation changed the output shape'
        assert np.array_equal(chunked, full), 'Chunked evaluation changed the result'
    # And through evaluate_fit, with the [-1,1] scaling
    assert np.array_equal(
        fitting.evaluate_fit(c, function, x, x2=y, minx=-1., maxx=1., minx2=-1., maxx2=1.),
        fitting.evaluate_fit(c, function, x.ravel(), x2=y.ravel(), minx=-1., maxx=1., minx2=-1.,
                             maxx2=1.).reshape(x.shape))


def test_evaluate_fit_2d_chunked_dtype():
    # Chunked evaluation must return the same type as numpy, e.g. float32
    rng = np.random.default_rng(42)
    c = rng.standard_normal((4, 5)).astype(np.float32)
    x, y = np.meshgrid(np.linspace(-1, 1, 101, dtype=np.float32),
                       np.linspace(-1, 1, 99, dtype=np.float32))
    full = np.polynomial.legendre.legval2d(x, y, c)
    chunked = fitting._chunked_val2d(np.polynomial.legendre.legval2d, x, y, c, chunk=1000)
    assert chunked.dtype == full.dtype, 'Chunked evaluation changed the output type'
    assert np.array_equal(chunked, full), 'Chunked evaluation changed the result'
