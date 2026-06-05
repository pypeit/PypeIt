"""
Controlled synthetic-data tests for :class:`pypeit.bspline.bspline.bspline`.

These tests verify the mathematical correctness of the fitting algorithms and
provide a numerical baseline for verifying the refactored code.  They are
intentionally independent of the existing golden-file integration tests in
``test_bspline.py``.
"""

import numpy as np
import pytest

from pypeit import bspline


def test_1d_exact_polynomial():
    """
    A cubic B-spline (nord=4) spans the full space of degree-3 polynomials,
    so a cubic polynomial must be represented exactly.  With noise-free,
    unit-weight data the RMS residual must be within a few times machine
    epsilon.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = x**3 - 2.*x**2 + 3.*x - 1.
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)
    err, yfit = sset.fit(x, ytrue, invvar)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'Cubic polynomial not fit to machine precision: '
        f'max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


def test_1d_smooth_function():
    """
    Fit a smooth sinusoidal function with breakpoint spacing much smaller
    than the oscillation period.  With no noise and dense knots the
    B-spline approximation error must be below 1e-3.

    The dominant fourth derivative is from the 0.5*cos(2πx/3) term
    (~2.4 in amplitude); the O(h^4) cubic B-spline error with h=0.5 is
    ~1e-3, so we use atol=2e-3 as a conservative bound.
    """
    n = 1000
    x = np.linspace(0., 10., n)
    ytrue = np.sin(2.*np.pi*x / 5.) + 0.5*np.cos(2.*np.pi*x / 3.)
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=0.5)
    err, yfit = sset.fit(x, ytrue, invvar)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=2e-3), (
        f'Smooth-function residuals too large: '
        f'max = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


def test_1d_fit_workit_equivalence():
    """
    ``fit`` and ``workit`` solve the same weighted normal equations.
    With unit weights the intermediate arrays are bitwise-identical, so the
    coefficients and model evaluations must be exactly equal.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = np.sin(2.*np.pi*x / 5.)
    invvar = np.ones(n, dtype=float)

    # Share an identical knot vector between the two objects.
    fullbkpt = bspline.bspline(x, nord=4, bkspace=1.0).breakpoints.copy()

    sset_fit = bspline.bspline(fullbkpt=fullbkpt, nord=4)
    err_fit, yfit_fit = sset_fit.fit(x, ytrue, invvar)

    sset_workit = bspline.bspline(fullbkpt=fullbkpt, nord=4)
    action, lower, upper = sset_workit.action(x)
    err_workit, yfit_workit = sset_workit.workit(x, ytrue, invvar, action, lower, upper)

    assert err_fit == 0, f'fit returned error code {err_fit}'
    assert err_workit == 0, f'workit returned error code {err_workit}'
    assert np.array_equal(sset_fit.coeff, sset_workit.coeff), (
        'fit and workit produced different coefficients'
    )
    assert np.array_equal(yfit_fit, yfit_workit), (
        'fit and workit produced different model evaluations'
    )


def test_1d_value_matches_fit():
    """
    After a successful fit, calling ``value()`` at the training x must
    return the same model as the ``yfit`` array returned by ``fit()``.
    All points within the fit domain must have ``mask=True``.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = np.sin(2.*np.pi*x / 5.)
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)
    err, yfit = sset.fit(x, ytrue, invvar)
    assert err == 0

    yval, mask = sset.value(x)

    assert mask.all(), 'mask should be True everywhere within the fit domain'
    assert np.array_equal(yval, yfit), (
        'value() at training x disagrees with fit() output'
    )


def test_1d_value_out_of_range_mask():
    """
    ``value()`` must return ``mask=False`` for any x strictly outside the
    valid breakpoint range and ``mask=True`` for points inside it.
    """
    n = 500
    x = np.linspace(2., 8., n)
    ytrue = np.sin(2.*np.pi*x / 5.)
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)
    sset.fit(x, ytrue, invvar)

    x_outside = np.array([0., 1., 9., 10.])     # definitely outside [2, 8]
    x_inside = np.linspace(3., 7., 100)          # safely inside [2, 8]
    x_eval = np.concatenate([x_outside, x_inside])

    _, mask = sset.value(x_eval)

    n_out = len(x_outside)
    assert not mask[:n_out].any(), (
        'Out-of-range points should have mask=False'
    )
    assert mask[n_out:].all(), (
        'In-range points should have mask=True'
    )


# ---------------------------------------------------------------------------
# 2D synthetic-data tests (x2 functionality)
# ---------------------------------------------------------------------------

def test_2d_exact_polynomial():
    """
    With npoly=2 (Legendre) the 2D model spans functions of the form
    h0(x) + h1(x)*x2norm, where x2norm = 2*x2 - 1 (with the default
    xmin=0, xmax=1).  When h0 and h1 are both representable by cubic
    B-splines the fit must be exact.

    x2 must be independent of x so that the P_0 and P_1 columns of the
    design matrix are linearly independent within each B-spline span.
    """
    n = 500
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)                 # independent of x
    x2norm = 2.*x2 - 1.                         # default xmin=0, xmax=1

    h0 = x**3 - 2.*x**2 + 3.*x - 1.
    h1 = x**2 - x + 1.
    ytrue = h0 + h1 * x2norm
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0, npoly=2)
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'2D polynomial not fit to machine precision: '
        f'max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


def test_2d_smooth_function():
    """
    Fit a smooth function with both x and x2 dependence.  The model has
    the form h0(x) + h1(x)*x2norm; with bkspace=0.5 the approximation
    error in each B-spline coefficient function is O(h^4), giving
    residuals below 2e-3.

    x2 must be independent of x so that the normal equations are
    non-degenerate.
    """
    n = 1000
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)                 # independent of x
    x2norm = 2.*x2 - 1.

    ytrue = np.sin(2.*np.pi*x / 5.) + 0.5*np.cos(2.*np.pi*x / 3.) * x2norm
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=0.5, npoly=2)
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=2e-3), (
        f'2D smooth-function residuals too large: '
        f'max = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


def test_2d_value_matches_fit():
    """
    After a 2D fit, calling ``value(x, x2=x2)`` at the training points
    must return the same model as the ``yfit`` returned by ``fit()``.
    All in-range points must have ``mask=True``.
    """
    n = 500
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)                 # independent of x
    x2norm = 2.*x2 - 1.

    ytrue = np.sin(2.*np.pi*x / 5.) + 0.3*np.cos(x) * x2norm
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0, npoly=2)
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)
    assert err == 0

    yval, mask = sset.value(x, x2=x2)

    assert mask.all(), 'mask should be True everywhere within the fit domain'
    assert np.array_equal(yval, yfit), (
        'value() with x2 at training points disagrees with fit() output'
    )
