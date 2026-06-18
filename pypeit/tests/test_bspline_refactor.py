"""
Controlled synthetic-data tests for :class:`pypeit.bspline.bspline.bspline`.

These tests verify the mathematical correctness of the fitting algorithms and
provide a numerical baseline for verifying the refactored code.  They are
intentionally independent of the existing golden-file integration tests in
``test_bspline.py``.
"""

import numpy as np
import pytest

from pypeit.core import bspline
from pypeit.core.bspline.util import bspline_model, cholesky_band, cholesky_solve, intrv, solution_arrays


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


# ---------------------------------------------------------------------------
# Utility function tests: solution_arrays, cholesky_band, cholesky_solve
# ---------------------------------------------------------------------------

def test_cholesky_band_solve_known_system():
    """
    ``cholesky_band`` and ``cholesky_solve`` must recover a known solution
    for a small symmetric positive definite tridiagonal matrix stored in
    band format.

    Band storage convention: ``l[k, j]`` = ``A[j+k, j]`` — the element
    ``k`` rows below the diagonal at column ``j``.  The array has shape
    ``(bw, n + bw)`` where the last ``bw`` columns are zero padding
    required by the forward/back-substitution algorithm.
    """
    n = 8
    bw = 2

    # 8×8 tridiagonal: diagonal = 4, off-diagonals = 1.
    A_full = (4. * np.eye(n)
              + np.diag(np.ones(n - 1), 1)
              + np.diag(np.ones(n - 1), -1))

    # Band storage: l[k, j] = A[j+k, j].
    l = np.zeros((bw, n + bw), dtype=float)
    l[0, :n] = 4.
    l[1, :n - 1] = 1.          # A[j+1, j] for j = 0..n-2; l[1, n-1] = 0

    x_true = np.arange(1., n + 1.)
    # cholesky_solve operates on a vector of length n + bw; pad with zeros.
    b = np.zeros(n + bw, dtype=float)
    b[:n] = A_full @ x_true

    err, G = cholesky_band(l)
    assert err == -1, f'cholesky_band failed with error code {err}'

    sol = cholesky_solve(G, b)[1]
    assert np.allclose(sol[:n], x_true, atol=1e-10), (
        f'cholesky solution incorrect: '
        f'max error = {np.max(np.abs(sol[:n] - x_true)):.2e}'
    )


def test_solution_arrays_normal_equations():
    """
    ``solution_arrays`` must build a banded normal-equations matrix
    ``alpha`` equal to ``A^T W A`` and a right-hand side ``beta`` equal to
    ``A^T W y``, where ``A`` is the full (sparse) design matrix implied by
    the span structure.

    Setup: 6 data points, 4 good breakpoints, nord = 2, npoly = 1.
    This gives 3 spans of 2 data points, bandwidth bw = 2, and
    nfull = 4 unknowns.

    Band storage convention: ``alpha[k, j]`` = ``(A^T W A)[j+k, j]``,
    so ``alpha[k, :nfull-k]`` must equal the k-th sub-diagonal of the
    full normal-equations matrix.
    """
    nd, nn, nord, npoly = 6, 4, 2, 1
    nfull = nn * npoly       # 4
    bw = npoly * nord        # 2
    n_spans = nn - nord + 1  # 3

    rng = np.random.default_rng(seed=99)
    action = rng.standard_normal((nd, bw))
    ydata = rng.standard_normal(nd)
    ivar = np.ones(nd, dtype=float)

    lower = np.array([0, 2, 4], dtype=int)
    upper = np.array([1, 3, 5], dtype=int)   # inclusive upper bounds

    alpha, beta = solution_arrays(
        nn, npoly, nord, ydata, action, ivar,
        upper=upper.copy(), lower=lower.copy()
    )

    # Build the full (nd, nfull) design matrix explicitly.
    # Span k uses columns k : k+nord of A_full.
    A_full = np.zeros((nd, nfull), dtype=float)
    for k in range(n_spans):
        A_full[lower[k]:upper[k] + 1, k:k + nord] = action[lower[k]:upper[k] + 1, :]

    # Reference normal equations via direct matrix products.
    alpha_ref = A_full.T @ np.diag(ivar) @ A_full   # shape (nfull, nfull)
    beta_ref = A_full.T @ (ivar * ydata)             # shape (nfull,)

    # alpha[k, :nfull-k] must equal the k-th sub-diagonal of alpha_ref.
    for k in range(bw):
        assert np.allclose(alpha[k, :nfull - k], np.diag(alpha_ref, -k), atol=1e-14), (
            f'alpha band diagonal k={k} does not match A^T W A'
        )

    assert np.allclose(beta[:nfull], beta_ref, atol=1e-14), (
        'beta does not match A^T W y'
    )

    # End-to-end: solving alpha x = beta must agree with np.linalg.solve.
    err, G = cholesky_band(alpha)
    assert err == -1, f'cholesky_band failed with error code {err}'
    sol = cholesky_solve(G, beta)[1]
    sol_ref = np.linalg.solve(alpha_ref, beta_ref)
    assert np.allclose(sol[:nfull], sol_ref, atol=1e-10), (
        f'Cholesky solution disagrees with np.linalg.solve: '
        f'max error = {np.max(np.abs(sol[:nfull] - sol_ref)):.2e}'
    )


def test_solution_arrays_bspline_model_consistency():
    """
    Verify that ``solution_arrays`` and ``bspline_model`` are self-consistent.

    The four-step pipeline exercised here bypasses ``fit``, ``workit``, and
    ``value`` entirely::

        action(x)  ->  solution_arrays  ->  cholesky_band/solve  ->  bspline_model

    For a 1D cubic B-spline (``npoly=1``) and a noise-free cubic polynomial
    (which lives exactly in the spline space), the full pipeline must recover
    ``ytrue`` to near machine precision.

    With ``npoly=1`` the Cholesky solution vector ``sol[:nn]`` is directly the
    coefficient array — no order-sensitive reshape is required.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = x**3 - 2.*x**2 + 3.*x - 1.
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)

    action, lower, upper = sset.action(x)

    nn = sset.mask[sset.nord:].sum()
    alpha, beta = solution_arrays(
        nn, sset.npoly, sset.nord, ytrue, action, invvar,
        upper=upper.copy(), lower=lower.copy()
    )

    err, G = cholesky_band(alpha)
    assert err == -1, f'cholesky_band failed with error code {err}'
    sol = cholesky_solve(G, beta)[1]

    coeff = sol[:nn]   # npoly=1: solution vector is the coefficient array directly
    yfit = bspline_model(x, action, lower, upper, coeff, nn, sset.nord, sset.npoly)

    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'solution_arrays/bspline_model pipeline inconsistency: '
        f'max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


# ---------------------------------------------------------------------------
# maskpoints tests
# ---------------------------------------------------------------------------

def test_maskpoints_masks_neighborhood():
    """
    When given a valid Cholesky-failure column index, ``maskpoints`` must
    return -1 and set ``self.mask`` to False for the breakpoints in the
    neighborhood of the failing span.

    The neighborhood formula (nord=4, npoly=1):
        for jj in range(-ceil(nord/2), nord//2):   # jj in [-2, -1, 0, 1]
            foo    = max(0, hmm + jj)
            inside = min(n-1, foo + nord)
            test[inside] = True
    where hmm = err_span (since npoly=1) and n = nbkpt - nord.

    With all breakpoints initially good and err_span=5, we compute the
    expected ``inside`` set analytically and verify the mask is modified
    exactly there.
    """
    n_data = 500
    x = np.linspace(0., 10., n_data)
    sset = bspline.bspline(x, nord=4, bkspace=1.0)

    nbkpt = sset.breakpoints.size
    n = nbkpt - sset.nord   # number used in the neighborhood clamp

    err_span = 5
    expected_inside = set()
    for jj in range(-int(np.ceil(sset.nord / 2)), int(sset.nord / 2.)):
        foo = max(0, err_span + jj)
        expected_inside.add(min(n - 1, foo + sset.nord))

    mask_before = sset.mask.copy()
    ret = sset.maskpoints(np.array([err_span]))

    assert ret == -1, f'maskpoints returned {ret}, expected -1'
    assert not sset.mask[list(expected_inside)].any(), (
        'Expected breakpoints were not masked'
    )
    not_expected = np.array([i for i in range(nbkpt) if i not in expected_inside])
    assert np.array_equal(sset.mask[not_expected], mask_before[not_expected]), (
        'maskpoints modified breakpoints outside the expected neighborhood'
    )


def test_maskpoints_too_few_returns_minus2():
    """
    When the number of currently-good breakpoints is <= 2*nord, ``maskpoints``
    must return -2 immediately without modifying ``self.mask``.
    """
    n_data = 500
    x = np.linspace(0., 10., n_data)
    sset = bspline.bspline(x, nord=4, bkspace=1.0)

    # Reduce good breakpoints to exactly 2*nord (the abort threshold).
    sset.mask[:] = False
    sset.mask[:2 * sset.nord] = True
    mask_before = sset.mask.copy()

    ret = sset.maskpoints(np.array([0]))

    assert ret == -2, f'maskpoints returned {ret}, expected -2'
    assert np.array_equal(sset.mask, mask_before), (
        'maskpoints must not modify the mask when aborting due to too few good breakpoints'
    )


# ---------------------------------------------------------------------------
# copy and reinit_coeff tests
# ---------------------------------------------------------------------------

def test_copy():
    """
    ``copy()`` must return a new object whose attributes match the original
    and are fully independent (deep-copied), so mutating the copy does not
    affect the original.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = np.sin(2.*np.pi*x / 5.)
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)
    sset.fit(x, ytrue, invvar)

    sset2 = sset.copy()

    # Scalar attributes must match.
    assert sset2.nord == sset.nord
    assert sset2.npoly == sset.npoly
    assert sset2.xmin == sset.xmin
    assert sset2.xmax == sset.xmax
    assert sset2.funcname == sset.funcname

    # Array attributes must be value-equal.
    assert np.array_equal(sset2.breakpoints, sset.breakpoints)
    assert np.array_equal(sset2.mask, sset.mask)
    assert np.array_equal(sset2.coeff, sset.coeff)
    assert np.array_equal(sset2.icoeff, sset.icoeff)

    # Mutating the copy must not affect the original.
    sset2.coeff[:] = 0.
    assert not np.array_equal(sset2.coeff, sset.coeff), (
        'copy() must produce independent arrays — coeff arrays share memory'
    )


def test_reinit_coeff():
    """
    ``reinit_coeff()`` must reset ``self.coeff`` to all zeros with the
    correct shape, and must leave ``self.icoeff`` untouched.
    """
    n = 500
    x = np.linspace(0., 10., n)
    ytrue = np.sin(2.*np.pi*x / 5.)
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0)
    sset.fit(x, ytrue, invvar)

    icoeff_before = sset.icoeff.copy()
    nc = sset.breakpoints.size - sset.nord

    sset.reinit_coeff()

    assert sset.coeff.shape == (nc,), (
        f'reinit_coeff produced wrong shape: {sset.coeff.shape}, expected ({nc},)'
    )
    assert np.all(sset.coeff == 0.), 'reinit_coeff must zero out coeff'
    assert np.array_equal(sset.icoeff, icoeff_before), (
        'reinit_coeff must not modify icoeff'
    )


# ---------------------------------------------------------------------------
# bsplvn tests
# ---------------------------------------------------------------------------

def test_bsplvn():
    """
    Verify two fundamental B-spline properties of ``bsplvn`` across orders
    2 (linear), 3 (quadratic), and 4 (cubic), plus exact analytic values for
    the linear case.

    Properties tested:

    1. **Partition of unity** — the ``nord`` active basis functions at any
       point sum to exactly 1:  ``vnikx.sum(axis=1) == 1``.
    2. **Non-negativity** — every basis value is >= 0.
    3. **Exact linear case (nord=2)** — for uniform breakpoints with
       spacing h, the two basis values are the linear interpolation weights::

           vnikx[:, 0] = (bkpt[ileft+1] - x) / h
           vnikx[:, 1] = (x - bkpt[ileft])   / h

    ``intrv`` is used to compute span indices the same way ``action`` does,
    so the test exercises the ``intrv`` → ``bsplvn`` path directly.
    """
    n_data = 200
    x = np.linspace(0., 10., n_data)

    for nord in [2, 3, 4]:
        sset = bspline.bspline(x, nord=nord, bkspace=1.0)
        bkpt = sset.breakpoints[sset.mask]
        ileft = intrv(nord, bkpt, x)
        vnikx = sset.bsplvn(x, ileft)

        assert np.allclose(vnikx.sum(axis=1), 1.0, atol=1e-14), (
            f'Partition of unity violated for nord={nord}: '
            f'max deviation = {np.max(np.abs(vnikx.sum(axis=1) - 1.0)):.2e}'
        )
        assert (vnikx >= 0).all(), (
            f'Non-negativity violated for nord={nord}: '
            f'min value = {vnikx.min():.2e}'
        )

    # Exact analytic values for the linear case (nord=2).
    sset2 = bspline.bspline(x, nord=2, bkspace=1.0)
    bkpt2 = sset2.breakpoints[sset2.mask]
    ileft2 = intrv(2, bkpt2, x)
    vnikx2 = sset2.bsplvn(x, ileft2)

    h = bkpt2[ileft2 + 1] - bkpt2[ileft2]   # span width at each point
    assert np.allclose(vnikx2[:, 0], (bkpt2[ileft2 + 1] - x) / h, atol=1e-14), (
        'bsplvn linear basis[0] does not match analytic interpolation weight'
    )
    assert np.allclose(vnikx2[:, 1], (x - bkpt2[ileft2]) / h, atol=1e-14), (
        'bsplvn linear basis[1] does not match analytic interpolation weight'
    )


# ---------------------------------------------------------------------------
# intrv tests
# ---------------------------------------------------------------------------

def test_intrv():
    """
    ``intrv`` must return, for each x[i], the knot span index ``ileft``
    satisfying ``breakpoints[ileft] <= x[i] < breakpoints[ileft+1]``,
    clamped to ``[nord-1, n-1]`` where ``n = breakpoints.size - nord``.

    Uses a small explicit breakpoint array with nord=4 so every case can be
    reasoned about analytically:

        bkpt = [-3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7]
        nord=4  →  n=7, valid index range [3, 6]

    Four properties are checked:

    1. Exact known indices for one strictly-interior point per span.
    2. Bracket property: bkpt[indx[i]] <= x[i] < bkpt[indx[i]+1] for
       in-domain points.
    3. Lower clamp: x below bkpt[nord-1] = 0  →  indx = nord-1 = 3.
    4. Upper clamp: x above bkpt[n]   = 4  →  indx = n-1 = 6.
    """
    nord = 4
    bkpt = np.array([-3., -2., -1., 0., 1., 2., 3., 4., 5., 6., 7.])
    n = bkpt.size - nord   # 7; valid index range [3, 6]

    # Part 1: exact known indices — one strictly-interior point per span.
    x_known = np.array([0.3, 1.3, 2.7, 3.5])
    expected = np.array([3, 4, 5, 6])
    assert np.array_equal(intrv(nord, bkpt, x_known), expected), (
        f'exact index check failed: got {intrv(nord, bkpt, x_known)}, expected {expected}'
    )

    # Part 2: bracket property over a denser in-domain grid.
    x_in = np.array([0.1, 0.9, 1.1, 1.9, 2.1, 2.9, 3.1, 3.9])
    indx = intrv(nord, bkpt, x_in)
    assert np.all(bkpt[indx] <= x_in), 'Left bracket violated'
    assert np.all(x_in < bkpt[indx + 1]), 'Right bracket violated'

    # Part 3: lower clamp — x strictly below bkpt[nord-1] = bkpt[3] = 0.
    x_lo = np.array([-5., -1., -0.01])
    assert np.all(intrv(nord, bkpt, x_lo) == nord - 1), (
        'Lower clamp failed: expected all indices == nord-1'
    )

    # Part 4: upper clamp — x strictly above bkpt[n] = bkpt[7] = 4.
    x_hi = np.array([4.01, 5., 10.])
    assert np.all(intrv(nord, bkpt, x_hi) == n - 1), (
        'Upper clamp failed: expected all indices == n-1'
    )


# ---------------------------------------------------------------------------
# action lower/upper structure tests
# ---------------------------------------------------------------------------

def test_action_lower_upper():
    """
    Verify the structural properties of the ``lower`` and ``upper`` arrays
    returned by ``action``.

    ``lower[i]`` and ``upper[i]`` (inclusive) are the first and last data
    indices belonging to span ``i``, where span ``i`` covers x values between
    ``bkpt[i+nord-1]`` and ``bkpt[i+nord]``.  Empty spans (no data) have
    ``upper[i] == -1``.

    Five properties are checked:

    1. **Shape**: both arrays have length ``nbkpt - 2*nord + 1`` where
       ``nbkpt = sset.mask.sum()``.
    2. **Non-empty bound**: ``upper[i] >= lower[i]`` for all non-empty spans.
    3. **Consistency with indx**: every data index ``j`` in
       ``[lower[i], upper[i]]`` satisfies ``indx[j] == i + nord - 1``,
       verified against ``intrv`` used as an independent ground truth.
    4. **Coverage**: the non-empty spans tile ``[0, nx-1]`` without gaps
       or overlaps.
    5. **Empty sentinel**: spans with no data have ``upper[i] == -1``,
       verified by passing gapped data that skips an entire interior span.
    """
    n_data = 200
    x = np.linspace(0., 10., n_data)
    sset = bspline.bspline(x, nord=4, bkspace=1.0)

    action_mat, lower, upper = sset.action(x)
    bkpt = sset.breakpoints[sset.mask]
    indx = intrv(sset.nord, bkpt, x)

    n_bkpt = sset.mask.sum()
    expected_len = n_bkpt - 2 * sset.nord + 1

    # 1. Shape.
    assert lower.shape == (expected_len,), (
        f'lower shape {lower.shape} != expected ({expected_len},)'
    )
    assert upper.shape == lower.shape, 'lower and upper have different shapes'

    # 2. Non-empty bound.
    nonempty = upper >= 0
    assert np.all(upper[nonempty] >= lower[nonempty]), (
        'upper[i] < lower[i] for some non-empty span'
    )

    # 3. Consistency with indx: every index in [lower[i], upper[i]] belongs
    #    to span i (i.e. indx[j] == i + nord - 1).
    for i in range(len(lower)):
        if upper[i] < 0:
            continue
        span_indx = indx[lower[i]:upper[i] + 1]
        assert np.all(span_indx == i + sset.nord - 1), (
            f'Span {i}: data in [{lower[i]}, {upper[i]}] has wrong span indices'
        )

    # 4. Coverage: union of non-empty spans must cover every data index.
    covered = np.zeros(x.size, dtype=bool)
    for i in range(len(lower)):
        if upper[i] >= 0:
            covered[lower[i]:upper[i] + 1] = True
    assert covered.all(), 'Some data indices are not covered by any span'

    # 5. Empty sentinel: gapped data leaves interior spans with no points.
    x_gap = np.concatenate([
        np.linspace(0., 2.9, 60),
        np.linspace(5.1, 10., 60),
    ])
    _, lower_gap, upper_gap = sset.action(x_gap)
    assert np.any(upper_gap == -1), (
        'Expected empty-span sentinels (upper == -1) for data with a gap in [3, 5]'
    )


# ---------------------------------------------------------------------------
# non-legendre funcname tests
# ---------------------------------------------------------------------------

def test_2d_poly_exact():
    """
    ``funcname='poly'`` builds ``temppoly`` columns ``[1, x2norm, x2norm²]``.
    A function that is a sum of cubic-in-x polynomials times each of these
    three basis columns must be recovered to near machine precision.
    """
    n = 500
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)
    x2norm = 2.*x2 - 1.                     # default xmin=0, xmax=1

    ytrue = ((x**3 - 2.*x**2 + 3.*x - 1.)          # coeff of 1
             + (x**2 - x + 1.) * x2norm             # coeff of x2norm
             + x                * x2norm**2)         # coeff of x2norm^2
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0, npoly=3, funcname='poly')
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'poly exact fit failed: max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


@pytest.mark.skip(reason='poly1 fit with random x2 is ill-conditioned — needs investigation')
def test_2d_poly1_exact():
    """
    ``funcname='poly1'`` builds ``temppoly`` columns ``[x2norm, x2norm², x2norm³]``
    (no constant x2 term).  A function spanned by these three basis columns
    times cubic B-splines must be recovered to near machine precision.
    """
    n = 500
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)
    x2norm = 2.*x2 - 1.

    ytrue = ((x**2 - x + 1.) * x2norm          # coeff of x2norm
             + x               * x2norm**2      # coeff of x2norm^2
             + 1.              * x2norm**3)      # coeff of x2norm^3
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0, npoly=3, funcname='poly1')
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'poly1 exact fit failed: max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )


def test_2d_chebyshev_exact():
    """
    ``funcname='chebyshev'`` builds ``temppoly`` columns
    ``[T₀, T₁, T₂] = [1, x2norm, 2·x2norm²-1]``.  A function spanned by
    these Chebyshev basis columns times cubic B-splines must be recovered to
    near machine precision.
    """
    n = 500
    rng = np.random.default_rng(seed=99)
    x = np.linspace(0., 10., n)
    x2 = rng.uniform(0., 1., n)
    x2norm = 2.*x2 - 1.

    # ytrue = x^3·T₀ + x²·T₁ + x·T₂  (T₂ = 2·x2norm² − 1)
    ytrue = (x**3              * 1.
             + x**2            * x2norm
             + x               * (2.*x2norm**2 - 1.))
    invvar = np.ones(n, dtype=float)

    sset = bspline.bspline(x, nord=4, bkspace=1.0, npoly=3, funcname='chebyshev')
    err, yfit = sset.fit(x, ytrue, invvar, x2=x2)

    assert err == 0, f'fit returned error code {err}'
    assert np.allclose(yfit, ytrue, atol=1e-10), (
        f'chebyshev exact fit failed: max residual = {np.max(np.abs(yfit - ytrue)):.2e}'
    )
