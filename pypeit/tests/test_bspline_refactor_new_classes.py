"""
Per-method unit tests for :class:`~pypeit.bspline.refactor.BSpline` and
:class:`~pypeit.bspline.refactor.BSpline2D`.

Each test function targets a single method or a tightly coupled pair.  The final
section contains cross-check tests that compare numerical output against the original
:class:`~pypeit.bspline.bspline.bspline` implementation.

All random-number generators use a fixed seed for reproducibility.
"""

import warnings

import numpy as np
import pytest

from pypeit.bspline.refactor import BSpline, BSpline2D, Knots
from pypeit.bspline.bspline import bspline


# ============================================================================
# Knots.build
# ============================================================================

def test_build_breakpoints_fullbkpt_returned_sorted():
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 10, 20)
    assert np.all(np.diff(Knots(full=pts, nord=4).breakpoints) >= 0)


def test_build_breakpoints_fullbkpt_padded_when_short():
    pts = np.array([0.0, 1.0, 2.0, 3.0])  # length 4 < 2*nord=8
    assert Knots(full=pts, nord=4).breakpoints.size >= 2 * 4


def test_build_breakpoints_bkspace_strategy():
    x = np.linspace(0, 10, 200)
    kv = Knots(spacing=1.0, x=x, nord=4).breakpoints
    assert kv.min() <= 0.0
    assert kv.max() >= 10.0


def test_build_breakpoints_nbkpts_strategy():
    x = np.linspace(0, 5, 100)
    assert Knots(count=10, x=x, nord=4).breakpoints.size >= 2 * 4


def test_build_breakpoints_everyn_strategy():
    x = np.linspace(0, 10, 300)
    assert Knots(stride=20, x=x, nord=4).breakpoints.size > 2 * 4


def test_build_breakpoints_bkpt_strategy():
    x = np.linspace(0, 10, 200)
    interior = np.array([2.0, 4.0, 6.0, 8.0])
    kv = Knots(interior=interior, x=x, nord=4).breakpoints
    assert kv[0] <= 0.0
    assert kv[-1] >= 10.0


def test_build_breakpoints_phantom_knots_at_each_end():
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=6, x=x, nord=nord).breakpoints
    interior_min = kv[nord - 1]
    interior_max = kv[-(nord)]
    assert np.all(kv[:nord - 1] < interior_min)
    assert np.all(kv[-(nord - 1):] > interior_max)


def test_build_breakpoints_raises_without_x_or_fullbkpt():
    with pytest.raises(ValueError):
        Knots(spacing=1.0).build(None, 4)



# ============================================================================
# BSpline._find_spans
# ============================================================================

def test_find_spans_all_indices_in_range():
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    n = kv.size - nord
    indx = BSpline._find_spans(x, kv, nord)
    assert np.all(indx >= nord - 1)
    assert np.all(indx <= n - 1)


def test_find_spans_bracketing():
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    indx = BSpline._find_spans(x, kv, nord)
    for i, (xi, il) in enumerate(zip(x, indx)):
        assert kv[il] <= xi
        if il < kv.size - 1:
            assert xi <= kv[il + 1] + 1e-12


def test_find_spans_clamped_at_lower_edge():
    kv = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    nord = 2
    x = np.array([0.0])
    indx = BSpline._find_spans(x, kv, nord)
    assert indx[0] >= nord - 1


def test_find_spans_monotone_input_gives_monotone_output():
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    indx = BSpline._find_spans(x, kv, nord)
    assert np.all(np.diff(indx) >= 0)


# ============================================================================
# BSpline._bspline_basis
# ============================================================================

def test_bspline_basis_output_shape():
    x = np.sort(np.random.default_rng(7).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.shape == (x.size, spl.nord)


def test_bspline_basis_c_order():
    x = np.sort(np.random.default_rng(8).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.flags['C_CONTIGUOUS']


def test_bspline_basis_partition_of_unity():
    x = np.sort(np.random.default_rng(9).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    basis = spl._bspline_basis(x, indx)
    np.testing.assert_allclose(basis.sum(axis=1), 1.0, atol=1e-12)


def test_bspline_basis_non_negative():
    x = np.sort(np.random.default_rng(10).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert np.all(basis >= -1e-14)


def test_bspline_basis_linear_case():
    kv = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)
    spl = BSpline(knots=kv, nord=2)
    x = np.array([1.5])
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], 2)
    basis = spl._bspline_basis(x, indx)
    np.testing.assert_allclose(basis.sum(axis=1), 1.0, atol=1e-14)


# ============================================================================
# BSpline._build_design_matrix
# ============================================================================

def test_build_design_matrix_shape():
    x = np.sort(np.random.default_rng(42).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.shape == (x.size, spl.nord)


def test_build_design_matrix_lower_upper_lengths():
    x = np.sort(np.random.default_rng(43).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    n = spl.mask.sum() - spl.nord
    assert lower.size == n - spl.nord + 1
    assert upper.size == n - spl.nord + 1


def test_build_design_matrix_full_data_coverage():
    x = np.sort(np.random.default_rng(44).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    covered = np.zeros(x.size, dtype=bool)
    for k in range(lower.size):
        if lower[k] <= upper[k]:
            covered[lower[k]:upper[k] + 1] = True
    assert covered.all()


def test_build_design_matrix_consistent_with_find_spans():
    x = np.sort(np.random.default_rng(45).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    for k in range(lower.size):
        mask_k = (indx == k + spl.nord - 1)
        if mask_k.any():
            assert lower[k] == np.where(mask_k)[0][0]
            assert upper[k] == np.where(mask_k)[0][-1]


# ============================================================================
# BSpline._assemble_normal_equations
# ============================================================================

def test_assemble_normal_equations_alpha_shape():
    rng = np.random.default_rng(13)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x) + 0.1 * rng.standard_normal(80)
    w = np.ones(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, w, lower, upper)
    bw = A.shape[1]
    nn = spl.mask[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert alpha.shape == (bw, nfull + bw)
    assert beta.shape == (nfull + bw,)


def test_assemble_normal_equations_diagonal_positive():
    rng = np.random.default_rng(14)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x) + 0.1 * rng.standard_normal(80)
    w = np.ones(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, w, lower, upper)
    nn = spl.mask[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert np.all(alpha[0, :nfull] > 0)


def test_assemble_normal_equations_solution_matches_lstsq():
    rng = np.random.default_rng(15)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x)
    invvar = np.ones(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)

    # Banded solve via fit()
    err, yfit_banded = spl.fit(x, y, invvar)
    assert err == 0

    # Build the full (N, nn) design matrix: A_full[i, k:k+nord] = A[i, :]
    # for rows i in span k, then solve with lstsq for reference.
    A, lower, upper = spl._build_design_matrix(x)
    nn = spl.mask[spl.nord:].sum()
    A_full = np.zeros((x.size, nn))
    for k in range(lower.size):
        if lower[k] <= upper[k]:
            sl = slice(lower[k], upper[k] + 1)
            A_full[sl, k:k + spl.nord] = A[sl, :]

    sol_lstsq, _, _, _ = np.linalg.lstsq(A_full, y, rcond=None)
    yfit_lstsq = A_full @ sol_lstsq

    np.testing.assert_allclose(yfit_banded, yfit_lstsq, atol=1e-10)


# ============================================================================
# BSpline._solve_banded
# ============================================================================

def _make_tridiagonal_alpha_beta(n=10):
    """Build a small SPD banded (tridiagonal) system in the format expected by
    `_solve_banded`: alpha shape (bw, nfull + bw) with bw=2.

    scipy lower banded format: alpha[0, j] = A[j,j] (diagonal),
    alpha[1, j] = A[j+1, j] (subdiagonal), so off-diagonal goes at alpha[1, :n-1].
    """
    d = 4.0 * np.ones(n)
    u = -1.0 * np.ones(n - 1)
    bw = 2
    alpha = np.zeros((bw, n + bw))
    alpha[0, :n] = d
    alpha[1, :n - 1] = u   # alpha[1, j] = A[j+1, j] for j = 0..n-2
    beta = np.zeros(n + bw)
    beta[:n] = np.arange(1, n + 1, dtype=float)
    return alpha, beta, n


def test_solve_banded_success_path():
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert bad_cols[0] == -1
    assert sol is not None
    assert chol is not None
    assert sol.shape == (n,)


def test_solve_banded_solution_correct():
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    # scipy lower: alpha[1, j] = A[j+1, j]; off-diagonals sit at alpha[1, :n-1]
    A_dense = (np.diag(alpha[0, :n])
               + np.diag(alpha[1, :n - 1], k=1)
               + np.diag(alpha[1, :n - 1], k=-1))
    ref = np.linalg.solve(A_dense, beta[:n])
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    np.testing.assert_allclose(sol, ref, atol=1e-10)


def test_solve_banded_failure_returns_bad_cols():
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    alpha[0, 3] = -1.0  # corrupt diagonal
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert sol is None
    assert 3 in bad_cols


def test_solve_banded_chol_diagonal_positive():
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert np.all(chol[0, :n] > 0)


# ============================================================================
# BSpline._evaluate_model
# ============================================================================

def test_evaluate_model_matches_explicit_matmul():
    rng = np.random.default_rng(99)
    x = np.sort(rng.uniform(0, 5, 60))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    spl.coeff[:] = rng.standard_normal(spl.coeff.size)

    yfit = spl._evaluate_model(A, lower, upper)

    # Verify one non-empty span manually
    coeffbk = spl.mask[spl.nord:].nonzero()[0]
    goodcoeff = spl.coeff[coeffbk]
    k = next(k for k in range(lower.size) if lower[k] <= upper[k])
    sl = slice(lower[k], upper[k] + 1)
    np.testing.assert_allclose(yfit[sl], A[sl, :] @ goodcoeff[k:k + spl.nord], atol=1e-14)


# ============================================================================
# BSpline._mask_breakpoints
# ============================================================================

def test_mask_breakpoints_masks_neighbourhood():
    x = np.sort(np.random.default_rng(5).uniform(0, 10, 200))
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    n_before = spl.mask.sum()
    result = spl._mask_breakpoints(np.array([8]))
    assert result == -1
    assert spl.mask.sum() < n_before


def test_mask_breakpoints_invalidates_cache():
    x = np.sort(np.random.default_rng(6).uniform(0, 10, 200))
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl._cached_design = spl._build_design_matrix(x)
    spl._mask_breakpoints(np.array([8]))
    assert spl._cached_design is None


def test_mask_breakpoints_too_few_returns_minus2():
    x = np.linspace(0, 5, 30)
    spl = BSpline(x=x, knots=Knots(count=3), nord=4)
    spl.mask[:] = False
    spl.mask[:2 * 4] = True  # leave only 2*nord active
    result = spl._mask_breakpoints(np.array([0]))
    assert result == -2


# ============================================================================
# BSpline.fit
# ============================================================================

def test_fit_cubic_polynomial_recovery():
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, 300))
    y = 1.0 + 2.0*x - 0.5*x**2 + 0.1*x**3
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, invvar)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-8)


def test_fit_smooth_function_residuals():
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 2 * np.pi, 500))
    y = np.sin(x)
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y, invvar)
    assert err == 0
    assert np.std(yfit - y) < 2e-3


def test_fit_returns_correct_length():
    rng = np.random.default_rng(2)
    x = np.sort(rng.uniform(0, 5, 100))
    y = rng.standard_normal(100)
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    err, yfit = spl.fit(x, y, invvar)
    assert yfit.shape == x.shape


def test_fit_design_matrix_cached():
    rng = np.random.default_rng(3)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, y, invvar)
    cache_id = id(spl._cached_design)
    spl.fit(x, y * 0.9, invvar)
    assert id(spl._cached_design) == cache_id


def test_fit_zero_invvar_points_ignored():
    rng = np.random.default_rng(4)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    invvar = np.ones_like(x)
    invvar[::5] = 0.0
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    err, yfit = spl.fit(x, y, invvar)
    assert err == 0
    assert np.isfinite(yfit).all()


def test_fit_reset_knots_rebuilds_breakpoints():
    rng = np.random.default_rng(50)
    x1 = np.sort(rng.uniform(0, 5, 100))
    invvar = np.ones(100)
    spl = BSpline(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1), invvar)
    bkpt_before = spl.breakpoints.copy()
    x2 = np.sort(rng.uniform(0, 10, 100))
    spl.fit(x2, np.sin(x2), invvar, reset_knots=True)
    assert spl.breakpoints.min() <= x2.min()
    assert spl.breakpoints.max() >= x2.max()
    assert not np.array_equal(spl.breakpoints, bkpt_before)


def test_fit_reset_knots_resets_mask():
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    invvar = np.ones(100)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y, invvar)
    spl.mask[5] = False
    assert not spl.mask.all()
    spl.fit(x, y, invvar, reset_knots=True)
    assert spl.mask.all()


# ============================================================================
# BSpline.value
# ============================================================================

def test_value_matches_fit_at_training_points():
    rng = np.random.default_rng(10)
    x = np.sort(rng.uniform(0, 10, 200))
    y = 2.0 + x - 0.3 * x**2
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=12), nord=4)
    err, yfit_fit = spl.fit(x, y, invvar)
    yfit_val, _ = spl.value(x)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_value_out_of_range_masked():
    x = np.linspace(1, 9, 100)
    y = np.sin(x)
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y, invvar)
    x_eval = np.array([0.0, 5.0, 10.0])
    yfit, mask = spl.value(x_eval)
    assert not mask[0]
    assert mask[1]
    assert not mask[2]


def test_value_unsorted_input():
    rng = np.random.default_rng(11)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.cos(x)
    invvar = np.ones_like(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, y, invvar)
    x_unsorted = rng.permutation(x)
    yfit, _ = spl.value(x_unsorted)
    yfit_sorted, _ = spl.value(np.sort(x_unsorted))
    np.testing.assert_allclose(yfit[np.argsort(x_unsorted)], yfit_sorted, atol=1e-12)


# ============================================================================
# BSpline2D._normalize_x2
# ============================================================================

def test_normalize_x2_maps_xmin_to_minus_one():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_x2(np.array([2.0])), [-1.0])


def test_normalize_x2_maps_xmax_to_plus_one():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_x2(np.array([8.0])), [1.0])


def test_normalize_x2_midpoint_is_zero():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 4.0
    np.testing.assert_allclose(spl._normalize_x2(np.array([2.0])), [0.0])


def test_normalize_x2_linear_mapping():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 10.0
    x2 = np.linspace(0, 10, 11)
    np.testing.assert_allclose(spl._normalize_x2(x2), np.linspace(-1, 1, 11), atol=1e-14)


# ============================================================================
# BSpline2D._poly_basis
# ============================================================================

def _make_poly_spl(funcname, npoly):
    spl = BSpline2D.__new__(BSpline2D)
    spl.npoly = npoly
    spl.funcname = funcname
    return spl


def test_poly_basis_legendre_shape():
    spl = _make_poly_spl('legendre', 4)
    P = spl._poly_basis(np.linspace(-1, 1, 20))
    assert P.shape == (20, 4)


def test_poly_basis_legendre_first_column_is_one():
    spl = _make_poly_spl('legendre', 3)
    P = spl._poly_basis(np.linspace(-1, 1, 10))
    np.testing.assert_allclose(P[:, 0], 1.0, atol=1e-14)


def test_poly_basis_chebyshev_shape():
    spl = _make_poly_spl('chebyshev', 3)
    P = spl._poly_basis(np.linspace(-1, 1, 15))
    assert P.shape == (15, 3)


def test_poly_basis_chebyshev_first_is_one_second_is_x():
    spl = _make_poly_spl('chebyshev', 2)
    x2norm = np.linspace(-1, 1, 10)
    P = spl._poly_basis(x2norm)
    np.testing.assert_allclose(P[:, 0], 1.0, atol=1e-14)
    np.testing.assert_allclose(P[:, 1], x2norm, atol=1e-14)


def test_poly_basis_poly_monomial_structure():
    spl = _make_poly_spl('poly', 4)
    x2norm = np.array([0.0, 0.5, 1.0])
    P = spl._poly_basis(x2norm)
    for k in range(4):
        np.testing.assert_allclose(P[:, k], x2norm**k, atol=1e-14)


def test_poly_basis_unknown_funcname_raises():
    spl = _make_poly_spl('invalid_basis', 2)
    with pytest.raises(ValueError):
        spl._poly_basis(np.array([0.0]))


# ============================================================================
# BSpline2D._build_design_matrix
# ============================================================================

def test_bspline2d_build_design_matrix_shape():
    rng = np.random.default_rng(20)
    x = np.sort(rng.uniform(0, 10, 100))
    x2 = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, npoly=3, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x, x2)
    assert A.shape == (x.size, spl.nord * spl.npoly)


def test_bspline2d_build_design_matrix_outer_product_structure():
    rng = np.random.default_rng(21)
    N = 50
    x = np.sort(rng.uniform(0, 5, N))
    x2 = rng.uniform(0, 1, N)
    npoly = 2
    spl = BSpline2D(x=x, npoly=npoly, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x, x2)

    indx = BSpline._find_spans(x, spl.breakpoints[spl.mask], spl.nord)
    B = spl._bspline_basis(x, indx)
    P = spl._poly_basis(spl._normalize_x2(x2))

    for ii in range(spl.nord):
        for jj in range(npoly):
            col = ii * npoly + jj
            np.testing.assert_allclose(A[:, col], B[:, ii] * P[:, jj], atol=1e-14,
                                       err_msg=f'Outer product mismatch at col {col}')


def test_bspline2d_build_design_matrix_x2_size_mismatch_raises():
    rng = np.random.default_rng(22)
    x = np.sort(rng.uniform(0, 10, 100))
    x2 = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, npoly=3, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=8), nord=4)
    with pytest.raises(ValueError):
        spl._build_design_matrix(x, x2[:-5])


def test_bspline2d_build_design_matrix_c_order():
    rng = np.random.default_rng(23)
    x = np.sort(rng.uniform(0, 10, 100))
    x2 = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, npoly=3, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x, x2)
    assert A.flags['C_CONTIGUOUS']


# ============================================================================
# BSpline2D.fit
# ============================================================================

def test_bspline2d_fit_exact_polynomial_legendre():
    """h0(x) + h1(x) * P1(x2) with Legendre basis should be recovered exactly."""
    rng = np.random.default_rng(30)
    N = 400
    x = np.sort(rng.uniform(0, 10, N))
    x2 = rng.uniform(0, 1, N)
    x2norm = 2 * x2 - 1
    y = 1.0 + 0.5 * x + 0.3 * x2norm
    invvar = np.ones(N)
    spl = BSpline2D(x=x, npoly=2, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, invvar, x2)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_smooth_function():
    rng = np.random.default_rng(31)
    N = 500
    x = np.sort(rng.uniform(0, 2 * np.pi, N))
    x2 = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.2 * x2)
    invvar = np.ones(N)
    spl = BSpline2D(x=x, npoly=3, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y, invvar, x2)
    assert err == 0
    assert np.std(yfit - y) < 2e-2


def test_bspline2d_fit_exact_polynomial_chebyshev():
    rng = np.random.default_rng(32)
    N = 400
    x = np.sort(rng.uniform(0, 8, N))
    x2 = rng.uniform(-1, 1, N)
    y = 2.0 + x2  # T0 + T1 in Chebyshev basis
    invvar = np.ones(N)
    spl = BSpline2D(x=x, npoly=2, xmin=-1.0, xmax=1.0, funcname='chebyshev',
                    knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, invvar, x2)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_exact_polynomial_poly():
    rng = np.random.default_rng(33)
    N = 400
    x = np.sort(rng.uniform(0, 5, N))
    x2 = rng.uniform(0, 2, N)
    x2norm = 2.0 * (x2 - 1.0) / 2.0 - 0.0   # xmin=0, xmax=2 → norm range [-1, 1]
    y = 1.0 + x2norm
    invvar = np.ones(N)
    spl = BSpline2D(x=x, npoly=2, xmin=0.0, xmax=2.0, funcname='poly',
                    knots=Knots(count=15), nord=4)
    err, yfit = spl.fit(x, y, invvar, x2)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_reset_knots_rebuilds_breakpoints():
    rng = np.random.default_rng(52)
    N = 200
    invvar = np.ones(N)
    x1 = np.sort(rng.uniform(0, 5, N))
    x2_1 = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x1, npoly=2, xmin=0.0, xmax=1.0, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1), invvar, x2_1)
    bkpt_before = spl.breakpoints.copy()
    x3 = np.sort(rng.uniform(0, 10, N))
    x2_2 = rng.uniform(0, 1, N)
    spl.fit(x3, np.sin(x3), invvar, x2_2, reset_knots=True)
    assert spl.breakpoints.min() <= x3.min()
    assert spl.breakpoints.max() >= x3.max()
    assert not np.array_equal(spl.breakpoints, bkpt_before)


def test_bspline2d_fit_reset_knots_resets_mask():
    rng = np.random.default_rng(53)
    N = 200
    invvar = np.ones(N)
    x = np.sort(rng.uniform(0, 5, N))
    x2 = rng.uniform(0, 1, N)
    y = np.sin(x)
    spl = BSpline2D(x=x, npoly=2, xmin=0.0, xmax=1.0, knots=Knots(count=10), nord=4)
    spl.fit(x, y, invvar, x2)
    spl.mask[5] = False
    assert not spl.mask.all()
    spl.fit(x, y, invvar, x2, reset_knots=True)
    assert spl.mask.all()


# ============================================================================
# BSpline2D.value
# ============================================================================

def test_bspline2d_value_matches_fit_at_training_points():
    rng = np.random.default_rng(40)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    x2 = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * x2
    invvar = np.ones(N)
    spl = BSpline2D(x=x, npoly=2, xmin=0.0, xmax=1.0, funcname='legendre',
                    knots=Knots(count=12), nord=4)
    err, yfit_fit = spl.fit(x, y, invvar, x2)
    yfit_val, _ = spl.value(x, x2)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_bspline2d_value_x2_required():
    rng = np.random.default_rng(41)
    x = np.sort(rng.uniform(0, 5, 50))
    x2 = rng.uniform(0, 1, 50)
    spl = BSpline2D(x=x, npoly=2, xmin=0.0, xmax=1.0,
                    knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), np.ones(50), x2)
    with pytest.raises(TypeError):
        spl.value(x)  # missing required x2


# ============================================================================
# Cross-check tests: new vs. original implementation
# ============================================================================

def test_crosscheck_1d_fit_yfit_match():
    rng = np.random.default_rng(50)
    x = np.sort(rng.uniform(0, 10, 200))
    y = np.sin(x) + 0.05 * rng.standard_normal(200)
    invvar = np.ones_like(x)
    nord = 4
    nbkpts = 15

    old = bspline(x=x, nord=nord, npoly=1, nbkpts=nbkpts)
    old_err, old_yfit = old.fit(x, y, invvar)

    new = BSpline(knots=old.breakpoints.copy(), nord=nord)
    new.mask[:] = old.mask
    new._cached_design = None
    new_err, new_yfit = new.fit(x, y, invvar)

    assert old_err == 0 and new_err == 0
    np.testing.assert_allclose(new_yfit, old_yfit, rtol=1e-10)


def test_crosscheck_1d_coeff_match():
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 10, 200))
    y = np.sin(x) + 0.05 * rng.standard_normal(200)
    invvar = np.ones_like(x)

    old = bspline(x=x, nord=4, npoly=1, nbkpts=15)
    old.fit(x, y, invvar)

    new = BSpline(knots=old.breakpoints.copy(), nord=4)
    new.mask[:] = old.mask
    new._cached_design = None
    new.fit(x, y, invvar)

    np.testing.assert_allclose(new.coeff, old.coeff, rtol=1e-10)


def test_crosscheck_2d_fit_yfit_match():
    rng = np.random.default_rng(52)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    x2 = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.3 * x2)
    invvar = np.ones(N)
    nord = 4
    npoly = 2
    nbkpts = 12

    old = bspline(x=x, nord=nord, npoly=npoly, nbkpts=nbkpts, funcname='legendre')
    old.xmin = 0.0
    old.xmax = 1.0
    old_err, old_yfit = old.fit(x, y, invvar, x2=x2)

    new = BSpline2D(knots=old.breakpoints.copy(), nord=nord, npoly=npoly,
                    xmin=0.0, xmax=1.0, funcname='legendre')
    new.mask[:] = old.mask
    new._cached_design = None
    new_err, new_yfit = new.fit(x, y, invvar, x2)

    assert old_err == 0 and new_err == 0
    np.testing.assert_allclose(new_yfit, old_yfit, rtol=1e-8)


def test_crosscheck_2d_coeff_transposed():
    """new.coeff should equal old.coeff.T (shape change from (npoly,nc) to (nc,npoly))."""
    rng = np.random.default_rng(53)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    x2 = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.3 * x2)
    invvar = np.ones(N)

    old = bspline(x=x, nord=4, npoly=2, nbkpts=12, funcname='legendre')
    old.xmin = 0.0
    old.xmax = 1.0
    old.fit(x, y, invvar, x2=x2)

    new = BSpline2D(knots=old.breakpoints.copy(), nord=4, npoly=2,
                    xmin=0.0, xmax=1.0, funcname='legendre')
    new.mask[:] = old.mask
    new._cached_design = None
    new.fit(x, y, invvar, x2)

    np.testing.assert_allclose(new.coeff, old.coeff.T, rtol=1e-8)


def test_crosscheck_1d_value_match():
    rng = np.random.default_rng(54)
    x = np.sort(rng.uniform(0, 10, 150))
    y = 1 + 2*x - 0.1*x**2
    invvar = np.ones_like(x)

    old = bspline(x=x, nord=4, npoly=1, nbkpts=12)
    old.fit(x, y, invvar)

    new = BSpline(knots=old.breakpoints.copy(), nord=4)
    new.mask[:] = old.mask
    new._cached_design = None
    new.fit(x, y, invvar)

    x_eval = np.sort(rng.uniform(0, 10, 50))
    old_val, _ = old.value(x_eval)
    new_val, _ = new.value(x_eval)
    np.testing.assert_allclose(new_val, old_val, rtol=1e-10)
