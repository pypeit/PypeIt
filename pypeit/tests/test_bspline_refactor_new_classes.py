"""
Per-method unit tests for :class:`~pypeit.bspline.refactor.BSpline` and
:class:`~pypeit.bspline.refactor.BSpline2D`.

Each test function targets a single method or a tightly coupled pair.  The final
section contains cross-check tests that compare numerical output against the original
:class:`~pypeit.bspline.bspline.bspline` implementation.

All random-number generators use a fixed seed for reproducibility.
"""

import warnings

from IPython import embed
import numpy as np
import pytest

from pypeit import dataPaths
from pypeit.bspline.refactor import BSpline, BSpline2D, Knots, bspline_profile_refactor
from pypeit.bspline.bspline import bspline
from pypeit.core.fitting import bspline_profile
from pypeit.core.basis import fchebyshev, flegendre


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


def test_build_breakpoints_spacing_too_large_fallback():
    x = np.linspace(0, 5, 100)
    kv = Knots(spacing=10.0, x=x, nord=4).breakpoints
    assert kv is not None
    assert kv.min() <= 0.0
    assert kv.max() >= 5.0


def test_build_breakpoints_stride_too_large_fallback():
    x = np.linspace(0, 5, 50)
    kv = Knots(stride=100, x=x, nord=4).breakpoints
    assert kv is not None
    assert kv.min() <= 0.0
    assert kv.max() >= 5.0


def test_build_breakpoints_interior_outside_range_fallback():
    x = np.linspace(0, 5, 100)
    interior = np.array([7.0, 8.0, 9.0])
    kv = Knots(interior=interior, x=x, nord=4).breakpoints
    assert kv is not None
    assert kv.min() <= 0.0
    assert kv.max() >= 5.0


def test_build_breakpoints_no_strategy_raises_valueerror():
    with pytest.raises(ValueError):
        Knots(stride=None).build(np.linspace(0, 5, 50), 4)


def test_knots_init_nord_without_x_warns():
    with pytest.warns(UserWarning):
        Knots(count=10, nord=4)


# ============================================================================
# Knots.copy
# ============================================================================

def test_knots_copy_strategy_params_preserved():
    orig = Knots(count=10, spread=1.5, stride=None)
    cp = orig.copy()
    assert cp.count == orig.count
    assert cp.spread == orig.spread
    assert cp.stride == orig.stride
    assert cp.spacing == orig.spacing


def test_knots_copy_breakpoints_deep_copied():
    orig = Knots(count=10, x=np.linspace(0, 5, 100), nord=4)
    cp = orig.copy()
    assert np.array_equal(cp.breakpoints, orig.breakpoints)
    orig_val = orig.breakpoints[0]
    cp._breakpoints[0] = 999.0
    assert orig.breakpoints[0] == orig_val


def test_knots_copy_unbuilt_returns_none_breakpoints():
    orig = Knots(count=10)
    cp = orig.copy()
    assert cp.breakpoints is None


# ============================================================================
# BSpline.__init__ — knots type validation
# ============================================================================

def test_bspline_init_invalid_knots_type_raises():
    with pytest.raises(TypeError):
        BSpline(knots='not_a_valid_knots_type')


# ============================================================================
# BSpline._uniq
# ============================================================================

def test_uniq_empty_array_raises():
    with pytest.raises(ValueError):
        BSpline._uniq(np.array([]))


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
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.shape == (x.size, spl.nord)


def test_bspline_basis_c_order():
    x = np.sort(np.random.default_rng(8).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.flags['C_CONTIGUOUS']


def test_bspline_basis_partition_of_unity():
    x = np.sort(np.random.default_rng(9).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    np.testing.assert_allclose(basis.sum(axis=1), 1.0, atol=1e-12)


def test_bspline_basis_non_negative():
    x = np.sort(np.random.default_rng(10).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert np.all(basis >= -1e-14)


def test_bspline_basis_linear_case():
    kv = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0], dtype=float)
    spl = BSpline(knots=kv, nord=2)
    x = np.array([1.5])
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], 2)
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
    n = spl.bkpt_gpm.sum() - spl.nord
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
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
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
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, lower, upper)
    bw = A.shape[1]
    nn = spl.bkpt_gpm[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert alpha.shape == (bw, nfull + bw)
    assert beta.shape == (nfull + bw,)


def test_assemble_normal_equations_diagonal_positive():
    rng = np.random.default_rng(14)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x) + 0.1 * rng.standard_normal(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, lower, upper)
    nn = spl.bkpt_gpm[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert np.all(alpha[0, :nfull] > 0)


def test_assemble_normal_equations_solution_matches_lstsq():
    rng = np.random.default_rng(15)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)

    # Banded solve via fit()
    err, yfit_banded = spl.fit(x, y)
    assert err == 0

    # Build the full (N, nn) design matrix: A_full[i, k:k+nord] = A[i, :]
    # for rows i in span k, then solve with lstsq for reference.
    A, lower, upper = spl._build_design_matrix(x)
    nn = spl.bkpt_gpm[spl.nord:].sum()
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


def test_solve_banded_linalg_error_branch():
    # [[1, 2], [2, 1]] is indefinite (det = -3); positive diagonal passes
    # the pre-check but Cholesky decomposition raises LinAlgError.
    bw, n = 2, 2
    alpha = np.zeros((bw, n + bw))
    alpha[0, :n] = [1.0, 1.0]  # positive diagonal
    alpha[1, 0] = 2.0           # large off-diagonal -> not positive definite
    beta = np.zeros(n + bw)
    beta[:n] = [1.0, 1.0]
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert sol is None
    assert bad_cols[0] != -1


# ============================================================================
# BSpline._evaluate_model
# ============================================================================

def test_evaluate_model_matches_explicit_matmul():
    rng = np.random.default_rng(99)
    x = np.sort(rng.uniform(0, 5, 60))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.reset_coeff()
    A, lower, upper = spl._build_design_matrix(x)
    spl.coeff[:] = rng.standard_normal(spl.coeff.size)

    yfit = spl._evaluate_model(A, lower, upper)

    # Verify one non-empty span manually
    coeffbk = spl.bkpt_gpm[spl.nord:].nonzero()[0]
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
    n_before = spl.bkpt_gpm.sum()
    result = spl._mask_breakpoints(np.array([8]))
    assert result == -1
    assert spl.bkpt_gpm.sum() < n_before


def test_mask_breakpoints_invalidates_cache():
    x = np.sort(np.random.default_rng(6).uniform(0, 10, 200))
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl._cached_design = spl._build_design_matrix(x)
    spl._mask_breakpoints(np.array([8]))
    assert spl._cached_design is None


def test_mask_breakpoints_too_few_returns_minus2():
    x = np.linspace(0, 5, 30)
    spl = BSpline(x=x, knots=Knots(count=3), nord=4)
    spl.bkpt_gpm[:] = False
    spl.bkpt_gpm[:2 * 4] = True  # leave only 2*nord active
    result = spl._mask_breakpoints(np.array([0]))
    assert result == -2


# ============================================================================
# BSpline.fit
# ============================================================================

def test_fit_cubic_polynomial_recovery():
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, 300))
    y = 1.0 + 2.0*x - 0.5*x**2 + 0.1*x**3
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-8)


def test_fit_smooth_function_residuals():
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 2 * np.pi, 500))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == 0
    assert np.std(yfit - y) < 2e-3


def test_fit_returns_correct_length():
    rng = np.random.default_rng(2)
    x = np.sort(rng.uniform(0, 5, 100))
    y = rng.standard_normal(100)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    err, yfit = spl.fit(x, y)
    assert yfit.shape == x.shape


def test_fit_design_matrix_cached():
    rng = np.random.default_rng(3)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, y)
    cache_id = id(spl._cached_design)
    spl.fit(x, y * 0.9)
    assert id(spl._cached_design) == cache_id


def test_fit_zero_invvar_points_ignored():
    rng = np.random.default_rng(4)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    ivar = np.ones_like(x)
    ivar[::5] = 0.0
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    err, yfit = spl.fit(x, y, ivar=ivar)
    assert err == 0
    assert np.isfinite(yfit).all()


def test_fit_reset_knots_rebuilds_breakpoints():
    rng = np.random.default_rng(50)
    x1 = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1))
    bkpt_before = spl.breakpoints.copy()
    x2 = np.sort(rng.uniform(0, 10, 100))
    spl.fit(x2, np.sin(x2), reset_knots=True)
    assert spl.breakpoints.min() <= x2.min()
    assert spl.breakpoints.max() >= x2.max()
    assert not np.array_equal(spl.breakpoints, bkpt_before)


def test_fit_reset_knots_resets_mask():
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y)
    spl.bkpt_gpm[5] = False
    assert not spl.bkpt_gpm.all()
    spl.fit(x, y, reset_knots=True)
    assert spl.bkpt_gpm.all()


def test_fit_returns_minus1_on_degenerate_cholesky():
    rng = np.random.default_rng(60)
    # Data only in [0, 2]; knots span [0, 10] -> spans in [2, 10] are empty,
    # leaving zero diagonal entries that trigger the degenerate pre-check.
    x = np.sort(rng.uniform(0, 2, 100))
    y = np.sin(x)
    full_bkpt = np.linspace(0, 10, 40)
    spl = BSpline(knots=Knots(full=full_bkpt), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == -1
    assert yfit.shape == x.shape


def test_fit_ivar_none_matches_unit_weights():
    rng = np.random.default_rng(70)
    x = np.sort(rng.uniform(0, 5, 150))
    y = np.sin(x)
    spl_w = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl_n = BSpline(x=x, knots=Knots(count=10), nord=4)
    err_w, yfit_w = spl_w.fit(x, y, ivar=np.ones(150))
    err_n, yfit_n = spl_n.fit(x, y)
    assert err_w == 0 and err_n == 0
    np.testing.assert_array_equal(yfit_w, yfit_n)


def test_fit_nonuniform_ivar_affects_result():
    rng = np.random.default_rng(72)
    N = 300
    x = np.sort(rng.uniform(0.5, 5.0, N))
    ivar = 1.0 / x**2  # noise amplitude proportional to x
    y = np.sin(x) + rng.standard_normal(N) / np.sqrt(ivar)
    spl_u = BSpline(x=x, knots=Knots(count=15), nord=4)
    spl_n = BSpline(x=x, knots=Knots(count=15), nord=4)
    err_u, yfit_u = spl_u.fit(x, y)
    err_n, yfit_n = spl_n.fit(x, y, ivar=ivar)
    assert err_u == 0 and err_n == 0
    assert np.isfinite(yfit_n).all()
    assert not np.allclose(yfit_u, yfit_n)


# ============================================================================
# BSpline.value
# ============================================================================

def test_value_matches_fit_at_training_points():
    rng = np.random.default_rng(10)
    x = np.sort(rng.uniform(0, 10, 200))
    y = 2.0 + x - 0.3 * x**2
    spl = BSpline(x=x, knots=Knots(count=12), nord=4)
    err, yfit_fit = spl.fit(x, y)
    yfit_val, _ = spl.value(x)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_value_out_of_range_masked():
    x = np.linspace(1, 9, 100)
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y)
    x_eval = np.array([0.0, 5.0, 10.0])
    yfit, mask = spl.value(x_eval)
    assert not mask[0]
    assert mask[1]
    assert not mask[2]


def test_value_unsorted_input():
    rng = np.random.default_rng(11)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.cos(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, y)
    x_unsorted = rng.permutation(x)
    yfit, _ = spl.value(x_unsorted)
    yfit_sorted, _ = spl.value(np.sort(x_unsorted))
    np.testing.assert_allclose(yfit[np.argsort(x_unsorted)], yfit_sorted, atol=1e-12)


def test_value_gap_masking():
    rng = np.random.default_rng(61)
    x = np.sort(rng.uniform(0, 10, 300))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=25), nord=4)
    spl.fit(x, y)
    # Mask three interior breakpoints to create a gap in the knot vector.
    # goodbk becomes [..., 9, 13, ...] so np.diff jumps by 4 > 2.
    spl.bkpt_gpm[10:13] = False
    spl._cached_design = None
    bkpts = spl.breakpoints
    x_gap = np.array([(bkpts[10] + bkpts[11]) / 2.0])
    _, vmask = spl.value(x_gap)
    assert not vmask[0]


# ============================================================================
# BSpline.copy
# ============================================================================

def test_bspline_copy_attributes():
    rng = np.random.default_rng(62)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    cp = spl.copy()
    assert cp.nord == spl.nord
    np.testing.assert_array_equal(cp.breakpoints, spl.breakpoints)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.icoeff, spl.icoeff)


def test_bspline_copy_arrays_are_independent():
    rng = np.random.default_rng(63)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    cp = spl.copy()
    orig_val = cp.coeff[0]
    spl.coeff[0] += 999.0
    assert cp.coeff[0] == orig_val


def test_bspline_copy_cache_is_cleared():
    rng = np.random.default_rng(64)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    assert spl._cached_design is not None
    cp = spl.copy()
    assert cp._cached_design is None


# ============================================================================
# BSpline2D._normalize_basis_x
# ============================================================================

def test_normalize_basis_x_maps_xmin_to_minus_one():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([2.0])), [-1.0])


def test_normalize_basis_x_maps_xmax_to_plus_one():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([8.0])), [1.0])


def test_normalize_basis_x_midpoint_is_zero():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 4.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([2.0])), [0.0])


def test_normalize_basis_x_linear_mapping():
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 10.0
    basis_x = np.linspace(0, 10, 11)
    np.testing.assert_allclose(spl._normalize_basis_x(basis_x), np.linspace(-1, 1, 11), atol=1e-14)


def test_normalize_basis_x_explicit_kwarg_overrides_stored():
    """xmin/xmax kwargs take precedence over self.xmin / self.xmax."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 10.0
    # Explicitly normalise to [2, 8]; points at 2 and 8 should map to ±1
    result = spl._normalize_basis_x(np.array([2.0, 5.0, 8.0]), xmin=2.0, xmax=8.0)
    np.testing.assert_allclose(result, [-1.0, 0.0, 1.0], atol=1e-14)


def test_normalize_basis_x_kwarg_none_falls_back_to_stored():
    """Passing xmin=None / xmax=None uses the stored attributes."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 4.0
    result = spl._normalize_basis_x(np.array([2.0]), xmin=None, xmax=None)
    np.testing.assert_allclose(result, [0.0], atol=1e-14)


def test_normalize_basis_x_no_mutation_with_explicit_kwarg():
    """Explicit xmin/xmax kwargs must not change self.xmin or self.xmax."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 10.0
    spl._normalize_basis_x(np.array([1.0, 5.0, 9.0]), xmin=1.0, xmax=9.0)
    assert spl.xmin == 0.0
    assert spl.xmax == 10.0


def test_normalize_basis_x_no_mutation_when_all_none():
    """When self.xmin/xmax and the kwargs are all None, bounds come from data
    and the stored attributes remain None."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = None
    spl.xmax = None
    basis_x = np.array([3.0, 5.0, 7.0])
    result = spl._normalize_basis_x(basis_x)
    # Data min/max span → endpoints map to ±1
    np.testing.assert_allclose(result[[0, -1]], [-1.0, 1.0], atol=1e-14)
    # Stored attributes must NOT be mutated
    assert spl.xmin is None
    assert spl.xmax is None


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
    basis_xnorm = np.linspace(-1, 1, 10)
    P = spl._poly_basis(basis_xnorm)
    np.testing.assert_allclose(P[:, 0], 1.0, atol=1e-14)
    np.testing.assert_allclose(P[:, 1], basis_xnorm, atol=1e-14)


def test_poly_basis_poly_monomial_structure():
    spl = _make_poly_spl('poly', 4)
    basis_xnorm = np.array([0.0, 0.5, 1.0])
    P = spl._poly_basis(basis_xnorm)
    for k in range(4):
        np.testing.assert_allclose(P[:, k], basis_xnorm**k, atol=1e-14)


def test_poly_basis_unknown_funcname_raises():
    spl = _make_poly_spl('invalid_basis', 2)
    with pytest.raises(ValueError):
        spl._poly_basis(np.array([0.0]))


def test_poly_basis_funcname_kwarg_overrides_stored():
    """funcname kwarg takes precedence over self.funcname."""
    spl = _make_poly_spl('legendre', 3)
    xnorm = np.linspace(-1, 1, 10)
    P_override = spl._poly_basis(xnorm, funcname='chebyshev')
    np.testing.assert_allclose(P_override, fchebyshev(xnorm, 3), atol=1e-14)


def test_poly_basis_npoly_kwarg_changes_column_count():
    """npoly kwarg changes the number of output columns without affecting self.npoly."""
    spl = _make_poly_spl('legendre', 3)
    P = spl._poly_basis(np.linspace(-1, 1, 10), npoly=5)
    assert P.shape == (10, 5)
    assert spl.npoly == 3


def test_poly_basis_funcname_kwarg_does_not_mutate_stored():
    """Passing funcname kwarg must not change self.funcname."""
    spl = _make_poly_spl('legendre', 2)
    spl._poly_basis(np.linspace(-1, 1, 5), funcname='chebyshev')
    assert spl.funcname == 'legendre'


# ============================================================================
# BSpline2D._build_design_matrix
# ============================================================================

def _setup_bspline2d(x, basis_x, npoly, xmin, xmax, funcname, knots, nord):
    """Helper: construct a BSpline2D and populate polynomial-basis state manually."""
    spl = BSpline2D(x=x, knots=knots, nord=nord)
    spl.npoly = npoly
    spl.xmin = xmin
    spl.xmax = xmax
    spl.funcname = funcname
    spl.P = spl._poly_basis(spl._normalize_basis_x(basis_x))
    return spl


def test_bspline2d_build_design_matrix_shape():
    rng = np.random.default_rng(20)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = _setup_bspline2d(x, basis_x, npoly=3, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.shape == (x.size, spl.nord * spl.npoly)


def test_bspline2d_build_design_matrix_outer_product_structure():
    rng = np.random.default_rng(21)
    N = 50
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    npoly = 2
    spl = _setup_bspline2d(x, basis_x, npoly=npoly, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)

    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    B = spl._bspline_basis(x, indx)
    P = spl._poly_basis(spl._normalize_basis_x(basis_x))

    for ii in range(spl.nord):
        for jj in range(npoly):
            col = ii * npoly + jj
            np.testing.assert_allclose(A[:, col], B[:, ii] * P[:, jj], atol=1e-14,
                                       err_msg=f'Outer product mismatch at col {col}')


def test_bspline2d_fit_basis_array_size_mismatch_raises():
    rng = np.random.default_rng(22)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    with pytest.raises(ValueError):
        spl.fit(x, np.sin(x), basis_x=basis_x, basis=np.ones((basis_x.size - 5, 3)))


def test_bspline2d_build_design_matrix_c_order():
    rng = np.random.default_rng(23)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = _setup_bspline2d(x, basis_x, npoly=3, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.flags['C_CONTIGUOUS']


def test_bspline2d_build_design_matrix_P_kwarg_used():
    """Explicit P kwarg changes the design matrix relative to self.P."""
    rng = np.random.default_rng(120)
    N = 60
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = _setup_bspline2d(x, basis_x, npoly=2, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=6), nord=4)
    A_stored, _, _ = spl._build_design_matrix(x)
    P_alt = spl._poly_basis(spl._normalize_basis_x(basis_x * 0.5 + 0.1))  # different basis_x
    A_alt, _, _ = spl._build_design_matrix(x, P=P_alt)
    assert not np.allclose(A_stored, A_alt)


def test_bspline2d_build_design_matrix_P_kwarg_no_mutation():
    """Passing P kwarg must not change self.P."""
    rng = np.random.default_rng(121)
    N = 60
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = _setup_bspline2d(x, basis_x, npoly=2, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=6), nord=4)
    P_orig = spl.P.copy()
    P_alt = rng.standard_normal((N, 2))
    spl._build_design_matrix(x, P=P_alt)
    np.testing.assert_array_equal(spl.P, P_orig)


def test_bspline2d_build_design_matrix_npoly_kwarg_gives_correct_shape():
    """npoly kwarg controls the column count of the returned design matrix."""
    rng = np.random.default_rng(122)
    N = 60
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = _setup_bspline2d(x, basis_x, npoly=2, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=6), nord=4)
    P3 = spl._poly_basis(spl._normalize_basis_x(basis_x), npoly=3)
    A3, _, _ = spl._build_design_matrix(x, P=P3, npoly=3)
    assert A3.shape == (N, spl.nord * 3)


# ============================================================================
# BSpline2D._evaluate_model
# ============================================================================

def test_bspline2d_evaluate_model_matches_explicit_einsum():
    rng = np.random.default_rng(65)
    N = 80
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    npoly = 3
    nord = 4
    spl = _setup_bspline2d(x, basis_x, npoly=npoly, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=nord)
    spl.reset_coeff()
    A, lower, upper = spl._build_design_matrix(x)
    spl.coeff[:] = rng.standard_normal(spl.coeff.shape)
    yfit = spl._evaluate_model(A, lower, upper)
    # Verify one non-empty span manually using the same einsum
    coeffbk = spl.bkpt_gpm[spl.nord:].nonzero()[0]
    goodcoeff = spl.coeff[coeffbk, :]
    k = next(k for k in range(lower.size) if lower[k] <= upper[k])
    sl = slice(lower[k], upper[k] + 1)
    expected = np.einsum(
        'nij,ij->n',
        A[sl, :].reshape(-1, nord, npoly),
        goodcoeff[k:k + nord, :],
    )
    np.testing.assert_allclose(yfit[sl], expected, atol=1e-14)


def test_bspline2d_evaluate_model_coeff_kwarg_differs_from_stored():
    """Explicit coeff kwarg is contracted against the design matrix instead of self.coeff."""
    rng = np.random.default_rng(130)
    N = 80
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    npoly = 3
    spl = _setup_bspline2d(x, basis_x, npoly=npoly, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    spl.reset_coeff()
    spl.coeff[:] = rng.standard_normal(spl.coeff.shape)
    A, lower, upper = spl._build_design_matrix(x)
    coeff_alt = rng.standard_normal(spl.coeff.shape)
    yfit_stored = spl._evaluate_model(A, lower, upper)
    yfit_alt = spl._evaluate_model(A, lower, upper, coeff=coeff_alt)
    assert not np.allclose(yfit_stored, yfit_alt)


def test_bspline2d_evaluate_model_coeff_kwarg_no_mutation():
    """Passing coeff kwarg must not change self.coeff."""
    rng = np.random.default_rng(131)
    N = 80
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    npoly = 3
    spl = _setup_bspline2d(x, basis_x, npoly=npoly, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    spl.reset_coeff()
    spl.coeff[:] = rng.standard_normal(spl.coeff.shape)
    coeff_before = spl.coeff.copy()
    A, lower, upper = spl._build_design_matrix(x)
    spl._evaluate_model(A, lower, upper, coeff=rng.standard_normal(spl.coeff.shape))
    np.testing.assert_array_equal(spl.coeff, coeff_before)


# ============================================================================
# BSpline2D.fit
# ============================================================================

def test_bspline2d_fit_exact_polynomial_legendre():
    """h0(x) + h1(x) * P1(basis_x) with Legendre basis should be recovered exactly."""
    rng = np.random.default_rng(30)
    N = 400
    x = np.sort(rng.uniform(0, 10, N))
    basis_x = rng.uniform(0, 1, N)
    basis_xnorm = 2 * basis_x - 1
    y = 1.0 + 0.5 * x + 0.3 * basis_xnorm
    spl = BSpline2D(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_smooth_function():
    rng = np.random.default_rng(31)
    N = 500
    x = np.sort(rng.uniform(0, 2 * np.pi, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    spl = BSpline2D(x=x, knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=3, xmin=0.0, xmax=1.0)
    assert err == 0
    assert np.std(yfit - y) < 2e-2


def test_bspline2d_fit_exact_polynomial_chebyshev():
    rng = np.random.default_rng(32)
    N = 400
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(-1, 1, N)
    y = 2.0 + basis_x  # T0 + T1 in Chebyshev basis
    spl = BSpline2D(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='chebyshev', npoly=2, xmin=-1.0, xmax=1.0)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_exact_polynomial_poly():
    rng = np.random.default_rng(33)
    N = 400
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 2, N)
    basis_xnorm = 2.0 * (basis_x - 1.0) / 2.0 - 0.0   # xmin=0, xmax=2 → norm range [-1, 1]
    y = 1.0 + basis_xnorm
    spl = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='poly', npoly=2, xmin=0.0, xmax=2.0)
    assert err == 0
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_reset_knots_rebuilds_breakpoints():
    rng = np.random.default_rng(52)
    N = 200
    x1 = np.sort(rng.uniform(0, 5, N))
    basis_x_1 = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1), basis_x=basis_x_1, npoly=2, xmin=0.0, xmax=1.0)
    bkpt_before = spl.breakpoints.copy()
    x3 = np.sort(rng.uniform(0, 10, N))
    basis_x_2 = rng.uniform(0, 1, N)
    spl.fit(x3, np.sin(x3), basis_x=basis_x_2, npoly=2, xmin=0.0, xmax=1.0, reset_knots=True)
    assert spl.breakpoints.min() <= x3.min()
    assert spl.breakpoints.max() >= x3.max()
    assert not np.array_equal(spl.breakpoints, bkpt_before)


def test_bspline2d_fit_reset_knots_resets_mask():
    rng = np.random.default_rng(53)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x)
    spl = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    spl.bkpt_gpm[5] = False
    assert not spl.bkpt_gpm.all()
    spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0, reset_knots=True)
    assert spl.bkpt_gpm.all()


def test_bspline2d_fit_returns_minus1_on_degenerate_cholesky():
    rng = np.random.default_rng(66)
    x = np.sort(rng.uniform(0, 2, 150))
    basis_x = rng.uniform(0, 1, 150)
    y = np.sin(x)
    full_bkpt = np.linspace(0, 10, 40)
    spl = BSpline2D(knots=Knots(full=full_bkpt), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    assert err == -1
    assert yfit.shape == x.shape


def test_bspline2d_fit_ivar_none_matches_unit_weights():
    rng = np.random.default_rng(71)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * basis_x
    spl_w = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl_n = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    err_w, yfit_w = spl_w.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0, ivar=np.ones(N))
    err_n, yfit_n = spl_n.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    assert err_w == 0 and err_n == 0
    np.testing.assert_array_equal(yfit_w, yfit_n)


def test_bspline2d_fit_nonuniform_ivar_affects_result():
    rng = np.random.default_rng(73)
    N = 400
    x = np.sort(rng.uniform(0.5, 5.0, N))
    basis_x = rng.uniform(0, 1, N)
    ivar = 1.0 / x**2  # noise amplitude proportional to x
    y = np.sin(x) + 0.1 * basis_x + rng.standard_normal(N) / np.sqrt(ivar)
    spl_u = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    spl_n = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err_u, yfit_u = spl_u.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    err_n, yfit_n = spl_n.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0, ivar=ivar)
    assert err_u == 0 and err_n == 0
    assert np.isfinite(yfit_n).all()
    assert not np.allclose(yfit_u, yfit_n)


# ============================================================================
# BSpline2D.value
# ============================================================================

def test_bspline2d_value_matches_fit_at_training_points():
    rng = np.random.default_rng(40)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * basis_x
    spl = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    err, yfit_fit = spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    yfit_val, _ = spl.value(x, basis_x=basis_x)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_bspline2d_value_basis_x_required():
    rng = np.random.default_rng(41)
    x = np.sort(rng.uniform(0, 5, 50))
    basis_x = rng.uniform(0, 1, 50)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    with pytest.raises(ValueError):
        spl.value(x)  # basis_x=None not permitted when funcname is set


def test_bspline2d_value_gap_masking():
    rng = np.random.default_rng(67)
    N = 300
    x = np.sort(rng.uniform(0, 10, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x)
    spl = BSpline2D(x=x, knots=Knots(count=25), nord=4)
    spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    spl.bkpt_gpm[10:13] = False
    spl._cached_design = None
    bkpts = spl.breakpoints
    x_gap = np.array([(bkpts[10] + bkpts[11]) / 2.0])
    basis_x_gap = np.array([0.5])
    _, vmask = spl.value(x_gap, basis_x=basis_x_gap)
    assert not vmask[0]


# --- value() override kwarg tests ---

def test_bspline2d_value_string_basis_matches_stored_funcname():
    """Passing basis='legendre' explicitly gives the same result as basis=None."""
    rng = np.random.default_rng(140)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * basis_x
    spl = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    yfit_none, _ = spl.value(x, basis_x=basis_x)
    yfit_str, _ = spl.value(x, basis='legendre', basis_x=basis_x)
    np.testing.assert_allclose(yfit_none, yfit_str, atol=1e-12)


def test_bspline2d_value_xmin_xmax_kwarg_affects_normalization():
    """Supplying different xmin/xmax bounds changes the evaluated polynomial basis."""
    rng = np.random.default_rng(141)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    # Strong unit-amplitude Legendre component guarantees coeff[:,1] ~ 1.
    x_norm = 2.0 * basis_x - 1.0
    y = np.sin(x) + x_norm
    spl = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    yfit_stored, _ = spl.value(x, basis='legendre', basis_x=basis_x)
    # xmax=2.0 maps x_norm_alt = basis_x - 1, which differs from x_norm by basis_x
    yfit_alt, _ = spl.value(x, basis='legendre', basis_x=basis_x, xmin=0.0, xmax=2.0)
    assert not np.allclose(yfit_stored, yfit_alt)


def test_bspline2d_value_coeff_kwarg_additivity():
    """Partial coeff evaluations sum to the full model (linearity)."""
    rng = np.random.default_rng(142)
    N = 300
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.2 * basis_x + 0.1 * basis_x**2
    spl = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=3, xmin=0.0, xmax=1.0)
    P = spl.P  # shape (N, 3)
    yfit_full, _ = spl.value(x, basis_x=basis_x)
    yfit_part0, _ = spl.value(x, basis=P[:, :1], coeff=spl.coeff[:, :1])
    yfit_part1, _ = spl.value(x, basis=P[:, 1:], coeff=spl.coeff[:, 1:])
    np.testing.assert_allclose(yfit_part0 + yfit_part1, yfit_full, atol=1e-12)


def test_bspline2d_value_overrides_do_not_mutate_instance():
    """Calling value() with override kwargs must not change any stored attributes."""
    rng = np.random.default_rng(143)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.2 * basis_x
    spl = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)

    npoly_before = spl.npoly
    funcname_before = spl.funcname
    xmin_before = spl.xmin
    xmax_before = spl.xmax
    coeff_before = spl.coeff.copy()
    P_before = spl.P.copy()

    coeff_alt = rng.standard_normal(spl.coeff.shape)
    spl.value(x, basis='chebyshev', basis_x=basis_x, xmin=0.1, xmax=0.9,
              coeff=coeff_alt)

    assert spl.npoly == npoly_before
    assert spl.funcname == funcname_before
    assert spl.xmin == xmin_before
    assert spl.xmax == xmax_before
    np.testing.assert_array_equal(spl.coeff, coeff_before)
    np.testing.assert_array_equal(spl.P, P_before)


# ============================================================================
# BSpline2D.copy
# ============================================================================

def test_bspline2d_copy_attributes():
    rng = np.random.default_rng(68)
    N = 100
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, basis='chebyshev', npoly=3, xmin=0.0, xmax=1.0)
    cp = spl.copy()
    assert cp.npoly == spl.npoly
    assert cp.xmin == spl.xmin
    assert cp.xmax == spl.xmax
    assert cp.funcname == spl.funcname
    assert cp.nord == spl.nord
    assert cp._cached_design is None
    assert cp.basis_x is spl.basis_x
    np.testing.assert_array_equal(cp.breakpoints, spl.breakpoints)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bspline2d_copy_coeff_independent():
    rng = np.random.default_rng(69)
    N = 100
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, npoly=3, xmin=0.0, xmax=1.0)
    cp = spl.copy()
    orig_val = cp.coeff[0, 0]
    spl.coeff[0, 0] += 999.0
    assert cp.coeff[0, 0] == orig_val


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
    new.bkpt_gpm[:] = old.mask
    new._cached_design = None
    new_err, new_yfit = new.fit(x, y, ivar=invvar)

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
    new.bkpt_gpm[:] = old.mask
    new._cached_design = None
    new.fit(x, y, ivar=invvar)

    np.testing.assert_allclose(new.coeff, old.coeff, rtol=1e-10)


def test_crosscheck_2d_fit_yfit_match():
    rng = np.random.default_rng(52)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.3 * basis_x)
    invvar = np.ones(N)
    nord = 4
    npoly = 2
    nbkpts = 12

    old = bspline(x=x, nord=nord, npoly=npoly, nbkpts=nbkpts, funcname='legendre')
    old.xmin = 0.0
    old.xmax = 1.0
    old_err, old_yfit = old.fit(x, y, invvar, x2=basis_x)

    new = BSpline2D(knots=old.breakpoints.copy(), nord=nord)
    new.bkpt_gpm[:] = old.mask
    new._cached_design = None
    new_err, new_yfit = new.fit(x, y, basis_x=basis_x, ivar=invvar,
                                basis='legendre', npoly=npoly, xmin=0.0, xmax=1.0)

    assert old_err == 0 and new_err == 0
    np.testing.assert_allclose(new_yfit, old_yfit, rtol=1e-8)


def test_crosscheck_2d_coeff_transposed():
    """new.coeff should equal old.coeff.T (shape change from (npoly,nc) to (nc,npoly))."""
    rng = np.random.default_rng(53)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.3 * basis_x)
    invvar = np.ones(N)

    old = bspline(x=x, nord=4, npoly=2, nbkpts=12, funcname='legendre')
    old.xmin = 0.0
    old.xmax = 1.0
    old.fit(x, y, invvar, x2=basis_x)

    new = BSpline2D(knots=old.breakpoints.copy(), nord=4)
    new.bkpt_gpm[:] = old.mask
    new._cached_design = None
    new.fit(x, y, basis_x=basis_x, ivar=invvar, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)

    np.testing.assert_allclose(new.coeff, old.coeff.T, rtol=1e-8)


def test_crosscheck_1d_value_match():
    rng = np.random.default_rng(54)
    x = np.sort(rng.uniform(0, 10, 150))
    y = 1 + 2*x - 0.1*x**2
    invvar = np.ones_like(x)

    old = bspline(x=x, nord=4, npoly=1, nbkpts=12)
    old.fit(x, y, invvar)

    new = BSpline(knots=old.breakpoints.copy(), nord=4)
    new.bkpt_gpm[:] = old.mask
    new._cached_design = None
    new.fit(x, y, ivar=invvar)

    x_eval = np.sort(rng.uniform(0, 10, 50))
    old_val, _ = old.value(x_eval)
    new_val, _ = new.value(x_eval)
    np.testing.assert_allclose(new_val, old_val, rtol=1e-10)


# ============================================================================
# BSpline2D.fit / value — array basis path
# ============================================================================

def test_bspline2d_fit_array_basis_matches_string_basis():
    """Passing a pre-built P array to fit() should give identical results to the string path."""
    rng = np.random.default_rng(80)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.2 * basis_x)

    spl_str = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err_s, yfit_s = spl_str.fit(x, y, basis_x=basis_x, basis='legendre', npoly=3, xmin=0.0, xmax=1.0)

    # Build the same P array that the string path would have produced
    P = spl_str._poly_basis(spl_str._normalize_basis_x(basis_x))

    spl_arr = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err_a, yfit_a = spl_arr.fit(x, y, basis_x=basis_x, basis=P)

    assert err_s == 0 and err_a == 0
    np.testing.assert_allclose(yfit_a, yfit_s, atol=1e-12)


def test_bspline2d_value_array_basis_raises_when_basis_none():
    """value() with no basis argument must raise when fit was done with an array basis
    and new evaluation points are requested."""
    rng = np.random.default_rng(81)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    # Build P directly so we don't need the object's internal state
    P = np.column_stack([np.ones(N), basis_x])  # shape (N, 2)
    spl = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, basis=P)
    # The fast path fires for the training arrays; use new arrays to reach the
    # funcname is None check.
    x_new = np.sort(rng.uniform(0, 5, N))
    basis_x_new = rng.uniform(0, 1, N)
    with pytest.raises(ValueError):
        spl.value(x_new, basis_x=basis_x_new)  # funcname is None → must raise


def test_bspline2d_value_array_basis_with_explicit_basis():
    """value() should succeed when a basis array is passed after an array-basis fit."""
    rng = np.random.default_rng(82)
    N = 300
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    P = np.column_stack([np.ones(N), basis_x])  # shape (N, 2)
    spl = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err, yfit_fit = spl.fit(x, y, basis_x=basis_x, basis=P)
    yfit_val, _ = spl.value(x, basis_x=basis_x, basis=P)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_bspline2d_fit_basis_array_1d_raises():
    """fit() must reject a 1-D basis array."""
    rng = np.random.default_rng(83)
    x = np.sort(rng.uniform(0, 5, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    with pytest.raises(ValueError):
        spl.fit(x, np.sin(x), basis_x=basis_x, basis=np.ones(100))  # 1-D, not 2-D


def test_bspline2d_value_restores_training_P_after_call():
    """value() must not clobber self.P with the evaluation basis."""
    rng = np.random.default_rng(84)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * basis_x

    spl = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    P_before = spl.P.copy()

    x_eval = np.sort(rng.uniform(0, 5, 50))
    basis_x_eval = rng.uniform(0, 1, 50)
    spl.value(x_eval, basis_x=basis_x_eval)

    np.testing.assert_array_equal(spl.P, P_before)


# ============================================================================
# BSpline.value / BSpline2D.value — interpolate keyword
# ============================================================================

def test_value_interpolate_equals_numpy_interp():
    """interpolate=True at new x must equal np.interp(x, self.x, self.yfit)."""
    rng = np.random.default_rng(90)
    x = np.sort(rng.uniform(0, 10, 200))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl.fit(x, y)

    x_eval = np.sort(rng.uniform(0.5, 9.5, 150))
    yfit_interp, _ = spl.value(x_eval, interpolate=True)

    np.testing.assert_array_equal(yfit_interp, np.interp(x_eval, spl.x, spl.yfit))


def test_value_interpolate_gpm_matches_full():
    """interpolate=True and interpolate=False must return identical gpm."""
    rng = np.random.default_rng(91)
    x = np.sort(rng.uniform(0, 10, 200))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl.fit(x, y)

    x_eval = np.sort(rng.uniform(0, 10, 100))
    _, gpm_full = spl.value(x_eval, interpolate=False)
    _, gpm_interp = spl.value(x_eval, interpolate=True)

    np.testing.assert_array_equal(gpm_full, gpm_interp)


def test_value_interpolate_not_identical_to_full():
    """interpolate=True approximates but does not reproduce full B-spline evaluation.

    With a well-sampled training grid (N=500, h≈0.02) the max linear-interpolation
    error is bounded above by 5e-3 (well below 1% of the signal amplitude), but
    exceeds 1e-4 near knot locations where the piecewise-polynomial model has
    non-negligible curvature.
    """
    rng = np.random.default_rng(92)
    x = np.sort(rng.uniform(0, 10, 500))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=50), nord=4)
    spl.fit(x, y)

    x_eval = np.sort(rng.uniform(0.5, 9.5, 80))
    yfit_full, gpm  = spl.value(x_eval, interpolate=False)
    yfit_interp, _  = spl.value(x_eval, interpolate=True)

    max_diff = np.max(np.abs(yfit_interp - yfit_full)[gpm])
    assert max_diff > 1e-4   # not identical — linear interpolation error is nonzero
    assert max_diff < 5e-3   # good approximation for a well-sampled model (h≈0.02)


def test_value_interpolate_training_x_fast_path():
    """When x is self.x the identity fast path fires before the interpolate check."""
    rng = np.random.default_rng(93)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.cos(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    err, yfit_fit = spl.fit(x, y)

    yfit_val, _ = spl.value(x, interpolate=True)
    assert yfit_val is yfit_fit


# ---------------------------------------------------------------------------
# bspline_profile_refactor — integration tests against reference data
# ---------------------------------------------------------------------------

def test_bspline_profile_refactor_spec():
    """
    bspline_profile_refactor reproduces the spectral flat-field reference fit
    from the gemini_gnirs_32 test data files.
    """
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit))
             for slit in [0, 1]]
    logrej = 0.5
    spec_samp_fine = 1.2
    for f in files:
        d = np.load(f)
        sset, _, spec_flat_fit, _, exit_status = bspline_profile_refactor(
            d['spec_coo_data'], d['spec_flat_data'],
            ivar=d['spec_ivar_data'], gpm=d['spec_gpm_data'],
            nord=4, upper=logrej, lower=logrej,
            kwargs_knots={'spacing': spec_samp_fine},
            kwargs_reject={'groupbadpix': True, 'maxrej': 5},
        )
        assert isinstance(sset, BSpline)
        assert np.allclose(d['spec_flat_fit'], spec_flat_fit), \
            'Bad spectral bspline_profile_refactor result'


def test_bspline_profile_refactor_spat():
    """
    bspline_profile_refactor reproduces the spatial flat-field reference fit
    from the gemini_gnirs_32 test data files.
    """
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spat_fit.npz'.format(slit))
             for slit in [0, 1]]
    for f in files:
        d = np.load(f)
        bkspace = np.fmax(
            1.0 / d['median_slit_width'] / 10.0,
            1.2 * np.median(np.diff(d['spat_coo_data'])),
        )
        sset, _, spat_flat_fit, _, exit_status = bspline_profile_refactor(
            d['spat_coo_data'], d['spat_flat_data'],
            nord=4, upper=5.0, lower=5.0,
            kwargs_knots={'spacing': bkspace},
        )
        assert isinstance(sset, BSpline)
        assert np.allclose(d['spat_flat_fit'], spat_flat_fit), \
            'Bad spatial bspline_profile_refactor result'


def test_bspline_profile_refactor_twod():
    """
    bspline_profile_refactor reproduces the 2D flat-field reference fit
    from the gemini_gnirs_32 test data files.
    """
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_twod_fit.npz'.format(slit))
             for slit in [0, 1]]
    spec_samp_coarse = 50.0
    twod_sigrej = 4.0
    for f in files:
        d = np.load(f)
        sset, _, twod_flat_fit, _, exit_status = bspline_profile_refactor(
            d['twod_spec_coo_data'], d['twod_flat_data'],
            basis=d['poly_basis'],
            ivar=d['twod_ivar_data'], gpm=d['twod_gpm_data'],
            nord=4, upper=twod_sigrej, lower=twod_sigrej,
            kwargs_knots={'spacing': spec_samp_coarse},
            kwargs_reject={'groupbadpix': True, 'maxrej': 10},
        )
        assert isinstance(sset, BSpline2D)
        assert np.allclose(d['twod_flat_fit'], twod_flat_fit), \
            'Bad 2D bspline_profile_refactor result'


# ============================================================================
# bspline_profile_refactor option coverage and agreement tests
# ============================================================================

def test_bspline_profile_refactor_string_basis():
    """basis='legendre' + basis_x activates the 2D BSpline2D path."""
    rng = np.random.default_rng(42)
    n = 300
    x = np.sort(rng.uniform(0, 10, n))
    basis_x = rng.uniform(-1, 1, n)
    npoly = 3
    y_true = (1.0 + 0.3 * basis_x) * (np.sin(x / 2) + 2.0)
    y = y_true + rng.normal(0, 0.05, n)
    ivar = np.full(n, 400.0)

    sset, outmask, yfit, reduced_chi, exit_status = bspline_profile_refactor(
        x, y, ivar=ivar, nord=4, npoly=npoly,
        basis='legendre', basis_x=basis_x, xmin=-1.0, xmax=1.0,
        upper=5.0, lower=5.0,
        kwargs_knots={'spacing': 1.0},
    )
    assert isinstance(sset, BSpline2D)
    assert exit_status == 0
    assert np.allclose(y_true, yfit, atol=0.15)


def test_bspline_profile_refactor_basis_array_shapes():
    """Flat 1D basis array (auto-reshaped) gives identical yfit to its 2D form."""
    rng = np.random.default_rng(99)
    n = 200
    x = np.sort(rng.uniform(0, 5, n))
    x2 = rng.uniform(-1, 1, n)
    npoly = 2
    basis_2d = flegendre(x2, npoly)
    basis_1d = basis_2d.flatten()
    y = np.sin(x) + 1.5 * x2 + rng.normal(0, 0.1, n)
    ivar = np.full(n, 100.0)

    _, _, yfit_2d, _, _ = bspline_profile_refactor(
        x, y, ivar=ivar, basis=basis_2d, kwargs_knots={'spacing': 0.5},
    )
    _, _, yfit_1d, _, _ = bspline_profile_refactor(
        x, y, ivar=ivar, basis=basis_1d, kwargs_knots={'spacing': 0.5},
    )
    assert np.allclose(yfit_2d, yfit_1d)


def test_bspline_profile_refactor_exit_statuses():
    """Exit status 4 (too few points) and 1 (maxiter exceeded)."""
    rng = np.random.default_rng(0)
    n = 300
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x) + rng.normal(0, 0.1, n)
    ivar = np.full(n, 100.0)

    # exit_status = 4: only 3 good points, fewer than nord=4
    gpm_few = np.zeros(n, dtype=bool)
    gpm_few[:3] = True
    _, _, _, _, es = bspline_profile_refactor(
        x, y, ivar=ivar, gpm=gpm_few, nord=4,
        kwargs_knots={'spacing': 1.0},
    )
    assert es == 4

    # exit_status = 1: maxiter=1 with outliers that need exactly one rejection
    # pass (converges in 2 iterations, but maxiter caps at 1)
    y_out = y.copy()
    y_out[rng.choice(n, 15, replace=False)] += 20.0
    _, _, _, _, es = bspline_profile_refactor(
        x, y_out, ivar=ivar, nord=4, maxiter=1, upper=3.0, lower=3.0,
        kwargs_knots={'spacing': 1.0},
    )
    assert es == 1


def test_bspline_profile_refactor_matches_bspline_profile_1d():
    """1D path (basis=None) agrees with bspline_profile(profile_basis=ones)."""
    rng = np.random.default_rng(77)
    n = 400
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x / 2) + 2.0 + rng.normal(0, 0.05, n)
    ivar = np.full(n, 400.0)
    profile_basis = np.ones((n, 1))
    spacing = 1.0

    _, _, yfit_old, _, _ = bspline_profile(
        x, y, ivar, profile_basis,
        upper=5, lower=5, nord=4,
        kwargs_bspline={'bkspace': spacing},
        quiet=True,
    )
    _, _, yfit_new, _, _ = bspline_profile_refactor(
        x, y, ivar=ivar, nord=4, upper=5, lower=5,
        kwargs_knots={'spacing': spacing},
    )
    assert np.allclose(yfit_old, yfit_new)


def test_bspline_profile_refactor_matches_bspline_profile_2d():
    """2D path agrees with bspline_profile when given the same pre-built basis."""
    rng = np.random.default_rng(13)
    n = 400
    x = np.sort(rng.uniform(0, 10, n))
    x2 = rng.uniform(-1, 1, n)
    npoly = 3
    profile_basis = flegendre(x2, npoly)
    y = (1.0 + 0.5 * x2) * (np.sin(x / 2) + 2.0) + rng.normal(0, 0.05, n)
    ivar = np.full(n, 400.0)
    spacing = 1.0

    _, _, yfit_old, _, _ = bspline_profile(
        x, y, ivar, profile_basis,
        upper=5, lower=5, nord=4,
        kwargs_bspline={'bkspace': spacing},
        quiet=True,
    )
    _, _, yfit_new, _, _ = bspline_profile_refactor(
        x, y, ivar=ivar, basis=profile_basis, nord=4, upper=5, lower=5,
        kwargs_knots={'spacing': spacing},
    )
    assert np.allclose(yfit_old, yfit_new)
