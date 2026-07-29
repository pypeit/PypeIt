"""
Per-method unit tests for :class:`~pypeit.core.bspline.BSpline`,
:class:`~pypeit.core.bspline.BSpline2D`, and :class:`~pypeit.core.bspline.Knots`.

Each test function targets a single method or a tightly coupled pair.  The
integration tests at the end verify that
:func:`~pypeit.core.fitting.iterative_bspline_fit` reproduces reference fits from
the PypeIt development-suite test data.

All random-number generators use a fixed seed for reproducibility.
"""

import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.axes
import matplotlib.pyplot

from IPython import embed
import numpy as np
import pytest

from pypeit import dataPaths
from pypeit.core.bspline import BSpline, BSpline2D, Knots, _cholesky_banded
from pypeit.core.fitting import iterative_bspline_fit
from pypeit.core.fitting import bspline_qa
from pypeit.core.basis import fchebyshev, flegendre


# ============================================================================
# Knots.build
# ============================================================================

def test_build_breakpoints_fullbkpt_returned_sorted():
    """full= knot vector is sorted on return."""
    rng = np.random.default_rng(0)
    pts = rng.uniform(0, 10, 20)
    assert np.all(np.diff(Knots(full=pts, nord=4).breakpoints) >= 0), 'full= breakpoints must be sorted ascending'


def test_build_breakpoints_fullbkpt_padded_when_short():
    """Short full= vector is padded to at least 2*nord entries."""
    pts = np.array([0.0, 1.0, 2.0, 3.0])  # length 4 < 2*nord=8
    assert Knots(full=pts, nord=4).breakpoints.size >= 2 * 4, 'short full= vector must be padded to at least 2*nord'


def test_build_breakpoints_bkspace_strategy():
    """spacing= strategy spans the full x range."""
    x = np.linspace(0, 10, 200)
    kv = Knots(spacing=1.0, x=x, nord=4).breakpoints
    assert kv.min() <= 0.0, 'breakpoint minimum must not exceed x.min()'
    assert kv.max() >= 10.0, 'breakpoint maximum must reach x.max()'


def test_build_breakpoints_nbkpts_strategy():
    """count= strategy produces a valid padded knot vector."""
    x = np.linspace(0, 5, 100)
    assert Knots(count=10, x=x, nord=4).breakpoints.size >= 2 * 4, 'count= must produce at least 2*nord breakpoints'


def test_build_breakpoints_everyn_strategy():
    """stride= strategy produces more than 2*nord breakpoints."""
    x = np.linspace(0, 10, 300)
    assert Knots(stride=20, x=x, nord=4).breakpoints.size > 2 * 4, 'stride= must produce more than 2*nord breakpoints'


def test_build_breakpoints_bkpt_strategy():
    """interior= strategy places knots at the x endpoints."""
    x = np.linspace(0, 10, 200)
    interior = np.array([2.0, 4.0, 6.0, 8.0])
    kv = Knots(interior=interior, x=x, nord=4).breakpoints
    assert kv[0] <= 0.0, 'first breakpoint must not exceed x.min()'
    assert kv[-1] >= 10.0, 'last breakpoint must reach x.max()'


def test_build_breakpoints_phantom_knots_at_each_end():
    """nord-1 phantom knots are placed before the first and after the last interior breakpoint."""
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=6, x=x, nord=nord).breakpoints
    interior_min = kv[nord - 1]
    interior_max = kv[-(nord)]
    assert np.all(kv[:nord - 1] < interior_min), 'phantom knots must lie before the first interior breakpoint'
    assert np.all(kv[-(nord - 1):] > interior_max), 'phantom knots must lie after the last interior breakpoint'


def test_build_breakpoints_raises_without_x_or_fullbkpt():
    """build() raises ValueError when neither x nor full= is provided."""
    with pytest.raises(ValueError):
        Knots(spacing=1.0).build(None, 4)


def test_build_breakpoints_spacing_too_large_fallback():
    """spacing larger than the x range falls back to two-endpoint breakpoints."""
    x = np.linspace(0, 5, 100)
    kv = Knots(spacing=10.0, x=x, nord=4).breakpoints
    assert kv is not None, 'large spacing must still produce a knot vector'
    assert kv.min() <= 0.0, 'fallback breakpoints must reach x.min()'
    assert kv.max() >= 5.0, 'fallback breakpoints must reach x.max()'


def test_build_breakpoints_stride_too_large_fallback():
    """stride >= x.size falls back to two-endpoint breakpoints."""
    x = np.linspace(0, 5, 50)
    kv = Knots(stride=100, x=x, nord=4).breakpoints
    assert kv is not None, 'large stride must still produce a knot vector'
    assert kv.min() <= 0.0, 'fallback breakpoints must reach x.min()'
    assert kv.max() >= 5.0, 'fallback breakpoints must reach x.max()'


def test_build_breakpoints_interior_outside_range_fallback():
    """interior= points outside the x range fall back to two-endpoint breakpoints."""
    x = np.linspace(0, 5, 100)
    interior = np.array([7.0, 8.0, 9.0])
    kv = Knots(interior=interior, x=x, nord=4).breakpoints
    assert kv is not None, 'out-of-range interior must still produce a knot vector'
    assert kv.min() <= 0.0, 'fallback breakpoints must reach x.min()'
    assert kv.max() >= 5.0, 'fallback breakpoints must reach x.max()'


def test_build_breakpoints_no_strategy_raises_valueerror():
    """build() raises ValueError when stride=None and no other strategy is given."""
    with pytest.raises(ValueError):
        Knots(stride=None).build(np.linspace(0, 5, 50), 4)


def test_knots_init_nord_without_x_warns():
    """Knots() with nord but no x emits a UserWarning."""
    with pytest.warns(UserWarning):
        Knots(count=10, nord=4)


def test_build_breakpoints_spread_doubles_phantom_spacing():
    """spread=2.0 must double the phantom-knot step at each end relative to spread=1.0."""
    x = np.arange(1000, dtype=float)
    k1 = Knots(interior=x[[0, -1]], x=x, nord=4)           # spread=1.0 (default)
    k2 = Knots(interior=x[[0, -1]], x=x, nord=4, spread=2.0)
    np.testing.assert_array_equal(
        np.diff(k2.breakpoints)[[0, -1]],
        2 * np.diff(k1.breakpoints)[[0, -1]],
    )


def test_build_breakpoints_spacing_count_equivalence():
    """spacing=10. and count=100 produce identical breakpoints for x = arange(1000)."""
    x = np.arange(1000, dtype=float)
    k_spacing = Knots(spacing=10., x=x, nord=4)
    k_count = Knots(count=100, x=x, nord=4)
    np.testing.assert_array_equal(k_spacing.breakpoints, k_count.breakpoints)


def test_build_breakpoints_single_interior_matches_count_one():
    """A single interior point collapses to [sx, ex], giving the same breakpoints as count=1."""
    x = np.arange(1000, dtype=float)
    k_interior = Knots(interior=np.array([500.]), x=x, nord=4)
    k_count = Knots(count=1, x=x, nord=4)
    np.testing.assert_array_equal(k_interior.breakpoints, k_count.breakpoints)


def test_build_breakpoints_stride_too_large_matches_count_one():
    """stride >= x.size falls back to [sx, ex], giving the same breakpoints as count=1."""
    x = np.arange(1000, dtype=float)
    k_stride = Knots(stride=x.size + 1, x=x, nord=4)
    k_count = Knots(count=1, x=x, nord=4)
    np.testing.assert_array_equal(k_stride.breakpoints, k_count.breakpoints)


def test_build_breakpoints_float_stride_uniform_spacing():
    """A float stride on a uniform grid produces uniformly spaced breakpoints."""
    x = np.arange(1000, dtype=float)
    k = Knots(stride=1.5, x=x, nord=4)
    diffs = np.diff(k.breakpoints)
    assert np.all(diffs == 1.5), (
        f'Expected all diffs == 1.5; got unique values {np.unique(diffs)}'
    )


def test_build_breakpoints_stride_irregular_grid():
    """stride on an irregular grid produces non-uniform breakpoint spacing."""
    x = np.arange(1000, dtype=float)
    rng = np.random.default_rng(99)
    x_irr = np.sort(x + 5. * rng.standard_normal(x.size))
    k = Knots(stride=10, x=x_irr, nord=4)
    assert np.std(np.diff(k.breakpoints)) > 1., (
        'Expected irregular breakpoint spacing for an irregular-grid stride'
    )


# ============================================================================
# Knots.copy
# ============================================================================

def test_knots_copy_strategy_params_preserved():
    """copy() preserves all strategy parameters."""
    orig = Knots(count=10, spread=1.5, stride=None)
    cp = orig.copy()
    assert cp.count == orig.count, 'copy must preserve count'
    assert cp.spread == orig.spread, 'copy must preserve spread'
    assert cp.stride == orig.stride, 'copy must preserve stride'
    assert cp.spacing == orig.spacing, 'copy must preserve spacing'


def test_knots_copy_breakpoints_deep_copied():
    """copy() deep-copies the breakpoint array so mutations do not propagate."""
    orig = Knots(count=10, x=np.linspace(0, 5, 100), nord=4)
    cp = orig.copy()
    assert np.array_equal(cp.breakpoints, orig.breakpoints), 'copy must duplicate breakpoint values'
    orig_val = orig.breakpoints[0]
    cp._breakpoints[0] = 999.0
    assert orig.breakpoints[0] == orig_val, 'mutating the copy must not change the original'


def test_knots_copy_unbuilt_returns_none_breakpoints():
    """copy() of an unbuilt Knots returns None breakpoints."""
    orig = Knots(count=10)
    cp = orig.copy()
    assert cp.breakpoints is None, 'copy of unbuilt Knots must have None breakpoints'


# ============================================================================
# BSpline.__init__ — knots type validation
# ============================================================================

def test_bspline_init_invalid_knots_type_raises():
    """BSpline() raises TypeError for an unsupported knots type."""
    with pytest.raises(TypeError):
        BSpline(knots='not_a_valid_knots_type')


# ============================================================================
# BSpline._uniq
# ============================================================================

def test_uniq_empty_array_raises():
    """_uniq() raises ValueError on an empty array."""
    with pytest.raises(ValueError):
        BSpline._uniq(np.array([]))


# ============================================================================
# BSpline._find_spans
# ============================================================================

def test_find_spans_all_indices_in_range():
    """_find_spans() returns indices clamped to [nord-1, n-1]."""
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    n = kv.size - nord
    indx = BSpline._find_spans(x, kv, nord)
    assert np.all(indx >= nord - 1), 'all span indices must be >= nord-1'
    assert np.all(indx <= n - 1), 'all span indices must be <= n-1'


def test_find_spans_bracketing():
    """_find_spans() returns the index of the enclosing knot span."""
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    indx = BSpline._find_spans(x, kv, nord)
    for i, (xi, il) in enumerate(zip(x, indx)):
        assert kv[il] <= xi, f'breakpoint[{il}] must be <= x[{i}]'
        if il < kv.size - 1:
            assert xi <= kv[il + 1] + 1e-12, f'x[{i}] must be <= breakpoint[{il+1}]'


def test_find_spans_clamped_at_lower_edge():
    """_find_spans() clamps x at or below the first interior knot to nord-1."""
    kv = np.array([0.0, 1.0, 2.0, 3.0, 4.0, 5.0])
    nord = 2
    x = np.array([0.0])
    indx = BSpline._find_spans(x, kv, nord)
    assert indx[0] >= nord - 1, 'lower-edge index must be clamped to nord-1'


def test_find_spans_monotone_input_gives_monotone_output():
    """_find_spans() produces non-decreasing span indices for sorted x."""
    x = np.linspace(0, 10, 100)
    nord = 4
    kv = Knots(count=8, x=x, nord=nord).breakpoints
    indx = BSpline._find_spans(x, kv, nord)
    assert np.all(np.diff(indx) >= 0), 'span indices must be non-decreasing for sorted x'


# ============================================================================
# BSpline._bspline_basis
# ============================================================================

def test_bspline_basis_output_shape():
    """_bspline_basis() returns shape (N, nord)."""
    x = np.sort(np.random.default_rng(7).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.shape == (x.size, spl.nord), f'expected basis shape {(x.size, spl.nord)}'


def test_bspline_basis_c_order():
    """_bspline_basis() returns a C-contiguous array."""
    x = np.sort(np.random.default_rng(8).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert basis.flags['C_CONTIGUOUS'], 'basis array must be C-contiguous'


def test_bspline_basis_partition_of_unity():
    """B-spline basis functions sum to 1 at every x."""
    x = np.sort(np.random.default_rng(9).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    np.testing.assert_allclose(basis.sum(axis=1), 1.0, atol=1e-12)


def test_bspline_basis_non_negative():
    """B-spline basis functions are non-negative everywhere."""
    x = np.sort(np.random.default_rng(10).uniform(0, 10, 50))
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    basis = spl._bspline_basis(x, indx)
    assert np.all(basis >= -1e-14), 'B-spline basis functions must be non-negative'


def test_bspline_basis_linear_case():
    """Partition-of-unity holds for a linear (nord=2) B-spline."""
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
    """_build_design_matrix() returns shape (N, nord)."""
    x = np.sort(np.random.default_rng(42).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.shape == (x.size, spl.nord), f'expected design matrix shape {(x.size, spl.nord)}'


def test_build_design_matrix_lower_upper_lengths():
    """lower and upper index vectors have length (nactive - nord + 1)."""
    x = np.sort(np.random.default_rng(43).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    n = spl.bkpt_gpm.sum() - spl.nord
    assert lower.size == n - spl.nord + 1, 'lower index vector has wrong length'
    assert upper.size == n - spl.nord + 1, 'upper index vector has wrong length'


def test_build_design_matrix_full_data_coverage():
    """Every data point falls in at least one span of the design matrix."""
    x = np.sort(np.random.default_rng(44).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    covered = np.zeros(x.size, dtype=bool)
    for k in range(lower.size):
        if lower[k] <= upper[k]:
            covered[lower[k]:upper[k] + 1] = True
    assert covered.all(), 'every data point must fall in at least one span'


def test_build_design_matrix_consistent_with_find_spans():
    """lower/upper bounds from _build_design_matrix match _find_spans span indices."""
    x = np.sort(np.random.default_rng(45).uniform(0, 10, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    indx = BSpline._find_spans(x, spl.breakpoints[spl.bkpt_gpm], spl.nord)
    for k in range(lower.size):
        mask_k = (indx == k + spl.nord - 1)
        if mask_k.any():
            assert lower[k] == np.where(mask_k)[0][0], f'lower[{k}] must equal first x index in span'
            assert upper[k] == np.where(mask_k)[0][-1], f'upper[{k}] must equal last x index in span'


# ============================================================================
# BSpline._assemble_normal_equations
# ============================================================================

def test_assemble_normal_equations_alpha_shape():
    """_assemble_normal_equations() returns alpha of shape (bw, nfull+bw)."""
    rng = np.random.default_rng(13)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x) + 0.1 * rng.standard_normal(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, lower, upper)
    bw = A.shape[1]
    nn = spl.bkpt_gpm[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert alpha.shape == (bw, nfull + bw), f'expected alpha shape {(bw, nfull + bw)}'
    assert beta.shape == (nfull + bw,), f'expected beta shape {(nfull + bw,)}'


def test_assemble_normal_equations_diagonal_positive():
    """Diagonal of the normal-equation matrix is positive for well-posed data."""
    rng = np.random.default_rng(14)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x) + 0.1 * rng.standard_normal(80)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    alpha, beta = spl._assemble_normal_equations(A, y, lower, upper)
    nn = spl.bkpt_gpm[spl.nord:].sum()
    nfull = spl._poly_scale(nn)
    assert np.all(alpha[0, :nfull] > 0), 'normal equation diagonal must be positive for well-posed data'


def test_assemble_normal_equations_solution_matches_lstsq():
    """Banded Cholesky solution matches numpy lstsq for a clean 1D fit."""
    rng = np.random.default_rng(15)
    x = np.sort(rng.uniform(0, 5, 80))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=6), nord=4)

    # Banded solve via fit()
    err, yfit_banded = spl.fit(x, y)
    assert err == 0, 'banded fit must succeed for clean polynomial data'

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
# _cholesky_banded
# ============================================================================

def test_cholesky_banded_success_known_output():
    """A diagonal (bw=1) SPD matrix has an analytically known Cholesky
    factor: the element-wise square root of the diagonal."""
    ab = np.array([[4.0, 9.0, 16.0, 25.0]])
    chol, info = _cholesky_banded(ab, lower=True)
    assert info == 0, '_cholesky_banded must succeed on a positive definite matrix'
    np.testing.assert_allclose(chol[0], [2.0, 3.0, 4.0, 5.0])


def test_cholesky_banded_failure_identifies_leading_minor_index():
    """A tridiagonal matrix with an all-positive diagonal can still be
    indefinite due to the off-diagonal coupling; `info` must report the
    1-indexed order of the leading minor that first fails, not merely flag
    a negative diagonal entry (there is none here).

    A_dense = [[4, -1,  0,  0],
               [-1, 4, -1,  0],
               [ 0,-1,  4, -5],
               [ 0, 0, -5,  4]]
    has leading-minor determinants 4, 15, 56, -151: positive definiteness
    first fails at order 4.
    """
    ab = np.array([
        [4.0, 4.0, 4.0, 4.0],
        [-1.0, -1.0, -5.0, 0.0],
    ])
    assert np.all(ab[0] > 0), 'diagonal must be all-positive for this to be a meaningful test'
    chol, info = _cholesky_banded(ab, lower=True)
    assert info == 4, 'positive-definiteness must first fail at leading minor order 4'


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
    """_solve_banded() returns a solution and Cholesky factor for an SPD system."""
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert bad_cols[0] == -1, 'bad_cols sentinel must be -1 for a successful solve'
    assert sol is not None, 'solution must be returned for an SPD system'
    assert chol is not None, 'Cholesky factor must be returned for an SPD system'
    assert sol.shape == (n,), f'solution must have length {n}'


def test_solve_banded_solution_correct():
    """_solve_banded() solution matches numpy.linalg.solve for a tridiagonal system."""
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
    """Corrupted diagonal: _solve_banded returns sol=None and reports bad columns."""
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    alpha[0, 3] = -1.0  # corrupt diagonal
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert sol is None, 'corrupted diagonal must cause solve to fail'
    assert 3 in bad_cols, 'column 3 must be reported as bad after diagonal corruption'


def test_solve_banded_chol_diagonal_positive():
    """Cholesky diagonal is positive for a successful SPD solve."""
    alpha, beta, n = _make_tridiagonal_alpha_beta()
    spl = BSpline.__new__(BSpline)
    sol, chol, bad_cols = spl._solve_banded(alpha, beta, mininf=0.0)
    assert np.all(chol[0, :n] > 0), 'Cholesky diagonal must be positive for an SPD system'


def test_solve_banded_linalg_error_branch():
    """Indefinite matrix with positive diagonal: _solve_banded returns sol=None."""
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
    assert sol is None, 'indefinite matrix must cause solve to fail'
    assert bad_cols[0] != -1, 'indefinite matrix must report at least one bad column'


# ============================================================================
# BSpline._evaluate_model
# ============================================================================

def test_evaluate_model_matches_explicit_matmul():
    """_evaluate_model() result matches an explicit span-by-span matrix multiply."""
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
    """_mask_breakpoints() reduces the active-breakpoint count."""
    x = np.sort(np.random.default_rng(5).uniform(0, 10, 200))
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    n_before = spl.bkpt_gpm.sum()
    result = spl._mask_breakpoints(np.array([8]))
    assert result == -1, '_mask_breakpoints must return -1 (retry) after masking'
    assert spl.bkpt_gpm.sum() < n_before, 'active breakpoint count must decrease after masking'


def test_mask_breakpoints_invalidates_cache():
    """_mask_breakpoints() clears _cached_design."""
    x = np.sort(np.random.default_rng(6).uniform(0, 10, 200))
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl._cached_design = spl._build_design_matrix(x)
    spl._mask_breakpoints(np.array([8]))
    assert spl._cached_design is None, 'design-matrix cache must be cleared after masking breakpoints'


def test_mask_breakpoints_too_few_returns_minus2():
    """_mask_breakpoints() returns -2 when masking would leave fewer than 2*nord active breakpoints."""
    x = np.linspace(0, 5, 30)
    spl = BSpline(x=x, knots=Knots(count=3), nord=4)
    spl.bkpt_gpm[:] = False
    spl.bkpt_gpm[:2 * 4] = True  # leave only 2*nord active
    result = spl._mask_breakpoints(np.array([0]))
    assert result == -2, '_mask_breakpoints must return -2 when too few active breakpoints remain'


# ============================================================================
# BSpline.reset_coeff / BSpline2D.reset_coeff
# ============================================================================

def test_reset_coeff_1d_shape_and_zeros():
    """reset_coeff() sets coeff and icoeff to zero arrays of shape (nc,)."""
    x = np.linspace(0., 10., 300)
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    spl.fit(x, np.sin(x))
    nc = spl.breakpoints.size - spl.nord
    spl.reset_coeff()
    assert spl.coeff.shape == (nc,), f'Wrong coeff shape: {spl.coeff.shape}'
    assert spl.icoeff.shape == (nc,), f'Wrong icoeff shape: {spl.icoeff.shape}'
    assert np.all(spl.coeff == 0.), 'reset_coeff must zero coeff'
    assert np.all(spl.icoeff == 0.), 'reset_coeff must zero icoeff'


def test_reset_coeff_1d_leaves_mask_and_cache_unchanged():
    """reset_coeff() must not modify bkpt_gpm or invalidate _cached_design."""
    rng = np.random.default_rng(42)
    x = np.sort(rng.uniform(0., 10., 200))
    spl = BSpline(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, np.sin(x))
    gpm_before = spl.bkpt_gpm.copy()
    cache_id = id(spl._cached_design)
    spl.reset_coeff()
    np.testing.assert_array_equal(spl.bkpt_gpm, gpm_before)
    assert id(spl._cached_design) == cache_id, 'reset_coeff must not invalidate the design-matrix cache'


def test_reset_coeff_1d_before_breakpoints_raises():
    """reset_coeff() raises ValueError when breakpoints have not been established."""
    spl = BSpline.__new__(BSpline)
    spl.knots = Knots()
    spl.nord = 4
    with pytest.raises(ValueError):
        spl.reset_coeff()


def test_reset_coeff_2d_shape_and_zeros():
    """BSpline2D.reset_coeff() sets coeff and icoeff to zero arrays of shape (nc, npoly)."""
    rng = np.random.default_rng(43)
    x = np.sort(rng.uniform(0., 10., 300))
    basis_x = rng.uniform(0., 1., 300)
    npoly = 3
    spl = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, basis='legendre', npoly=npoly, xmin=0., xmax=1.)
    nc = spl.breakpoints.size - spl.nord
    spl.reset_coeff()
    assert spl.coeff.shape == (nc, npoly), f'Wrong coeff shape: {spl.coeff.shape}'
    assert spl.icoeff.shape == (nc, npoly), f'Wrong icoeff shape: {spl.icoeff.shape}'
    assert np.all(spl.coeff == 0.), 'reset_coeff must zero coeff'
    assert np.all(spl.icoeff == 0.), 'reset_coeff must zero icoeff'


def test_reset_coeff_2d_before_npoly_raises():
    """BSpline2D.reset_coeff() raises ValueError when npoly has not been set by fit()."""
    x = np.linspace(0., 5., 100)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    with pytest.raises(ValueError):
        spl.reset_coeff()


# ============================================================================
# BSpline.fit
# ============================================================================

def test_fit_cubic_polynomial_recovery():
    """fit() recovers a degree-3 polynomial exactly (error < 1e-8)."""
    rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, 300))
    y = 1.0 + 2.0*x - 0.5*x**2 + 0.1*x**3
    spl = BSpline(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == 0, 'fit must succeed on a cubic polynomial'
    np.testing.assert_allclose(yfit, y, atol=1e-8)


def test_fit_smooth_function_residuals():
    """fit() on sin(x) gives residual std < 2e-3 with 30 interior knots."""
    rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 2 * np.pi, 500))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == 0, 'fit must succeed on sin(x)'
    assert np.std(yfit - y) < 2e-3, 'residual std must be below 2e-3 for a smooth function'


def test_fit_returns_correct_length():
    """fit() returns yfit with the same length as x."""
    rng = np.random.default_rng(2)
    x = np.sort(rng.uniform(0, 5, 100))
    y = rng.standard_normal(100)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    err, yfit = spl.fit(x, y)
    assert yfit.shape == x.shape, 'yfit must have the same shape as x'


def test_fit_design_matrix_cached():
    """Design matrix is reused across fit() calls with the same x."""
    rng = np.random.default_rng(3)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, y)
    cache_id = id(spl._cached_design)
    spl.fit(x, y * 0.9)
    assert id(spl._cached_design) == cache_id, 'second fit() call must reuse the cached design matrix'


def test_fit_zero_invvar_points_ignored():
    """Points with ivar=0 are excluded and fit() still converges."""
    rng = np.random.default_rng(4)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    ivar = np.ones_like(x)
    ivar[::5] = 0.0
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    err, yfit = spl.fit(x, y, ivar=ivar)
    assert err == 0, 'fit must succeed when some ivar values are zero'
    assert np.isfinite(yfit).all(), 'yfit must be finite even with zero-ivar points'


def test_fit_reset_knots_rebuilds_breakpoints():
    """reset_knots=True rebuilds breakpoints spanning the new x range."""
    rng = np.random.default_rng(50)
    x1 = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1))
    bkpt_before = spl.breakpoints.copy()
    x2 = np.sort(rng.uniform(0, 10, 100))
    spl.fit(x2, np.sin(x2), reset_knots=True)
    assert spl.breakpoints.min() <= x2.min(), 'new breakpoints must span x2.min()'
    assert spl.breakpoints.max() >= x2.max(), 'new breakpoints must span x2.max()'
    assert not np.array_equal(spl.breakpoints, bkpt_before), 'breakpoints must change after reset_knots on a new x range'


def test_fit_reset_knots_resets_mask():
    """reset_knots=True restores a fully-True bkpt_gpm."""
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y)
    spl.bkpt_gpm[5] = False
    assert not spl.bkpt_gpm.all(), 'manually masked breakpoint must not be all-True'
    spl.fit(x, y, reset_knots=True)
    assert spl.bkpt_gpm.all(), 'reset_knots=True must restore a fully-True bkpt_gpm'


def test_fit_returns_minus1_on_degenerate_cholesky():
    """Data confined to a subset of the knot span returns err=-1."""
    rng = np.random.default_rng(60)
    # Data only in [0, 2]; knots span [0, 10] -> spans in [2, 10] are empty,
    # leaving zero diagonal entries that trigger the degenerate pre-check.
    x = np.sort(rng.uniform(0, 2, 100))
    y = np.sin(x)
    full_bkpt = np.linspace(0, 10, 40)
    spl = BSpline(knots=Knots(full=full_bkpt), nord=4)
    err, yfit = spl.fit(x, y)
    assert err == -1, 'sparse data must trigger degenerate Cholesky (err=-1)'
    assert yfit.shape == x.shape, 'yfit must have correct shape even on degenerate fit'


def test_fit_ivar_none_matches_unit_weights():
    """ivar=None and ivar=ones give identical fit results."""
    rng = np.random.default_rng(70)
    x = np.sort(rng.uniform(0, 5, 150))
    y = np.sin(x)
    spl_w = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl_n = BSpline(x=x, knots=Knots(count=10), nord=4)
    err_w, yfit_w = spl_w.fit(x, y, ivar=np.ones(150))
    err_n, yfit_n = spl_n.fit(x, y)
    assert err_w == 0 and err_n == 0, 'both unit-ivar and no-ivar fits must succeed'
    np.testing.assert_array_equal(yfit_w, yfit_n)


def test_fit_nonuniform_ivar_affects_result():
    """Non-uniform ivar changes the fitted values relative to unit weights."""
    rng = np.random.default_rng(72)
    N = 300
    x = np.sort(rng.uniform(0.5, 5.0, N))
    ivar = 1.0 / x**2  # noise amplitude proportional to x
    y = np.sin(x) + rng.standard_normal(N) / np.sqrt(ivar)
    spl_u = BSpline(x=x, knots=Knots(count=15), nord=4)
    spl_n = BSpline(x=x, knots=Knots(count=15), nord=4)
    err_u, yfit_u = spl_u.fit(x, y)
    err_n, yfit_n = spl_n.fit(x, y, ivar=ivar)
    assert err_u == 0 and err_n == 0, 'both weighted and unweighted fits must succeed'
    assert np.isfinite(yfit_n).all(), 'weighted fit must produce finite values'
    assert not np.allclose(yfit_u, yfit_n), 'weighted and unweighted fits must differ under non-uniform noise'


# ============================================================================
# BSpline.value
# ============================================================================

def test_value_matches_fit_at_training_points():
    """value() at training x matches the yfit returned by fit()."""
    rng = np.random.default_rng(10)
    x = np.sort(rng.uniform(0, 10, 200))
    y = 2.0 + x - 0.3 * x**2
    spl = BSpline(x=x, knots=Knots(count=12), nord=4)
    err, yfit_fit = spl.fit(x, y)
    yfit_val, _ = spl.value(x)
    np.testing.assert_allclose(yfit_fit, yfit_val, atol=1e-12)


def test_value_out_of_range_masked():
    """value() masks points outside the breakpoint range."""
    x = np.linspace(1, 9, 100)
    y = np.sin(x)
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y)
    x_eval = np.array([0.0, 5.0, 10.0])
    yfit, mask = spl.value(x_eval)
    assert not mask[0], 'x below breakpoint range must be masked'
    assert mask[1], 'x inside breakpoint range must not be masked'
    assert not mask[2], 'x above breakpoint range must be masked'


def test_value_unsorted_input():
    """value() on unsorted x gives results consistent with sorted evaluation."""
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
    """value() masks x points that fall in a manually removed knot-span gap."""
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
    assert not vmask[0], 'x in a masked knot-span gap must be masked by value()'


# ============================================================================
# BSpline.copy
# ============================================================================

def test_bspline_copy_attributes():
    """copy() preserves nord, breakpoints, mask, and coefficients."""
    rng = np.random.default_rng(62)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    cp = spl.copy()
    assert cp.nord == spl.nord, 'copy must preserve nord'
    np.testing.assert_array_equal(cp.breakpoints, spl.breakpoints)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.icoeff, spl.icoeff)


def test_bspline_copy_arrays_are_independent():
    """copy() deep-copies arrays; mutating the source does not affect the copy."""
    rng = np.random.default_rng(63)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    cp = spl.copy()
    orig_val = cp.coeff[0]
    spl.coeff[0] += 999.0
    assert cp.coeff[0] == orig_val, 'mutating source coeff must not affect the copy'


def test_bspline_copy_cache_is_cleared():
    """copy() starts with a cold design-matrix cache."""
    rng = np.random.default_rng(64)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSpline(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, np.sin(x))
    assert spl._cached_design is not None, 'design-matrix cache must be populated after fit()'
    cp = spl.copy()
    assert cp._cached_design is None, 'copy must start with a cold design-matrix cache'


# ============================================================================
# BSpline2D._normalize_basis_x
# ============================================================================

def test_normalize_basis_x_maps_xmin_to_minus_one():
    """_normalize_basis_x() maps xmin to -1."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([2.0])), [-1.0])


def test_normalize_basis_x_maps_xmax_to_plus_one():
    """_normalize_basis_x() maps xmax to +1."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 2.0
    spl.xmax = 8.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([8.0])), [1.0])


def test_normalize_basis_x_midpoint_is_zero():
    """_normalize_basis_x() maps the midpoint to 0."""
    spl = BSpline2D.__new__(BSpline2D)
    spl.xmin = 0.0
    spl.xmax = 4.0
    np.testing.assert_allclose(spl._normalize_basis_x(np.array([2.0])), [0.0])


def test_normalize_basis_x_linear_mapping():
    """_normalize_basis_x() maps [xmin, xmax] linearly to [-1, 1]."""
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
    assert spl.xmin == 0.0, 'explicit xmin kwarg must not mutate self.xmin'
    assert spl.xmax == 10.0, 'explicit xmax kwarg must not mutate self.xmax'


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
    assert spl.xmin is None, 'stored xmin must remain None when all bounds were None'
    assert spl.xmax is None, 'stored xmax must remain None when all bounds were None'


# ============================================================================
# BSpline2D._poly_basis
# ============================================================================

def _make_poly_spl(funcname, npoly):
    spl = BSpline2D.__new__(BSpline2D)
    spl.npoly = npoly
    spl.funcname = funcname
    return spl


def test_poly_basis_legendre_shape():
    """_poly_basis() with 'legendre' returns shape (N, npoly)."""
    spl = _make_poly_spl('legendre', 4)
    P = spl._poly_basis(np.linspace(-1, 1, 20))
    assert P.shape == (20, 4), 'Legendre basis must have shape (N, npoly)'


def test_poly_basis_legendre_first_column_is_one():
    """First Legendre basis column is identically 1."""
    spl = _make_poly_spl('legendre', 3)
    P = spl._poly_basis(np.linspace(-1, 1, 10))
    np.testing.assert_allclose(P[:, 0], 1.0, atol=1e-14)


def test_poly_basis_chebyshev_shape():
    """_poly_basis() with 'chebyshev' returns shape (N, npoly)."""
    spl = _make_poly_spl('chebyshev', 3)
    P = spl._poly_basis(np.linspace(-1, 1, 15))
    assert P.shape == (15, 3), 'Chebyshev basis must have shape (N, npoly)'


def test_poly_basis_chebyshev_first_is_one_second_is_x():
    """First Chebyshev basis column is 1; second is x."""
    spl = _make_poly_spl('chebyshev', 2)
    basis_xnorm = np.linspace(-1, 1, 10)
    P = spl._poly_basis(basis_xnorm)
    np.testing.assert_allclose(P[:, 0], 1.0, atol=1e-14)
    np.testing.assert_allclose(P[:, 1], basis_xnorm, atol=1e-14)


def test_poly_basis_poly_monomial_structure():
    """_poly_basis() with 'poly' gives columns x^0, x^1, ..., x^(npoly-1)."""
    spl = _make_poly_spl('poly', 4)
    basis_xnorm = np.array([0.0, 0.5, 1.0])
    P = spl._poly_basis(basis_xnorm)
    for k in range(4):
        np.testing.assert_allclose(P[:, k], basis_xnorm**k, atol=1e-14)


def test_poly_basis_unknown_funcname_raises():
    """_poly_basis() raises ValueError for an unrecognised function name."""
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
    assert P.shape == (10, 5), 'npoly kwarg must override number of basis columns'
    assert spl.npoly == 3, 'npoly kwarg must not mutate self.npoly'


def test_poly_basis_funcname_kwarg_does_not_mutate_stored():
    """Passing funcname kwarg must not change self.funcname."""
    spl = _make_poly_spl('legendre', 2)
    spl._poly_basis(np.linspace(-1, 1, 5), funcname='chebyshev')
    assert spl.funcname == 'legendre', 'funcname kwarg must not mutate self.funcname'


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
    """_build_design_matrix() returns shape (N, nord*npoly) for BSpline2D."""
    rng = np.random.default_rng(20)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = _setup_bspline2d(x, basis_x, npoly=3, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.shape == (x.size, spl.nord * spl.npoly), f'2D design matrix must have shape (N, nord*npoly)'


def test_bspline2d_build_design_matrix_outer_product_structure():
    """Each column of the 2D design matrix equals B_i * P_j for the corresponding B-spline and polynomial basis."""
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
    """fit() raises ValueError when the basis array row count does not match x.size."""
    rng = np.random.default_rng(22)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    with pytest.raises(ValueError):
        spl.fit(x, np.sin(x), basis_x=basis_x, basis=np.ones((basis_x.size - 5, 3)))


def test_bspline2d_build_design_matrix_c_order():
    """BSpline2D._build_design_matrix() returns a C-contiguous array."""
    rng = np.random.default_rng(23)
    x = np.sort(rng.uniform(0, 10, 100))
    basis_x = rng.uniform(0, 1, 100)
    spl = _setup_bspline2d(x, basis_x, npoly=3, xmin=0.0, xmax=1.0,
                           funcname='legendre', knots=Knots(count=8), nord=4)
    A, lower, upper = spl._build_design_matrix(x)
    assert A.flags['C_CONTIGUOUS'], '2D design matrix must be C-contiguous'


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
    assert not np.allclose(A_stored, A_alt), 'different P kwarg must change the design matrix'


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
    assert A3.shape == (N, spl.nord * 3), 'npoly kwarg must change the design matrix column count'


# ============================================================================
# BSpline2D._evaluate_model
# ============================================================================

def test_bspline2d_evaluate_model_matches_explicit_einsum():
    """_evaluate_model() result matches an explicit einsum span-by-span contraction."""
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
    assert not np.allclose(yfit_stored, yfit_alt), 'different coeff kwarg must change the evaluated model'


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
    assert err == 0, 'Legendre 2D fit must succeed for a linear polynomial signal'
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_smooth_function():
    """BSpline2D.fit() on a sinusoidal modulated signal gives residual std < 2e-2."""
    rng = np.random.default_rng(31)
    N = 500
    x = np.sort(rng.uniform(0, 2 * np.pi, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    spl = BSpline2D(x=x, knots=Knots(count=30), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=3, xmin=0.0, xmax=1.0)
    assert err == 0, '2D fit on a sinusoidal signal must succeed'
    assert np.std(yfit - y) < 2e-2, 'residual std must be below 2e-2 for the sinusoidal signal'


def test_bspline2d_fit_exact_polynomial_chebyshev():
    """BSpline2D.fit() with Chebyshev basis recovers a linear polynomial exactly."""
    rng = np.random.default_rng(32)
    N = 400
    x = np.sort(rng.uniform(0, 8, N))
    basis_x = rng.uniform(-1, 1, N)
    y = 2.0 + basis_x  # T0 + T1 in Chebyshev basis
    spl = BSpline2D(x=x, knots=Knots(count=20), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='chebyshev', npoly=2, xmin=-1.0, xmax=1.0)
    assert err == 0, 'Chebyshev 2D fit must succeed for a linear polynomial signal'
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_exact_polynomial_poly():
    """BSpline2D.fit() with 'poly' basis recovers a linear polynomial exactly."""
    rng = np.random.default_rng(33)
    N = 400
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 2, N)
    basis_xnorm = 2.0 * (basis_x - 1.0) / 2.0 - 0.0   # xmin=0, xmax=2 → norm range [-1, 1]
    y = 1.0 + basis_xnorm
    spl = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, basis='poly', npoly=2, xmin=0.0, xmax=2.0)
    assert err == 0, "'poly' basis 2D fit must succeed for a linear polynomial signal"
    np.testing.assert_allclose(yfit, y, atol=1e-6)


def test_bspline2d_fit_reset_knots_rebuilds_breakpoints():
    """reset_knots=True rebuilds 2D-fit breakpoints spanning the new x range."""
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
    assert spl.breakpoints.min() <= x3.min(), 'new breakpoints must span x3.min()'
    assert spl.breakpoints.max() >= x3.max(), 'new breakpoints must span x3.max()'
    assert not np.array_equal(spl.breakpoints, bkpt_before), 'breakpoints must change after 2D reset_knots on a new x range'


def test_bspline2d_fit_reset_knots_resets_mask():
    """reset_knots=True restores a fully-True bkpt_gpm for a 2D fit."""
    rng = np.random.default_rng(53)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x)
    spl = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    spl.bkpt_gpm[5] = False
    assert not spl.bkpt_gpm.all(), 'manually masked breakpoint must not be all-True'
    spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0, reset_knots=True)
    assert spl.bkpt_gpm.all(), '2D reset_knots=True must restore a fully-True bkpt_gpm'


def test_bspline2d_fit_returns_minus1_on_degenerate_cholesky():
    """2D fit with data confined to a subset of the knot span returns err=-1."""
    rng = np.random.default_rng(66)
    x = np.sort(rng.uniform(0, 2, 150))
    basis_x = rng.uniform(0, 1, 150)
    y = np.sin(x)
    full_bkpt = np.linspace(0, 10, 40)
    spl = BSpline2D(knots=Knots(full=full_bkpt), nord=4)
    err, yfit = spl.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    assert err == -1, '2D fit with sparse data must trigger degenerate Cholesky (err=-1)'
    assert yfit.shape == x.shape, '2D yfit must have correct shape even on degenerate fit'


def test_bspline2d_fit_ivar_none_matches_unit_weights():
    """ivar=None and ivar=ones give identical 2D fit results."""
    rng = np.random.default_rng(71)
    N = 200
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    y = np.sin(x) + 0.1 * basis_x
    spl_w = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    spl_n = BSpline2D(x=x, knots=Knots(count=10), nord=4)
    err_w, yfit_w = spl_w.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0, ivar=np.ones(N))
    err_n, yfit_n = spl_n.fit(x, y, basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    assert err_w == 0 and err_n == 0, 'both unit-ivar and no-ivar 2D fits must succeed'
    np.testing.assert_array_equal(yfit_w, yfit_n)


def test_bspline2d_fit_nonuniform_ivar_affects_result():
    """Non-uniform ivar changes the 2D fitted values relative to unit weights."""
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
    assert err_u == 0 and err_n == 0, 'both weighted and unweighted 2D fits must succeed'
    assert np.isfinite(yfit_n).all(), 'weighted 2D fit must produce finite values'
    assert not np.allclose(yfit_u, yfit_n), 'weighted and unweighted 2D fits must differ under non-uniform noise'


# ============================================================================
# BSpline2D.value
# ============================================================================

def test_bspline2d_value_matches_fit_at_training_points():
    """BSpline2D.value() at training x matches the yfit returned by fit()."""
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
    """BSpline2D.value() raises ValueError when basis_x is omitted for a string-basis fit."""
    rng = np.random.default_rng(41)
    x = np.sort(rng.uniform(0, 5, 50))
    basis_x = rng.uniform(0, 1, 50)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, npoly=2, xmin=0.0, xmax=1.0)
    with pytest.raises(ValueError):
        spl.value(x)  # basis_x=None not permitted when funcname is set


def test_bspline2d_value_gap_masking():
    """BSpline2D.value() masks x points in a manually removed knot-span gap."""
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
    assert not vmask[0], '2D value() must mask x in a manually removed knot-span gap'


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
    assert not np.allclose(yfit_stored, yfit_alt), 'different xmin/xmax bounds must change the evaluated model'


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

    assert spl.npoly == npoly_before, 'value() override kwargs must not mutate npoly'
    assert spl.funcname == funcname_before, 'value() override kwargs must not mutate funcname'
    assert spl.xmin == xmin_before, 'value() override kwargs must not mutate xmin'
    assert spl.xmax == xmax_before, 'value() override kwargs must not mutate xmax'
    np.testing.assert_array_equal(spl.coeff, coeff_before)
    np.testing.assert_array_equal(spl.P, P_before)


# ============================================================================
# BSpline2D.copy
# ============================================================================

def test_bspline2d_copy_attributes():
    """BSpline2D.copy() preserves npoly, xmin, xmax, funcname, nord, and breakpoints."""
    rng = np.random.default_rng(68)
    N = 100
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, basis='chebyshev', npoly=3, xmin=0.0, xmax=1.0)
    cp = spl.copy()
    assert cp.npoly == spl.npoly, 'copy must preserve npoly'
    assert cp.xmin == spl.xmin, 'copy must preserve xmin'
    assert cp.xmax == spl.xmax, 'copy must preserve xmax'
    assert cp.funcname == spl.funcname, 'copy must preserve funcname'
    assert cp.nord == spl.nord, 'copy must preserve nord'
    assert cp._cached_design is None, 'copy must start with a cold design-matrix cache'
    assert cp.basis_x is spl.basis_x, 'copy must share the basis_x reference'
    np.testing.assert_array_equal(cp.breakpoints, spl.breakpoints)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bspline2d_copy_coeff_independent():
    """BSpline2D.copy() deep-copies coeff; mutating the original does not affect the copy."""
    rng = np.random.default_rng(69)
    N = 100
    x = np.sort(rng.uniform(0, 5, N))
    basis_x = rng.uniform(0, 1, N)
    spl = BSpline2D(x=x, knots=Knots(count=8), nord=4)
    spl.fit(x, np.sin(x), basis_x=basis_x, npoly=3, xmin=0.0, xmax=1.0)
    cp = spl.copy()
    orig_val = cp.coeff[0, 0]
    spl.coeff[0, 0] += 999.0
    assert cp.coeff[0, 0] == orig_val, 'mutating source coeff must not affect the 2D copy'


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

    assert err_s == 0 and err_a == 0, 'both string-basis and array-basis 2D fits must succeed'
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
    assert max_diff > 1e-4, 'interpolated values must differ from full evaluation (linear error must be nonzero)'
    assert max_diff < 5e-3, 'interpolated values must closely approximate full evaluation (max diff < 5e-3)'


def test_value_interpolate_training_x_fast_path():
    """When x is self.x the identity fast path fires before the interpolate check."""
    rng = np.random.default_rng(93)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.cos(x)
    spl = BSpline(x=x, knots=Knots(count=8), nord=4)
    err, yfit_fit = spl.fit(x, y)

    yfit_val, _ = spl.value(x, interpolate=True)
    assert yfit_val is yfit_fit, 'interpolate=True at training x must return the cached yfit'


# ---------------------------------------------------------------------------
# iterative_bspline_fit — integration tests against reference data
# ---------------------------------------------------------------------------

def test_iterative_bspline_fit_spec():
    """
    iterative_bspline_fit reproduces the spectral flat-field reference fit
    from the gemini_gnirs_32 test data files.
    """
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_spec_fit.npz'.format(slit))
             for slit in [0, 1]]
    logrej = 0.5
    spec_samp_fine = 1.2
    for f in files:
        d = np.load(f)
        sset, _, spec_flat_fit, _, exit_status = iterative_bspline_fit(
            d['spec_coo_data'], d['spec_flat_data'],
            ivar=d['spec_ivar_data'], gpm=d['spec_gpm_data'],
            nord=4, upper=logrej, lower=logrej,
            kwargs_knots={'spacing': spec_samp_fine},
            kwargs_reject={'groupbadpix': True, 'maxrej': 5},
        )
        assert isinstance(sset, BSpline), 'spectral fit must return a BSpline'
        assert np.allclose(d['spec_flat_fit'], spec_flat_fit), \
            'Bad spectral iterative_bspline_fit result'


def test_iterative_bspline_fit_spat():
    """
    iterative_bspline_fit reproduces the spatial flat-field reference fit
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
        sset, _, spat_flat_fit, _, exit_status = iterative_bspline_fit(
            d['spat_coo_data'], d['spat_flat_data'],
            nord=4, upper=5.0, lower=5.0,
            kwargs_knots={'spacing': bkspace},
        )
        assert isinstance(sset, BSpline), 'spatial fit must return a BSpline'
        assert np.allclose(d['spat_flat_fit'], spat_flat_fit), \
            'Bad spatial iterative_bspline_fit result'


def test_iterative_bspline_fit_twod():
    """
    iterative_bspline_fit reproduces the 2D flat-field reference fit
    from the gemini_gnirs_32 test data files.
    """
    files = [dataPaths.tests.get_file_path('gemini_gnirs_32_{0}_twod_fit.npz'.format(slit))
             for slit in [0, 1]]
    spec_samp_coarse = 50.0
    twod_sigrej = 4.0
    for f in files:
        d = np.load(f)
        sset, _, twod_flat_fit, _, exit_status = iterative_bspline_fit(
            d['twod_spec_coo_data'], d['twod_flat_data'],
            basis=d['poly_basis'],
            ivar=d['twod_ivar_data'], gpm=d['twod_gpm_data'],
            nord=4, upper=twod_sigrej, lower=twod_sigrej,
            kwargs_knots={'spacing': spec_samp_coarse},
            kwargs_reject={'groupbadpix': True, 'maxrej': 10},
        )
        assert isinstance(sset, BSpline2D), '2D fit must return a BSpline2D'
        assert np.allclose(d['twod_flat_fit'], twod_flat_fit), \
            'Bad 2D iterative_bspline_fit result'


# ============================================================================
# iterative_bspline_fit option coverage and agreement tests
# ============================================================================

def test_iterative_bspline_fit_string_basis():
    """basis='legendre' + basis_x activates the 2D BSpline2D path."""
    rng = np.random.default_rng(42)
    n = 300
    x = np.sort(rng.uniform(0, 10, n))
    basis_x = rng.uniform(-1, 1, n)
    npoly = 3
    y_true = (1.0 + 0.3 * basis_x) * (np.sin(x / 2) + 2.0)
    y = y_true + rng.normal(0, 0.05, n)
    ivar = np.full(n, 400.0)

    sset, outmask, yfit, reduced_chi, exit_status = iterative_bspline_fit(
        x, y, ivar=ivar, nord=4, npoly=npoly,
        basis='legendre', basis_x=basis_x, xmin=-1.0, xmax=1.0,
        upper=5.0, lower=5.0,
        kwargs_knots={'spacing': 1.0},
    )
    assert isinstance(sset, BSpline2D), 'string-basis fit must return a BSpline2D'
    assert exit_status == 0, 'string-basis fit must converge normally'
    assert np.allclose(y_true, yfit, atol=0.15), 'fitted values must approximate the true signal within 0.15'


def test_iterative_bspline_fit_basis_array_shapes():
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

    _, _, yfit_2d, _, _ = iterative_bspline_fit(
        x, y, ivar=ivar, basis=basis_2d, kwargs_knots={'spacing': 0.5},
    )
    _, _, yfit_1d, _, _ = iterative_bspline_fit(
        x, y, ivar=ivar, basis=basis_1d, kwargs_knots={'spacing': 0.5},
    )
    assert np.allclose(yfit_2d, yfit_1d), '2D and auto-reshaped 1D basis arrays must give identical fit'


def test_iterative_bspline_fit_exit_statuses():
    """Exit status 4 (too few points) and 1 (maxiter exceeded)."""
    rng = np.random.default_rng(0)
    n = 300
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x) + rng.normal(0, 0.1, n)
    ivar = np.full(n, 100.0)

    # exit_status = 4: only 3 good points, fewer than nord=4
    gpm_few = np.zeros(n, dtype=bool)
    gpm_few[:3] = True
    _, _, _, _, es = iterative_bspline_fit(
        x, y, ivar=ivar, gpm=gpm_few, nord=4,
        kwargs_knots={'spacing': 1.0},
    )
    assert es == 4, 'too few good points must give exit status 4'

    # exit_status = 1: maxiter=1 with outliers that need exactly one rejection
    # pass (converges in 2 iterations, but maxiter caps at 1)
    y_out = y.copy()
    y_out[rng.choice(n, 15, replace=False)] += 20.0
    _, _, _, _, es = iterative_bspline_fit(
        x, y_out, ivar=ivar, nord=4, maxiter=1, upper=3.0, lower=3.0,
        kwargs_knots={'spacing': 1.0},
    )
    assert es == 1, 'maxiter=1 with outliers must give exit status 1'


# ============================================================================
# bspline_qa
# ============================================================================

def _make_bspline_fit(rng=None, n=300, spacing=1.0):
    """Return (sset, gpm, xdata, ydata, yfit) from iterative_bspline_fit."""
    if rng is None:
        rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x) + rng.normal(0, 0.1, n)
    ivar = np.full(n, 100.0)
    sset, gpm, yfit, _, _ = iterative_bspline_fit(
        x, y, ivar=ivar, nord=4, upper=5, lower=5,
        kwargs_knots={'spacing': spacing},
    )
    return sset, gpm, x, y, yfit


def test_bspline_qa_returns_axes_when_show_false():
    """bspline_qa returns a matplotlib Axes when show=False."""
    sset, gpm, xdata, ydata, yfit = _make_bspline_fit()
    matplotlib.pyplot.figure()
    ax = bspline_qa(xdata, ydata, sset, gpm, yfit, show=False)
    assert isinstance(ax, matplotlib.axes.Axes), 'show=False must return a matplotlib Axes'
    matplotlib.pyplot.close('all')


def test_bspline_qa_breakpoints_on_axes():
    """bspline_qa plots breakpoints at the correct x positions on the axes."""
    sset, gpm, xdata, ydata, yfit = _make_bspline_fit()
    matplotlib.pyplot.figure()
    ax = bspline_qa(xdata, ydata, sset, gpm, yfit, show=False)

    # The breakpoint line is the 4th line added (index 3).
    bkpt_line = ax.lines[3]
    plotted_bkx = bkpt_line.get_xdata()
    expected_bkx = sset.breakpoints[sset.bkpt_gpm]
    assert np.allclose(plotted_bkx, expected_bkx), 'breakpoint x positions on axes must match bkpt_gpm-selected breakpoints'
    matplotlib.pyplot.close('all')


def test_bspline_qa_breakpoint_y_from_value():
    """bspline_qa evaluates breakpoint y-values via BSpline.value."""
    sset, gpm, xdata, ydata, yfit = _make_bspline_fit()
    matplotlib.pyplot.figure()
    ax = bspline_qa(xdata, ydata, sset, gpm, yfit, show=False)

    bkpt_line = ax.lines[3]
    plotted_bky = bkpt_line.get_ydata()
    expected_bky, _ = sset.value(sset.breakpoints[sset.bkpt_gpm])
    assert np.allclose(plotted_bky, expected_bky), 'breakpoint y positions on axes must equal BSpline.value at breakpoints'
    matplotlib.pyplot.close('all')
