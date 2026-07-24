"""
Unit tests for :class:`~pypeit.containers.bspline.BSplineContainer` and
:class:`~pypeit.containers.bspline.BSpline2DContainer`.

All random-number generators use a fixed seed for reproducibility.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from pypeit.core.bspline import BSpline, BSpline2D, Knots
from pypeit.containers.bspline import BSplineContainer, BSpline2DContainer


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_1d_fitted(rng=None, n=300):
    """Return a fitted BSplineContainer and the (x, y) data used to fit it."""
    if rng is None:
        rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x) + 0.05 * rng.standard_normal(n)
    spl = BSplineContainer(x=x, knots=Knots(count=20), nord=4)
    spl.fit(x, y)
    return spl, x, y


def _make_2d_fitted_string(rng=None, n=300):
    """Return a fitted BSpline2DContainer using a string (Legendre) basis."""
    if rng is None:
        rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 8, n))
    basis_x = rng.uniform(0, 1, n)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    spl = BSpline2DContainer(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    return spl, x, basis_x, y


def _make_2d_fitted_array(rng=None, n=300):
    """Return a fitted BSpline2DContainer using a pre-built array basis."""
    if rng is None:
        rng = np.random.default_rng(2)
    x = np.sort(rng.uniform(0, 8, n))
    basis_x = rng.uniform(0, 1, n)
    P = np.column_stack([np.ones(n), basis_x])
    y = np.sin(x) * (1 + 0.2 * basis_x)
    spl = BSpline2DContainer(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis=P)
    return spl, x, basis_x, P, y


# ===========================================================================
# BSplineContainer — construction and basic behaviour
# ===========================================================================

def test_bsplinecontainer_is_bspline():
    """BSplineContainer is also a BSpline instance."""
    rng = np.random.default_rng(10)
    x = np.sort(rng.uniform(0, 10, 200))
    spl = BSplineContainer(x=x, knots=Knots(count=15), nord=4)
    assert isinstance(spl, BSpline), 'BSplineContainer must be a BSpline instance'
    assert isinstance(spl, BSplineContainer), 'BSplineContainer must be a BSplineContainer instance'


def test_bsplinecontainer_fit_value_match_bspline():
    """fit() and value() results are numerically identical to BSpline."""
    rng = np.random.default_rng(11)
    x = np.sort(rng.uniform(0, 10, 300))
    y = np.sin(x)

    ref = BSpline(x=x, knots=Knots(count=20), nord=4)
    err_r, yfit_r = ref.fit(x, y)

    cnt = BSplineContainer(x=x, knots=Knots(count=20), nord=4)
    err_c, yfit_c = cnt.fit(x, y)

    assert err_r == 0 and err_c == 0, 'both fits must succeed'
    np.testing.assert_array_equal(yfit_r, yfit_c)


def test_bsplinecontainer_bkpt_full_synced_after_init():
    """bkpt_full is populated from knots.breakpoints at construction."""
    rng = np.random.default_rng(12)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSplineContainer(x=x, knots=Knots(count=10), nord=4)
    np.testing.assert_array_equal(spl.bkpt_full, spl.knots.breakpoints)


def test_bsplinecontainer_bkpt_full_synced_after_reset_knots():
    """bkpt_full stays in sync after fit() with reset_knots=True."""
    rng = np.random.default_rng(13)
    x1 = np.sort(rng.uniform(0, 5, 100))
    x2 = np.sort(rng.uniform(0, 10, 100))
    spl = BSplineContainer(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1))
    spl.fit(x2, np.sin(x2), reset_knots=True)
    np.testing.assert_array_equal(spl.bkpt_full, spl.knots.breakpoints)


# ===========================================================================
# BSplineContainer — from_bspline classmethod
# ===========================================================================

def test_bsplinecontainer_from_bspline_attributes():
    """from_bspline copies nord, breakpoints, mask, and coefficients."""
    rng = np.random.default_rng(20)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    ref = BSpline(x=x, knots=Knots(count=10), nord=4)
    ref.fit(x, y)

    cnt = BSplineContainer.from_bspline(ref)
    assert isinstance(cnt, BSplineContainer), 'from_bspline must return a BSplineContainer'
    assert cnt.nord == ref.nord, 'from_bspline must preserve nord'
    np.testing.assert_array_equal(cnt.bkpt_full, ref.knots.breakpoints)
    np.testing.assert_array_equal(cnt.bkpt_gpm, ref.bkpt_gpm)
    np.testing.assert_array_equal(cnt.coeff, ref.coeff)
    np.testing.assert_array_equal(cnt.icoeff, ref.icoeff)


def test_bsplinecontainer_from_bspline_arrays_independent():
    """from_bspline deep-copies arrays; mutating the source does not affect the container."""
    rng = np.random.default_rng(21)
    x = np.sort(rng.uniform(0, 5, 100))
    ref = BSpline(x=x, knots=Knots(count=10), nord=4)
    ref.fit(x, np.sin(x))

    cnt = BSplineContainer.from_bspline(ref)
    orig_coeff = cnt.coeff[0]
    ref.coeff[0] += 999.0
    assert cnt.coeff[0] == orig_coeff, 'mutating the source must not affect the container coeff'


def test_bsplinecontainer_from_bspline_value_matches():
    """Container created by from_bspline evaluates identically to the source BSpline."""
    rng = np.random.default_rng(22)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    ref = BSpline(x=x, knots=Knots(count=10), nord=4)
    ref.fit(x, y)

    cnt = BSplineContainer.from_bspline(ref)
    x_eval = np.sort(rng.uniform(0.5, 4.5, 50))
    yfit_ref, _ = ref.value(x_eval)
    yfit_cnt, _ = cnt.value(x_eval)
    np.testing.assert_allclose(yfit_cnt, yfit_ref, atol=1e-14)


# ===========================================================================
# BSplineContainer — copy
# ===========================================================================

def test_bsplinecontainer_copy_returns_container():
    """copy() returns a BSplineContainer."""
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert isinstance(cp, BSplineContainer), 'copy must return a BSplineContainer'


def test_bsplinecontainer_copy_attributes_match():
    """copy() preserves nord, coeff, and bkpt_gpm."""
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert cp.nord == spl.nord, 'copy must preserve nord'
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bsplinecontainer_copy_arrays_independent():
    """copy() deep-copies arrays; mutating the original does not affect the copy."""
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    orig = cp.coeff[0]
    spl.coeff[0] += 999.0
    assert cp.coeff[0] == orig, 'mutating the source must not affect the copy coeff'


def test_bsplinecontainer_copy_yfit_copied():
    """copy() carries over yfit as an independent copy."""
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert cp.yfit is not None, 'copy must carry over yfit'
    np.testing.assert_array_equal(cp.yfit, spl.yfit)
    spl.yfit[0] += 999.0
    assert cp.yfit[0] != spl.yfit[0], 'mutating the original yfit must not affect the copy'


def test_bsplinecontainer_copy_cache_cleared():
    """copy() starts with a cold design-matrix cache."""
    spl, x, y = _make_1d_fitted()
    assert spl._cached_design is not None, 'design-matrix cache must be populated after fit()'
    cp = spl.copy()
    assert cp._cached_design is None, 'copy must start with a cold design-matrix cache'


# ===========================================================================
# BSplineContainer — round-trip FITS serialisation
# ===========================================================================

def test_bsplinecontainer_roundtrip_value():
    """Round-trip FITS I/O preserves the spline evaluation."""
    spl, x, y = _make_1d_fitted()
    x_eval = np.sort(np.random.default_rng(30).uniform(0.5, 9.5, 80))
    yfit_orig, _ = spl.value(x_eval)

    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    yfit_loaded, _ = loaded.value(x_eval)
    np.testing.assert_allclose(yfit_loaded, yfit_orig, atol=1e-14)


def test_bsplinecontainer_roundtrip_attributes():
    """Round-trip FITS I/O preserves nord, breakpoints, mask, and coefficients."""
    spl, x, y = _make_1d_fitted()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    assert loaded.nord == spl.nord, 'round-trip must preserve nord'
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)
    np.testing.assert_array_equal(loaded.bkpt_gpm, spl.bkpt_gpm)
    np.testing.assert_allclose(loaded.coeff, spl.coeff, atol=1e-14)
    np.testing.assert_allclose(loaded.icoeff, spl.icoeff, atol=1e-14)


def test_bsplinecontainer_roundtrip_knots_rebuilt():
    """_validate() reconstructs the Knots instance from bkpt_full on load."""
    spl, x, y = _make_1d_fitted()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    assert loaded.knots is not None, 'loaded container must have a reconstructed Knots instance'
    np.testing.assert_array_equal(loaded.knots.breakpoints, spl.knots.breakpoints)


def test_bsplinecontainer_roundtrip_after_reset_knots():
    """Round-trip is correct after fitting a second data range with reset_knots."""
    rng = np.random.default_rng(31)
    x1 = np.sort(rng.uniform(0, 5, 100))
    spl = BSplineContainer(x=x1, knots=Knots(count=10), nord=4)
    spl.fit(x1, np.sin(x1))
    x2 = np.sort(rng.uniform(0, 10, 100))
    spl.fit(x2, np.sin(x2), reset_knots=True)

    x_eval = np.sort(rng.uniform(0.5, 9.5, 60))
    yfit_orig, _ = spl.value(x_eval)

    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    yfit_loaded, _ = loaded.value(x_eval)
    np.testing.assert_allclose(yfit_loaded, yfit_orig, atol=1e-14)


def test_bsplinecontainer_roundtrip_overwrite():
    """to_file with overwrite=True must not raise when the file already exists."""
    spl, x, y = _make_1d_fitted()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        spl.to_file(str(path), overwrite=True)
        loaded = BSplineContainer.from_file(str(path))
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)


def test_bsplinecontainer_empty_write():
    """An unfitted BSplineContainer must serialise and load without error,
    with all array attributes remaining None."""
    spl = BSplineContainer()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline_empty.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))
    assert loaded.bkpt_full is None, 'empty container bkpt_full must be None after round-trip'
    assert loaded.coeff is None, 'empty container coeff must be None after round-trip'


# ===========================================================================
# BSpline2DContainer — construction and basic behaviour
# ===========================================================================

def test_bspline2dcontainer_is_bspline2d():
    """BSpline2DContainer is also a BSpline2D instance."""
    rng = np.random.default_rng(40)
    x = np.sort(rng.uniform(0, 8, 200))
    spl = BSpline2DContainer(x=x, knots=Knots(count=12), nord=4)
    assert isinstance(spl, BSpline2D), 'BSpline2DContainer must be a BSpline2D instance'
    assert isinstance(spl, BSpline2DContainer), 'BSpline2DContainer must be a BSpline2DContainer instance'


def test_bspline2dcontainer_fit_value_match_bspline2d():
    """fit() and value() results are numerically identical to BSpline2D."""
    rng = np.random.default_rng(41)
    x = np.sort(rng.uniform(0, 8, 300))
    basis_x = rng.uniform(0, 1, 300)
    y = np.sin(x) * (1 + 0.2 * basis_x)

    ref = BSpline2D(x=x, knots=Knots(count=15), nord=4)
    err_r, yfit_r = ref.fit(x, y, basis_x=basis_x, basis='legendre',
                            npoly=2, xmin=0.0, xmax=1.0)

    cnt = BSpline2DContainer(x=x, knots=Knots(count=15), nord=4)
    err_c, yfit_c = cnt.fit(x, y, basis_x=basis_x, basis='legendre',
                             npoly=2, xmin=0.0, xmax=1.0)

    assert err_r == 0 and err_c == 0, 'both 2D fits must succeed'
    np.testing.assert_array_equal(yfit_r, yfit_c)


# ===========================================================================
# BSpline2DContainer — from_bspline2d classmethod
# ===========================================================================

def test_bspline2dcontainer_from_bspline2d_string_basis():
    """from_bspline2d copies string-basis metadata; basis array is not stored."""
    rng = np.random.default_rng(50)
    x = np.sort(rng.uniform(0, 8, 200))
    basis_x = rng.uniform(0, 1, 200)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    ref = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    ref.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)

    cnt = BSpline2DContainer.from_bspline2d(ref)
    assert isinstance(cnt, BSpline2DContainer), 'from_bspline2d must return a BSpline2DContainer'
    assert cnt.funcname == 'legendre', "from_bspline2d must preserve funcname as 'legendre'"
    assert cnt.npoly == 2, 'from_bspline2d must preserve npoly'
    assert cnt.xmin == 0.0, 'from_bspline2d must preserve xmin'
    assert cnt.xmax == 1.0, 'from_bspline2d must preserve xmax'
    assert cnt.basis is None, 'string-basis fit must not store the basis array'


def test_bspline2dcontainer_from_bspline2d_array_basis():
    """from_bspline2d stores the basis array when the fit used an array basis."""
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 8, 200))
    basis_x = rng.uniform(0, 1, 200)
    P = np.column_stack([np.ones(200), basis_x])
    y = np.sin(x) * (1 + 0.2 * basis_x)
    ref = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    ref.fit(x, y, basis_x=basis_x, basis=P)

    cnt = BSpline2DContainer.from_bspline2d(ref)
    assert cnt.funcname is None, 'array-basis fit must have funcname=None'
    assert cnt.basis is not None, 'array-basis fit must store the basis array'
    np.testing.assert_array_equal(cnt.basis, ref.P)


def test_bspline2dcontainer_from_bspline2d_value_matches():
    """Container created by from_bspline2d evaluates identically to the source BSpline2D."""
    rng = np.random.default_rng(52)
    x = np.sort(rng.uniform(0, 8, 200))
    basis_x = rng.uniform(0, 1, 200)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    ref = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    ref.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)

    cnt = BSpline2DContainer.from_bspline2d(ref)
    x_eval = np.sort(rng.uniform(0.5, 7.5, 80))
    basis_x_eval = rng.uniform(0, 1, 80)
    yfit_ref, _ = ref.value(x_eval, basis_x=basis_x_eval)
    yfit_cnt, _ = cnt.value(x_eval, basis_x=basis_x_eval)
    np.testing.assert_allclose(yfit_cnt, yfit_ref, atol=1e-14)


# ===========================================================================
# BSpline2DContainer — copy
# ===========================================================================

def test_bspline2dcontainer_copy_returns_container():
    """copy() returns a BSpline2DContainer."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    assert isinstance(cp, BSpline2DContainer), 'copy must return a BSpline2DContainer'


def test_bspline2dcontainer_copy_string_basis_attributes():
    """copy() preserves funcname, npoly, coeff, and bkpt_gpm for a string-basis fit."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    assert cp.funcname == spl.funcname, 'copy must preserve funcname'
    assert cp.npoly == spl.npoly, 'copy must preserve npoly'
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bspline2dcontainer_copy_string_basis_P_copied():
    """copy() deep-copies P; mutating the original does not affect the copy."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    assert spl.P is not None, 'P must be populated after a string-basis fit'
    cp = spl.copy()
    assert cp.P is not None, 'copy must carry over P'
    np.testing.assert_array_equal(cp.P, spl.P)
    orig = cp.P[0, 0]
    spl.P[0, 0] += 999.0
    assert cp.P[0, 0] == orig, 'mutating the source P must not affect the copy'


def test_bspline2dcontainer_copy_array_basis_P_restored():
    """copy() restores P from the stored basis array for array-basis fits."""
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    cp = spl.copy()
    assert cp.P is not None, 'copy of an array-basis fit must carry over P'
    np.testing.assert_array_equal(cp.P, spl.P)


def test_bspline2dcontainer_copy_arrays_independent():
    """copy() deep-copies coeff; mutating the original does not affect the copy."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    orig = cp.coeff[0, 0]
    spl.coeff[0, 0] += 999.0
    assert cp.coeff[0, 0] == orig, 'mutating source coeff must not affect the 2D copy'


# ===========================================================================
# BSpline2DContainer — round-trip FITS serialisation (string basis)
# ===========================================================================

def test_bspline2dcontainer_roundtrip_string_basis_value():
    """Round-trip FITS I/O preserves string-basis evaluation."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    rng = np.random.default_rng(60)
    x_eval = np.sort(rng.uniform(0.5, 7.5, 80))
    basis_x_eval = rng.uniform(0, 1, 80)
    yfit_orig, _ = spl.value(x_eval, basis_x=basis_x_eval)

    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    yfit_loaded, _ = loaded.value(x_eval, basis_x=basis_x_eval)
    np.testing.assert_allclose(yfit_loaded, yfit_orig, atol=1e-14)


def test_bspline2dcontainer_roundtrip_string_basis_attributes():
    """Round-trip FITS I/O preserves string-basis metadata and coefficients."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.nord == spl.nord, 'round-trip must preserve nord'
    assert loaded.npoly == spl.npoly, 'round-trip must preserve npoly'
    assert loaded.funcname == spl.funcname, 'round-trip must preserve funcname'
    assert abs(loaded.xmin - spl.xmin) < 1e-14, 'round-trip must preserve xmin'
    assert abs(loaded.xmax - spl.xmax) < 1e-14, 'round-trip must preserve xmax'
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)
    np.testing.assert_allclose(loaded.coeff, spl.coeff, atol=1e-14)


def test_bspline2dcontainer_roundtrip_string_basis_has_no_stored_basis():
    """String-basis fit: basis array is None after round-trip."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.basis is None, 'string-basis round-trip must not restore basis array'


# ===========================================================================
# BSpline2DContainer — round-trip FITS serialisation (array basis)
# ===========================================================================

def test_bspline2dcontainer_roundtrip_array_basis_basis_stored():
    """Array-basis fit: basis array is stored and survives round-trip."""
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_arr.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.basis is not None, 'array-basis round-trip must store the basis array'
    np.testing.assert_allclose(loaded.basis, spl.P, atol=1e-14)


def test_bspline2dcontainer_roundtrip_array_basis_P_restored():
    """Array-basis fit: P is restored from the stored basis on load."""
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_arr.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.P is not None, 'array-basis round-trip must restore P'
    np.testing.assert_allclose(loaded.P, spl.P, atol=1e-14)


def test_bspline2dcontainer_roundtrip_array_basis_value():
    """Round-trip FITS I/O preserves array-basis evaluation."""
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    rng = np.random.default_rng(70)
    x_eval = np.sort(rng.uniform(0.5, 7.5, 80))
    basis_x_eval = rng.uniform(0, 1, 80)
    P_eval = np.column_stack([np.ones(80), basis_x_eval])
    yfit_orig, _ = spl.value(x_eval, basis=P_eval)

    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_arr.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    yfit_loaded, _ = loaded.value(x_eval, basis=P_eval)
    np.testing.assert_allclose(yfit_loaded, yfit_orig, atol=1e-14)


def test_bspline2dcontainer_roundtrip_overwrite():
    """to_file with overwrite=True must not raise when the file already exists."""
    spl, x, basis_x, y = _make_2d_fitted_string()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        spl.to_file(str(path), overwrite=True)
        loaded = BSpline2DContainer.from_file(str(path))
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)


def test_bspline2dcontainer_empty_write():
    """An unfitted BSpline2DContainer must serialise and load without error,
    with all array attributes remaining None."""
    spl = BSpline2DContainer()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_empty.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))
    assert loaded.bkpt_full is None, 'empty 2D container bkpt_full must be None after round-trip'
    assert loaded.coeff is None, 'empty 2D container coeff must be None after round-trip'


# ===========================================================================
# BSpline2DContainer — version and datamodel
# ===========================================================================

def test_bsplinecontainer_version():
    """BSplineContainer version string is '2.0.0'."""
    assert BSplineContainer.version == '2.0.0', "BSplineContainer version must be '2.0.0'"


def test_bspline2dcontainer_version():
    """BSpline2DContainer version string is '2.0.0'."""
    assert BSpline2DContainer.version == '2.0.0', "BSpline2DContainer version must be '2.0.0'"


def test_bspline2dcontainer_datamodel_has_basis():
    """BSpline2DContainer datamodel includes a 'basis' field."""
    assert 'basis' in BSpline2DContainer.datamodel, "BSpline2DContainer datamodel must include a 'basis' field"
