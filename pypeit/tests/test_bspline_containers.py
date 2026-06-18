"""
Unit tests for :class:`~pypeit.bspline.containers.BSplineContainer` and
:class:`~pypeit.bspline.containers.BSpline2DContainer`.

All random-number generators use a fixed seed for reproducibility.
"""

import tempfile
from pathlib import Path

import numpy as np
import pytest

from pypeit.core.bspline.refactor import BSpline, BSpline2D, Knots
from pypeit.core.bspline.containers import BSplineContainer, BSpline2DContainer


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_1d_fitted(rng=None, n=300):
    if rng is None:
        rng = np.random.default_rng(0)
    x = np.sort(rng.uniform(0, 10, n))
    y = np.sin(x) + 0.05 * rng.standard_normal(n)
    spl = BSplineContainer(x=x, knots=Knots(count=20), nord=4)
    spl.fit(x, y)
    return spl, x, y


def _make_2d_fitted_string(rng=None, n=300):
    if rng is None:
        rng = np.random.default_rng(1)
    x = np.sort(rng.uniform(0, 8, n))
    basis_x = rng.uniform(0, 1, n)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    spl = BSpline2DContainer(x=x, knots=Knots(count=15), nord=4)
    spl.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)
    return spl, x, basis_x, y


def _make_2d_fitted_array(rng=None, n=300):
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
    rng = np.random.default_rng(10)
    x = np.sort(rng.uniform(0, 10, 200))
    spl = BSplineContainer(x=x, knots=Knots(count=15), nord=4)
    assert isinstance(spl, BSpline)
    assert isinstance(spl, BSplineContainer)


def test_bsplinecontainer_fit_value_match_bspline():
    rng = np.random.default_rng(11)
    x = np.sort(rng.uniform(0, 10, 300))
    y = np.sin(x)

    ref = BSpline(x=x, knots=Knots(count=20), nord=4)
    err_r, yfit_r = ref.fit(x, y)

    cnt = BSplineContainer(x=x, knots=Knots(count=20), nord=4)
    err_c, yfit_c = cnt.fit(x, y)

    assert err_r == 0 and err_c == 0
    np.testing.assert_array_equal(yfit_r, yfit_c)


def test_bsplinecontainer_bkpt_full_synced_after_init():
    rng = np.random.default_rng(12)
    x = np.sort(rng.uniform(0, 5, 100))
    spl = BSplineContainer(x=x, knots=Knots(count=10), nord=4)
    np.testing.assert_array_equal(spl.bkpt_full, spl.knots.breakpoints)


def test_bsplinecontainer_bkpt_full_synced_after_reset_knots():
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
    rng = np.random.default_rng(20)
    x = np.sort(rng.uniform(0, 5, 100))
    y = np.sin(x)
    ref = BSpline(x=x, knots=Knots(count=10), nord=4)
    ref.fit(x, y)

    cnt = BSplineContainer.from_bspline(ref)
    assert isinstance(cnt, BSplineContainer)
    assert cnt.nord == ref.nord
    np.testing.assert_array_equal(cnt.bkpt_full, ref.knots.breakpoints)
    np.testing.assert_array_equal(cnt.bkpt_gpm, ref.bkpt_gpm)
    np.testing.assert_array_equal(cnt.coeff, ref.coeff)
    np.testing.assert_array_equal(cnt.icoeff, ref.icoeff)


def test_bsplinecontainer_from_bspline_arrays_independent():
    rng = np.random.default_rng(21)
    x = np.sort(rng.uniform(0, 5, 100))
    ref = BSpline(x=x, knots=Knots(count=10), nord=4)
    ref.fit(x, np.sin(x))

    cnt = BSplineContainer.from_bspline(ref)
    orig_coeff = cnt.coeff[0]
    ref.coeff[0] += 999.0
    assert cnt.coeff[0] == orig_coeff


def test_bsplinecontainer_from_bspline_value_matches():
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
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert isinstance(cp, BSplineContainer)


def test_bsplinecontainer_copy_attributes_match():
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert cp.nord == spl.nord
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bsplinecontainer_copy_arrays_independent():
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    orig = cp.coeff[0]
    spl.coeff[0] += 999.0
    assert cp.coeff[0] == orig


def test_bsplinecontainer_copy_yfit_copied():
    spl, x, y = _make_1d_fitted()
    cp = spl.copy()
    assert cp.yfit is not None
    np.testing.assert_array_equal(cp.yfit, spl.yfit)
    spl.yfit[0] += 999.0
    assert cp.yfit[0] != spl.yfit[0]


def test_bsplinecontainer_copy_cache_cleared():
    spl, x, y = _make_1d_fitted()
    assert spl._cached_design is not None
    cp = spl.copy()
    assert cp._cached_design is None


# ===========================================================================
# BSplineContainer — round-trip FITS serialisation
# ===========================================================================

def test_bsplinecontainer_roundtrip_value():
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
    spl, x, y = _make_1d_fitted()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    assert loaded.nord == spl.nord
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)
    np.testing.assert_array_equal(loaded.bkpt_gpm, spl.bkpt_gpm)
    np.testing.assert_allclose(loaded.coeff, spl.coeff, atol=1e-14)
    np.testing.assert_allclose(loaded.icoeff, spl.icoeff, atol=1e-14)


def test_bsplinecontainer_roundtrip_knots_rebuilt():
    spl, x, y = _make_1d_fitted()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline.fits'
        spl.to_file(str(path))
        loaded = BSplineContainer.from_file(str(path))

    assert loaded.knots is not None
    np.testing.assert_array_equal(loaded.knots.breakpoints, spl.knots.breakpoints)


def test_bsplinecontainer_roundtrip_after_reset_knots():
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
    assert loaded.bkpt_full is None
    assert loaded.coeff is None


# ===========================================================================
# BSpline2DContainer — construction and basic behaviour
# ===========================================================================

def test_bspline2dcontainer_is_bspline2d():
    rng = np.random.default_rng(40)
    x = np.sort(rng.uniform(0, 8, 200))
    spl = BSpline2DContainer(x=x, knots=Knots(count=12), nord=4)
    assert isinstance(spl, BSpline2D)
    assert isinstance(spl, BSpline2DContainer)


def test_bspline2dcontainer_fit_value_match_bspline2d():
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

    assert err_r == 0 and err_c == 0
    np.testing.assert_array_equal(yfit_r, yfit_c)


# ===========================================================================
# BSpline2DContainer — from_bspline2d classmethod
# ===========================================================================

def test_bspline2dcontainer_from_bspline2d_string_basis():
    rng = np.random.default_rng(50)
    x = np.sort(rng.uniform(0, 8, 200))
    basis_x = rng.uniform(0, 1, 200)
    y = np.sin(x) * (1 + 0.2 * basis_x)
    ref = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    ref.fit(x, y, basis_x=basis_x, basis='legendre', npoly=2, xmin=0.0, xmax=1.0)

    cnt = BSpline2DContainer.from_bspline2d(ref)
    assert isinstance(cnt, BSpline2DContainer)
    assert cnt.funcname == 'legendre'
    assert cnt.npoly == 2
    assert cnt.xmin == 0.0
    assert cnt.xmax == 1.0
    assert cnt.basis is None  # string basis: basis not stored


def test_bspline2dcontainer_from_bspline2d_array_basis():
    rng = np.random.default_rng(51)
    x = np.sort(rng.uniform(0, 8, 200))
    basis_x = rng.uniform(0, 1, 200)
    P = np.column_stack([np.ones(200), basis_x])
    y = np.sin(x) * (1 + 0.2 * basis_x)
    ref = BSpline2D(x=x, knots=Knots(count=12), nord=4)
    ref.fit(x, y, basis_x=basis_x, basis=P)

    cnt = BSpline2DContainer.from_bspline2d(ref)
    assert cnt.funcname is None
    assert cnt.basis is not None
    np.testing.assert_array_equal(cnt.basis, ref.P)


def test_bspline2dcontainer_from_bspline2d_value_matches():
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
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    assert isinstance(cp, BSpline2DContainer)


def test_bspline2dcontainer_copy_string_basis_attributes():
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    assert cp.funcname == spl.funcname
    assert cp.npoly == spl.npoly
    np.testing.assert_array_equal(cp.coeff, spl.coeff)
    np.testing.assert_array_equal(cp.bkpt_gpm, spl.bkpt_gpm)


def test_bspline2dcontainer_copy_string_basis_P_copied():
    spl, x, basis_x, y = _make_2d_fitted_string()
    assert spl.P is not None
    cp = spl.copy()
    assert cp.P is not None
    np.testing.assert_array_equal(cp.P, spl.P)
    orig = cp.P[0, 0]
    spl.P[0, 0] += 999.0
    assert cp.P[0, 0] == orig


def test_bspline2dcontainer_copy_array_basis_P_restored():
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    cp = spl.copy()
    assert cp.P is not None
    np.testing.assert_array_equal(cp.P, spl.P)


def test_bspline2dcontainer_copy_arrays_independent():
    spl, x, basis_x, y = _make_2d_fitted_string()
    cp = spl.copy()
    orig = cp.coeff[0, 0]
    spl.coeff[0, 0] += 999.0
    assert cp.coeff[0, 0] == orig


# ===========================================================================
# BSpline2DContainer — round-trip FITS serialisation (string basis)
# ===========================================================================

def test_bspline2dcontainer_roundtrip_string_basis_value():
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
    spl, x, basis_x, y = _make_2d_fitted_string()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.nord == spl.nord
    assert loaded.npoly == spl.npoly
    assert loaded.funcname == spl.funcname
    assert abs(loaded.xmin - spl.xmin) < 1e-14
    assert abs(loaded.xmax - spl.xmax) < 1e-14
    np.testing.assert_array_equal(loaded.bkpt_full, spl.bkpt_full)
    np.testing.assert_allclose(loaded.coeff, spl.coeff, atol=1e-14)


def test_bspline2dcontainer_roundtrip_string_basis_has_no_stored_basis():
    spl, x, basis_x, y = _make_2d_fitted_string()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.basis is None


# ===========================================================================
# BSpline2DContainer — round-trip FITS serialisation (array basis)
# ===========================================================================

def test_bspline2dcontainer_roundtrip_array_basis_basis_stored():
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_arr.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.basis is not None
    np.testing.assert_allclose(loaded.basis, spl.P, atol=1e-14)


def test_bspline2dcontainer_roundtrip_array_basis_P_restored():
    spl, x, basis_x, P, y = _make_2d_fitted_array()
    with tempfile.TemporaryDirectory() as td:
        path = Path(td) / 'bspline2d_arr.fits'
        spl.to_file(str(path))
        loaded = BSpline2DContainer.from_file(str(path))

    assert loaded.P is not None
    np.testing.assert_allclose(loaded.P, spl.P, atol=1e-14)


def test_bspline2dcontainer_roundtrip_array_basis_value():
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
    assert loaded.bkpt_full is None
    assert loaded.coeff is None


# ===========================================================================
# BSpline2DContainer — version and datamodel
# ===========================================================================

def test_bsplinecontainer_version():
    assert BSplineContainer.version == '2.0.0'


def test_bspline2dcontainer_version():
    assert BSpline2DContainer.version == '2.0.0'


def test_bspline2dcontainer_datamodel_has_basis():
    assert 'basis' in BSpline2DContainer.datamodel
