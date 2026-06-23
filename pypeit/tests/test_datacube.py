"""Unit tests for pypeit.core.datacube helpers."""
import numpy as np

from pypeit import utils
from pypeit.core.datacube import resample_spec_to_grid


def test_resample_spec_to_grid_identity():
    """Resampling onto the native grid returns flux/ivar unchanged inside it."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.sin(wave / 100.0) + 5.0
    ivar = np.full_like(wave, 4.0)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, wave)
    assert np.all(cov)
    np.testing.assert_allclose(fg, flux)
    np.testing.assert_allclose(ig, ivar)


def test_resample_spec_to_grid_interpolates_variance_not_ivar():
    """Variance (1/ivar) must be interpolated, not the inverse variance.

    With ivar varying linearly between two native pixels, interpolating ivar
    directly and interpolating variance then inverting give different answers
    at the midpoint.  The helper must use the variance-space result.
    """
    wave = np.array([4000.0, 4002.0])
    flux = np.array([1.0, 1.0])
    ivar = np.array([1.0, 9.0])  # var = 1.0 and 1/9
    mid = np.array([4001.0])
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, mid)
    # Variance-space: var = (1.0 + 1/9)/2 = 5/9 -> ivar = 9/5 = 1.8
    np.testing.assert_allclose(ig, 1.8)
    # The naive (wrong) ivar-space interpolation would give (1+9)/2 = 5.0
    assert not np.isclose(ig[0], 5.0)


def test_resample_spec_to_grid_zero_outside_coverage():
    """Grid points outside the native wavelength range stay zero."""
    wave = np.linspace(4500.0, 4600.0, 51)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave)
    grid = np.linspace(4000.0, 5000.0, 101)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, grid)
    outside = (grid < 4500.0) | (grid > 4600.0)
    assert np.all(cov == ~outside)
    assert np.all(fg[outside] == 0.0)
    assert np.all(ig[outside] == 0.0)


def test_resample_spec_to_grid_ivar_none():
    """A None inverse variance yields all-zero ivar but valid flux."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    grid = np.linspace(4100.0, 4900.0, 51)
    fg, ig, cov = resample_spec_to_grid(wave, flux, None, grid)
    assert np.all(cov)
    np.testing.assert_allclose(fg, 1.0)
    assert np.all(ig == 0.0)


def test_resample_spec_to_grid_below_min_good():
    """Too few valid samples returns all-zero arrays and an empty mask."""
    wave = np.array([4000.0, -1.0, 0.0])  # only one positive wavelength
    flux = np.array([1.0, 1.0, 1.0])
    ivar = np.array([1.0, 1.0, 1.0])
    grid = np.linspace(4000.0, 5000.0, 11)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, grid, min_good=2)
    assert not np.any(cov)
    assert np.all(fg == 0.0)
    assert np.all(ig == 0.0)


def test_resample_spec_to_grid_nonpositive_ivar_drops_out():
    """Native pixels with ivar <= 0 contribute zero inverse variance."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    ivar = np.zeros_like(wave)  # no weight anywhere
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, wave)
    assert np.all(ig == 0.0)
    # utils.inverse semantics: zero ivar -> zero variance contribution
    assert np.all(utils.inverse(ig) == 0.0)
