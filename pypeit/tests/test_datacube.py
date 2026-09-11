"""Unit tests for pypeit.core.datacube helpers."""
import numpy as np

from pypeit import utils
from pypeit.core.datacube import resample_spec_to_grid


def test_resample_spec_to_grid_identity():
    """Resampling onto the native grid reproduces flux/ivar on interior pixels.

    The flux-conserving resampler maps each interior output pixel one-to-one
    onto its native pixel, so the values are reproduced exactly there.  The two
    edge pixels can be flagged as partially covered (output borders float just
    outside the native coverage), so they are excluded from the comparison.
    """
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.sin(wave / 100.0) + 5.0
    ivar = np.full_like(wave, 4.0)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, wave)
    interior = slice(1, -1)
    assert np.all(cov[interior]), \
        "interior output pixels should all be flagged as covered"
    np.testing.assert_allclose(
        fg[interior], flux[interior], rtol=1e-6,
        err_msg="resampled flux should match the input on interior pixels")
    np.testing.assert_allclose(
        ig[interior], ivar[interior], rtol=1e-6,
        err_msg="resampled ivar should match the input on interior pixels")


def test_resample_spec_to_grid_propagates_variance_in_quadrature():
    """Error is propagated in quadrature (variance space), not ivar space.

    Resampling onto the native grid leaves the per-pixel variance unchanged on
    interior pixels, so the returned inverse variance matches the input there.
    A naive average of the inverse variance would not preserve this.
    """
    wave = np.linspace(4000.0, 4100.0, 51)
    flux = np.ones_like(wave)
    # Inverse variance varying across the spectrum.
    ivar = np.linspace(1.0, 9.0, wave.size)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, wave)
    interior = slice(1, -1)
    np.testing.assert_allclose(
        ig[interior], ivar[interior], rtol=1e-6,
        err_msg="variance must be propagated in quadrature; interior ivar "
                "should match the input")


def test_resample_spec_to_grid_zero_outside_coverage():
    """Grid points outside the native wavelength range stay zero."""
    wave = np.linspace(4500.0, 4600.0, 51)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave)
    grid = np.linspace(4000.0, 5000.0, 101)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, grid)
    outside = (grid < 4500.0) | (grid > 4600.0)
    assert np.all(fg[outside] == 0.0), \
        "flux outside the native coverage should be zero"
    assert np.all(ig[outside] == 0.0), \
        "ivar outside the native coverage should be zero"
    assert not np.any(cov[outside]), \
        "no output pixel outside the native coverage should be flagged covered"
    # Interior of the covered range carries real flux.
    assert np.any(cov), \
        "at least some output pixels within the native range should be covered"


def test_resample_spec_to_grid_ivar_none():
    """A None inverse variance yields all-zero ivar but valid flux."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    grid = np.linspace(4100.0, 4900.0, 51)
    fg, ig, cov = resample_spec_to_grid(wave, flux, None, grid)
    assert np.all(cov), \
        "all in-range output pixels should be covered when resampling flat flux"
    np.testing.assert_allclose(
        fg[cov], 1.0, rtol=1e-6,
        err_msg="resampled flat (unit) flux should remain 1.0 where covered")
    assert np.all(ig == 0.0), \
        "ivar should be all zero when no input inverse variance is supplied"


def test_resample_spec_to_grid_below_min_good():
    """Too few valid samples returns all-zero arrays and an empty mask."""
    wave = np.array([4000.0, -1.0, 0.0])  # only one positive wavelength
    flux = np.array([1.0, 1.0, 1.0])
    ivar = np.array([1.0, 1.0, 1.0])
    grid = np.linspace(4000.0, 5000.0, 11)
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, grid, min_good=2)
    assert not np.any(cov), \
        "too few valid samples should leave every output pixel uncovered"
    assert np.all(fg == 0.0), \
        "too few valid samples should give all-zero flux"
    assert np.all(ig == 0.0), \
        "too few valid samples should give all-zero ivar"


def test_resample_spec_to_grid_masked_native_pixels_drop_out():
    """Native pixels with ivar <= 0 are masked and carry no inverse variance."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    ivar = np.zeros_like(wave)  # no weight anywhere
    fg, ig, cov = resample_spec_to_grid(wave, flux, ivar, wave)
    assert np.all(ig == 0.0), \
        "masked native pixels (ivar <= 0) should carry no inverse variance"
    # utils.inverse semantics: zero ivar -> zero variance contribution
    assert np.all(utils.inverse(ig) == 0.0), \
        "utils.inverse of zero ivar should be zero"
