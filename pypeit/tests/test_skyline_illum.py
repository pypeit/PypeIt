"""Tests for sky-line-based fiber illumination correction."""
import numpy as np

from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph


spec = MMTBINOSPECIFUSpectrograph()


def _make_mock_data(nfibers=40, nspec=2000, wmin=5000.0, wmax=9000.0,
                    sky_fiber_boost=1.3, n_sky_fibers=4):
    """Create synthetic 2D data with sky lines and fiber throughput
    variation.

    Uses ~2 Ang/pixel dispersion (realistic for Binospec G270) so
    that the 4-Ang line measurement window spans multiple pixels.

    Returns sciimg, waveimg, slitmask, spat_ids, sky_indices.
    """
    nspat_per_fiber = 3
    nspat = nfibers * nspat_per_fiber
    sciimg = np.zeros((nspec, nspat))
    waveimg = np.zeros((nspec, nspat))
    slitmask = np.zeros((nspec, nspat), dtype=int)

    wave_arr = np.linspace(wmin, wmax, nspec)
    spat_ids = np.arange(1, nfibers + 1)
    sky_indices = np.arange(n_sky_fibers)

    # Throughput varies smoothly across fibers, sky fibers are brighter
    throughput = np.ones(nfibers)
    throughput[sky_indices] = sky_fiber_boost
    # Add some variation to science fibers
    rng = np.random.default_rng(42)
    throughput[n_sky_fibers:] += 0.05 * rng.standard_normal(
        nfibers - n_sky_fibers)

    # Sky continuum + emission lines
    sky_lines_ang = np.array([5577.34, 6300.304, 6863.951, 7340.881,
                              7993.327, 8465.353])
    line_sigma_ang = 2.0  # Angstroms (FWHM ~ 4.7 Ang)

    for i in range(nfibers):
        cols = slice(i * nspat_per_fiber, (i + 1) * nspat_per_fiber)
        slitmask[:, cols] = spat_ids[i]
        waveimg[:, cols] = wave_arr[:, None]

        # Sky spectrum: continuum + lines (in wavelength space)
        sky_spec = np.full(nspec, 100.0)  # continuum
        for wl in sky_lines_ang:
            sky_spec += 500.0 * np.exp(
                -0.5 * ((wave_arr - wl) / line_sigma_ang) ** 2)

        # Apply throughput
        sky_spec *= throughput[i]
        sciimg[:, cols] = sky_spec[:, None]

    return sciimg, waveimg, slitmask, spat_ids, sky_indices


def test_uniform_throughput_gives_unit_correction():
    """If all fibers have same throughput, correction should be ~1."""
    # Use n_sky_fibers=0 so no fibers get a boost, and the random
    # variation still applies to all fibers equally
    sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
        sky_fiber_boost=1.0, n_sky_fibers=0)
    corr = spec.compute_skyline_illum(
        sciimg, waveimg, slitmask, spat_ids)
    # Correction should track the random throughput variation (~5%),
    # so individual fibers may deviate from 1.0, but the median
    # correction across all fibers should be ~1.0
    assert np.abs(np.median(corr[slitmask > 0]) - 1.0) < 0.02, \
        "uniform throughput should give a median correction of ~1.0, got " \
        f"{np.median(corr[slitmask > 0]):.4f}"


def test_sky_fibers_get_higher_correction():
    """Sky fibers with 30% more throughput should get corr > 1."""
    sciimg, waveimg, slitmask, spat_ids, sky_idx = _make_mock_data(
        sky_fiber_boost=1.3)
    corr = spec.compute_skyline_illum(
        sciimg, waveimg, slitmask, spat_ids)
    # Sky fiber pixels should have correction > 1
    for i in sky_idx:
        sky_pix = slitmask == spat_ids[i]
        assert np.median(corr[sky_pix]) > 1.1, \
            f"sky fiber {spat_ids[i]} (30% throughput boost) should get a " \
            f"correction > 1.1, got {np.median(corr[sky_pix]):.3f}"


def test_correction_normalizes_sky_level():
    """After applying correction, sky line flux should be
    equalized across fibers."""
    sciimg, waveimg, slitmask, spat_ids, sky_idx = _make_mock_data(
        sky_fiber_boost=1.3)
    corr = spec.compute_skyline_illum(
        sciimg, waveimg, slitmask, spat_ids)
    corrected = sciimg / np.where(corr > 0.1, corr, 1.0)

    # Measure flux in a sky line region for sky vs science fibers
    wave_arr = waveimg[:, 1]
    line_pix = np.abs(wave_arr - 5577.34) < 10
    sky_flux = []
    sci_flux = []
    nspat_per = 3
    for i in range(len(spat_ids)):
        cols = slice(i * nspat_per, (i + 1) * nspat_per)
        f = np.sum(corrected[line_pix, cols])
        if i in sky_idx:
            sky_flux.append(f)
        else:
            sci_flux.append(f)
    # After correction, sky and science fibers should have similar
    # flux (within 5%)
    ratio = np.mean(sky_flux) / np.mean(sci_flux)
    assert abs(ratio - 1.0) < 0.05, \
        "after correction, sky/science sky-line flux ratio should be ~1.0 " \
        f"(within 5%), got {ratio:.4f}"


def test_no_lines_in_range():
    """If no sky lines fall in wavelength range, return ones."""
    sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
        wmin=3000.0, wmax=4000.0)
    corr = spec.compute_skyline_illum(
        sciimg, waveimg, slitmask, spat_ids)
    assert np.allclose(corr, 1.0), \
        "no sky lines in the wavelength range should yield a unit " \
        "(all-ones) correction"


def test_extreme_corrections_clipped():
    """Corrections should be clipped to [0.3, 3.0]."""
    sciimg, waveimg, slitmask, spat_ids, _ = _make_mock_data(
        sky_fiber_boost=5.0)
    corr = spec.compute_skyline_illum(
        sciimg, waveimg, slitmask, spat_ids)
    assert np.all(corr >= 0.3), \
        f"corrections must be clipped to >= 0.3, min was {corr.min():.3f}"
    assert np.all(corr <= 3.0), \
        f"corrections must be clipped to <= 3.0, max was {corr.max():.3f}"
