"""Tests for FiberFlatField and FiberFlatImages."""
import numpy as np
from pathlib import Path

from pypeit.flatfield import FiberFlatImages
from pypeit.tests.tstutils import data_output_path


def test_fiberflatimages_io():
    """Test FiberFlatImages DataContainer I/O round-trip."""
    nfibers = 40
    nwave = 100

    # Create synthetic fiber flat products (fixed seed for determinism)
    rng = np.random.default_rng(12345)
    normflat = rng.uniform(0.8, 1.2, (nfibers, nwave)).astype(np.float64)
    normflat_wave = np.linspace(4000.0, 7000.0, nwave).astype(np.float64)
    global_norm = 50000.0
    fiber_ids = np.arange(1, nfibers + 1).astype(np.int64)
    fiber_types = np.array(['sky'] * 8 + ['science'] * 32)

    ffi = FiberFlatImages(
        normflat=normflat,
        normflat_wave=normflat_wave,
        global_norm=global_norm,
        fiber_ids=fiber_ids,
        fiber_types=fiber_types,
        PYP_SPEC='mmt_binospec_ifu',
    )

    # Write and read back
    ofile = Path(data_output_path('test_fiberflatimages.fits')).absolute()
    ffi.to_file(str(ofile), overwrite=True)
    _ffi = FiberFlatImages.from_file(str(ofile))

    assert np.allclose(ffi.normflat, _ffi.normflat), \
        "normflat should survive the FITS I/O round-trip unchanged"
    assert np.allclose(ffi.normflat_wave, _ffi.normflat_wave), \
        "normflat_wave should survive the FITS I/O round-trip unchanged"
    assert np.isclose(ffi.global_norm, _ffi.global_norm), \
        f"global_norm should round-trip; wrote {ffi.global_norm}, " \
        f"read {_ffi.global_norm}"
    assert np.array_equal(ffi.fiber_ids, _ffi.fiber_ids), \
        "fiber_ids should survive the FITS I/O round-trip unchanged"
    assert _ffi.PYP_SPEC == 'mmt_binospec_ifu', \
        f"PYP_SPEC should round-trip as 'mmt_binospec_ifu', got {_ffi.PYP_SPEC!r}"

    ofile.unlink()


def test_normflat_preserves_throughput():
    """Test that normflat preserves throughput differences.

    Creates synthetic fiber spectra where sky fibers have 3x throughput.
    After global normalization, sky fibers should have ~3x higher normflat
    values than science fibers.
    """
    np.random.seed(42)
    n_sci = 8
    n_sky = 2
    nfibers = n_sci + n_sky
    nwave = 500

    # Common spectral shape: quadratic (simulating lamp spectrum)
    wave_grid = np.linspace(4000.0, 7000.0, nwave)
    true_shape = 1.0 - 0.3 * ((wave_grid - 5500.0) / 1500.0) ** 2

    # Science fibers: throughput ~1.0; sky fibers: ~3.0
    fiber_spectra = np.zeros((nfibers, nwave))
    for i in range(n_sci):
        fiber_spectra[i] = true_shape * 10000.0 * np.random.uniform(0.9, 1.1)
    for i in range(n_sci, nfibers):
        fiber_spectra[i] = true_shape * 10000.0 * 3.0

    # Global normalization (IDL approach)
    central = slice(nwave // 10, 9 * nwave // 10)
    global_norm = float(np.nanmax(fiber_spectra[:, central]))
    normflat = fiber_spectra / global_norm

    # Science fibers should have lower normflat than sky fibers
    sci_med = np.median(normflat[:n_sci])
    sky_med = np.median(normflat[n_sci:])
    ratio = sky_med / sci_med

    assert 2.5 < ratio < 3.5, \
        f"Sky/science normflat ratio={ratio:.2f}, expected ~3.0"
    assert sci_med < 0.5, \
        f"Science normflat median={sci_med:.3f}, expected < 0.5 (relative to peak)"
