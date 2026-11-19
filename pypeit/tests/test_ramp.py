"""Tests for the instrument-independent up-the-ramp core (:mod:`pypeit.core.ramp`).

These build ramps directly from numpy -- no spectrograph, no raw-file layout --
to exercise the generic fitting, noise calibration, and sidecar IO in isolation
from any instrument.
"""
import os

import numpy as np
import pytest
from astropy.io import fits

from pypeit.core import ramp
from pypeit.ext.fitramp import fitramp


def synth_ramp(ngroups, ny=64, nx=64, rate=20., sig=8., grptime=2., seed=42):
    """
    Build a synthetic ramp in electrons plus its covariance.

    The interior accumulates ``rate`` e-/s with Poisson shot noise; every read
    carries Gaussian read noise of ``sig`` electrons.  Returns the reads
    (electrons, shape ``(ngroups, ny, nx)``) and the matching Covar.
    """
    rng = np.random.default_rng(seed)
    lam = np.full((ny, nx), rate * grptime)
    increments = rng.poisson(lam, size=(ngroups, ny, nx)).astype(float)
    reads = np.cumsum(increments, axis=0)
    reads += rng.normal(0., sig, size=(ngroups, ny, nx))
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    return reads, covar


def test_ramp_diffs_shape_and_scaling():
    reads, covar = synth_ramp(6, ny=8, nx=8, seed=1)
    diffs = ramp.ramp_diffs(reads, covar)
    assert diffs.shape == (5, 8, 8), \
        'ramp_diffs must drop one group and preserve the spatial shape'
    # First difference equals (reads[1]-reads[0]) / delta_t[0]; assert_allclose
    # raises on mismatch, verifying each difference is scaled by its delta_t.
    expected0 = (reads[1] - reads[0]) / covar.delta_t[0]
    np.testing.assert_allclose(diffs[0], expected0, rtol=1e-12)


def test_calibrate_sigma_recovers_injected_noise():
    sig = 8.
    reads, covar = synth_ramp(15, rate=2., sig=sig, seed=7)
    diffs = ramp.ramp_diffs(reads, covar)
    sig_cal = ramp.calibrate_sigma(diffs, covar)
    assert np.abs(sig_cal - sig) < 1.5, \
        f'calibrated read noise {sig_cal} must recover the injected {sig}'


def test_calibrate_sigma_uncertainty_deterministic():
    reads, covar = synth_ramp(15, rate=2., sig=8., seed=7)
    diffs = ramp.ramp_diffs(reads, covar)
    sig1, err1 = ramp.calibrate_sigma(diffs, covar, return_err=True)
    sig2, err2 = ramp.calibrate_sigma(diffs, covar, return_err=True)
    assert np.isfinite(err1) and err1 > 0., \
        'the bootstrap uncertainty must be finite and positive'
    assert (sig1, err1) == (sig2, err2), \
        'the calibration must be deterministic for a fixed seed'


def test_fit_ramp_recovers_rate():
    rate = 20.
    reads, covar = synth_ramp(8, rate=rate, sig=8., seed=3)
    diffs = ramp.ramp_diffs(reads, covar)
    countrate = ramp.fit_ramp(diffs, covar, 8.)
    interior = countrate[10:-10, 10:-10]
    assert np.abs(np.median(interior) - rate) < 1.0, \
        f'fitted count rate {np.median(interior)} must recover the injected {rate}'


def test_fit_ramp_threaded_matches_serial():
    """Threaded, chunked fitting must be numerically identical to serial."""
    reads, covar = synth_ramp(10, ny=96, nx=96, seed=11)
    diffs = ramp.ramp_diffs(reads, covar)
    serial = ramp.fit_ramp(diffs, covar, 8., workers=1)
    threaded = ramp.fit_ramp(diffs, covar, 8., workers=4)
    assert np.array_equal(serial, threaded, equal_nan=True), \
        'threaded, chunked ramp fitting must be numerically identical to serial'
    s_serial = ramp.calibrate_sigma(diffs, covar, workers=1, nrows=40)
    s_threaded = ramp.calibrate_sigma(diffs, covar, workers=4, nrows=40)
    assert s_serial == s_threaded, \
        'threaded noise calibration must match the serial result exactly'


def test_effective_ronoise_matches_montecarlo():
    """The Brandt (2024a) formula must match the total-count noise of a fitted
    pure-noise ramp."""
    rng = np.random.default_rng(99)
    ngroups, npix, sig, grptime = 20, 2000, 10., 2.
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    reads = rng.normal(0., sig, size=(ngroups, npix))
    diffs = np.diff(reads, axis=0) / np.asarray(covar.delta_t)[:, None]
    result = fitramp.fit_ramps(diffs, covar, np.full(npix, sig),
                               countrateguess=np.zeros(npix))
    measured = np.std(result.countrate * grptime * (ngroups - 1))
    expected = ramp.effective_ronoise(sig, ngroups)
    assert np.abs(measured / expected - 1.) < 0.06, \
        f'the effective read-noise formula ({expected}) must match the ' \
        f'Monte-Carlo total-count noise ({measured})'


def test_rampfit_path_default_and_custom():
    from pathlib import Path
    p = ramp.rampfit_path('/data/raw/sci.0001.fits', '/data/rdx')
    assert p == Path('/data/rdx/RampFit/sci.0001_rampfit.fits'), \
        'the sidecar must default to <redux>/RampFit/<raw stem>_rampfit<ext>'
    p = ramp.rampfit_path('/data/raw/sci.0001.fits', '/data/rdx', 'Ramps')
    assert p == Path('/data/rdx/Ramps/sci.0001_rampfit.fits'), \
        'a custom rampfit_dir must replace the RampFit subdirectory name'


def _minimal_raw(tmp_path, name='sci.fits'):
    """A minimal 2-HDU raw cube (primary + one data ext with metadata)."""
    prihead = fits.Header()
    prihead['INSTRUME'] = 'synth'
    dhead = fits.Header()
    dhead['IMAGETYP'] = 'object'
    dhead['EXPTIME'] = 10.
    hdu = fits.HDUList([fits.PrimaryHDU(header=prihead),
                        fits.ImageHDU(np.zeros((8, 8), dtype='float32'),
                                      header=dhead)])
    path = tmp_path / name
    hdu.writeto(path)
    return path, hdu


def test_write_rampfit_roundtrip_and_freshness(tmp_path):
    raw, hdu = _minimal_raw(tmp_path)
    rate = np.full((8, 8), 3.0)
    sidecar = ramp.rampfit_path(raw, tmp_path)
    ramp.write_rampfit(sidecar, rate, hdu, sig=7.0, eff_ronoise=2.5, ngroups=6,
                       raw_mtime=raw.stat().st_mtime, raw_file=raw)
    hdu.close()
    assert sidecar == tmp_path / 'RampFit' / 'sci_rampfit.fits', \
        'the sidecar must land in the RampFit subdir under the redux path'
    with fits.open(sidecar) as shdu:
        assert shdu[0].header['RAMPFIT'], 'the sidecar must be flagged RAMPFIT'
        assert shdu[0].header['NGROUPS'] == 6, \
            'the sidecar must record the source ramp group count'
        assert np.isclose(shdu[0].header['RAMPSIG'], 7.0)
        assert np.isclose(shdu[0].header['RAMPRON'], 2.5)
        assert shdu[1].header['BUNIT'] == 'e-/s'
        assert shdu[1].header['IMAGETYP'] == 'object', \
            'metadata from the raw data header must be preserved'
        np.testing.assert_allclose(shdu[1].data, rate, rtol=1e-5, atol=1e-6)
    assert ramp.rampfit_fresh(sidecar, raw), \
        'a just-written sidecar must be fresh for its raw source'
    st = raw.stat()
    os.utime(raw, (st.st_atime, st.st_mtime + 10.))
    assert not ramp.rampfit_fresh(sidecar, raw), \
        'bumping the raw mtime must make the sidecar stale'


def test_write_rampfit_atomic_on_failure(tmp_path, monkeypatch):
    """A crash/full-disk mid-writeto must not leave a truncated sidecar."""
    raw, hdu = _minimal_raw(tmp_path)
    sidecar = ramp.rampfit_path(raw, tmp_path)

    def boom(self, *args, **kwargs):
        raise OSError('disk full')
    monkeypatch.setattr(fits.HDUList, 'writeto', boom)
    with pytest.raises(OSError):
        ramp.write_rampfit(sidecar, np.zeros((8, 8)), hdu, sig=7.0,
                           eff_ronoise=2.5, ngroups=6,
                           raw_mtime=raw.stat().st_mtime)
    hdu.close()
    assert not sidecar.exists(), 'a failed write must not leave a partial sidecar'
    assert not list(sidecar.parent.glob('*.tmp*')), \
        'a failed atomic write must leave no temporary file behind'


def test_rampfit_fresh_rejects_mismatched_source(tmp_path):
    """RAWPATH must distinguish a same-named cube from a different directory."""
    (tmp_path / 'n1').mkdir()
    (tmp_path / 'n2').mkdir()
    raw1, hdu1 = _minimal_raw(tmp_path / 'n1', 'sci.0001.fits')
    raw2, _ = _minimal_raw(tmp_path / 'n2', 'sci.0001.fits')
    os.utime(raw2, (raw1.stat().st_atime, raw1.stat().st_mtime))  # identical mtimes
    sidecar = ramp.rampfit_path(raw1, tmp_path)
    ramp.write_rampfit(sidecar, np.zeros((8, 8)), hdu1, sig=7.0, eff_ronoise=2.5,
                       ngroups=6, raw_mtime=raw1.stat().st_mtime, raw_file=raw1)
    hdu1.close()
    assert ramp.rampfit_fresh(sidecar, raw1), \
        'the sidecar must be fresh for the cube it was built from'
    assert not ramp.rampfit_fresh(sidecar, raw2), \
        'RAWPATH must distinguish a same-named cube from another directory'
