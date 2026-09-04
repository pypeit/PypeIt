"""Tests for the vendored up-the-ramp fitting module (pypeit.ext.fitramp.fitramp)."""
import numpy as np

from pypeit.ext.fitramp import fitramp


def make_ramp(ngroups, npix, rate, sig, grptime, rng):
    """Synthetic ramp reads (electrons): cumulative Poisson signal + Gaussian read noise."""
    increments = rng.poisson(rate * grptime, size=(ngroups, npix)).astype(float)
    reads = np.cumsum(increments, axis=0)
    reads += rng.normal(0., sig, size=(ngroups, npix))
    return reads


def test_slope_recovery():
    rng = np.random.default_rng(42)
    ngroups, npix, rate, sig, grptime = 20, 500, 50., 10., 2.
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    reads = make_ramp(ngroups, npix, rate, sig, grptime, rng)
    diffs = np.diff(reads, axis=0) / np.asarray(covar.delta_t)[:, None]
    sig_arr = np.full(npix, sig)
    result = fitramp.fit_ramps(diffs, covar, sig_arr)
    assert np.abs(np.mean(result.countrate) - rate) < 1.0
    # chisq should be ~ndiffs-1 per pixel
    assert np.abs(np.median(result.chisq) - (ngroups - 2)) < 3.0


def test_jump_detection_and_masked_fit():
    rng = np.random.default_rng(12345)
    ngroups, npix, rate, sig, grptime = 20, 200, 5., 10., 2.
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    reads = make_ramp(ngroups, npix, rate, sig, grptime, rng)
    # Inject a 5000 e- cosmic-ray hit at read index 10 in the first 50 pixels
    reads[10:, :50] += 5000.
    diffs = np.diff(reads, axis=0) / np.asarray(covar.delta_t)[:, None]
    sig_arr = np.full(npix, sig)
    diffs2use, countrates = fitramp.mask_jumps(diffs, covar, sig_arr)
    # The jumped diff (index 9) must be masked for all affected pixels
    assert np.all(diffs2use[9, :50] == 0)
    # Unaffected pixels should be almost entirely unmasked
    assert np.mean(diffs2use[:, 50:]) > 0.95
    # Refit with the mask recovers the true rate for the hit pixels
    result = fitramp.fit_ramps(diffs, covar, sig_arr, diffs2use=diffs2use,
                               countrateguess=countrates * (countrates > 0))
    assert np.abs(np.median(result.countrate[:50]) - rate) < 2.0
