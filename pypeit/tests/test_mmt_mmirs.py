"""Tests for MMT/MMIRS up-the-ramp fitting support."""
import numpy as np
from astropy.io import fits

from pypeit.core import fitramp
from pypeit.spectrographs import mmt_mmirs


def synth_ramp_hdulist(ngroups, ny=64, nx=64, rate=20., sig=8., gain=0.95,
                       grptime=2., seed=42, imagetyp='object'):
    """
    Build an MMIRS-like raw HDUList.

    Mimics the real file layout: ext 1 is the FINAL read (highest EXTVER) and
    carries the metadata; exts 2..ngroups hold EXTVER 1..ngroups-1.  The outer
    4-pixel border is signal-free (reference pixels); the interior accumulates
    ``rate`` e-/s.  Reads are in ADU (divided by ``gain``).
    """
    rng = np.random.default_rng(seed)
    lam = np.zeros((ny, nx))
    lam[4:-4, 4:-4] = rate * grptime          # e- per group interval
    increments = rng.poisson(lam, size=(ngroups, ny, nx)).astype(float)
    reads = np.cumsum(increments, axis=0)
    reads += rng.normal(0., sig, size=(ngroups, ny, nx))
    reads /= gain                              # electrons -> ADU
    exptime = grptime * (ngroups - 1)

    hdus = [fits.PrimaryHDU()]
    order = [ngroups - 1] + list(range(ngroups - 1))   # final read first
    for i in order:
        h = fits.ImageHDU(reads[i].astype('float32'))
        h.header['EXTVER'] = i + 1
        h.header['DATASEC'] = f'[5:{ny-4},5:{nx-4}]'
        h.header['CCDSUM'] = '1 1'
        h.header['FILTER'] = 'zJ'
        h.header['DISPERSE'] = 'J'
        h.header['GAIN'] = gain
        h.header['GRPTIME'] = grptime
        h.header['EXPTIME'] = exptime
        h.header['IMAGETYP'] = imagetyp
        h.header['HXRGGRUP'] = ngroups
        h.header['INSTRUME'] = 'mmirs'
        h.header['RA'] = 180.0
        h.header['DEC'] = 30.0
        h.header['OBJECT'] = 'synthetic'
        h.header['APERTURE'] = '1pixel-long'
        h.header['AIRMASS'] = 1.2
        h.header['DATE-OBS'] = '2017-04-09T06:36:25'
        hdus.append(h)
    return fits.HDUList(hdus)


def test_load_ramp_orders_and_trims():
    ngroups, ny, nx = 5, 64, 64
    hdu = synth_ramp_hdulist(ngroups, ny=ny, nx=nx, sig=0.5, seed=1)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    assert reads.shape == (ngroups, ny - 8, nx - 8)
    assert head1['HXRGGRUP'] == ngroups
    # Reads must be in time order: interior signal strictly accumulates
    interior_means = reads[:, 10:-10, 10:-10].mean(axis=(1, 2))
    assert np.all(np.diff(interior_means) > 0)


def test_calibrate_sigma_recovers_truth():
    ngroups, sig, gain, grptime = 15, 8., 0.95, 2.
    hdu = synth_ramp_hdulist(ngroups, rate=2., sig=sig, gain=gain,
                             grptime=grptime, seed=7)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    reads *= gain                              # ADU -> electrons
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads, covar)
    sig_cal = mmt_mmirs.mmirs_calibrate_sigma(diffs, covar)
    assert np.abs(sig_cal - sig) < 1.5


def test_fit_ramp_recovers_rate():
    ngroups, rate, sig, gain, grptime = 8, 20., 8., 0.95, 2.
    hdu = synth_ramp_hdulist(ngroups, rate=rate, sig=sig, gain=gain,
                             grptime=grptime, seed=3)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    reads *= gain
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads, covar)
    countrate = mmt_mmirs.mmirs_fit_ramp(diffs, covar, sig)
    interior = countrate[10:-10, 10:-10]
    assert np.abs(np.median(interior) - rate) < 1.0


def test_effective_ronoise_formula():
    """Monte-Carlo check: total-count noise of a fitted pure-noise ramp."""
    rng = np.random.default_rng(99)
    ngroups, npix, sig, grptime = 20, 2000, 10., 2.
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    reads = rng.normal(0., sig, size=(ngroups, npix))
    diffs = np.diff(reads, axis=0) / np.asarray(covar.delta_t)[:, None]
    sig_arr = np.full(npix, sig)
    result = fitramp.fit_ramps(diffs, covar, sig_arr,
                               countrateguess=np.zeros(npix))
    total_exp = grptime * (ngroups - 1)
    measured = np.std(result.countrate * total_exp)
    expected = mmt_mmirs.mmirs_effective_ronoise(sig, ngroups)
    assert np.abs(measured / expected - 1.) < 0.06
