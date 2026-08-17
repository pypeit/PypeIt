"""Tests for MMT/MMIRS up-the-ramp fitting support."""
import os
from pathlib import Path

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from pypeit import log
from pypeit.ext.fitramp import fitramp
from pypeit.metadata import PypeItMetaData
from pypeit.pypeitsetup import PypeItSetup
from pypeit.scripts.fit_ramp import FitRamp
from pypeit.spectrographs import mmt_mmirs
from pypeit.spectrographs.util import load_spectrograph


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
        h.header['RDNOISE'] = 3.14
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


def test_calibrate_sigma_returns_uncertainty():
    """With return_err=True the calibration also returns a finite, positive
    bootstrap uncertainty, and the point estimate still recovers the truth."""
    ngroups, sig, gain, grptime = 15, 8., 0.95, 2.
    hdu = synth_ramp_hdulist(ngroups, rate=2., sig=sig, gain=gain,
                             grptime=grptime, seed=7)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    reads *= gain
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads, covar)
    sig_cal, sig_err = mmt_mmirs.mmirs_calibrate_sigma(diffs, covar,
                                                       return_err=True)
    assert np.abs(sig_cal - sig) < 1.5
    assert np.isfinite(sig_err) and sig_err > 0.
    # Deterministic: same seed -> identical uncertainty.
    _, sig_err2 = mmt_mmirs.mmirs_calibrate_sigma(diffs, covar, return_err=True)
    assert sig_err == sig_err2


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


def test_fit_ramp_threaded_matches_serial():
    """Threaded, chunked fitting must be numerically identical to the serial
    row-by-row fit for both the count-rate fit and the noise calibration."""
    ngroups, gain, grptime = 10, 0.95, 2.
    # Detector well beyond one chunk (nb=16) so multiple blocks are dispatched
    # (kept square: mmirs_read_amp assumes a square frame).
    hdu = synth_ramp_hdulist(ngroups, ny=96, nx=96, rate=20., sig=8.,
                             gain=gain, grptime=grptime, seed=11)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    reads *= gain
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads, covar)

    serial = mmt_mmirs.mmirs_fit_ramp(diffs, covar, 8., workers=1)
    threaded = mmt_mmirs.mmirs_fit_ramp(diffs, covar, 8., workers=4)
    assert np.array_equal(serial, threaded, equal_nan=True)

    sig_serial = mmt_mmirs.mmirs_calibrate_sigma(diffs, covar, workers=1, nrows=40)
    sig_threaded = mmt_mmirs.mmirs_calibrate_sigma(diffs, covar, workers=4, nrows=40)
    assert sig_serial == sig_threaded


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


def _write_synth(hdulist, path):
    hdulist.writeto(path, overwrite=True)
    return path


def _metadata_for(files, idnames):
    """Minimal PypeItMetaData carrying directory/filename/idname columns."""
    spec = load_spectrograph('mmt_mmirs')
    par = spec.default_pypeit_par()
    data = Table({'filename': [f.name for f in files],
                  'directory': [str(f.parent) for f in files],
                  'idname': idnames})
    fitstbl = PypeItMetaData(spec, par=par, data=data)
    return spec, fitstbl


def test_setup_warns_and_skips_single_hdu_file(tmp_path, monkeypatch):
    raw = _write_synth(synth_ramp_hdulist(4), tmp_path / 'sci.fits')
    bad_files = [tmp_path / 'extra.fits', tmp_path / 'intra.fits']
    for bad_file in bad_files:
        fits.PrimaryHDU(np.zeros((80, 80))).writeto(bad_file)

    warnings = []
    monkeypatch.setattr(log, 'warning',
                        lambda message, *args, **kwargs: warnings.append(message))
    setup = PypeItSetup.from_file_root(tmp_path, 'mmt_mmirs')
    _, _, fitstbl = setup.run(setup_only=True)

    assert fitstbl['filename'].tolist() == [raw.name]
    assert any('Skipping this file' in message for message in warnings)

    with pytest.raises(IndexError):
        PypeItMetaData(setup.spectrograph, par=setup.par,
                       files=[bad_files[0]], strict=True)


def test_cache_metadata_records_darks(tmp_path):
    sci = _write_synth(synth_ramp_hdulist(4, seed=11), tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, seed=12,
                                           imagetyp='dark'),
                        tmp_path / 'dark.fits')
    spec, fitstbl = _metadata_for([sci, dark], ['object', 'dark'])
    # cache_metadata ran inside PypeItMetaData construction (Task 2 hook)
    assert spec._ramp_dark_files == [dark]
    # ... and recorded the reduction directory for RampFit output
    assert spec._ramp_output_dir == Path(fitstbl.par['rdx']['redux_path'])


def test_cache_metadata_applies_rdx_ramp_overrides(tmp_path):
    """[rdx] ramp_fit_cores / ramp_fit_chunk_rows override the class defaults."""
    sci = _write_synth(synth_ramp_hdulist(4, seed=31), tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    default_workers = spec.ramp_fit_workers
    default_chunk = spec.ramp_fit_chunk_rows

    par = spec.default_pypeit_par()
    par['rdx']['ramp_fit_cores'] = 3
    par['rdx']['ramp_fit_chunk_rows'] = 40
    data = Table({'filename': [sci.name], 'directory': [str(sci.parent)],
                  'idname': ['object']})
    # Construction triggers the cache_metadata hook.
    PypeItMetaData(spec, par=par, data=data)

    assert spec.ramp_fit_workers == 3
    assert spec.ramp_fit_chunk_rows == 40
    assert mmt_mmirs.MMTMMIRSSpectrograph.ramp_fit_workers == default_workers
    assert mmt_mmirs.MMTMMIRSSpectrograph.ramp_fit_chunk_rows == default_chunk


def test_cache_metadata_keeps_ramp_defaults_when_unset(tmp_path):
    """With no [rdx] override, the ramp-fit class defaults are left untouched."""
    sci = _write_synth(synth_ramp_hdulist(4, seed=32), tmp_path / 'sci.fits')
    spec, _ = _metadata_for([sci], ['object'])
    assert spec.ramp_fit_workers == mmt_mmirs.MMTMMIRSSpectrograph.ramp_fit_workers
    assert spec.ramp_fit_chunk_rows == mmt_mmirs.MMTMMIRSSpectrograph.ramp_fit_chunk_rows


def test_ramp_sigma_from_dark_and_cached(tmp_path):
    sig_true = 8.
    sci = _write_synth(synth_ramp_hdulist(6, rate=20., sig=sig_true, seed=21),
                       tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true, seed=22,
                                           imagetyp='dark'),
                        tmp_path / 'dark.fits')
    spec, _ = _metadata_for([sci, dark], ['object', 'dark'])

    gain, grptime, ngroups = 0.95, 2., 6
    with fits.open(sci) as hdu:
        reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    sig = spec.get_ramp_sigma(diffs, covar)
    assert np.abs(sig - sig_true) < 1.5
    assert spec._ramp_sigma == sig
    # Cached: works even after the dark file disappears
    dark.unlink()
    assert spec.get_ramp_sigma(diffs, covar) == sig


def test_ramp_sigma_weighted_mean_over_darks(tmp_path, monkeypatch):
    """Read noise is the inverse-variance weighted mean over all darks with
    enough reads, weighted by each dark's calibrated-noise uncertainty."""
    sci = _write_synth(synth_ramp_hdulist(6, seed=201), tmp_path / 'sci.fits')
    d1 = _write_synth(synth_ramp_hdulist(12, imagetyp='dark', seed=202),
                      tmp_path / 'd1.fits')
    d2 = _write_synth(synth_ramp_hdulist(15, imagetyp='dark', seed=203),
                      tmp_path / 'd2.fits')
    spec, _ = _metadata_for([sci, d1, d2], ['object', 'dark', 'dark'])

    # Controlled (sigma, err) per dark, keyed by ndiffs (= nreads - 1).
    controlled = {11: (8.0, 1.0), 14: (6.0, 2.0)}

    def fake_calib(diffs, covar, **kwargs):
        return controlled[diffs.shape[0]]
    monkeypatch.setattr(mmt_mmirs, 'mmirs_calibrate_sigma', fake_calib)

    covar = fitramp.Covar([2. * (i + 1) for i in range(6)])
    sig = spec.get_ramp_sigma(np.zeros((5, 4, 4)), covar)
    w1, w2 = 1. / 1.0 ** 2, 1. / 2.0 ** 2
    expected = (w1 * 8.0 + w2 * 6.0) / (w1 + w2)
    assert np.isclose(sig, expected)
    assert spec._ramp_sigma == sig            # cached


def test_ramp_sigma_reports_max_of_ivar_err_and_sem(tmp_path, monkeypatch):
    """The logged combined uncertainty is max(inverse-variance err, SEM), so
    the larger between-dark scatter is not hidden by tiny within-dark errors."""
    sci = _write_synth(synth_ramp_hdulist(6, seed=221), tmp_path / 'sci.fits')
    d1 = _write_synth(synth_ramp_hdulist(12, imagetyp='dark', seed=222),
                      tmp_path / 'd1.fits')
    d2 = _write_synth(synth_ramp_hdulist(15, imagetyp='dark', seed=223),
                      tmp_path / 'd2.fits')
    spec, _ = _metadata_for([sci, d1, d2], ['object', 'dark', 'dark'])

    controlled = {11: (8.0, 1.0), 14: (6.0, 2.0)}

    def fake_calib(diffs, covar, **kwargs):
        return controlled[diffs.shape[0]]
    monkeypatch.setattr(mmt_mmirs, 'mmirs_calibrate_sigma', fake_calib)

    infos = []
    monkeypatch.setattr(log, 'info',
                        lambda message, *a, **k: infos.append(message))

    covar = fitramp.Covar([2. * (i + 1) for i in range(6)])
    spec.get_ramp_sigma(np.zeros((5, 4, 4)), covar)

    # inverse-variance err = sqrt(1/(1+0.25)) = 0.894; SEM (ddof=1) over
    # {8,6} = sqrt(2)/sqrt(2) = 1.000 -> the SEM wins.
    assert any('+/- 1.00 e-' in m for m in infos)


def test_ramp_sigma_ignores_darks_with_too_few_reads(tmp_path, monkeypatch):
    """A dark with fewer than ramp_min_cal_groups reads is not calibrated
    and does not contribute to the combined read noise."""
    sci = _write_synth(synth_ramp_hdulist(6, seed=211), tmp_path / 'sci.fits')
    good = _write_synth(synth_ramp_hdulist(12, imagetyp='dark', seed=212),
                        tmp_path / 'good.fits')
    few = load_spectrograph('mmt_mmirs').ramp_min_cal_groups - 1
    shallow = _write_synth(synth_ramp_hdulist(few, imagetyp='dark', seed=213),
                           tmp_path / 'shallow.fits')
    spec, _ = _metadata_for([sci, good, shallow], ['object', 'dark', 'dark'])

    seen = []

    def fake_calib(diffs, covar, **kwargs):
        seen.append(diffs.shape[0])
        return (7.0, 1.0)
    monkeypatch.setattr(mmt_mmirs, 'mmirs_calibrate_sigma', fake_calib)

    covar = fitramp.Covar([2. * (i + 1) for i in range(6)])
    sig = spec.get_ramp_sigma(np.zeros((5, 4, 4)), covar)
    assert seen == [11]                       # only the 12-read dark calibrated
    assert np.isclose(sig, 7.0)


def test_ramp_sigma_matches_dark_exptime_to_science(tmp_path, monkeypatch):
    """Only darks whose EXPTIME matches the science frame are calibrated, so
    the effective read noise is measured at the science ramp depth."""
    # synth EXPTIME = grptime*(ngroups-1); grptime defaults to 2 s.
    sci = _write_synth(synth_ramp_hdulist(20, seed=241), tmp_path / 'sci.fits')
    match = _write_synth(synth_ramp_hdulist(20, imagetyp='dark', seed=242),
                         tmp_path / 'match.fits')          # EXPTIME 38 s
    other = _write_synth(synth_ramp_hdulist(12, imagetyp='dark', seed=243),
                         tmp_path / 'other.fits')          # EXPTIME 22 s
    spec, _ = _metadata_for([sci, match, other], ['object', 'dark', 'dark'])

    seen = []

    def fake_calib(diffs, covar, **kwargs):
        seen.append(diffs.shape[0])
        return (9.0, 1.0)
    monkeypatch.setattr(mmt_mmirs, 'mmirs_calibrate_sigma', fake_calib)

    covar = fitramp.Covar([2. * (i + 1) for i in range(20)])
    sig = spec.get_ramp_sigma(np.zeros((19, 4, 4)), covar, exptime=2. * 19)
    assert seen == [19]                       # only the 20-read (matched) dark
    assert np.isclose(sig, 9.0)


def test_ramp_sigma_no_exptime_match_falls_back(tmp_path):
    """Darks are recorded but none match the science EXPTIME -> self-calibrate
    from the frame rather than use a mismatched-depth dark."""
    sig_true = 8.
    gain, grptime = 0.95, 2.
    sci = _write_synth(synth_ramp_hdulist(15, rate=2., sig=sig_true, gain=gain,
                                          grptime=grptime, seed=251),
                       tmp_path / 'sci.fits')              # EXPTIME 28 s
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true,
                                           imagetyp='dark', seed=252),
                        tmp_path / 'dark.fits')            # EXPTIME 22 s
    spec, _ = _metadata_for([sci, dark], ['object', 'dark'])

    reads, head1 = mmt_mmirs.mmirs_load_ramp(fits.open(sci))
    covar = fitramp.Covar([grptime * (i + 1) for i in range(15)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    sig = spec.get_ramp_sigma(diffs, covar, exptime=grptime * 14)
    assert np.abs(sig - sig_true) < 1.5       # self-caled from the science frame
    assert spec._ramp_sigma is None           # self-cal not cached


def test_ramp_sigma_floored_at_header_rdnoise(tmp_path, monkeypatch):
    """A derived noise below the header RDNOISE floor is clamped up to it, since
    values below the instantaneous read noise are unphysical."""
    sci = _write_synth(synth_ramp_hdulist(12, seed=261), tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, imagetyp='dark', seed=262),
                        tmp_path / 'dark.fits')
    spec, _ = _metadata_for([sci, dark], ['object', 'dark'])
    monkeypatch.setattr(mmt_mmirs, 'mmirs_calibrate_sigma',
                        lambda d, c, **k: (2.0, 1.0))     # below the 3.14 floor

    covar = fitramp.Covar([2. * (i + 1) for i in range(12)])
    sig = spec.get_ramp_sigma(np.zeros((11, 4, 4)), covar, exptime=2. * 11,
                              ron_floor=3.14)
    assert np.isclose(sig, 3.14)


def test_ramp_sigma_selfcal_fallback(tmp_path):
    """No darks recorded, but enough groups -> self-calibrate from the frame."""
    sig_true = 8.
    gain, grptime, ngroups = 0.95, 2., 15
    hdu = synth_ramp_hdulist(ngroups, rate=2., sig=sig_true, gain=gain,
                             grptime=grptime, seed=31)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    spec = load_spectrograph('mmt_mmirs')          # fresh: no darks recorded
    sig = spec.get_ramp_sigma(diffs, covar)
    assert np.abs(sig - sig_true) < 1.5
    assert spec._ramp_sigma is None                # self-cal is not cached


def test_ramp_sigma_few_groups_uses_published_guess(tmp_path):
    """Too few reads to self-calibrate and no dark -> the published guess."""
    gain, grptime = 0.95, 2.
    ngroups = load_spectrograph('mmt_mmirs').ramp_min_cal_groups - 1
    hdu = synth_ramp_hdulist(ngroups, rate=2., sig=8., gain=gain,
                             grptime=grptime, seed=131)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    spec = load_spectrograph('mmt_mmirs')          # fresh: no darks recorded
    sig = spec.get_ramp_sigma(diffs, covar)
    # No fit is attempted; the published per-read noise is returned verbatim.
    assert sig == spec.ramp_sig_guess
    assert spec._ramp_sigma is None                # not cached


def test_ramp_sigma_missing_dark_falls_back(tmp_path):
    """A recorded dark that vanishes before first use -> self-cal fallback."""
    sig_true = 8.
    gain, grptime, ngroups = 0.95, 2., 15
    sci = _write_synth(synth_ramp_hdulist(ngroups, rate=2., sig=sig_true,
                                          gain=gain, grptime=grptime, seed=61),
                       tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true, seed=62,
                                           imagetyp='dark'),
                        tmp_path / 'dark.fits')
    spec, _ = _metadata_for([sci, dark], ['object', 'dark'])
    dark.unlink()                     # dark disappears before any sigma call

    with fits.open(sci) as hdu:
        reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    sig = spec.get_ramp_sigma(diffs, covar)
    assert np.abs(sig - sig_true) < 1.5
    assert spec._ramp_sigma is None   # fell back to self-cal, not cached


def test_get_rawimage_ramp_path(tmp_path):
    ngroups, ny, nx, rate, gain, grptime = 6, 64, 64, 20., 0.95, 2.
    exptime = grptime * (ngroups - 1)
    path = _write_synth(synth_ramp_hdulist(ngroups, ny=ny, nx=nx, rate=rate,
                                           gain=gain, grptime=grptime,
                                           seed=41),
                        tmp_path / 'sci_ramp.fits')
    spec = load_spectrograph('mmt_mmirs')
    spec._ramp_output_dir = tmp_path
    detpar, img, hdu, exp, rdsec, oscan = spec.get_rawimage(str(path), 1)
    assert img.shape == (ny - 8, nx - 8)
    assert exp == exptime
    # Effective read noise replaced the CDS default of 3.14
    assert detpar['ronoise'][0] != 3.14
    assert detpar['ronoise'][0] < 15.
    # Rate recovery: image is ADU; convert back to e-/s
    rate_meas = np.median(img) * gain / exptime
    assert np.abs(rate_meas - rate) < 1.5


@pytest.mark.parametrize('ngroups', [2, 4])
def test_get_rawimage_cds_path(tmp_path, ngroups):
    """Frames below ramp_min_reads (5) must use the unchanged CDS path; 4 reads
    is above the old threshold of 3, so this guards the raised minimum."""
    spec = load_spectrograph('mmt_mmirs')
    assert ngroups < spec.ramp_min_reads
    path = _write_synth(synth_ramp_hdulist(ngroups, ny=64, nx=64, rate=20.,
                                           seed=51),
                        tmp_path / 'sci_cds.fits')
    spec._ramp_output_dir = tmp_path
    detpar, img, hdu, exp, rdsec, oscan = spec.get_rawimage(str(path), 1)
    assert img.shape == (56, 56)
    assert detpar['ronoise'][0] == 3.14
    # CDS frames never get a preprocessed sidecar
    assert not mmt_mmirs.mmirs_rampfit_path(path, tmp_path).exists()


def test_findobj_trace_defaults():
    par = load_spectrograph('mmt_mmirs').default_pypeit_par()
    assert par['reduce']['findobj']['trace_npoly'] == 1
    assert par['reduce']['findobj']['trace_maxdev'] == 50.
    assert par['reduce']['findobj']['trace_maxshift'] == 20.


def test_rampfit_path():
    p = mmt_mmirs.mmirs_rampfit_path('/data/raw/sci.0001.fits', '/data/rdx')
    assert p == Path('/data/rdx/RampFit/sci.0001.fits')


def test_write_rampfit_roundtrip(tmp_path):
    ngroups = 6
    raw = _write_synth(synth_ramp_hdulist(ngroups, rate=20., seed=71),
                       tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    with fits.open(raw) as hdu:
        detpar = spec.get_detector_par(1, hdu=hdu)
        rate, sig, eff = spec._ramp_fit_image(hdu, detpar)
        sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
        mmt_mmirs.mmirs_write_rampfit(sidecar, rate, hdu, sig, eff,
                                      raw.stat().st_mtime)
    assert sidecar == tmp_path / 'RampFit' / 'sci.fits'
    assert sidecar.exists()
    with fits.open(sidecar) as shdu:
        assert shdu[0].header['RAMPFIT']
        assert shdu[0].header['NGROUPS'] == ngroups
        assert np.isclose(shdu[0].header['RAMPSIG'], sig)
        assert np.isclose(shdu[0].header['RAMPRON'], eff)
        assert shdu[1].header['BUNIT'] == 'e-/s'
        assert shdu[1].header['DATASEC'] == \
                f'[1:{rate.shape[0]},1:{rate.shape[1]}]'
        # Metadata cards preserved for pypeit_setup on RampFit dirs
        assert shdu[1].header['IMAGETYP'] == 'object'
        assert shdu[1].header['EXPTIME'] == 10.
        np.testing.assert_allclose(shdu[1].data, rate, rtol=1e-5, atol=1e-3)
    assert mmt_mmirs.mmirs_rampfit_fresh(sidecar, raw)
    # Bumping the raw mtime makes the sidecar stale
    st = raw.stat()
    os.utime(raw, (st.st_atime, st.st_mtime + 10.))
    assert not mmt_mmirs.mmirs_rampfit_fresh(sidecar, raw)


def test_write_rampfit_atomic_on_failure(tmp_path, monkeypatch):
    """A crash/full-disk mid-writeto must not leave a truncated sidecar."""
    ngroups = 6
    raw = _write_synth(synth_ramp_hdulist(ngroups, rate=20., seed=75),
                       tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    with fits.open(raw) as hdu:
        detpar = spec.get_detector_par(1, hdu=hdu)
        rate, sig, eff = spec._ramp_fit_image(hdu, detpar)
        sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)

        def boom(self, *args, **kwargs):
            raise OSError('disk full')
        monkeypatch.setattr(fits.HDUList, 'writeto', boom)

        with pytest.raises(OSError):
            mmt_mmirs.mmirs_write_rampfit(sidecar, rate, hdu, sig, eff,
                                          raw.stat().st_mtime)
    assert not sidecar.exists()
    assert sidecar.parent.exists()
    assert not list(sidecar.parent.glob('*.tmp*'))


def test_rampfit_fresh_missing_or_invalid(tmp_path):
    raw = _write_synth(synth_ramp_hdulist(4, seed=72), tmp_path / 'sci.fits')
    sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
    # No sidecar
    assert not mmt_mmirs.mmirs_rampfit_fresh(sidecar, raw)
    # Sidecar without RAWMTIME card
    sidecar.parent.mkdir()
    fits.PrimaryHDU().writeto(sidecar)
    assert not mmt_mmirs.mmirs_rampfit_fresh(sidecar, raw)


def test_rampfit_fresh_missing_raw_file(tmp_path):
    """A sidecar that outlives its raw file must be treated as not fresh."""
    ngroups = 6
    raw = _write_synth(synth_ramp_hdulist(ngroups, rate=20., seed=76),
                       tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    with fits.open(raw) as hdu:
        detpar = spec.get_detector_par(1, hdu=hdu)
        rate, sig, eff = spec._ramp_fit_image(hdu, detpar)
        sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
        mmt_mmirs.mmirs_write_rampfit(sidecar, rate, hdu, sig, eff,
                                      raw.stat().st_mtime)
    raw.unlink()
    assert mmt_mmirs.mmirs_rampfit_fresh(sidecar, raw) is False


def test_count_reads():
    hdu = synth_ramp_hdulist(5, seed=73)
    assert mmt_mmirs.mmirs_count_reads(hdu) == 5


def test_get_rawimage_writes_and_reuses_sidecar(tmp_path, monkeypatch):
    # No _ramp_output_dir recorded: exercises the cwd fallback
    monkeypatch.chdir(tmp_path)
    path = _write_synth(synth_ramp_hdulist(6, rate=20., seed=81),
                        tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    detpar1, img1, hdu1, exp1, _, _ = spec.get_rawimage(str(path), 1)
    sidecar = mmt_mmirs.mmirs_rampfit_path(path, tmp_path)
    assert sidecar.exists()
    assert mmt_mmirs.mmirs_rampfit_fresh(sidecar, path)

    # Second load must come from the sidecar: make refitting impossible
    def boom(*args, **kwargs):
        raise AssertionError('ramp was re-fit despite a fresh sidecar')
    monkeypatch.setattr(mmt_mmirs, 'mmirs_fit_ramp', boom)
    spec2 = load_spectrograph('mmt_mmirs')
    detpar2, img2, hdu2, exp2, _, _ = spec2.get_rawimage(str(path), 1)
    assert exp2 == exp1
    np.testing.assert_allclose(img2, img1, rtol=1e-4, atol=1e-3)
    assert np.isclose(detpar2['ronoise'][0], detpar1['ronoise'][0])
    # Returned hdu is the preprocessed file, with metadata intact
    assert hdu2[0].header['RAMPFIT']
    assert hdu2[1].header['IMAGETYP'] == 'object'


def test_get_rawimage_direct_preprocessed(tmp_path, monkeypatch):
    """A preprocessed file itself (e.g. listed by setup on RampFit) loads."""
    path = _write_synth(synth_ramp_hdulist(6, rate=20., seed=82),
                        tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    spec._ramp_output_dir = tmp_path
    _, img1, *_ = spec.get_rawimage(str(path), 1)
    sidecar = mmt_mmirs.mmirs_rampfit_path(path, tmp_path)

    def boom(*args, **kwargs):
        raise AssertionError('ramp was re-fit for a preprocessed input')
    monkeypatch.setattr(mmt_mmirs, 'mmirs_fit_ramp', boom)
    spec2 = load_spectrograph('mmt_mmirs')
    detpar2, img2, hdu2, exp2, _, _ = spec2.get_rawimage(str(sidecar), 1)
    np.testing.assert_allclose(img2, img1, rtol=1e-4, atol=1e-3)
    assert np.isclose(detpar2['ronoise'][0],
                      fits.getval(sidecar, 'RAMPRON'))


def test_get_rawimage_stale_sidecar_refits(tmp_path):
    path = _write_synth(synth_ramp_hdulist(6, rate=20., seed=83),
                        tmp_path / 'sci.fits')
    spec = load_spectrograph('mmt_mmirs')
    spec._ramp_output_dir = tmp_path
    spec.get_rawimage(str(path), 1)
    sidecar = mmt_mmirs.mmirs_rampfit_path(path, tmp_path)
    old_mtime = fits.getval(sidecar, 'RAWMTIME')

    # Raw cube "changes": sidecar must be refit and rewritten
    st = path.stat()
    os.utime(path, (st.st_atime, st.st_mtime + 10.))
    spec2 = load_spectrograph('mmt_mmirs')
    spec2._ramp_output_dir = tmp_path
    spec2.get_rawimage(str(path), 1)
    assert fits.getval(sidecar, 'RAWMTIME') != old_mtime
    assert mmt_mmirs.mmirs_rampfit_fresh(sidecar, path)


def test_get_rawimage_unwritable_fallback(tmp_path, monkeypatch):
    """Sidecar write failure must warn and continue, not raise."""
    path = _write_synth(synth_ramp_hdulist(6, rate=20., seed=84),
                        tmp_path / 'sci.fits')

    def read_only(*args, **kwargs):
        raise OSError('read-only filesystem')
    monkeypatch.setattr(mmt_mmirs, 'mmirs_write_rampfit', read_only)
    spec = load_spectrograph('mmt_mmirs')
    spec._ramp_output_dir = tmp_path
    detpar, img, *_ = spec.get_rawimage(str(path), 1)
    assert img.shape == (56, 56)
    assert detpar['ronoise'][0] != 3.14
    assert not mmt_mmirs.mmirs_rampfit_path(path, tmp_path).exists()


def test_mmirs_ramp_script(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)          # keep the log file out of the repo
    raw = _write_synth(synth_ramp_hdulist(6, rate=20., seed=91),
                       tmp_path / 'sci.fits')
    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw), '--sig', '8.0']))
    sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
    assert sidecar.exists()
    assert np.isclose(fits.getval(sidecar, 'RAMPSIG'), 8.0)

    # Re-run without --force: fresh sidecar is skipped, file untouched
    mtime_ns = sidecar.stat().st_mtime_ns
    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw), '--sig', '8.0']))
    assert sidecar.stat().st_mtime_ns == mtime_ns

    # --force refits and overwrites
    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw), '--sig', '9.0', '--force']))
    assert np.isclose(fits.getval(sidecar, 'RAMPSIG'), 9.0)


def test_mmirs_ramp_script_skips_few_reads(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    raw = _write_synth(synth_ramp_hdulist(2, seed=92), tmp_path / 'cds.fits')
    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw)]))    # must not raise
    assert not mmt_mmirs.mmirs_rampfit_path(raw, tmp_path).exists()


def test_mmirs_ramp_script_continues_after_write_failure(tmp_path, monkeypatch):
    """A write failure for one file must not abort the rest of the batch."""
    monkeypatch.chdir(tmp_path)
    raw1 = _write_synth(synth_ramp_hdulist(6, rate=20., seed=95),
                        tmp_path / 'sci1.fits')
    raw2 = _write_synth(synth_ramp_hdulist(6, rate=20., seed=96),
                        tmp_path / 'sci2.fits')
    bad_sidecar = mmt_mmirs.mmirs_rampfit_path(raw1, tmp_path)
    real_write_rampfit = mmt_mmirs.mmirs_write_rampfit

    def flaky_write_rampfit(rampfit_file, *args, **kwargs):
        if Path(rampfit_file) == bad_sidecar:
            raise OSError('disk full')
        return real_write_rampfit(rampfit_file, *args, **kwargs)
    monkeypatch.setattr(mmt_mmirs, 'mmirs_write_rampfit', flaky_write_rampfit)

    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw1), str(raw2), '--sig', '8.0',
                                         '--odir', str(tmp_path)]))
    assert not bad_sidecar.exists()
    assert mmt_mmirs.mmirs_rampfit_path(raw2, tmp_path).exists()


def test_mmirs_ramp_script_dark(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    sig_true = 8.
    raw = _write_synth(synth_ramp_hdulist(6, rate=20., sig=sig_true, seed=93),
                       tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true,
                                           seed=94, imagetyp='dark'),
                        tmp_path / 'dark.fits')
    FitRamp.main(FitRamp.parse_args(['mmt_mmirs', str(raw), '--dark', str(dark)]))
    sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
    assert sidecar.exists()
    assert np.abs(fits.getval(sidecar, 'RAMPSIG') - sig_true) < 1.5


def test_fit_ramp_script_rejects_unsupported_spectrograph(tmp_path, monkeypatch):
    """A spectrograph without up-the-ramp support must be rejected up front."""
    from pypeit import PypeItError
    monkeypatch.chdir(tmp_path)
    # The guard fires before any file is opened, so the path need not exist.
    with pytest.raises(PypeItError):
        FitRamp.main(FitRamp.parse_args(['shane_kast_blue',
                                         str(tmp_path / 'nonexistent.fits')]))


# ---------------------------------------------------------------------------
# A-B nod dither support
# ---------------------------------------------------------------------------

def _dither_header(ra_deg, dec_deg, catra_deg, catdec_deg, posang, frameno):
    """Build a 2-element headarr with MMIRS dither cards in ext 1."""
    from astropy.coordinates import Angle
    from astropy import units
    h1 = fits.Header()
    h1['RA'] = str(ra_deg / 15.0)                       # stored in hours
    h1['DEC'] = str(dec_deg)
    h1['CAT-RA'] = Angle(catra_deg, unit=units.deg).to_string(unit=units.deg, sep=':')
    h1['CAT-DEC'] = Angle(catdec_deg, unit=units.deg).to_string(unit=units.deg, sep=':')
    h1['POSANGLE'] = posang
    h1['FILENAME'] = f'MMIRS/2019.0913/nep.as1_mos.{frameno}'
    return [fits.Header(), h1]


def test_dithoff_recovers_abab_pattern():
    spec = load_spectrograph('mmt_mmirs')
    # Catalog target at the nod midpoint; slit PA 166.382 deg (nep.as1 geometry).
    catra, catdec, pa = 260.614, 65.8828, 166.382
    # Four positions offset ALONG the slit by {+1.8, -1.4, +1.4, -1.8} arcsec.
    expected = [1.8, -1.4, 1.4, -1.8]
    par = np.radians(pa)
    cosd = np.cos(np.radians(catdec))
    offs = []
    for off in expected:
        # place pointing at catalog + off arcsec along the slit axis
        dra_as = off * np.sin(par)
        ddec_as = off * np.cos(par)
        ra_deg = catra + (dra_as / 3600.0) / cosd
        dec_deg = catdec + ddec_as / 3600.0
        headarr = _dither_header(ra_deg, dec_deg, catra, catdec, pa, 1822)
        offs.append(spec.compound_meta(headarr, 'dithoff'))
    assert np.allclose(offs, expected, atol=0.02)


def test_frameno_parsed_from_filename():
    spec = load_spectrograph('mmt_mmirs')
    headarr = _dither_header(260.614, 65.8828, 260.614, 65.8828, 166.382, 1823)
    assert spec.compound_meta(headarr, 'frameno') == 1823


def test_dithoff_missing_card_returns_zero():
    spec = load_spectrograph('mmt_mmirs')
    headarr = [fits.Header(), fits.Header()]   # no dither cards
    assert spec.compound_meta(headarr, 'dithoff') == 0.0


def _nod_table(dithoffs, frametype='science', setup='A'):
    """Minimal table mimicking what set_combination_groups passes to
    get_comb_group: unique comb_id per frame, bkg_id=-1, increasing mjd."""
    n = len(dithoffs)
    return Table({
        'frametype': [frametype] * n,
        'setup': [setup] * n,
        'dithoff': np.array(dithoffs, dtype=float),
        'mjd': 59000.0 + np.arange(n) / 1440.0,   # 1-min cadence, time-ordered
        'comb_id': np.arange(n, dtype=int) + 1,
        'bkg_id': np.full(n, -1, dtype=int),
        # init_meta seeds these from the 4-char default 'None', so the real
        # metadata table carries them as width-4 unicode; mirror that here so
        # the tests catch truncation of longer labels (e.g. "ABA'B'").
        'dithpat': np.full(n, 'None'),
        'dithpos': np.full(n, 'None'),
    })


def test_get_comb_group_pairs_ababprime():
    spec = load_spectrograph('mmt_mmirs')
    tbl = _nod_table([1.8, -1.4, 1.4, -1.8, 1.8, -1.4, 1.4, -1.8])
    out = spec.get_comb_group(tbl)
    # Greedy sequential pairing across the midpoint (0.0): (0,1),(2,3),(4,5),(6,7)
    assert out['bkg_id'].tolist() == [2, 1, 4, 3, 6, 5, 8, 7]
    assert out['dithpos'].tolist() == ["A", "B", "A'", "B'", "A", "B", "A'", "B'"]
    assert set(out['dithpat']) == {"ABA'B'"}


def test_get_comb_group_stare_no_pairing():
    spec = load_spectrograph('mmt_mmirs')
    # Peak-to-peak 0.3" < nod_min_offset (1.0") -> treated as a stare.
    tbl = _nod_table([0.0, 0.1, -0.1, 0.2])
    out = spec.get_comb_group(tbl)
    assert out['bkg_id'].tolist() == [-1, -1, -1, -1]


def test_get_comb_group_odd_frame_unpaired():
    spec = load_spectrograph('mmt_mmirs')
    tbl = _nod_table([1.8, -1.4, 1.8])
    out = spec.get_comb_group(tbl)
    # (0,1) pair; frame 2 has no remaining opposite-nod partner.
    assert out['bkg_id'].tolist() == [2, 1, -1]


def test_get_comb_group_longslit_simple_ab():
    spec = load_spectrograph('mmt_mmirs')
    # No mask/decker columns present -> logic is slit-count agnostic.
    tbl = _nod_table([2.5, -2.5])
    out = spec.get_comb_group(tbl)
    assert out['bkg_id'].tolist() == [2, 1]
    assert out['dithpos'].tolist() == ["A", "B"]


def test_pypeit_file_keys_include_dither_columns():
    spec = load_spectrograph('mmt_mmirs')
    keys = spec.pypeit_file_keys()
    for col in ['dithpat', 'dithpos', 'dithoff', 'frameno']:
        assert col in keys
