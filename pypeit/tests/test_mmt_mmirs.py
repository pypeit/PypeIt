"""Tests for MMT/MMIRS up-the-ramp fitting support."""
import os
from pathlib import Path

import numpy as np
import pytest
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


def _write_synth(hdulist, path):
    hdulist.writeto(path, overwrite=True)
    return path


def _metadata_for(files, idnames):
    """Minimal PypeItMetaData carrying directory/filename/idname columns."""
    from astropy.table import Table
    from pypeit.metadata import PypeItMetaData
    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('mmt_mmirs')
    par = spec.default_pypeit_par()
    data = Table({'filename': [f.name for f in files],
                  'directory': [str(f.parent) for f in files],
                  'idname': idnames})
    fitstbl = PypeItMetaData(spec, par=par, data=data)
    return spec, fitstbl


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


def test_ramp_sigma_from_dark_and_cached(tmp_path):
    sig_true = 8.
    sci = _write_synth(synth_ramp_hdulist(6, rate=20., sig=sig_true, seed=21),
                       tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true, seed=22,
                                           imagetyp='dark'),
                        tmp_path / 'dark.fits')
    spec, _ = _metadata_for([sci, dark], ['object', 'dark'])

    gain, grptime, ngroups = 0.95, 2., 6
    from astropy.io import fits as _fits
    with _fits.open(sci) as hdu:
        reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    sig = spec.get_ramp_sigma(diffs, covar)
    assert np.abs(sig - sig_true) < 1.5
    assert spec._ramp_sigma == sig
    # Cached: works even after the dark file disappears
    dark.unlink()
    assert spec.get_ramp_sigma(diffs, covar) == sig


def test_ramp_sigma_selfcal_fallback(tmp_path):
    """No darks recorded (or too few groups) -> self-calibrate from the frame."""
    sig_true = 8.
    gain, grptime, ngroups = 0.95, 2., 15
    hdu = synth_ramp_hdulist(ngroups, rate=2., sig=sig_true, gain=gain,
                             grptime=grptime, seed=31)
    reads, head1 = mmt_mmirs.mmirs_load_ramp(hdu)
    covar = fitramp.Covar([grptime * (i + 1) for i in range(ngroups)])
    diffs = mmt_mmirs.mmirs_ramp_diffs(reads * gain, covar)

    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('mmt_mmirs')          # fresh: no darks recorded
    sig = spec.get_ramp_sigma(diffs, covar)
    assert np.abs(sig - sig_true) < 1.5
    assert spec._ramp_sigma is None                # self-cal is not cached


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

    from astropy.io import fits as _fits
    with _fits.open(sci) as hdu:
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
    from pypeit.spectrographs.util import load_spectrograph
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


def test_get_rawimage_cds_path(tmp_path):
    """2-read files must still use the unchanged CDS path."""
    path = _write_synth(synth_ramp_hdulist(2, ny=64, nx=64, rate=20., seed=51),
                        tmp_path / 'sci_cds.fits')
    from pypeit.spectrographs.util import load_spectrograph
    spec = load_spectrograph('mmt_mmirs')
    spec._ramp_output_dir = tmp_path
    detpar, img, hdu, exp, rdsec, oscan = spec.get_rawimage(str(path), 1)
    assert img.shape == (56, 56)
    assert detpar['ronoise'][0] == 3.14
    # CDS frames never get a preprocessed sidecar
    assert not mmt_mmirs.mmirs_rampfit_path(path, tmp_path).exists()


def test_findobj_trace_defaults():
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.spectrographs.util import load_spectrograph
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
    from pypeit.scripts.mmirs_ramp import MMIRSRamp
    monkeypatch.chdir(tmp_path)          # keep the log file out of the repo
    raw = _write_synth(synth_ramp_hdulist(6, rate=20., seed=91),
                       tmp_path / 'sci.fits')
    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw), '--sig', '8.0']))
    sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
    assert sidecar.exists()
    assert np.isclose(fits.getval(sidecar, 'RAMPSIG'), 8.0)

    # Re-run without --force: fresh sidecar is skipped, file untouched
    mtime_ns = sidecar.stat().st_mtime_ns
    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw), '--sig', '8.0']))
    assert sidecar.stat().st_mtime_ns == mtime_ns

    # --force refits and overwrites
    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw), '--sig', '9.0', '--force']))
    assert np.isclose(fits.getval(sidecar, 'RAMPSIG'), 9.0)


def test_mmirs_ramp_script_skips_few_reads(tmp_path, monkeypatch):
    from pypeit.scripts.mmirs_ramp import MMIRSRamp
    monkeypatch.chdir(tmp_path)
    raw = _write_synth(synth_ramp_hdulist(2, seed=92), tmp_path / 'cds.fits')
    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw)]))    # must not raise
    assert not mmt_mmirs.mmirs_rampfit_path(raw, tmp_path).exists()


def test_mmirs_ramp_script_continues_after_write_failure(tmp_path, monkeypatch):
    """A write failure for one file must not abort the rest of the batch."""
    from pypeit.scripts.mmirs_ramp import MMIRSRamp
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

    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw1), str(raw2), '--sig', '8.0',
                                         '--odir', str(tmp_path)]))
    assert not bad_sidecar.exists()
    assert mmt_mmirs.mmirs_rampfit_path(raw2, tmp_path).exists()


def test_mmirs_ramp_script_dark(tmp_path, monkeypatch):
    from pypeit.scripts.mmirs_ramp import MMIRSRamp
    monkeypatch.chdir(tmp_path)
    sig_true = 8.
    raw = _write_synth(synth_ramp_hdulist(6, rate=20., sig=sig_true, seed=93),
                       tmp_path / 'sci.fits')
    dark = _write_synth(synth_ramp_hdulist(12, rate=0.1, sig=sig_true,
                                           seed=94, imagetyp='dark'),
                        tmp_path / 'dark.fits')
    MMIRSRamp.main(MMIRSRamp.parse_args([str(raw), '--dark', str(dark)]))
    sidecar = mmt_mmirs.mmirs_rampfit_path(raw, tmp_path)
    assert sidecar.exists()
    assert np.abs(fits.getval(sidecar, 'RAMPSIG') - sig_true) < 1.5
