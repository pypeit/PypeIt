"""
Unit tests for the Keck/OSIRIS spectrograph skeleton.

Synthetic headers reproduce the card sets of the four eras/flavours of
OSIRIS raw files recorded in the development-suite reports
(``pypeitdev/keck_osiris/Reports/01*`` and ``02*``):

- 2015 (Hawaii-2 detector, ITIME in s)
- 2016 truncated header (ISSKY, PA_SPEC, DTMP7, weather cards missing)
- 2016 full header (no INSTRUME, no weather cards)
- 2024 (complete header, ITIME in ms)
"""
import numpy as np
import pytest

from astropy.io import fits
from astropy.table import Table

from pypeit import PypeItError
from pypeit.spectrographs.keck_osiris import KeckOSIRISSpectrograph
from pypeit.spectrographs.util import load_spectrograph


# Header cards common to all synthetic headers
COMMON = dict(RA=187.27791667, DEC=2.05241667, AIRMASS=1.06, COADDS=1,
              ROTPOSN=312.0, INSTANGL=312.0, ROTREFAN=0.0, PARANG=-10.0,
              RAOFF=0.0, DECOFF=0.0)

HDR_2015 = dict(COMMON, INSTRUME='OSIRIS', CURRINST='OSIRIS',
                OBSTYPE='astro', ISSKY=1, SFILTER='Kn5', SSCALE='0.035',
                ITIME=600.0, SAMPMODE=1, NUMREADS=32, MJD_OBS=57173.465035,
                TARGNAME='v5558 sgr', PA_SPEC=0.0, DTMP7=59.2,
                WXPRESS=618.3, WXOUTTMP=3.0, WXOUTHUM=8.4)

HDR_2016_TRUNC = dict(COMMON, CURRINST='OSIRIS', OBSTYPE='a',
                      SFILTER='Hn3', SSCALE='0.100', ITIME=240000,
                      SAMPMODE=3, NUMREADS=162, MJD_OBS=57468.428963,
                      TARGNAME='3C 273')

HDR_2016_FULL = dict(HDR_2016_TRUNC, ISSKY=0, PA_SPEC=0, DTMP7=55.3,
                     MJD_OBS=57468.443104)

HDR_2024 = dict(COMMON, INSTRUME='OSIRIS', CURRINST='OSIRIS', OBSTYPE='a',
                ISSKY=0, SFILTER='Kbb', SSCALE='0.050', ITIME=2000,
                TRUITIME=1.475, SAMPMODE=3, NUMREADS=1,
                MJD_OBS=60343.480002,
                TARGNAME='HD89239(A0V)', ROTPOSN=487.0,
                INSTANGL=311.99999999, PA_SPEC=175, DTMP7=54.19,
                WXPRESS=617.57, WXOUTTMP=0.5, WXOUTHUM=30.0)


def make_header(cards):
    """Build a FITS header from a dict (MJD_OBS -> MJD-OBS)."""
    hdr = fits.Header()
    for key, value in cards.items():
        hdr[key.replace('MJD_OBS', 'MJD-OBS')] = value
    return hdr


def write_raw(path, cards, shape=(8, 8), value=2.0):
    """Write a tiny OSIRIS-like raw file (SCI, VAR, DQ)."""
    hdr = make_header(cards)
    sci = fits.PrimaryHDU(data=np.full(shape, value, dtype=np.float32),
                          header=hdr)
    sci.header['EXTNAME'] = 'SCI'
    var = fits.ImageHDU(data=np.ones(shape, dtype=np.float32), name='VAR')
    dq = np.full(shape, 9, dtype=np.uint8)
    dq[0, 0] = 0
    dqhdu = fits.ImageHDU(data=dq, name='DQ')
    fits.HDUList([sci, var, dqhdu]).writeto(path, overwrite=True)
    return path


@pytest.fixture(scope='module')
def spec():
    return load_spectrograph('keck_osiris')


def test_load_spectrograph(spec):
    assert isinstance(spec, KeckOSIRISSpectrograph)
    assert spec.name == 'keck_osiris'
    assert spec.pypeline == 'Fiber'
    assert spec.supported
    assert spec.camera == 'OSIRIS'
    assert spec.configuration_keys() == ['dispname', 'decker']
    for key in ['issky', 'idname', 'comb_id', 'bkg_id']:
        assert key in spec.pypeit_file_keys()


def test_exptime_units(spec):
    # 2015: ITIME in seconds
    assert spec.get_meta_value(make_header(HDR_2015), 'exptime') == 600.
    # 2016: ITIME in milliseconds
    for cards in [HDR_2016_TRUNC, HDR_2016_FULL]:
        assert spec.get_meta_value(make_header(cards), 'exptime') == 240.
    # 2024: TRUITIME (s per coadd) takes precedence over ITIME
    assert spec.get_meta_value(make_header(HDR_2024), 'exptime') == 1.475
    # COADDS multiplies
    assert np.isclose(spec.get_meta_value(make_header(dict(HDR_2024,
                                                           COADDS=3)),
                                          'exptime'), 3 * 1.475)
    # Without TRUITIME, a 2 s exposure has ITIME = 2000 ms; the era
    # decides the units
    notru = {k: v for k, v in HDR_2024.items() if k != 'TRUITIME'}
    assert spec.get_meta_value(make_header(notru), 'exptime') == 2.
    # Without MJD-OBS, ITIME > 2000 flags milliseconds
    nomjd = {k: v for k, v in HDR_2016_FULL.items() if k != 'MJD_OBS'}
    assert KeckOSIRISSpectrograph.exptime_from_header(nomjd) == 240.


def test_metadata_tolerance(spec):
    trunc = make_header(HDR_2016_TRUNC)
    # Missing INSTRUME -> CURRINST; missing ISSKY -> 0
    assert spec.get_meta_value(trunc, 'instrument') == 'OSIRIS'
    assert spec.get_meta_value(trunc, 'issky') == 0
    # Weather defaults (Maunakea)
    assert spec.get_meta_value(trunc, 'pressure') == 611.
    assert spec.get_meta_value(trunc, 'temperature') == 1.5
    assert spec.get_meta_value(trunc, 'humidity') == 20.
    # Missing DTMP7 is tolerated
    assert spec.get_meta_value(trunc, 'tmatemp') is None
    # OBSTYPE first letter
    assert spec.get_meta_value(make_header(HDR_2015), 'idname') == 'a'
    # Position angle from the rotator (matches PA_SPEC)
    assert np.isclose(spec.get_meta_value(make_header(HDR_2024), 'posang'),
                      175., atol=1e-3)
    assert np.isclose(spec.get_meta_value(trunc, 'posang'), 0.)
    # Configuration values are strings
    assert spec.get_meta_value(trunc, 'decker') == '0.100'
    assert spec.get_meta_value(make_header(HDR_2024), 'dispname') == 'Kbb'
    # Lenslet size in degrees
    assert np.isclose(spec.get_meta_value(trunc, 'slitwid'), 0.1/3600)
    # Full 2024 header uses the recorded weather
    assert spec.get_meta_value(make_header(HDR_2024), 'pressure') == 617.57
    assert spec.get_meta_value(make_header(HDR_2024), 'tmatemp') == 54.19


def test_frame_typing(spec):
    # Build a metadata-like table from synthetic headers
    cards = [
        HDR_2016_FULL,                                         # science
        dict(HDR_2016_FULL, ISSKY=1),                          # sky
        dict(HDR_2024, SFILTER='Drk', TARGNAME='unknown',
             SSCALE='0.020'),                                  # dark
        HDR_2024,                                              # standard
        dict(HDR_2024, ISSKY=1),                               # std sky
        dict(HDR_2024, TARGNAME='NGC3245', ITIME=600000,
             TRUITIME=598.963),                                # science
        dict(HDR_2016_FULL, OBSTYPE='s', TARGNAME='star'),     # standard
        dict(HDR_2016_FULL, OBSTYPE='c'),                      # calib
    ]
    keys = ['target', 'dispname', 'exptime', 'idname', 'issky', 'mjd']
    rows = [spec.get_meta_value(make_header(c), keys) for c in cards]
    tbl = Table(rows=rows, names=keys)

    sci = spec.check_frame_type('science', tbl)
    assert np.array_equal(sci, [1, 1, 0, 0, 0, 1, 0, 0])
    for ftype in ['trace', 'arc', 'tilt']:
        assert np.array_equal(spec.check_frame_type(ftype, tbl), sci)
    assert np.array_equal(spec.check_frame_type('sky', tbl),
                          [0, 1, 0, 0, 1, 0, 0, 0])
    assert np.array_equal(spec.check_frame_type('dark', tbl),
                          [0, 0, 1, 0, 0, 0, 0, 0])
    assert np.array_equal(spec.check_frame_type('standard', tbl),
                          [0, 0, 0, 1, 1, 0, 1, 0])
    for ftype in ['bias', 'pixelflat', 'illumflat', 'align', 'pinhole',
                  'scattlight']:
        assert not np.any(spec.check_frame_type(ftype, tbl))


def test_standard_names(spec):
    names = ['HD101060(A0V)', 'HD 89239', 'NGC3245', '3C 273',
             'NGC3245-PSF1', 'HIP 12345', 'v5558 sgr']
    texp = [10., 2., 600., 240., 30., 600., 600.]
    assert np.array_equal(spec.is_standard_target(names, texp),
                          [1, 1, 0, 0, 0, 0, 0])


def test_read_noise(spec):
    assert spec.read_noise(2, 1) == 22.
    assert spec.read_noise(3, 1) == 22.
    assert np.isclose(spec.read_noise(3, 4), 12.2)
    assert np.isclose(spec.read_noise(3, 16), 8.1)
    assert np.isclose(spec.read_noise(3, 32), 7.8)
    assert np.isclose(spec.read_noise(3, 162), 7.8)
    assert spec.read_noise(None, None) == 8.


def test_era_rejection(spec, tmp_path):
    # Pre-2016 frames are refused
    old = write_raw(tmp_path / 's150531_a025001.fits', HDR_2015)
    with pytest.raises(PypeItError, match='H2RG'):
        spec.get_rawimage(old, 1)
    with pytest.raises(PypeItError):
        spec.check_era(None)


def test_get_rawimage(spec, tmp_path):
    raw = write_raw(tmp_path / 's160321_a002011.fits', HDR_2016_FULL)
    detpar, img, hdu, exptime, datasec, oscan = spec.get_rawimage(raw, 1)
    # ADU/s/coadd -> ADU
    assert exptime == 240.
    assert np.allclose(img, 2.0 * 240.)
    assert np.all(datasec == 1) and np.all(oscan == 0)
    assert detpar['specaxis'] == 1
    assert np.isclose(detpar['ronoise'][0], 7.8)
    assert detpar['platescale'] == 0.1
    # VAR and DQ remain available
    var, dq = spec.read_var_dq(hdu)
    assert var.shape == img.shape and dq[0, 0] == 0
    hdu.close()


def test_static_masks(spec):
    hot, dead = spec.load_static_masks()
    assert hot.shape == (2048, 2048) and dead.shape == (2048, 2048)
    assert hot.sum() == 4002
    assert dead.sum() == 1045
    # BPM from the shape alone: static masks, re-oriented
    bpm = spec.bpm(None, 1, shape=(2048, 2048))
    assert bpm.sum() == np.sum(hot | dead)
    detpar = spec.get_detector_par(1)
    assert np.array_equal(bpm.astype(bool),
                          spec.orient_image(detpar, hot | dead))


def test_wave_coeffs(spec):
    epochs = spec.wave_coeff_epochs()
    assert len(epochs) == 11
    # 2016-03-21 (Hn3 test data)
    assert spec.wave_coeff_file(57468.43).name \
            == 'osiris_wave_coeffs_160101_170508.fits'
    # Inclusive end date: late on 2017-05-08 is still the 2016 epoch
    assert spec.wave_coeff_file(57881.9).name \
            == 'osiris_wave_coeffs_160101_170508.fits'
    # 2024 (KOA data) -> open-ended latest epoch
    assert spec.wave_coeff_file(60343.5).name == 'osiris_wave_coeffs.fits'
    with pytest.raises(PypeItError):
        spec.wave_coeff_file(50000.)
    coeffs = spec.load_wave_coeffs(60343.5)
    assert coeffs.shape == (66, 51, 4)
    # Every shipped epoch file exists
    for f in epochs['file']:
        assert spec.spec_data_path().get_file_path(f).is_file()


def test_grating_expansion(spec):
    cte = spec.load_al_cte()
    assert len(cte) == 46
    assert np.isclose(spec.grating_expansion(spec.TEMP_REF_WAVE), 1.)
    assert spec.grating_expansion(None) == 1.
    # Warmer grating -> larger expansion
    assert spec.grating_expansion(60.) > 1.


def test_lenslet_maps(spec):
    s2l = spec.load_slice_to_lenslet()
    assert s2l.shape == (1216, 3)
    assert np.array_equal(np.sort(s2l[:, 0]), np.arange(1216))
    lmap = spec.load_lenslet_mapping()
    assert lmap.shape == (64, 19)
    # The two maps agree (rows in the slice map start at 2)
    assert np.all(lmap[s2l[:, 1] - 2, s2l[:, 2]] == s2l[:, 0])


def test_default_par(spec):
    par = spec.default_pypeit_par()
    proc = par['scienceframe']['process']
    assert not proc['use_darkimage']
    assert not proc['use_biasimage']
    assert not proc['use_overscan']
    assert not proc['use_pixelflat']
    assert not proc['mask_cr']
    assert par['calibrations']['flatfield']['method'] == 'skip'
    # M3: Keck-coefficient wavelengths, raw-frame cosmic-ray masking
    assert par['calibrations']['wavelengths']['method'] == 'fixed'
    assert par['calibrations']['osiris_rect_crmask'] == 'lacosmic'
    assert par['calibrations']['osiris_rect_method'] == 'lstsq'
    assert not par['reduce']['skysub']['joint_fit']
    assert spec.has_rectification
    assert par['flexure']['spec_method'] == 'skip'
    assert par['sensfunc']['algorithm'] == 'IR'
    # Boxcar radius of half a lenslet
    cpar = spec.config_specific_par(make_header(HDR_2024))
    assert np.isclose(cpar['reduce']['extraction']['boxcar_radius'], 0.025)
    # M4: per-lenslet OH refinement, OH list per band, no DTMP7 term,
    # Davies OH scaling of A-B residuals
    assert par['calibrations']['wavelengths']['fiber_sky_refine'] \
        == 'zeropoint'
    assert cpar['calibrations']['wavelengths']['lamps'] == ['OH_OSIRIS']
    cpar = spec.config_specific_par(make_header(HDR_2016_FULL))
    assert cpar['calibrations']['wavelengths']['lamps'] == ['OH_NIRES']
    assert not par['calibrations']['osiris_grating_temp']
    assert par['reduce']['skysub']['osiris_oh_scaling'] == 'global'
    assert spec.frame_expansion(make_header(HDR_2024)) == 1.0
    assert spec.frame_expansion(make_header(HDR_2024), apply=True) > 1.0


# ---------------------------------------------------------------------------
# Milestone M3: lenslet layout, fiber hooks, identity tilts, Keck wavelengths
# (Generated by JXP and Claude)
# ---------------------------------------------------------------------------
MJD_2016 = 57468.43
MJD_2024 = 60343.48


def test_wave_inversion_roundtrip(spec):
    # Inverting the Keck pixel(lambda) relation of the shipped coefficients
    # and re-evaluating it returns the detector columns
    from pypeit.core import lenslet_rectify as lr
    for mjd, filt in [(MJD_2016, 'Hn3'), (MJD_2024, 'Kbb')]:
        coeffs = spec.load_wave_coeffs(mjd)
        lay = lr.lenslet_layout(filt, coeffs)
        order = lr.keck_order(filt)
        sel = np.linspace(0, len(lay) - 1, 25).astype(int)
        cl = coeffs[lay['crow'][sel], lay['ccol'][sel]]
        for frac in [1.0, spec.grating_expansion(55.3)]:
            lam = lr.keck_wavelengths(cl, order, frac=frac)
            x = lr.keck_pixels(lam, cl, order, frac=frac)
            # Within the lenslet windows (narrowband) / whole row
            for k, i in enumerate(sel):
                w = slice(lay['xlo'][i], lay['xhi'][i] + 1)
                assert np.max(np.abs(x[k, w] - np.arange(2048)[w])) < 1e-3
            assert np.all(np.diff(lam, axis=1) < 0)


def test_lenslet_layout_broadband(spec):
    from pypeit.core import lenslet_rectify as lr
    lay = lr.lenslet_layout('Kbb', spec.load_wave_coeffs(MJD_2024))
    assert lay.meta['broadband'] and lay.meta['nblock'] == 19
    assert lay.meta['nspat'] == 1216
    # Identity: spatial pixel = slice; lenslet (row+1, col+16) of the grid
    assert np.array_equal(lay['pos'], lay['slice'])
    assert np.array_equal(lay['crow'], lay['slice'] % 64 + 1)
    assert np.array_equal(lay['ccol'], lay['slice'] // 64 + 16)
    assert np.array_equal(lay['lens_id'], lay['crow'] * 51 + lay['ccol'])
    # Broadband: the window is the filter band (DRP cube range)
    width = lay['xhi'] - lay['xlo'] + 1
    assert np.all(lay['xhi'] - lay['xlo'] < 2047)
    assert 1400 < np.median(width) < 1700
    assert len(lay) > 1150


def test_lenslet_layout_narrowband(spec):
    from pypeit.core import lenslet_rectify as lr
    lay = lr.lenslet_layout('Hn3', spec.load_wave_coeffs(MJD_2016))
    assert not lay.meta['broadband']
    assert lay.meta['nblock'] == 51 and lay.meta['nspat'] == 64 * 51
    assert len(np.unique(lay['pos'])) == len(lay)
    assert len(np.unique(lay['lens_id'])) == len(lay)
    # One Keck lenslet column per block
    for b in range(lay.meta['nblock']):
        assert len(np.unique(lay['ccol'][lay['bidx'] == b])) == 1
    # Spectra sharing a slice have disjoint raw-column windows of about
    # one narrowband spectrum (~390 px) plus the gaps
    for sp in [0, 500, 900]:
        idx = np.where(lay['slice'] == sp)[0]
        idx = idx[np.argsort(lay['xlo'][idx])]
        assert len(idx) >= 3
        assert np.all(lay['xlo'][idx[1:]] > lay['xhi'][idx[:-1]])
    width = lay['xhi'] - lay['xlo'] + 1
    assert 400 < np.median(width) < 520
    # No lenslet clipped onto row 65 twice
    assert np.all(lay['crow'] <= 65)
    # Reference lenslets sit in their blocks
    ref = lr.block_reference(lay)
    assert np.array_equal(lay['bidx'][ref], np.arange(lay.meta['nblock']))


def test_fiber_hooks(spec):
    spec.set_lenslet_state('Hn3', MJD_2016, scale=0.1)
    blocks = spec.get_fiber_blocks(1)
    assert len(blocks) == 51
    assert all(b['box_radius'] == 0.5 for b in blocks)
    pos = np.concatenate([b['fiber_positions'] for b in blocks])
    ids = np.concatenate([b['fiber_ids'] for b in blocks])
    assert len(pos) == len(spec.lenslet_layout())
    # Positions -> ids round trip, also via the cube builder interface
    meta = spec.get_fiber_metadata(1, pos.astype(int), slit_centers=pos)
    assert np.array_equal(meta['fiber_id'], ids)
    assert np.all(meta['fiber_type'] == 'SCI')
    res = spec.identify_fibers_in_block(1, 7, blocks[7]['fiber_positions'])
    assert np.array_equal(res['fiber_id'], blocks[7]['fiber_ids'])
    # A position of another block is not matched
    res = spec.identify_fibers_in_block(1, 7, blocks[8]['fiber_positions'][:2])
    assert np.all(res['fiber_id'] == -1)
    lidx = spec.get_science_fiber_layout_indices(1, meta['fiber_id'],
                                                 meta['fiber_type'])
    x, y = spec.load_sky_layout()
    assert x.size == 66 * 51 and np.all(lidx < x.size)
    # Lenslet scale 0.1": neighbouring columns are 0.1" apart in x
    assert np.isclose(x[1] - x[0], -0.1)
    # Block-slit edges
    left, right = spec.get_block_slit_edges(np.zeros((10, 64 * 51)), 1)
    assert left.shape == (10, 51)
    assert np.allclose(right[0, :-1] + 1, left[0, 1:] + 1)
    assert np.isclose(right[0, 0] - left[0, 0], 64.)
    assert left[0, 0] == -0.5
    with pytest.raises(PypeItError):
        spec.get_block_slit_edges(np.zeros((10, 100)), 1)


def _block_slits(spec, nspec=2048):
    """Block-slits of the current lenslet layout (helper)."""
    from pypeit import slittrace
    left, right = spec.get_block_slit_edges(
        np.zeros((nspec, spec.lenslet_layout().meta['nspat'])), 1)
    return slittrace.SlitTraceSet(left, right, 'Fiber', nspec=nspec,
                                  nspat=spec.lenslet_layout().meta['nspat'],
                                  PYP_SPEC=spec.name)


def test_fixed_tilts(spec):
    # Keck-derived tilts: one 2D Legendre per block, finite and close to
    # the identity (the lenslet stagger is a few tens of pixels at most)
    spec.set_lenslet_state('Kbb', MJD_2024, scale=0.05)
    slits = _block_slits(spec)
    tilts = spec.get_fixed_tilts(slits, raw_header=make_header(HDR_2024))
    assert tilts.nslit == slits.nslits
    assert tilts.coeffs.shape == (spec.TILT_ORDER[0] + 1,
                                  spec.TILT_ORDER[1] + 1, slits.nslits)
    assert np.any(tilts.coeffs[:, 1:, :] != 0)
    slitmask = slits.slit_img()
    img = tilts.fit2tiltimg(slitmask)
    on = slitmask > -1
    ident = np.tile(np.arange(2048)[:, None] / 2047., (1, slitmask.shape[1]))
    assert np.all(np.isfinite(img[on]))
    assert np.max(np.abs(img[on] - ident[on])) * 2047 < 80.
    # The base class provides no fixed tilts
    assert load_spectrograph('shane_kast_blue').get_fixed_tilts(slits) is None


@pytest.mark.parametrize('cards,scale', [(HDR_2024, 0.05),
                                         (HDR_2016_FULL, 0.1)])
def test_tilt_waveimg_vs_keck(spec, cards, scale):
    # The wavelength image built from the Keck-derived tilts and the block
    # solutions reproduces every lenslet's Keck solution inside the band
    # to < 0.05 pixel rms (M4 task 8)
    from pypeit.core import lenslet_rectify as lr
    hdr = make_header(cards)
    spec.set_lenslet_state(hdr['SFILTER'], hdr['MJD-OBS'], scale=scale)
    slits = _block_slits(spec)
    par = spec.default_pypeit_par()['calibrations']['wavelengths']
    wv = spec.get_fixed_wavecalib(slits, hdr, par)
    tilts = spec.get_fixed_tilts(slits, raw_header=hdr)
    slitmask = slits.slit_img()
    waveimg = wv.build_waveimg(tilts.fit2tiltimg(slitmask), slits)
    lay = spec.lenslet_layout()
    lam = spec.lenslet_wavelengths(spec.frame_expansion(hdr))
    band = lr.lenslet_band_mask(lam, lay)[:, ::-1]
    exact = lam[:, ::-1] * 10.
    model = waveimg[:, np.asarray(lay['pos'])].T
    disp = np.abs(np.gradient(exact, axis=1))
    dpix = ((model - exact) / disp)[band]
    rms = np.sqrt(np.mean(dpix**2))
    assert rms < 0.05
    assert 1.4826 * np.median(np.abs(dpix)) < 0.03


def test_fixed_wavecalib(spec):
    # One WaveFit per block that reproduces the Keck solution of the
    # block's reference lenslet; per-lenslet wavelengths from the hook
    from pypeit.core import lenslet_rectify as lr
    from pypeit import specobjs, specobj
    hdr = make_header(HDR_2016_FULL)
    spec.set_lenslet_state('Hn3', hdr['MJD-OBS'], scale=0.1)
    slits = _block_slits(spec)
    par = spec.default_pypeit_par()['calibrations']['wavelengths']
    assert par['method'] == 'fixed'
    wv = spec.get_fixed_wavecalib(slits, hdr, par)
    lay = spec.lenslet_layout()
    ref = lr.block_reference(lay)
    lam = spec.lenslet_wavelengths(spec.frame_expansion(hdr))
    for b in [0, 20, 50]:
        w = lam[ref[b]][::-1] * 10.
        assert np.max(np.abs(wv.wv_fits[b].wave_soln - w)
                       / np.abs(np.gradient(w))) < 0.1
    assert np.all([f.rms < 0.05 for f in wv.wv_fits])
    # Hn3 band (vacuum, Angstrom) at the reference lenslet window centre
    i = ref[25]
    xc = int(lay['xcen'][i])
    assert 15900. < lam[i][xc] * 10 < 16800.
    # Per-lenslet wavelengths for two lenslets
    sobjs = specobjs.SpecObjs()
    for k in [0, 1000]:
        so = specobj.SpecObj(PYPELINE='Fiber', DET='DET01', OBJTYPE='science',
                             SLITID=1)
        so.MASKDEF_ID = int(lay['lens_id'][k])
        sobjs.add_sobj(so)
    waves = spec.get_fiber_wavelengths(1, sobjs, [hdr])
    assert waves.shape == (2, 2048)
    assert np.allclose(waves[1], lam[1000][::-1] * 10.)


def test_rectify_step(spec, tmp_path, monkeypatch):
    # Generic RawImage.rectify step with a stand-in transformation: the
    # image, variance components and BPM are replaced, and the variance of
    # the result is the returned (flux-independent) variance only
    from pypeit.images import rawimage
    raw = write_raw(tmp_path / 'raw.fits', HDR_2016_FULL, shape=(64, 64),
                    value=2.0)

    def fake_rectify(image, ivar, bpm, headarr, par, base_var=None,
                     rn2img=None, filename=None):
        assert image.shape == (64, 64) and ivar.shape == (64, 64)
        assert np.all(base_var > 0) and rn2img is not None
        out = image[:, ::2] + image[:, 1::2]
        return dict(image=out, rn2img=np.full(out.shape, 4.),
                    proc_var=np.full(out.shape, 12.),
                    bpm=np.zeros(out.shape, dtype=bool))

    monkeypatch.setattr(spec, 'rectify', fake_rectify)
    par = spec.default_pypeit_par()['scienceframe']['process']
    img = rawimage.RawImage(str(raw), spec, 1).process(par)
    assert 'rectify' in img.process_steps
    assert img.image.shape == (64, 32)
    assert np.allclose(img.ivar, 1 / 16.)
    assert not img.shot_noise
    # Two raw pixels of 2 ADU/s x 240 s x 2.15 e-/ADU, minus the tabulated
    # dark current (90 e-/hr)
    assert np.allclose(img.image, 2 * (2.0 * 240. * 2.15 - 90. * 240 / 3600))
    # Darks (and biases) are not rectified: they are subtracted from the
    # raw frames before the rectification (M4 task 2)
    img = rawimage.RawImage(str(raw), spec, 1).process(par, rectify=False)
    assert 'rectify' not in img.process_steps
    assert img.image.shape == (64, 64)


# ---------------------------------------------------------------------------
# Milestone M4: fiber cube hooks (Generated by JXP and Claude)
# ---------------------------------------------------------------------------
def test_fiber_cube_setup(spec):
    # The cube hook sets the lenslet state from a spec1d-like header
    hdr = fits.Header({'SFILTER': 'Kbb', 'MJD': MJD_2024, 'SSCALE': '0.050',
                       'DECKER': '0.05'})
    out = spec.fiber_cube_setup(hdr)
    assert out['wcs_ref'] == 'origin'
    assert np.isclose(out['spatial_scale'], 0.05)
    assert spec.dispname == 'Kbb'
    assert spec.lenslet_layout().meta['nblock'] == 19
    # Base class: no state, no defaults
    assert load_spectrograph('shane_kast_blue').fiber_cube_setup(hdr) == {}


def test_fiber_cube_grid(spec, tmp_path):
    # A cube built from lenslet spectra on the regular lenslet grid holds
    # each lenslet's spectrum at its own spaxel, and the WCS reference pixel
    # (wcs_ref='origin') is the DRP reference spaxel (layout origin)
    from pypeit.core import datacube
    hdr = fits.Header({'SFILTER': 'Kbb', 'MJD': MJD_2024, 'SSCALE': '0.050',
                       'RA': 150., 'DEC': 20., 'ROTPOSN': 312.,
                       'INSTANGL': 312.})
    spec.fiber_cube_setup(hdr)
    x, y = spec.load_sky_layout()
    lay = spec.lenslet_layout()
    ids = np.asarray(lay['lens_id'])
    nspec = 50
    wave = np.tile(20000. + 2.5 * np.arange(nspec), (ids.size, 1))
    flux = np.tile(ids[:, None].astype(float), (1, nspec))
    meta = dict(fiber_id=ids, fiber_name=np.asarray(lay['name']),
                fiber_type=np.full(ids.size, 'SCI'))
    data = {'DET01': dict(flux=flux, ivar=np.ones_like(flux), wave=wave,
                          fiber_meta=meta)}
    out = tmp_path / 'cube.fits'
    datacube.build_cube_common(data, spec, x, y, hdr, 'spec1d_test.fits',
                               spatial_scale=0.05, method='nearest',
                               output=str(out), wcs_ref='origin')
    with fits.open(out) as hdu:
        cube = hdu['FLUX'].data
        h = hdu['FLUX'].header
    # 19 lenslet columns + 2 padding pixels, 64 rows + 2
    assert cube.shape[1:] == (66, 21)
    # Lenslet id at the reference pixel: column 24, flipped row 31
    ix, iy = int(round(h['CRPIX1'])) - 1, int(round(h['CRPIX2'])) - 1
    assert np.isclose(h['CRPIX1'] - 1 - ix, 0.) and h['WCSREF'] == 'origin'
    crow = 66 - 1 - 31
    assert np.allclose(cube[:, iy, ix], crow * 51 + 24)


# ---------------------------------------------------------------------------
# Milestone M5: robustness (truncated headers, header eras, matrices) and
# the cached Keck-derived tilts (Generated by JXP and Claude)
# ---------------------------------------------------------------------------
def test_truncated_header_hooks(spec):
    # The 2016 truncated header (no ISSKY, PA_SPEC, DTMP7, INSTRUME or
    # weather cards) goes through every header-dependent hook
    hdr = make_header(HDR_2016_TRUNC)
    # No DTMP7: no grating expansion, whatever the parameter says
    assert spec.frame_expansion(hdr, apply=True) == 1.0
    assert spec.frame_expansion(hdr, apply=False) == 1.0
    spec.set_lenslet_state('Hn3', hdr['MJD-OBS'], scale=0.1)
    slits = _block_slits(spec)
    par = spec.default_pypeit_par()['calibrations']['wavelengths']
    wv = spec.get_fixed_wavecalib(slits, hdr, par)
    assert len(wv.wv_fits) == slits.nslits == 51
    tilts = spec.get_fixed_tilts(slits, raw_header=hdr)
    assert tilts.nslit == 51
    # Cube defaults from the reduced-file header cards
    out = spec.fiber_cube_setup(fits.Header({'SFILTER': 'Hn3',
                                             'MJD': hdr['MJD-OBS']}))
    assert np.isclose(out['spatial_scale'], 0.05)    # SSCALE default
    # All metadata keys of the pypeit file parse
    for key in spec.pypeit_file_keys():
        if key in ['filename', 'frametype', 'calib', 'comb_id', 'bkg_id',
                   'setup', 'directory']:
            continue
        spec.get_meta_value(hdr, key)


def test_header_missing_cards(spec):
    # Headers without the filter or the date give a clear error naming the
    # missing cards, not a KeyError
    for drop in ['SFILTER', 'MJD_OBS']:
        cards = {k: v for k, v in HDR_2024.items() if k != drop}
        with pytest.raises(PypeItError, match=drop.replace('_', '-')):
            spec._set_state_from_header(make_header(cards))


def test_era_boundary(spec):
    # The H2RG era starts on 2016-01-01 (MJD 57388) inclusive
    spec.check_era(KeckOSIRISSpectrograph.MJD_H2RG)
    spec.check_era(60343.5)
    with pytest.raises(PypeItError, match='s150531_a025001'):
        spec.check_era(KeckOSIRISSpectrograph.MJD_H2RG - 0.01,
                       filename='/raw/s150531_a025001.fits')
    with pytest.raises(PypeItError, match='MJD-OBS'):
        spec.check_era(np.nan)


def test_missing_matrix_errors(spec, tmp_path):
    # No matrix directory: the error says how to set it
    hdr = make_header(HDR_2024)
    cpar = spec.default_pypeit_par()['calibrations']
    cpar['osiris_rectmat_dir'] = str(tmp_path / 'missing')
    with pytest.raises(PypeItError, match='osiris_rectmat_dir'):
        spec.rect_matrix(hdr, cpar)
    # A directory with the wrong matrices and a misnamed file: the error
    # lists both
    (tmp_path / 's160318_c005___infl_Hn3_100.fits').touch()
    (tmp_path / 'Kbb_050_matrix.fits').touch()
    cpar['osiris_rectmat_dir'] = str(tmp_path)
    with pytest.raises(PypeItError) as err:
        spec.rect_matrix(hdr, cpar)
    msg = str(err.value)
    assert 's160318_c005___infl_Hn3_100.fits' in msg
    assert 'Kbb_050_matrix.fits' in msg
    assert '230701-250310' in msg and 'tkserver' in msg


def test_fixed_tilts_cache(spec):
    # The Keck-derived tilts are computed once per configuration and the
    # returned objects do not share their coefficients with the cache
    hdr = make_header(HDR_2024)
    spec.set_lenslet_state('Kbb', hdr['MJD-OBS'], scale=0.05)
    slits = _block_slits(spec)
    t1 = spec.get_fixed_tilts(slits, raw_header=hdr)
    n = len(KeckOSIRISSpectrograph._tilt_cache)
    t1.coeffs[:] = 0.
    t2 = load_spectrograph('keck_osiris').get_fixed_tilts(slits,
                                                          raw_header=hdr)
    assert len(KeckOSIRISSpectrograph._tilt_cache) == n
    assert np.any(t2.coeffs != 0.)
