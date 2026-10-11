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
    assert not spec.supported
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
    assert par['flexure']['spec_method'] == 'skip'
    assert par['sensfunc']['algorithm'] == 'IR'
    # Boxcar radius of half a lenslet
    cpar = spec.config_specific_par(make_header(HDR_2024))
    assert np.isclose(cpar['reduce']['extraction']['boxcar_radius'], 0.025)
