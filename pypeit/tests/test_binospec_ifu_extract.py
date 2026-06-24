"""Unit tests for pypeit_binospec_ifu_extract helpers."""
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io.fits import Header

from pypeit import PypeItError
from pypeit.onespec import OneSpec
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph
from pypeit.scripts.binospec_ifu_extract import BinospecIFUExtract
from pypeit.core.datacube import (
    _SKY_LINE_MASK_HALFWIDTH,
    _SKY_LINE_WAVELENGTHS,
    compute_fiber_fluxes,
    load_fibers,
    project_to_sky,
    resample_and_combine,
    sky_line_mask,
    write_onespec,
)

# Real IFU spectrograph supplies the shared TAN/POSANG WCS convention used by
# project_to_sky.
_IFU_SPEC = MMTBINOSPECIFUSpectrograph()


def test_sky_line_mask_marks_each_line():
    wave = np.array(_SKY_LINE_WAVELENGTHS, dtype=float)
    mask = sky_line_mask(wave)
    assert mask.shape == wave.shape
    assert np.all(mask), 'Each sky-line wavelength itself should be masked'


def test_sky_line_mask_excludes_far_wavelengths():
    # 5000 Angstrom is well clear of every line in the table
    wave = np.array([5000.0, 5577.34, 5000.0], dtype=float)
    mask = sky_line_mask(wave)
    assert mask.tolist() == [False, True, False]


def test_sky_line_mask_halfwidth_boundary():
    line = _SKY_LINE_WAVELENGTHS[0]
    # Just inside ± halfwidth → masked; just outside → not masked.
    inside = line + _SKY_LINE_MASK_HALFWIDTH - 0.01
    outside = line + _SKY_LINE_MASK_HALFWIDTH + 0.01
    mask = sky_line_mask(np.array([inside, outside]))
    assert mask.tolist() == [True, False]


def _fake_header(ra=180.0, dec=30.0, posang=0.0):
    h = Header()
    h['RA'] = ra
    h['DEC'] = dec
    h['POSANG'] = posang
    return h


def test_project_to_sky_origin():
    """Fiber at offset (0,0) should land on the pointing."""
    hdr = _fake_header(ra=180.0, dec=30.0, posang=0.0)
    ra, dec = project_to_sky(np.array([0.0]), np.array([0.0]), hdr, _IFU_SPEC)
    assert ra.shape == (1,)
    assert dec.shape == (1,)
    np.testing.assert_allclose(ra, 180.0, atol=1e-9)
    np.testing.assert_allclose(dec, 30.0, atol=1e-9)


def test_project_to_sky_posang_zero():
    """With POSANG=0, +y should move Dec north (positive)."""
    hdr = _fake_header(ra=180.0, dec=30.0, posang=0.0)
    ra, dec = project_to_sky(np.array([0.0]), np.array([1.0]), hdr, _IFU_SPEC)
    # 1 arcsec north of the pointing
    np.testing.assert_allclose(dec, 30.0 + 1.0 / 3600.0, atol=1e-7)
    np.testing.assert_allclose(ra, 180.0, atol=1e-7)


def test_project_to_sky_posang_90():
    """POSANG=90 rotates the layout: +y in instrument frame → +RA (east)."""
    hdr = _fake_header(ra=180.0, dec=0.0, posang=90.0)
    ra, dec = project_to_sky(np.array([0.0]), np.array([1.0]), hdr, _IFU_SPEC)
    # PA=90 means +y points East; at dec=0, +1 arcsec east → +1/3600 deg in RA.
    np.testing.assert_allclose(ra - 180.0, 1.0 / 3600.0, atol=1e-7)
    np.testing.assert_allclose(dec, 0.0, atol=1e-7)


def test_project_to_sky_handles_string_ra_dec():
    """Match binospec_ifu_cube: header may have RA in hourangle or degrees."""
    hdr = Header()
    hdr['RA'] = '12:00:00'   # 180 deg
    hdr['DEC'] = '+30:00:00'
    hdr['POSANG'] = 0.0
    ra, dec = project_to_sky(np.array([0.0]), np.array([0.0]), hdr, _IFU_SPEC)
    np.testing.assert_allclose(ra, 180.0, atol=1e-6)
    np.testing.assert_allclose(dec, 30.0, atol=1e-6)


def test_project_to_sky_missing_posang_defaults_to_zero():
    hdr = Header()
    hdr['RA'] = 180.0
    hdr['DEC'] = 30.0
    # No POSANG key.
    ra, dec = project_to_sky(np.array([0.0]), np.array([1.0]), hdr, _IFU_SPEC)
    np.testing.assert_allclose(dec, 30.0 + 1.0 / 3600.0, atol=1e-7)


def test_compute_fiber_fluxes_basic():
    """Sum over each fiber's wavelength range, ignoring sky lines."""
    # Two fibers, one wavelength grid, identical flat flux of 1.0 -> integral
    # equals the count of unmasked pixels in [wave_min, wave_max].
    wave = np.linspace(4000.0, 5000.0, 1001)
    flux = np.ones((2, wave.size))
    waves = [wave, wave]
    fluxes = [flux[0], flux[1]]

    out = compute_fiber_fluxes(waves, fluxes,
                                wave_min=4000.0, wave_max=5000.0)
    assert out.shape == (2,)
    # No sky lines fall in 4000-5000 Angstrom, so every pixel is summed.
    np.testing.assert_allclose(out, [wave.size, wave.size])


def test_compute_fiber_fluxes_excludes_sky_line():
    """A fiber with a tall spike at 5577 should have most flux masked out."""
    wave = np.linspace(5500.0, 5700.0, 2001)
    flat = np.ones_like(wave)
    spiky = flat.copy()
    near_5577 = np.abs(wave - 5577.34) < 0.5
    spiky[near_5577] = 1000.0  # huge sky-line residual

    out = compute_fiber_fluxes([wave, wave], [flat, spiky],
                                wave_min=5500.0, wave_max=5700.0)
    # The masked region knocks out the spike entirely; the two fibers' integrals
    # should be equal (within floating-point) because the spike sits inside the
    # ±10 Angstrom mask of the 5577 line.
    np.testing.assert_allclose(out[0], out[1])


def test_compute_fiber_fluxes_handles_nan():
    """nansum skips NaN pixels; result equals the count of finite unmasked pixels."""
    wave = np.linspace(4000.0, 4100.0, 11)
    flux = np.ones_like(wave)
    flux[5] = np.nan
    out = compute_fiber_fluxes([wave], [flux],
                                wave_min=4000.0, wave_max=4100.0)
    np.testing.assert_allclose(out, [10.0])  # 11 pixels - 1 NaN


def test_compute_fiber_fluxes_per_fiber_wave_grids():
    """Each fiber may have its own native wavelength grid."""
    w0 = np.linspace(4000.0, 4100.0, 11)
    w1 = np.linspace(4050.0, 4150.0, 11)
    f0 = np.ones_like(w0)
    f1 = np.ones_like(w1) * 2.0
    # Slider range is narrower than both grids; each fiber clips independently
    # to its own pixel positions.
    out = compute_fiber_fluxes([w0, w1], [f0, f1],
                                wave_min=4060.0, wave_max=4090.0)
    # Both fibers fully inside [4060, 4090]: w0 has 4 pixels, w1 has 4 pixels
    # (4060,4070,4080,4090).  Integrals: 4 and 8.
    np.testing.assert_allclose(out, [4.0, 8.0])


def test_resample_and_combine_two_identical_fibers():
    """Two identical fibers → flux doubles, ivar doubles."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave) * 4.0  # var = 0.25
    out_w, out_f, out_iv = resample_and_combine([wave, wave],
                                                  [flux, flux],
                                                  [ivar, ivar])
    np.testing.assert_allclose(out_f, 2.0 * np.ones_like(out_w))
    # var_out = sum(var_i) = 0.25 + 0.25 = 0.5 → ivar_out = 2.0
    np.testing.assert_allclose(out_iv, 2.0 * np.ones_like(out_w))


def test_resample_and_combine_offset_grids():
    """Two fibers with offset native grids combine over their overlap."""
    w0 = np.linspace(4000.0, 4100.0, 21)
    w1 = np.linspace(4050.0, 4150.0, 21)
    f0 = np.ones_like(w0)
    f1 = np.ones_like(w1)
    iv0 = np.ones_like(w0)
    iv1 = np.ones_like(w1)
    out_w, out_f, out_iv = resample_and_combine([w0, w1], [f0, f1],
                                                  [iv0, iv1])
    # The output grid is restricted to the overlap [4050, 4100].
    assert out_w[0] >= 4050.0
    assert out_w[-1] <= 4100.0
    # Both fibers contribute on the overlap → flux ≈ 2
    np.testing.assert_allclose(out_f, 2.0 * np.ones_like(out_w), atol=1e-6)


def test_resample_and_combine_no_overlap_raises():
    w0 = np.linspace(4000.0, 4100.0, 11)
    w1 = np.linspace(5000.0, 5100.0, 11)
    f = np.ones(11)
    iv = np.ones(11)
    with pytest.raises(PypeItError):
        resample_and_combine([w0, w1], [f, f], [iv, iv])


def test_resample_and_combine_single_fiber():
    wave = np.linspace(4000.0, 5000.0, 51)
    flux = 3.0 * np.ones_like(wave)
    ivar = 4.0 * np.ones_like(wave)
    out_w, out_f, out_iv = resample_and_combine([wave], [flux], [ivar])
    np.testing.assert_allclose(out_f, flux)
    np.testing.assert_allclose(out_iv, ivar)


def test_resample_and_combine_zero_ivar_pixels():
    """Pixels with ivar=0 should be excluded from variance combine."""
    wave = np.linspace(4000.0, 4100.0, 11)
    flux = np.ones_like(wave)
    iv0 = np.ones_like(wave)
    iv1 = np.ones_like(wave)
    iv1[5] = 0.0  # bad pixel in fiber 1
    _, out_f, out_iv = resample_and_combine([wave, wave], [flux, flux],
                                              [iv0, iv1])
    # At pixel 5 only fiber 0 contributes ivar.  var_out = 1 → ivar_out = 1.
    np.testing.assert_allclose(out_iv[5], 1.0)
    # Elsewhere both contribute: var = 1+1 = 2, ivar = 0.5
    np.testing.assert_allclose(out_iv[0], 0.5)


def test_resample_and_combine_all_zero_wavelengths_raises():
    """Fibers with no valid wavelength data raise PypeItError."""
    with pytest.raises(PypeItError, match="No valid wavelength data"):
        resample_and_combine([np.zeros(10)], [np.ones(10)],
                              [np.ones(10)])


# ---------------------------------------------------------------------------
# Fixtures for load_fibers tests
# ---------------------------------------------------------------------------

class _FakeSpec(SimpleNamespace):
    """Just enough of the SpecObj API for load_fibers to read."""


class _FakeSpecObjs(list):
    @property
    def DET(self):
        return np.array([s.DET for s in self], dtype=object)

    @property
    def nobj(self):
        return len(self)

    def __getitem__(self, key):
        if isinstance(key, np.ndarray) and key.dtype == bool:
            return _FakeSpecObjs([s for s, k in zip(self, key) if k])
        return super().__getitem__(key)


def _make_spec(det, slitid, spat, wave, flux, ivar):
    return _FakeSpec(DET=det, SLITID=slitid, SPAT_PIXPOS=spat,
                     OPT_WAVE=wave, OPT_COUNTS=flux, OPT_COUNTS_IVAR=ivar,
                     BOX_WAVE=wave, BOX_COUNTS=flux, BOX_COUNTS_IVAR=ivar)


class _FakeSpectrograph:
    def __init__(self, fiber_types_by_det):
        self._types = fiber_types_by_det

    def get_fiber_metadata(self, det, slit_spat_ids, slit_centers=None):
        types = np.array(self._types[det])
        assert len(types) == len(slit_spat_ids), (
            f'Fixture misconfiguration: {len(types)} types for '
            f'{len(slit_spat_ids)} fibers on det {det}')
        return {'fiber_id': np.arange(len(slit_spat_ids)),
                'fiber_type': types}

    def get_science_fiber_layout_indices(self, det, fiber_ids, fiber_types):
        # Sci fibers map to consecutive indices; sky → -1.
        out = np.full_like(fiber_ids, -1)
        sci = fiber_types != 'SKY'
        # Fake distinct layout indices per detector.
        offset = 0 if det == 1 else 100
        out[sci] = offset + np.arange(np.sum(sci))
        return out


# ---------------------------------------------------------------------------
# load_fibers tests
# ---------------------------------------------------------------------------

def test_load_fibers_drops_sky_and_returns_native_arrays():
    wave = np.linspace(4000.0, 5000.0, 11)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave) * 4.0
    sobjs = _FakeSpecObjs([
        _make_spec('DET01', 1, 100.0, wave, flux, ivar),
        _make_spec('DET01', 2, 200.0, wave, flux, ivar),
        _make_spec('DET02', 1, 100.0, wave, flux, ivar),
    ])
    spec = _FakeSpectrograph(fiber_types_by_det={1: ['SCI', 'SKY'],
                                                  2: ['SCI']})
    # Two distinct layouts.  Index 0 = (0,0), index 100 = (1,0) arcsec.
    targetx = np.zeros(200)
    targety = np.zeros(200)
    targetx[0] = 0.0
    targety[0] = 0.0
    targetx[100] = 1.0
    targety[100] = 0.0

    fibers = load_fibers(sobjs, spec, targetx, targety, prefix='OPT')
    # Two science fibers survive: DET01 #0 and DET02 #0.
    assert len(fibers) == 2
    assert fibers[0]['fiber_type'] != 'SKY'
    np.testing.assert_array_equal(fibers[0]['wave'], wave)
    np.testing.assert_array_equal(fibers[0]['flux'], flux)
    np.testing.assert_array_equal(fibers[0]['ivar'], ivar)
    # Spatial coords picked up from layout.
    xs = sorted(f['x'] for f in fibers)
    assert xs == [0.0, 1.0]


def test_load_fibers_boxcar_prefix():
    wave = np.linspace(4000.0, 5000.0, 11)
    flux = np.full_like(wave, 7.0)
    ivar = np.ones_like(wave)
    sobjs = _FakeSpecObjs([_make_spec('DET01', 1, 100.0, wave, flux, ivar)])
    spec = _FakeSpectrograph(fiber_types_by_det={1: ['SCI']})
    targetx = np.zeros(200)
    targety = np.zeros(200)

    fibers = load_fibers(sobjs, spec, targetx, targety, prefix='BOX')
    np.testing.assert_array_equal(fibers[0]['flux'], flux)


def test_load_fibers_no_science_raises():
    wave = np.linspace(4000.0, 5000.0, 11)
    sobjs = _FakeSpecObjs([_make_spec('DET01', 1, 100.0, wave,
                                       np.ones(11), np.ones(11))])
    spec = _FakeSpectrograph(fiber_types_by_det={1: ['SKY']})
    targetx = np.zeros(200)
    targety = np.zeros(200)

    with pytest.raises(PypeItError):
        load_fibers(sobjs, spec, targetx, targety, prefix='OPT')


def test_load_fibers_drops_unmatched_science_fibers():
    """Science fibers with layout == -1 (e.g. dead/unmatched) are dropped."""

    class _UnmatchedSpec(_FakeSpectrograph):
        def get_science_fiber_layout_indices(self, det, fiber_ids,
                                              fiber_types):
            # Mark every fiber as unmatched (layout = -1).
            return np.full_like(fiber_ids, -1)

    wave = np.linspace(4000.0, 5000.0, 11)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave)
    sobjs = _FakeSpecObjs([_make_spec('DET01', 1, 100.0, wave, flux, ivar)])
    spec = _UnmatchedSpec(fiber_types_by_det={1: ['SCI']})
    targetx = np.zeros(200)
    targety = np.zeros(200)

    with pytest.raises(PypeItError):
        load_fibers(sobjs, spec, targetx, targety, prefix='OPT')


def test_write_onespec_round_trip(tmp_path):
    wave = np.linspace(4000.0, 5000.0, 51)
    flux = np.ones_like(wave) * 3.0
    ivar = np.ones_like(wave) * 4.0
    hdr = Header()
    hdr['RA'] = 180.0
    hdr['DEC'] = 30.0
    hdr['OBJECT'] = 'TEST_OBJ'
    out = tmp_path / 'extract1d.fits'
    write_onespec(wave, flux, ivar, hdr, 'mmt_binospec_ifu', str(out))
    assert out.exists()
    # Round-trip via OneSpec.from_file
    one = OneSpec.from_file(str(out))
    np.testing.assert_allclose(one.wave, wave)
    np.testing.assert_allclose(one.flux, flux)
    np.testing.assert_allclose(one.ivar, ivar)
    assert one.PYP_SPEC == 'mmt_binospec_ifu'
    # wave_grid_mid is not populated for a non-coadd write
    assert one.wave_grid_mid is None
    # Raw header keys propagate through to the saved primary HDU
    assert one.head0['RA'] == pytest.approx(180.0)
    assert one.head0['DEC'] == pytest.approx(30.0)
    assert one.head0['OBJECT'] == 'TEST_OBJ'


def test_parser_defaults():
    parser = BinospecIFUExtract.get_parser()
    args = parser.parse_args(['spec1d_xx.fits'])
    assert args.spec1d_file == 'spec1d_xx.fits'
    assert args.output is None
    assert args.boxcar is False


def test_parser_boxcar_and_output():
    parser = BinospecIFUExtract.get_parser()
    args = parser.parse_args(['spec1d_xx.fits', '-o', 'out.fits',
                              '--boxcar'])
    assert args.spec1d_file == 'spec1d_xx.fits'
    assert args.output == 'out.fits'
    assert args.boxcar is True
