"""
Unit tests for the MMT/Binospec spectrographs (MOS/longslit and IFU).

Consolidates the Binospec unit tests into a single module, grouped by topic:

- Setup grouping and frame typing (issue #2137)
- Static bad pixel masks
- Overscan subtraction and nonlinearity correction
- Automatic combination-group assignment on IFU data
- ``pypeit_binospec_ifu_extract`` helpers
"""
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.io.fits import Header
from astropy.table import Table

from pypeit import dataPaths
from pypeit import PypeItError
from pypeit.onespec import OneSpec
from pypeit.core.procimg import clean_overscan_vector
from pypeit.spectrographs.mmt_binospec import (binospec_read_amp,
                                               MMTBINOSPECSpectrograph,
                                               MMTBINOSPECIFUSpectrograph)
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts.binospec_ifu_extract import BinospecIFUExtract
from pypeit.core.datacube import (
    load_fibers,
    project_to_sky,
    resample_and_combine,
    write_onespec,
)
from pypeit.gui.binospec_ifu_extract import (
    _SKY_LINE_MASK_HALFWIDTH,
    _SKY_LINE_WAVELENGTHS,
    compute_fiber_fluxes,
    sky_line_mask,
)


# ---------------------------------------------------------------------------
# Setup grouping (configuration_keys) and frame typing (issue #2137)
# ---------------------------------------------------------------------------

def test_configuration_keys_includes_decker():
    """Setups must be split by slit mask (decker), not just grating."""
    spec = MMTBINOSPECSpectrograph()
    keys = spec.configuration_keys()
    assert 'dispname' in keys, \
        f"dispname (grating) must be a configuration key; got {keys}"
    assert 'decker' in keys, \
        "decker (MASK) must be a configuration key so different slit masks " \
        "get separate setups"


def _frame_table(**columns):
    """Build a minimal metadata table for frame-typing tests.

    Every column is broadcast to the same length.
    """
    n = max(len(np.atleast_1d(v)) for v in columns.values())
    data = {k: np.broadcast_to(np.atleast_1d(v), n).copy()
            for k, v in columns.items()}
    return Table(data)


def test_incan_on_is_flat():
    """INCAN=on with screen deployed and arc lamp off is a flat/trace."""
    spec = MMTBINOSPECSpectrograph()
    fitstbl = _frame_table(exptime=10.0, lampstat01='off',
                           lampstat02='deployed', lampstat03='on')
    for ftype in ['pixelflat', 'trace', 'illumflat']:
        assert spec.check_frame_type(ftype, fitstbl)[0], \
            f"INCAN=on frame should be typed as {ftype}"


def test_incan_off_is_not_flat():
    """INCAN=off must NOT be typed as a flat/trace even if screen deployed."""
    spec = MMTBINOSPECSpectrograph()
    fitstbl = _frame_table(exptime=10.0, lampstat01='off',
                           lampstat02='deployed', lampstat03='off')
    for ftype in ['pixelflat', 'trace', 'illumflat']:
        assert not spec.check_frame_type(ftype, fitstbl)[0], \
            f"INCAN=off frame must not be typed as {ftype}"


# ---------------------------------------------------------------------------
# Static bad pixel masks
# ---------------------------------------------------------------------------

@pytest.mark.parametrize('det', [1, 2])
def test_bpm_file_exists(det):
    """Both static BPM files exist in static_calibs."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    assert bpm_path.exists(), f'BPM file not found: {bpm_path}'


@pytest.mark.parametrize('det', [1, 2])
def test_bpm_shape(det):
    """Each BPM file has the correct NumPy shape and dtype."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    bpm = fits.getdata(bpm_path)
    assert bpm.shape == (4096, 4112)
    assert bpm.dtype == np.int8


@pytest.mark.parametrize('det,min_bad,max_bad', [
    (1, 30000, 80000),
    (2, 15000, 50000),
])
def test_bpm_nonzero_count(det, min_bad, max_bad):
    """Each BPM has a reasonable number of masked pixels."""
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    bpm = fits.getdata(bpm_path)
    n_bad = np.sum(bpm > 0)
    assert min_bad < n_bad < max_bad, \
        f'det={det}: {n_bad} bad pixels outside expected range [{min_bad}, {max_bad}]'


@pytest.mark.parametrize('det', [1, 2])
def test_bpm_applied(det):
    """bpm() returns a mask that includes the static BPM pixels."""
    spec = MMTBINOSPECSpectrograph()
    # Use shape parameter to avoid needing a raw data file
    bpm_img = spec.bpm(filename=None, det=det, shape=(4096, 4112))

    # Load the static file directly for comparison
    bpm_path = dataPaths.static_calibs.get_file_path(
        f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
    static_bpm = fits.getdata(bpm_path)

    # Every pixel marked in the static file should be marked in bpm_img
    assert np.all(bpm_img[static_bpm > 0] > 0)


# ---------------------------------------------------------------------------
# Overscan subtraction and nonlinearity correction
# ---------------------------------------------------------------------------

def test_clean_no_outliers():
    """A smooth vector should be returned unchanged."""
    rng = np.random.default_rng(42)
    vec = 1000.0 + rng.normal(scale=1.0, size=100)
    cleaned = clean_overscan_vector(vec)
    np.testing.assert_array_equal(cleaned, vec)


def test_clean_single_outlier():
    """A single large outlier should be interpolated over."""
    vec = np.full(100, 1000.0)
    vec[50] = 2000.0  # outlier: deviation = 1000 >> nsig*rdnoise = 4.0
    cleaned = clean_overscan_vector(vec)
    # Outlier should be replaced with interpolated value (~1000)
    assert abs(cleaned[50] - 1000.0) < 1.0
    # Non-outlier pixels should be unchanged
    assert cleaned[0] == 1000.0
    assert cleaned[99] == 1000.0


def test_clean_edge_outlier():
    """Outliers at edges should be handled (extrapolation)."""
    vec = np.full(50, 500.0)
    vec[0] = 1500.0  # outlier at left edge
    vec[49] = 1500.0  # outlier at right edge
    cleaned = clean_overscan_vector(vec)
    assert abs(cleaned[0] - 500.0) < 1.0
    assert abs(cleaned[49] - 500.0) < 1.0


def test_clean_constant_vector():
    """A constant vector with strict threshold should be unchanged."""
    vec = np.full(20, 500.0)
    # Even with nsig=0, constant vector has zero deviation
    # from its median-filtered version, so nothing is rejected
    cleaned = clean_overscan_vector(vec, nsig=0.0)
    np.testing.assert_array_equal(cleaned, vec)


def test_clean_custom_params():
    """Custom window and nsig should be respected."""
    vec = np.full(20, 100.0)
    vec[10] = 200.0
    # With nsig=100, threshold = 100*4 = 400, so 200-100=100 < 400
    # Outlier should NOT be cleaned
    cleaned = clean_overscan_vector(vec, w=5, nsig=100.0)
    assert cleaned[10] == 200.0


def _make_fake_amp_hdu(bias_level=1000.0, signal=500.0):
    """Create a fake Binospec amplifier HDU with known overscan.

    Layout matches real data: NAXIS1=2114, NAXIS2=2072,
    DATASEC=[51:2098,1:2056].  The image is stored in standard
    FITS orientation (not transposed).

    The data section is filled with bias_level + signal.
    The overscan regions contain only bias_level.
    """
    nx, ny = 2114, 2072
    img = np.full((ny, nx), bias_level, dtype=np.float32)
    # Data section: cols 50:2098, rows 0:2056 (0-indexed)
    img[0:2056, 50:2098] += signal

    hdr = fits.Header()
    hdr['DATASEC'] = '[51:2098,1:2056]'
    hdr['DETSEC'] = '[1:2048,1:2056]'
    hdr['NAXIS1'] = nx
    hdr['NAXIS2'] = ny

    hdu_primary = fits.PrimaryHDU()
    hdu_ext = fits.ImageHDU(data=img, header=hdr)
    hdulist = fits.HDUList([hdu_primary, hdu_ext])
    return hdulist


def test_bias_subtracted():
    """Overscan subtraction should remove the bias level."""
    bias = 1000.0
    signal = 500.0
    hdulist = _make_fake_amp_hdu(bias_level=bias, signal=signal)
    data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)

    # After overscan subtraction, data section should be close to
    # signal only (bias removed)
    med_data = np.median(data)
    assert abs(med_data - signal) < 5.0, \
        f"Bias not removed: median={med_data}, expected ~{signal}"


def test_output_shape():
    """Output data should have datasec dimensions (2048 x 2056)."""
    hdulist = _make_fake_amp_hdu()
    data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)
    # Note: binospec_read_amp transposes the image, so shape is
    # (x, y) = (2048, 2056) after cropping to datasec
    assert data.shape == (2048, 2056), f"Unexpected shape: {data.shape}"


def test_zero_fake_overscan():
    """Returned overscan should be all zeros (fake)."""
    hdulist = _make_fake_amp_hdu()
    data, overscan, datasec, biassec = binospec_read_amp(hdulist, 1)
    assert np.all(overscan == 0), "Overscan should be fake zeros"


def test_row_dependent_bias_removed():
    """Row-dependent bias structure should be removed by overscan."""
    nx, ny = 2114, 2072
    # Create a bias pattern that varies along FITS rows (axis 0)
    row_bias = np.linspace(990, 1010, ny).astype(np.float32)
    img = np.broadcast_to(row_bias[:, None], (ny, nx)).copy()
    # Add signal to data section
    img[0:2056, 50:2098] += 500.0

    hdr = fits.Header()
    hdr['DATASEC'] = '[51:2098,1:2056]'
    hdr['DETSEC'] = '[1:2048,1:2056]'
    hdr['NAXIS1'] = nx
    hdr['NAXIS2'] = ny
    hdu_primary = fits.PrimaryHDU()
    hdu_ext = fits.ImageHDU(data=img, header=hdr)
    hdulist = fits.HDUList([hdu_primary, hdu_ext])

    data, _, _, _ = binospec_read_amp(hdulist, 1)

    # After overscan subtraction, the row-dependent bias pattern
    # should be mostly removed. Check that the row-wise variation
    # in the data is much less than the original 20 ADU range.
    col_medians = np.median(data, axis=0)  # median along x for each y
    assert np.ptp(col_medians) < 5.0, \
        f"Row-dependent bias not removed: range={np.ptp(col_medians)}"


def test_coefficients_shape():
    """Nonlinearity coefficients should be 8x5 (8 amps, degree 4)."""
    coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
    assert coeffs.shape == (8, 5)


def test_coefficients_zero_constant():
    """All constant terms should be zero."""
    coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
    np.testing.assert_array_equal(coeffs[:, 0], 0.0)


def test_coefficients_near_unity_linear():
    """Linear terms should be close to 1.0 (small correction)."""
    coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
    assert np.all(np.abs(coeffs[:, 1] - 1.0) < 0.01)


def test_correction_applied_in_read_amp():
    """binospec_read_amp should apply nonlinearity correction."""
    nx, ny = 2114, 2072
    bias_level = 1000.0
    signal = 10000.0
    img = np.full((ny, nx), bias_level, dtype=np.float32)
    # Data section gets bias + signal; overscan has only bias
    img[0:2056, 50:2098] += signal

    hdr = fits.Header()
    hdr['DATASEC'] = '[51:2098,1:2056]'
    hdr['DETSEC'] = '[1:2048,1:2056]'
    hdr['NAXIS1'] = nx
    hdr['NAXIS2'] = ny
    hdu_primary = fits.PrimaryHDU()
    hdu_ext = fits.ImageHDU(data=img, header=hdr)
    hdulist = fits.HDUList([hdu_primary, hdu_ext])

    data, _, _, _ = binospec_read_amp(hdulist, 1)

    # After overscan subtraction, data ~ signal. Nonlinearity
    # correction then maps signal -> polyval(signal, coeffs).
    coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs[0]
    expected = np.polynomial.polynomial.polyval(signal, coeffs)
    med_data = np.median(data)
    assert abs(med_data - expected) < 5.0, \
        f"Nonlinearity not applied: got {med_data}, expected ~{expected}"


def test_correction_is_small():
    """At typical science levels (~1000 ADU), correction < 1%."""
    coeffs = MMTBINOSPECSpectrograph.nonlinearity_coeffs
    test_counts = 1000.0
    for i in range(8):
        corrected = np.polynomial.polynomial.polyval(test_counts, coeffs[i])
        ratio = corrected / test_counts
        assert 0.99 < ratio < 1.01, \
            f"Amp {i+1}: correction too large: {ratio}"


# ---------------------------------------------------------------------------
# Automatic combination-group assignment on IFU data
#
# The Binospec IFU has no nod partner (sky comes from dedicated sky fibers via
# a joint fit), so get_comb_group never sets bkg_id; it only groups
# science/standard frames that share a pointing onto a common comb_id so
# same-pointing frames coadd while distinct dither positions each reduce to
# their own spec1d/spec2d.
# ---------------------------------------------------------------------------

def _comb_table(rows, setups=None):
    """
    Build a minimal metadata table shaped like the one
    :func:`~pypeit.metadata.PypeItMetaData.set_combination_groups` hands to
    ``get_comb_group``: every science/standard frame already carries a unique
    ``comb_id`` (1..k in table order) and ``bkg_id`` is -1 everywhere.  The
    ``dithoff``/``dithpos`` columns are seeded with their card defaults, as
    ``PypeItMetaData._build`` does from the ``init_meta`` definitions.

    Args:
        rows (list): sequence of ``(frametype, ra_deg, dec_deg, mjd)`` tuples.
        setups (list, optional): per-row setup label; defaults to all ``'A'``.
    """
    n = len(rows)
    ftype = [r[0] for r in rows]
    comb = np.full(n, -1, dtype=int)
    sci = [i for i, f in enumerate(ftype) if 'science' in f or 'standard' in f]
    for k, i in enumerate(sci):
        comb[i] = k + 1
    t = Table()
    t['filename'] = [f'f{i}.fits' for i in range(n)]
    t['frametype'] = ftype
    t['setup'] = setups if setups is not None else ['A'] * n
    t['ra'] = [r[1] for r in rows]
    t['dec'] = [r[2] for r in rows]
    t['mjd'] = [r[3] for r in rows]
    t['comb_id'] = comb
    t['bkg_id'] = np.full(n, -1, dtype=int)
    t['dithoff'] = np.zeros(n, dtype=float)
    t['dithpos'] = np.full(n, 'None')
    return t


def _offset(ra0, dec0, east_as, north_as):
    """Offset a pointing by (east, north) arcsec, returning (ra_deg, dec_deg)."""
    ra = ra0 + (east_as / 3600.0) / np.cos(np.radians(dec0))
    dec = dec0 + north_as / 3600.0
    return ra, dec


# JADES field pointing used throughout (12:36:52.6 +62:07:56.6).
RA0, DEC0 = 189.21916666666664, 62.13238888888889


def test_same_pointing_shares_comb_id():
    # Eight science frames at one identical pointing must coadd: a single
    # shared comb_id, no background pairing, ~zero dither offset.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(8)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb == comb[0]), \
        'all same-pointing science frames must share one comb_id so they coadd'
    assert comb[0] > 0, 'the shared comb_id must be a real (positive) group id'
    assert np.all(np.asarray(t['bkg_id']) == -1), \
        'bkg_id must stay unset: IFU sky comes from sky fibers, not a nod partner'
    assert np.allclose(np.asarray(t['dithoff'], dtype=float), 0.0, atol=1e-6), \
        'a single pointing has zero dither offset'
    assert len(set(str(x) for x in t['dithpos'])) == 1, \
        'a single pointing yields exactly one dithpos label'


def test_two_position_dither_splits_comb_id():
    # A two-position dither (0.5" apart, well above the 0.1" tolerance) taken
    # ABAB: frames at the same position share a comb_id, the two positions get
    # different comb_ids, and nothing is background-paired.
    spec = load_spectrograph('mmt_binospec_ifu')
    posA = (RA0, DEC0)
    posB = _offset(RA0, DEC0, 0.5, 0.0)
    seq = [posA, posB, posA, posB]
    rows = [('science', p[0], p[1], 60819.14 + 0.01 * i)
            for i, p in enumerate(seq)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[2], 'frames at position A must share a comb_id'
    assert comb[1] == comb[3], 'frames at position B must share a comb_id'
    assert comb[0] != comb[1], \
        'the two dither positions must get different comb_ids'
    assert len(set(comb)) == 2, 'a two-position dither yields exactly two comb_ids'
    assert np.all(np.asarray(t['bkg_id']) == -1), \
        'bkg_id must stay unset for the IFU'
    assert len(set(str(x) for x in t['dithpos'])) == 2, \
        'two distinct positions yield two dithpos labels'


def test_subtolerance_jitter_stays_one_group():
    # Pointing wander below the tolerance (0.02") must not split the group.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = []
    for i in range(4):
        ra, dec = _offset(RA0, DEC0, 0.02 * (i % 2), 0.0)
        rows.append(('science', ra, dec, 60819.14 + 0.01 * i))
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb == comb[0]), \
        'pointing wander below the tolerance must not split the group'


def test_large_offset_gets_own_comb_id():
    # A large sky-acquisition offset (60") lands in its own cluster: the six
    # on-object frames share one comb_id, the offset frame gets its own. No
    # special sky handling -- it is just another cluster.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(6)]
    sky = _offset(RA0, DEC0, 60.0, 0.0)
    rows.append(('science', sky[0], sky[1], 60819.14 + 0.07))
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert np.all(comb[:6] == comb[0]), \
        'the six on-object frames must share one comb_id'
    assert comb[6] != comb[0], \
        'a large sky-acquisition offset must get its own comb_id'
    assert len(set(comb)) == 2, 'on-object + offset gives exactly two comb_ids'


def test_two_setups_not_merged():
    # Identical pointing but two different configurations must never share a
    # comb_id -- combination cannot cross a setup.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('science', RA0, DEC0, 60819.14 + 0.01 * i) for i in range(4)]
    setups = ['A', 'A', 'B', 'B']
    t = spec.get_comb_group(_comb_table(rows, setups=setups))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[1], 'same-setup same-pointing frames share a comb_id'
    assert comb[2] == comb[3], 'same-setup same-pointing frames share a comb_id'
    assert comb[0] != comb[2], \
        'identical pointing in different setups must never share a comb_id'


def test_standard_at_own_pointing_separate():
    # A standard star observed at a different pointing than the science field
    # gets its own comb_id; the science frames still coadd together.
    spec = load_spectrograph('mmt_binospec_ifu')
    std = _offset(RA0, DEC0, 200.0, 120.0)
    rows = [('science', RA0, DEC0, 60819.14),
            ('science', RA0, DEC0, 60819.15),
            ('science', RA0, DEC0, 60819.16),
            ('standard', std[0], std[1], 60819.20)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == comb[1] == comb[2], 'the three science frames must coadd'
    assert comb[3] != comb[0], \
        'a standard at a different pointing gets its own comb_id'
    assert np.all(np.asarray(t['bkg_id']) == -1), \
        'bkg_id must stay unset for the IFU'


def test_calibration_frames_untouched():
    # Non science/standard frames keep comb_id/bkg_id == -1.
    spec = load_spectrograph('mmt_binospec_ifu')
    rows = [('arc,tilt', RA0, DEC0, 60819.10),
            ('pixelflat,illumflat,trace', RA0, DEC0, 60819.11),
            ('science', RA0, DEC0, 60819.14),
            ('science', RA0, DEC0, 60819.15)]
    t = spec.get_comb_group(_comb_table(rows))
    comb = np.asarray(t['comb_id'])
    assert comb[0] == -1 and comb[1] == -1, \
        'non science/standard frames must keep comb_id == -1'
    assert comb[2] == comb[3] and comb[2] > 0, \
        'the two science frames still share a real comb_id'


# ---------------------------------------------------------------------------
# pypeit_binospec_ifu_extract helpers
# ---------------------------------------------------------------------------

# Real IFU spectrograph supplies the shared TAN/POSANG WCS convention used by
# project_to_sky.
_IFU_SPEC = MMTBINOSPECIFUSpectrograph()


def test_sky_line_mask_marks_each_line():
    wave = np.array(_SKY_LINE_WAVELENGTHS, dtype=float)
    mask = sky_line_mask(wave)
    assert mask.shape == wave.shape, \
        f"mask shape {mask.shape} should match wave shape {wave.shape}"
    assert np.all(mask), 'Each sky-line wavelength itself should be masked'


def test_sky_line_mask_excludes_far_wavelengths():
    # 5000 Angstrom is well clear of every line in the table
    wave = np.array([5000.0, 5577.34, 5000.0], dtype=float)
    mask = sky_line_mask(wave)
    assert mask.tolist() == [False, True, False], \
        f"only the 5577.34 sky line should be masked, got {mask.tolist()}"


def test_sky_line_mask_halfwidth_boundary():
    line = _SKY_LINE_WAVELENGTHS[0]
    # Just inside ± halfwidth → masked; just outside → not masked.
    inside = line + _SKY_LINE_MASK_HALFWIDTH - 0.01
    outside = line + _SKY_LINE_MASK_HALFWIDTH + 0.01
    mask = sky_line_mask(np.array([inside, outside]))
    assert mask.tolist() == [True, False], \
        "point just inside the halfwidth should be masked and just outside " \
        f"not, got {mask.tolist()}"


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
    assert ra.shape == (1,), f"ra should have shape (1,), got {ra.shape}"
    assert dec.shape == (1,), f"dec should have shape (1,), got {dec.shape}"
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
    assert out.shape == (2,), \
        f"expected one summed flux per fiber, shape (2,), got {out.shape}"
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
    assert out_w[0] >= 4050.0, \
        f"output grid should start at the overlap (>= 4050), got {out_w[0]:.2f}"
    assert out_w[-1] <= 4100.0, \
        f"output grid should end at the overlap (<= 4100), got {out_w[-1]:.2f}"
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
    assert len(fibers) == 2, \
        f"two science fibers should survive (sky dropped), got {len(fibers)}"
    assert fibers[0]['fiber_type'] != 'SKY', \
        f"surviving fibers must be science, got {fibers[0]['fiber_type']!r}"
    np.testing.assert_array_equal(fibers[0]['wave'], wave)
    np.testing.assert_array_equal(fibers[0]['flux'], flux)
    np.testing.assert_array_equal(fibers[0]['ivar'], ivar)
    # Spatial coords picked up from layout.
    xs = sorted(f['x'] for f in fibers)
    assert xs == [0.0, 1.0], \
        f"spatial x coords should be picked up from the layout, got {xs}"


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
    assert out.exists(), f"write_onespec should create the output file {out}"
    # Round-trip via OneSpec.from_file
    one = OneSpec.from_file(str(out))
    np.testing.assert_allclose(one.wave, wave)
    np.testing.assert_allclose(one.flux, flux)
    np.testing.assert_allclose(one.ivar, ivar)
    assert one.PYP_SPEC == 'mmt_binospec_ifu', \
        f"PYP_SPEC should round-trip, got {one.PYP_SPEC!r}"
    # wave_grid_mid is not populated for a non-coadd write
    assert one.wave_grid_mid is None, \
        f"wave_grid_mid should be None for a non-coadd write, got {one.wave_grid_mid}"
    # Raw header keys propagate through to the saved primary HDU
    assert one.head0['RA'] == pytest.approx(180.0), \
        f"RA header should propagate to the primary HDU, got {one.head0['RA']}"
    assert one.head0['DEC'] == pytest.approx(30.0), \
        f"DEC header should propagate to the primary HDU, got {one.head0['DEC']}"
    assert one.head0['OBJECT'] == 'TEST_OBJ', \
        f"OBJECT header should propagate to the primary HDU, got {one.head0['OBJECT']!r}"


def test_parser_defaults():
    parser = BinospecIFUExtract.get_parser()
    args = parser.parse_args(['spec1d_xx.fits'])
    assert args.spec1d_file == 'spec1d_xx.fits', \
        f"positional spec1d_file misparsed: {args.spec1d_file!r}"
    assert args.output is None, \
        f"output should default to None, got {args.output!r}"
    assert args.boxcar is False, \
        f"boxcar should default to False, got {args.boxcar}"


def test_parser_boxcar_and_output():
    parser = BinospecIFUExtract.get_parser()
    args = parser.parse_args(['spec1d_xx.fits', '-o', 'out.fits',
                              '--boxcar'])
    assert args.spec1d_file == 'spec1d_xx.fits', \
        f"positional spec1d_file misparsed: {args.spec1d_file!r}"
    assert args.output == 'out.fits', \
        f"-o should set output to 'out.fits', got {args.output!r}"
    assert args.boxcar is True, \
        f"--boxcar should set boxcar True, got {args.boxcar}"
