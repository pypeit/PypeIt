"""Tests for fiber block-slit extraction logic."""
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.io import fits
from astropy.table import Table

from pypeit.find_objects import FiberFindObjects
from pypeit.slittrace import SlitTraceSet
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph


spec = MMTBINOSPECIFUSpectrograph()


def test_get_fiber_blocks_det01():
    """DET01 should have 21 blocks: 5 sky (8 fibers) + 16 science (20)."""
    blocks = spec.get_fiber_blocks(1)
    assert len(blocks) == 21, f"DET01 should have 21 blocks, got {len(blocks)}"
    # Sky blocks: 1, 6, 11, 16, 21
    for b in [0, 5, 10, 15, 20]:  # 0-indexed
        assert blocks[b]['type'] == 'sky', \
            f"block {b} should be a sky block, got {blocks[b]['type']}"
        assert blocks[b]['nfibers'] == 8, \
            f"sky block {b} should have 8 fibers, got {blocks[b]['nfibers']}"
    # Science blocks
    for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
        assert blocks[b]['type'] == 'science', \
            f"block {b} should be a science block, got {blocks[b]['type']}"
        assert blocks[b]['nfibers'] == 20, \
            f"science block {b} should have 20 fibers, got {blocks[b]['nfibers']}"


def test_get_fiber_blocks_det02():
    """DET02 should have 21 live blocks (dead fibers in block -1 are excluded).

    DET02 has 4 dead fibers assigned to block_id -1 in the reference profile.
    These are excluded from the returned block list, leaving 21 standard blocks.
    Some science blocks have 19 fibers instead of 20 because of those dead fibers.
    """
    blocks = spec.get_fiber_blocks(2)
    assert len(blocks) == 21, \
        f"DET02 should have 21 live blocks, got {len(blocks)}"
    # Sky blocks are at the same positions as DET01
    for b in [0, 5, 10, 15, 20]:  # 0-indexed
        assert blocks[b]['type'] == 'sky', \
            f"block {b} should be a sky block, got {blocks[b]['type']}"
        assert blocks[b]['nfibers'] == 8, \
            f"sky block {b} should have 8 fibers, got {blocks[b]['nfibers']}"
    # Science blocks have 19 or 20 fibers (some have dead fibers replaced)
    for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
        assert blocks[b]['type'] == 'science', \
            f"block {b} should be a science block, got {blocks[b]['type']}"
        assert blocks[b]['nfibers'] in (19, 20), \
            f"DET02 science block {b} should have 19 or 20 fibers, got " \
            f"{blocks[b]['nfibers']}"


def test_det02_reference_is_mirrored_onto_pypeit_frame():
    """DET02 reference rows/blocks should increase with PypeIt pixel."""
    with fits.open(spec._ifu_calib_path() / 'fiber_ref_profile.fits') as hdu:
        raw = hdu[2].data.copy()
    ref = spec.load_fiber_ref_profile(2)
    live = ref['FIB_BLOCK'] >= 0

    assert np.all(np.diff(ref['TR_PIX'][live]) > 0), \
        "DET02 live TR_PIX should increase with row index"
    assert np.all(np.diff(ref['FIB_BLOCK'][live]) >= 0), \
        "DET02 FIB_BLOCK should increase with TR_PIX"
    np.testing.assert_array_equal(
        ref['FIB_ID'], raw['FIB_ID'][::-1],
        err_msg="DET02 rows should be the stored IDL rows reversed")
    raw_live = raw['FIB_BLOCK'] >= 0
    np.testing.assert_allclose(
        ref['TR_PIX'][live],
        spec.ref_mirror_pix_b - raw['TR_PIX'][raw_live][::-1],
        err_msg="DET02 TR_PIX should be mirrored about ref_mirror_pix_b")
    np.testing.assert_allclose(
        ref['TR_H3'], -raw['TR_H3'][::-1],
        err_msg="Odd Gauss-Hermite terms should change sign under the mirror")


def test_det02_dead_fibers_match_idl():
    """Issue #2220: DET02 dead fibers must be the IDL ones (B27/B79/B159/B199)."""
    ref = spec.load_fiber_ref_profile(2)
    names = np.char.strip(ref['FIB_NAME'])
    dead = set(names[ref['FIB_DEAD_FLAG'] == 1])
    assert dead == {'B27_DEAD', 'B79_DEAD', 'B159_DEAD', 'B199_DEAD'}, \
        f"Unexpected DET02 dead fibers: {sorted(dead)}"
    for name in ('B122', 'B162', 'B242', 'B294'):
        assert name in names[ref['FIB_BLOCK'] >= 0], \
            f"{name} should be a live DET02 fiber"


def test_det02_reference_profile_is_mirrored():
    """DET02 PROF_REF should be the stored profile mirrored onto PypeIt pixels."""
    with fits.open(spec._ifu_calib_path() / 'fiber_ref_profile.fits') as hdu:
        raw = hdu[3].data['PROF_REF'].astype(float)
    np.testing.assert_allclose(spec._load_ref_spatial_profile(1), raw[0],
                               err_msg="DET01 PROF_REF should be unchanged")
    prof = spec._load_ref_spatial_profile(2)
    x = np.arange(raw.shape[1], dtype=float)
    pix = np.array([500, 2000, 3500])
    np.testing.assert_allclose(prof[pix],
                               np.interp(spec.ref_mirror_pix_b - pix, x, raw[1]),
                               err_msg="DET02 PROF_REF should be mirrored")


def test_det02_live_reference_matches_mirrored_fibers():
    """Matching the mirrored DET02 positions should recover every live fiber by name."""
    ref = spec.load_fiber_ref_profile(2)
    live = ref['FIB_BLOCK'] >= 0
    meta = spec.get_fiber_metadata(2, np.round(ref['TR_PIX'][live]).astype(int),
                                   ref['TR_PIX'][live].astype(float))
    np.testing.assert_array_equal(meta['fiber_id'], ref['FIB_ID'][live],
                                  err_msg="Every live DET02 fiber should match itself")


@pytest.mark.parametrize('det', [1, 2])
def test_science_fiber_layout_indices_match_by_name(det):
    """Each science fiber should map to the layout entry with its own name."""
    ref = spec.load_fiber_ref_profile(det)
    names = np.char.strip(ref['FIB_NAME'])
    types = np.where(np.char.startswith(names, 'SKY'), 'SKY', 'SCI')
    idx = spec.get_science_fiber_layout_indices(det, ref['FIB_ID'], types)
    with fits.open(spec._ifu_calib_path() / 'bino_IFU_sky_layout.fits') as hdu:
        layout_names = np.char.strip(hdu[1].data[0]['TARGET_NAME'])
    mapped = idx >= 0
    assert mapped.sum() == (320 if det == 1 else 316), \
        f"DET{det:02d}: unexpected number of mapped science fibers {mapped.sum()}"
    np.testing.assert_array_equal(layout_names[idx[mapped]], names[mapped],
                                  err_msg="Layout entry name should match fiber name")
    assert all('_DEAD' in n for n in names[~mapped & (types == 'SCI')]), \
        "Only dead science fibers should be unmapped"


def test_block_fiber_positions():
    """Fiber positions within blocks should be monotonically ordered."""
    blocks = spec.get_fiber_blocks(1)
    for block in blocks:
        positions = block['fiber_positions']
        assert np.all(np.diff(positions) > 0), \
            f"Block {block['block_id']} positions not monotonic"


def test_block_gaps():
    """Inter-block gaps should be > 50 pixels."""
    blocks = spec.get_fiber_blocks(1)
    for i in range(len(blocks) - 1):
        gap = blocks[i+1]['min_pix'] - blocks[i]['max_pix']
        assert gap > 50, f"Gap between blocks {i} and {i+1} is only {gap:.1f} px"


def test_fiber_position_shift_from_slits():
    """The slit-derived fiber shift should recover the observed bulk offset."""
    blocks = spec.get_fiber_blocks(1)
    nspec = 11
    shift = 3.25
    margin = 7.0
    left = np.zeros((nspec, len(blocks)))
    right = np.zeros_like(left)
    for i, block in enumerate(blocks):
        left[:, i] = block['min_pix'] + shift - margin
        right[:, i] = block['max_pix'] + shift + margin

    slits = SlitTraceSet(left, right, 'Fiber', nspat=4096,
                         PYP_SPEC='mmt_binospec_ifu')
    recovered = spec.get_fiber_position_shift(slits, 1)
    assert np.isclose(recovered, shift), \
        f"recovered fiber shift should be {shift}, got {recovered:.4f}"


def test_ifu_config_specific_keeps_standard_tilts():
    """IFU config overrides should keep the intended tilt fit orders."""
    fitstbl_row = Table({'dispname': ['x270'], 'decker': ['IFU']})
    par = spec.config_specific_par(fitstbl_row)

    assert par['calibrations']['tilts']['spat_order'] == 3, \
        "IFU tilt spat_order should be 3, got " \
        f"{par['calibrations']['tilts']['spat_order']}"
    assert par['calibrations']['tilts']['spec_order'] == 5, \
        "IFU tilt spec_order should be 5, got " \
        f"{par['calibrations']['tilts']['spec_order']}"
    assert par['calibrations']['tilts']['maxdev2d'] == 0.5, \
        "IFU tilt maxdev2d should be 0.5, got " \
        f"{par['calibrations']['tilts']['maxdev2d']}"
    assert par['calibrations']['wavelengths']['n_final'] == 5, \
        "IFU wavelengths n_final should be 5, got " \
        f"{par['calibrations']['wavelengths']['n_final']}"
    assert par['reduce']['skysub']['bspline_spacing'] == 1.05, \
        "x270 bspline_spacing should be 1.05, got " \
        f"{par['reduce']['skysub']['bspline_spacing']}"


def test_det02_static_illumination_is_mirror_mapped():
    """DET02 illumination should be reversed to follow the reversed reference rows."""
    illum_file = spec._ifu_calib_path() / 'fiber_illumination.fits'
    with fits.open(illum_file) as hdu:
        raw = hdu[1].data['F_ILLUM'].copy()

    np.testing.assert_allclose(
        spec.load_fiber_illumination(1), raw[0],
        err_msg="DET01 illumination should match the stored side-A row")
    np.testing.assert_allclose(
        spec.load_fiber_illumination(2), raw[1][::-1],
        err_msg="DET02 illumination should be the side-B row reversed onto "
                "PypeIt's reference order")


def test_attach_fiber_throughput_maps_by_id():
    """``fiber_throughput`` should be looked up by MASKDEF_ID."""
    sobjs = [
        SimpleNamespace(MASKDEF_ID=42),
        SimpleNamespace(MASKDEF_ID=99),
        SimpleNamespace(MASKDEF_ID=None),
    ]
    fiber_flatimages = SimpleNamespace(
        fiber_ids=np.array([42, 99]),
        fiber_throughput=np.array([1.05, 1.35]))
    spectrograph = SimpleNamespace()  # no load_fiber_illumination

    FiberFindObjects._attach_fiber_throughput(
        sobjs, fiber_flatimages, spectrograph, det=1)

    assert sobjs[0].fiber_throughput == 1.05, \
        f"MASKDEF_ID=42 should map to throughput 1.05, got {sobjs[0].fiber_throughput}"
    assert sobjs[1].fiber_throughput == 1.35, \
        f"MASKDEF_ID=99 should map to throughput 1.35, got {sobjs[1].fiber_throughput}"
    assert sobjs[2].fiber_throughput == 1.0, \
        f"MASKDEF_ID=None should fall back to unity, got {sobjs[2].fiber_throughput}"


def test_attach_fiber_throughput_falls_back_to_unity():
    """Missing flat / throughput should yield unity scalars."""
    sobjs = [SimpleNamespace(MASKDEF_ID=42)]
    spectrograph = SimpleNamespace()
    FiberFindObjects._attach_fiber_throughput(
        sobjs, None, spectrograph, det=1)
    assert sobjs[0].fiber_throughput == 1.0, \
        "missing flat/throughput should yield a unity scalar, got " \
        f"{sobjs[0].fiber_throughput}"


def test_attach_fiber_throughput_multiplies_in_static_illumination():
    """Static fiber_illumination should multiply the flat-derived scalar."""
    sobjs = [
        SimpleNamespace(MASKDEF_ID=42),
        SimpleNamespace(MASKDEF_ID=99),
    ]
    fiber_flatimages = SimpleNamespace(
        fiber_ids=np.array([42, 99]),
        fiber_throughput=np.array([1.0, 1.30]))
    spectrograph = SimpleNamespace(
        load_fiber_illumination=lambda _det: np.array([0.95, 1.05]),
        load_fiber_ref_profile=lambda _det: {'FIB_ID': np.array([42, 99])})

    FiberFindObjects._attach_fiber_throughput(
        sobjs, fiber_flatimages, spectrograph, det=1)

    np.testing.assert_allclose(
        sobjs[0].fiber_throughput, 1.0 * 0.95,
        err_msg="static illumination (0.95) should multiply the flat scalar (1.0)")
    np.testing.assert_allclose(
        sobjs[1].fiber_throughput, 1.30 * 1.05,
        err_msg="static illumination (1.05) should multiply the flat scalar (1.30)")


def test_identify_fibers_in_block():
    """Given fiber peak positions within a block, identify each fiber."""
    blocks = spec.get_fiber_blocks(1)
    # Use the reference positions as "detected" positions for block 2 (science, 20 fibers)
    block = blocks[1]  # block 2 (0-indexed = 1)
    detected_positions = block['fiber_positions']
    result = spec.identify_fibers_in_block(
        det=1, block_idx=1, detected_positions=detected_positions)
    assert len(result['fiber_id']) == 20, \
        f"block 2 should yield 20 fiber ids, got {len(result['fiber_id'])}"
    assert len(result['fiber_name']) == 20, \
        f"block 2 should yield 20 fiber names, got {len(result['fiber_name'])}"
    # All should match since we used reference positions
    assert np.all(result['fiber_id'] > 0), \
        "all fibers should match when reference positions are used as detections; " \
        f"unmatched: {np.sum(result['fiber_id'] <= 0)}"


def test_identify_fibers_with_offset():
    """Fibers should still match with a small position offset."""
    blocks = spec.get_fiber_blocks(1)
    block = blocks[1]
    # Add 2-pixel offset (simulating flexure)
    detected_positions = block['fiber_positions'] + 2.0
    result = spec.identify_fibers_in_block(
        det=1, block_idx=1, detected_positions=detected_positions)
    assert np.all(result['fiber_id'] > 0), \
        "all fibers should still match with a 2-pixel offset; unmatched: " \
        f"{np.sum(result['fiber_id'] <= 0)}"
