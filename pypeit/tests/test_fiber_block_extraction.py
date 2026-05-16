"""Tests for fiber block-slit extraction logic."""
from types import SimpleNamespace

import numpy as np
from astropy.io import fits
from astropy.table import Table

from pypeit.extraction import FiberExtract
from pypeit.find_objects import FiberFindObjects
from pypeit.slittrace import SlitTraceSet
from pypeit.spectrographs.mmt_binospec import MMTBINOSPECIFUSpectrograph


spec = MMTBINOSPECIFUSpectrograph()


def test_get_fiber_blocks_det01():
    """DET01 should have 21 blocks: 5 sky (8 fibers) + 16 science (20)."""
    blocks = spec.get_fiber_blocks(1)
    assert len(blocks) == 21
    # Sky blocks: 1, 6, 11, 16, 21
    for b in [0, 5, 10, 15, 20]:  # 0-indexed
        assert blocks[b]['type'] == 'sky'
        assert blocks[b]['nfibers'] == 8
    # Science blocks
    for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
        assert blocks[b]['type'] == 'science'
        assert blocks[b]['nfibers'] == 20


def test_get_fiber_blocks_det02():
    """DET02 should have 21 live blocks (dead fibers in block -1 are excluded).

    DET02 has 4 dead fibers assigned to block_id -1 in the reference profile.
    These are excluded from the returned block list, leaving 21 standard blocks.
    Some science blocks have 19 fibers instead of 20 because of those dead fibers.
    """
    blocks = spec.get_fiber_blocks(2)
    assert len(blocks) == 21
    # Sky blocks are at the same positions as DET01
    for b in [0, 5, 10, 15, 20]:  # 0-indexed
        assert blocks[b]['type'] == 'sky'
        assert blocks[b]['nfibers'] == 8
    # Science blocks have 19 or 20 fibers (some have dead fibers replaced)
    for b in [1, 2, 3, 4, 6, 7, 8, 9, 11, 12, 13, 14, 16, 17, 18, 19]:
        assert blocks[b]['type'] == 'science'
        assert blocks[b]['nfibers'] in (19, 20)


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
    assert np.isclose(spec.get_fiber_position_shift(slits, 1), shift)


def test_ifu_config_specific_keeps_standard_tilts():
    """IFU config overrides should keep the intended tilt fit orders."""
    fitstbl_row = Table({'dispname': ['x270'], 'decker': ['IFU']})
    par = spec.config_specific_par(fitstbl_row)

    assert par['calibrations']['tilts']['spat_order'] == 3
    assert par['calibrations']['tilts']['spec_order'] == 5
    assert par['reduce']['skysub']['bspline_spacing'] == 1.05


def test_det02_static_illumination_is_mirror_mapped():
    """DET02 illumination should be reversed onto PypeIt's reference order."""
    illum_file = spec._ifu_calib_path() / 'fiber_illumination.fits'
    with fits.open(illum_file) as hdu:
        raw = hdu[1].data['F_ILLUM'].copy()

    np.testing.assert_allclose(spec.load_fiber_illumination(1), raw[0])
    np.testing.assert_allclose(spec.load_fiber_illumination(2), raw[1][::-1])


def test_reconstruct_fiber_skymodel_is_per_pixel():
    """Diagnostic SKYMODEL should be in SCIIMG-like per-pixel units."""
    nspec, nspat = 3, 6
    left = np.zeros((nspec, 1))
    right = np.full((nspec, 1), 5.0)
    slits = SlitTraceSet(left, right, 'Fiber', nspat=nspat,
                         PYP_SPEC='mmt_binospec_ifu')

    findobj = FiberFindObjects.__new__(FiberFindObjects)
    findobj.sciImg = SimpleNamespace(image=np.zeros((nspec, nspat)))
    findobj.slits = slits
    findobj.spat_flexure_shift = None

    sobj = SimpleNamespace(
        BOX_COUNTS_SKY=np.full(nspec, 12.0),
        flat_corr=np.array([1.0, 2.0, 1.0]),
        TRACE_SPAT=np.full(nspec, 2.0),
        BOX_R_PIX=1.0,
        SLITID=slits.spat_id[0])

    sky_2d = findobj._reconstruct_fiber_skymodel([sobj])

    np.testing.assert_allclose(sky_2d[0, 1:4], 4.0)
    np.testing.assert_allclose(sky_2d[1, 1:4], 8.0)
    np.testing.assert_allclose(sky_2d[2, 1:4], 4.0)
    assert np.all(sky_2d[:, [0, 4, 5]] == 0.0)


def test_apply_flat_correction_does_not_double_apply_static_illumination():
    """Extracted normflat should be the default fiber throughput correction."""
    findobj = FiberFindObjects.__new__(FiberFindObjects)
    findobj.sciImg = SimpleNamespace(detector=SimpleNamespace(det=1))
    findobj.spectrograph = SimpleNamespace(
        use_static_fiber_illumination=False,
        load_fiber_illumination=lambda _det: (_ for _ in ()).throw(
            AssertionError('static illumination should not be loaded')))

    sobj = SimpleNamespace(
        BOX_COUNTS=np.array([10.0, 20.0]),
        BOX_WAVE=np.array([5000.0, 5001.0]),
        BOX_COUNTS_IVAR=np.ones(2),
        MASKDEF_ID=42)
    fiber_flatimages = SimpleNamespace(
        normflat=np.array([[0.5, 0.5]]),
        normflat_wave=np.array([5000.0, 5001.0]),
        fiber_ids=np.array([42]))

    findobj._apply_flat_correction([sobj], fiber_flatimages)

    np.testing.assert_allclose(sobj.BOX_COUNTS, [20.0, 40.0])
    np.testing.assert_allclose(sobj.BOX_COUNTS_IVAR, [0.25, 0.25])
    np.testing.assert_allclose(sobj.flat_corr, [0.5, 0.5])


def test_apply_flat_correction_can_apply_static_illumination():
    """Static illumination should multiply the normalized flat divisor."""
    findobj = FiberFindObjects.__new__(FiberFindObjects)
    findobj.sciImg = SimpleNamespace(detector=SimpleNamespace(det=1))
    findobj.spectrograph = SimpleNamespace(
        use_static_fiber_illumination=True,
        load_fiber_illumination=lambda _det: np.array([0.8]),
        load_fiber_ref_profile=lambda _det: {'FIB_ID': np.array([42])})

    sobj = SimpleNamespace(
        BOX_COUNTS=np.array([10.0, 20.0]),
        BOX_WAVE=np.array([5000.0, 5001.0]),
        BOX_COUNTS_IVAR=np.ones(2),
        MASKDEF_ID=42)
    fiber_flatimages = SimpleNamespace(
        normflat=np.array([[0.5, 0.5]]),
        normflat_wave=np.array([5000.0, 5001.0]),
        fiber_ids=np.array([42]))

    findobj._apply_flat_correction([sobj], fiber_flatimages)

    np.testing.assert_allclose(sobj.BOX_COUNTS, [25.0, 50.0])
    np.testing.assert_allclose(sobj.BOX_COUNTS_IVAR, [0.16, 0.16])
    np.testing.assert_allclose(sobj.flat_corr, [0.4, 0.4])


def test_profile_extraction_can_drive_box_sky_fit():
    """Profile-fit spectra should be usable as the 1D sky-fit input."""
    sobj = SimpleNamespace(
        OPT_WAVE=np.array([5000.0, 5001.0]),
        OPT_COUNTS=np.array([11.0, 12.0]),
        OPT_COUNTS_IVAR=np.array([0.5, 0.25]),
        OPT_COUNTS_SIG=np.array([1.4, 2.0]),
        OPT_COUNTS_NIVAR=np.array([0.4, 0.2]),
        OPT_MASK=np.array([True, False]),
        OPT_FWHM=np.array([3.0, 3.1]),
        OPT_FLAT=np.array([100.0, 101.0]),
        OPT_COUNTS_SKY=np.array([1.0, 1.5]),
        OPT_COUNTS_SIG_DET=np.array([0.1, 0.2]))

    FiberFindObjects._copy_optimal_to_boxcar(sobj)

    np.testing.assert_allclose(sobj.BOX_WAVE, sobj.OPT_WAVE)
    np.testing.assert_allclose(sobj.BOX_COUNTS, sobj.OPT_COUNTS)
    np.testing.assert_allclose(sobj.BOX_COUNTS_SKY, sobj.OPT_COUNTS_SKY)
    assert np.array_equal(sobj.BOX_MASK, sobj.OPT_MASK)
    assert sobj.BOX_COUNTS is not sobj.OPT_COUNTS


def test_profiled_fiber_skymodel_uses_spatial_profile():
    """Profiled SKYMODEL should distribute 1D sky using flat profiles."""
    profile = np.zeros((2, 5))
    profile[0, 1:4] = [0.2, 0.6, 0.2]
    profile[1, 1:4] = [0.1, 0.8, 0.1]
    sobj = SimpleNamespace(
        SLITID=101,
        OBJID=7,
        BOX_COUNTS_SKY=np.array([10.0, 10.0]),
        flat_corr=np.array([1.0, 2.0]))

    skymodel = FiberExtract._reconstruct_profiled_fiber_skymodel(
        [sobj], {(101, 7): profile})

    np.testing.assert_allclose(skymodel[0, 1:4], [2.0, 6.0, 2.0])
    np.testing.assert_allclose(skymodel[1, 1:4], [2.0, 16.0, 2.0])
    assert np.all(skymodel[:, [0, 4]] == 0.0)


def test_clean_calibration_image_residual_cr():
    """clean_calibration_image must catch a synthetic diagonal CR
    without damaging synthetic arc lines."""
    rng = np.random.default_rng(7)
    image = 10.0 + rng.normal(scale=1.0, size=(120, 160))
    image[35, :] += 500.0
    image[85, :] += 300.0
    xs = np.arange(20, 135)
    ys = np.rint(15 + 0.48 * xs).astype(int)
    image[ys, xs] += 1000.0

    arc_line_pixels = image[[35, 85]].copy()

    bpm_2d = np.zeros(image.shape, dtype=bool)
    fullmask = SimpleNamespace(flagged=lambda flag: bpm_2d.copy())
    calib_image = SimpleNamespace(
        image=image,
        ivar=None,
        fullmask=fullmask,
        update_mask_cr=lambda mask: None,
    )
    process_par = {
        'cr_median_width': 51,
        'sigclip': 10.0,
        'sigfrac': 0.3,
        'objlim': 0.0,
        'lamaxiter': 2,
        'grow': 2.0,
        'rmcompact': False,
    }

    spec.clean_calibration_image(calib_image, "arc", det=2, process_par=process_par)

    assert np.nanmedian(calib_image.image[ys, xs]) < 50.0
    # Arc-line rows are preserved everywhere except where the trail
    # crosses them (at most a handful of pixels per row).
    for row, original in zip([35, 85], arc_line_pixels):
        unchanged = np.isclose(calib_image.image[row, :], original, atol=2.0)
        assert unchanged.sum() >= original.size - 5


def test_identify_fibers_in_block():
    """Given fiber peak positions within a block, identify each fiber."""
    blocks = spec.get_fiber_blocks(1)
    # Use the reference positions as "detected" positions for block 2 (science, 20 fibers)
    block = blocks[1]  # block 2 (0-indexed = 1)
    detected_positions = block['fiber_positions']
    result = spec.identify_fibers_in_block(
        det=1, block_idx=1, detected_positions=detected_positions)
    assert len(result['fiber_id']) == 20
    assert len(result['fiber_name']) == 20
    # All should match since we used reference positions
    assert np.all(result['fiber_id'] > 0)


def test_identify_fibers_with_offset():
    """Fibers should still match with a small position offset."""
    blocks = spec.get_fiber_blocks(1)
    block = blocks[1]
    # Add 2-pixel offset (simulating flexure)
    detected_positions = block['fiber_positions'] + 2.0
    result = spec.identify_fibers_in_block(
        det=1, block_idx=1, detected_positions=detected_positions)
    assert np.all(result['fiber_id'] > 0)
