"""
Module to test functions in pypeit.coadd2d
"""
from pathlib import Path

import numpy as np
import pytest

from astropy.io import fits
from astropy.table import Table

from pypeit import PypeItError
from pypeit.coadd2d import CoAdd2dStack, MultiSlitCoAdd2D
from pypeit.core import coadd as core_coadd
from pypeit.images.imagebitmask import ImageBitMaskArray
from pypeit.slittrace import SlitTraceSet
from pypeit.spec2dobj import Spec2DObj
from pypeit.spectrographs.util import load_spectrograph

# Spectrograph used to provide the detector and the parameters of the mock coadds
SPECTROGRAPH = load_spectrograph('shane_kast_blue')


@pytest.fixture
def coadd():
    # Multislit coadd of two exposures.  Skip __init__, which needs real spec2d
    # files, and set nexp directly; the tests set any other attribute they need
    # (e.g., coadd2d_stack, see _coadd2d_stack)
    _coadd = MultiSlitCoAdd2D.__new__(MultiSlitCoAdd2D)
    _coadd.nexp = 2
    return _coadd


def test_check_input_list(coadd):
    weights = coadd.check_input([0.3, np.float32(0.7)], 'weights')
    assert isinstance(weights, list), f'weights should be returned as a list, got {type(weights)}'
    assert np.allclose(weights, [0.3, 0.7]), f'weights values changed: {weights}'
    offsets = coadd.check_input(np.array([0, -5]), 'offsets')
    assert isinstance(offsets, np.ndarray), f'offsets should be returned as an ndarray, got {type(offsets)}'
    assert np.array_equal(offsets, [0, -5]), f'offsets values changed: {offsets}'

    # Wrong number of elements
    with pytest.raises(PypeItError, match='same number of elements'):
        coadd.check_input([1.0, 2.0, 3.0], 'offsets')

    # Lists with the right length but non-numerical (or boolean) elements are rejected
    for bad in [['[1.0', '2.0'], [True, False], [1.0, None]]:
        with pytest.raises(PypeItError, match='Unrecognized format'):
            coadd.check_input(bad, 'weights')

    # Unknown type
    with pytest.raises(PypeItError, match='Unrecognized type'):
        coadd.check_input([1.0, 2.0], 'shifts')


def _coadd2d_stack(platescales, slits_list=None, maskdef_designtab=None):
    # Build a CoAdd2dStack with only the detectors (one per exposure, each
    # with its own platescale), slits, and slitmask design needed by the tests
    nexp = len(platescales)
    detectors = []
    for platescale in platescales:
        detector = SPECTROGRAPH.get_detector_par(1)
        detector.platescale = platescale
        detectors.append(detector)
    return CoAdd2dStack(specobjs_list=[], slits_list=[] if slits_list is None else slits_list,
                        slitmask_stack=[], sciimg_stack=[], sciivar_stack=[], skymodel_stack=[],
                        mask_stack=[], waveimg_stack=[], exptime_stack=[], exptime_coadd=None,
                        redux_path=Path.cwd(), detectors=detectors,
                        spectrograph=SPECTROGRAPH.name, pypeline=SPECTROGRAPH.pypeline,
                        maskdef_designtab_list=[maskdef_designtab]*nexp,
                        spat_flexure_list=[None]*nexp)


def test_offsets_same_pixscale(coadd):
    # With equal pixel scales, the offsets are the usual difference of positions
    coadd.coadd2d_stack = _coadd2d_stack([0.1, 0.1])
    assert np.allclose(coadd.pixscale_ratios, 1.)
    positions = np.array([3.0, -2.5])
    assert np.allclose(coadd.offsets_from_positions(positions), positions[0] - positions)


def test_offsets_and_spat_bins_different_pixscale(coadd):
    platescales = np.array([0.100, 0.103])
    coadd.coadd2d_stack = _coadd2d_stack(platescales)
    assert np.allclose(coadd.pixscales, platescales)
    assert np.allclose(coadd.pixscale_ratios, platescales / platescales[0])

    # Slit center (pixels) and dither offset (arcsec) of each exposure; the
    # source is at -dithoff arcsec from the slit center
    slitcen = np.array([50.3, 47.8])
    dithoff = np.array([-1.0, 0.5])
    srcpos = slitcen - dithoff / platescales

    # Offsets in pixels of each exposure; they reproduce the dither in arcsec
    offsets = coadd.offsets_from_positions(-dithoff / platescales)
    assert np.allclose(offsets * platescales, dithoff - dithoff[0])

    nspec, nspat = 5, 100
    thismask_stack = [np.ones((nspec, nspat), dtype=bool)] * 2
    ref_trace_stack = [np.full(nspec, c - o) for c, o in zip(slitcen, offsets)]
    _, dspat_stack = core_coadd.get_spat_bins(thismask_stack, ref_trace_stack,
                                              spat_scale=coadd.pixscale_ratios)

    # Every pixel must land on the coadd grid (in pixels of the first exposure)
    # at a position set only by its angular distance from the source
    xpix = np.arange(nspat)
    for i in range(2):
        sky = (xpix - srcpos[i]) * platescales[i] / platescales[0]
        assert np.allclose(dspat_stack[i][0] - sky, srcpos[0] - slitcen[0])

    # Without accounting for the pixel scales the alignment away from the
    # reference position degrades
    _, dspat_stack = core_coadd.get_spat_bins(thismask_stack, ref_trace_stack)
    sky = (xpix - srcpos[1]) * platescales[1] / platescales[0]
    assert not np.allclose(dspat_stack[1][0] - sky, srcpos[0] - slitcen[0])

    # Wrong number of scale factors
    with pytest.raises(PypeItError, match='spat_scale'):
        core_coadd.get_spat_bins(thismask_stack, ref_trace_stack, spat_scale=[1.])


def test_maskdef_offsets_different_pixscale(coadd):
    # Exercise the `offsets = maskdef_offsets` path of a multislit coadd with
    # two exposures having different platescales, with the slitmask-design
    # quantities built as in SlitTraceSet.get_maskdef_objpos
    nspec, nspat = 5, 120
    platescales = np.array([0.100, 0.104])
    # design distances (arcsec) of the object from the slit edges
    topdist, botdist = 3.0, 2.0
    # measured slit center, edge loss on each side, and slitmask offset (pixels)
    slitcen = np.array([60.2, 55.7])
    edgeloss = np.array([0.5, 1.2])
    maskdef_offset = np.array([1.3, -2.1])

    slitlen = (topdist + botdist) / platescales - 2 * edgeloss
    left = slitcen - slitlen / 2
    # expected position of the object from the left edge
    objpos = topdist / platescales - edgeloss
    # measured position of the object
    srcpos = left + objpos + maskdef_offset

    slits_list = [SlitTraceSet(np.full((nspec, 1), left[i]), np.full((nspec, 1), left[i] + slitlen[i]),
                               SPECTROGRAPH.pypeline, nspat=nspat, maskdef_id=np.array([7]),
                               maskdef_objpos=np.array([objpos[i]]),
                               maskdef_slitcen=np.full((nspec, 1), slitcen[i]),
                               maskdef_offset=maskdef_offset[i])
                  for i in range(2)]
    spat_ids = slits_list[0].spat_id
    designtab = Table({'SPAT_ID': spat_ids, 'TRACEID': [0], 'TRACESROW': [0], 'TRACELPIX': [0.],
                       'TRACERPIX': [0.], 'SLITLMASKDEF': [0.], 'SLITRMASKDEF': [0.]})

    coadd.coadd2d_stack = _coadd2d_stack(platescales, slits_list=slits_list,
                                         maskdef_designtab=designtab)
    coadd.obj_id_bri = None
    coadd.spat_ids = spat_ids
    coadd.nspec_array = [nspec]
    coadd.par = SPECTROGRAPH.default_pypeit_par()
    coadd.par['calibrations']['slitedges']['use_maskdesign'] = True
    coadd.par['coadd2d']['offsets'] = 'maskdef_offsets'

    coadd.compute_offsets()
    # offsets, converted to arcsec, are the difference of the slitmask offsets
    assert np.allclose(coadd.offsets * platescales,
                       maskdef_offset[0] * platescales[0] - maskdef_offset * platescales)

    # the object lands at the same position of the coadd spatial grid
    ref_trace_stack = coadd.reference_trace_stack(0, offsets=coadd.offsets)
    thismask_stack = [np.ones((nspec, nspat), dtype=bool)] * 2
    _, dspat_stack = core_coadd.get_spat_bins(thismask_stack, ref_trace_stack,
                                              spat_scale=coadd.pixscale_ratios)
    xpix = np.arange(nspat)
    obj_dspat = [np.interp(srcpos[i], xpix, dspat_stack[i][0]) for i in range(2)]
    assert np.isclose(obj_dspat[0], obj_dspat[1])

    # the expected object position and slit center propagated to the coadd are
    # consistent among the exposures, i.e. their mean equals the value of the
    # first exposure
    maskdef_dict = coadd.get_maskdef_dict(0, ref_trace_stack)
    assert np.isclose(maskdef_dict['maskdef_objpos'], obj_dspat[0])
    assert np.isclose(maskdef_dict['maskdef_slitcen'], maskdef_offset[0])


def _spec2dobj(exptime, units):
    # A minimal Spec2DObj with one slit, a given exposure time, and given units
    shape = (100, 50)
    left = np.full((shape[0], 1), 10., dtype=float)
    right = np.full((shape[0], 1), 40., dtype=float)
    slits = SlitTraceSet(left, right, 'MultiSlit', nspat=shape[1], PYP_SPEC=SPECTROGRAPH.name)
    flex = Table()
    flex['spat_id'] = slits.spat_id
    flex['sci_spec_flexure'] = np.zeros(1)
    s2dobj = Spec2DObj(sciimg=np.full(shape, 10.), ivarraw=np.ones(shape), skymodel=np.ones(shape),
                       bkg_redux_skymodel=None, objmodel=np.zeros(shape), ivarmodel=np.ones(shape),
                       scaleimg=np.ones(shape), waveimg=np.full(shape, 5000.),
                       bpmmask=ImageBitMaskArray(shape), detector=SPECTROGRAPH.get_detector_par(1),
                       sci_spat_flexure=None, sci_spec_flexure=flex, vel_type=None, vel_corr=None,
                       slits=slits, wavesol=None, tilts=np.ones(shape), maskdef_designtab=None)
    s2dobj.head0 = fits.Header({'EXPTIME': exptime})
    s2dobj.units = units
    return s2dobj


def test_exptime_rescaling(coadd):
    # Images in detector counts with different exposure times are rescaled to
    # the median exposure time; flux-calibrated images are not
    coadd.detname = 'DET01'
    coadd.spectrograph = SPECTROGRAPH
    for units, rescaled in [('e-', True), (None, True), ('MJy', False)]:
        stack = coadd.load_coadd2d_stacks([_spec2dobj(100., units), _spec2dobj(200., units)])
        assert stack.exptime_coadd == 200., 'The reference exposure time should be the median.'
        assert stack.units == units, 'The stack should have the units of the input images.'
        expected = 20. if rescaled else 10.
        assert np.allclose(stack.sciimg_stack[0], expected), \
            f'Images in {units} should {"" if rescaled else "not "}be rescaled to the ' \
            'reference exposure time.'
        assert np.allclose(stack.sciimg_stack[1], 10.), \
            'The image with the reference exposure time should not be rescaled.'

    with pytest.raises(PypeItError, match='different units'):
        coadd.load_coadd2d_stacks([_spec2dobj(100., 'e-'), _spec2dobj(100., 'MJy')])
