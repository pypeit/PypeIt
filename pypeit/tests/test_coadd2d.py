"""
Module to run tests on the 2D coadd utilities that do not require data files.
"""
from types import SimpleNamespace

import numpy as np
import pytest
from astropy.table import Table

from pypeit import PypeItError
from pypeit.coadd2d import CoAdd2D, MultiSlitCoAdd2D
from pypeit.slittrace import SlitTraceSet


@pytest.fixture
def coadd():
    # check_input only depends on the number of exposures, so skip __init__, which
    # needs real spec2d files, and set nexp directly
    _coadd = CoAdd2D.__new__(CoAdd2D)
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

def _slits(nspec=100):
    left = np.stack([np.linspace(20., 23., nspec), np.linspace(80., 83., nspec)], axis=1)
    right = left + 40.
    return SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit', nspat=200,
                        PYP_SPEC='dummy')


def test_slit_centers():
    slits = _slits()
    assert np.array_equal(CoAdd2D.slit_centers(slits), slits.center)
    # The slit centers follow the spatial flexure of each slit edge
    spat_flexure = np.array([[1.0, 2.0], [-0.5, 0.5]])
    assert np.allclose(CoAdd2D.slit_centers(slits, spat_flexure),
                       slits.center + np.mean(spat_flexure, axis=1)[None, :])
    # The global padding of the slit edges does not change the slit centers
    slits.pad = np.array([5., 1.])
    assert np.allclose(CoAdd2D.slit_centers(slits), slits.center)


def test_offset_slit_cen_spat_flexure():
    # Two exposures with different spatial flexure and a dither.  The offsets
    # (e.g. from the header, the user, or the slitmask design) only describe
    # the dither, so the reference traces must include the spatial flexure of
    # each exposure for the object to be at the same position relative to the
    # reference trace in both exposures.
    slits = _slits()
    spat_flexure = [np.zeros((2, 2)), np.full((2, 2), 2.5)]
    dither = np.array([0., -3.])
    offsets = dither[0] - dither
    obj_pos = 5.
    coadd = SimpleNamespace(coadd2d_stack=SimpleNamespace(slits_list=[slits, slits],
                                                          spat_flexure_list=spat_flexure),
                            slit_centers=CoAdd2D.slit_centers)
    for slit_idx in range(slits.nslits):
        ref_traces = CoAdd2D.offset_slit_cen(coadd, slit_idx, offsets)
        obj_traces = [slits.center[:, slit_idx] + spat_flexure[iexp][slit_idx, 0] + dither[iexp]
                      + obj_pos for iexp in range(2)]
        rel = [obj - ref for obj, ref in zip(obj_traces, ref_traces)]
        assert np.allclose(rel[0], rel[1]), 'Exposures are misregistered by the spatial flexure'
    # Exposures without spatial flexure
    coadd.coadd2d_stack.spat_flexure_list = [None, None]
    ref_traces = CoAdd2D.offset_slit_cen(coadd, 0, offsets)
    assert np.allclose(ref_traces[1], slits.center[:, 0] - offsets[1])


def _maskdef_coadd(spat_flexure_list, maskdef_offset):
    # A MultiSlitCoAdd2D with only the attributes used by get_maskdef_dict and
    # offset_slit_cen (two exposures of the same two slits)
    slits = _slits()
    slits.maskdef_id = np.array([101, 102])
    slits.maskdef_objpos = np.array([12., 15.])
    slits.maskdef_slitcen = slits.center + np.array([0.3, -0.4])[None, :]
    designtab = Table({'SPAT_ID': slits.spat_id, 'MASKDEF_ID': slits.maskdef_id})
    for col in ['TRACEID', 'TRACESROW', 'TRACELPIX', 'TRACERPIX', 'SLITLMASKDEF', 'SLITRMASKDEF']:
        designtab[col] = np.zeros(slits.nslits)
    coadd = MultiSlitCoAdd2D.__new__(MultiSlitCoAdd2D)
    coadd.par = {'calibrations': {'slitedges': {'use_maskdesign': True}},
                 'coadd2d': {'offsets': 'maskdef_offsets', 'spat_samp_fact': 1.0}}
    coadd.coadd2d_stack = SimpleNamespace(slits_list=[slits, slits],
                                          spat_flexure_list=spat_flexure_list,
                                          maskdef_designtab_list=[designtab, designtab.copy()])
    coadd.spat_ids = slits.spat_id
    coadd.nexp = 2
    coadd.nspec_array = [slits.nspec, slits.nspec]
    coadd.maskdef_offset = np.asarray(maskdef_offset)
    # Offsets as computed by compute_offsets for offsets = 'maskdef_offsets'
    offsets = coadd.maskdef_offset[0] - coadd.maskdef_offset
    return coadd, offsets


def test_maskdef_dict_spat_flexure():
    # The slitmask-design slit center and target position in the coadded frame
    # must not depend on the spatial flexure of the individual exposures: the
    # maskdef offsets do not include the flexure, while the reference traces
    # (offset_slit_cen) do.
    maskdef_offset = [0., -3.]
    ref = {}
    for flexure in [[None, None], [np.zeros((2, 2)), np.zeros((2, 2))]]:
        coadd, offsets = _maskdef_coadd(flexure, maskdef_offset)
        for slit_idx in range(2):
            d = coadd.get_maskdef_dict(slit_idx, coadd.offset_slit_cen(slit_idx, offsets))
            if slit_idx not in ref:
                ref[slit_idx] = d
            # No flexure (None) and zero flexure give the same result
            assert np.isclose(d['maskdef_slitcen'], ref[slit_idx]['maskdef_slitcen'])
            assert np.isclose(d['maskdef_objpos'], ref[slit_idx]['maskdef_objpos'])
            assert d['maskdef_id'] == 101 + slit_idx

    # Different flexure for each exposure and slit (same for both edges of a slit)
    flexure = [np.array([[0.5, 0.5], [-1.0, -1.0]]), np.array([[2.5, 2.5], [1.5, 1.5]])]
    coadd, offsets = _maskdef_coadd(flexure, maskdef_offset)
    for slit_idx in range(2):
        d = coadd.get_maskdef_dict(slit_idx, coadd.offset_slit_cen(slit_idx, offsets))
        assert np.isclose(d['maskdef_slitcen'], ref[slit_idx]['maskdef_slitcen']), \
            'The slitmask slit center in the coadd depends on the spatial flexure'
        assert np.isclose(d['maskdef_objpos'], ref[slit_idx]['maskdef_objpos']), \
            'The slitmask target position in the coadd depends on the spatial flexure'

    # Different flexure for the left and right edges: the slit center is still
    # independent of the flexure
    flexure = [np.array([[0.5, 1.5], [-1.0, 0.0]]), np.array([[2.5, 1.0], [1.5, 2.0]])]
    coadd, offsets = _maskdef_coadd(flexure, maskdef_offset)
    for slit_idx in range(2):
        d = coadd.get_maskdef_dict(slit_idx, coadd.offset_slit_cen(slit_idx, offsets))
        assert np.isclose(d['maskdef_slitcen'], ref[slit_idx]['maskdef_slitcen'])
