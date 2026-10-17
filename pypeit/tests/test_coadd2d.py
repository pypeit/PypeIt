"""
Module to run tests on the 2D coadd utilities that do not require data files.
"""
from types import SimpleNamespace

import numpy as np

from pypeit.coadd2d import CoAdd2D
from pypeit.slittrace import SlitTraceSet


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
