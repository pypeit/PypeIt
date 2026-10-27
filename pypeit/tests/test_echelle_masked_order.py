"""Regression coverage for echelle orders masked after wavelength calibration."""

import numpy as np
from types import SimpleNamespace

from pypeit.core import findobj_skymask
from pypeit.find_objects import EchelleFindObjects


def test_ech_objfind_with_masked_middle_order():
    nspec, nspat = 80, 60
    spat = np.arange(nspat)
    image = np.zeros((nspec, nspat), dtype=float)
    slitmask = np.full((nspec, nspat), -1, dtype=int)
    left = np.array([2., 22., 42.])
    right = left + 12.
    spat_ids = np.array([8, 28, 48])
    for iord in range(3):
        slitmask[:, int(left[iord]):int(right[iord])] = spat_ids[iord]
        image += 100. * np.exp(-0.5*((spat - spat_ids[iord])/1.5)**2)

    # Keep the masked order in the slitmask image, as the reduction does, but
    # pass only the two usable orders to object finding.
    keep = np.array([True, False, True])
    sobjs = findobj_skymask.ech_objfind(
        image, np.ones_like(image), slitmask,
        np.tile(left[keep], (nspec, 1)),
        np.tile(right[keep], (nspec, 1)),
        spat_ids[keep], np.array([55, 53]),
        np.tile(np.array([[0], [nspec-1]]), (1, 2)),
        plate_scale=np.full(2, 0.2), ncoeff=2, npca=1,
        snr_thresh=3., min_snr=0., nabove_min_snr=1,
        box_radius=0.6, trim_edg=(1, 1))

    assert len(sobjs) == 2
    assert set(sobjs.ECH_ORDER) == {55, 53}
    assert set(sobjs.SLITID) == {8, 48}


def test_echelle_findobjects_excludes_bad_order_pixels():
    nspec, nspat = 30, 60
    slitmask = np.full((nspec, nspat), -1, dtype=int)
    left = np.array([2., 22., 42.])
    right = left + 12.
    spat_ids = np.array([8, 28, 48])
    for iord in range(3):
        slitmask[:, int(left[iord]):int(right[iord])] = spat_ids[iord]

    finder = EchelleFindObjects.__new__(EchelleFindObjects)
    finder.reduce_bpm = np.array([False, True, False])
    finder.slitmask = slitmask
    finder.slits_left = np.tile(left, (nspec, 1))
    finder.slits_right = np.tile(right, (nspec, 1))
    finder.slits = SimpleNamespace(spat_id=spat_ids)
    finder.order_vec = np.array([55, 54, 53])
    finder.spectrograph = SimpleNamespace(order_platescale=lambda *args, **kwargs: np.full(3, 0.2))
    finder.binning = '1,1'
    finder.sciImg = SimpleNamespace(
        detector=SimpleNamespace(name='DET01'),
        select_flag=lambda **kwargs: np.ones_like(slitmask, dtype=bool))
    finder.objtype = 'science'
    finder.pypeline = 'Echelle'
    finder.par = {'reduce': {'findobj': {'force_center_obj': True}}}
    finder.det = 'DET01'
    finder.steps = []

    sobjs, nobj = finder.find_objects_pypeline(
        np.ones_like(slitmask, dtype=float), np.ones_like(slitmask, dtype=float))
    assert nobj == len(sobjs) == 2
    assert set(sobjs.ECH_ORDER) == {55, 53}
    assert set(sobjs.SLITID) == {8, 48}

    finder.reduce_bpm[:] = True
    sobjs, nobj = finder.find_objects_pypeline(
        np.ones_like(slitmask, dtype=float), np.ones_like(slitmask, dtype=float))
    assert nobj == len(sobjs) == 0
