"""
Module to run tests on skysub routines (mainly for IFU)
"""
from pathlib import Path
from IPython import embed

import numpy as np

from pypeit.core import skysub
from pypeit.images.buildimage import SkyRegions
from pypeit.slittrace import SlitTraceSet
from pypeit.tests.tstutils import data_output_path


def test_userregions():
    """ Run the parameter setup script
    """
    regions = [":10",
               "30:60",
               "80:",
               ":10,35:65,80:",
               "10:20;",
               "10:20,"]
    result = [[np.array([], dtype=int), np.array([99], dtype=int)],
              [np.array([299], dtype=int), np.array([598], dtype=int)],
              [np.array([798], dtype=int), np.array([], dtype=int)],
              [np.array([349, 798], dtype=int), np.array([99, 648], dtype=int)],
              [],
              [np.array([99], dtype=int), np.array([199], dtype=int)]
              ]
    resstat = [0, 0, 0, 0, 1, 2]
    nslits = 2
    maxsl = 100
    for rr, reg in enumerate(regions):
        status, regs = skysub.read_userregions(reg, nslits, maxslitlength=maxsl)
        if status != 1:
            assert len(regs) == nslits, 'Bad number of slits'
            assert np.array_equal(np.where(regs[0][1:] & ~regs[0][:-1])[0], result[rr][0]), \
                    'Bad result for slit 1'
            assert np.array_equal(np.where(~regs[0][1:] & regs[0][:-1])[0], result[rr][1]), \
                    'Bad result for slit 2'
        assert status == resstat[rr], 'Unexpected status'


def test_generatemask():
    maxsl = 1000
    nslits = 2
    reg = "80:"
    tstmsk = np.zeros((1000, 1000), dtype=bool)
    tstmsk[:, 744:901] = True
    status, regs = skysub.read_userregions(reg, nslits, maxslitlength=maxsl)
    assert status == 0, 'Bad generation of user regions'
    slits = SlitTraceSet(left_init=np.full((1000, 1), 120, dtype=float),
                         right_init=np.full((1000, 1), 900, dtype=float), binspec=1, binspat=1,
                         pypeline='IFU', nspat=1000, PYP_SPEC='dummy')
    skymask = skysub.generate_mask("IFU", regs, slits, slits.left_init, slits.right_init)
    assert np.array_equal(skymask, tstmsk), 'Sky mask mismatch'


def test_skyregions_io():
    tstmsk = np.zeros((1000, 1000), dtype=bool)
    tstmsk[:, 744:901] = True
    ofile = Path(data_output_path('test_skyregions.fits')).absolute()
    if ofile.exists():
        ofile.unlink()
    SkyRegions(image=tstmsk.astype(float), PYP_SPEC='dummy').to_file(file_path=ofile)
    assert ofile.exists(), 'SkyRegions file not written'
    reg = SkyRegions.from_file(ofile)
    assert np.array_equal(tstmsk, reg.image.astype(bool)), 'Bad read'

    # Clean-up
    ofile.unlink()


def _make_skyoptimal_inputs(npoly=1, use_spatial=False, seed=1234):
    """
    Build small, deterministic, synthetic 2d inputs for
    :func:`~pypeit.core.skysub.skyoptimal`.

    The ``localmask`` is deliberately irregular (not a simple rectangular
    sub-image) and some pixels inside it have ``ivar <= 0``, so that both the
    boolean-index extraction/scatter and the internal good-pixel
    sub-selection are exercised.
    """
    rng = np.random.default_rng(seed)
    nspec, nspat, nobj = 50, 16, 2

    piximg = np.broadcast_to(np.arange(nspec, dtype=float)[:, None], (nspec, nspat)).copy()
    spatial_img = np.broadcast_to(np.arange(nspat, dtype=float)[None, :], (nspec, nspat)).copy()

    spat = np.arange(nspat, dtype=float)
    centers = [5.0, 10.0]
    sigma = 1.3
    oprof = np.zeros((nspec, nspat, nobj), dtype=float)
    for iobj, cen in enumerate(centers):
        profile = np.exp(-0.5 * ((spat - cen) / sigma) ** 2)
        profile /= profile.sum()
        oprof[:, :, iobj] = np.broadcast_to(profile[None, :], (nspec, nspat))

    sky_true = 4.0 + 0.03 * piximg
    flux_true = [8.0, 5.0]
    noise_sigma = 0.3
    data = sky_true + flux_true[0] * oprof[:, :, 0] + flux_true[1] * oprof[:, :, 1] \
            + rng.normal(scale=noise_sigma, size=(nspec, nspat))
    ivar = np.full((nspec, nspat), 1.0 / noise_sigma ** 2)

    # Irregular local mask: a band of central spatial columns, minus a
    # scattered set of pixels within that band.
    localmask = np.zeros((nspec, nspat), dtype=bool)
    localmask[:, 2:14] = True
    iband = np.where(localmask)
    nband = iband[0].size
    idrop = rng.choice(nband, size=nband // 10, replace=False)
    localmask[iband[0][idrop], iband[1][idrop]] = False

    # Zero out ivar for a separate scattered subset of pixels within the
    # (already trimmed) local mask, to exercise the internal ivar>0
    # sub-selection.
    iband = np.where(localmask)
    nband = iband[0].size
    izero = rng.choice(nband, size=nband // 10, replace=False)
    ivar[iband[0][izero], iband[1][izero]] = 0.0

    fullbkpt = np.linspace(0.0, float(nspec - 1), 11)

    return piximg, data, ivar, oprof, (spatial_img if use_spatial else None), localmask, fullbkpt


def test_skyoptimal():
    """
    Regression/sanity test for :func:`~pypeit.core.skysub.skyoptimal` on
    synthetic 2d data with an irregular ``localmask`` and some ``ivar <= 0``
    pixels inside it.
    """
    for npoly, use_spatial in [(1, False), (3, True)]:
        piximg, data, ivar, oprof, spatial_img, localmask, fullbkpt = \
                _make_skyoptimal_inputs(npoly=npoly, use_spatial=use_spatial)

        sky_bmodel, obj_bmodel, gpm = skysub.skyoptimal(
                piximg, data, ivar, oprof, localmask, spatial_img=spatial_img,
                fullbkpt=fullbkpt, sigrej=3.0, npoly=npoly)

        assert sky_bmodel.shape == piximg.shape
        assert obj_bmodel.shape == piximg.shape
        assert gpm.shape == piximg.shape

        # Outputs are only ever populated within localmask.
        assert not np.any(sky_bmodel[~localmask]), 'sky_bmodel populated outside localmask'
        assert not np.any(obj_bmodel[~localmask]), 'obj_bmodel populated outside localmask'
        assert not np.any(gpm[~localmask]), 'gpm populated outside localmask'

        # gpm can only be True where ivar>0 within localmask.
        assert not np.any(gpm & ~(localmask & (ivar > 0))), 'gpm set for a masked/bad-ivar pixel'

        # The recovered sky model should track the smooth, noise-free sky
        # trend used to build the synthetic data.
        sky_true = 4.0 + 0.03 * piximg
        assert np.allclose(sky_bmodel[localmask], sky_true[localmask], atol=0.5)


def test_skyoptimal_all_masked():
    """
    skyoptimal must degrade to all-zero models and an all-False gpm when
    there are no good pixels to fit.
    """
    piximg, data, ivar, oprof, spatial_img, localmask, fullbkpt = \
            _make_skyoptimal_inputs(npoly=1, use_spatial=False)
    ivar_all_bad = np.zeros_like(ivar)

    sky_bmodel, obj_bmodel, gpm = skysub.skyoptimal(
            piximg, data, ivar_all_bad, oprof, localmask, fullbkpt=fullbkpt, sigrej=3.0, npoly=1)

    assert not np.any(sky_bmodel)
    assert not np.any(obj_bmodel)
    assert not np.any(gpm)
    assert sky_bmodel.shape == piximg.shape
    assert gpm.shape == piximg.shape


