"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path

import numpy as np
import pytest

from pypeit.tests.tstutils import data_output_path
from pypeit import PypeItError
from pypeit import wavetilts
from pypeit.slittrace import SlitTraceSet
from pypeit.core import tracewave


# Test WaveTilts
def test_wavetilts():

    instant_dict = dict(coeffs=np.ones((6,4,1)),
                        nslit=1,
                        spat_order=np.array([3]),
                        spec_order=np.array([5]),
                        spat_id=np.array([150]),
                        func2d='legendre2d')

    wvtilts = wavetilts.WaveTilts(**instant_dict)
    wvtilts.set_paths(data_output_path(''), 'A', '1', 'DET01')
    # I/O
    ofile = Path(wvtilts.get_path()).absolute()
    wvtilts.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    _wvtilts = wavetilts.WaveTilts.from_file(ofile)

    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(wvtilts[key],_wvtilts[key])
        else:
            assert wvtilts[key] == _wvtilts[key]
    # Write again
    wvtilts.to_file(overwrite=True)
    # Clean-up
    ofile.unlink()


def test_fit2tiltimg_negative_spat_id(tmp_path):
    """
    Regression test: fit2tiltimg must fill slits whose spat_id is negative.

    Echelle edge orders can have a negative spat_id (their spatial position
    extrapolates below the detector edge).  SlitTraceSet.slit_img uses -1 as
    the ONLY off-slit sentinel, so the selection of good slits must test
    against -1 -- not `>= 0`, which would silently drop the negative-spat_id
    order, leave its tilts all zero, and later crash global sky subtraction
    with "Infinities in action matrix".
    """
    nspec, nspat = 100, 40
    # Two slits: one with a negative spat_id (edge order), one positive.
    spat_id = np.array([-12, 30])
    # fit2tiltimg reads the slit edges from the Slits calibration file, so
    # build the slits and write them to the calibrations directory.
    left = np.tile(np.array([-20., 25.]), (nspec, 1))
    right = np.tile(np.array([-4., 35.]), (nspec, 1))
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit',
                         nspat=nspat, PYP_SPEC='dummy')
    assert np.array_equal(slits.spat_id, spat_id), 'Unexpected slit spat_id values'
    slits.set_paths(tmp_path, 'A', '1', 'DET01')
    slits.to_file()

    # A gentle spectral ramp so the evaluated tilts vary across the slit
    # (tilts = 0.2 + 0.6 * spec_fraction, within the fit2tilts clamp).
    coeffs = np.zeros((6, 4, 2))
    coeffs[0, 0, :] = 0.5
    coeffs[1, 0, :] = 0.3

    wvtilts = wavetilts.WaveTilts(coeffs=coeffs, nslit=2,
                                  spat_order=np.array([3, 3]),
                                  spec_order=np.array([5, 5]),
                                  spat_id=spat_id, func2d='legendre2d',
                                  slits_filename=Path(slits.get_path()).name)
    wvtilts.set_paths(tmp_path, 'A', '1', 'DET01')

    # Build a slitmask: -1 off-slit, spat_id on-slit (mimics slit_img output).
    slitmask = np.full((nspec, nspat), -1, dtype=int)
    slitmask[:, 5:15] = -12
    slitmask[:, 25:35] = 30

    tilts = wvtilts.fit2tiltimg(slitmask)

    neg = slitmask == -12
    pos = slitmask == 30
    # The negative-spat_id slit must be filled with varying (non-zero) tilts,
    # exactly like the positive one.  Under the `>= 0` bug it stays all zero.
    assert np.unique(tilts[neg]).size > 1, \
        'Negative-spat_id slit was not filled (regressed to all-zero tilts)'
    assert not np.allclose(tilts[neg], 0.), \
        'Negative-spat_id slit tilts are all zero'
    assert np.unique(tilts[pos]).size > 1, 'Positive-spat_id slit was not filled'
    # Off-slit pixels remain zero.
    assert np.all(tilts[slitmask == -1] == 0.), 'Off-slit pixels should be zero'


# Test the spatial flexure with the tilts
def test_fit2tilts():
    """ Test fit2tilts function
    """
    lflex = -0.3  # Spatial flexure of left slit
    rflex = +0.7  # Spatial flexure of right slit
    spat_flexure = np.array([[lflex, rflex]])
    nspec = 40    # Number of spectral pixels
    nspat = 10    # Number of spatial pixels
    left_init = 2.0
    right_init = 7.0
    slits_left = np.full(nspec, left_init)
    slits_right = np.full(nspec, right_init)
    slits = SlitTraceSet(left_init=slits_left,
                         right_init=slits_right,
                         pypeline='MultiSlit',
                         nspat=nspat, PYP_SPEC='dummy')

    left, right, _ = slits.select_edges(initial=True, spat_flexure=spat_flexure)
    slitmask = slits.slit_img(initial=True, spat_flexure=spat_flexure)
    # Create a slit mask
    thismask = (slitmask == slits.spat_id[0])
    _spec_eval, _spat_eval = tracewave.fit2tilts_prepareSlit(slits_left, slits_right, thismask, relative_spat_flexure=spat_flexure[0, :])

    # The spectral coordinate is unchanged by the spatial flexure
    spec_pix, spat_pix = np.where(thismask)
    assert np.allclose(_spec_eval, spec_pix / (nspec - 1)), 'Spectral coordinate changed'
    # Each pixel of the shifted slit must be evaluated at the same fractional position along
    # the slit in the tilts frame; i.e., the left/right edges of the shifted slit map onto the
    # left/right edges of the tilts slit.
    frac = (spat_pix - left[spec_pix, 0]) / (right[spec_pix, 0] - left[spec_pix, 0])
    expected = left_init + frac * (right_init - left_init)
    assert np.allclose(_spat_eval * (nspat - 1), expected), \
        'Spatial coordinate not mapped to the same fractional position along the slit'

    # A uniform flexure of +dx pixels means that a pixel at x must be evaluated at x - dx
    dx = 1.5
    slitmask = slits.slit_img(initial=True, spat_flexure=np.full((1, 2), dx))
    thismask = slitmask == slits.spat_id[0]
    _, _spat_eval = tracewave.fit2tilts_prepareSlit(slits_left, slits_right, thismask,
                                                    relative_spat_flexure=np.full(2, dx))
    assert np.allclose(_spat_eval * (nspat - 1), np.where(thismask)[1] - dx), \
        'Uniform spatial flexure applied with the wrong sign'

    # No flexure means no change of the spatial coordinate
    _, _spat_eval = tracewave.fit2tilts_prepareSlit(slits_left, slits_right, thismask)
    assert np.allclose(_spat_eval * (nspat - 1), np.where(thismask)[1]), \
        'Spatial coordinate changed without spatial flexure'


def test_fit2tilts_prepareSlit_no_flexure():
    # Without a relative spatial flexure, the evaluation coordinates must be
    # *exactly* the pixel coordinates (as used before the spatial flexure was
    # introduced), even for edges that are not at integer pixel positions
    nspec, nspat = 50, 30
    spec = np.arange(nspec)
    slits_left = 3.37 + 0.011 * spec
    slits_right = 21.83 + 0.007 * spec
    spat = np.arange(nspat)[None, :]
    thismask = (spat > slits_left[:, None]) & (spat < slits_right[:, None])
    spec_pix, spat_pix = np.where(thismask)
    for flex in [None, np.zeros(2)]:
        _spec_eval, _spat_eval = tracewave.fit2tilts_prepareSlit(slits_left, slits_right, thismask,
                                                                 relative_spat_flexure=flex)
        assert np.array_equal(_spec_eval, spec_pix / (nspec - 1)), 'Spectral coordinate changed'
        assert np.array_equal(_spat_eval, spat_pix / (nspat - 1)), 'Spatial coordinate changed'


def test_fit2tiltimg_spat_flexure_shift(tmp_path):
    """
    A frame whose spatial flexure is dx pixels larger than that of the tilts frame
    is the tilts frame shifted by dx pixels (the slits and the tilts move together),
    so its tilts image must be the tilts image of the tilts frame, rolled by dx pixels.
    """
    nspec, nspat = 100, 60
    left = np.tile([10., 35.], (nspec, 1))
    right = left + 15.
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit', nspat=nspat,
                         PYP_SPEC='dummy')
    slits.set_paths(tmp_path, 'A', '1', 'DET01')
    slits.to_file()

    # Tilts that vary with spatial position, measured on a tilts frame with its own flexure
    coeffs = np.zeros((4, 4, 2))
    coeffs[0, 0, :] = 0.5
    coeffs[1, 0, :] = 0.4
    coeffs[0, 1, :] = 0.05
    coeffs[1, 1, :] = 0.05
    coeffs[0, 2, :] = 0.03
    tilts_flex = np.full((slits.nslits, 2), 0.5)
    wvtilts = wavetilts.WaveTilts(coeffs=coeffs, nslit=slits.nslits, spat_id=slits.spat_id,
                                  spat_order=np.array([3, 3]), spec_order=np.array([3, 3]),
                                  func2d='legendre2d', spat_flexure=tilts_flex,
                                  slits_filename=Path(slits.get_path()).name)
    wvtilts.set_paths(tmp_path, 'A', '1', 'DET01')

    # Tilts image of the tilts frame
    mask0 = slits.slit_img(initial=True, spat_flexure=tilts_flex)
    t0 = wvtilts.fit2tiltimg(mask0, spat_flexure=tilts_flex)
    # Tilts image of a frame shifted by dx pixels relative to the tilts frame
    dx = 3
    mask1 = slits.slit_img(initial=True, spat_flexure=tilts_flex + dx)
    t1 = wvtilts.fit2tiltimg(mask1, spat_flexure=tilts_flex + dx)

    on = mask1 != -1
    assert np.array_equal(on, np.roll(mask0, dx, axis=1) != -1), 'Slit images inconsistent'
    expected = np.roll(t0, dx, axis=1)
    assert np.allclose(t1[on], expected[on]), 'The tilts must move with the slits'
    # The test is only meaningful if the tilts vary across the slits
    assert not np.allclose(t1[on], np.roll(t0, -dx, axis=1)[on]), 'Test is not sensitive'



def test_fit2tilts_prepareSlit_bad_input():
    nspec, nspat = 20, 10
    thismask = np.zeros((nspec, nspat), dtype=bool)
    thismask[:, 3:6] = True
    left, right = np.full(nspec, 2.), np.full(nspec, 7.)
    with pytest.raises(PypeItError):
        tracewave.fit2tilts_prepareSlit(left, right, thismask, relative_spat_flexure=np.zeros(3))
    with pytest.raises(PypeItError):
        tracewave.fit2tilts_prepareSlit(left[:-1], right, thismask)
    with pytest.raises(PypeItError):
        tracewave.fit2tilts_prepareSlit(left[:-1], right[:-1], thismask)
