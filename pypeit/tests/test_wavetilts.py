"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path

import numpy as np

from pypeit.tests.tstutils import data_output_path
from pypeit import wavetilts


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


def test_fit2tiltimg_negative_spat_id():
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
    spat_id = np.array([-12, 150])
    # A gentle spectral ramp so the evaluated tilts vary across the slit
    # (tilts = 0.2 + 0.6 * spec_fraction, within the fit2tilts clamp).
    coeffs = np.zeros((6, 4, 2))
    coeffs[0, 0, :] = 0.5
    coeffs[1, 0, :] = 0.3

    wvtilts = wavetilts.WaveTilts(coeffs=coeffs, nslit=2,
                                  spat_order=np.array([3, 3]),
                                  spec_order=np.array([5, 5]),
                                  spat_id=spat_id, func2d='legendre2d')

    # Build a slitmask: -1 off-slit, spat_id on-slit (mimics slit_img output).
    slitmask = np.full((nspec, nspat), -1, dtype=int)
    slitmask[:, 5:15] = -12
    slitmask[:, 25:35] = 150

    tilts = wvtilts.fit2tiltimg(slitmask)

    neg = slitmask == -12
    pos = slitmask == 150
    # The negative-spat_id slit must be filled with varying (non-zero) tilts,
    # exactly like the positive one.  Under the `>= 0` bug it stays all zero.
    assert np.unique(tilts[neg]).size > 1, \
        'Negative-spat_id slit was not filled (regressed to all-zero tilts)'
    assert not np.allclose(tilts[neg], 0.), \
        'Negative-spat_id slit tilts are all zero'
    assert np.unique(tilts[pos]).size > 1, 'Positive-spat_id slit was not filled'
    # Off-slit pixels remain zero.
    assert np.all(tilts[slitmask == -1] == 0.), 'Off-slit pixels should be zero'

