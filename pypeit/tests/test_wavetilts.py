"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path

import numpy as np

from pypeit.tests.tstutils import data_output_path
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

