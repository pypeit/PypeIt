"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import numpy as np


from pypeit.tests.tstutils import data_path
from pypeit import wavetilts


instant_dict = dict(coeffs=np.ones((6,4,1)),
                    nslit=1,
                    spat_order=np.array([3]),
                    spec_order=np.array([5]),
                    spat_id=np.array([150]),
                    func2d='legendre2d')

# Test WaveTilts
def test_wavetilts():
    #
    wvtilts = wavetilts.WaveTilts(**instant_dict)
    # I/O
    outfile = data_path('tst_wavetilts.fits')
    wvtilts.to_file(outfile, overwrite=True)
    _wvtilts = wavetilts.WaveTilts.from_file(outfile)

    # Test
    for key in instant_dict.keys():
        if isinstance(instant_dict[key], np.ndarray):
            assert np.array_equal(wvtilts[key],_wvtilts[key])
        else:
            assert wvtilts[key] == _wvtilts[key]
    # Write again
    wvtilts.to_file(outfile, overwrite=True)
    os.remove(outfile)

