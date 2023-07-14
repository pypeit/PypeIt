"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
from pathlib import Path

import numpy as np

from pypeit.tests.tstutils import data_path
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
    wvtilts.set_paths(data_path(''), 'A', '1', 'DET01')
    # I/O
    ofile = Path(wvtilts.get_path()).resolve()
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

