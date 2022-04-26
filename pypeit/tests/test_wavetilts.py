"""
Module to run tests on WaveTilts and BuildWaveTilts classes
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import numpy as np


from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters, cooked_required
from pypeit import wavetilts
from pypeit import slittrace
from pypeit.core import tracewave, pixels
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@pytest.fixture
@cooked_required
def master_dir():
    return os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'shane_kast_blue')

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

