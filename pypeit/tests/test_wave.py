
from IPython import embed

import numpy as np

from astropy import units

from pypeit.core.wave import airtovac

def test_airtovac():
    air = 5000.
    vac = airtovac(air * units.AA).value
    assert isinstance(vac, (float, np.floating)), 'Should return scalar'
    assert np.isclose(vac, 5001.4), 'Bad value'

    air = [5000., 5100.]
    _vac = airtovac(air * units.AA).value
    assert isinstance(_vac, np.ndarray), 'Should return array'
    assert np.isclose(vac, _vac[0]), 'Bad value'

    air = [1950.]
    vac = airtovac(air * units.AA).value
    assert isinstance(vac, np.ndarray), 'Should return array'
    assert np.isclose(vac[0], air[0]), 'Should not alter wave < 2000'
