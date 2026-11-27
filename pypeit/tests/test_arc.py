"""
Module to run tests on ararclines
"""
import pytest
import numpy as np

from pypeit.core import arc
from pypeit.core import skyspec


def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = skyspec.load_sky_spectrum('paranal_sky.fits')
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux)
    assert (len(arx_w) > 3275)

# Many more functions in pypeit.core.arc that need tests!

def test_fit2darc_requires_two_unique_orders():
    import numpy as np
    import pytest
    from pypeit import PypeItWavelengthCalibrationError
    from pypeit.core import arc

    for orders in (np.array([]), np.array([30., 30.])):
        with pytest.raises(PypeItWavelengthCalibrationError,
                           match='needs >= 2 unique orders'):
            arc.fit2darc(np.ones(orders.size), np.arange(orders.size),
                         orders, nspec=2000)
