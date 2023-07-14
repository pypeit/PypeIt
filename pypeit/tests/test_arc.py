"""
Module to run tests on ararclines
"""
import pytest

from pypeit.core import arc
from pypeit import data

def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = data.load_sky_spectrum('paranal_sky.fits')
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux.value)
    assert (len(arx_w) > 3275)

# Many more functions in pypeit.core.arc that need tests!

