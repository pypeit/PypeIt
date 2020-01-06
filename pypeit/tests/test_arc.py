"""
Module to run tests on ararclines
"""
import os
import numpy as np
import pytest

from linetools.spectra import xspectrum1d

import pypeit
from pypeit.core import arc

import pkg_resources

def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    sky_file = pkg_resources.resource_filename('pypeit', 'data/sky_spec/paranal_sky.fits')
    arx_sky = xspectrum1d.XSpectrum1D.from_file(sky_file)
    arx_amp_true, arx_amp, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux.value)
    assert (len(arx_w[0]) > 3275)

# Many more functions in pypeit.core.arc that need tests!

