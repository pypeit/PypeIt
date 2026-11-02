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



def test_fit_arcspec_tiny_fwhm():
    """A fit window below three pixels (measured FWHM < 1 pixel on a slit without arc light) is skipped
    instead of raising TypeError from curve_fit (keck_mosfire long2pos_specphot 2016-01-08)."""
    x = np.arange(50, dtype=float)
    y = np.exp(-0.5 * ((x - 25.0) / 1.5) ** 2) * 100. + 1.
    pixt = np.array([25])
    ampl, cent, widt, centerr = arc.fit_arcspec(x, y, pixt, 0)     # fitp = round(fwhm) = 0
    assert cent[0] == -999. and ampl[0] == -999.
    ampl, cent, widt, centerr = arc.fit_arcspec(x, y, pixt, 6)
    assert abs(cent[0] - 25.0) < 1e-3
