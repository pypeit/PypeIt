"""
Module to run tests on simple fitting routines for arrays
"""
import os

import pytest

import numpy as np

from linetools.spectra.io import readspec

import pypeit
from pypeit.core import wave


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_flex_shift():
    # Dummy slf
    # Read spectra
    obj_spec = readspec(data_path('obj_lrisb_600_sky.fits'))
    arx_file = pypeit.__path__[0]+'/data/sky_spec/sky_LRISb_600.fits'
    arx_spec = readspec(arx_file)
    # Call
    flex_dict = wave.flex_shift(obj_spec, arx_spec, mxshft=60)

#    # Apply
#    from scipy import interpolate
#    print(flex_dict['shift'])
#    npix = len(obj_spec.wavelength)
#    x = np.linspace(0., 1., npix)
#    f = interpolate.interp1d(x, obj_spec.wavelength.value, bounds_error=False,
#                             fill_value="extrapolate")
#    new_wave = f(x+flex_dict['shift']/(npix-1))
#
#    from matplotlib import pyplot
#    pyplot.plot(arx_spec.wavelength, arx_spec.flux)
#    pyplot.plot(obj_spec.wavelength, obj_spec.flux)
#    pyplot.plot(new_wave, obj_spec.flux)
#    pyplot.show()
    assert np.abs(flex_dict['shift'] - 43.7) < 0.1

