"""
Module to run tests on simple fitting routines for arrays
"""
import os

import pytest

import numpy as np

from linetools.spectra.io import readspec

from pypeit.core import flexure, arc
from pypeit import data

from pypeit.tests.tstutils import data_path

from IPython import embed


def test_flex_shift():
    # Dummy slf
    # Read spectra
    obj_spec = readspec(data_path('obj_lrisb_600_sky.fits'))
    arx_file = os.path.join(data.Paths.sky_spec, 'sky_LRISb_600.fits')
    arx_spec = readspec(arx_file)
    arx_lines = arc.detect_lines(arx_spec.flux.value)

    # Call
    flex_dict = flexure.spec_flex_shift(obj_spec, arx_spec, arx_lines, mxshft=60)

    assert np.abs(flex_dict['shift'] - 43.7) < 0.1
