"""
Module to run tests on WaveCalib and WaveFit
"""
import sys
import io
import os
import shutil
import inspect

import pytest

import numpy as np

from pypeit.core.wavecal import wv_fitting
from pypeit.core import fitting

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_wavefit():
    out_file = data_path('test_wavefit.fits')
    if os.path.isfile(out_file):
        os.remove(out_file)
    pypeitFit = fitting.PypeItFit(fitc=np.arange(100).astype(float))
    waveFit = wv_fitting.WaveFit(fit=pypeitFit, pixel_fit=np.arange(100).astype(float))
    # Write
    waveFit.to_file(out_file)
    # Read
    waveFit2 = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.fit.fitc, waveFit2.fit.fitc)
    waveFit2.to_file(out_file, overwrite=True)
    # Finish
    os.remove(out_file)

