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
from pypeit import wavecalib
from pypeit import slittrace

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_wavefit():
    "Fuss with the WaveFit DataContainer"
    out_file = data_path('test_wavefit.fits')
    if os.path.isfile(out_file):
        os.remove(out_file)
    pypeitFit = fitting.PypeItFit(fitc=np.arange(100).astype(float))
    waveFit = wv_fitting.WaveFit(pypeitfit=pypeitFit, pixel_fit=np.arange(100).astype(float))

    # Write
    waveFit.to_file(out_file)
    # Read
    waveFit2 = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.pypeitfit.fitc, waveFit2.pypeitfit.fitc)
    waveFit2.to_file(out_file, overwrite=True)
    # Finish
    os.remove(out_file)


def test_wavecalib():
    "Fuss with the WaveCalib DataContainer"
    out_file = data_path('test_wavecalib.fits')
    if os.path.isfile(out_file):
        os.remove(out_file)
    # Piecese
    pypeitFit = fitting.PypeItFit(fitc=np.arange(100).astype(float))
    waveFit = wv_fitting.WaveFit(pypeitfit=pypeitFit, pixel_fit=np.arange(100).astype(float))

    waveCalib = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit]), spat_id=np.asarray([232]))

    # Write
    waveCalib.to_file(out_file)

    # Read
    waveCalib2 = wavecalib.WaveCalib.from_file(out_file)

    # Test
    assert np.array_equal(waveCalib.spat_id, waveCalib2.spat_id), 'Bad spat_id'
    assert np.array_equal(waveCalib.wv_fits[0].pypeitfit.fitc,
                          waveCalib2.wv_fits[0].pypeitfit.fitc), 'Bad fitc'

    # Finish
    os.remove(out_file)
