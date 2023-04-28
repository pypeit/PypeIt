"""
Module to run tests on WaveCalib and WaveFit
"""
from pathlib import Path
import os

from IPython import embed

import pytest

import numpy as np

from pypeit.core.wavecal import wv_fitting
from pypeit.core import fitting
from pypeit import wavecalib
from pypeit import slittrace
from pypeit.tests.tstutils import data_path


def test_wavefit():
    "Fuss with the WaveFit DataContainer"
    out_file = Path(data_path('test_wavefit.fits')).resolve()
    if out_file.exists():
        out_file.unlink()
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float))
    #
    ions = np.asarray(['CdI', 'HgI'])
    ion_bits = np.zeros(len(ions), dtype=wv_fitting.WaveFit.bitmask.minimum_dtype())
    for kk,ion in enumerate(ions):
        ion_bits[kk] = wv_fitting.WaveFit.bitmask.turn_on(ion_bits[kk], ion)
    # Do it
    waveFit = wv_fitting.WaveFit(99, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10), sigrej=3.,
                                 ion_bits=ion_bits)

    # Write
    waveFit.to_file(out_file)
    # Read
    waveFit2 = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.pypeitfit.fitc, waveFit2.pypeitfit.fitc)
    # Write again
    waveFit2.to_file(out_file, overwrite=True)
    # And one more read
    waveFit2b = wv_fitting.WaveFit.from_file(out_file)
    assert np.array_equal(waveFit.pypeitfit.fitc, waveFit2b.pypeitfit.fitc)
    # Finish
    out_file.unlink()

    # No fit
    waveFit3 = wv_fitting.WaveFit(99, pypeitfit=None, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10), sigrej=3.,
                                 ion_bits=ion_bits)
    waveFit3.to_file(out_file)
    waveFit4 = wv_fitting.WaveFit.from_file(out_file)
    assert waveFit4.pypeitfit is None

    # Clean up
    out_file.unlink()


def test_wavecalib():
    "Fuss with the WaveCalib DataContainer"
    # Pieces
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float), xval=np.linspace(1,100., 100))
    # 2D fit
    pypeitFit2 = fitting.PypeItFit(fitc=np.linspace((1,2),(10,20),10),
                                   xval=np.linspace(1,100., 100),
                                   x2=np.linspace(1, 100., 100))
    waveFit = wv_fitting.WaveFit(232, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10))

    waveCalib = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit]),
                                    nslits=1, spat_ids=np.asarray([232]),
                                    wv_fit2d=np.array([pypeitFit2]))
    waveCalib.set_paths(data_path(''), 'A', '1', 'DET01')

    ofile = Path(waveCalib.get_path()).resolve()

    # Write
    waveCalib.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    # Read
    waveCalib2 = wavecalib.WaveCalib.from_file(ofile)

    # Test
    assert np.array_equal(waveCalib.spat_ids, waveCalib2.spat_ids), 'Bad spat_ids'
    assert np.array_equal(waveCalib.wv_fits[0].pypeitfit.fitc,
                          waveCalib2.wv_fits[0].pypeitfit.fitc), 'Bad fitc'
    assert np.array_equal(waveCalib.wv_fit2d[0].xval, waveCalib2.wv_fit2d[0].xval)

    # Write again!
    waveCalib2.to_file(overwrite=True)

    # Finish
    ofile.unlink()

    # With None (failed wave)
    spat_ids = np.asarray([232, 949])
    waveCalib3 = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit, wv_fitting.WaveFit(949)]),
                                    nslits=2, spat_ids=spat_ids,
                                    wv_fit2d=np.array([pypeitFit2, pypeitFit2]))
    waveCalib3.set_paths(data_path(''), 'A', '1', 'DET01')
    waveCalib3.to_file(overwrite=True)
    waveCalib4 = wavecalib.WaveCalib.from_file(ofile)

    # Check masking
    slits = slittrace.SlitTraceSet(left_init=np.full((1000,2), 2, dtype=float),
                         right_init=np.full((1000,2), 8, dtype=float),
                         pypeline='MultiSlit', spat_id=spat_ids,
                         nspat=2, PYP_SPEC='dummy')
    slits.mask_wvcalib(waveCalib3)
    assert slits.bitmask.flagged(slits.mask[1], flag='BADWVCALIB')

    # Finish
    ofile.unlink()


def test_wvcalib_no2d():
    """ Test WaveCalib without a 2D fit (not echelle) """
    # Pieces
    pypeitFit = fitting.PypeItFit(fitc=np.arange(5).astype(float), xval=np.linspace(1,100., 100))
    waveFit = wv_fitting.WaveFit(232, pypeitfit=pypeitFit, pixel_fit=np.arange(10).astype(float),
                                 wave_fit=np.linspace(1.,10.,10))

    # With None (failed wave)
    spat_ids = np.asarray([232, 949])
    waveCalib = wavecalib.WaveCalib(wv_fits=np.asarray([waveFit, wv_fitting.WaveFit(949)]),
                                    nslits=2, spat_ids=spat_ids)
    waveCalib.set_paths(data_path(''), 'A', '1', 'DET01')
    ofile = Path(waveCalib.get_path()).resolve()

    waveCalib.to_file(overwrite=True)
    assert ofile.exists(), 'File not written'

    waveCalib2 = wavecalib.WaveCalib.from_file(ofile)
    assert len(waveCalib2.wv_fits) == 2, 'Should be two wavelength fits'
    assert waveCalib2.wv_fits[0].spat_id == 232, 'Wrong spat_id'
    assert np.array_equal(waveCalib2.wv_fits[0].wave_fit, waveCalib.wv_fits[0].wave_fit), \
            'Bad wave_fit read'

    # Finish
    ofile.unlink()

