"""
Module to run tests on SpecObj
"""
import sys
import os

from IPython import embed

import pytest

import numpy as np

from astropy.io import fits

from pypeit import specobj
from pypeit.tests.tstutils import data_path
from pypeit import msgs


def test_init():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    # Test
    assert sobj.PYPELINE == 'MultiSlit'
    assert sobj['PYPELINE'] == 'MultiSlit'
    assert sobj.NAME == 'SPAT-----SLIT0000-DET01'


def test_assignment():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    #
    sobj.PYPELINE = 'Blah'
    # Quick test on datamodel
    with pytest.raises(TypeError):
        sobj.PYPELINE = 2
    #
    sobj.SPAT_PIXPOS = 523.0
    sobj.PYPELINE = 'MultiSlit'
    sobj.set_name()
    assert sobj.NAME == 'SPAT0523-SLIT0000-DET01'


def test_hdu():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    #
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    sobj['BOX_COUNTS'] = np.ones_like(sobj.BOX_WAVE)
    sobj['TRACE_SPAT'] = np.arange(100) * 2.
    # Test
    hdul = sobj.to_hdu()#force_dict_bintbl=True)
    assert len(hdul) == 1  # Should be one BinTableHDU
    assert isinstance(hdul[0], fits.hdu.table.BinTableHDU)
    assert len(hdul[0].data) == 100
    assert 'TRACE_SPAT' in hdul[0].data.dtype.names
    assert 'SLITID' in hdul[0].header.keys()

def test_io():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    # Can we handle 1 array?
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    ofile = data_path('tmp.fits')
    sobj.to_file(ofile, overwrite=True)
    _sobj = specobj.SpecObj.from_file(ofile)
    assert np.array_equal(sobj.BOX_WAVE, _sobj.BOX_WAVE)

    # Cleanup
    os.remove(ofile)

def test_iotwo():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    #
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    sobj['BOX_COUNTS'] = np.ones_like(sobj.BOX_WAVE)
    sobj['TRACE_SPAT'] = np.arange(100) * 2.
    sobj['BOX_MASK'] = np.arange(100).astype(bool)

    # Write table
    ofile = data_path('tmp.fits')
    sobj.to_file(ofile, overwrite=True)
    sobj2 = specobj.SpecObj.from_file(ofile)
    #
    assert isinstance(sobj2, specobj.SpecObj)
    assert np.array_equal(sobj.BOX_WAVE, sobj2.BOX_WAVE)
    assert sobj2.PYPELINE == 'MultiSlit'

    # Cleanup
    os.remove(ofile)

def test_copy():
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    #
    sobj['BOX_WAVE'] = np.arange(100).astype(float)
    sobj.smash_snr = 1.
    # Copy
    sobj2 = specobj.SpecObj.copy(sobj)
    assert np.isclose(sobj2.smash_snr, 1.)
    # Check
    assert np.array_equal(sobj.BOX_WAVE, sobj2.BOX_WAVE)


def test_from_arrays():
    wave = np.linspace(5000., 6000, 1000)
    flux = np.ones_like(wave)
    ivar = 0.1*np.ones_like(wave)
    sobj = specobj.SpecObj.from_arrays('MultiSlit', wave, flux, ivar)

    assert sobj.OPT_WAVE[0] == 5000.


