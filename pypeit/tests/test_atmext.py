"""
Module to run tests on the atmospheric extinction class.
"""
from IPython import embed

import numpy as np
import pytest

from pypeit import telescopes
from pypeit import PypeItError
from pypeit.core.atmextinction import AtmosphericExtinction


def test_instantiate():
    wave = np.array([3200., 3250., 3300., 3350., 3400.])
    mag_ext = np.array([1.084, 0.948, 0.858, 0.794, 0.745])

    atmext = AtmosphericExtinction(wave, mag_ext)
    assert wave[0] == atmext.wave[0], 'Bad wave assignment'
    assert mag_ext[0] == atmext.mag_ext[0], 'Bad extinction assignment'
    assert atmext.size == wave.size, 'Bad length'

    indx = np.arange(wave.size)
    np.random.shuffle(indx)
    atmext = AtmosphericExtinction(wave[indx], mag_ext[indx])
    assert wave[indx[0]] == atmext.wave[0], 'Bad wave assignment'

    atmext = AtmosphericExtinction(wave[indx], mag_ext[indx], assume_sorted=False)
    assert wave[0] == atmext.wave[0], 'Bad wave assignment'


def test_find_closest():
    mtham = telescopes.ShaneTelescopePar()
    _file = AtmosphericExtinction.closest_extinction_file(mtham['longitude'], mtham['latitude'])
    assert _file == 'mthamextinct.dat', 'Found the wrong file'

    # Fail on location in the western Mediterranean Sea
    with pytest.raises(PypeItError):
        _file = AtmosphericExtinction.closest_extinction_file(0., 37.3413889)


def test_from_coordinates():

    mtham = telescopes.ShaneTelescopePar()
    atmext = AtmosphericExtinction.from_coordinates(mtham['longitude'], mtham['latitude'])
    assert atmext.file == 'mthamextinct.dat', 'Found the wrong file'

    np.testing.assert_allclose(atmext.wave[0], 3200.)
    np.testing.assert_allclose(atmext.mag_ext[0], 1.084)

    # Fail on location in the western Mediterranean Sea
    with pytest.raises(PypeItError):
        atmext = AtmosphericExtinction.from_coordinates(0., 37.3413889)


def test_from_file():

    mtham = telescopes.ShaneTelescopePar()
    _file = AtmosphericExtinction.closest_extinction_file(mtham['longitude'], mtham['latitude'])
    assert _file == 'mthamextinct.dat', 'Found the wrong file'
    atmext = AtmosphericExtinction.from_file(_file)
    assert atmext.file == _file, 'File name not recorded correctly'

    np.testing.assert_allclose(atmext.wave[0], 3200.)
    np.testing.assert_allclose(atmext.mag_ext[0], 1.084)

    # Fail on nonexistant extinction filename
    with pytest.raises(PypeItError):
        extinct = AtmosphericExtinction.from_file('northpoleextinct.dat')


def test_correction():
    # Load
    atmext = AtmosphericExtinction.from_file('mthamextinct.dat')

    # Wavelength vector
    wave = np.arange(3000.,10000.,5.)

    # Test
    corr = atmext.correction_factor(wave, airmass=1.5)
    assert np.round(corr[0], decimals=3) == 4.471, 'Bad correction factor'

    flux = np.ones(wave.size, dtype=float)
    cflux = atmext.correct(flux, corr)
    assert np.round(cflux[0], decimals=3) == 4.471, 'Bad correction'

    cflux = atmext.correct(2*flux, corr)
    assert np.round(cflux[0], decimals=3) == 2*4.471, 'Bad correction'

    ivar = np.full(wave.size, 0.5, dtype=float)
    cflux, civar = atmext.correct(2*flux, corr, ivar=ivar)
    assert np.round(cflux[0], decimals=3) == 2*4.471, 'Bad correction'
    assert np.round(civar[0], decimals=3) == 0.025, 'Bad correction'

    # Test for failures for shape mismatches
    with pytest.raises(PypeItError):
        cflux = atmext.correct(flux, corr[1:])
    with pytest.raises(PypeItError):
        cflux = atmext.correct(flux[1:], corr)
    with pytest.raises(PypeItError):
        cflux = atmext.correct(flux, corr, ivar=ivar[1:])

    # Try 2D
    wave = wave.reshape(2,-1)
    flux = flux.reshape(2,-1)
    ivar = ivar.reshape(2,-1)
    corr = atmext.correction_factor(wave, airmass=1.5)
    assert corr.ndim == 2, 'Wrong dimensionality for correction factor'

    cflux, civar = atmext.correct(2*flux, corr, ivar=ivar)
    assert np.round(cflux[0,0], decimals=3) == 2*4.471, 'Bad correction'
    assert np.round(civar[0,0], decimals=3) == 0.025, 'Bad correction'


