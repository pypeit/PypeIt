"""
Module to run tests on a few flux routines
"""

import numpy as np

from pypeit.core import flux_calib
from pypeit.par.pypeitpar import Coadd1DPar


def test_filter_scale():
    # Test scale_in_filter() method which is called in coadding
    wave = np.arange(3000.,10000.)
    flux = np.ones_like(wave)
    gdm = np.ones_like(wave, dtype=bool)
    #
    par = Coadd1DPar()
    par['filter'] = 'DECAM-R'
    par['filter_mag'] = 17.
    # Run
    scale = flux_calib.scale_in_filter(wave, flux, gdm, par)
    assert np.isclose(scale, 41.698475048180406, rtol=1e-3)


def test_fnu_to_flam():
    # f_lambda = f_nu c / lambda^2, in units of 1e-17 erg/s/cm^2/Angstrom
    assert flux_calib.is_flux_density('MJy'), 'MJy is a flux density per unit frequency.'
    for unit in ['MJy/sr', 'e-', 'ADU', 'unknown']:
        assert not flux_calib.is_flux_density(unit), f'{unit} is not a flux density.'
    wave = np.array([0., 5000., 30000.])
    factor = flux_calib.fnu_to_flam(wave, 'Jy')
    expected = 1e-23 * 2.99792458e18 / wave[1:]**2 / 1e-17
    assert np.allclose(factor[1:], expected), 'Wrong conversion from Jy to f_lambda.'
    assert factor[0] == 0., 'The conversion should be 0 where the wavelength is not positive.'
    assert np.allclose(flux_calib.fnu_to_flam(wave, 'MJy')[1:], 1e6 * expected), \
        'Wrong conversion from MJy to f_lambda.'
