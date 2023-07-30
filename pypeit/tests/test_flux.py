"""
Module to run tests on a few flux routines
"""
import os

import numpy as np
import pytest

from astropy import units
from astropy.coordinates import SkyCoord

from pypeit.core import flux_calib
from pypeit import telescopes
from pypeit.par.pypeitpar import Coadd1DPar


from pypeit.pypmsgs import PypeItError



def test_blackbody():
    a, teff = 2.65, 10086  # Parameter of J1245+4238
    wave, flam = flux_calib.blackbody_func(a, teff)
    flam_scl = flam*flux_calib.BB_SCALE_FACTOR  # In units 10^-17 erg/s/cm2/A
    res = np.interp(4000.0, wave, flam_scl)
    # The following value is close to the value shown in the Figure 15 of Suzuki & Fukugita (2018).
    assert(np.isclose(res, 89.6419630016348))


def test_find_standard():
    # G191b2b
    coord = SkyCoord('J050630.6+524951.0', unit=(units.hourangle, units.deg))  #
    # Grab
    std_dict = flux_calib.find_standard_file(coord.ra.value, coord.dec.value)
    # Test
    assert std_dict['name'] == 'G191B2B'
    assert os.path.split(std_dict['cal_file'])[1] == 'g191b2b_stisnic_002.fits'
    assert std_dict['std_source'] == 'calspec'
    # Fail to find
    # near G191b2b
    coord = SkyCoord('J050630.6+522201.0', unit=(units.hourangle, units.deg))  #
    with pytest.raises(PypeItError):
        std_dict = flux_calib.find_standard_file(coord.ra.value, coord.dec.value)


def test_load_extinction():
    # Load
    mtham = telescopes.ShaneTelescopePar()
    lon = mtham['longitude']
    lat = mtham['latitude']
    extinct = flux_calib.load_extinction_data(lon, lat, 'closest')
    np.testing.assert_allclose(extinct['wave'][0], 3200.)
    assert extinct['wave'].unit == units.AA
    np.testing.assert_allclose(extinct['mag_ext'][0], 1.084)

    # Fail on location in the western Mediterranean Sea
    with pytest.raises(PypeItError):
        extinct = flux_calib.load_extinction_data(0., 37.3413889, 'closest')

    # Fail on nonexistant extinction filename
    with pytest.raises(PypeItError):
        extinct = flux_calib.load_extinction_data(lon, lat, 'northpoleextinct.dat')


def test_extinction_correction():
    # Load
    mtham = telescopes.ShaneTelescopePar()
    lon = mtham['longitude']
    lat = mtham['latitude']
    extinct = flux_calib.load_extinction_data(lon, lat, 'closest')
    # Correction
    wave = np.arange(3000.,10000.)*units.AA
    AM=1.5
    flux_corr = flux_calib.extinction_correction(wave, AM, extinct)
    # Test
    np.testing.assert_allclose(flux_corr[0], 4.47095192)


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
