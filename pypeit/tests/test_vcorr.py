"""
Module to run tests on arvcorr
"""
import numpy as np
import pytest

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units

from linetools import utils as ltu

from pypeit.core import wave
from pypeit import specobj
from pypeit import specobjs
from pypeit.tests.tstutils import dummy_fitstbl
from pypeit import telescopes

mjd = 57783.269661
RA = '07:06:23.45'
DEC = '+30:20:50.5'
hdr_equ = 2000.
keck = telescopes.KeckTelescopePar()
lon = keck['longitude']
lat = keck['latitude']
alt = keck['elevation']


def test_geovelocity():
    """ Test the full geomotion velocity calculation
    """
    loc = (lon * units.deg, lat * units.deg, alt * units.m,)

    radec = SkyCoord(RA, DEC, unit=(units.hourangle, units.deg), frame='icrs')
    obstime = Time(mjd, format='mjd', scale='utc', location=loc)

    corrhelio = wave.geomotion_velocity(obstime, radec, frame="heliocentric")
    corrbary = wave.geomotion_velocity(obstime, radec, frame="barycentric")

    # IDL
    # vhel = x_keckhelio(106.59770833333332, 30.34736111111111, 2000., jd=2457783.769661)
    # print, vhel ===> 12.621846
    #    vrotate = -0.25490532
    assert np.isclose(corrhelio, -12.654140253363275, rtol=1e-5)
    assert np.isclose(corrbary, -12.666516016238132, rtol=1e-5)


def test_geocorrect():
    """
    """

    fitstbl = dummy_fitstbl()
    # Specobj (wrap in a list to mimic a slit)
    scidx = 5
    obstime = Time(fitstbl['mjd'][scidx], format='mjd')#'%Y-%m-%dT%H:%M:%S.%f')
    radec = ltu.radec_to_coord((fitstbl["ra"][scidx], fitstbl["dec"][scidx]))

    helio, hel_corr = wave.geomotion_correct(radec, obstime, lon, lat, alt, 'heliocentric')
    # IDL
    # vhel = x_keckhelio(106.59770833333332, 30.34736111111111, 2000., jd=2457045.5036)
    # print, vhel ===> 8.8838761
    assert np.isclose(helio, -8.873666027875592, rtol=1e-5)  # Checked against x_keckhelio
    assert np.isclose(1-hel_corr, 2.959892574860845e-05, rtol=1e-5)

    # Now apply to a specobj
    npix = 1000
    sobj = specobj.SpecObj('MultiSlit', 'DET01', SLITID=0)
    sobj.BOX_WAVE = np.linspace(4000., 6000., npix)
    sobj.BOX_COUNTS = 50.*(sobj.BOX_WAVE/5000.)**-1.
    sobj.BOX_COUNTS_IVAR = 1./sobj.BOX_COUNTS.copy()
    sobj.apply_helio(hel_corr, 'heliocentric')
    assert np.isclose(sobj.BOX_WAVE[0], 3999.8816042970057, rtol=1e-8)


