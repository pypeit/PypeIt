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

mjd = 57783.269661
RA = '07:06:23.45'
DEC = '+30:20:50.5'
hdr_equ = 2000.
lon = 155.47833            # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
lat = 19.82833             # Latitude of the telescope
alt = 4160.0               # Elevation of the telescope (in m)


@pytest.fixture
def fitstbl():
    return dummy_fitstbl()


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
    #    vrotate = -0.25490532
    assert np.isclose(corrhelio, -12.49764005490221, rtol=1e-5)
    assert np.isclose(corrbary, -12.510015817405023, rtol=1e-5)


def test_geocorrect(fitstbl):
    """
    """

    # Specobj (wrap in a list to mimic a slit)
    scidx = 5
    obstime = Time(fitstbl['mjd'][scidx], format='mjd')#'%Y-%m-%dT%H:%M:%S.%f')
    radec = ltu.radec_to_coord((fitstbl["ra"][scidx], fitstbl["dec"][scidx]))

    helio, hel_corr = wave.geomotion_correct(radec, obstime, lon, lat, alt, 'heliocentric')
    assert np.isclose(helio, -9.17461338, rtol=1e-5)  # Checked against x_keckhelio
    #assert np.isclose(helio, -9.3344957, rtol=1e-5)  # Original
    assert np.isclose(1-hel_corr, 3.060273748e-05, rtol=1e-5)

    # Now apply to a specobj
    npix = 1000
    sobj = specobj.SpecObj('MultiSlit', 1, SLITID=0)
    sobj.BOX_WAVE = np.linspace(4000., 6000., npix)
    sobj.BOX_COUNTS = 50.*(sobj.BOX_WAVE/5000.)**-1.
    sobj.BOX_COUNTS_IVAR = 1./sobj.BOX_COUNTS.copy()
    sobj.apply_helio(hel_corr, 'heliocentric')
    assert np.isclose(sobj.BOX_WAVE[0], 3999.877589008, rtol=1e-8)
