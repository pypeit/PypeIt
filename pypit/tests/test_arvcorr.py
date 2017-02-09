# Module to run tests on arvcorr

import numpy as np
import pytest

from pypit import pyputils
msgs = pyputils.get_dummy_logger()

from pypit import arvcorr as py_arvcorr
from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy import units as u

mjd = 57783.269661
RA = '07:06:23.45'
DEC = '+30:20:50.5'
hdr_equ = 2000.
lon = 155.47833            # Longitude of the telescope (NOTE: West should correspond to positive longitudes)
lat = 19.82833              # Latitude of the telescope
alt = 4160.0               # Elevation of the telescope (in m)


def test_vhelio():
    """ Test the full helio-centric calculation
    """
    loc = (lon * u.deg, lat * u.deg, alt * u.m,)

    radec = SkyCoord(RA, DEC, unit=(u.hourangle, u.deg), frame='icrs')
    obstime = Time(mjd, format='mjd', scale='utc', location=loc)

    corrhelio = py_arvcorr.correction(obstime, radec, frame="heliocentric")
    corrbary  = py_arvcorr.correction(obstime, radec, frame="barycentric")

    # IDL
    # vhel = x_keckhelio(106.59770833333332, 30.34736111111111, 2000., jd=2457783.769661)
    #    vrotate = -0.25490532
    pytest.set_trace()
    assert np.isclose(corrhelio, -12.49764005490221)
    assert np.isclose(corrbary, -12.510015817405023)

