from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypit import arparse as settings
from pypit import armsgs
from astropy.time import Time
from astropy.coordinates import SkyCoord, solar_system, EarthLocation, ICRS, UnitSphericalRepresentation, CartesianRepresentation
from astropy import units as u

# Logging
msgs = armsgs.get_logger()


def calculate(slf, fitsdict, idx):
    """
    Correct the wavelength calibration solution to the desired reference frame
    """

    frame = settings.argflag["reduce"]["calibrate"]["refframe"]
    mjd = float(fitsdict["time"][idx])
    lat = settings.spect['mosaic']['latitude']
    lon = settings.spect['mosaic']['longitude']
    alt = settings.spect['mosaic']['elevation']
    loc = (lon * u.deg, lat * u.deg, alt * u.m,)

    radec = SkyCoord(slf._fitsdict["ra"][idx], slf._fitsdict["dec"][idx], unit=(u.hourangle, u.deg), frame='fk5')
    obstime = Time(mjd, format='mjd', scale='utc', location=loc)

    vcorr = correction(obstime, radec, frame=frame)

    msgs.info("{0:s} velocity correction = {1:+.4f} km/s for file:".format(frame.title(), vcorr) + msgs.newline() +
              fitsdict["filename"][idx])
    w = np.where(slf._waveids == -999999.9)
    slf._waveids += slf._waveids*vcorr/299792.458
    slf._waveids[w] = -999999.9
    return slf._waveids


def correction(time, skycoord, frame="heliocentric"):
    """ Perform a barycentric/heliocentric velocity correction.

    For the correciton, this routine uses the ephemeris:  astropy.coordinates.solar_system_ephemeris.set
    For more information see `~astropy.coordinates.solar_system_ephemeris`.

    Parameters
    ----------
    time : astropy.time.Time
      The time of observation, including the location.
    skycoord: astropy.coordinates.SkyCoord
      The RA and DEC of the pointing, as a SkyCoord quantity.

    Returns
    -------
    vcorr : float
      The velocity correction that should be added to the original velocity.
    """

    # Check that the RA/DEC of the object is ICRS compatible
    if not skycoord.is_transformable_to(ICRS()):
        msgs.error("Cannot transform RA/DEC of object to the ICRS")

    # Calculate ICRS position and velocity of Earth's geocenter
    ep, ev = solar_system.get_body_barycentric_posvel('earth', time)
    # Calculate GCRS position and velocity of observatory
    op, ov = time.location.get_gcrs_posvel(time)
    # ICRS and GCRS are axes-aligned. Can add the velocities
    velocity = ev + ov
    if frame == "heliocentric":
        # ICRS position and velocity of the Sun
        sp, sv = solar_system.get_body_barycentric_posvel('sun', time)
        velocity += sv

    # Get unit ICRS vector in direction of SkyCoord
    sc_cartesian = skycoord.icrs.represent_as(UnitSphericalRepresentation).represent_as(CartesianRepresentation)
    return sc_cartesian.dot(velocity).to(u.km / u.s).value
