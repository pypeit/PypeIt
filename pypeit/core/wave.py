""" Routines related to flexure, air2vac, etc.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
 """

import numpy as np


from astropy import units
from astropy.coordinates import solar_system, ICRS
from astropy.coordinates import UnitSphericalRepresentation, CartesianRepresentation
from astropy.time import Time


from pypeit import msgs

from IPython import embed

def geomotion_calculate(radec, time, longitude, latitude, elevation, refframe):
    """
    Correct the wavelength calibration solution to the desired reference frame
    """

    # Time
    loc = (longitude * units.deg, latitude * units.deg, elevation * units.m,)
    obstime = Time(time.value, format=time.format, scale='utc', location=loc)
    return geomotion_velocity(obstime, radec, frame=refframe)


def geomotion_correct(radec, time, longitude, latitude, elevation, refframe):
    """
    Correct the wavelength of every pixel to a barycentric/heliocentric frame.

    Args:
        radec (astropy.coordiantes.SkyCoord):
        time (:obj:`astropy.time.Time`):
        gd_slitord (`numpy.ndarray`_):
            Array of good slit/order IDs
        fitstbl : Table/PypeItMetaData
            Containing the properties of every fits file
        longitude (float): deg
        latitude (float): deg
        elevation (float): m
        refframe (str):

    Returns:
        tuple: Two objects are returned:

            - float: The velocity correction that should be applied to
              the wavelength array.
            - float: The relativistic velocity correction that should be
              multiplied by the wavelength array to convert each
              wavelength into the user-specified reference frame.

    """
    # Calculate
    vel = geomotion_calculate(radec, time, longitude, latitude, elevation, refframe)
    vel_corr = np.sqrt((1. + vel/299792.458) / (1. - vel/299792.458))

    # Return
    return vel, vel_corr


def geomotion_velocity(time, skycoord, frame="heliocentric"):
    """ Perform a barycentric/heliocentric velocity correction.

    For the correciton, this routine uses the ephemeris:  astropy.coordinates.solar_system_ephemeris.set
    For more information see `~astropy.coordinates.solar_system_ephemeris`.

    Parameters
    ----------
    time : astropy.time.Time
        The time of observation, including the location.
    skycoord: astropy.coordinates.SkyCoord
        The RA and DEC of the pointing, as a SkyCoord quantity.
    frame : str
        The reference frame that should be used for the calculation.

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
    return sc_cartesian.dot(velocity).to(units.km / units.s).value


def airtovac(wave):
    """ Convert air-based wavelengths to vacuum

    Parameters
    ----------
    wave: Quantity array
        Wavelengths 

    Returns
    -------
    wave: Quantity array
        Wavelength array corrected to vacuum wavelengths
    """
    # Convert to AA
    wave = wave.to(units.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength*factor
    # Units
    new_wave = wavelength*units.AA
    new_wave.to(wave.unit)

    return new_wave


def vactoair(wave):
    """Convert to air-based wavelengths from vacuum

    Parameters
    ----------
    wave: Quantity array
        Wavelengths 

    Returns
    -------
    wave: Quantity array
        Wavelength array corrected to air

    """
    # Convert to AA
    wave = wave.to(units.AA)
    wavelength = wave.value

    # Standard conversion format
    sigma_sq = (1.e4/wavelength)**2. #wavenumber squared
    factor = 1 + (5.792105e-2/(238.0185-sigma_sq)) + (1.67918e-3/(57.362-sigma_sq))
    factor = factor*(wavelength>=2000.) + 1.*(wavelength<2000.) #only modify above 2000A

    # Convert
    wavelength = wavelength/factor
    new_wave = wavelength*units.AA
    new_wave.to(wave.unit)

    return new_wave
