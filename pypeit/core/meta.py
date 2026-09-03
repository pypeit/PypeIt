"""
Provides methods common to :class:`~pypeit.metadata.PypeItMetaData` and
:class:`~pypeit.spectrographs.spectrograph.Spectrograph` that define the
common metadata used for all spectrographs.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import numpy as np

from astropy import units, coordinates
from astropy.time import Time

from IPython import embed


def convert_radec(ra, dec):
    """
    Handle multiple ra,dec inputs and return decimal degrees

    If ra, dec are str but do *not* have J or ':' in the RA term,
    then they will be converted to floats

    Args:
        ra (str or float or `numpy.ndarray`_):
            RA as decimal deg (float) or  hh:mm:ss.s (str)
        dec (str or float or `numpy.ndarray`_):
            DEC as decimal deg (float) or  +dd:mm:ss.s (str)
            Must be the same format as ra

    Returns:
        tuple:
           float,float of ra,dec in decimal deg if input is str or float
           np.ndarray, np.ndarray of ra,dec in decimal deg if input is np.ndarray

    """
    if isinstance(ra, str):
        if (('J' in ra) or (':' in ra)) or (' ' in ra.strip()):
            coord = coordinates.SkyCoord(ra, dec, unit=(units.hourangle, units.deg))
            return coord.ra.value, coord.dec.value
        else:
            return float(ra), float(dec)
    elif isinstance(ra, np.ndarray):
        if isinstance(ra[0], str):
            if (('J' in ra[0]) or (':' in ra[0])) or (' ' in ra[0].strip()):
                coords = coordinates.SkyCoord(ra,dec, unit=(units.hourangle, units.deg))
                return coords.ra.value, coords.dec.value
            else:
                return ra.astype(float), dec.astype(float)
        elif isinstance(ra[0], float):
            return ra.astype(float), dec.astype(float)
        else:
            raise IOError("Bad ra, dec format!!")
    else:
        return ra, dec


def airmass(ra, dec, obstime, longitude, latitude, elevation, exptime=None):
    """
    Compute the airmass of an observation from its pointing, time, and the
    observatory location, using `astropy`.

    This is intended for instruments whose raw headers do not record an
    airmass (or elevation) card directly, so it can be reused by any
    spectrograph's ``compound_meta`` method.  It follows the `astropy`
    observation-planning example:
    https://docs.astropy.org/en/latest/coordinates/example_gallery_plot_obs_planning.html

    The airmass returned is the secant of the zenith distance
    (:math:`\\sec z`) computed from the target altitude in the local
    horizontal (alt-az) frame.

    Args:
        ra (:obj:`float`):
            Target right ascension in decimal degrees (ICRS).
        dec (:obj:`float`):
            Target declination in decimal degrees (ICRS).
        obstime (:obj:`float`, :obj:`str`, or `astropy.time.Time`_):
            UT of the observation.  A :obj:`float` is interpreted as an
            MJD, a :obj:`str` is parsed by `astropy.time.Time`_ (e.g. an
            ISOT string), and an `astropy.time.Time`_ is used directly.
        longitude (:obj:`float`):
            East longitude of the observatory in decimal degrees.
        latitude (:obj:`float`):
            Latitude of the observatory in decimal degrees.
        elevation (:obj:`float`):
            Elevation of the observatory in meters.
        exptime (:obj:`float`, optional):
            Exposure time in seconds.  If provided, ``obstime`` is taken to
            be the *start* of the exposure and the airmass is computed at
            the exposure midpoint (``obstime + exptime/2``).  If None, the
            airmass is computed at ``obstime``.

    Returns:
        :obj:`float`: The airmass (:math:`\\sec z`) of the observation.
        Returns `numpy.nan`_ if the target is below the horizon.
    """
    # Observatory location on the Earth
    location = coordinates.EarthLocation.from_geodetic(
        lon=longitude*units.deg, lat=latitude*units.deg,
        height=elevation*units.m)
    # Interpret the observation time: MJD float, parseable string, or Time
    if isinstance(obstime, Time):
        time = obstime
    elif isinstance(obstime, (int, float, np.integer, np.floating)):
        time = Time(obstime, format='mjd')
    else:
        time = Time(obstime)
    # Compute at the midpoint of the exposure, if the exposure time is given
    if exptime is not None:
        time = time + (exptime/2.)*units.s
    # Target coordinate and its altitude/azimuth as seen from the site
    coord = coordinates.SkyCoord(ra*units.deg, dec*units.deg)
    altaz = coord.transform_to(
        coordinates.AltAz(obstime=time, location=location))
    # Below the horizon -> airmass is undefined
    if altaz.alt < 0*units.deg:
        return np.nan
    # AltAz.secz is the airmass (sec of the zenith distance)
    return float(altaz.secz)


def define_core_meta():
    """
    Define the core set of meta data that must be defined
    to run PypeIt.

    See the metadata.rst file for further discussion

    .. warning::

        The keys should all be <= 8 length as they are all written
        to the Header.

    Each meta entry is a dict with the following keys:
       - dtype: str, float, int
       - comment: str
       - rtol: float, optional -- Sets the relative tolerance for
         float meta when used to set a configuration

    Each meta dtype must be scalar or str.  No tuple, list, ndarray, etc.

    Returns:
        dict: core_meta

    """
    # Mainly to format output to PypeIt file
    core_meta = {}

    # Target
    core_meta['ra'] = dict(dtype=float, comment='(J2000) RA in decimal degrees')
    core_meta['dec'] = dict(dtype=float, comment='(J2000) DEC in decimal degrees')
    core_meta['target'] = dict(dtype=str, comment='Name of the target')

    # Instrument related
    core_meta['dispname'] = dict(dtype=str, comment='Disperser name')
    core_meta['decker'] = dict(dtype=str, comment='Slit/mask/decker name')
    core_meta['binning'] = dict(dtype=str, comment='"spectral,spatial" binning')

    # Obs
    core_meta['mjd'] = dict(dtype=float, comment='Observation MJD; Read by astropy.time.Time format=mjd')
    core_meta['airmass'] = dict(dtype=float, comment='Airmass')
    core_meta['exptime'] = dict(dtype=float, comment='Exposure time (s)')

    # Test me
    # TODO: Do we need this?
    for key in core_meta.keys():
        assert len(key) <= 8

    # Return
    return core_meta


def define_additional_meta(nlamps=20):
    """
    Defines meta that tends to be instrument-specific and not used as
    widely in the code.

    See :func:`define_core_meta` for additional details

    For meta used to define configurations, the rtol key specifies
    the relative tolerance for a match

    Args:
        nlamps (:obj:`int`, optional):
            Number of calibrations lamps for this instrument.

    Returns:
        :obj:`dict`: Describes the additional meta data used in
        pypeit.
    """
    # TODO: Can we consolidate dither (only read by vlt_xshooter, and never
    # used in the code) with the MOSFIRE dith* keywords?
    additional_meta = {'amp': dict(dtype=str, comment='Amplifier used'),
                       'arm': dict(dtype=str, comment='Name of arm (e.g. NIR for X-Shooter)'),
                       'calpos': dict(dtype=str, comment='Position of calibration system (KCWI)'),
                       'datasec': dict(dtype=str, comment='Data section (windowing)'),
                       'dateobs': dict(dtype=str, comment='Observation date'),
                       'decker_secondary': dict(dtype=str, comment='Partial Slitmask/decker name. '
                                                                  'It differs from decker. This is currently '
                                                                  'only needed for the reduction of some '
                                                                  'MOSFIRE masks, which use calibrations '
                                                                  'taken with a partially different decker '
                                                                  'name than the one used for the associated '
                                                                  'science frames.'),
                       'detector': dict(dtype=str, comment='Name of detector'),
                       'dichroic': dict(dtype=str, comment='Beam splitter'),
                       'dispangle': dict(dtype=float, comment='Angle of the disperser', rtol=0.),
                       'cenwave': dict(dtype=float, comment='Central wavelength of the disperser', rtol=0.),
                       # TODO what is the difference between dither and dithoff? Also, we should rename these to be
                       # more clearly the offset along the slit to distinguish from the IFU case of RA_off and DEC_off
                       'dither': dict(dtype=float, comment='Dither amount in arcsec'),
                       'dithoff': dict(dtype=float, comment='Dither offset'),
                       'dithpat': dict(dtype=str, comment='Dither pattern'),
                       'dithpos': dict(dtype=str, comment='Dither position'),
                       'posang': dict(dtype=float, comment='Position angle of the observation (degrees, positive is East from North)'),
                       'ra_off': dict(dtype=float, comment='Dither offset in RA'),
                       'dec_off': dict(dtype=float, comment='Dither offset in DEC'),
                       'echangle':dict(dtype=float, comment='Echelle angle'),
                       'filter1': dict(dtype=str, comment='First filter in optical path'),
                       'filter2': dict(dtype=str, comment='Second filter in optical path'),
                       'frameno': dict(dtype=str, comment='Frame number provided by instrument software'),
                       'hatch': dict(dtype=str, comment='Position of instrument hatch'),
                       'humidity': dict(dtype=float, comment='Humidity at observation time (as a percentage, not a fraction)'),
                       'idname': dict(dtype=str, comment='Instrument supplied frametype (e.g. bias)'),
                       'instrument': dict(dtype=str, comment='Header supplied instrument name'),
                       'mode': dict(dtype=str, comment='Observing mode'),
                       'nodpix': dict(dtype=int, comment='Number of rows shuffled when nodding, used for Gemini/GMOS'),
                       'object': dict(dtype=str, comment='Alternative object name (cf. target)'),
                       'obstime': dict(dtype=str, comment='Observation time'),
                       'oscansec': dict(dtype=str, comment='Overscan section (windowing)'),
                       'parangle': dict(dtype=float, comment='Parallactic angle (units.radian)'),
                       'pressure': dict(dtype=float, comment='Pressure (units.pascal) at observation time'),
                       'seq_expno': dict(dtype=int, comment='Number of exposure in observing sequence'),
                       'slitwid': dict(dtype=float, comment='Slit width, sometimes distinct from decker'),
                       'slitlength': dict(dtype=float, comment='Slit length, used only for long slits'),
                       'slitrange': dict(dtype=str, comment='Slit pixel range as colon-separated values, used for MOSFIRE longslit'),
                       'temperature': dict(dtype=float, comment='Temperature (units.K) at observation time'),
                       'utc': dict(dtype=str, comment='UTC of observation'),
                       'mirror': dict(dtype=str, comment='Position of an instrument mirror (e.g. IN or OUT)'),
                       'xd': dict(dtype=float, comment='Cross disperser (e.g. red or blue for HIRES)'),
                       'xdangle':dict(dtype=float, comment='Cross disperser angle'),
                       'readmode': dict(dtype=str,
                                        comment='Read mode of the image'),
                       'savemode': dict(dtype=str,
                                        comment='Save mode of the image'),
                       'dit': dict(dtype=float,
                                   comment='Detector integration time'),
                       'ndit': dict(dtype=int,
                                    comment='Number of integrations'),
                       'camera': dict(dtype=str,
                                      comment='Camera (e.g. N1.8/N3.75/N30 for LUCI)'),
                       'camera_pos':dict(dtype=str, comment='Camera position (e.g. LongRed, ShortBlue, '
                                                            'etc for Gemini GNIRS)'),
                       }

    for kk in range(nlamps):
        additional_meta['lampstat{:02d}'.format(kk+1)] \
                = dict(dtype=str, comment='Status of a given lamp (e.g off/on)')
        additional_meta['lampshst{:02d}'.format(kk + 1)] \
            = dict(dtype=str, comment='Status of a lamp shutter (e.g closed/open)')
    return additional_meta


def get_meta_data_model(nlamps=20):
    """
    Construct full metadata model general to all spectrographs.

    This is a wrapper for :func:`define_core_meta` and
    :func:`define_additional_meta` that checks that the keys defined
    by both are unique (a coding issue) and returns a single combined
    dictionary.

    Args:
        nlamps (:obj:`int`, optional):
            Number of calibrations lamps for this instrument, passed
            directly to :func:`define_additional_meta`.

    Returns:
        :obj:`dict`: Dictionary with the full metadata model common to all
        spectrographs.

    Raises:
        ValueError:
            Raised if the coding of func:`define_core_meta` and
            :func:`define_additional_meta` do not produce unique
            keys. This should never be raised in the released
            version.
    """
    core = define_core_meta()
    add = define_additional_meta(nlamps=nlamps)
    if np.any(np.isin(list(core.keys()), list(add.keys()))):
        raise ValueError('CODING ERROR: Keys in core and additional meta data are not unique!')
    core.update(add)
    return core

