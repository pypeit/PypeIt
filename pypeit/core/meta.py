"""
Provides methods common to :class:`pypeit.metadata.PypeItMetaData` and
:class:`pypeit.spectographs.spectrograph.Spectrograph` that define the
common metadata used for all spectrographs.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import numpy as np

from astropy import units, coordinates

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
                       'dither': dict(dtype=float, comment='Dither amount in arcsec'),
                       'dithpat': dict(dtype=str, comment='Dither pattern'),
                       'dithpos': dict(dtype=str, comment='Dither position'),
                       'dithoff': dict(dtype=float, comment='Dither offset'),
                       'echangle':dict(dtype=float, comment='Echelle angle'),
                       'filter1': dict(dtype=str, comment='First filter in optical path'),
                       'frameno': dict(dtype=str, comment='Frame number provided by instrument software'),
                       'hatch': dict(dtype=str, comment='Position of instrument hatch'),
                       'humidity': dict(dtype=float, comment='Relative humidity (0 to 1) at observation time'),
                       'idname': dict(dtype=str, comment='Instrument supplied frametype (e.g. bias)'),
                       'instrument': dict(dtype=str, comment='Header supplied instrument name'),
                       'mode': dict(dtype=str, comment='Observing mode'),
                       'object': dict(dtype=str, comment='Alternative object name (cf. target)'),
                       'obstime': dict(dtype=str, comment='Observation time'),
                       'oscansec': dict(dtype=str, comment='Overscan section (windowing)'),
                       'pressure': dict(dtype=float, comment='Pressure (units.bar) at observation time'),
                       'seq_expno': dict(dtype=int, comment='Number of exposure in observing sequence'),
                       'slitwid': dict(dtype=float, comment='Slit width, sometimes distinct from decker'),
                       'slitlength': dict(dtype=float, comment='Slit length, used only for long slits'),
                       'temperature': dict(dtype=float, comment='Temperature (units.K) at observation time'),
                       'utc': dict(dtype=str, comment='UTC of observation'),
                       'xd': dict(dtype=float, comment='Cross disperser (e.g. red or blue for HIRES)'),
                       'xdangle':dict(dtype=float, comment='Cross disperser angle'),
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

