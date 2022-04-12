# -*- coding: utf-8 -*-
"""
Data locations for built-in PypeIt data files

.. include:: ../include/links.rst
"""
import os
from pkg_resources import resource_filename

from linetools.spectra import xspectrum1d

from pypeit import io
from pypeit import msgs

__all__ = ['Paths', 'load_telluric_grid', 'load_thar_spec',
           'load_sky_spectrum']


class Paths_meta(type):
    """Paths_meta MetaClass for Paths; only needed until python>=3.9

    The use of this metaclass is necessary until PypeIt sets python>=3.9, at
    which time, the methods shown below can be in the base Paths() class with
    the dual decorators:
        @classmethod
        @property

    This entire machinery is in place because we do not instantiate the Paths()
    class, but rather use it only as a container for the hardwired paths with
    error checking.

    TODO: Upon upgrade to python>=3.9, turn the metaclass into the Paths()
          class, and add the @classmethod decorator to each method.
    """
    def __init__(cls, *args, **kwargs):

        # Class Attributes -- Hardwired Paths
        cls._data = resource_filename('pypeit', 'data')

        # Telluric Corrections
        cls._telgrid = os.path.join(cls._data, 'telluric', 'atm_grids')
        cls._tel_model = os.path.join(cls._data, 'telluric', 'models')

        # Wavelength Calibrations
        cls._arclines = os.path.join(cls._data, 'arc_lines')
        cls._reid_arxiv = os.path.join(cls._data, 'arc_lines', 'reid_arxiv')
        cls._linelist = os.path.join(cls._data, 'arc_lines', 'lists')
        cls._nist = os.path.join(cls._data, 'arc_lines', 'NIST')
        cls._arc_plot = os.path.join(cls._data, 'arc_lines', 'plots')

        # Flux Calibrations
        cls._standards = os.path.join(cls._data, 'standards')
        cls._extinction = os.path.join(cls._data, 'extinction')
        cls._skisim = os.path.join(cls._data, 'skisim')
        cls._filters = os.path.join(cls._data, 'filters')
        cls._sensfuncs = os.path.join(cls._data, 'sensfuncs')

        # Other
        cls._sky_spec = os.path.join(cls._data, 'sky_spec')
        cls._static_calibs = os.path.join(cls._data, 'static_calibs')

    @property
    def data(cls):
        return check_isdir(cls._data)

    # Telluric Corrections
    @property
    def telgrid(cls):
        return check_isdir(cls._telgrid)
    @property
    def tel_model(cls):
        return check_isdir(cls._tel_model)

    # Wavelength Calibrations
    @property
    def arclines(cls):
        return check_isdir(cls._arclines)
    @property
    def reid_arxiv(cls):
        return check_isdir(cls._reid_arxiv)
    @property
    def linelist(cls):
        return check_isdir(cls._linelist)
    @property
    def nist(cls):
        return check_isdir(cls._nist)
    @property
    def arc_plot(cls):
        return check_isdir(cls._arc_plot)

    # Flux Calibrations
    @property
    def standards(cls):
        return check_isdir(cls._standards)
    @property
    def extinction(cls):
        return check_isdir(cls._extinction)
    @property
    def skisim(cls):
        return check_isdir(cls._skisim)
    @property
    def filters(cls):
        return check_isdir(cls._filters)
    @property
    def sensfuncs(cls):
        return check_isdir(cls._sensfuncs)

    # Other
    @property
    def sky_spec(cls):
        return check_isdir(cls._sky_spec)
    @property
    def static_calibs(cls):
        return check_isdir(cls._static_calibs)


class Paths(metaclass=Paths_meta):
    """Paths List of hardwired paths within the pypeit.data module

    [extended_summary]
    """


# Loading Functions for Particular File Types ================================#
def load_telluric_grid(filename):
    """
    Load a telluric atmospheric grid

    NOTE: This is where the path to the data directory is added!

    Args:
        filename: str
          The filename (NO PATH) of the telluric atmospheric grid to use.

    Returns:
        telgrid: `astropy.io.fits.HDUList`
          Telluric Grid FITS HDU list
    """
    # Check for existance of file parameter
    if not filename:
        msgs.error("No file specified for telluric correction.  "
                   "See https://pypeit.readthedocs.io/en/latest/telluric.html")

    # Add standard data path to the filename, as contained in default pypeit pars
    file_with_path = os.path.join(Paths.telgrid, filename)

    # Check for existance of file
    if not os.path.isfile(file_with_path):
        msgs.error(f"File {file_with_path} is not on your disk.  "
                   "You likely need to download the Telluric files.  "
                   "See https://pypeit.readthedocs.io/en/release/installing.html#atmospheric-model-grids")

    return io.fits_open(file_with_path)


def load_thar_spec():
    """
    Load the archived ThAr spectrum

    NOTE: This is where the path to the data directory is added!

    Args:
        filename: str
          The filename (NO PATH) of the telluric atmospheric grid to use.

    Returns:
        thar_spec: `astropy.io.fits.HDUList`
          ThAr Spectrum FITS HDU list
    """
    return io.fits_open(os.path.join(Paths.arclines, 'thar_spec_MM201006.fits'))


def load_sky_spectrum(sky_file):
    """
    Load a sky spectrum into an XSpectrum1D object

    NOTE: This is where the path to the data directory is added!

    .. todo::

        Try to eliminate the XSpectrum1D dependancy

    Args:
        sky_file: str
          The filename (NO PATH) of the sky file to use.

    Returns:
        sky_spec: XSpectrum1D
          spectrum
    """
    return xspectrum1d.XSpectrum1D.from_file(os.path.join(Paths.sky_spec, sky_file))


# Utility Function ===========================================================#
def check_isdir(path):
    """_check_isdir Check that the hardwired directory exists

    If yes, return the directory path, else raise NotADirectoryError
    """
    if not os.path.isdir(path):
        raise NotADirectoryError(f"Unable to find {path}.  "
                                    "Check your installation.")
    return path
