# -*- coding: utf-8 -*-
"""
Data utilities for built-in PypeIt data files


NOTE: The remote package data only changes with a new release of PypeIt, and
      the installation of a new version from PyPI or conda will obliterate
      the existing package data.  Therefore, the latest versions of a file
      will be downloaded following an upgrade.

      There is the issue, however, of an older version of PypeIt downloading a
      file associated with a newer version of the pipeline.  Perhaps on the
      remote server, the package data files can be sorted in subdirectories
      by PypeIt version number, so the proper files are always downloaded for
      the pipeline version in use.

.. include:: ../include/links.rst
"""
import os
from pkg_resources import resource_filename

from astropy.utils import data as astropy_data
from linetools.spectra import xspectrum1d

from pypeit import io
from pypeit import msgs

__all__ = ['Paths', 'load_telluric_grid', 'load_thar_spec',
           'load_sky_spectrum', 'get_reid_arxiv_filepath']


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


# Remote-fetch functions for package data not distributed via PyPI ===========#
def get_reid_arxiv_filepath(arxiv_file, use_local=False):
    """get_reid_arxiv_filepath Return the full path to the `reid_arxiv` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the `reid_arxiv` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the `reid_arxiv`
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the `redi_arxiv` file, and
    downloads it from the remote server into the proper `Paths` location for
    future use.  As most users will need only a small number of `reid_arxiv`
    files for thier particular reductions, the remote fetch will only occur
    once per file (per version on PypeIt installed via PyPI or conda).

    Args:
        arxiv_file: str
          The base filename of the `reid_arxiv` file to be located
        use_local: bool, optional
          [STUB FOR FUTURE FUNCTIONALITY]  If the WavelengthSolutionPar
          parameter `use_local` is set to True, look for the `reid_arxiv`
          file on the local filesystem rather than in the PypeIt package
          data.

    Returns:
        calibfile: str
          The full path to the `reid_arxiv` file
    """
    # Full path within the package data structure:
    reid_path = os.path.join(Paths.reid_arxiv, arxiv_file)

    # Check if the file does not exist in the package directory
    if not os.path.isfile(reid_path):
        msgs.warn(f"reid_arxiv file {arxiv_file} does not exist in the package directory.")

        github_url = f"https://github.com/pypeit/PypeIt/blob/release/pypeit/data/arc_lines/reid_arxiv/{arxiv_file}?raw=true"

        reid_path = astropy_data.download_file(github_url, cache=True, timeout=10, pkgname='pypeit')

    # Return current functionality for now
    return reid_path


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



def main():
    pass
if __name__ == '__main__':
    main()