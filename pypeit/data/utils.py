# -*- coding: utf-8 -*-
"""
Data utilities for built-in PypeIt data files

.. include:: ../include/links.rst
"""
import os
import shutil
import urllib

from astropy.utils import data as astropy_data
from linetools.spectra import xspectrum1d
from pkg_resources import resource_filename

from pypeit import io
from pypeit import msgs
from pypeit import __version__ as pypeit_version

__all__ = ['Paths', 'load_telluric_grid', 'load_thar_spec',
           'load_sky_spectrum', 'get_reid_arxiv_filepath',
           'get_skisim_filepath', 'get_sensfunc_filepath',
           'fetch_remote_file', 'write_file_to_cache']


# Package-Data Paths =========================================================#
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
def get_reid_arxiv_filepath(arxiv_file):
    """get_reid_arxiv_filepath Return the full path to the `reid_arxiv` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the `reid_arxiv` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the `reid_arxiv`
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the `reid_arxiv` file, and
    downloads it from the remote server using AstroPy's `download_file`
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally `~/.pypeit/cache`) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the `reid_arxiv` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of `reid_arxiv` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        arxiv_file: str
          The base filename of the `reid_arxiv` file to be located

    Returns:
        calibfile: str
          The full path to the `reid_arxiv` file
    """
    # Full path within the package data structure:
    reid_path = os.path.join(Paths.reid_arxiv, arxiv_file)

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not os.path.isfile(reid_path):

        # Output a brief warning for now -- makes it easier to find in the output
        msgs.warn(f"reid_arxiv file {arxiv_file} does not exist in the package directory.")

        reid_path = fetch_remote_file(arxiv_file, "arc_lines/reid_arxiv")

        # If a development version, copy into the package directory, point path there
        if ".dev" in pypeit_version:
            shutil.copy(reid_path, os.path.join(Paths.reid_arxiv, arxiv_file))
            reid_path = os.path.join(Paths.reid_arxiv, arxiv_file)

    # Return the path to the `reid_arxiv` file
    return reid_path


def get_skisim_filepath(skisim_file):
    """get_skisim_filepath Return the full path to the `skisim` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the `skisim` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the `skisim`
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the `skisim` file, and
    downloads it from the remote server using AstroPy's `download_file`
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally `~/.pypeit/cache`) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the `skisim` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of `skisim` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        skisim_file: str
          The base filename of the `skisim` file to be located

    Returns:
        calibfile: str
          The full path to the `skisim` file
    """
    # Full path within the package data structure:
    skisim_path = os.path.join(Paths.skisim, skisim_file)

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not os.path.isfile(skisim_path):

        # Output a brief warning for now -- makes it easier to find in the output
        msgs.warn(f"skisim file {skisim_file} does not exist in the package directory.")

        skisim_path = fetch_remote_file(skisim_file, "skisim")

        # If a development version, copy into the package directory, point path there
        if ".dev" in pypeit_version:
            shutil.copy(skisim_path, os.path.join(Paths.skisim, skisim_file))
            skisim_path = os.path.join(Paths.skisim, skisim_file)

    # Return the path to the `skisim` file
    return skisim_path


def get_sensfunc_filepath(sensfunc_file):
    """get_sensfunc_filepath Return the full path to the `sensfunc` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the `sensfunc` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the `sensfunc`
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the `sensfunc` file, and
    downloads it from the remote server using AstroPy's `download_file`
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally `~/.pypeit/cache`) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the `sensfunc` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of `sensfunc` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        sensfunc_file: str
          The base filename of the `sensfunc` file to be located

    Returns:
        calibfile: str
          The full path to the `sensfunc` file
    """
    # Full path within the package data structure:
    sensfunc_path = os.path.join(Paths.sensfuncs, sensfunc_file)

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not os.path.isfile(sensfunc_path):

        # Output a brief warning for now -- makes it easier to find in the output
        msgs.warn(f"sensfunc file {sensfunc_file} does not exist in the package directory.")

        sensfunc_path = fetch_remote_file(sensfunc_file, "sensfuncs")

        # If a development version, copy into the package directory, point path there
        if ".dev" in pypeit_version:
            shutil.copy(sensfunc_path, os.path.join(Paths.sensfuncs, sensfunc_file))
            sensfunc_path = os.path.join(Paths.sensfuncs, sensfunc_file)

    # Return the path to the `sensfunc` file
    return sensfunc_path


def get_telgrid_filepath(telgrid_file):
    """get_sensfunc_filepath Return the full path to the `sensfunc` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the `sensfunc` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the `sensfunc`
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the `sensfunc` file, and
    downloads it from the remote server using AstroPy's `download_file`
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally `~/.pypeit/cache`) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the `sensfunc` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of `sensfunc` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        sensfunc_file: str
          The base filename of the `sensfunc` file to be located

    Returns:
        calibfile: str
          The full path to the `sensfunc` file
    """
    # Full path within the package data structure:
    telgrid_path = os.path.join(Paths.telgrid, telgrid_file)

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for most installations
    if not os.path.isfile(telgrid_path):

        # Output a brief warning for now -- makes it easier to find in the output
        msgs.warn(f"telgrid file {telgrid_file} does not exist in the package directory.")

        telgrid_path = fetch_remote_file(telgrid_file, 'telgrid', remote_host='s3_cloud')

        # If a development version, MOVE into the package directory, point path there
        if ".dev" in pypeit_version:
            shutil.move(telgrid_path, os.path.join(Paths.telgrid, telgrid_file))
            telgrid_path = os.path.join(Paths.telgrid, telgrid_file)

    # Return the path to the `telgrid` file
    return telgrid_path


def fetch_remote_file(filename, filetype, remote_host='github', install_script=False,
                      force_update=False):
    """fetch_remote_file Use `astropy.utils.data` to fetch file from remote or cache

    The function `download_file` will first look in the local cache (the option
    `cache=True` is used with this function to retrieve downloaded files from
    the cache, as needed) before downloading the file from the remote server.

    Args:
        filename: str
          The base filename to search for
        filetype: str
          The subdirectory of `pypeit/data/` in which to find the file
          (e.g., `arc_lines/reid_arxiv` or `sensfuncs`)
        remote_host: str, optional
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  [Default: 'github']
        install_script: bool, optional
          This function is being called from an install script (i.e.,
          `pypeit_install_telluric`) -- relates to warnings displayed.
          [Default: False]
        force_update: bool, optional
          Force `astropy_data.download_file()` to update the cache by downloading
          the latest version.  [Default: False]

    Returns:
        path_to_file: str
          The local path to the desired file in the cache
    """
    remote_url = _build_remote_url(filename, filetype, remote_host=remote_host)

    if remote_host == "s3_cloud" and not install_script:
        # Display a warning that this may take a while, and the user
        #   may wish to download using the `pypeit_install_telluric` script
        msgs.warn(f"You may wish to download {filename}{msgs.newline()}"
                    f"independently from your reduction by using the{msgs.newline()}"
                    "`pypeit_install_telluric` script.")

    # Get the file from cache, if available, or download from the remote server
    cache = 'update' if force_update else True
    try:
        return astropy_data.download_file(remote_url, cache=cache, timeout=10, pkgname="pypeit")
    except urllib.error.URLError as error:
        msgs.error(f"Error downloading {filename}: {error}")


def write_file_to_cache(filename, cachename, filetype, remote_host="github"):
    """write_file_to_cache Use `astropy.utils.data` to save local file to cache

    This function writes a local file to the PypeIt cache as if it came from a
    remote server.  This is useful for being able to use lcoally created files
    in place of PypeIt-distributed versions.

    Args:
        filename: str
          The filename of the local file to save
        cachename: str
          The name of the cached version of the file
        filetype: str
          The subdirectory of `pypeit/data/` in which to find the file
          (e.g., `arc_lines/reid_arxiv` or `sensfuncs`)
        remote_host: str, optional
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  [Default: 'github']
    """
   # Build the `url_key` as if this file were in the remote location
    url_key = _build_remote_url(cachename, filetype, remote_host=remote_host)

    # Use `import file_to_cache()` to place the `filename` into the cache
    astropy_data.import_file_to_cache(url_key, filename, pkgname="pypeit")


def _build_remote_url(f_name, f_type, remote_host=""):
    """build_remote_url Build the remote URL for the `astropy.utils.data` functions

    This function keeps the URL-creation in one place.  In the event that files
    are moved from GitHub or S3_Cloud, this is the only place that would need
    to be changed.

    Args:
        f_name: str
          The base filename to search for
        f_type: str
          The subdirectory of `pypeit/data/` in which to find the file
          (e.g., `arc_lines/reid_arxiv` or `sensfuncs`)
        remote_host: str, optional
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  [Default: '']

    Returns:
        url: str
          The URL of the `f_name` of `f_type` on server `remote_host`
    """
    if remote_host == "github":
        # Build up the remote_url for GitHub
        # Look in the current `develop` branch if the code is not a tagged release
        tag = "develop" if ".dev" in pypeit_version else pypeit_version
        # TODO: If we host these files elsewhere, need to change this hard-code
        return (f"https://github.com/pypeit/PypeIt/blob/{tag}/pypeit/"
                   f"data/{f_type}/{f_name}?raw=true")

    if remote_host == "s3_cloud":
        # Build up the remote_url for S3 Cloud
        # TODO: Put the correct path here once we get it from @profxj
        return f"https://s3/{f_type}/{f_name}"

    msgs.error(f"Remote host type {remote_host} is not supported for package data caching.")


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

    # Get the data path for the filename, whether in the package directory or cache
    file_with_path = get_telgrid_filepath(filename)

    # Check for existance of file
    # NOTE: With the use of `get_telgrid_filepath()`, this should never run
    if not os.path.isfile(file_with_path):
        msgs.error(f"File {file_with_path} is not on your disk.  "
                   "You likely need to download the Telluric files.  "
                   "See https://pypeit.readthedocs.io/en/release/installing.html"
                   "#atmospheric-model-grids")

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
