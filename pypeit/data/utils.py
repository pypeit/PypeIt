# -*- coding: utf-8 -*-
"""
Data utilities for built-in ``PypeIt`` data files

.. note::

    If the hostname URL for the telluric atmospheric grids on S3 changes, the
    only place that needs to change is the file ``s3_url.txt``.

----

Implementation Documentation
----------------------------

This module contains the organization scheme for the ``pypeit/data`` files
needed by the ``PypeIt`` package.  Any routine in the package that needs to load
a data file stored in this directory should use the paths supplied by this
module and not call `resource_filename
<https://setuptools.pypa.io/en/latest/pkg_resources.html#resource-extraction>`__
or attempt to otherwise directly access the package directory structure.  In
this way, if structural changes to this directory are needed, only this module
need be modified and the remainder of the package can remain ignorant of those
changes and continue to call the paths supplied by this module.

Furthermore, all paths returned by this module are :obj:`pathlib.Path` objects
rather than pure strings, with all of the functionality therein contained.

Most (by number) of the package data files here are distributed with the
``PypeIt`` package and are accessed via the :class:`~pypeit.data.Paths`
class.  For instance, the NIR spectrophotometry for Vega is accessed via:

.. code-block:: python

    vega_file = data.Paths.standards / 'vega_tspectool_vacuum.dat'

For some directories, however, the size of the included files is large enough
that it was beginning to cause problems with distributing the package via PyPI.
For these specific directories, the data is still stored in the GitHub
repository but is not distributed with the PyPI package.  In order to access
and use these files, we use the AstroPy download/cache system, and specific
functions (``get_*_filepath()``) are required to interact with these files.
Currently, the directories covered by the AstroPy download/cache system are:

    * arc_lines/reid_arxiv
    * skisim
    * sensfuncs

From time to time, it may be necessary to add additional files/directories to
the AstroPy download/cache system.  In this case, there is a particular
sequence of steps required.  The caching routines look for remote-hosted data
files in either the ``develop`` tree or a tagged version tree (e.g., ``1.8.0``)
of the repository, any new files must be already present in the repo before
testing a new ``get_*_filepath()`` routine.  Order of operations is:

    #. Add any new remote-hosted files to the GitHub repo via a separate PR
       that also modifies ``MANIFEST.in`` to exclude these files from the
       distributed package.

    #. Create a new ``get_*_filepath()`` function in this module, following the
       example of one of the existing functions.  Elsewhere in ``PypeIt``, load
       the needed file by invoking the new ``get_*_filepath()`` function.  An
       example of this can be found in ``pypeit/core/flux_calib.py`` where
       ``get_skisim_filepath()`` is called to locate sky transmission files.

If new package-included data are added that are not very large (total directory
size < a few MB), it is not necessary to use the AstroPy cache/download system.
In this case, simply add the directory path to the
:class:`~pypeit.data.Paths` class and access the enclosed files similarly
to the Vega example above.

.. include:: ../include/links.rst
"""
import pathlib
import shutil
import urllib.error

import astropy.utils.data
import github
from linetools.spectra import xspectrum1d
from pkg_resources import resource_filename
import requests

from pypeit import io
from pypeit import msgs
from pypeit import __version__

__all__ = ['Paths', 'load_telluric_grid', 'load_thar_spec',
           'load_sky_spectrum', 'get_reid_arxiv_filepath',
           'get_skisim_filepath', 'get_sensfunc_filepath',
           'get_telgrid_filepath', 'get_linelist_filepath',
           'get_extinctfile_filepath',
           'fetch_remote_file', 'search_cache',
           'write_file_to_cache']


# Package-Data Paths =========================================================#
class Paths:
    """List of hardwired paths within the pypeit.data module

    Each `@property` method returns a :obj:`pathlib.Path` object
    """

    # Class Attributes -- Hardwired Paths
    _data = pathlib.Path(resource_filename('pypeit', 'data'))

    # Telluric Corrections
    _telgrid = _data / 'telluric' / 'atm_grids'
    _tel_model = _data / 'telluric' / 'models'

    # Wavelength Calibrations
    _arclines = _data / 'arc_lines'
    _reid_arxiv = _arclines / 'reid_arxiv'
    _linelist = _arclines / 'lists'
    _nist = _arclines / 'NIST'
    _arc_plot = _arclines /'plots'

    # Flux Calibrations
    _standards = _data / 'standards'
    _extinction = _data / 'extinction'
    _skisim = _data / 'skisim'
    _filters = _data / 'filters'
    _sensfuncs = _data / 'sensfuncs'

    # Other
    _sky_spec = _data / 'sky_spec'
    _static_calibs = _data / 'static_calibs'
    _spectrographs = _data / 'spectrographs'

    @classmethod
    @property
    def data(cls):
        return cls.check_isdir(cls._data)

    # Telluric Corrections
    @classmethod
    @property
    def telgrid(cls):
        return cls.check_isdir(cls._telgrid)
    @classmethod
    @property
    def tel_model(cls):
        return cls.check_isdir(cls._tel_model)

    # Wavelength Calibrations
    @classmethod
    @property
    def arclines(cls):
        return cls.check_isdir(cls._arclines)
    @classmethod
    @property
    def reid_arxiv(cls):
        return cls.check_isdir(cls._reid_arxiv)
    @classmethod
    @property
    def linelist(cls):
        return cls.check_isdir(cls._linelist)
    @classmethod
    @property
    def nist(cls):
        return cls.check_isdir(cls._nist)
    @classmethod
    @property
    def arc_plot(cls):
        return cls.check_isdir(cls._arc_plot)

    # Flux Calibrations
    @classmethod
    @property
    def standards(cls):
        return cls.check_isdir(cls._standards)
    @classmethod
    @property
    def extinction(cls):
        return cls.check_isdir(cls._extinction)
    @classmethod
    @property
    def skisim(cls):
        return cls.check_isdir(cls._skisim)
    @classmethod
    @property
    def filters(cls):
        return cls.check_isdir(cls._filters)
    @classmethod
    @property
    def sensfuncs(cls):
        return cls.check_isdir(cls._sensfuncs)

    # Other
    @classmethod
    @property
    def sky_spec(cls):
        return cls.check_isdir(cls._sky_spec)
    @classmethod
    @property
    def static_calibs(cls):
        return cls.check_isdir(cls._static_calibs)
    @classmethod
    @property
    def spectrographs(cls):
        return cls.check_isdir(cls._spectrographs)

    @staticmethod
    def check_isdir(path):
        """Check that the hardwired directory exists

        If yes, return the directory path, else raise an error message
        """
        if not path.is_dir():
            msgs.error(f"Unable to find {path}.  "
                        "Check your installation.")
        return path


# Remote-fetch functions for package data not distributed via PyPI ===========#
def get_reid_arxiv_filepath(arxiv_file: str) -> tuple[pathlib.Path, str]:
    """Return the full path to the ``reid_arxiv`` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the ``reid_arxiv`` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the ``reid_arxiv``
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the ``reid_arxiv`` file, and
    downloads it from the remote server using AstroPy's ``download_file()``
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally ``~/.pypeit/cache``) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the ``reid_arxiv`` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of ``reid_arxiv`` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        arxiv_file (str):
          The base filename of the ``reid_arxiv`` file to be located

    Returns:
        tuple: The full path and whether the path is in the cache:

           * reid_path (:obj:`pathlib.Path`): The full path to the ``reid_arxiv`` file
           * arxiv_fmt (:obj:`str`): The extension of the ``reid_arxiv`` file (format)
    """
    # Full path within the package data structure:
    reid_path = Paths.reid_arxiv / arxiv_file
    arxiv_fmt = arxiv_file.split(".")[-1].lower()

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not reid_path.is_file():

        # Output an informational message
        msgs.info(f"reid_arxiv file {arxiv_file} does not exist in{msgs.newline()}"
                  "the package directory.  Checking cache or downloading the file now.")

        reid_path = fetch_remote_file(arxiv_file, "arc_lines/reid_arxiv")

    # Return the path to the `reid_arxiv` file, and the file format
    return reid_path, arxiv_fmt


def get_skisim_filepath(skisim_file: str) -> pathlib.Path:
    """Return the full path to the ``skisim`` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the ``skisim`` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the ``skisim``
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the ``skisim`` file, and
    downloads it from the remote server using AstroPy's ``download_file()``
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally ``~/.pypeit/cache``) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the ``skisim`` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of ``skisim`` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        skisim_file (str):
          The base filename of the ``skisim`` file to be located

    Returns:
        :obj:`pathlib.Path`: The full path to the ``skisim`` file
    """
    # Full path within the package data structure:
    skisim_path = Paths.skisim / skisim_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not skisim_path.is_file():

        # Output an informational message
        msgs.info(f"skisim file {skisim_file} does not exist in{msgs.newline()}"
                  "the package directory.  Checking cache or downloading the file now.")

        skisim_path = fetch_remote_file(skisim_file, "skisim")

    # Return the path to the `skisim` file
    return skisim_path


def get_sensfunc_filepath(sensfunc_file: str, symlink_in_pkgdir=False) -> pathlib.Path:
    """Return the full path to the ``sensfunc`` file

    In an attempt to reduce the size of the PypeIt package as distributed on
    PyPI, the ``sensfunc`` files are not longer distributed with the package.
    The collection of files are hosted remotely, and only the ``sensfunc``
    files needed by a particular user are downloaded to the local machine.

    This function checks for the local existance of the ``sensfunc`` file, and
    downloads it from the remote server using AstroPy's ``download_file()``
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally ``~/.pypeit/cache``) and is not placed into the package
    directory itself.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the ``sensfunc`` files will be
    downloaded again (matching the new version #) to catch any changes.

    As most users will need only a small number of ``sensfunc`` files for thier
    particular reductions, the remote fetch will only occur once per file (per
    version of PypeIt).

    Args:
        sensfunc_file (str):
          The base filename of the ``sensfunc`` file to be located
        symlink_in_pkgdir (:obj:`bool`, optional):
          Create a symlink (with the canonical filename) in the package directory
          pointing to the cached downloaded file.  Defaults to False.

    Returns:
        :obj:`pthlib.Path`: The full path to the ``sensfunc`` file
    """
    # Full path within the package data structure:
    sensfunc_path = Paths.sensfuncs / sensfunc_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not sensfunc_path.is_file():

        # Output an informational message
        msgs.info(f"sensfunc file {sensfunc_file} does not exist in{msgs.newline()}"
                  "the package directory.  Checking cache or downloading the file now.")

        sensfunc_path = fetch_remote_file(sensfunc_file, "sensfuncs")

        # If requested, copy to package data directory and point the path there
        if symlink_in_pkgdir:
            path_in_pkgdir = Paths.sensfuncs / sensfunc_file
            # Create the symlink
            path_in_pkgdir.symlink_to(sensfunc_path)
            # Return the path to the symlink in the package directory
            sensfunc_path = path_in_pkgdir

    # Return the path to the `sensfunc` file
    return sensfunc_path


def get_telgrid_filepath(telgrid_file: str) -> pathlib.Path:
    """Return the full path to the ``telgrid`` file

    Atmospheric Telluric Grid files are not part of the PypeIt package itself
    due to their large (~4-8GB) size.  These files are hosted remotely (see
    the PyepIt documentation), and only the ``telgrid`` files needed by a
    particular user are downloaded to the local machine.

    This function checks for the local existance of the ``telgrid`` file, and
    downloads it from the remote server using AstroPy's ``download_file()``
    function.  The file downloaded in this fashion is kept in the PypeIt
    cache (nominally ``~/.pypeit/cache``) and is not placed into the package
    directory itself.

    As most users will need only a small number of ``telgrid`` files for thier
    particular reductions, the remote fetch will only occur once per file.

    Args:
        telgrid_file (str):
          The base filename of the ``telgrid`` file to be located

    Returns:
        :obj:`pathlib.Path`: The full path to the ``telgrid`` file
    """
    # Full path within the package data structure:
    telgrid_path = Paths.telgrid / telgrid_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for most installations
    if not telgrid_path.is_file():

        # Output a brief warning for now -- makes it easier to find in the output
        msgs.info(f"telgrid file {telgrid_file} does not exist in{msgs.newline()}"
                  "the package directory.  Checking cache or downloading the file now.")

        telgrid_path = fetch_remote_file(telgrid_file, 'telluric/atm_grids', remote_host='s3_cloud')

        # If a development version, MOVE into the package directory, point path there
        if ".dev" in __version__:
            shutil.move(telgrid_path, Paths.telgrid / telgrid_file)
            telgrid_path = Paths.telgrid / telgrid_file

    # Return the path to the `telgrid` file
    return telgrid_path


def get_linelist_filepath(linelist_file: str) -> pathlib.Path:
    """Return the full path to the ``linelist`` file

    It is desired to allow users to utilize their own arc line lists for
    wavelength calibration without modifying the distributed version of the
    package.  We can utilize the ``astropy`` download/cache system added
    previously to include this functionality.

    Using the script ``pypeit_install_linelist``, custom arc line lists can be
    installed into the PypeIt cache (nominally ``~/.pypeit/cache``), and are
    not placed into the package directory itself.

    Given the line list filename, this function checks first for the existance
    of the file in the package directory, then checks the PypeIt cache.  For
    all built-in line lists, this function returns the file location within the
    package directory.  For user-supplied lists that were installed using the
    script, this function returns the location within the cache.

    The cache keeps a hash of the file URL, which contains the PypeIt version
    number.  As users update to newer versions, the ``linelist`` files must be
    reinstalled using the included script.

    Args:
        linelist_file (str):
          The base filename of the ``linelist`` file to be located

    Returns:
        :obj:`pathlib.Path`: The full path to the ``linelist`` file
    """
    # Full path within the package data structure:
    linelist_path = Paths.linelist / linelist_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should only be the case for user-installed line lists
    if not linelist_path.is_file():

        linelist_path = fetch_remote_file(linelist_file, "arc_lines/lists")

        # Output an informational message
        msgs.info(f"Using line list file {linelist_file}{msgs.newline()}"
                  "that was found in the cache.")

    # Return the path to the `linelist` file
    return linelist_path


def get_extinctfile_filepath(extinction_file: str) -> pathlib.Path:
    """Return the full path to the ``extinction`` file

    Unlike other get_*_filepath() functions, the extinction files are included
    with the PyPI distribution since they are small text files.  The purpose
    of this function is to be able to load in user-installed extinction files
    from observatories not already included.  Users may self-install such
    files using the ``pypeit_install_extinctfile`` script.

    Args:
        extinction_file (str):
          The base filename of the ``extinction`` file to be located

    Returns:
        :obj:`pathlib.Path`: The full path to the ``extinction`` file
    """
    # Full path within the package data structure:
    extinction_path = Paths.extinction / extinction_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case only for user-installed extinction files
    if not extinction_path.is_file():

        extinction_path = fetch_remote_file(extinction_file, "extinction")

        # Output an informational message
        msgs.info(f"Using extinction file {extinction_file}{msgs.newline()}"
                  "that was found in the cache.")

    # Return the path to the `extinction` file
    return extinction_path


# AstroPy download/cache infrastructure ======================================#
def fetch_remote_file(
    filename: str,
    filetype: str,
    remote_host='github',
    install_script=False,
    force_update=False,
    full_url=None
) -> pathlib.Path:
    """Use ``astropy.utils.data`` to fetch file from remote or cache

    The function ``download_file()`` will first look in the local cache (the option
    ``cache=True`` is used with this function to retrieve downloaded files from
    the cache, as needed) before downloading the file from the remote server.

    The remote file can be forcibly downloaded through the use of ``force_update``.

    Args:
        filename (str):
          The base filename to search for
        filetype (str):
          The subdirectory of ``pypeit/data/`` in which to find the file
          (e.g., ``arc_lines/reid_arxiv`` or ``sensfuncs``)
        remote_host (:obj:`str`, optional):
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  Defaults to 'github'].
        install_script (:obj:`bool`, optional):
          This function is being called from an install script (i.e.,
          ``pypeit_install_telluric``) -- relates to warnings displayed.
          Defaults to False.
        force_update (:obj:`bool`, optional):
          Force ``astropy.utils.data.download_file()`` to update the cache by
          downloading the latest version.  Defaults to False.
        full_url (:obj:`str`, optional):
          The full url (i.e., skip _build_remote_url())  Defaults to None.

    Returns:
        :obj:`pathlib.Path`: The local path to the desired file in the cache
    """
    # In some cases, we have the full URL already, but most of the time not
    if full_url:
        remote_url, sources = full_url, None
    else:
        remote_url, sources = _build_remote_url(filename, filetype, remote_host=remote_host)

    if remote_host == "s3_cloud" and not install_script:
        # Display a warning that this may take a while, and the user
        #   may wish to download using the `pypeit_install_telluric` script
        msgs.warn(f"You may wish to download {filename}{msgs.newline()}"
                    f"independently from your reduction by using the{msgs.newline()}"
                    "`pypeit_install_telluric` script.")

    # Get the file from cache, if available, or download from the remote server
    try:
        cache_fn = astropy.utils.data.download_file(
            remote_url,
            sources=sources,
            timeout=10,
            cache="update" if force_update else True,
            pkgname="pypeit"
        )

    except urllib.error.URLError as error:
        if remote_host == "s3_cloud" and (requests.head(sources[0]).status_code in
                                         [requests.codes.forbidden, requests.codes.not_found]):

            err_msg = (
                f"The file {filename}{msgs.newline()}"
                f"is not hosted in the cloud.  Please download this file from{msgs.newline()}"
                f"the PypeIt Google Drive and install it using the script{msgs.newline()}"
                f"pypeit_install_telluric --local.  See instructions at{msgs.newline()}"
                "https://pypeit.readthedocs.io/en/latest/installing.html#additional-data"
            )

        elif filetype == "arc_lines/lists":
            err_msg = (
                f"Cannot find local arc line list {filename}{msgs.newline()}"
                f"Use the script `pypeit_install_linelist` to install{msgs.newline()}"
                f"your custom line list into the cache.  See instructions at{msgs.newline()}"
                "https://pypeit.readthedocs.io/en/latest/wave_calib.html#line-lists"
            )

        elif filetype == "extinction":
            err_msg = (
                f"Cannot find local extinction file {filename}{msgs.newline()}"
                f"Use the script `pypeit_install_extinctfile` to install{msgs.newline()}"
                f"your custom extinction file into the cache.  See instructions at{msgs.newline()}"
                "https://pypeit.readthedocs.io/en/latest/fluxing.html#extinction-correction"
            )

        else:
            err_msg = (
                f"Error downloading {filename}: {error}{msgs.newline()}"
                f"URL attempted: {remote_url}{msgs.newline()}"
                f"If the error relates to the server not being found,{msgs.newline()}"
                f"check your internet connection.  If the remote server{msgs.newline()}"
                f"name has changed, please contact the PypeIt development{msgs.newline()}"
                "team."
            )

        # Raise the appropriate error message
        msgs.error(err_msg)

    # If no error, return the pathlib object
    return pathlib.Path(cache_fn)


def search_cache(pattern_str: str) -> list[pathlib.Path]:
    """Search the cache for items matching a pattern string

    This function searches the PypeIt cache for files whose URL keys contain
    the input ``pattern_str``, and returns the local filesystem path to those
    files.

    Args:
        pattern_str (str):
          The filename pattern to match

    Returns:
        list: The list of local paths for the objects whose normal filenames
        match the ``pattern_str``.
    """
    # Retreive a dictionary of the cache contents
    cache_dict = astropy.utils.data.cache_contents(pkgname="pypeit")

    # Return just the local filenames' Paths for items matching the `pattern_str`
    return [pathlib.Path(cache_dict[url]) for url in cache_dict if pattern_str in url]


def write_file_to_cache(filename, cachename, filetype, remote_host="github"):
    """Use ``astropy.utils.data`` to save local file to cache

    This function writes a local file to the PypeIt cache as if it came from a
    remote server.  This is useful for being able to use locally created or
    separately downloaded files in place of PypeIt-distributed versions.

    Args:
        filename (str):
          The filename of the local file to save
        cachename (str):
          The name of the cached version of the file
        filetype (str):
          The subdirectory of ``pypeit/data/`` in which to find the file
          (e.g., ``arc_lines/reid_arxiv`` or ``sensfuncs``)
        remote_host (:obj:`str`, optional):
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  Defaults to 'github'.
    """
   # Build the `url_key` as if this file were in the remote location
    url_key, _ = _build_remote_url(cachename, filetype, remote_host=remote_host)

    # Use `import file_to_cache()` to place the `filename` into the cache
    astropy.utils.data.import_file_to_cache(url_key, filename, pkgname="pypeit")


def _build_remote_url(f_name, f_type, remote_host=""):
    """Build the remote URL for the ``astropy.utils.data`` functions

    This function keeps the URL-creation in one place.  In the event that files
    are moved from GitHub or S3_Cloud, this is the only place that would need
    to be changed.

    Args:
        f_name (str):
          The base filename to search for
        f_type (str):
          The subdirectory of ``pypeit/data/`` in which to find the file
          (e.g., ``arc_lines/reid_arxiv`` or ``sensfuncs``)
        remote_host (:obj:`str`, optional):
          The remote host scheme.  Currently only 'github' and 's3_cloud' are
          supported.  Defaults to ''.

    Returns:
        tuple: URLs for both cache storage name and current download link:

           * url (str): The URL of the ``f_name`` of ``f_type`` on server ``remote_host``
           * sources (:obj:`list` or :obj:`None`): For 's3_cloud', the list of URLs to
                actually try, passed to ``download_file()``, used in the event that the
                S3 location changes.  We maintain the static URL for the name to prevent
                re-downloading of large data files in the event the S3 location changes
                (but the file itself is unchanged).  If None (e.g. for 'github'), then
                ``download_file()`` is unaffected, and the ``url`` (above) is what
                controls the download.
    """
    if remote_host == "github":
        # Hard-wire the URL based on PypeIt Version
        data_url = (
            "https://raw.githubusercontent.com/pypeit/PypeIt/"
            f"{'develop' if '.dev' in __version__ else __version__}/pypeit/data"
        )
        return f"{data_url}/{f_type}/{f_name}", None

    if remote_host == "s3_cloud":
        # Build up the (permanent, fake) `remote_url` and (fluid, real) `sources` for S3 Cloud
        return (f"https://s3.cloud.com/pypeit/{f_type}/{f_name}",
               [f"https://{_get_s3_hostname()}/pypeit/{f_type}/{f_name}"])

    msgs.error(f"Remote host type {remote_host} is not supported for package data caching.")


def _get_s3_hostname():
    """Get the current S3 hostname from the package file

    Since the S3 server hostname used to hold package data such as telluric
    atmospheric grids may change periodically, we keep the current hostname
    in a separate file (s3_url.txt), and pull the current version from the
    PypeIt ``release`` branch whenever needed.

    NOTE: When/if the S3 URL changes, the ``release`` branch version of
          ``s3_url.txt`` can be updated easily with a hotfix PR, and this
          routine will pull it.

    If GitHub cannot be reached, the routine uses the version of ``s3_url.txt``
    included with the package distribution.

    Returns:
        str: The current hostname URL of the S3 server holding package data
    """
    # Try getting the latest version from the server, else use what's included
    try:
        remote_url = (
            github.Github()
            .get_repo("pypeit/PypeIt")
            .get_contents("pypeit/data/s3_url.txt", "release")
            .download_url
        )
        filepath = astropy.utils.data.download_file(
            remote_url, cache="update", timeout=10, pkgname="pypeit"
        )
    except (
        requests.exceptions.ConnectionError,
        requests.exceptions.RequestException,
        urllib.error.URLError,
        github.GithubException
    ):
        filepath = Paths.data / "s3_url.txt"

    # Open the file and return the URL
    with open(filepath, "r", encoding="utf-8") as fileobj:
        return fileobj.read().strip()


# Loading Functions for Particular File Types ================================#
def load_telluric_grid(filename: str):
    """Load a telluric atmospheric grid

    NOTE: This is where the path to the data directory is added!

    Args:
        filename (str):
          The filename (NO PATH) of the telluric atmospheric grid to use.

    Returns:
        (:obj:`astropy.io.fits.HDUList`): Telluric Grid FITS HDU list
    """
    # Check for existance of file parameter
    if not filename:
        msgs.error("No file specified for telluric correction.  "
                   "See https://pypeit.readthedocs.io/en/latest/telluric.html")

    # Get the data path for the filename, whether in the package directory or cache
    file_with_path = get_telgrid_filepath(filename)

    # Check for existance of file
    # NOTE: With the use of `get_telgrid_filepath()`, this should never run
    if not file_with_path.is_file():
        msgs.error(f"File {file_with_path} is not on your disk.  "
                   "You likely need to download the Telluric files.  "
                   "See https://pypeit.readthedocs.io/en/release/installing.html"
                   "#atmospheric-model-grids")

    return io.fits_open(file_with_path)


def load_thar_spec():
    """Load the archived ThAr spectrum

    NOTE: This is where the path to the data directory is added!

    Args:
        filename (str):
          The filename (NO PATH) of the telluric atmospheric grid to use.

    Returns:
        (:obj:`astropy.io.fits.HDUList`): ThAr Spectrum FITS HDU list
    """
    return io.fits_open(Paths.arclines / 'thar_spec_MM201006.fits')


def load_sky_spectrum(sky_file):
    """Load a sky spectrum into an XSpectrum1D object

    NOTE: This is where the path to the data directory is added!

    .. todo::

        Try to eliminate the XSpectrum1D dependancy

    Args:
        sky_file (str):
          The filename (NO PATH) of the sky file to use.

    Returns:
        (:obj:`XSpectrum1D`): Sky spectrum
    """
    return xspectrum1d.XSpectrum1D.from_file(str(Paths.sky_spec / sky_file))
