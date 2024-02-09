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
module and not call, e.g. `importlib.resources.files
<https://docs.python.org/3/library/importlib.resources.html#importlib.resources.files>`__
or attempt to otherwise directly access the package directory structure.  In
this way, if structural changes to this directory are needed, only this module
need be modified and the remainder of the package can remain ignorant of those
changes and continue to call the paths supplied by this module.

Furthermore, all paths returned by this module are :obj:`pathlib.Path` objects
rather than pure strings, with all of the functionality therein contained.

Most (by number) of the package data files here are distributed with the
``PypeIt`` package and are accessed via the :class:`~pypeit.data.utils.Paths`
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
:class:`~pypeit.data.utils.Paths` class and access the enclosed files similarly
to the Vega example above.

.. include:: ../include/links.rst
"""
from importlib import resources
import pathlib
import shutil
import urllib.error

import astropy.utils.data
import github
from linetools.spectra import xspectrum1d
import requests

# NOTE: To avoid circular imports, avoid (if possible) importing anything from
# pypeit!  msgs is the possible exception for now.
from pypeit import msgs

# AstroPy download/cache infrastructure ======================================#
def fetch_remote_file(
    filename: str,
    filetype: str,
    remote_host: str='github',
    install_script: bool=False,
    force_update: bool=False,
    full_url: str=None,
) -> pathlib.Path:
    """
    Use `astropy.utils.data`_ to fetch file from remote or cache

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
            supported.  Defaults to 'github'.
        install_script (:obj:`bool`, optional):
            This function is being called from an install script (i.e.,
            ``pypeit_install_telluric``) -- relates to warnings displayed.
            Defaults to False.
        force_update (:obj:`bool`, optional):
            Force `astropy.utils.data.download_file`_ to update the cache by
            downloading the latest version.  Defaults to False.
        full_url (:obj:`str`, optional):
            The full url.  If None, use :func:`_build_remote_url`).  Defaults to
            None.

    Returns:
        `pathlib.Path`_: The local path to the desired file in the cache
    """
    # In some cases, we have the full URL already, but most of the time not
    if full_url:
        remote_url, sources = full_url, None
    else:
        remote_url, sources = _build_remote_url(filename, filetype, remote_host=remote_host)

    if remote_host == "s3_cloud" and not install_script:
        # Display a warning that this may take a while, and the user
        #   may wish to download using the `pypeit_install_telluric` script
        # TODO: Is this only true for telluric files?
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
            pkgname="pypeit",
        )

    except urllib.error.URLError as error:
        if remote_host == "s3_cloud" and (
            requests.head(sources[0]).status_code
            in [requests.codes.forbidden, requests.codes.not_found]
        ):
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

    except TimeoutError as error:
        msgs.error(f"Timeout Error encountered: {error}")

    # If no error, return the pathlib object
    return pathlib.Path(cache_fn).resolve()


def search_cache(pattern_str: str) -> list[pathlib.Path]:
    """
    Search the cache for items matching a pattern string

    This function searches the PypeIt cache for files whose URL keys contain
    the input ``pattern_str``, and returns the local filesystem path to those
    files.

    Args:
        pattern_str (:obj:`str`):
            The filename pattern to match

    Returns:
        :obj:`list`: The list of local paths for the objects whose normal
        filenames match the ``pattern_str``.
    """
    # Retreive a dictionary of the cache contents
    cache_dict = astropy.utils.data.cache_contents(pkgname="pypeit")

    # Return just the local filenames' Paths for items matching the `pattern_str`
    return [pathlib.Path(cache_dict[url]) for url in cache_dict if pattern_str in url]


def write_file_to_cache(filename: str, cachename: str, filetype: str, remote_host: str="github"):
    """
    Use `astropy.utils.data`_ to save local file to cache

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


def _build_remote_url(f_name: str, f_type: str, remote_host: str=None):
    """
    Build the remote URL for the `astropy.utils.data`_ functions

    This function keeps the URL-creation in one place.  In the event that files
    are moved from GitHub or S3_Cloud, this is the only place that would need
    to be changed.

    Parameters
    ----------
    f_name : str
        The base filename to search for
    f_type : str
        The subdirectory of ``pypeit/data/`` in which to find the file
        (e.g., ``arc_lines/reid_arxiv`` or ``sensfuncs``)
    remote_host : :obj:`str`, optional
        The remote host scheme.  Currently only 'github' and 's3_cloud' are
        supported.  Defaults to None.

    Returns
    -------
    url : str
        The URL of the ``f_name`` of ``f_type`` on server ``remote_host``
    sources : :obj:`list` or :obj:`None`
        For 's3_cloud', the list of URLs to actually try, passed to
        `astropy.utils.data.download_file`_, used in the event that the S3
        location changes.  We maintain the static URL for the name to prevent
        re-downloading large data files in the event the S3 location changes
        (but the file itself is unchanged).  If None (e.g. for 'github'), then
        `astropy.utils.data.download_file`_ is unaffected, and the ``url``
        (above) is what controls the download.
    """
    if remote_host == "github":
        # Hard-wire the URL based on PypeIt Version
        data_url = (
            "https://raw.githubusercontent.com/pypeit/PypeIt/"
            f"{'develop' if '.dev' in __version__ else __version__}"
            f"/pypeit/data/{f_type}/{f_name}"
        )
        return data_url, None

    if remote_host == "s3_cloud":
        # Build up the (permanent, fake) `remote_url` and (fluid, real) `sources` for S3 Cloud
        return (f"https://s3.cloud.com/pypeit/{f_type}/{f_name}",
               [f"https://{_get_s3_hostname()}/pypeit/{f_type}/{f_name}"])

    msgs.error(f"Remote host type {remote_host} is not supported for package data caching.")


def _get_s3_hostname() -> str:
    """
    Get the current S3 hostname from the package file

    Since the S3 server hostname used to hold package data such as telluric
    atmospheric grids may change periodically, we keep the current hostname
    in a separate file (s3_url.txt), and pull the current version from the
    PypeIt ``release`` branch whenever needed.

    .. note::

        When/if the S3 URL changes, the ``release`` branch version of
        ``s3_url.txt`` can be updated easily with a hotfix PR, and this routine
        will pull it.

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
        github.GithubException,
        TimeoutError,
    ):
        filepath = Paths.data / "s3_url.txt"

    # Open the file and return the URL
    with open(filepath, "r", encoding="utf-8") as fileobj:
        return fileobj.read().strip()

