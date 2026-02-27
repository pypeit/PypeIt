# -*- coding: utf-8 -*-
"""
PypeIt uses the `astropy.utils.data`_ caching system to limit the size of its
package distribution in PyPI by enabling on-demand downloading of reference
files needed for specific data-reduction steps.  This module provides the
low-level utility functions that interface with the cache.

Access to the data files are handled in the code base using the
:class:`~pypeit.pypeitdata.PypeItDataPaths` object instantiated every time
PypeIt is imported.

To get the location of your pypeit cache (by default ``~/.pypeit/cache``) you
can run:

.. code-block:: python

    import astropy.config.paths
    print(astropy.config.paths.get_cache_dir('pypeit'))

.. note::

    If the hostname URL for the telluric atmospheric grids on S3 changes, the
    only place that needs to change is the file ``pypeit/data/s3_url.txt``.

.. include:: ../include/links.rst
"""
import shutil

from functools import reduce
from importlib import resources
import pathlib
import urllib.error
from urllib.parse import urljoin, urlparse
from datetime import datetime
import warnings

import packaging

from IPython import embed

import numpy as np

import astropy.utils.data
import github
import requests


# NOTE: pygit2 is only used for testing purposes.  It is not a requirement for a
# general user.  Hence the try block below.
try:
    from pygit2 import Repository
except ImportError:
    Repository = None
    GitError = None
else:
    from pygit2 import GitError


# NOTE: To avoid circular imports, avoid (if possible) importing anything from
# pypeit into this module!  Objects created or available in pypeit/__init__.py
# are the exceptions, for now.
#from pypeit import log
from .exceptions import PypeItError, PypeItPathError
from .version import version as __version__
#from pypeit import PypeItError, PypeItPathError
#from pypeit import __version__


__PYPEIT_DATA__ = resources.files('pypeit') / 'data'
__PYPEIT_REPO_PATH__ = 'pypeit/PypeIt'


def git_repo():
    """
    Get a reference to the local repository, if possible.
    """
    if Repository is None:
        # pygit2 not available
        return None
    try:
        return Repository(resources.files('pypeit'))
    except GitError:
        # PypeIt not in a git repo
        return None


# For development versions, try to get the branch name
def git_branch():
    """
    Return the name/hash of the currently checked out branch
    
    Returns:
        :obj:`str`: Branch name or hash. Defaults to "develop" if PypeIt is not
        currently in a repository or pygit2 is not installed.
    """
    repo = git_repo()
    if repo is None:
        return 'develop' if '.dev' in __version__ else __version__
    return str(repo.head.target) if repo.head_is_detached else str(repo.head.shorthand)


def git_remote_path():
    """
    The main path to the GitHub repository.

    This defaults to the main repository if the repository cannot be defined
    (see :func:`git_repo`) or if the "origin" remote URL cannot be determined.

    Returns:
        :obj:`str`: Remote path
    """
    repo = git_repo()
    if repo is None:
        return __PYPEIT_REPO_PATH__
    try:
        url = repo.remotes['origin'].url
    except KeyError:
        return __PYPEIT_REPO_PATH__
    return urlparse(url).path.replace('.git','').removeprefix('/')


def github_contents(repo, branch, path, recursive=True):
    """
    (Recursively) Acquire a listing of the contents of a repository directory.

    Args:
        repo (`github.Repository`_):
            Repository to search
        branch (:obj:`str`):
            Name of the branch or commit hash
        path (:obj:`str`):
            Path relative to the top-level directory of the repository to
            search.
        recursive (:obj:`bool`, optional):
            Flag to search the directory recursively.  If False, subdirectory
            names are included in the list of returned objects.  If True,
            subdirectories are removed from the listing and replaced by their
            contents; in this case the list of all objects should only include
            repository files.

    Returns:
        :obj:`list`: A list of `github.ContentFile`_ objects with the repo
        contents.
    """
    try:
        # Collect the contents
        contents = repo.get_contents(path, branch)
    except github.GithubException as e:
        raise PypeItPathError(f'{path} not found in the {branch} of the GitHub tree.') from e

    # If not searching recursively, we're done
    if not recursive:
        return contents

    # Check if any of the contents are directories 
    is_dir = [c.type == 'dir' for c in contents]
    # If not, we're done
    if not any(is_dir):
        return contents

    # For each directory, append the directory contents recursively
    is_dir = np.where(is_dir)[0]
    for indx in is_dir:
        contents.extend(github_contents(repo, branch, contents[indx].path))

    # Remove the directories from the list
    return [c for i,c in enumerate(contents) if i not in is_dir] 


def git_most_recent_tag():
    """
    Return the version number for the most recent tag and the date of its last
    commit.

    Returns:
        :obj:`tuple`: The version number and a ISO format string with the date
        of the last commit included in the tag.  If ``pygit2`` is not installed
        or no tags are found, the returned version is the same as
        ``pypeit.__version__`` and the date is None.
    """
    if Repository is None:
        return __version__, None
    repo = Repository(resources.files('pypeit'))
    tags = [packaging.version.parse(ref.split('/')[-1]) \
                for ref in repo.references if 'refs/tags' in ref]
    if len(tags) == 0:
        warnings.warn('Unable to find any tags in pypeit repository.')
        return __version__, None
    latest_version = str(sorted(tags)[-1])
    timestamp = repo.resolve_refish(f'refs/tags/{latest_version}')[0].author.time
    return latest_version, datetime.fromtimestamp(timestamp).isoformat()


# AstroPy download/cache infrastructure ======================================#
def fetch_remote_file(
    filename: str,
    filetype: str,
    remote_host: str='github',
    install_script: bool=False,
    force_update: bool=False,
    full_url: str=None,
    return_none: bool=False,
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
        return_none (:obj:`bool`, optional):
            Return None if the file is not found.  Defaults to False.

    Returns:
        `Path`_: The local path to the desired file in the cache
    """
    # In some cases, we have the full URL already, but most of the time not
    if full_url:
        remote_url, sources = full_url, None
    else:
        remote_url, sources = _build_remote_url(filename, filetype, remote_host=remote_host)

    if remote_host == "s3_cloud" and not install_script:
        # Display a warning that this may take a while, and the user may wish to
        # download use an install script
        warnings.warn(f'Note: If this file takes a while to download, you may wish to used one of '
                  'the install scripts (e.g., pypeit_install_telluric) to install the file '
                  'independent of this processing script.')

    # Get the file from cache, if available, or download from the remote server
    # TODO: Make timeout a function argument?
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
                f"The file {filename}\n"
                f"is not hosted in the cloud.  Please download this file from"
                f"the PypeIt Google Drive and install it using the script"
                f"pypeit_install_telluric --local.  See instructions at"
                "https://pypeit.readthedocs.io/en/latest/installing.html#additional-data"
            )

        elif filetype == "arc_lines/lists":
            err_msg = (
                f"Cannot find local arc line list {filename}\n"
                f"Use the script `pypeit_install_linelist` to install"
                f"your custom line list into the cache.  See instructions at"
                "https://pypeit.readthedocs.io/en/latest/wave_calib.html#line-lists"
            )

        elif filetype == "extinction":
            err_msg = (
                f"Cannot find local extinction file {filename}\n"
                f"Use the script `pypeit_install_extinctfile` to install"
                f"your custom extinction file into the cache.  See instructions at"
                "https://pypeit.readthedocs.io/en/latest/fluxing.html#extinction-correction"
            )

        elif return_none:
            return None

        else:
            err_msg = (
                f"Error downloading {filename}: {error}\n"
                f"URL attempted: {remote_url}\n"
                f"If the error relates to the server not being found,"
                f"check your internet connection.  If the remote server"
                f"name has changed, please contact the PypeIt development"
                "team."
            )

        # Raise the appropriate error message
        raise PypeItError(err_msg)

    except TimeoutError as error:
        raise PypeItError(f"Timeout Error encountered: {error}")

    # If no error, return the pathlib object
    return pathlib.Path(cache_fn).resolve()


def search_cache(pattern: str, path_only=True):
    """
    Search the cache for items matching a pattern string.

    This function searches the PypeIt cache for files whose URL keys contain the
    input ``pattern``, and returns the local filesystem path to those files.

    Args:
        pattern (:obj:`str`):
            The pattern to match within the file name of the source url.  This
            can be None, meaning that the full contents of the cache is
            returned.  However, note that setting ``pattern`` to None and
            ``path_only=True`` may not be very useful given the abstraction of
            the file names.
        path_only (:obj:`bool`, optional):
            Only return the path(s) to the files found in the cache.  If False,
            a dictionary is returned where each key is the source url, and the
            value is the local path.

    Returns:
        :obj:`list`, :obj:`dict`: If ``path_only`` is True, this is a
        :obj:`list` of local paths for the objects whose normal filenames match
        the ``pattern``.  Otherwise, this is a dictionary with keys matching the
        original source url, and the value set to the local path.
    """
    # Retrieve a dictionary of the cache contents
    contents = astropy.utils.data.cache_contents(pkgname="pypeit")
    contents = {k:pathlib.Path(v) for k, v in contents.items() if pattern is None or pattern in k}
    return list(contents.values()) if path_only else contents
    

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

    # Use `import_file_to_cache()` to place the `filename` into the cache
    astropy.utils.data.import_file_to_cache(url_key, filename, pkgname="pypeit")


def remove_from_cache(cache_url=None, pattern=None, allow_multiple=False):
    """
    Remove a previously downloaded file from the pypeit-specific
    `astropy.utils.data`_ cache.

    To specify the file, the full URL can be provided or a name used in a cache
    search.

    Args:
        cache_url (:obj:`list`, :obj:`str`, optional):
            One or more URLs in the cache to be deleted (if they exist in the
            cache).  If ``allow_multiple`` is False, this must be a single
            string.
        pattern (:obj:`str`, optional):
            A pattern to use when searching the cache for the relevant file(s).
            If ``allow_mulitple`` is False, this must return a single file,
            otherwise the function will issue a warning and nothing will be
            deleted.
        allow_multiple (:obj:`bool`, optional):
            If the search pattern yields multiple results, remove them all.
    """
    if cache_url is None:
        _url = search_cache(pattern, path_only=False)
        if len(_url) == 0:
            warnings.warn(f'Cache does not include a file matching the pattern {pattern}.')
            return
        _url = list(_url.keys())
    elif not isinstance(cache_url, list):
        _url = [cache_url]
    else:
        _url = cache_url

    if len(_url) > 1 and not allow_multiple:
        warnings.warn('Function found or was provided with multiple entries to be removed.  Either '
                  'set allow_multiple=True, or try again with a single url or more specific '
                  'pattern.  URLs passed/found are:\n' + '\n'.join(_url))
        return

    # Use `clear_download_cache` to remove the file
    for u in _url:
        astropy.utils.data.clear_download_cache(hashorurl=u, pkgname='pypeit')


def parse_cache_url(url):
    """
    Parse a URL from the cache into its relevant components.

    Parameters
    ----------
    url : :obj:`str`
        URL of a file in the pypeit cache. A valid cache URL must include either
        ``'github'`` or ``'s3.cloud'`` in its address.

    Returns
    -------
    host : :obj:`str`
        Host name, either ``'github'`` or ``'s3_cloud'``.  None if the ``url``
        is not valid.
    fork : :obj:`str`
        Fork name.  None if the ``url`` is not valid or if the host is
        ``'s3_cloud'``.
    branch : :obj:`str`
        Branch name.  None if the ``url`` is not valid or if the host is
        ``'s3_cloud'``.
    dir : :obj:`str`
        Directory name.  None if the ``url`` is not valid.
    file : :obj:`str`
        File name.  None if the ``url`` is not valid.
    """
    url_parts = urlparse(url)

    # Get the host
    if 'github' in url_parts.netloc:
        _path = pathlib.PurePosixPath(url_parts.path)
        root_tuple = _path.parts[:3] if _path.is_absolute() else ('/', *_path.parts[:2])
        fork = pathlib.PurePosixPath(*root_tuple)
        sub_path = pathlib.PurePosixPath(url_parts.path).relative_to(fork)
        branch = sub_path.parts[0]
        f_type = str(sub_path.parent.relative_to(pathlib.PurePosixPath(f'{branch}/pypeit/data')))
        return 'github', str(fork), branch, f_type, sub_path.name

    elif 's3.cloud' in url_parts.netloc:
        # NOTE: I'm assuming "s3.cloud" will always be in the url ...
        s3_root = pathlib.PurePosixPath('/pypeit')
        sub_path = pathlib.PurePosixPath(url_parts.path).relative_to(s3_root)
        return 's3_cloud', None, None, str(sub_path.parent), sub_path.name

    # Unknown host
    warnings.warn(f'URL not recognized as a pypeit cache url:\n\t{url}')
    return None, None, None, None, None


def list_cache_contents(contents):
    """
    Print the list of cache contents

    Parameters
    ----------
    contents : :obj:`dict`
        A dictionary with key-value pairs that provide the original source url
        (key) and the path to the local file (value).  This can be generated
        using :func:`search_cache`.
    """
    print(f' {"HOST":>10} {"FORK":>20} {"BRANCH":>15} {"SUBDIR":>15} {"FILE":<20}')
    for url in contents.keys():
        head, fork, branch, subdir, f = parse_cache_url(url)
        print(f' {head:>10} {"..." if fork is None else fork:>20}'
              f'{"..." if branch is None else branch:>20}'
              f' {subdir:>20} {f:<30}')


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
        parts = ['https://raw.githubusercontent.com', f'/{git_remote_path()}/', f'{git_branch()}/',
                 'pypeit/', 'data/'] + [f'{p}/' for p in pathlib.Path(f_type).parts] + [f'{f_name}']
        return reduce(lambda a, b: urljoin(a, b), parts), None

    if remote_host == "s3_cloud":
        # Build up the (permanent, fake) `remote_url` and (fluid, real) `sources` for S3 Cloud
        parts = [f'{p}/' for p in pathlib.Path(f_type).parts] + [f'{f_name}']
        parts_perm = ['https://s3.cloud.com/pypeit/'] + parts
        parts_fake = [f'https://{_get_s3_hostname()}/pypeit/'] + parts
        return reduce(lambda a, b: urljoin(a, b), parts_perm), \
                [reduce(lambda a, b: urljoin(a, b), parts_fake)]

    raise PypeItError(f"Remote host type {remote_host} is not supported for package data caching.")


def _get_s3_hostname() -> str:
    """
    Get the current S3 hostname from the package file

    Since the S3 server hostname used to hold package data such as telluric
    atmospheric grids may change periodically, we keep the current hostname in a
    separate file (``pypeit/data/s3_url.txt``), and pull the current version
    from the PypeIt ``release`` branch whenever needed.

    .. note::

        When/if the S3 URL changes, the ``release`` branch version of
        ``pypeit/data/s3_url.txt`` can be updated easily with a hotfix PR, and
        this routine will pull it.

    If GitHub cannot be reached, the routine uses the version of
    ``pypeit/data/s3_url.txt`` included with the package distribution.

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
        filepath = __PYPEIT_DATA__ / 's3_url.txt'

    # Open the file and return the URL
    with open(filepath, "r", encoding="utf-8") as fileobj:
        return fileobj.read().strip()

# NOTE: A better approach may be to subclass from Path.  I briefly tried that,
# but quickly realized it was going to be more complicated than I'd hoped.  This
# is a clean and sufficient solution for now.  As of python 3.12, the
# pathlib.Path class supports subclassing.  Return to this as soon as we require
# python>=3.12.
class PypeItDataPath:
    """
    Convenience class that enables a general interface between a pypeit data
    directory and the rest of the code, regardless of whether or not the
    directory is expected to leverage the cache system for package-installed (as
    opposed to GitHub/source-installed) versions.

    Args:
        subdirs (:obj:`str`, `Path`_):
            The subdirectories within the main ``pypeit/data`` directory that
            contains the data.
        remote_host (:obj:`str`, optional):
            The remote host for the data.  By definition, all files in this data
            path must have the *same* host.  Currently must be None,
            ``'github'``, or ``'s3_cloud'``.  If None, all the files in this
            path are expected to be local to *any* pypeit installation.

    Attributes:
        host (:obj:`str`):
            String representing the remote host
        subdirs (:obj:`str`):
            The subdirectory path within the ``pypeit/data`` directory.
        data (`Path`_):
            Path to the top-level data directory on the user's system.
        path (`Path`_):
            Path to the specific data directory.
    """

    def __init__(self, subdirs, remote_host=None):
        if remote_host not in [None, 's3_cloud', 'github']:
            raise PypeItError(f'Remote host not recognized: {self.host}')
        self.host = remote_host
        self.subdirs = subdirs
#        self.data = self.check_isdir(cache.__PYPEIT_DATA__)
        self.data = self.check_isdir(__PYPEIT_DATA__)
        self.path = self.check_isdir(self.data / self.subdirs)

    def glob(self, pattern):
        """
        Search for all contents of :attr:`path` that match the provided search
        string; see `Path.glob`_.

        .. important::

            This method *only* works for files that are on-disk in the correct
            directory.  It does *not* return any files that are in the cache or
            that have not yet been downloaded into the cache.

        Args:
            pattern (:obj:`str`):
                Search string with wild cards.

        Returns:
            generator: Generator object that provides contents matching the
            search string.
        """
        return self.path.glob(pattern)

    def __repr__(self):
        """Provide a string representation of the path.  Mimics pathlib."""
        return f"{self.__class__.__name__}('{str(self.path)}')"

    def __truediv__(self, p):
        """
        Instantiate a new path object that points to a subdirectory using the
        (true) division operator, ``/``.

        This operation should only be used for contents that *exist* on the
        users distribution.  I.e., any file that is distributed using the cache
        should *not* use this; use :func:`get_file_path` instead.

        Args:
            p (:obj:`str`, `Path`_):
                A sub-directory or file within the :attr:`path` directory.  If
                this is a sub-directory, a new :class:`PypeItDataPath` object is
                returned.  If a file, a `Path`_ object is returned.

        Returns:
            :class:`PypeItDataPath`, `Path`_: Path to contents of :attr:`path`,
            where the output type is based on the type of ``p``.

        Raises:
            :class:`~pypeit.exceptions.PypeItPathError`:
                Raised if the requested contents *do not exist*.
        """
        if (self.path / p).is_dir():
            # Create a new PypeItData object; inherit the host from the parent
            # directory
            return PypeItDataPath(str((self.path / p).relative_to(self.data)),
                                  remote_host=self.host)
        if (self.path / p).is_file():
            return self.path / p
        raise PypeItPathError(
            f'{str(self.path / p)} is not a valid PypeIt data path or is a file that does not '
            'exist.'
        )

    @staticmethod
    def check_isdir(path:pathlib.Path) -> pathlib.Path:
        """
        Check that the hardwired directory exists.

        Args:
            path (`Path`_):
                The path to check.  This *must* be a directory (not a file).

        Returns:
            `Path`_: The input path is returned if it is valid.

        Raises:
            :class:`~pypeit.exceptions.PypeItPathError`:
                Raised if the path does not exist or is not a directory.
        """
        if not path.is_dir():
            raise PypeItPathError(f'Unable to find {path}.  Check your installation.')
        return path

    @staticmethod
    def _get_file_path_return(f, return_format, format=None):
        """
        Convience function that formats the return of :func:`get_file_path`,
        depending on whether or not the "format" of the file is requested.

        Args:
            f (`Path`_):
                The file path to return.
            return_format (:obj:`bool`):
                If True, parse and return the file suffix (e.g., 'fits').  If
                False, only the file path is returned.
            format (:obj:`str`):
                If ``return_format`` is True, override (i.e. do not parse) the
                format from the file name, and use this string instead.  Ignored
                if None or ``return_format`` is False.

        Returns:
            `Path`_, tuple: The file path and, if requested, the file format.
            See the arguments list.
        """
        if return_format:
            _format = PypeItDataPath._parse_format(f) if format is None else format
            return f, _format
        return f

    @staticmethod
    def _parse_format(f):
        """
        Parse the file format, ignoring gz extensions.

        Args:
            f (`Path`_):
                File path to parse.

        Returns:
            :obj:`str`: The extension of the file that should indicate its
            format.  Any ``'.gz'`` extensions are stripped.
        """
        _f = f.with_suffix('') if f.suffix == '.gz' else f
        return _f.suffix.replace('.','').lower()

    def get_file_path(self, data_file, force_update=False, to_pkg=None, return_format=False,
                      return_none=False):
        """
        Return the path to a file.

        The file must either exist locally or be downloadable from the
        :attr:`host`.  To *define* a path to a file that does not meet these
        criteria, use ``self.path / data_file``.

        If ``data_file`` is a valid path to a file or is a file within
        :attr:`path`, the full path is returned.  Otherwise, it is assumed that
        the file is accessible remotely in the GitHub repository and can be
        downloaded using :func:`~pypeit.cache.fetch_remote_file`.  Note,
        ``data_file`` *must* be a file, not a subdirectory within :attr:`path`.

        Throughout the code base, this is the main function that should be used
        to obtain paths to files within :attr:`path`.  I.e., for *any* data in
        the ``pypeit/data`` directory, developers should be using ``p =
        path.get_file_path(file)`` instead of ``p = path / file`` to define the
        path to a data file.

        Args:
            data_file (:obj:`str`, `Path`_):
                File name or path.  See above.
            force_update (:obj:`bool`, optional):
                If the file is in the cache, force
                `astropy.utils.data.download_file`_ to update the cache by
                downloading the latest version.
            to_pkg (:obj:`str`, optional):
                If the file is in the cache, this argument affects how the
                cached file is connected to the package installation.  If
                ``'symlink'``, a symbolic link is created in the package
                directory tree that points to the cached file.  If ``'move'``,
                the cached file is *moved* (not copied) from the cache into the
                package directory tree.  If anything else (including None), no
                operation is performed; no warning is issued if the value of
                ``to_pkg`` is not one of these three options (None,
                ``'symlink'``, or ``'move'``).  This argument is ignored if
                ``data_file`` is a value path or a file within :attr:`path`.
            return_format (:obj:`bool`, optional):
                If True, the returned object is a :obj:`tuple` that includes the
                file path and its format (e.g., ``'fits'``).  If False, only the
                file path is returned.
            return_none (:obj:`bool`, optional):
                If True, return None if the file does not exist.  If False, an
                error is raised if the file does not exist.

        Returns:
            `Path`_, tuple: The file path and, if requested, the file format;
            see ``return_format``.
        """
        # Make sure the file is a Path object
        _data_file = pathlib.Path(data_file).absolute()

        # Check if the file exists on disk, as provided
        if _data_file.is_file():
            # If so, assume this points directly to the file to be read
            return self._get_file_path_return(_data_file, return_format)

        # Otherwise, construct the file name given the root path:
        _data_file = self.path / data_file

        # If the file exists, return it
        if _data_file.is_file():
            # Return the full path and the file format
            return self._get_file_path_return(_data_file, return_format)

        # Get the path to the cached file
        # NOTE: fetch_remote_file will only return the name of the cached file
        # if the file exists in the cache and force_update is False.
        subdir = str(self.path.relative_to(self.data))
#        _cached_file = cache.fetch_remote_file(data_file, subdir, remote_host=self.host,
        _cached_file = fetch_remote_file(data_file, subdir, remote_host=self.host,
                                               force_update=force_update, return_none=return_none)
        if _cached_file is None:
            warnings.warn(f'File {data_file} not found in the cache.')
            return None

        # If we've made it this far, the file is being pulled from the cache.
        if to_pkg is None:
            # The cache file is not symlinked or moved, meaning the filename is
            # always ``contents``.  That means, if the format is returned, we
            # should use the *expected* file name.
            format = PypeItDataPath._parse_format(_data_file) if return_format else None
            # Return the cached file path
            return self._get_file_path_return(_cached_file, return_format, format=format)

        # Create a symlink to the cached file or move it into the package
        # data directory
        if to_pkg == 'symlink':
            _data_file.symlink_to(_cached_file)
        elif to_pkg == 'move':
            # Move the file
            shutil.move(_cached_file, _data_file)
            # ... and delete it from the cache
#            cache.remove_from_cache(pattern=data_file)
            remove_from_cache(pattern=data_file)
        # Return the symlinked or moved file
        return self._get_file_path_return(_data_file, return_format)


class PypeItDataPaths:
    """
    List of hardwired data path objects, primarily for developers.

    The top-level directory for all attributes is ``pypeit/data``.  All of these
    directories should, at minimum, include a README file that is
    version-controlled and hosted by GitHub.  I.e., the code assumes these paths
    exist, and maintaining a version-controlled README ensures that is true,
    even if the directory is empty otherwise.
    """

    defined_paths = {
                     # Class attribute name
                     'tests':
                        # Subdirectory in pypeit/data
                        {'path': 'tests',
                        # String name for the remote host; None means the data
                        # should be in *all* installations.
                         'host': 'github'
                        },
                     # Telluric
                     'telgrid': {'path': 'telluric/atm_grids', 'host': 's3_cloud'},
                     'tel_model': {'path': 'telluric/models', 'host': None},
                     # Wavelength Calibrations
                     'arclines': {'path': 'arc_lines', 'host': None},
                     'reid_arxiv': {'path': 'arc_lines/reid_arxiv', 'host': 'github'},
                     'linelist': {'path': 'arc_lines/lists', 'host': None},
                     'nist': {'path': 'arc_lines/NIST', 'host': 'github'},
                     'arc_plot': {'path': 'arc_lines/plots', 'host': None},
                     # Flux Calibrations
                     'standards': {'path': 'standards', 'host': 'github'},
                     'extinction': {'path': 'extinction', 'host': None},
                     'skisim': {'path': 'skisim', 'host': 'github'},
                     'filters': {'path': 'filters', 'host': None},
                     'sensfunc': {'path': 'sensfuncs', 'host': 'github'},
                     # Pixel Flats
                     'pixelflat': {'path': 'pixelflats', 'host': 'github'},
                     # Other
                     'sky_spec': {'path': 'sky_spec', 'host': None},
                     'static_calibs': {'path': 'static_calibs', 'host': None},
                     'spectrographs': {'path': 'spectrographs', 'host': None},
                     # Known lists of interesting lines for data analysis
                     'line_lists': {'path': 'line_lists', 'host': None},
                    }
    """
    Dictionary providing the metadata for all the paths defined by the class.
    """

    def __init__(self):
        for key, a in PypeItDataPaths.defined_paths.items():
            setattr(self, key, PypeItDataPath(a['path'], remote_host=a['host']))


    @classmethod
    def github_paths(cls):
        """
        Return the subset paths hosted on GitHub.

        Returns:
            :obj:`dict`: A dictionary with the same format as
            :attr:`defined_paths`, but only includes those paths hosted on
            GitHub.
        """
        paths = {}
        for name, meta in cls.defined_paths.items():
            if meta['host'] != 'github':
                continue
            paths[name] = meta
        return paths

    @classmethod
    def s3_paths(cls):
        """
        Return the subset paths hosted on aws s3.

        Returns:
            :obj:`dict`: A dictionary with the same format as
            :attr:`defined_paths`, but only includes those paths hosted on
            aws s3.
        """
        paths = {}
        for name, meta in cls.defined_paths.items():
            if meta['host'] != 's3_cloud':
                continue
            paths[name] = meta
        return paths

    @classmethod
    def remote_paths(cls):
        """
        Return the subset paths with data hosted remotely.

        Returns:
            :obj:`dict`: A dictionary with the same format as
            :attr:`defined_paths`, but only includes those paths with data
            hosted remotely.
        """
        paths = cls.github_paths()
        paths.update(cls.s3_paths())
        return paths
