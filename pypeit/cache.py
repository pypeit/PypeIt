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
from functools import reduce
from importlib import resources
import pathlib
import urllib.error
from urllib.parse import urljoin, urlparse
from datetime import datetime

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
from pypeit.pypmsgs import PypeItPathError
from pypeit import msgs
from pypeit import __version__


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
        msgs.warn('Unable to find any tags in pypeit repository.')
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
        msgs.warn(f'Note: If this file takes a while to download, you may wish to used one of '
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

        elif return_none:
            return None

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
            msgs.warn(f'Cache does not include a file matching the pattern {pattern}.')
            return
        _url = list(_url.keys())
    elif not isinstance(cache_url, list):
        _url = [cache_url]
    else:
        _url = cache_url

    if len(_url) > 1 and not allow_multiple:
        msgs.warn('Function found or was provided with multiple entries to be removed.  Either '
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
    msgs.warn(f'URL not recognized as a pypeit cache url:\n\t{url}')
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

    msgs.error(f"Remote host type {remote_host} is not supported for package data caching.")


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

