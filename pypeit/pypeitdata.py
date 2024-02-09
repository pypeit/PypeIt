# -*- coding: utf-8 -*-
"""
Basic class used to facilitate the cache system used for data files.

.. include:: ../include/links.rst
"""
from importlib import resources
import pathlib

from pypeit import msgs

class PypeItDataPaths:
    """
    List of hardwired paths within the pypeit.data module

    Each `@property` method returns a :obj:`pathlib.Path` object
    """
    def __init__(self):
        # Determine the location of the data directory on this system
        self.pkg = resources.files('pypeit')
        self.data = self.check_isdir(self.pkg / 'data')

        # Telluric
        self.telluric = self.check_isdir(self.data / 'telluric')
        self.telgrid  = self.check_isdir(self.telluric / 'atm_grids')
        self.tel_model = self.check_isdir(self.telluric / 'models')

        # Wavelength Calibrations
        self.arclines = self.check_isdir(self.data / 'arc_lines')
        self.reid_arxiv = self.check_isdir(self.arc_lines / 'reid_arxiv')
        self.linelist = self.check_isdir(self.arc_lines / 'lists')
        self.nist = self.check_isdir(self.arc_lines / 'NIST')
        self.arc_plot = self.check_isdir(self.arc_lines / 'plots')

        # Flux Calibrations
        self.standards = self.check_isdir(self.data / 'standards')
        self.extinction = self.check_isdir(self.data / 'extinction')
        self.skisim = self.check_isdir(self.data / 'skisim')
        self.filters = self.check_isdir(self.data / 'filters')
        self.sensfuncs = self.check_isdir(self.data / 'sensfuncs')

        # Other
        self.sky_spec = self.check_isdir(self.data / 'sky_spec')
        self.static_calibs = self.check_isdir(self.data / 'static_calibs')
        self.spectrographs = self.check_isdir(self.data / 'spectrographs')

    @staticmethod
    def check_isdir(path:pathlib.Path) -> pathlib.Path:
        """Check that the hardwired directory exists

        If yes, return the directory path, else raise an error message
        """
        if not path.is_dir():
            msgs.error(f"Unable to find {path}.  Check your installation.")
        return path
    
    @staticmethod
    def _get_file_path_return(f, return_format):
        if return_format:
            return f, f.suffix.replace('.','').lower()
        return f

    def get_file_path(self, data_root, data_file, symlink_in_pkdir=False, return_format=False):
        """
        Test
        """
        _package_root = resources.files('pypeit')

        # Make sure the file is a Path object
        _data_file = pathlib.Path(data_file).resolve()

        # Check if the file exists on disk, as provided 
        if _data_file.is_file():
            # If so, assume this points directly to the file to be read
            return _get_file_path_return(_data_file, return_format)

        # Otherwise, construct the file name given the root path:
        _data_file = _package_root / data_root / data_file

        # If the file exists, return it 
        if _data_file.is_file():
            # Return the full path and the file format
            return _get_file_path_return(_data_file, return_format)

        # If it does not, inform the user and download it into the cache.
        # NOTE: This should not be required for from-source (dev) installations.
        msgs.info(f'{data_file} does not exist in the expected package directory ({data_root}).  '
                    'Checking cache or downloading the file now.')

        _cached_file = fetch_remote_file(data_file, data_root)

        # If requested, create a symlink to the cached file in the package data
        # directory
        if symlink_in_pkgdir:
            # Create and return values for the symlink
            _data_file.symlink_to(_cached_file)
            return _get_file_path_return(_data_file, return_format)
    
        return _get_file_path_return(_cached_file, return_format)

def _get_reid_arxiv_filepath(arxiv_file):
    return get_data_file(Paths.reid_arxiv, arxiv_file)

def _get_skisim_filepath(skisim_file):
    return get_data_file(Paths.skisim, skisim_file)

def _get_sensfunc_filepath(sensfunc_file, symlink_in_pkgdir=False):
    return get_data_file(Paths.sensfuncs, sensfunc_file, symlink_in_pkgdir=symlink_in_pkgdir)

def _get_telgrid_filepath(telgrid_file):
    return get_data_file(Paths.telgrid, telgrid_file)

def _get_linelist_filepath(linelist_file):
    return get_data_file(Paths.linelist, linelist_file)

def _get_extinctfile_filepath(extinction_file):
    return get_data_file(Paths.extinction, extinction_file)

def get_test_filepath(test_file):
    return get_data_file(Paths.tests, test_file)

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
        arxiv_file (:obj:`str`):
          The base filename of the ``reid_arxiv`` file to be located

    Returns:
        tuple: The full path and whether the path is in the cache:

           * reid_path (:obj:`~pathlib.Path`): The full path to the ``reid_arxiv`` file
           * arxiv_fmt (:obj:`str`): The extension of the ``reid_arxiv`` file (format)
    """
    # First check that what is passed in works
    if not isinstance(arxiv_file, (str, pathlib.Path)):
        msgs.error(f"Incorrect or non-existent arxiv file specified: {arxiv_file}")
    if isinstance(arxiv_file, pathlib.Path):
        arxiv_file = str(arxiv_file)

    # Check if the `arxiv_file` already comes with a full path
    if (reid_path := pathlib.Path(arxiv_file)).name != arxiv_file:
        # Check existence, return it with format
        if not reid_path.is_file():
            msgs.error(f"Incorrect or non-existent arxiv file specified: {reid_path}")
        return reid_path, reid_path.suffix.replace('.','').lower()

    # Else, full path within the package data structure:
    reid_path = Paths.reid_arxiv / arxiv_file

    # Check if the file does NOT exist in the package directory
    # NOTE: This should be the case for all but from-source installations
    if not reid_path.is_file():

        # Output an informational message
        msgs.info(f"reid_arxiv file {arxiv_file} does not exist in{msgs.newline()}"
                  "the package directory.  Checking cache or downloading the file now.")

        reid_path = fetch_remote_file(arxiv_file, "arc_lines/reid_arxiv")

    # Return the path to the `reid_arxiv` file, and the file format
    return reid_path, arxiv_file.split('.')[-1].lower()


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


def get_sensfunc_filepath(sensfunc_file: str, symlink_in_pkgdir: bool=False) -> pathlib.Path:
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
        :obj:`pathlib.Path`: The full path to the ``sensfunc`` file
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
    remote_host: str='github',
    install_script: bool=False,
    force_update: bool=False,
    full_url: str=None,
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
        msgs.error(f"Timeout Error enountered: {error}")

    # If no error, return the pathlib object
    return pathlib.Path(cache_fn).resolve()


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


def write_file_to_cache(filename: str, cachename: str, filetype: str, remote_host: str="github"):
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


def _build_remote_url(f_name: str, f_type: str, remote_host: str=""):
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
        github.GithubException,
        TimeoutError,
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
        (`astropy.io.fits.HDUList`_): Telluric Grid FITS HDU list
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
        (`astropy.io.fits.HDUList`_): ThAr Spectrum FITS HDU list
    """
    return io.fits_open(Paths.arclines / 'thar_spec_MM201006.fits')


def load_sky_spectrum(sky_file: str) -> xspectrum1d.XSpectrum1D:
    """Load a sky spectrum into an XSpectrum1D object

    NOTE: This is where the path to the data directory is added!

    .. todo::

        Try to eliminate the XSpectrum1D dependancy

    Args:
        sky_file (str):
          The filename (NO PATH) of the sky file to use.

    Returns:
        (`linetools.spectra.xspectrum1d.XSpectrum1D`_): Sky spectrum
    """
    return xspectrum1d.XSpectrum1D.from_file(str(Paths.sky_spec / sky_file))

