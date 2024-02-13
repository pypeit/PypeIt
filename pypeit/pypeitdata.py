# -*- coding: utf-8 -*-
"""
PypeIt uses the `astropy.utils.data`_ caching system to limit the size of its
package distribution in PyPI by enabling on-demand downloading of reference
files needed for specific data-reduction steps.  This module provides the class
used to access to the data files in the code base.

The :mod:`~pypeit.cache` module implements the low-level function used to
interface with the PypeIt cache.  To get the location of your pypeit cache (by
default ``~/.pypeit/cache``) you can run:

.. code-block:: python

    import astropy.config.paths
    print(astropy.config.paths.get_cache_dir('pypeit'))

Every time PypeIt is imported, a new instance of
:class:`~pypeit.pypeitdata.PypeItDataPaths` is created, and this instance is
used to set paths to PypeIt data files.  "Data files" in this context
essentially refers to anything in the ``pypeit/data`` directory tree; see the
attributes of :class:`~pypeit.pypeitdata.PypeItDataPaths` for the list of
directories that can be directly accessed.

Some of files in these directories are included in the package distribution, but
most are not.  Regardless, the :class:`~pypeit.pypeitdata.PypeItDataPaths`
object should be used to define the relevant file paths.  For example, to access
a given NIST line list, one would run:

.. code-block:: python

    from pypeit import dataPaths
    thar = dataPaths.nist.get_file_path('ThAr_vacuum.ascii')

All of the attributes of :class:`~pypeit.pypeitdata.PypeItDataPaths` are
:class:`~pypeit.pypeitdata.PypeItDataPath` objects, such that the code above is
effectively equivalent to:

.. code-block:: python

    from pypeit.pypeitdata import PypeItDataPath
    thar = PypeItDataPath('arc_lines/NIST').get_file_path('ThAr_vacuum.ascii')

Although :class:`~pypeit.pypeitdata.PypeItDataPath` objects can be treated
similarly to `Path`_ objects, you should always use the
:func:`~pypeit.pypeitdata.PypeItDataPath.get_file_path` function to access the
relevant file path.  Behind the scenes, this function looks for the requested
file in your package distribution and/or downloads the file to your cache before
returning the appropriate path.

Data directories that *MUST* exist as part of the package distribution are:

    - ``pypeit/data/spectrographs/keck_deimos/gain_ronoise``

.. include:: ../include/links.rst
"""
from importlib import resources
import pathlib
import shutil

from IPython import embed

from pypeit import msgs
from pypeit import cache

# NOTE: A better approach may be to subclass from Path.  I briefly tried that,
# but quickly realized it was going to be more complicated than I'd hoped.  This
# is a clean and sufficient solution for now.
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

    Attributes:
        data (`Path`_):
            Path to the top-level data directory on the user's system.
        path (`Path`_):
            Path to the specific data directory.
    """

    _data = resources.files('pypeit') / 'data'
    """
    Package-level data path.  This is defined without checking that it exists!
    """

    def __init__(self, subdirs):
        self.data = self.check_isdir(self._data)
        self.path = self.check_isdir(self.data / subdirs)

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
            :class:`~pypeit.pypmsgs.PypeItPathError`:
                Raised if the requested contents *do not exist*.
        """
        if (self.path / p).is_dir():
            return PypeItDataPath(str((self.path / p).relative_to(self.data)))
        if (self.path / p).is_file(): 
            return self.path / p
        msgs.error(f'{str(self.path / p)} is not a valid PypeIt data path or is a file '
                   'that does not exist.', cls='PypeItPathError')

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
            :class:`~pypeit.pypmsgs.PypeItPathError`:
                Raised if the path does not exist or is not a directory.
        """
        if not path.is_dir():
            msgs.error(f"Unable to find {path}.  Check your installation.", cls='PypeItPathError')
        return path
    
    @staticmethod
    def _get_file_path_return(f, return_format):
        """
        Convience function that formats the return of :func:`get_file_path`,
        depending on whether or not the "format" of the file is requested.

        Args:
            f (`Path`_):
                The file path to return.
            return_format (:obj:`bool`):
                If True, parse and return the file suffix (e.g., 'fits').  If
                False, only the file path is returned.

        Returns:
            `Path`_, tuple: The file path and, if requested, the file format.
            See the arguments list.
        """
        if return_format:
            return f, PypeItDataPath._parse_format(f)
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

    def get_file_path(self, data_file, to_pkg=None, return_format=False):
        """
        Return the path to a file.

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

        Returns:
            `Path`_, tuple: The file path and, if requested, the file format;
            see ``return_format``.
        """
        # Make sure the file is a Path object
        _data_file = pathlib.Path(data_file).resolve()

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

        # If it does not, inform the user and download it into the cache.
        # NOTE: This should not be required for from-source (dev) installations.
        msgs.info(f'{data_file} does not exist in the expected package directory ({self.path}).  '
                    'Checking cache or downloading the file now.')

        # TODO: Fix self.path
        _cached_file = cache.fetch_remote_file(data_file, self.path)

        # If requested, either create a symlink to the cached file or move it
        # into the package data directory
        if to_pkg is not None:
            # Create and return values for the symlink
            if to_pkg == 'symlink':
                _data_file.symlink_to(_cached_file)
            elif to_pkg == 'move':
                shutil.move(_cached_file, _data_file)
            return self._get_file_path_return(_data_file, return_format)
        # Return the cached file path 
        return self._get_file_path_return(_cached_file, return_format)
    

class PypeItDataPaths:
    """
    List of hardwired data path objects, primarily for developers.

    The top-level directory for all attributes is ``pypeit/data``.  All of these
    directories should, at minimum, include a README file that is
    version-controlled and hosted by GitHub.

    Attributes:
        tests (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``tests/``
        telgrid (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``telluric/atm_grids``
        tel_model (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``telluric/models``
        arclines (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``arc_lines``
        reid_arxiv (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``arc_lines/reid_arxiv``
        linelist (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``arc_lines/lists``
        nist (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``arc_lines/NIST``
        arc_plot (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``arc_lines/plots``
        standards (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``standards``
        extinction (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``extinction``
        skisim (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``skisim``
        filters (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``filters``
        sensfunc (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``sensfuncs``
        sky_spec (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``sky_spec``
        static_calibs (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``static_calibs``
        spectrographs (:class:`~pypeit.pypeitdata.PypeItDataPath`):
            Path to ``spectrographs``
    """
    def __init__(self):
        # Tests
        self.tests = PypeItDataPath('tests')

        # Telluric
        self.telgrid  = PypeItDataPath('telluric/atm_grids')
        self.tel_model = PypeItDataPath('telluric/models')

        # Wavelength Calibrations
        self.arclines = PypeItDataPath('arc_lines')
        self.reid_arxiv = self.arclines / 'reid_arxiv'
        self.linelist = self.arclines / 'lists'
        self.nist = self.arclines / 'NIST'
        self.arc_plot = self.arclines / 'plots'

        # Flux Calibrations
        self.standards = PypeItDataPath('standards')
        self.extinction = PypeItDataPath('extinction')
        self.skisim = PypeItDataPath('skisim')
        self.filters = PypeItDataPath('filters')
        self.sensfunc = PypeItDataPath('sensfuncs')

        # Other
        self.sky_spec = PypeItDataPath('sky_spec')
        self.static_calibs = PypeItDataPath('static_calibs')
        self.spectrographs = PypeItDataPath('spectrographs')


