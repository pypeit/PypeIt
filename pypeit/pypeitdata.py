# -*- coding: utf-8 -*-
"""
Basic class used to facilitate the cache system used for data files.

.. include:: ../include/links.rst
"""
from importlib import resources
import pathlib

from IPython import embed

from pypeit import msgs
from pypeit import cache


class PypeItDataPath:
    """
    """
    def __init__(self, subdirs):
        self.data = self.check_isdir(resources.files('pypeit') / 'data')
        self.path = self.check_isdir(self.data / subdirs)

    def glob(self, search_str):
        return self.path.glob(search_str)

    def __repr__(self):
        return str(self.path)
    
    def __truediv__(self, subdir):
        return self.path / subdir

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

    def get_file_path(self, data_file, symlink_in_pkgdir=False, return_format=False):
        """
        Test
        """
        # Make sure the file is a Path object
        _data_file = pathlib.Path(data_file).resolve()

        # Check if the file exists on disk, as provided 
        if _data_file.is_file():
            # If so, assume this points directly to the file to be read
            return _get_file_path_return(_data_file, return_format)

        # Otherwise, construct the file name given the root path:
        _data_file = self.path / data_file

        # If the file exists, return it 
        if _data_file.is_file():
            # Return the full path and the file format
            return _get_file_path_return(_data_file, return_format)

        # If it does not, inform the user and download it into the cache.
        # NOTE: This should not be required for from-source (dev) installations.
        msgs.info(f'{data_file} does not exist in the expected package directory ({data_root}).  '
                    'Checking cache or downloading the file now.')

        # TODO: Fix self.path
        _cached_file = cache.fetch_remote_file(data_file, self.path)

        # If requested, create a symlink to the cached file in the package data
        # directory
        if symlink_in_pkgdir:
            # Create and return values for the symlink
            _data_file.symlink_to(_cached_file)
            return _get_file_path_return(_data_file, return_format)
    
        return _get_file_path_return(_cached_file, return_format)
    

class PypeItDataPaths:
    """
    List of hardwired paths within the pypeit.data module

    Each `@property` method returns a :obj:`pathlib.Path` object
    """
    def __init__(self):
        # Tests
        self.tests = PypeItDataPath('tests')

        # Telluric
        self.telgrid  = PypeItDataPath('telluric/atm_grids')
        self.tel_model = PypeItDataPath('telluric/models')

        # Wavelength Calibrations
        self.reid_arxiv = PypeItDataPath('arc_lines/reid_arxiv')
        self.linelist = PypeItDataPath('arc_lines/lists')
        self.nist = PypeItDataPath('arc_lines/NIST')
        self.arc_plot = PypeItDataPath('arc_lines/plots')

        # Flux Calibrations
        self.standards = PypeItDataPath('standards')
        self.extinction = PypeItDataPath('extinction')
        self.skisim = PypeItDataPath('skisim')
        self.filters = PypeItDataPath('filters')
        self.sensfuncs = PypeItDataPath('sensfuncs')

        # Other
        self.sky_spec = PypeItDataPath('sky_spec')
        self.static_calibs = PypeItDataPath('static_calibs')
        self.spectrographs = PypeItDataPath('spectrographs')


