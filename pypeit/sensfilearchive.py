"""
Provides a class that handles archived sensfunc files.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
from abc import ABC, abstractmethod
import os

from astropy.io import fits

from pypeit import msgs
from pypeit import data

class SensFileArchive(ABC):
    """Class for managing archived SensFunc files. This is an abstract class that instantitates 
    child classes based on spectrograph.
    """

    @abstractmethod
    def get_archived_sensfile(self, fitsfile):
        """Get the full path name of the archived sens file that can be used to flux calibrate a given fitsfile
        
        Args:
            fitsfile (str): The fitsfile to find an archived SensFunc file for.

        Return:
            str: The full pathname of the archived SensFunc.

        Raises:
            PypeItError: Raised an archived SensFunc file can't be found for the given fits file.
        """
        pass

    @classmethod
    def get_instance(cls, spectrograph_name):
        """Return a SensFuncArchive instance that will find archived SensFuncs for a specific spectrograph.
        
        Args:
            spectrograph_name (str): 
                The spectrograph name for the SensFuncArchive instance to return.

        Return:
            pypeit.sensfilearchive.SensFileArchive: 
                A SensFuncArchive object to find archived sensfuncs for a specific spectrograph.

        Raises:
            ValueError: Raised if the passed in spectrograph is not supported.
        """

        for child in cls.__subclasses__():
            if child.spec_name == spectrograph_name:
                return child()

        raise ValueError(f"No SensFileArchive found for {spectrograph_name}")

    @classmethod
    def supported_spectrographs(cls):
        """Return which spectrograph names support Archived SensFuncs.
        
        Return: list of str
        """
        return [child.spec_name for child in cls.__subclasses__()]


class DEIMOSSensFileArchive(SensFileArchive):
    """SensFileArchive subclass specifically for keck_deimos SensFuncs."""
    spec_name = "keck_deimos"

    def get_archived_sensfile(self, fitsfile, symlink_in_pkgdir=False):
        """Get the full path name of the archived sens file that can be used to flux calibrate a given fitsfile
        
        Args:
            fitsfile (str): The fitsfile to find an archived SensFunc file for.
            symlink_in_pkgdir (bool): Create a symlink to the the cached file in the package directory (default False)

        Return:
            str: The full pathname of the archived SensFunc.

        Raises:
            PypeItError: Raised an archived SensFunc file can't be found for the given fits file.
        """
        header = fits.getheader(fitsfile)
        
        grating = header['DISPNAME']
        if grating not in ["600ZD", "830G", "900ZD", "1200B", "1200G"]:
            msgs.error(f"There are no archived SensFuncFiles for keck_deimos grating {grating}.")
        
        archived_file = data.get_sensfunc_filepath(f"keck_deimos_{grating}_sensfunc.fits",
                                                   symlink_in_pkgdir=symlink_in_pkgdir)
        msgs.info(f"Found archived sensfile '{archived_file}'")
        return archived_file

