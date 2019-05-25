"""
Implements the master frame base class.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""
import os
import sys
import warnings
import time

from abc import ABCMeta

import numpy as np

from astropy.io import fits

from pypeit import msgs

# These imports are largely just to make the versions available for
# writing to the header. See `initialize_header`
import scipy
import astropy
import sklearn
import pypeit

class MasterFrame(object):
    """
    Base class for all PypeIt calibration master frames.

    Main methods of the class are to implement the default paths for the
    master frames and their I/O methods.

    Args:
        master_type (str):
            Master frame type, used when constructing the output file
            name.  See :attr:`file_name`.
        master_dir (str or None):
            Name of the MasterFrame folder. If None, set to current
            working directory.
        master_key (str or None):
            Root of the MasterFrame names, e.g. 'A_1_01'.  If None, set
            to 'master'
        file_format (:obj:`str`, optional):
            File format for the master frame.  Typically 'fits' or 'json'.
        reuse_masters (bool, optional):
            Reuse already created master files from disk.

    Attributes:
        master_type (:obj:`str`):
            See initialization arguments.
        master_dir (:obj:`str`):
            See initialization arguments.
        master_key (:obj:`str`):
            See initialization arguments.
        file_format (:obj:`str`):
            See initialization arguments.
        reuse_masters (:obj:`bool`):
            See initialization arguments.
    """
    __metaclass__ = ABCMeta

    def __init__(self, master_type, master_dir=None, master_key=None, file_format='fits',
                 reuse_masters=False):

        # Output names
        self.master_type = master_type
        self.master_dir = os.getcwd() if master_dir is None else master_dir
        self.master_key = 'master' if master_key is None else master_key
        self.file_format = file_format

        # Other parameters
        self.reuse_masters = reuse_masters

    # TODO: Kludge to allow for a one-liner that gets the master frame
    # file name
    @staticmethod
    def construct_file_name(master_type, master_key, file_format='fits'):
        """
        Name of the MasterFrame file.
        """
        return 'Master{0}_{1}.{2}'.format(master_type, master_key, file_format)

    @property
    def file_name(self):
        """
        Name of the MasterFrame file.
        """
        return MasterFrame.construct_file_name(self.master_type, self.master_key,
                                               self.file_format)

    @property
    def file_path(self):
        """
        Full path to MasterFrame file.
        """
        return os.path.join(self.master_dir, self.file_name)

    def save(self, data, extnames, outfile=None, overwrite=True, raw_files=None, steps=None,
             checksum=True):
        """
        Base interface for saving master frame image data to a fits
        file.

        Data is always placed in extensions, with the PRIMARY extension
        only containing the header data.  More complicated master frame
        data models should overload this function.

        Args:
            data (:obj:`list`, `numpy.ndarray`_):
                One or more data arrays to save.
            extnames (:obj:`list`, :obj:`str`):
                The names for the data extensions.
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
            raw_files (:obj:`list`, optional):
                List of processed raw files used to construct the master
                frame.
            steps (:obj:`list`, optional):
                The list of steps executed by the derived class to
                construct the master frame.
            checksum (:obj:`bool`, optional):
                Passed to `astropy.io.fits.HDUList.writeto` to add
                the DATASUM and CHECKSUM keywords fits header(s).
        """
        _outfile = self.file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Log
        msgs.info('Saving master frame to {0}'.format(_outfile))

        # Format the output
        ext = extnames if isinstance(extnames, list) else [extnames]
        if len(ext) > 1 and not isinstance(data, list):
            msgs.error('Input data type should be list, one numpy.ndarray per extension.')
        _data = data if isinstance(data, list) else [data]

        # Build the header
        hdr = self.initialize_header()
        #   - List the completed steps
        if steps is not None:
            hdr['STEPS'] = (','.join(steps), 'Completed reduction steps')
        #   - Provide the file names
        if raw_files is not None:
            nfiles = len(raw_files)
            ndig = int(np.log10(nfiles))+1
            for i in range(nfiles):
                hdr['F{0}'.format(i+1).zfill(ndig)] = (raw_files[i], 'PypeIt: Processed raw file')
        
        # Write the fits file
        fits.HDUList([fits.PrimaryHDU(header=hdr)]
                        + [ fits.ImageHDU(data=d, name=n) for d,n in zip(_data, ext)]
                     ).writeto(_outfile, overwrite=True, checksum=checksum)

        msgs.info('Master frame written to {0}'.format(_outfile))

    # TODO: Add a base-level staticmethod one-liner?
    # TODO: include checksum keyword, used to validate data when
    # loading?
    def load_master(self, ext, ifile=None, return_header=False):
        """
        Generic master file reader.

        This generic reader assumes the file is in fits format.
        
        Args:
            ext (:obj:`str`, :obj:`int`, :obj:`list`):
                One or more file extensions with data to return.  The
                extension can be designated by its 0-indexed integer
                number or its name.
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            tuple: Returns the image data from each provided extension.
            If return_header is true, the primary header is also
            returned.  If nothing is loaded, either because
            :attr:`reuse_masters` is `False` or the file does not exist,
            everything is returned as None (one per expected extension
            and one if the header is requested).
        """
        # Format the input and set the tuple for an empty return
        _ifile = self.file_path if ifile is None else ifile
        _ext = ext if isinstance(ext, list) else [ext]
        n_ext = len(_ext)
        empty_return = ((None,)*(n_ext+1) if return_header else (None,)*n_ext) \
                            if n_ext > 1 or return_header else None
        if not self.reuse_masters:
            # User does not want to load masters
            msgs.warn('PypeIt will not reuse masters!')
            return empty_return
        
        if not os.path.isfile(_ifile):
            # Master file doesn't exist
            msgs.warn('No Master {0} frame found: {1}'.format(self.master_type, self.file_path))
            return empty_return

        # Read and return
        msgs.info('Loading Master {0} frame: {1}'.format(self.master_type, _ifile))
        return MasterFrame.load_from_file(_ifile, _ext, return_header=return_header)

    # TODO: include checksum keyword, used to validate data when
    # loading?
    @staticmethod
    def load_from_file(filename, ext, return_header=False):
        """
        Generic static method to read a master file.

        Args:
            filename (:obj:`str`):
                Name of the master frame file.
            ext (:obj:`str`, :obj:`int`, :obj:`list`):
                One or more file extensions with data to return.  The
                extension can be designated by its 0-indexed integer
                number or its name.
            return_header (:obj:`bool`, optional):
                Return the header.

        Returns:
            tuple: Returns the image data from each provided extension.
            If return_header is true, the primary header is also
            returned.
        """
        # Check file exists
        if not os.path.isfile(filename):
            msgs.error('File does not exist: {0}'.format(filename))
        
        # Format the input and set the tuple for an empty return
        _ext = ext if isinstance(ext, list) else [ext]
        n_ext = len(_ext)

        # Open the file
        hdu = fits.open(filename)

        # Only one extension
        if n_ext == 1:
            data = hdu[_ext[0]].data.astype(np.float)
            return (data, hdu[0].header) if return_header else data
        # Multiple extensions
        data = tuple([None if hdu[k].data is None else hdu[k].data.astype(np.float) for k in _ext ])
        return data + (hdu[0].header,) if return_header else data


    def initialize_header(self, hdr=None):
        """
        Initialize the master frame header.

        The function writes information generic to all PypeIt master
        frame headers with basic information.

        Args:
            hdr (`astropy.io.fits.Header`, optional):
                Header object to update with basic summary
                information. The object is modified in-place and also
                returned. If None, an empty header is instantiated,
                edited, and returned.

        Returns:
            `astropy.io.fits.Header`: The initialized (or edited)
            fits header.
        """
        # Add versioning; this hits the highlights but should it add
        # the versions of all packages included in the requirements.txt
        # file?
        if hdr is None:
            hdr = fits.Header()
        hdr['VERSPYT'] = ('.'.join([ str(v) for v in sys.version_info[:3]]), 'Python version')
        hdr['VERSNPY'] = (np.__version__, 'Numpy version')
        hdr['VERSSCI'] = (scipy.__version__, 'Scipy version')
        hdr['VERSAST'] = (astropy.__version__, 'Astropy version')
        hdr['VERSSKL'] = (sklearn.__version__, 'Scikit-learn version')
        hdr['VERSPYP'] = (pypeit.__version__, 'PypeIt version')

        # Save the master frame type and key, in case the file name is
        # changed.
        hdr['MSTRTYP'] = (self.master_type, 'PypeIt: Master frame type')
        hdr['MSTRDIR'] = (self.master_dir, 'PypeIt: Master directory')
        hdr['MSTRKEY'] = (self.master_key, 'PypeIt: Calibration key')
        hdr['MSTRREU'] = (self.reuse_masters, 'PypeIt: Reuse existing masters')

        # Save the date of the reduction
        hdr['DATE'] = (time.strftime('%Y-%m-%d',time.gmtime()), 'UTC date created')

        # TODO: Anything else?

        return hdr
