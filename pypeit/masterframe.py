"""
Implements the master frame base class.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html
"""
import os
import warnings

from abc import ABCMeta

import numpy as np

from astropy.io import fits

from pypeit import msgs

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
            Name of the MasterFrame folder, e.g. MF_keck_deimos. If
            None, set to './'
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

    @property
    def file_name(self):
        """
        Name of the MasterFrame file.
        """
        return 'Master{0}_{1}.{2}'.format(self.master_type, self.master_key, self.file_format)

    @property
    def file_path(self):
        """
        Full path to MasterFrame file.
        """
        return os.path.join(self.master_dir, self.file_name)

    def save(self, data, extnames, outfile=None, overwrite=True, raw_files=None, steps=None):
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
        hdr = fits.Header()
        #   - Set the master frame type
        hdr['FRAMETYP'] = (self.master_type, 'PypeIt: Master calibration frame type')
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
                     ).writeto(_outfile, overwrite=True)

        msgs.info('Master frame written to {0}'.format(_outfile))

    # TODO: have ext default to provide all extensions?
    # TODO: Add a base-level staticmethod one-liner?
    def load(self, ext, ifile=None, return_header=False):
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
        hdu = fits.open(_ifile)
        # Only one extension
        if n_ext == 1:
            data = hdu[_ext[0]].data.astype(np.float)
            return (data, hdu[0].header) if return_header else data
        # Multiple extensions
        data = tuple([ hdu[k].data.astype(np.float) for k in _ext ])
        return data + (hdu[0].header,) if return_header else data

    @staticmethod
    def parse_hdr_raw_files(hdr):
        """
        Parse the file names written to the header by :func:`save`.

        Args:
            hdr (fits.Header):
                Astropy Header object
        
        Returns:
            list: The list of files ordered as they were listed in the
            header.
        """
        raw_files = {}
        prefix = 'F'
        for k, v in hdr.items():
            if k[:len(prefix)] == prefix:
                try:
                    i = int(k[len(prefix):])-1
                except ValueError:
                    continue
                raw_files[i] = v
        # Convert from dictionary with integer keys to an appropriately
        # sorted list
        return [ raw_files[i] for i in range(max(raw_files.keys())+1) ]
    
