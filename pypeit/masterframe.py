"""
Implements the master frame base class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
from IPython import embed

from abc import ABCMeta

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit.images import pypeitimage
from pypeit.io import initialize_header
from pypeit.spectrographs import util

import astropy

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
    def master_file_path(self):
        """
        Full path to MasterFrame file.
        """
        return os.path.join(self.master_dir, self.file_name)

    def chk_load_master(self, ifile):
        """
        Simple check to see if reuse_masters is set
        and whether the master file exists.

        If both of these are true, it returns the filename
        otherwise it returns None.

        Args:
            ifile (:obj:`str` or None):
                Name of the master frame file one is considering.
                Defaults to :attr:`master_file_path`.

        Returns:
            None or str:
        """
        # Format the input and set the tuple for an empty return
        _ifile = self.master_file_path if ifile is None else ifile
        if not self.reuse_masters:
            # User does not want to load masters
            msgs.warn('PypeIt will not reuse masters!')
            return
        if not os.path.isfile(_ifile):
            # Master file doesn't exist
            msgs.warn('No Master {0} frame found: {1}'.format(self.master_type, self.master_file_path))
            return
        # Else return the file
        return _ifile


    def build_master_header(self, hdr=None, steps=None, raw_files=None):
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
            steps (:obj:`list`, optional):
                The list of steps executed by the derived class to
                construct the master frame.
            raw_files (:obj:`list`, optional):
                List of processed raw files used to construct the master
                frame.

        Returns:
            `astropy.io.fits.Header`: The initialized (or edited)
            fits header.
        """
        # Standard init
        hdr = initialize_header(hdr)

        # Save the master frame type and key, in case the file name is
        # changed.
        hdr['MSTRTYP'] = (self.master_type, 'PypeIt: Master frame type')
        hdr['MSTRDIR'] = (self.master_dir, 'PypeIt: Master directory')
        hdr['MSTRKEY'] = (self.master_key, 'PypeIt: Calibration key')
        hdr['MSTRREU'] = (self.reuse_masters, 'PypeIt: Reuse existing masters')

        # Other info, pulled from the Child
        if hasattr(self, 'spectrograph'):
            hdr['PYP_SPEC'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph name')

        #   - List the completed steps
        if steps is not None:
            hdr['STEPS'] = (','.join(steps), 'Completed reduction steps')
        #   - Provide the file names
        if raw_files is not None:
            nfiles = len(raw_files)
            ndig = int(np.log10(nfiles)) + 1
            for i in range(nfiles):
                hdr['F{0}'.format(i + 1).zfill(ndig)] = (raw_files[i], 'PypeIt: Processed raw file')

        return hdr


def items_from_master_file(master_file):
    """
    Grab items from the Master file.  In particular, generate the
    spectrograph object

    Either the header or some other part of the object

    Args:
        master_file (str):
            Full path to the file

    Returns:
        tuple: Returns a :class:`pypeit.spectrograph.Spectrograph`
        instance generated from information in the the master object and
        a :obj:`list` of extra objects, currently the header from the
        master file (developer note: avoid having extras this grow out
        of control!)
    """
    ext = master_file.split('.')[-1]

    if ext == 'fits':
        # Spectrograph from header
        hdu = fits.open(master_file)
        head0 = hdu[0].header
        spec_name = head0['PYP_SPEC']
        spectrograph = util.load_spectrograph(spec_name)
        extras = [head0]
    else:
        msgs.error("Not read for this type of master file")
    # Return
    return spectrograph, extras


