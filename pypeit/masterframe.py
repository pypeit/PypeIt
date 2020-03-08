"""
Implements the master frame base class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import inspect
from IPython import embed

from abc import ABCMeta

import numpy as np

from astropy.io import fits

from pypeit import msgs
from pypeit.images import pypeitimage
from pypeit.io import initialize_header
from pypeit.spectrographs import util

import astropy

def construct_file_name(master_obj, master_key, master_dir=None):
    basefile = 'Master{0}_{1}.{2}'.format(master_obj.master_type, master_key,
                                          master_obj.file_format)
    if master_dir is not None:
        basefile =  os.path.join(master_dir, basefile)
    return basefile


def build_master_header(master_obj, master_key, master_dir, spectrograph,
                        hdr=None, steps=None, raw_files=None):
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
    _hdr = initialize_header(hdr)

    # Save the master frame type and key and version, in case the file name is
    # changed.
    _hdr['MSTRTYP'] = (master_obj.master_type, 'PypeIt: Master frame type')
    _hdr['MSTRDIR'] = (master_dir, 'PypeIt: Master directory')
    _hdr['MSTRKEY'] = (master_key, 'PypeIt: Calibration key')
    _hdr['MSTRVER'] = (master_obj.version, 'PypeIt: Master datamodel version')
    #_hdr['MSTRREU'] = (self.reuse_masters, 'PypeIt: Reuse existing masters')

    # Other info, pulled from the Child
    _hdr['PYP_SPEC'] = (spectrograph, 'PypeIt: Spectrograph name')

    #   - List the completed steps
    if steps is not None:
        _hdr['STEPS'] = (','.join(steps), 'Completed reduction steps')
    #   - Provide the file names
    if raw_files is not None:
        nfiles = len(raw_files)
        ndig = int(np.log10(nfiles)) + 1
        for i in range(nfiles):
            _hdr['F{0}'.format(i + 1).zfill(ndig)] = (raw_files[i], 'PypeIt: Processed raw file')
    # Return
    return _hdr


class NewMasterFrame(object):
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
        master_dir (:obj:`str` or None):
            See initialization arguments.
        reuse_masters (:obj:`bool`):
            See initialization arguments.
    """
    def __init__(self, master_dir, reuse_masters, save_masters=True):

        # Output names
        #self.master_dir = os.getcwd() if master_dir is None else master_dir
        self.master_dir = master_dir
        #self.master_key = 'master' if master_key is None else master_key

        # Pulled from master
        #self.master_type = None
        #self.file_format = None

        # Other parameters
        self.reuse_masters = reuse_masters
        self.save_masters = save_masters

        # Other attrib
        self.key_dict = {}

    def construct_file_name(self, master_obj, master_key, add_path=True):
        """
        Name of the MasterFrame file.
        """
        basefile = 'Master{0}_{1}.{2}'.format(master_obj.master_type, master_key,
                                          master_obj.file_format)
        #
        if add_path:
            return os.path.join(self.master_dir, basefile)
        else:
            return basefile

#    @property
#    def file_name(self):
#        """
#        Name of the MasterFrame file.
#        """
#        return MasterFrame.construct_file_name(self.master_type, self.master_key,
#                                               self.file_format)

#    @property
#    def master_file_path(self):
#        """
#        Full path to MasterFrame file.
#        """
#        return os.path.join(self.master_dir, self.file_name)

    @property
    def exists(self):
        """
        Check if the output file already exists.

        Returns:
            bool
        """
        return os.path.isfile(self.master_file_path)

    # TODO: Why doesn't ifile default to None instead of having to provide None?
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
            msgs.warn('You requested that PypeIt not reuse existing masters!')
            return
        if not self.exists:
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
        _hdr = initialize_header(hdr)

        # Save the master frame type and key and version, in case the file name is
        # changed.
        _hdr['MSTRTYP'] = (self.master_type, 'PypeIt: Master frame type')
        _hdr['MSTRDIR'] = (self.master_dir, 'PypeIt: Master directory')
        _hdr['MSTRKEY'] = (self.master_key, 'PypeIt: Calibration key')
        _hdr['MSTRVER'] = (self.master_version, 'PypeIt: Master datamodel version')
        _hdr['MSTRREU'] = (self.reuse_masters, 'PypeIt: Reuse existing masters')

        # Other info, pulled from the Child
        # TODO: I think this should be put in the relevant derived
        # classes, instead of the kludge I had to put in for
        # SlitTraceSet.
        if hasattr(self, 'spectrograph'):
            if isinstance(self.spectrograph, str):
                _hdr['PYP_SPEC'] = (self.spectrograph, 'PypeIt: Spectrograph name')
            else:
                # Assume it's a Spectrograph object
                _hdr['PYP_SPEC'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph name')

        #   - List the completed steps
        if steps is not None:
            _hdr['STEPS'] = (','.join(steps), 'Completed reduction steps')
        #   - Provide the file names
        if raw_files is not None:
            nfiles = len(raw_files)
            ndig = int(np.log10(nfiles)) + 1
            for i in range(nfiles):
                _hdr['F{0}'.format(i + 1).zfill(ndig)] = (raw_files[i], 'PypeIt: Processed raw file')

        return _hdr


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

    # Version
    master_version = None

    # TODO: set master_type and file_format to be class attributes
    # (instead of instance attributes) that each derived class has to
    # define.
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

    @property
    def exists(self):
        """
        Check if the output file already exists.

        Returns:
            bool
        """
        return os.path.isfile(self.master_file_path)

    # TODO: Why doesn't ifile default to None instead of having to provide None?
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
            msgs.warn('You requested that PypeIt not reuse existing masters!')
            return
        if not self.exists:
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
        _hdr = initialize_header(hdr)

        # Save the master frame type and key and version, in case the file name is
        # changed.
        _hdr['MSTRTYP'] = (self.master_type, 'PypeIt: Master frame type')
        _hdr['MSTRDIR'] = (self.master_dir, 'PypeIt: Master directory')
        _hdr['MSTRKEY'] = (self.master_key, 'PypeIt: Calibration key')
        _hdr['MSTRVER'] = (self.master_version, 'PypeIt: Master datamodel version')
        _hdr['MSTRREU'] = (self.reuse_masters, 'PypeIt: Reuse existing masters')

        # Other info, pulled from the Child
        # TODO: I think this should be put in the relevant derived
        # classes, instead of the kludge I had to put in for
        # SlitTraceSet.
        if hasattr(self, 'spectrograph'):
            if isinstance(self.spectrograph, str):
                _hdr['PYP_SPEC'] = (self.spectrograph, 'PypeIt: Spectrograph name')
            else:
                # Assume it's a Spectrograph object
                _hdr['PYP_SPEC'] = (self.spectrograph.spectrograph, 'PypeIt: Spectrograph name')

        #   - List the completed steps
        if steps is not None:
            _hdr['STEPS'] = (','.join(steps), 'Completed reduction steps')
        #   - Provide the file names
        if raw_files is not None:
            nfiles = len(raw_files)
            ndig = int(np.log10(nfiles)) + 1
            for i in range(nfiles):
                _hdr['F{0}'.format(i + 1).zfill(ndig)] = (raw_files[i], 'PypeIt: Processed raw file')

        return _hdr


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


