"""
Implements the master frame base class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
from IPython import embed
from abc import ABCMeta

import numpy as np

from pypeit import msgs
from pypeit.io import initialize_header

from astropy.io import fits

# DEFINE HERE AS THE DATAMODEL
sep1 = '_'  # Separation between master type and key
sep2 = '.'  # Separation between master key and extension

def construct_file_name(master_obj, master_key, master_dir=None):
    """
    Generate a MasterFrame filename
    Args:
        master_obj (object):
            MasterFrame object to be named.  This provides the master_type and file_format
        master_key (str):
            Designation
        master_dir (str, optional):
            Path to the master frame folder
    Returns:
        str:
    """
    basefile = 'Master{0}{1}{2}{3}{4}'.format(master_obj.master_type, sep1,
                                              master_key, sep2,
                                              master_obj.master_file_format)
    filename = os.path.join(master_dir, basefile) if master_dir is not None else basefile
    # Return
    return filename


def grab_key_mdir(inp, from_filename=False):
    """
    Grab master_key and master_dir by parsing a filename or inspecting a header
    Args:
        inp (:obj:`str` or astropy.io.fits.Header):
            Either a filename or a Header of a FITS file
        from_filename (bool, optional):
            If true, parse the input filename using the naming model
    Returns:
        tuple:  str, str of master_key and master_dir
    """
    if from_filename:
        # Grab the last folder of the path
        master_dir = os.path.dirname(inp)
        # Parse
        base = os.path.basename(inp)
        pos1 = base.find(sep1)
        pos2 = base.find(sep2)
        master_key = base[pos1+1:pos2]
    else:
        if isinstance(inp, str):
            head0 = fits.getheader(inp)
        elif isinstance(inp, fits.Header):
            head0 = inp
        # Grab it
        master_key = head0['MSTRKEY'] if 'MSTRKEY' in head0.keys() else None
        master_dir = head0['MSTRDIR'] if 'MSTRDIR' in head0.keys() else None
    # Return
    return master_key, master_dir


def build_master_header(master_obj, master_key, master_dir,
                        hdr=None, steps=None, raw_files=None):
    """
    Initialize the master frame header.

    This builds a generic header that is written to all PypeIt master
    frames.

    Args:
        master_obj (object):
            MasterFrame object to be named. This provides the
            master_type and file_format
        master_key (str):
            Designation
        master_dir (str):
            Path to the master frame folder
        spectrograph (str):
            Name of the spectrograph
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

    # Spectrograph
    if master_obj.PYP_SPEC is None:
        msgs.error("The object needs to include PYP_SPEC this so that it was written to the Header")
    _hdr['PYP_SPEC'] = (master_obj.PYP_SPEC, 'PypeIt: Spectrograph name')  # This may be over-written by itself

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

