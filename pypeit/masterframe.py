"""
Implements the master frame base class.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit.io import initialize_header

import astropy

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
    basefile = 'Master{0}_{1}.{2}'.format(master_obj.master_type, master_key,
                                          master_obj.file_format)
    filename = os.path.join(master_dir, basefile) if master_dir is not None else basefile
    # Return
    return filename


def build_master_header(master_obj, master_key, master_dir, spectrograph,
                        hdr=None, steps=None, raw_files=None):
    """
    Initialize the master frame header.

    The function writes information generic to all PypeIt master
    frame headers with basic information.

    Args:
        master_obj (object):
            MasterFrame object to be named.  This provides the master_type and file_format
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
