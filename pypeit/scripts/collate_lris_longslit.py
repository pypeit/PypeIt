"""
This script collates LRIS long slit 1d spectra by filename and slit position, 
runs flux calibration on them, and then coadds them.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from datetime import datetime
import os.path
from functools import partial
import traceback
import sys
import shutil

import numpy as np
from astropy.coordinates import Angle
from astropy.time import Time
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import coadd1d
from pypeit import msgs
from pypeit import par
from pypeit.utils import is_float
from pypeit.core import wave
from pypeit.archive import ArchiveMetadata, ArchiveDir
from pypeit.core.collate import collate_spectra_by_source, SourceObject
from pypeit.scripts import scriptbase
from pypeit.slittrace import SlitTraceBitMask
from pypeit.spec2dobj import AllSpec2DObj
from pypeit.sensfilearchive import SensFileArchive
from pypeit import fluxcalibrate
from pypeit.specobjs import SpecObjs
from pypeit import inputfiles


def find_slits_to_exclude(spec2d_files, par):
    """
    Find slits that should be excluded according to the input parameters.

    The slit mask ids are returned in a map alongside the text labels for the
    flags that caused the slit to be excluded.

    Args:
        spec2d_files (:obj:`list`): 
            List of spec2d files to build the map from.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Parameters from a ``.collate1d`` file

    Returns:
        :obj:`dict`: Mapping of slit mask ids to the flags that caused the slit
        to be excluded.
    """

    # Get the types of slits to exclude from our parameters
    exclude_flags = par['collate1d']['exclude_slit_trace_bm']
    if isinstance(exclude_flags, str):
        exclude_flags = [exclude_flags]

    # Go through the slit_info of all spec2d files and find
    # which slits should be excluded based on their flags
    bit_mask = SlitTraceBitMask()
    exclude_map = dict()
    for spec2d_file in spec2d_files:

        allspec2d = AllSpec2DObj.from_fits(spec2d_file, chk_version=par['rdx']['chk_version'])
        for sobj2d in [allspec2d[det] for det in allspec2d.detectors]:
            for (slit_id, mask, slit_mask_id) in sobj2d['slits'].slit_info:
                for flag in exclude_flags:
                    if bit_mask.flagged(mask, flag=flag):
                        if slit_mask_id not in exclude_map:
                            exclude_map[slit_mask_id] = {flag}
                        else:
                            exclude_map[slit_mask_id].add(flag)

    return exclude_map

def read_spec1d_files(par, spec1d_files, failure_msgs):
    """
    Read spec1d files.

    Args:
        par (`obj`:pypeit.par.pypeitpar.PypeItPar): 
            Parameters for collating, fluxing, and coadding.
        spec1d_files (list of str):
            List of spec1d files to read.
        failure_msgs(list of str):
            Return parameter describing any failures that occurred when reading.

    Returns:
        list of str: The SpecObjs objects that were successfully read.
        list of str: The spec1d files that were successfully read.
    """

    specobjs_list = []
    good_spec1d_files = []
    for spec1d_file in spec1d_files:
        try:
            sobjs = SpecObjs.from_fitsfile(spec1d_file, chk_version=par['rdx']['chk_version'])
            specobjs_list.append(sobjs)
            good_spec1d_files.append(spec1d_file)
        except Exception as e:
            formatted_exception = traceback.format_exc()
            msgs.warn(formatted_exception)
            msgs.warn(f"Failed to read {spec1d_file}, skipping it.")
            failure_msgs.append(f"Failed to read {spec1d_file}, skipping it.")
            failure_msgs.append(formatted_exception)

    return specobjs_list, good_spec1d_files
