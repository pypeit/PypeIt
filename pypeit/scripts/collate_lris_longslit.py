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

def build_parameters(args):
    """
    Read the command-line arguments and the input ``.collate1d`` file (if any), 
    to build the parameters needed by ``collate_1d``.

    Args:
        args (`argparse.Namespace`_):
            The parsed command line as returned by the ``argparse`` module.

    Returns:
        :obj:`tuple`: Returns three objects: a
        :class:`~pypeit.par.pypeitpar.PypeItPar` instance with the parameters
        for collate_1d, a
        :class:`~pypeit.spectrographs.spectrograph.Spectrograph` instance with
        the spectrograph parameters used to take the data, and a :obj:`list`
        with the spec1d files read from the command line or ``.collate1d`` file.
    """
    # First we need to get the list of spec1d files
    if args.input_file is not None:
        collateFile = inputfiles.Collate1DFile.from_file(args.input_file)
        cfg_lines, spec1d_files = collateFile.cfg_lines, collateFile.filenames

        # Look for a coadd1d file
        input_file_root, input_file_ext = os.path.splitext(args.input_file)
        coadd1d_config_name = input_file_root + ".coadd1d"
        if os.path.exists(coadd1d_config_name):
            coadd1DFile = inputfiles.Coadd1DFile.from_file(coadd1d_config_name)
            cfg_lines += coadd1DFile.cfg_lines

    else:
        cfg_lines = None
        spec1d_files = []

    if args.spec1d_files is not None and len(args.spec1d_files) > 0:
        spec1d_files = args.spec1d_files

    if spec1d_files is None or len(spec1d_files) == 0:
        parser = Collate1D.get_parser()
        print("Missing arguments: A list of spec1d files must be specified via command line or config file.")
        parser.print_usage()
        sys.exit(1)

    # Get the spectrograph for these files and then create a ParSet. 
    spectrograph = load_spectrograph(spec1d_files[0])
    spectrograph_def_par = spectrograph.default_pypeit_par()

    if cfg_lines is not None:
        # Build using config file
        params = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), 
                                                    merge_with=(cfg_lines,))
    else:
        # No config file, use the defaults and supplement with command line args
        params = spectrograph_def_par
        params['collate1d'] = pypeitpar.Collate1DPar()

    # command line arguments take precedence over config file parameters
    if args.tolerance is not None:
        params['collate1d']['tolerance'] = args.tolerance

    if args.match_using is not None:
        params['collate1d']['match_using'] = args.match_using

    if args.exclude_slit_trace_bm is not None and len(args.exclude_slit_trace_bm) > 0:
        params['collate1d']['exclude_slit_trace_bm'] = args.exclude_slit_trace_bm.split(',')

    if args.exclude_serendip:
        params['collate1d']['exclude_serendip'] = True

    if args.wv_rms_thresh is not None:
        params['collate1d']['wv_rms_thresh'] = args.wv_rms_thresh

    if args.dry_run:
        params['collate1d']['dry_run'] = True

    if args.ignore_flux:
        params['collate1d']['ignore_flux'] = True

    if args.flux:
        params['collate1d']['flux'] = True

    if args.outdir is not None:
        params['collate1d']['outdir'] = args.outdir

    if args.spec1d_outdir is not None:
        params['collate1d']['spec1d_outdir'] = args.spec1d_outdir

    if args.refframe is not None:
        params['collate1d']['refframe'] = args.refframe

    if args.chk_version is True:
        params['rdx']['chk_version'] = True

    return params, spectrograph, spec1d_files
