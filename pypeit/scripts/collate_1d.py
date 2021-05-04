#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script collates multiple 1d spectra in multiple files by object, 
runs flux calibration/coadding on them, and produces files suitable
for KOA archiving.
"""

from datetime import datetime
import argparse
from glob import glob
import os.path
from functools import partial
import re

import numpy as np
from astropy.coordinates import Angle
from astropy.io import fits

from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import coadd1d
from pypeit import msgs
from pypeit import par
from pypeit.utils import is_float
from pypeit.archive.archive_dir import ArchiveDir
from pypeit.archive.archive_metadata import ArchiveMetadata
from pypeit.core.collate import group_spectra_by_source, find_slits_to_exclude, SourceObject

# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)



def extract_id(header):
    """
    Pull an id from a file's header.

    This will give preference to a KOAID, but will return an id based on the file name if a KOAID can't be found.
    A KOAID is of the format II.YYYYMMDD.xxxxx.fits See the `KOA FAQ <https://www2.keck.hawaii.edu/koa/public/faq/koa_faq.php>`_ 
    for more information.

    Args:
        header (str):   A fits file header.

    Returns:
        str: The an id extracted from the header.
    """

    # First check for the KOAID keyword

    if 'KOAID' in header:
        return header['KOAID']
    else:
        # Attempt to pull KOAID from file name
        filename = header['FILENAME']
        if len(filename) >= 17:
            koaid = filename[0:17]
            if re.match(r'..\.\d{8}\.\d{5}$', koaid) is not None:
                # KOA seems to append .fits to the ID
                return koaid + ".fits"

        # For non KOA products, we use the filename
        return filename

def get_metadata_by_id(header_keys, file_info):
    """
    Gets the metadata from a FITS header used for the by id portion
    of the archive and appends it to by_id_metadata

    Args:
        filename (str): A filename from a file originating in the KOA.
    
    """
    if isinstance(file_info, SourceObject):
        return (None, None, None)


    filename = os.path.basename(file_info)

    header = fits.getheader(file_info)

    # Extract koa id from source image filename in header
    id = extract_id(header)

    # Build data row, which starts with koaid and filename + the metadata
    data_row = [id, filename] + [None if x not in header else header[x] for x in header_keys]

    return ([data_row], file_info, filename)

def get_object_based_metadata(object_header_keys, spec_obj_keys, source_header_keys, file_info):
    """
    Gets the metadata from a SourceObject instance used for the by object
    portion of the archive and appends it to by_object_metadata.

    Args:
    source (:obj:`pypeit.scripts.collate_1d.SourceObject`)): The source object containing the
        headers, filenames and SpecObj information for a coadd output file.
    """

    if not isinstance(file_info, SourceObject):
        return (None, None, None)

    header_data = [file_info.spec1d_header_list[0][x] for x in object_header_keys]
    spec_obj_data = [file_info.spec_obj_list[0][x] for x in spec_obj_keys]
    shared_data = [os.path.basename(file_info.coaddfile)] + spec_obj_data + header_data
    result_rows = []
    for header in file_info.spec1d_header_list:
        id = extract_id(header)
        unique_data = [header[x] if x in header else None for x in source_header_keys]
        result_rows.append(shared_data + [id] + unique_data)

    return (result_rows, file_info.coaddfile, file_info.coaddfile)

def coadd(par, source):
    """coadd the spectra for a given source.

    Args:
        par (`obj`:Collate1DPar): Paramters for the coadding
        source (`obj`:SourceObject): The SourceObject with information on
            which files and spectra to coadd.
    """
    par['coadd1d']['coaddfile'] = source.coaddfile
    par['coadd1d']['flux_value'] = False
    spectrograph = load_spectrograph(par['rdx']['spectrograph'])

    # Instantiate
    coAdd1d = coadd1d.CoAdd1D.get_instance(source.spec1d_file_list,
                                           [x.NAME for x in source.spec_obj_list],
                                           spectrograph=spectrograph, par=par['coadd1d'])

    # Run
    coAdd1d.run()
    # Save to file
    coAdd1d.save(source.coaddfile)

def read_file_list(lines, block):
    """
    Read a list of lines from the "read" portion of a PypeIt config file.

    Reads from a block of file names in a configuration file such as:

    spec1d read
    Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits
    spec1d end


    Args:
        lines (:obj:`numpy.ndarray`): List of lines to read.
        block (str): Name of the block to read from (e.g. 'spec1d')

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values, i.e. lines that did not
          contain the "read" list.

        files (list):
          Contains the list of lines read.

    Raises:
        ValueError if there was no list to read or there was a syntax error with the "read" portion.
    """

    files = []
    is_config = np.ones(len(lines), dtype=bool)

    s, e = par.util._find_pypeit_block(lines, block)
    if s >= 0 and e < 0:
        raise ValueError(f"Missing '{block} end' in collate1d config file")
    elif (s < 0) or (s==e):
        raise ValueError(f"Missing {block} block in collate1d config file")
    else:
        files = lines[s:e]

    is_config[s-1:e+1] = False

    return is_config, files


def read_coadd1d_config(coadd1d_file):
    """Read coadd1d configuration from a file.

    This will read any configuration keys, but skip the
    coadd1d section that contains files and objects to coadd.
    That information is generated by group_spectra_by_source
    instead.

    Args:
        coadd1d_file (str): 
            Path name of coadd1d file to read

    Return:
        cfg_lines (list): 
            Config lines to modify ParSet values

    """

    lines = par.util._read_pypeit_file_lines(coadd1d_file)

    # Strip out the coadd list
    try:
        is_config, coadd_lines = read_file_list(lines, 'coadd1d')
        cfg_lines = list(lines[is_config])
    except ValueError as e:
        # The coadd1d files were missing, which is okay because they
        # are ignored when collating
        cfg_lines = list(lines)

    return cfg_lines


def find_spec2d_from_spec1d(spec1d_files):
    """
    Find the spec2d files corresponding to the given list of spec1d files.
    This looks for the spec2d files in  the same directory as the spec1d files.
    It will exit with an error if a spec2d file cannot be found.

    Args:
    spec1d_files (list of str): List of spec1d files generated by PypeIt.

    Returns:
    list of str: List of the matching spec2d files.
    """

    spec2d_files = []
    for spec1d_file in spec1d_files:
        # Check for a corresponding 2d file
        (path, filename) = os.path.split(spec1d_file)
        spec2d_file = os.path.join(path, filename.replace('spec1d', 'spec2d', 1))

        if not os.path.exists(spec2d_file):
            msgs.error(f'Could not find matching spec2d file for {spec1d_file}')

        spec2d_files.append(spec2d_file)

    return spec2d_files

def build_parameters(args):
    """
    Read the command line arguments and the input .collate1d file (if any), 
    to build the parameters needed by collate_1d.

    Args:
    args (:obj:`argparse.Namespace`): The parsed command line as returned
        by the argparse module.

    Returns:
    :obj:`pypeit.par.pypeitpar.PypeItPar`: 
        The parameters for collate_1d.

    :obj:`pypeit.spectrographs.spectrograph.Spectrograph`:
        The spectrograph for the given spec1d files.

    list of 'str': The spec1d files read from the command line or .collate1d file.
    """
    # First we need to get the list of spec1d files
    if args.input_file is not None:
        (cfg_lines, spec1d_files) = par.util.parse_tool_config(args.input_file, 'spec1d', check_files=True)

        # Look for a coadd1d file
        (input_file_root, input_file_ext) = os.path.splitext(args.input_file)
        coadd1d_config_name = input_file_root + ".coadd1d"
        if os.path.exists(coadd1d_config_name):
            cfg_lines += par.util.parse_tool_config(coadd1d_config_name, 'coadd1d')[0]

    else:
        cfg_lines = None
        spec1d_files = []

    if args.spec1d_files is not None and len(args.spec1d_files) > 0:
        spec1d_files = args.spec1d_files

    if spec1d_files is None or len(spec1d_files) == 0:
        msgs.error("A list of spec1d files must be specified via command line or config file.")

    # Get the spectrograph for these files and then create a ParSet. 
    spectrograph = load_spectrograph(spec1d_files[0])
    spectrograph_def_par = spectrograph.default_pypeit_par()

    if cfg_lines is not None:
        # Build using config file
        params = pypeitpar.PypeItPar.from_cfg_lines(cfg_lines=spectrograph_def_par.to_config(), merge_with=cfg_lines)
    else:
        # No config file, use the defaults and supplement with command line args
        params = spectrograph_def_par
        params['collate1d'] = pypeitpar.Collate1DPar()

    # command line arguments take precedence over config file parameters
    if args.tolerance is not None:
        params['collate1d']['tolerance'] = args.tolerance

    if args.match is not None:
        params['collate1d']['match_using'] = args.match

    if args.exclude_slit is not None and len(args.exclude_slit) > 0:
        params['collate1d']['slit_exclude_flags'] = args.exclude_slit

    if args.dry_run:
        params['collate1d']['dry_run'] = True

    if args.archive_dir is not None:
        params['collate1d']['archive_root'] = args.archive_dir

    return params, spectrograph, spec1d_files

def parse_args(options=None, return_parser=False):
    """Parse the command line arguments"""

    # A blank Colate1DPar to avoid duplicating the help text.
    blank_par = pypeitpar.Collate1DPar()

    parser = argparse.ArgumentParser(description='Flux/Coadd multiple 1d spectra from multiple nights and prepare a directory for the KOA.',
                                     formatter_class=SmartFormatter)

    parser.add_argument('input_file', type=str,
                        help='R|(Optional) File for guiding the collate process.\n'
                             'Parameters in this file are overidden by the command\n'
                             'line. The file must have the following format:\n'
                             '\n'
                             '[collate1d]\n'
                             '  tolerance          <tolerance>\n'
                             '  archive_root       <directory for archive files>\n'
                             '  slit_exclude_flags <slit types to exclude>\n'
                             '  match_using        Whether to match using "pixel" or\n'
                             '                     "ra/dec"\n'
                             '  dry_run            If set the matches are displayed\n'
                             '                     without any processing\n'
                             '\n'
                             'spec1d read\n'
                             '<path to spec1d files, wildcards allowed>\n'
                             '...\n'
                             'end\n',                        
                        nargs='?')
    parser.add_argument('--spec1d_files', type=str, nargs='*', help='One or more spec1d files to ' \
                        'flux/coadd/archive. Can contain wildcards')
    parser.add_argument('--par_outfile', default='collate1d.par', type=str,
                        help='Output to save the parameters')
    parser.add_argument('--tolerance', type=str, help=blank_par.descr['tolerance'])
    parser.add_argument('--match', type=str, choices=blank_par.options['match_using'], help=blank_par.descr['match_using'])
    parser.add_argument('--dry_run', action='store_true', help=blank_par.descr['dry_run'])
    parser.add_argument('--archive_dir', type=str, help=blank_par.descr['archive_root'])
    parser.add_argument('--exclude_slit', type=str, nargs='*', help=blank_par.descr['slit_exclude_flags'])

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):

    start_time = datetime.now()
    (par, spectrograph, spec1d_files) = build_parameters(args)

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)

    # Make sure archive dir, if specified, exists
    if par['collate1d']['archive_root'] is not None:
        os.makedirs(par['collate1d']['archive_root'], exist_ok=True)

    if len(par['collate1d']['slit_exclude_flags']) > 0:
        spec2d_files = find_spec2d_from_spec1d(spec1d_files)
        exclude_map = find_slits_to_exclude(spec2d_files, par)
    else:
        spec2d_files = []
        exclude_map = dict()

    if par['collate1d']['match_using'] == 'pixel':
        tolerance = float(par['collate1d']['tolerance'])
    else:
        # For ra/dec matching, the default unit is arcseconds. We check for
        # this case by seeing if the passed in tolerance is a floating point number
        if is_float(par['collate1d']['tolerance']):
            tolerance =  float(par['collate1d']['tolerance'])
        else:
            tolerance = Angle(par['collate1d']['tolerance']).arcsec

    source_list = group_spectra_by_source(spec1d_files, exclude_map, par['collate1d']['match_using'], tolerance)

    #sensfunc, how to identify standard file

    # fluxing etc goes here

    for source in source_list:

        msgs.info(f'Creating {source.coaddfile} from the following sources:')
        for i in range(len(source.spec_obj_list)):
            msgs.info(f'    {source.spec1d_file_list[i]}: {source.spec_obj_list[i].NAME} ({source.spec_obj_list[i].MASKDEF_OBJNAME})')

        if not args.dry_run:
            coadd(par, source)

    if not args.dry_run:

        if par['collate1d']['archive_root'] is not None:
            metadata_root = par['collate1d']['archive_root']
            copy = True
        else:
            metadata_root = os.getcwd()
            copy = False

        archive = create_archive(metadata_root, copy)
        archive.add(spec1d_files)
        archive.add(spec2d_files)
        archive.add(source_list)
        archive.save()

    total_time = datetime.now() - start_time

    msgs.info(f'Total duration: {total_time}')

    return 0

def create_archive(archive_root, copy_to_archive):


    ID_BASED_HEADER_KEYS  = ['RA', 'DEC', 'TARGET', 'PJROGPI', 'SEMESTER', 'PROGID', 'DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
    OBJECT_BASED_HEADER_KEYS = ['DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
    OBJECT_BASED_SPEC_KEYS   = ['MASKDEF_OBJNAME', 'MASKDEF_ID', 'DET', 'RA', 'DEC']
    OBJECT_BASED_SOURCE_KEYS = ['GUIDFWHM', 'PJROGPI', 'SEMESTER', 'PROGID']

    by_id_names = ['id', 'filename'] + [x.lower() for x in ID_BASED_HEADER_KEYS]
    by_id_metadata = ArchiveMetadata(os.path.join(archive_root, "by_id_meta.dat"), 
                                        by_id_names, 
                                        partial(get_metadata_by_id, ID_BASED_HEADER_KEYS),
                                        append=True)

    by_object_names = ['filename'] + \
                        [x.lower() for x in OBJECT_BASED_SPEC_KEYS] + \
                        [x.lower() for x in OBJECT_BASED_HEADER_KEYS] + \
                        ['source_id'] + \
                        [x.lower() for x in OBJECT_BASED_SOURCE_KEYS]

    by_object_metadata = ArchiveMetadata(os.path.join(archive_root, "by_object_meta.dat"),
                                            by_object_names,
                                            partial(get_object_based_metadata, 
                                                    OBJECT_BASED_HEADER_KEYS,
                                                    OBJECT_BASED_SPEC_KEYS,
                                                    OBJECT_BASED_SOURCE_KEYS),
                                            append=True)

    # metadatas in archive object
    return ArchiveDir(archive_root, [by_id_metadata, by_object_metadata],
                      copy_to_archive=copy_to_archive)



def entry_point():
    main(parse_args())

if __name__ == '__main__':
    entry_point()
