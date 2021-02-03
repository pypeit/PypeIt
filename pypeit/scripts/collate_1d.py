import sys
from datetime import datetime, timezone
import argparse
from glob import glob
import os.path
import shutil
from functools import partial

import numpy as np
from astropy.time import Time
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits, ascii
from astropy.table import QTable

from pypeit import specobjs
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import coadd1d
from pypeit.spec2dobj import AllSpec2DObj, Spec2DObj
from pypeit import msgs
from pypeit import par
from pypeit.slittrace import SlitTraceBitMask

#import yappi

# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

class ADAPArchive():
    def __init__(self, archive_root):
        self.archive_root = archive_root

    def write_metadata(self, file, metadata):
        (root, ext) = os.path.splitext(file)
        output_file = root + '.dat'
        with open(output_file, "w") as f:
            ascii.write(metadata, f, format='ipac')

    def get_file_based_metadata(self, header):
        metadata = ['FILENAME', 'RA', 'DEC', 'TARGET', 'DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
        data = [[header[x]] for x in metadata]
        names = [x.lower() for x in metadata]
        return QTable(data, names=names)

    def get_source_based_metadata(self, source):
        header_metadata = ['DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
        spec_obj_metadata = ['MASKDEF_OBJNAME', 'RA', 'DEC']
        header_data = [[source.spec1d_header_list[0][x]] for x in header_metadata]
        spec_obj_data = [[source.spec_obj_list[0][x]] for x in spec_obj_metadata]
        names = [x.lower() for x in spec_obj_metadata] + [x.lower() for x in header_metadata]
        return QTable(spec_obj_data + header_data, names=names)

    def add_spec1d_files(self, spec1d_files):

        for i in range(len(spec1d_files)):
            sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_files[i])

            spec1d_file = self.archive_file(spec1d_files[i], Time(sobjs.header['MJD'], format='mjd'))
            metadata = self.get_file_based_metadata(sobjs.header)
            self.write_metadata(spec1d_file, metadata)

    def add_spec2d_files(self, spec2d_files):

        for i in range(len(spec2d_files)):
            allspec2d = AllSpec2DObj.from_fits(spec2d_files[i])

            spec2d_file = self.archive_file(spec2d_files[i], Time(allspec2d['meta']['head0']['MJD'], format='mjd'))
            metadata = self.get_file_based_metadata(allspec2d['meta']['head0'])
            self.write_metadata(spec2d_file, metadata)

    def add_coadd_sources(self, sources):
        for source in sources:
            coadd_file = self.archive_file(source.coaddfile, Time(source.spec1d_header_list[0]['MJD'], format='mjd'))
            metadata = self.get_source_based_metadata(source)
            self.write_metadata(coadd_file, metadata)

    def get_archive_path(self, time):
        return os.path.join(self.archive_root, time.strftime('%Y/%m/%d'))

    def archive_file(self, orig_file, time):

        if not os.path.exists(orig_file):
            msgs.error(f'File {orig_file} does not exist')

        path = self.get_archive_path(time)
        os.makedirs(path, exist_ok=True)

        dest_file = os.path.join(path, os.path.basename(orig_file))

        msgs.info(f'Copying {orig_file} to archive root {self.archive_root}')
        try:
            shutil.copy2(orig_file, dest_file)
        except:
            msgs.error(f'Failed to copy {orig_file} to {dest_file}')

        return dest_file


class SourceObject:
    def __init__(self, spec1d_obj, spec1d_header, spec1d_file, spectrograph, match_type):
        self.spec1d_file_list = [spec1d_file]
        self.spec1d_header_list = [spec1d_header]
        self.spec_camera = spectrograph.camera
        self.spec_obj_list = [spec1d_obj]
        self.match_type = match_type

        if (match_type == 'ra/dec'):
            self.coord = SkyCoord(spec1d_obj.RA, spec1d_obj.DEC, unit='deg')
        else:
            self.coord = spec1d_obj['SPAT_PIXPOS']

        self.coaddfile = self.build_coadd_file_name()

    def build_coadd_file_name(self):
        """Build the output file name for coadding.
        The filename convention is J<hmsdms+dms>_DEIMOS_<YYYYMMDD>.fits.

        Args:
            source (`obj`:SourceObject): The Source object being coaadded. The firsts entry in spec1d_header_list is used
            to build the filename.

        Return: str  The name of the coadd output file.
        """
        time_portion = Time(self.spec1d_header_list[0]['MJD'], format="mjd").strftime('%Y%m%d')
        if self.match_type == 'ra/dec':
            coord_portion = 'J' + self.coord.to_string('hmsdms', sep='', precision=2).replace(' ', '')
        else:
            coord_portion = self.spec_obj_list[0]['NAME'].split('_')[0]

        return f'{coord_portion}_{self.spec_camera}_{time_portion}.fits'

    def match(self, spec_obj, thresh):
        if self.match_type == 'ra/dec':
            coord2 = SkyCoord(ra=spec_obj.RA, dec=spec_obj.DEC, unit='deg')
            return self.coord.separation(coord2) <= thresh
        else:
            coord2 =spec_obj['SPAT_PIXPOS'] 
            return np.fabs(coord2 - self.coord) <= thresh

def find_slits_to_exclude(spec2d_files, par):
    """Find slits that should be excluded according to the input parameters. The slit mask ids are returned in a map
    alongside the text labels for the flags that caused the slit to be excluded.

    Args:
        spec2d_files (:obj:`list` of str): List of spec2d files to build the map from.
        par (:obj:`Collate1DPar): Parameters from .collate1d file

    Return:
        :obj:`dict` Mapping of slit mask ids to the flags that caused the slit to be excluded.
    """

    # Get the types of slits to exclude from our parameters
    exclude_flags = par['collate1d']['slit_exclude_flags']
    if isinstance(exclude_flags, str):
        exclude_flags = [exclude_flags]

    bit_mask = SlitTraceBitMask()
    exclude_map = dict()
    for spec2d_file in spec2d_files:

        allspec2d = AllSpec2DObj.from_fits(spec2d_file)
        for sobj2d in [allspec2d[det] for det in allspec2d.detectors]:
            for (slit_id, mask, slit_mask_id) in sobj2d['slits'].slit_info:
                for flag in exclude_flags:
                    if bit_mask.flagged(mask, flag):
                        if slit_mask_id not in exclude_map:
                            exclude_map[slit_mask_id] = {flag}
                        else:
                            exclude_map[slit_mask_id].add(flag)

    return exclude_map

def config_key_match(spectrograph, header1, header2):
    for key in spectrograph.configuration_keys():
        # Ignore "decker" because it's valid to coadd spectra with different slit masks
        if key != "decker":
            if key not in header1 and key not in header2:
                # Both are missing the key, this is ok
                continue
            elif key not in header1 or key not in header2:
                # One has a value and the other doesn't so they don't match
                return False

            if header1[key] != header2[key]:
                return False

    # No mismatches were found
    return True

def group_spectra_by_source(spec1d_files, spectrograph, exclude_map, match_type, thresh):
    """Given a list of spec1d files from PypeIt, group the spectra within the files by their source object.
    The grouping is done by comparing the RA/DEC of each spectra using a given threshold.

    Args:
        spec1d_files (list of str): A list of spec1d files created by PypeIt
        exclude_map: (dict): Mapping of excluded slit mask ids to the reason they were excluded.
                                   Any spectra for these slits will be ignored.
        thresh: (`obj`:astropy.coordinates.Angle): Maximum angular distance that two spectra can be from each other
                                                   to be considered to be from the same source.

    Returns: list of `obj`:SourceObject The spectra and spec1d files grouped into SourceObjects

    """
    source_list = []

    for spec1d_file in spec1d_files:
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)

        for sobj in sobjs.specobjs:

            if sobj.MASKDEF_ID in exclude_map:
                msgs.info(f'Excluding {sobj.MASKDEF_ID} in {spec1d_file} because of flags {exclude_map[sobj.MASKDEF_ID]}')
                continue
            if sobj.OPT_COUNTS is None:
                msgs.info(f'Excluding {sobj.NAME} in {spec1d_file} because of missing OPT_COUNTS')
                continue

            found = False
            for source in source_list:
                if config_key_match(spectrograph, source.spec1d_header_list[0], sobjs.header) and \
                    source.match(sobj, thresh):
                    source.spec_obj_list.append(sobj)
                    source.spec1d_file_list.append(spec1d_file)
                    source.spec1d_header_list.append(sobjs.header)
                    found = True

            if not found:
                source_list.append(SourceObject(sobj, sobjs.header, spec1d_file, spectrograph, match_type))

    return source_list

def coadd(par, source):
    """coadd the spectra for a given source.

    Args:
        par (`obj`:Collate1DPar): Paramters for the coadding
        source (`obj`:SourceObject): The SourceObject with information on which files and spectra to coadd.
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
    Read a list of files from a PypeIt .collate1d file, akin to a standard PypeIt file.

    Reads from a block of file names in a configuration file such as:

    spec1d read
    Science/spec1d_DE.20130409.20629-S13A-SDF-z6clus_DEIMOS_2013Apr09T054342.730.fits
    spec1d end


    Args:
        lines (:obj:`numpy.ndarray`): List of lines to read.
        block (str): Name of the block to read files from (e.g. 'spec1d', 'spec2d')

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values, i.e. lines that did not contain the file list.

        files (list):
          Contains the lsit of files read.


    """

    files = []
    is_config = np.ones(len(lines), dtype=bool)

    s, e = par.util._find_pypeit_block(lines, block)
    if s >= 0 and e < 0:
        msgs.error(f"Missing '{block} end' in collate1d config file")
    elif (s < 0) or (s==e):
        msgs.error(f"Missing {block} block in collate1d config file")
    else:
        files = lines[s:e]

    is_config[s-1:e+1] = False

    return is_config, files

def read_collate_file(collate_file):
    """
    Read the file lists from a PypeIt .collate1d file, akin to a standard PypeIt file.

    The top is a config block that sets ParSet parameters

    After that are the list of spec1d and spec 2d files.

    Args:
        collate_file (str):
          Name of the .collate1d file

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values
        spec1dfiles (list):
          Contains spec1dfiles to be coadded
        spec2dfiles (list):
          Contains spec2dfiles that slit metadata is pulled from

    """

    # Read in the config file
    msgs.info(f'Loading the {collate_file} config file')
    lines = par.util._read_pypeit_file_lines(collate_file)

    # Parse the file list
    is_config, spec1d_files = read_file_list(lines, 'spec1d')

    # Filter out file lists from configuration
    cfg_lines = list(lines[is_config])

    return  cfg_lines, spec1d_files

def find_spec2d_from_spec1d(spec1d_files):
    """
    Find the spec2d files corresponding to the given list of spec1d files.
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
    Read the command line arguments and the input .collate1d file (if any), to build the parameters needed by collate_1d.
    """

    # First we need to get the list of spec1d files
    if args.input_file is not None:
        (cfg_lines, spec1d_patterns) = read_collate_file(args.input_file)
    else:
        cfg_lines = None
        spec1d_patterns = []

    if args.spec1d_files is not None and len(args.spec1d_files) > 0:
        spec1d_patterns = args.spec1d_files

    # Expand any wild cards in the spec1d files
    spec1d_files = []
    for pattern in spec1d_patterns:
        spec1d_files += glob(pattern)

    if spec1d_files is None or len(spec1d_files) == 0:
        msgs.error("A list of spec1d files must be specified via command line or config file.")

    # Expand any wildcards in the spec2d files

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
    if args.thresh is not None:
        params['collate1d']['threshold'] = args.thresh

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

    # A blank Colate1DPar to avoid duplicating the help text.
    blank_par = pypeitpar.Collate1DPar()

    parser = argparse.ArgumentParser(description='Flux/Coadd multiple 1d spectra from multiple nights and archive the results.',
                                     formatter_class=SmartFormatter)

    parser.add_argument('input_file', type=str,
                        help='File with stuff and things (TODO real documentation)',
                        nargs='?')
    parser.add_argument('--spec1d_files', type=str, nargs='*', help='One or more spec1d files to ' \
                        'flux/coadd/archive. Can contain wildcards')
    parser.add_argument('--par_outfile', default='collate1d.par', type=str,
                        help='Output to save the parameters')
    parser.add_argument('--thresh', type=str, help=blank_par.descr['threshold'])
    parser.add_argument('--match', type=str, choices=blank_par.options['match_using'], help=blank_par.descr['match_using'])
    parser.add_argument('--dry_run', action='store_true', help=blank_par.descr['dry_run'])
    parser.add_argument('--archive_dir', type=str, help=blank_par.descr['archive_root'])
    parser.add_argument('--exclude_slit', type=str, nargs='*', help=blank_par.descr['slit_exclude_flags'])

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

def main(args):

    #yappi.set_clock_type("wall")
    #yappi.start(builtins=False)
    (par, spectrograph, spec1d_files) = build_parameters(args)

    spec2d_files = find_spec2d_from_spec1d(spec1d_files)

    # Write the par to disk
    print("Writing the parameters to {}".format(args.par_outfile))
    par.to_config(args.par_outfile)

    exclude_map = find_slits_to_exclude(spec2d_files, par)

    if par['collate1d']['match_using'] == 'pixel':
        threshold = np.float(par['collate1d']['threshold'])
    else:
        threshold = Angle(par['collate1d']['threshold'])

    source_list = group_spectra_by_source(spec1d_files, spectrograph, exclude_map, par['collate1d']['match_using'], threshold)

    #sensfunc, how to identify standard file

    # fluxing etc goes here

    for source in source_list:

        msgs.info(f'Creating {source.coaddfile} from the following sources:')
        for i in range(len(source.spec_obj_list)):
            msgs.info(f'    {source.spec1d_file_list[i]}: {source.spec_obj_list[i].NAME} ({source.spec_obj_list[i].MASKDEF_OBJNAME})')

        if not args.dry_run:
            coadd(par, source)

    if not args.dry_run and par['collate1d']['archive_root'] is not None:
        archive = ADAPArchive(par['collate1d']['archive_root'])
        archive.add_spec1d_files(spec1d_files)
        archive.add_spec2d_files(spec2d_files)
        archive.add_coadd_sources(source_list)

    #yappi.stop()
    #all_thread_func_stats = yappi.get_func_stats()

    #with open("yappi_results.txt", "w") as text_file:
    #    all_thread_func_stats.print_all(out=text_file, columns={0: ("name", 60),
    #                                                            1: ("ncall", 15),
    #                                                            2: ("tsub", 8),
    #                                                            3: ("ttot", 8),
    #                                                            4: ("tavg", 8)})
    #all_thread_func_stats.save("yappi_results.pstat", type="pstat")
    #all_thread_func_stats.save("yappi_results.callgrind", type="callgrind")
    return 0