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

import sys
from datetime import datetime, timezone
import argparse
from glob import glob
import os.path
import shutil
from functools import partial
import re

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle
from astropy.io import fits, ascii
from astropy.table import Table

from pypeit import specobjs
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit import coadd1d
from pypeit.spec2dobj import AllSpec2DObj, Spec2DObj
from pypeit import msgs
from pypeit import par
from pypeit.slittrace import SlitTraceBitMask
from pypeit.utils import is_float

# A trick from stackoverflow to allow multi-line output in the help:
#https://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text
class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)

class ArchiveDir():
    """
    Copies files and metadata to a directory for archival purposes.

    Files are all copied to the top level directory in the archive. 

    If a file originates from KOA the KOAID will be extracted either from the
    ``KOAID`` header keyword or from the ``FILENAME`` header keyword.

    A KOAID has the format: ``II.YYYYMMDD.xxxxx``
        See the `KOA FAQ <https://www2.keck.hawaii.edu/koa/public/faq/koa_faq.php>`_ 
        for more information.

    Metadata is written to two files in the 
    `ipac <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_ format. 

    ``by_id_meta.dat`` contains metadata for the spec1d and spec2d files in
    the archive. It is organzied by the id (either KOAID, or file name) of the 
    original science image.

    ``by_object_meta.dat`` contains metadata for the coadded output files.
    This may have multiple rows for each file depending on how many science
    images were coadded. The primary key is a combined key of the source 
    object name, filename, and koaid columns.

    Args:
        archive_root (str): Root directory where the files and metadata will be
        placed. This will be created if needed.

    Attributes:
    by_id_metadata (:obj:`list of :obj:`list of str``): List of metadata rows
        for the spec1d and spec2d files in the archive diredctory. This will be
        loaded from any existing by_id_meta.dat file in the archive_root on 
        object initialization.

    by_object_metadata (:obj:`list of :obj:`list of str``): List of metadata rows
        for the coadd output files in the archive directory. This will be
        loaded from any existing by_id_meta.dat file in the archive_root on 
        object initialization.

    """

    # Header and SpecObj keys for metadata needed for the IPAC files
    _ID_BASED_HEADER_KEYS  = ['RA', 'DEC', 'TARGET', 'PJROGPI', 'SEMESTER', 'PROGID', 'DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
    _OBJECT_BASED_HEADER_KEYS = ['DISPNAME', 'DECKER', 'BINNING', 'MJD', 'AIRMASS', 'EXPTIME']
    _OBJECT_BASED_SPEC_KEYS   = ['MASKDEF_OBJNAME', 'MASKDEF_ID', 'DET', 'RA', 'DEC']
    _OBJECT_BASED_SOURCE_KEYS = ['GUIDFWHM', 'PJROGPI', 'SEMESTER', 'PROGID']

    def __init__(self, archive_root, copy=True):
        self.archive_root = archive_root
        
        # Load metadata from any pre-existing metadata files.
        # Because astropy Tables are slow at adding rows, we convert 
        # the metadata to a list of lists for performance as adding rows is
        # the primary feature of this class.
        self._by_id_file = os.path.join(archive_root, 'by_id_meta.dat')
        if os.path.exists(self._by_id_file):
            by_id_table = ascii.read(self._by_id_file)
            self.by_id_metadata = [list(row) for row in by_id_table]
        else:
            self.by_id_metadata = []

        self._by_object_file = os.path.join(archive_root, 'by_object_meta.dat')
        if os.path.exists(self._by_object_file):
            by_object_table = ascii.read(self._by_object_file)
            self.by_object_metadata = [list(row) for row in by_object_table]
        else:
            self.by_object_metadata = []

        self._copy_files = copy

    def save(self):
        """
        Saves the metadata in this class to IPAC files in the archive directory
        """
        by_id_names = ['id', 'filename'] + [x.lower() for x in self._ID_BASED_HEADER_KEYS]
        with open(self._by_id_file, 'w') as f:
            ascii.write(Table(rows=self.by_id_metadata, names=by_id_names), f, format='ipac')

        by_object_names = ['filename'] + \
                          [x.lower() for x in self._OBJECT_BASED_SPEC_KEYS] + \
                          [x.lower() for x in self._OBJECT_BASED_HEADER_KEYS] + \
                          ['source_id'] + \
                          [x.lower() for x in self._OBJECT_BASED_SOURCE_KEYS]

        with open(self._by_object_file, 'w') as f:
            ascii.write(Table(rows=self.by_object_metadata, names=by_object_names), f, format='ipac')

    def _extract_id(self, header):
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

    def _get_id_based_metadata(self, filename, header):
        """
        Gets the metadata from a FITS header used for the by id portion
        of the archive and appends it to by_id_metadata

        Args:
            filename (str): A filename from a file originating in the KOA.
        
        """
        # Extract koa id from source image filename in header
        id =  self._extract_id(header)

        # Build data row, which starts with koaid and filename + the metadata
        data_row = [id, filename] + [None if x not in header else header[x] for x in self._ID_BASED_HEADER_KEYS]

        self.by_id_metadata.append(data_row)

    def _get_object_based_metadata(self, source):
        """
        Gets the metadata from a SourceObject instance used for the by object
        portion of the archive and appends it to by_object_metadata.

        Args:
        source (:obj:`pypeit.scripts.collate_1d.SourceObject`)): The source object containing the
            headers, filenames and SpecObj information for a coadd output file.
        """

        header_data = [source.spec1d_header_list[0][x] for x in self._OBJECT_BASED_HEADER_KEYS]
        spec_obj_data = [source.spec_obj_list[0][x] for x in self._OBJECT_BASED_SPEC_KEYS]
        shared_data = [os.path.basename(source.coaddfile)] + spec_obj_data + header_data
        
        for header in source.spec1d_header_list:
            id = self._extract_id(header)
            unique_data = [header[x] if x in header else None for x in self._OBJECT_BASED_SOURCE_KEYS]
            self.by_object_metadata.append(shared_data + [id] + unique_data)

    def add_files(self, files):
        """
        Copy fits files to the archive directory and add their metadata to the
        by_id_metadata member variable. 

        Args:
        files (list of str): List of full pathnames for tbe fits files
        """

        for file in files:
            header = fits.getheader(file)

            self._archive_file(file)
            self._get_id_based_metadata(os.path.basename(file), header)

    def add_coadd_sources(self, sources):
        """
        Copy coadd output files to the archive directory and add their metadata to the
        by_object_metadata member variable. 

        Args:
        sources (list of :obj:`pypeit.scripts.collate_1d.SourceObject`): 
            List of full pathnames for tbe spec2d files
        """
        for source in sources:
            self._archive_file(source.coaddfile)
            self._get_object_based_metadata(source)
            

    def _archive_file(self, orig_file):
        """
        Copies a file to the archive directory, if copying
        is enable.

        Args:
        orig_file (str): Path to the file to copy.

        Returns:
        str: The full path to the new copy in the archive.
        """

        if self._copy_files is False:
            return orig_file

        if not os.path.exists(orig_file):
            msgs.error(f'File {orig_file} does not exist')

        os.makedirs(self.archive_root, exist_ok=True)
        dest_file = os.path.join(self.archive_root, os.path.basename(orig_file))

        msgs.info(f'Copying {orig_file} to archive root {self.archive_root}')
        try:
            shutil.copy2(orig_file, dest_file)
        except:
            msgs.error(f'Failed to copy {orig_file} to {dest_file}')

        return dest_file


class SourceObject:
    """ A group of reduced spectra from the same source object. This contains
    the information needed to coadd the spectra and archive the metadata.

    An instance is initiated with the first spectra of the group. Additional
    spectra can be compared with this object to see if it matches using the
    match method, and are added to it if they do.

    Args:
    spec1d_obj (:obj:`pypeit.specobj.SpecObj`):
        The initial spectra of the group as a SpecObj.
    spec1d_header (:obj:`astropy.io.fits.Header`): 
        The header for the first spec1d file in the group.
    spec1d_file (str): Filename of the first spec1d file in the group.
    spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`): 
        The spectrograph that was used to take the data.
    match_type (str): How spectra should be compared. 'ra/dec' means the
        spectra should be compared using the sky coordinates in RA and DEC.
        'pixel' means the spectra should be compared by the spatial pixel
        coordinates in the image.

    Attributes:
    spec_obj_list (list of :obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
        The list of spectra in the group as SpecObj objects.
    spec1d_file_list (list of str): 
        The pathnames of the spec1d files in the group.
    spec1d_header_list: (list of :obj:`astropy.io.fits.Header`):
        The headers of the spec1d files in the group
    """
    def __init__(self, spec1d_obj, spec1d_header, spec1d_file, spectrograph, match_type):
        self.spec_obj_list = [spec1d_obj]
        self.spec1d_file_list = [spec1d_file]
        self.spec1d_header_list = [spec1d_header]
        self._spectrograph = spectrograph
        self.match_type = match_type

        if (match_type == 'ra/dec'):
            self.coord = SkyCoord(spec1d_obj.RA, spec1d_obj.DEC, unit='deg')
        else:
            self.coord = spec1d_obj['SPAT_PIXPOS']

        self.coaddfile = self.build_coadd_file_name()

    def build_coadd_file_name(self):
        """Build the output file name for coadding.
        The filename convention is J<hmsdms+dms>_<instrument name>_<YYYYMMDD>.fits
        when matching by RA/DEC and SPAT_<spatial position>_<instrument name>_<YYYYMMDD>.fits
        when matching by pixel position..

        Currently instrument_name is taken from spectrograph.camera

        Return: str  The name of the coadd output file.
        """
        time_portion = Time(self.spec1d_header_list[0]['MJD'], format="mjd").strftime('%Y%m%d')
        if self.match_type == 'ra/dec':
            coord_portion = 'J' + self.coord.to_string('hmsdms', sep='', precision=2).replace(' ', '')
        else:
            coord_portion = self.spec_obj_list[0]['NAME'].split('_')[0]
        instrument_name = self._spectrograph.camera

        return f'{coord_portion}_{instrument_name}_{time_portion}.fits'

    def _config_key_match(self, header):
        """
        Check to see if the configuration keys from a spec1d file match the
        ones for this SourceObject.

        Args:
        header (:obj:`astropy.io.fits.Header`):
            Header from a spec1d file.

        Returns (bool): True if the configuration keys match, 
            false if they do not.
        """
        # Make sure the spectrograph matches
        if 'PYP_SPEC' not in header or header['PYP_SPEC'] != self._spectrograph.name:
            return False

        first_header = self.spec1d_header_list[0]

        for key in self._spectrograph.configuration_keys():
            # Ignore "decker" because it's valid to coadd spectra with different slit masks
            if key != "decker":
                if key not in first_header and key not in header:
                    # Both are missing the key, this is ok
                    continue
                elif key not in first_header or key not in header:
                    # One has a value and the other doesn't so they don't match
                    return False

                if first_header[key] != header[key]:
                    return False

        # No mismatches were found
        return True

    def match(self, spec_obj, spec1d_header, tolerance, unit = u.arcsec):
        """Determine if a SpecObj matches this group within the given tolerance.
        This will also compare the configuration keys to make sure the SpecObj
        is compatible with the ones in this SourceObject.

        Args:
        spec_obj (:obj:`pypeit.specobj.SpecObj`): 
            The SpecObj to compare with this SourceObject.

        spec1d_header (:obj:`astropy.io.fits.Header`):
            The header from the spec1d that dontains the SpecObj.
            
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit1`` argument).

        unit (:obj:`astropy.units.Unit`):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Igored if match_type is 'pixel'.

        Returns (bool): True if the SpecObj matches this group,
            False otherwise.
        """

        if not self._config_key_match(spec1d_header):
            return False

        if self.match_type == 'ra/dec':
            coord2 = SkyCoord(ra=spec_obj.RA, dec=spec_obj.DEC, unit='deg')
            return self.coord.separation(coord2) <= Angle(tolerance, unit=unit)
        else:
            coord2 =spec_obj['SPAT_PIXPOS'] 
            return np.fabs(coord2 - self.coord) <= tolerance

def find_slits_to_exclude(spec2d_files, par):
    """Find slits that should be excluded according to the input parameters.
    The slit mask ids are returned in a map alongside the text labels for the
    flags that caused the slit to be excluded.

    Args:
        spec2d_files (:obj:`list` of str): 
            List of spec2d files to build the map from.
        par (:obj:`Collate1DPar):
            Parameters from .collate1d file

    Return:
        :obj:`dict` Mapping of slit mask ids to the flags that caused the
        slit to be excluded.
    """

    # Get the types of slits to exclude from our parameters
    exclude_flags = par['collate1d']['slit_exclude_flags']
    if isinstance(exclude_flags, str):
        exclude_flags = [exclude_flags]

    # Go through the slit_info of all spec2d files and find
    # which slits should be excluded based on their flags
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

def group_spectra_by_source(spec1d_files, exclude_map, match_type, tolerance, unit=u.arcsec):
    """Given a list of spec1d files from PypeIt, group the spectra within the
    files by their source object. The grouping is done by comparing the 
    position of each spectra (using either pixel or RA/DEC) using a given tolerance.

    Args:
        spec1d_files (list of str): A list of spec1d files created by PypeIt.
        exclude_map (dict): Mapping of excluded slit mask ids to the reason
            they were excluded.  Any spectra for these slits will be ignored.
        match_type (str): How to match the positions of spectra. 'pixel' for
            matching by spatial pixel distance or 'ra/dec' for matching by
            angular distance.
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit1`` argument).
        unit (:obj:`astropy.units.Unit`):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Igored if match_type is 'pixel'.

    Returns (list of `obj`:SourceObject): The grouped spectra as SourceObjects.

    """
    source_list = []

    for spec1d_file in spec1d_files:
        sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)

        for sobj in sobjs.specobjs:

            if sobj.MASKDEF_ID in exclude_map:
                msgs.info(f'Excluding {sobj.MASKDEF_ID} in {spec1d_file} because of flags {exclude_map[sobj.MASKDEF_ID]}')
                continue
            if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
                msgs.info(f'Excluding {sobj.NAME} in {spec1d_file} because of missing both OPT_COUNTS and BOX_COUNTS')
                continue

            # Search for a SourceObject that matches this SpecObj.
            # If one can't be found, trat this as a new SourceObject.
            found = False
            for source in source_list:
                if  source.match(sobj, sobjs.header, tolerance, unit):
                    source.spec_obj_list.append(sobj)
                    source.spec1d_file_list.append(spec1d_file)
                    source.spec1d_header_list.append(sobjs.header)
                    found = True

            if not found:
                spectrograph = load_spectrograph(sobjs.header['PYP_SPEC'])
                source_list.append(SourceObject(sobj, sobjs.header, spec1d_file, spectrograph, match_type))

    return source_list

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
    Read the file lists from a PypeIt .collate1d file, akin to a standard
    PypeIt file.

    The top is a config block that sets ParSet parameters

    After that are the list of spec1d files.

    Args:
        collate_file (str):
          Name of the .collate1d file

    Returns:
        cfg_lines (list):
          Config lines to modify ParSet values
        spec1dfiles (list):
          Contains spec1dfiles to be coadded
    """

    # Read in the config file
    msgs.info(f'Loading the {collate_file} config file')
    lines = par.util._read_pypeit_file_lines(collate_file)

    # Parse the file list
    is_config, spec1d_files = read_file_list(lines, 'spec1d')

    # Filter out file lists from configuration
    cfg_lines = list(lines[is_config])

    return  cfg_lines, spec1d_files

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
    is_config, coadd_lines = read_file_list(lines, 'coadd1d')

    cfg_lines = list(lines[is_config])

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
        (cfg_lines, spec1d_patterns) = read_collate_file(args.input_file)

        # Look for a coadd1d file
        (input_file_root, input_file_ext) = os.path.splitext(args.input_file)
        coadd1d_config_name = input_file_root + ".coadd1d"
        if os.path.exists(coadd1d_config_name):
            cfg_lines += read_coadd1d_config(coadd1d_config_name)

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
            archive = ArchiveDir(par['collate1d']['archive_root'])
        else:
            archive = ArchiveDir(os.getcwd(), copy=False)

        archive.add_files(spec1d_files)
        archive.add_files(spec2d_files)
        archive.add_coadd_sources(source_list)
        archive.save()

    total_time = datetime.now() - start_time

    msgs.info(f'Total duration: {total_time}')

    return 0

def entry_point():
    main(parse_args())

if __name__ == '__main__':
    entry_point()
