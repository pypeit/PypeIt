#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This module contains code for collating multiple 1d spectra by source object.

.. include:: ../include/links.rst
"""


import copy
import os.path
import shutil
import traceback
from functools import partial

import numpy as np
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import SkyCoord, Angle

from pypeit.spectrographs.util import load_spectrograph
from pypeit import log
from pypeit import outputPaths
from pypeit import outputfiles
from pypeit import PypeItError
from pypeit import coadd1d
from pypeit.core import wave
from pypeit.utils import radec_to_coord
from pypeit.slittrace import SlitTraceBitMask
from pypeit.spec2dobj import AllSpec2DObj
from pypeit.sensfilearchive import SensFileArchive
from pypeit import fluxcalibrate
from pypeit.specobjs import SpecObjs
from pypeit.archive import ArchiveMetadata, ArchiveDir


class SourceObject:

    """A group of reduced spectra from the same source object. This contains
    the information needed to coadd the spectra and archive the metadata.

    An instance is initiated with the first spectra of the group. Additional
    spectra can be compared with this object to see if it matches using the
    match method, and are added to it if they do.

    Args:
        spec1d_obj (:obj:`pypeit.specobj.SpecObj`):
            The initial spectra of the group as a SpecObj.
        spec1d_header (`astropy.io.fits.Header`_): 
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
        spec1d_header_list: (list of `astropy.io.fits.Header`_):
            The headers of the spec1d files in the group
    """
    def __init__(self, spec1d_obj, spec1d_header, spec1d_file, spectrograph, match_type):
        self.spec_obj_list = [spec1d_obj]
        self.spec1d_file_list = [spec1d_file]
        self.spec1d_header_list = [spec1d_header]
        self._spectrograph = spectrograph
        self.match_type = match_type

        if (match_type == 'ra/dec'):
            try:
                self.coord = SkyCoord(spec1d_obj.RA, spec1d_obj.DEC, unit='deg')
            except Exception as e:
                raise PypeItError(f"Cannot do ra/dec matching on {spec1d_obj.NAME}, could not read RA/DEC.")
        else:
            self.coord = spec1d_obj['SPAT_PIXPOS']

    @classmethod
    def build_source_objects(cls, specobjs_list, spec1d_files, match_type):
        """Build a list of SourceObjects from a list of spec1d files. There will be one SourceObject per
        SpecObj in the resulting list (i.e. no combining or collating is done by this method).

        Args:
            specobjs_list (list of :obj:`pypeit.specobjs.SpecObjs`): List of SpecObjs objects to build from.

            spec1d_files (list of str): List of spec1d filenames corresponding to each SpecObjs object.

            match_type (str):           What type of matching the SourceObjects will be configured for.
                                        Must be either 'ra/dec' or 'pixel'
        Returns: 
            list of :obj:`SourceObject`: A list of uncollated SourceObjects with one SpecObj per SourceObject.
        """
        result = []
        for i, sobjs in enumerate(specobjs_list):
            spectrograph = load_spectrograph(sobjs.header['PYP_SPEC'])
            for sobj in sobjs:
                result.append(SourceObject(sobj, sobjs.header, spec1d_files[i], spectrograph, match_type))
    
        return result

    def _config_key_match(self, header):
        """
        Check to see if the configuration keys from a spec1d file match the
        ones for this SourceObject.

        Args:
            header (`astropy.io.fits.Header`_):
                Header from a spec1d file.

        Returns:
            bool: True if the configuration keys match, 
                  false if they do not.
        """
        # Make sure the spectrograph matches
        if 'PYP_SPEC' not in header or header['PYP_SPEC'] != self._spectrograph.name:
            return False

        # Build the config to compare against from the first header, 
        # ignoring "decker" because it's valid to coadd spectra with different slit masks
        first_config = {key: self.spec1d_header_list[0][key] 
                            for key in self._spectrograph.configuration_keys() if key != 'decker'}

        second_config = {key: header[key] 
                            for key in self._spectrograph.configuration_keys() if key != 'decker'}

        # Use spectrograph.same_configuration to compare the configurations, with check_keys set to False
        # so only the keys in first_config are checked
        return self._spectrograph.same_configuration([first_config, second_config],check_keys=False)

    def match(self, spec_obj, spec1d_header, tolerance, unit = u.arcsec):
        """Determine if a SpecObj matches this group within the given tolerance.
        This will also compare the configuration keys to make sure the SpecObj
        is compatible with the ones in this SourceObject.

        Args:
            spec_obj (:obj:`pypeit.specobj.SpecObj`): 
                The SpecObj to compare with this SourceObject.
            spec1d_header (`astropy.io.fits.Header`_):
                The header from the spec1d that dontains the SpecObj.
            tolerance (float): 
                Maximum distance that two spectra can be from each other to be 
                considered to be from the same source. Measured in floating
                point pixels or as an angular distance (see ``unit1`` argument).
            unit (`astropy.units.Unit`_):
                Units of ``tolerance`` argument if match_type is 'ra/dec'. 
                Defaults to arcseconds. Igored if match_type is 'pixel'.

        Returns:
            bool: True if the SpecObj matches this group, False otherwise.
        """

        if not self._config_key_match(spec1d_header):
            return False

        if self.match_type == 'ra/dec':
            coord2 = SkyCoord(ra=spec_obj.RA, dec=spec_obj.DEC, unit='deg')
            return self.coord.separation(coord2) <= Angle(tolerance, unit=unit)
        else:
            coord2 =spec_obj['SPAT_PIXPOS'] 
            return np.fabs(coord2 - self.coord) <= tolerance

    def combine(self, other_source_object):
        """Combine this SourceObject with another. The two objects must be from the
        same spectrograph and use the same match type.

        Args:
            other_source_object (:obj:`SourceObject`): The other object to combine with.

        Returns:
            (:obj:`SourceObject`): This SourceObject, now combined with other_source_object.
        """
        
        if other_source_object._spectrograph.name != self._spectrograph.name or \
           other_source_object.match_type != self.match_type:
           raise PypeItError(f"Can't append incompatible source objects. {self._spectrograph.name}/{self.match_type} does not match {other_source_object._spectrograph.name}/{other_source_object.match_type}")

        self.spec_obj_list += other_source_object.spec_obj_list
        self.spec1d_file_list += other_source_object.spec1d_file_list
        self.spec1d_header_list += other_source_object.spec1d_header_list

        return self

def create_report_archive():
    """
    Create an report archive with the desired metadata information.

    Metadata is written to three files in the `ipac
    <https://irsa.ipac.caltech.edu/applications/DDGEN/Doc/ipac_tbl.html>`_
    format:

        - ``collate_report.dat`` contains metadata to report on the coadded output files
          from the collate process. Like ``coadded_files.dat`` it may have more
          than one row per output file.  This file is always written to the current directory.     

    Returns:
        :class:`~pypeit.archive.ArchiveDir`: Object for archiving files and/or
        metadata.
    """
    archive_metadata_list = []

    COADDED_SPEC1D_HEADER_KEYS  = ['DISPNAME', 'DECKER',   'BINNING', 'MJD', 'AIRMASS', 'EXPTIME','GUIDFWHM', 'PROGPI', 'SEMESTER', 'PROGID']
    COADDED_SPEC1D_COLUMN_NAMES = ['dispname', 'slmsknam', 'binning', 'mjd', 'airmass', 'exptime','guidfwhm', 'progpi', 'semester', 'progid']

    COADDED_SOBJ_KEYS  =        ['MASKDEF_OBJNAME', 'MASKDEF_ID', 'NAME',        'DET', 'RA',    'DEC',    'MASKDEF_OBJMAG', 'MASKDEF_OBJMAG_BAND', 'S2N', 'MASKDEF_EXTRACT', 'WAVE_RMS']
    COADDED_SOBJ_COLUMN_NAMES = ['maskdef_objname', 'maskdef_id', 'pypeit_name', 'det', 'objra', 'objdec', 'maskdef_objmag', 'maskdef_objmag_band', 's2n', 'maskdef_extract', 'wave_rms']

    report_names = ['filename'] + \
                   COADDED_SOBJ_COLUMN_NAMES + \
                   ['spec1d_filename'] + \
                   COADDED_SPEC1D_COLUMN_NAMES

    report_formats = {'s2n':      '%.2f',
                      'wave_rms': '%.3f'}

    report_metadata = ArchiveMetadata(str(outputPaths.collate / "collate_report.dat"),
                                      report_names,
                                      partial(get_report_metadata,
                                              COADDED_SPEC1D_HEADER_KEYS,
                                              COADDED_SOBJ_KEYS),
                                      append=True,
                                      formats= report_formats)
    archive_metadata_list.append(report_metadata)

    # metadatas in archive object
    return ArchiveDir(str(outputPaths.collate), archive_metadata_list, copy_to_archive=False)

def get_report_metadata(object_header_keys, spec_obj_keys, file_info):
    """
    Gets the metadata from a :class:`~pypeit.core.collate.SourceObject` instance
    used building a report on the results of collation.  It is intended to be
    wrapped in by functools partial object that passes in object_header_keys and
    spec_obj_keys. file_info is then passed as in by the
    :class:`~pypeit.archive.ArchiveMetadata` object.  Unlike the other
    get_*_metadata functions, this is not used for archiving; it is used for
    reporting on the results of collating.

    If another type of file is added to the ArchiveMetadata object, the
    file_info argument will not be a :class:`~pypeit.core.collate.SourceObject`,
    In this case, a list of ``None`` values are returned.

    Args:
        object_header_keys (list of str):
            The keys to read fom the spec1d headers from the SourceObject.

        spec_obj_keys (list of str):
            The keys to read from the (:class:`~pypeit.specobj.SpecObj`) objects in
            the SourceObject.

        file_info (:class:`~pypeit.core.collate.SourceObject`): 
            The source object containing the headers, filenames and SpecObj
            information for a coadd output file.

    Returns:
        tuple: A tuple of two lists:.

            - **data_rows** (:obj:`list` of :obj:`list`): The metadata rows
              built from the source object.

            - **files_to_copy** (iterable): An list of tuples of files to copy.
              Because this function is not used for archving data, this is
              always ``None``.

    """

    if not isinstance(file_info, SourceObject):
        return (None, None)

    coaddfile = build_coadd_file_name(file_info)
    result_rows = []
    for i in range(len(file_info.spec1d_header_list)):

        # Get the spec_obj metadata needed for the report
        spec_obj = file_info.spec_obj_list[i]
        header = file_info.spec1d_header_list[i]

        # Get the spec1d header metadata needed for the report
        # Use getattr for the spec_obj data
        spec_obj_data = [getattr(spec_obj, x) for x in spec_obj_keys]
        spec1d_filename =  os.path.basename(file_info.spec1d_file_list[i])
        header_data = [header[x] if x in header else None for x in object_header_keys]
        result_rows.append([coaddfile] + spec_obj_data + [spec1d_filename] + header_data)

    return (result_rows, None)


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

def exclude_source_objects(source_objects, exclude_map, par):
    """
    Exclude :class:`~pypeit.core.collate.SourceObject` objects based on a slit
    exclude map and the user's parameters.

    Args:
        source_objects (:obj:`list`): 
            List of uncollated :class:`~pypeit.core.collate.SourceObject`
            objects to filter. There should only be one
            :class:`~pypeit.specobj.SpecObj` per
            :class:`~pypeit.core.collate.SourceObject`.
        exclude_map (:obj:`dict`): 
            Mapping of excluded slit ids to the reasons they should be excluded.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`): 
            Configuration parameters from the command line or a configuration
            file.

    Returns:
        tuple: Tuple containing two lists:

            - **filtered_objects** (:obj:`list`): A list of
              :class:`~pypeit.core.collate.SourceObject` with any excluded ones
              removed.

            - **excluded_messages** (:obj:`list`): A list of messages explaining why some source objects were excluded.
    """
    filtered_objects = []
    excluded_messages= []
    for source_object in source_objects:

        sobj = source_object.spec_obj_list[0]
        spec1d_file = source_object.spec1d_file_list[0]

        if par['collate1d']['exclude_serendip'] and sobj.MASKDEF_OBJNAME == 'SERENDIP':
            msg = f'Excluding SERENDIP object from {sobj.NAME} in {spec1d_file}'
            log.info(msg)
            excluded_messages.append(msg)
            continue

        if par['collate1d']['wv_rms_thresh'] is not None and sobj.WAVE_RMS > par['collate1d']['wv_rms_thresh']:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} due to wave_rms {sobj.WAVE_RMS} > threshold {par["collate1d"]["wv_rms_thresh"]}'
            log.info(msg)
            excluded_messages.append(msg)
            continue

        if sobj.MASKDEF_ID in exclude_map:
            msg = f'Excluding {sobj.NAME} with mask id: {sobj.MASKDEF_ID} in {spec1d_file} because of flags {exclude_map[sobj.MASKDEF_ID]}'
            log.info(msg)
            excluded_messages.append(msg)
            continue

        if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing both OPT_COUNTS and BOX_COUNTS'
            log.warning(msg)
            excluded_messages.append(msg)
            continue

        if par['coadd1d']['ex_value'] == 'OPT':
            msg = None
            if sobj.OPT_COUNTS is None:
                msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing OPT_COUNTS. Consider changing ex_value to "BOX".'
            elif sobj.OPT_MASK is None:
                msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing OPT_MASK. Consider changing ex_value to "BOX".'
            else:
                if len(sobj.OPT_COUNTS[sobj.OPT_MASK]) == 0:
                    msg = f'Excluding {sobj.NAME} in {spec1d_file} because all of OPT_COUNTS was masked out. Consider changing ex_value to "BOX".'
            
            if msg is not None:
                log.warning(msg)
                excluded_messages.append(msg)
                continue

        if par['coadd1d']['ex_value'] == 'BOX':
            msg = None
            if sobj.BOX_COUNTS is None:
                msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing BOX_COUNTS. Consider changing ex_value to "OPT".'
            elif sobj.BOX_MASK is None:
                msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing BOX_MASK. Consider changing ex_value to "OPT".'
            else:
                if len(sobj.BOX_COUNTS[sobj.BOX_MASK]) == 0:
                    msg = f'Excluding {sobj.NAME} in {spec1d_file} because all of BOX_COUNTS was masked out. Consider changing ex_value to "OPT".'

            if msg is not None:
                log.warning(msg)
                excluded_messages.append(msg)
                continue

        filtered_objects.append(source_object)
    return (filtered_objects, excluded_messages)

def read_spec1d_files(par, spec1d_files, failure_log):
    """
    Read spec1d files.

    Args:
        par (`obj`:pypeit.par.pypeitpar.PypeItPar): 
            Parameters for collating, fluxing, and coadding.
        spec1d_files (list of str):
            List of spec1d files to read.
        failure_log(list of str):
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
            log.warning(formatted_exception)
            log.warning(f"Failed to read {spec1d_file}, skipping it.")
            failure_log.append(f"Failed to read {spec1d_file}, skipping it.")
            failure_log.append(formatted_exception)

    return specobjs_list, good_spec1d_files

def flux(par, spectrograph, spec1d_files, failed_fluxing_log):
    """
    Flux calibrate spec1d files using archived sens func files.

    Args:
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`): 
            Parameters for collating, fluxing, and coadding.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph for the files to flux.
        spec1d_files (list of str):
            List of spec1d files to flux calibrate.
        failed_fluxing_log(list of str):
            Return parameter describing any failures that occurred when fluxing.

    Returns:
        list of str: The spec1d files that were successfully flux calibrated.
    """

    # Make sure fluxing from archive is supported for this spectrograph
    if spectrograph.name not in SensFileArchive.supported_spectrographs():
        raise PypeItError(f"Flux calibrating {spectrograph.name} with an archived sensfunc is not supported.")

    par['fluxcalib']['extrap_sens'] = True

    sf_archive = SensFileArchive.get_instance(spectrograph.name)
    flux_calibrated_files = []
    for spec1d_file in spec1d_files:

        # Get the archived sens file to use
        try:
            sens_file = sf_archive.get_archived_sensfile(spec1d_file)
        except Exception:
            formatted_exception = traceback.format_exc()
            log.warning(formatted_exception)
            log.warning(f"Could not find archived sensfunc to flux {spec1d_file}, skipping it.")
            failed_fluxing_log.append(f"Could not find archived sensfunc to flux {spec1d_file}, skipping it.")
            failed_fluxing_log.append(formatted_exception)
            continue
            
        # Flux calibrate the spec1d file
        try:
            log.info(f"Running flux calibrate on {spec1d_file}")
            FxCalib = fluxcalibrate.flux_calibrate([spec1d_file], [sens_file], par=par['fluxcalib'],
                                                   chk_version=par['rdx']['chk_version'])
            flux_calibrated_files.append(spec1d_file)

        except Exception:
            formatted_exception = traceback.format_exc()
            log.warning(formatted_exception)
            log.warning(f"Failed to flux calibrate {spec1d_file}, skipping it.")
            failed_fluxing_log.append(f"Failed to flux calibrate {spec1d_file}, skipping it.")
            failed_fluxing_log.append(formatted_exception)
            continue

    # Return the succesfully fluxed files
    return flux_calibrated_files

def build_coadd_file_name(source_object):
    """Build the output file name for coadding.
    The filename convention is J<hmsdms+dms>_<instrument name>_<YYYYMMDD>.fits
    when matching by RA/DEC and SPAT_<spatial position>_<instrument name>_<YYYYMMDD>.fits
    when matching by pixel position. The date portion may be <YYYYMMDD-YYYMMDD> if the
    files coadded span more than one date.

    Currently instrument_name is taken from spectrograph.camera

    Returns: 
        str:  The name of the coadd output file.
    """
    mjd_list = [float(h['MJD']) for h in source_object.spec1d_header_list]
    start_mjd = np.min(mjd_list)
    end_mjd = np.max(mjd_list)

    start_date_portion = Time(start_mjd, format="mjd").strftime('%Y%m%d')
    end_date_portion = Time(end_mjd, format="mjd").strftime('%Y%m%d')

    date_portion = f"{start_date_portion}_{end_date_portion}"

    if source_object.match_type == 'ra/dec':
        coord_portion = 'J' + source_object.coord.to_string('hmsdms', sep='', precision=2).replace(' ', '')
    else:
        coord_portion = source_object.spec_obj_list[0]['NAME'].split('_')[0]
    instrument_name = source_object._spectrograph.camera
    
    return f'{coord_portion}_{instrument_name}_{date_portion}.fits'

def refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_log):

    refframe = par['collate1d']['refframe']
    log.info(f"Performing a {refframe} correction")

    for spec1d in spec1d_files:
        # Get values from the fits header needed to calculate the correction
        try:
            sobjs = SpecObjs.from_fitsfile(spec1d, chk_version=par['rdx']['chk_version'])
            hdr_ra = sobjs.header['RA']
            hdr_dec = sobjs.header['DEC']
            hdr_radec = radec_to_coord((hdr_ra, hdr_dec))
            obstime = Time(sobjs.header['MJD'], format='mjd')
        except Exception as e:
            msg = f'Failed to perform {refframe} correction on {spec1d}: {e}'
            log.info(msg)
            spec1d_failure_log.append(msg)
            continue

        corrected_at_least_one = False
        for sobj in sobjs:
            if sobj['VEL_CORR'] is not None:
                # Don't double correct
                msg = f"Not performing {refframe} correction for {spec1d} object {sobj['NAME']} because it has already been corrected."
                log.info(msg)
                spec1d_failure_log.append(msg)
                continue
            # Use the SpecObj RA/DEC if it's available, otherwise use the value from the header
            if sobj['RA'] is not None and sobj['DEC'] is not None:
                radec = radec_to_coord((sobj['RA'], sobj['DEC']))
            else:
                radec = hdr_radec
            # Calculate the correction
            vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                   spectrograph.telescope['longitude'],
                                                   spectrograph.telescope['latitude'],
                                                   spectrograph.telescope['elevation'],
                                                   refframe)
            # Apply correction to objects
            log.info(f'Applying {refframe} correction to {spec1d} object {sobj["NAME"]} = {vel} km/s, {vel_corr}')
            sobj.apply_helio(vel_corr, refframe)
            corrected_at_least_one = True
        if corrected_at_least_one:
            sobjs.write_to_fits(subheader = sobjs.header, outfile=spec1d, overwrite=True)

def copy_spec1d_to_outdir(spec1d_files, outdir):
    """Copy the spec1d files to the requested outdir, preserving the originals

    Args:
        spec1d_files (list of str): List of spec1d files generated by PypeIt.
        outdir (str): Directory to copy the spec1d files.

    Return:
        list of str: The pathnames of the newly copied files.
    """

    # Make sure the spec1d output directory exists
    os.makedirs(outdir, exist_ok=True)
    
    copied_files = []
    for spec1d_file in spec1d_files:
        new_file = os.path.join(outdir, os.path.basename(spec1d_file))
        shutil.copy2(spec1d_file, new_file)
        copied_files.append(new_file)

    return copied_files


def coadd(par, coaddfile, source):
    """coadd the spectra for a given source.

    Args:
        par (:class:`~pypeit.par.pypeitpar.Collate1DPar`):
            Parameters for the coadding
        source (:class:`~pypeit.core.collate.SourceObject`):
            The SourceObject with information on which files and spectra to
            coadd.
    """
    # Set destination file for coadding
    par['coadd1d']['coaddfile'] = coaddfile
    
    # Determine if we should coadd flux calibrated data
    flux_key = par['coadd1d']['ex_value'] + "_FLAM"

    if par['collate1d']['ignore_flux'] is True:
        # Use non fluxed if asked to
        log.info(f"Ignoring flux for {coaddfile}.")
        par['coadd1d']['flux_value'] = False

    elif False in [x[flux_key] is not None for x in source.spec_obj_list]:               
        # Do not use fluxed data if one or more objects have not been flux calibrated 
        log.info(f"Not all spec1ds for {coaddfile} are flux calibrated, using counts instead.")
        par['coadd1d']['flux_value'] = False
    
    else:
        # Use fluxed data
        log.info(f"Using flux for {coaddfile}.")
        par['coadd1d']['flux_value'] = True


    # Instantiate
    spectrograph = load_spectrograph(par['rdx']['spectrograph'], pypeit_fits=True)
    coAdd1d = coadd1d.CoAdd1D.get_instance(source.spec1d_file_list,
                                           [x.NAME for x in source.spec_obj_list],
                                           spectrograph=spectrograph, par=par['coadd1d'])

    # Run
    coAdd1d.run()
    # Save to file
    coAdd1d.save(coaddfile)

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
        basename, ext = os.path.splitext(filename[len('spec1d_'):])
        spec2d_file = str(outputfiles.coadd_output_file(path, basename, twod=True, ext=ext))

        if not os.path.exists(spec2d_file):
            raise PypeItError(f'Could not find matching spec2d file for {spec1d_file}')

        spec2d_files.append(spec2d_file)

    return spec2d_files


def write_warnings(excluded_obj_log, failed_source_log, spec1d_failure_log):
    """
    Write gathered warning messages to a `collate_warnings.txt` file.

    Args:
        excluded_obj_log (:obj:`list` of :obj:`str`): 
            Messages about which objects were excluded from collating and why.

        failed_source_log (:obj:`list` of :obj:`str`): 
            Messages about which objects failed coadding and why.

        spec1d_failure_log (:obj:`list` of :obj:`str`): 
            Messages about failures with spec1d files and why.

    """
    report_filename = str(outputPaths.collate / "collate_warnings.txt")

    with open(report_filename, "w") as f:
        print("pypeit_collate_1d warnings", file=f)

        if len(spec1d_failure_log) > 0:
            print("\nspec1d_* failures\n", file=f)
            for msg in spec1d_failure_log:
                print(msg, file=f)

        print("\nExcluded Objects:\n", file=f)
        for msg in excluded_obj_log:
            print(msg, file=f)

        print("\nFailed to Coadd:\n", file=f)
        for msg in failed_source_log:
            print(msg, file=f)


def collate_spectra_by_source(source_list, tolerance, unit=u.arcsec):
    """Given a list of spec1d files from PypeIt, group the spectra within the
    files by their source object. The grouping is done by comparing the 
    position of each spectra (using either pixel or RA/DEC) using a given tolerance.

    Args:
        source_list (list of :obj:`SourceObject`): A list of source objects, one
            SpecObj per object, ready for collation.
        tolerance (float): 
            Maximum distance that two spectra can be from each other to be 
            considered to be from the same source. Measured in floating
            point pixels or as an angular distance (see ``unit`` argument).
        unit (`astropy.units.Unit`_):
            Units of ``tolerance`` argument if match_type is 'ra/dec'. 
            Defaults to arcseconds. Ignored if match_type is 'pixel'.

    Returns:
        list: The collated spectra as SourceObjects.

    """

    collated_list = []
    for source in source_list:

        # Search for a collated SourceObject that matches this one.
        # If one can't be found, treat this as a new collated SourceObject.
        found = False
        for collated_source in collated_list:
            if  collated_source.match(source.spec_obj_list[0], 
                                      source.spec1d_header_list[0],
                                      tolerance, unit):

                collated_source.combine(source)
                found = True

        if not found:
            collated_list.append(copy.deepcopy(source))

    return collated_list

def collate_1d(par, spectrograph, tolerance, spec1d_files):

    # Filter out unwanted source objects based on our parameters.
    # First filter them out based on the exclude_slit_trace_bm parameter
    if len(par['collate1d']['exclude_slit_trace_bm']) > 0:
        spec2d_files = find_spec2d_from_spec1d(spec1d_files)
        exclude_map = find_slits_to_exclude(spec2d_files, par)
    else:
        spec2d_files = []
        exclude_map = dict()

    # Flux the spec1ds based on a archived sensfunc
    spec1d_failure_log = []
    copied_spec1d = False
    if par['collate1d']['flux'] and not par['collate1d']['dry_run']:
        if par['collate1d']['spec1d_outdir'] is not None:
            # Fluxing modifies the spec1d files, copy them to a new output directory
            # if requested
            spec1d_files = copy_spec1d_to_outdir(spec1d_files, par['collate1d']['spec1d_outdir'])
            copied_spec1d = True                
        spec1d_files = flux(par, spectrograph, spec1d_files, spec1d_failure_log)
    

    # Perform reference frame correction
    if par['collate1d']['refframe'] in  ['heliocentric', 'barycentric'] and not par['collate1d']['dry_run']:
        if not copied_spec1d and par['collate1d']['spec1d_outdir'] is not None:
            # Refframe correction modifies the spec1d files, copy them to a new output directory
            # if requested and fluxing hasn't already done so
            spec1d_files = copy_spec1d_to_outdir(spec1d_files, par['collate1d']['spec1d_outdir'])

        refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_log)

    # Read in the spec1d files        
    specobjs_to_coadd, spec1d_files = read_spec1d_files(par, spec1d_files, spec1d_failure_log)

    # Build source objects from spec1d file, this list is not collated 
    source_objects = SourceObject.build_source_objects(specobjs_to_coadd, spec1d_files,
                                                        par['collate1d']['match_using'])

    # Filter out unwanted SpecObj objects based on parameters 
    (objects_to_coadd, excluded_obj_log) = exclude_source_objects(source_objects, exclude_map, par)

    # Collate the spectra
    source_list = collate_spectra_by_source(objects_to_coadd, tolerance)

    # Coadd the spectra
    successful_source_list = []
    failed_source_log = []
    for source in source_list:

        coaddfile = str(outputPaths.collate / build_coadd_file_name(source))
        log.info(f'Creating {coaddfile} from the following sources:')
        for i in range(len(source.spec_obj_list)):
            log.info(f'    {source.spec1d_file_list[i]}: {source.spec_obj_list[i].NAME} '
                        f'({source.spec_obj_list[i].MASKDEF_OBJNAME})')

        # Exclude sources with a single object to coadd
        if len(source.spec_obj_list) == 1:
            excluded_obj_log.append(f"Excluding {source.spec_obj_list[0].NAME} in {source.spec1d_file_list[0]} because there's no other SpecObj to coadd with.")
            continue

        if not par['collate1d']['dry_run']:
            try:
                coadd(par, coaddfile, source)
                successful_source_list.append(source)
            except Exception:
                formatted_exception = traceback.format_exc()
                log.warning(formatted_exception)
                log.warning(f"Failed to coadd {coaddfile}, skipping")
                failed_source_log.append(f"Failed to coadd {coaddfile}:")
                failed_source_log.append(formatted_exception)

    # Create collate_report.dat
    archive = create_report_archive()
    archive.add(successful_source_list)
    archive.save()

    write_warnings(excluded_obj_log, failed_source_log,
                    spec1d_failure_log)


