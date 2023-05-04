"""
This script collates multiple 1d spectra in multiple files by object, 
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
from linetools import utils as ltu
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



def get_report_metadata(object_header_keys, spec_obj_keys, file_info):
    """
    Gets the metadata from a SourceObject instance used building a report
    on the results of collation.  It is intended to be wrapped in by functools
    partial object that passes in object_header_keys and spec_obj_keys. file_info
    is then passed as in by the :obj:`pypeit.archive.ArchiveMetadata` object.
    Unlike the other get_*_metadata functions, this is not used for archiving; it is
    used for reporting on the results of collating.

    If another type of file is added to the ArchiveMetadata object, the file_info
    argument will not be a SourceObject, In this case, a list of ``None`` values are 
    returned.

    Args:
        object_header_keys (list of str):
            The keys to read fom the spec1d headers from the SourceObject.

        spec_obj_keys (list of str):
            The keys to read from the (:obj:`pypeit.specobj.SpecObj`) objects in the SourceObject.

        file_info (:obj:`pypeit.scripts.collate_1d.SourceObject`)): 
            The source object containing the headers, filenames and SpecObj information for a coadd output file.

    Returns:
        tuple: A tuple of two lists:.

               **data_rows** (:obj:`list` of :obj:`list`): The metadata rows built from the source object.

               **files_to_copy** (iterable):
               An list of tuples of files to copy. Because this function is not used for
               archving data, this is always ``None``.
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
        par (:class:`~pypeit.par.pypeitpar.Collate1DPar`):
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

               **filtered_objects** (:obj:`list`): A list of :class:`~pypeit.core.collate.SourceObject` 
               with any excluded ones removed.

               **missing_archive_msgs** (:obj:`list`): A list of messages explaining why some source 
               objects were excluded.
    """
    filtered_objects = []
    excluded_messages= []
    for source_object in source_objects:

        sobj = source_object.spec_obj_list[0]
        spec1d_file = source_object.spec1d_file_list[0]

        if par['collate1d']['exclude_serendip'] and sobj.MASKDEF_OBJNAME == 'SERENDIP':
            msg = f'Excluding SERENDIP object from {sobj.NAME} in {spec1d_file}'
            msgs.info(msg)
            excluded_messages.append(msg)
            continue

        if par['collate1d']['wv_rms_thresh'] is not None and sobj.WAVE_RMS > par['collate1d']['wv_rms_thresh']:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} due to wave_rms {sobj.WAVE_RMS} > threshold {par["collate1d"]["wv_rms_thresh"]}'
            msgs.info(msg)
            excluded_messages.append(msg)
            continue

        if sobj.MASKDEF_ID in exclude_map:
            msg = f'Excluding {sobj.NAME} with mask id: {sobj.MASKDEF_ID} in {spec1d_file} because of flags {exclude_map[sobj.MASKDEF_ID]}'
            msgs.info(msg)
            excluded_messages.append(msg)
            continue

        if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing both OPT_COUNTS and BOX_COUNTS'
            msgs.warn(msg)
            excluded_messages.append(msg)
            continue

        if par['coadd1d']['ex_value'] == 'OPT' and sobj.OPT_COUNTS is None:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing OPT_COUNTS. Consider changing ex_value to "BOX".'
            msgs.warn(msg)
            excluded_messages.append(msg)
            continue

        if par['coadd1d']['ex_value'] == 'BOX' and sobj.BOX_COUNTS is None:
            msg = f'Excluding {sobj.NAME} in {spec1d_file} because of missing BOX_COUNTS. Consider changing ex_value to "OPT".'
            msgs.warn(msg)
            excluded_messages.append(msg)
            continue

        filtered_objects.append(source_object)
    return (filtered_objects, excluded_messages)

def flux(par, spectrograph, spec1d_files, failed_fluxing_msgs):
    """
    Flux calibrate spec1d files using archived sens func files.

    Args:
        par (`obj`:pypeit.par.pypeitpar.PypeItPar): 
            Parameters for collating, fluxing, and coadding.
        spectrograph (`obj`:pypeit.spectrographs.spectrograph):
            Spectrograph for the files to flux.
        spec1d_files (list of str):
            List of spec1d files to flux calibrate.
        failed_fluxing_msgs(list of str):
            Return parameter describing any failures that occurred when fluxing.

    Returns:
        list of str: The spec1d files that were successfully flux calibrated.
    """

    # Make sure fluxing from archive is supported for this spectrograph
    if spectrograph.name not in SensFileArchive.supported_spectrographs():
        msgs.error(f"Flux calibrating {spectrograph.name} with an archived sensfunc is not supported.")

    par['fluxcalib']['extrap_sens'] = True

    sf_archive = SensFileArchive.get_instance(spectrograph.name)
    flux_calibrated_files = []
    for spec1d_file in spec1d_files:

        # Get the archived sens file to use
        try:
            sens_file = sf_archive.get_archived_sensfile(spec1d_file)
        except Exception:
            formatted_exception = traceback.format_exc()
            msgs.warn(formatted_exception)
            msgs.warn(f"Could not find archived sensfunc to flux {spec1d_file}, skipping it.")
            failed_fluxing_msgs.append(f"Could not find archived sensfunc to flux {spec1d_file}, skipping it.")
            failed_fluxing_msgs.append(formatted_exception)
            
        # Flux calibrate the spec1d file
        try:
            msgs.info(f"Running flux calibrate on {spec1d_file}")
            FxCalib = fluxcalibrate.FluxCalibrate.get_instance([spec1d_file], [sens_file],
                                                                par=par['fluxcalib'])
            flux_calibrated_files.append(spec1d_file)

        except Exception:
            formatted_exception = traceback.format_exc()
            msgs.warn(formatted_exception)
            msgs.warn(f"Failed to flux calibrate {spec1d_file}, skipping it.")
            failed_fluxing_msgs.append(f"Failed to flux calibrate {spec1d_file}, skipping it.")
            failed_fluxing_msgs.append(formatted_exception)

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

def refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_msgs):

    refframe = par['collate1d']['refframe']
    msgs.info(f"Performing a {refframe} correction")

    for spec1d in spec1d_files:
        # Get values from the fits header needed to calculate the correction
        try:
            sobjs = SpecObjs.from_fitsfile(spec1d)
            hdr_ra = sobjs.header['RA']
            hdr_dec = sobjs.header['DEC']
            hdr_radec = ltu.radec_to_coord((hdr_ra, hdr_dec))
            obstime = Time(sobjs.header['MJD'], format='mjd')
        except Exception as e:
            msg = f'Failed to perform {refframe} correction on {spec1d}: {e}'
            msgs.info(msg)
            spec1d_failure_msgs.append(msg)
            continue

        corrected_at_least_one = False
        for sobj in sobjs:
            if sobj['VEL_CORR'] is not None:
                # Don't double correct
                msg = f"Not performing {refframe} correction for {spec1d} object {sobj['NAME']} because it has already been corrected."
                msgs.info(msg)
                spec1d_failure_msgs.append(msg)
                continue
            # Use the SpecObj RA/DEC if it's available, otherwise use the value from the header
            if sobj['RA'] is not None and sobj['DEC'] is not None:
                radec = ltu.radec_to_coord((sobj['RA'], sobj['DEC']))
            else:
                radec = hdr_radec
            # Calculate the correction
            vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                   spectrograph.telescope['longitude'],
                                                   spectrograph.telescope['latitude'],
                                                   spectrograph.telescope['elevation'],
                                                   refframe)
            # Apply correction to objects
            msgs.info(f'Applying {refframe} correction to {spec1d} object {sobj["NAME"]} = {vel} km/s, {vel_corr}')
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
        par (`obj`:Collate1DPar): Paramters for the coadding
        source (`obj`:SourceObject): The SourceObject with information on
            which files and spectra to coadd.
    """
    # Set destination file for coadding
    par['coadd1d']['coaddfile'] = coaddfile
    
    # Determine if we should coadd flux calibrated data
    flux_key = par['coadd1d']['ex_value'] + "_FLAM"

    if par['collate1d']['ignore_flux'] is True:
        # Use non fluxed if asked to
        msgs.info(f"Ignoring flux for {coaddfile}.")
        par['coadd1d']['flux_value'] = False

    elif False in [x[flux_key] is not None for x in source.spec_obj_list]:               
        # Do not use fluxed data if one or more objects have not been flux calibrated 
        msgs.info(f"Not all spec1ds for {coaddfile} are flux calibrated, using counts instead.")
        par['coadd1d']['flux_value'] = False
    
    else:
        # Use fluxed data
        msgs.info(f"Using flux for {coaddfile}.")
        par['coadd1d']['flux_value'] = True


    # Instantiate
    spectrograph = load_spectrograph(par['rdx']['spectrograph'])
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
        spec2d_file = os.path.join(path, filename.replace('spec1d', 'spec2d', 1))

        if not os.path.exists(spec2d_file):
            msgs.error(f'Could not find matching spec2d file for {spec1d_file}')

        spec2d_files.append(spec2d_file)

    return spec2d_files


def write_warnings(par, excluded_obj_msgs, failed_source_msgs, spec1d_failure_msgs, start_time, total_time):
    """
    Write gathered warning messages to a `collate_warnings.txt` file.

    Args:
        excluded_obj_msgs (:obj:`list` of :obj:`str`): 
            Messages about which objects were excluded from collating and why.

        failed_source_msgs (:obj:`list` of :obj:`str`): 
            Messages about which objects failed coadding and why.

        spec1d_failure_msgs (:obj:)`list` of :obj:`str`): 
            Messages about failures with spec1d files and why.

    """
    report_filename = os.path.join(par['collate1d']['outdir'], "collate_warnings.txt")

    with open(report_filename, "w") as f:
        print("pypeit_collate_1d warnings", file=f)
        print(f"\nStarted {start_time.isoformat(sep=' ')}", file=f)
        print(f"Duration: {total_time}", file=f)

        if len(spec1d_failure_msgs) > 0:
            print("\nspec1d_* failures\n", file=f)
            for msg in spec1d_failure_msgs:
                print(msg, file=f)

        print("\nExcluded Objects:\n", file=f)
        for msg in excluded_obj_msgs:
            print(msg, file=f)

        print("\nFailed to Coadd:\n", file=f)
        for msg in failed_source_msgs:
            print(msg, file=f)

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

    if args.exclude_slit_bm is not None and len(args.exclude_slit_bm) > 0:
        params['collate1d']['exclude_slit_trace_bm'] = args.exclude_slit_bm.split(',')

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

    return params, spectrograph, spec1d_files

def create_report_archive(par):
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

    COADDED_SOBJ_KEYS  =        ['MASKDEF_OBJNAME', 'MASKDEF_ID', 'NAME',        'DET', 'RA',    'DEC',    'S2N', 'MASKDEF_EXTRACT', 'WAVE_RMS']
    COADDED_SOBJ_COLUMN_NAMES = ['maskdef_objname', 'maskdef_id', 'pypeit_name', 'det', 'objra', 'objdec', 's2n',     'maskdef_extract', 'wave_rms']

    report_names = ['filename'] + \
                   COADDED_SOBJ_COLUMN_NAMES + \
                   ['spec1d_filename'] + \
                   COADDED_SPEC1D_COLUMN_NAMES

    report_formats = {'s2n':      '%.2f',
                      'wave_rms': '%.3f'}

    report_metadata = ArchiveMetadata(os.path.join(par['collate1d']['outdir'], "collate_report.dat"),
                                      report_names,
                                      partial(get_report_metadata,
                                              COADDED_SPEC1D_HEADER_KEYS,
                                              COADDED_SOBJ_KEYS),
                                      append=True,
                                      formats= report_formats)
    archive_metadata_list.append(report_metadata)

    # metadatas in archive object
    return ArchiveDir(par['collate1d']['outdir'], archive_metadata_list, copy_to_archive=False)


class Collate1D(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        # A blank Colate1DPar to avoid duplicating the help text.
        blank_par = pypeitpar.Collate1DPar()

        parser = super().get_parser(description='Flux/Coadd multiple 1d spectra from multiple '
                                                'nights and prepare a directory for the KOA.',
                                    width=width, formatter=scriptbase.SmartFormatter)

        # TODO: Is the file optional?  If so, shouldn't the first argument start
        # with '--'?
        parser.add_argument('input_file', type=str,
                            help='R|(Optional) File for guiding the collate process.  '
                                 'Parameters in this file are overidden by the command line. The '
                                 'file must have the following format:\n'
                                 '\n'
                                 'F|[collate1d]\n'
                                 'F|  tolerance             <tolerance>\n'
                                 'F|  outdir                <directory to place output files>\n'
                                 'F|  spec1d_outdir         <directory to place modified spec1ds, if any>\n'
                                 'F|  exclude_slit_trace_bm <slit types to exclude>\n'
                                 'F|  exclude_serendip      If set serendipitous objects are skipped.\n'  
                                 'F|  match_using           Whether to match using "pixel" or\n'
                                 'F|                        "ra/dec"\n'
                                 'F|  dry_run               If set the matches are displayed\n'
                                 'F|                        without any processing\n'
                                 'F|  flux                  Flux calibrate using archived sensfuncs.\n'
                                 'F|  ignore_flux           Ignore any flux calibration information in\n'
                                 'F|                        spec1d files.\n'
                                 'F|  wv_rms_thresh         If set, any objects with a wavelength rms > than the input\n'
                                 'F|                        value are skipped, else all wavelength rms values are accepted.\n'
                                 'F|  refframe              Perform reference frame correction prior to coadding.\n'
                                f'F|                        Options are {pypeitpar.WavelengthSolutionPar.valid_reference_frames()}. Defaults to None.\n'
                                 '\n'
                                 'F|spec1d read\n'
                                 'F|<path to spec1d files, wildcards allowed>\n'
                                 'F|...\n'
                                 'F|end\n',                        
                            nargs='?')
        parser.add_argument('--spec1d_files', type=str, nargs='*',
                            help='One or more spec1d files to flux/coadd/archive. '
                                 'Can contain wildcards')
        parser.add_argument('--par_outfile', default=None, type=str,
                            help='Output to save the parameters')
        parser.add_argument('--outdir', type=str, help=blank_par.descr['outdir'] + " Defaults to the current directory.")
        parser.add_argument('--spec1d_outdir', type=str, help=blank_par.descr['spec1d_outdir'] + " Defaults to overwriting existing spec1ds.")
        parser.add_argument('--tolerance', type=str, help=blank_par.descr['tolerance'])
        parser.add_argument('--match_using', type=str, choices=blank_par.options['match_using'],
                            help=blank_par.descr['match_using'])
        parser.add_argument('--dry_run', action='store_true', help=blank_par.descr['dry_run'])
        parser.add_argument('--ignore_flux', default=False, action='store_true', help=blank_par.descr['ignore_flux'])
        parser.add_argument('--flux', default=False, action = 'store_true', help=blank_par.descr['flux'])
        parser.add_argument('--exclude_slit_bm', type=str, # nargs='*',
                            help=blank_par.descr['exclude_slit_trace_bm']+' Comma separated.')
        parser.add_argument('--exclude_serendip', action='store_true',
                            help=blank_par.descr['exclude_serendip'])
        parser.add_argument("--wv_rms_thresh", type=float, default = None, help=blank_par.descr['wv_rms_thresh'])
        parser.add_argument("--refframe", type=str, default = None, choices = pypeitpar.WavelengthSolutionPar.valid_reference_frames(),
                            help=blank_par.descr['refframe'])
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename collate_1d_YYYYMMDD-HHMM.log')
        return parser

    @staticmethod
    def main(args):

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('collate_1d', args.verbosity)

        start_time = datetime.now()
        (par, spectrograph, spec1d_files) = build_parameters(args)

        outdir = par['collate1d']['outdir'] 
        os.makedirs(outdir, exist_ok=True)

        # Write the par to disk
        if args.par_outfile is None:
            args.par_outfile = os.path.join(outdir, 'collate1d.par')
        print("Writing the parameters to {}".format(args.par_outfile))
        # Gather up config lines for the sections relevant to collate_1d
        config_lines = par['collate1d'].to_config(section_name='collate1d',include_descr=False) + ['']
        config_lines += par['coadd1d'].to_config(section_name='coadd1d',include_descr=False)
        if par['collate1d']['flux']:
            config_lines += [''] + par['fluxcalib'].to_config(section_name='fluxcalib',include_descr=False)
        with open(args.par_outfile, "w") as f:
            for line in config_lines:
                print (line, file=f)

        # Parse the tolerance based on the match type
        if par['collate1d']['match_using'] == 'pixel':
            tolerance = float(par['collate1d']['tolerance'])
        else:
            # For ra/dec matching, the default unit is arcseconds. We check for
            # this case by seeing if the passed in tolerance is a floating point number
            if is_float(par['collate1d']['tolerance']):
                tolerance =  float(par['collate1d']['tolerance'])
            else:
                tolerance = Angle(par['collate1d']['tolerance']).arcsec

        # Filter out unwanted source objects based on our parameters.
        # First filter them out based on the exclude_slit_trace_bm parameter
        if len(par['collate1d']['exclude_slit_trace_bm']) > 0:
            spec2d_files = find_spec2d_from_spec1d(spec1d_files)
            exclude_map = find_slits_to_exclude(spec2d_files, par)
        else:
            spec2d_files = []
            exclude_map = dict()

        # Flux the spec1ds based on a archived sensfunc
        spec1d_failure_msgs = []
        copied_spec1d = False
        if par['collate1d']['flux']:
            if par['collate1d']['spec1d_outdir'] is not None:
                # Fluxing modifies the spec1d files, copy them to a new output directory
                # if requested
                spec1d_files = copy_spec1d_to_outdir(spec1d_files, par['collate1d']['spec1d_outdir'])
                copied_spec1d = True                
            spec1d_files = flux(par, spectrograph, spec1d_files, spec1d_failure_msgs)

        # Perform reference frame correction
        if par['collate1d']['refframe'] in  ['heliocentric', 'barycentric']:
            if not copied_spec1d and par['collate1d']['spec1d_outdir'] is not None:
                # Refframe correction modifies the spec1d files, copy them to a new output directory
                # if requested and fluxing hasn't already done so
                spec1d_files = copy_spec1d_to_outdir(spec1d_files, par['collate1d']['spec1d_outdir'])

            refframe_correction(par, spectrograph, spec1d_files, spec1d_failure_msgs)

        # Build source objects from spec1d file, this list is not collated 
        source_objects = SourceObject.build_source_objects(spec1d_files,
                                                           par['collate1d']['match_using'])

        # Filter based on the coadding ex_value, and the exclude_serendip 
        # boolean
        (objects_to_coadd, excluded_obj_msgs) = exclude_source_objects(source_objects, exclude_map, par)

        # Collate the spectra
        source_list = collate_spectra_by_source(objects_to_coadd, tolerance)

        # Coadd the spectra
        successful_source_list = []
        failed_source_msgs = []
        for source in source_list:
            coaddfile = os.path.join(par['collate1d']['outdir'], build_coadd_file_name(source))
            msgs.info(f'Creating {coaddfile} from the following sources:')
            for i in range(len(source.spec_obj_list)):
                msgs.info(f'    {source.spec1d_file_list[i]}: {source.spec_obj_list[i].NAME} '
                          f'({source.spec_obj_list[i].MASKDEF_OBJNAME})')

            if not args.dry_run:
                try:
                    coadd(par, coaddfile, source)
                    successful_source_list.append(source)
                except Exception:
                    formatted_exception = traceback.format_exc()
                    msgs.warn(formatted_exception)
                    msgs.warn(f"Failed to coadd {coaddfile}, skipping")
                    failed_source_msgs.append(f"Failed to coadd {coaddfile}:")
                    failed_source_msgs.append(formatted_exception)

        # Create collate_report.dat
        archive = create_report_archive(par)
        archive.add(successful_source_list)
        archive.save()

        total_time = datetime.now() - start_time

        write_warnings(par, excluded_obj_msgs, failed_source_msgs,
                       spec1d_failure_msgs, start_time, total_time)

        msgs.info(f'Total duration: {total_time}')

        return 0
