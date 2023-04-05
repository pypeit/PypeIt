"""
Script for quick-look PypeIt reductions.

Use cases:

  1. User inputs N files: arc, flat, science(s)
  2. User inputs 1 or more science files for a fixed-format instrument (e.g. NIRES)
  3. User inputs 1 folder of files
  4. User inputs 1 folder of files including 1 new science frame
  5. User inputs an ASCII file of files
  6. User inputs 2 science files with A-B [and calibs or uses defaults]
  7. User inputs N science files with A and B (only), stacks all files at A and B independently, A-B, add pos+neg
  8. User inputs N science files with an arbitrary set of dither patterns that are encoded in the headers (e.g. MOSFIRE, currently this works for just one dither pattern, and that may be all we need). Total stack is computed

Notes with JFH:
  1. Label B images as "sky" for A-B redux
  2. Write spec2D A images to disk with a minus sign and call B
  3. Consider not writing out but return instead

.. include:: ../include/links.rst
"""
from pathlib import Path
import time
import datetime
import shutil

from IPython import embed

import numpy as np

import configobj

from astropy.table import Table

from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import io
from pypeit import pypeit
from pypeit.calibframe import CalibFrame
from pypeit.core.parse import parse_binning
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs
from pypeit.slittrace import SlitTraceSet 
from pypeit import inputfiles 


def folder_name_from_scifiles(sci_files:list):
    """
    Construct the directory name for QL output.

    If one file, use the filename without its suffix (see
    :func:`~pypeit.io.remove_suffix`).  If multiple files, use a dash-separated
    string of the first and last filenames (also without their suffixes).

    Args:
        sci_files (:obj:`list`):
            List of :obj:`str` or `Path`_ objects with the science files.

    Returns:
        :obj:`str`: Directory name
    """
    if len(sci_files) == 1:
        return io.remove_suffix(sci_files[0])
    return '-'.join([io.remove_suffix(f) for f in [sci_files[0], sci_files[-1]]])


def generate_sci_pypeitfile(redux_path:str, 
                            ref_calib_dir:Path,
                            ps_sci, 
                            det:str=None,
                            remove_sci_dir:bool=True, 
                            slitspatnum:str=None,
                            maskID:str=None,
                            boxcar_radius:float=None,
                            bkg_redux:bool=False,
                            stack:bool=True):
    """
    Prepare to reduce the science frames by:

        - Correcting the setup and calibration group for the science frames to
          be the same as the associated calibration files.

        - Creating the path for the science reductions, and including a symlink
          to the pre-processed (reference) calibration frames.

        - Writing the pypeit file with the requested parameter adjustments.
    
    Args:
        redux_path (:obj:`str`):
            Path to the redux folder
        ref_calib_dir (`Path`_):
            Path with the pre-processed calibration frames.  A symlink will be
            created to this directory from within ``redux_path`` to mimic the
            location of the calibrations expected by :class:`~pypeit.PypeIt`.
        ps_sci (:class:`~pypeit.pypeitsetup.PypeItSetup`):
            Setup object for the science frame(s) only.
        det (:obj:`str`, optional):
            Detector/mosaic identifier.  If None, all detectors are reduced.
        remove_sci_dir (:obj:`bool`, optional):
            Remove the science directory if it exists.
        slitspatnum (:obj:`str`, optional):  
            Used to identify the slit that should be reduced; see
            :ref:`reduxpar`.  If None, all slits are reduced.
        maskID (:obj:`str`, optional):
            Slit identifier from the mask design, used to select a single slit
            to reduce.  If None, all slits are reduced.
        boxcar_radius (:obj:`float`, optional):
            Boxcar radius in arcsec used for extraction.  
        bkg_redux (:obj:`bool`, optional):
            Setup for dithered, difference-imaging reduction.
        stack (:obj:`bool`, optional):
            Reduce all of the science frames by stacking them all into a single
            image.

    Returns:
        :obj:`str`:  The name of the pypeit file.
    """

    # Check the directory with the reference calibrations exists
    if not ref_calib_dir.exists():
        msgs.error(f'Reference calibration directory does not exist: {ref_calib_dir}')

    # Get the setup and calibration group to use for the science frame(s)
    setup, calib = get_setup_calib(ref_calib_dir)

    # Force them in the fitstbl
    # TODO: Make this series of operations a function in PypeItMetaData?
    ps_sci.fitstbl['setup'] = setup
    ps_sci.fitstbl.unique_configurations(force=True)
    ps_sci.fitstbl['calib'] = calib
    ps_sci.fitstbl._set_calib_group_bits()

    # If stacking, make sure the combination group is added and set to the same
    # value for all science frames
    if len(ps_sci.fitstbl) > 1 and stack and not bkg_redux:
        ps_sci.fitstbl.set_combination_groups(assign_objects=False)
        ps_sci.fitstbl['comb_id'] = 1

    # Set the I/O directories
    _redux_path = Path(redux_path).resolve()
    sci_dir = _redux_path / folder_name_from_scifiles(ps_sci.fitstbl['filename'].data)
    calib_dir = sci_dir / ps_sci.par['calibrations']['calib_dir']

    # Science reduction folder
    if sci_dir.exists() and remove_sci_dir:
        shutil.rmtree(sci_dir)
    # NOTE: This is done here (as opposed to waiting for the pypeit.PypeIt call
    # to create it) so that the symlink to the calibrations directory can be
    # created.
    if not sci_dir.exists():
        sci_dir.mkdir(parents=True)

    # If the science directory wasn't removed and the new calibration directory
    # already exists, check that it points to the right directory.  If not,
    # raise an error.
    if calib_dir.exists() and calib_dir.is_link() and calib_dir.readlink() != ref_calib_dir:
        msgs.error(f'Symlink to calibrations directory ({calib_dir}) already exists and points '
                   f'to {calib_dir.readlink()} instead of {ref_calib_dir}.  Re-run quicklook '
                   f'forcing the existing reductions in {sci_dir} to be removed.')
    # Create the symlink if it doesn't already exist
    if not calib_dir.exists():
        calib_dir.symlink_to(ref_calib_dir, target_is_directory=True)
        
    # Configure
    cfg = {}
    cfg['rdx'] = {}
    cfg['rdx']['spectrograph'] = ps_sci.spectrograph.name
    cfg['rdx']['redux_path'] = str(sci_dir)
    cfg['rdx']['quicklook'] = True
    if det is not None:
        cfg['rdx']['detnum'] = det

    # TODO: Allow for input configurations?
#    if input_cfg_dict is not None:
#        full_cfg.merge(configobj.ConfigObj(input_cfg_dict))

    # Select to reduce a specific slit selected by its maskID
    if maskID is not None:
        # Get all the slit trace files in the calibrations directory
        slittrace_files = np.array(sorted(calib_dir.glob('Slits*')))
        # Filter based on the setup and calibration group
        keep = np.ones(slittrace_files.size, dtype=bool)
        for i, f in enumerate(slittrace_files):
            _setup, _calib, _detname = CalibFrame.parse_calib_key(
                    CalibFrame.parse_key_dir(str(f), from_filename=True)[0])
            keep[i] = _setup == setup and _calib in ['all', calib]
        if not any(keep):
            msgs.error('Could not find valid Slits calibration frame!')
        slittrace_files = slittrace_files[keep]

        # Iterate through each file to find the one with the relevant mask ID.
        detname = None
        for sliittrace_file in slittrace_files:
            slitTrace = SlitTraceSet.from_file(sliittrace_file)
            if maskID in slitTrace.maskdef_id:
                detname = slitTrace.detname
                # Mosaic?
                mosaic = True if detname[0:3] == 'MSC' else False
                det_id = np.where(ps_sci.spectrograph.list_detectors(mosaic=mosaic) == detname)[0]
                if len(det_id) == 0:
                    detname = None  # Reset
                    continue        # TODO: or fault?
                # Set det for reduction
                detnum = [ps_sci.spectrograph.allowed_mosaics[det_id[0]]] if mosaic else det_id[0]+1
                break
        if detname is None:
            msgs.error(f'Could not find a SlitTrace file with maskID={maskID}')

        # Add to config
        cfg['rdx']['detnum'] = detnum
        cfg['rdx']['maskIDs'] = maskID

    # Select to reduce a slit based on its spatial position on the detector?
    # TODO: Are slitspatnum and maskID mutually exclusive?  I.e., should this
    # use `elif` instead of `if`?
    if slitspatnum is not None:
        cfg['rdx']['slitspatnum'] = slitspatnum

    # Add the boxcar radius if specified
    if boxcar_radius is not None:
        cfg['redux'] = {'extraction': {'boxcar_radius': boxcar_radius}}

    # Write the pypeit file (note this returns the filename, not the list)
    return ps_sci.fitstbl.write_pypeit(output_path=sci_dir,
                                       cfg_lines=configobj.ConfigObj(cfg).write(),
                                       write_bkg_pairs=True, configs=setup)[0]


def match_science_to_calibs(ps_sci:pypeitsetup.PypeItSetup, 
                            calib_dir:str):
    """
    Match a set of science frames to the set of pre-made calibrations in the
    specified reduction folder, if any exists.
    
    Args:
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
            Object providing metadata and parameters necessary to execute PypeIt
            data reduction.  It is expected that this *only* includes science
            frames.
        calib_dir (:obj:`str`):
            Full path to the calibration directory.  This directory must exist.

    Returns:
        :obj:`tuple`: Provides the PypeIt file used to build the calibrations
        and the full path to the calibration directory.
    """
    # Check on one setup
    if len(ps_sci.fitstbl.configs.keys()) > 1:
        msgs.error('Your science files come from more than one setup. Please reduce them separately.')

    # Check the calibration directory exists
    _calib_dir = Path(calib_dir).resolve()
    if not _calib_dir.exists():
        msgs.error(f'Calibration directory does not exist: {_calib_dir}')

    # Check against existing calibration PypeIt files
    pypeit_files = sorted(_calib_dir.glob(
                        f'{ps_sci.spectrograph.name}_*/{ps_sci.spectrograph.name}_calib_*.pypeit'))
    mtch = []
    for pypeit_file in pypeit_files:
        # Read
        pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)

        # Check for a match
        match = True
        for key in ps_sci.spectrograph.configuration_keys():
            if ps_sci.fitstbl.configs['A'][key] != pypeitFile.setup[key]:
                match = False
                break

        # Found one, so add it
        if match:
            mtch.append(pypeit_file)

    # Are we ok?
    if len(mtch) == 0:
        msgs.error('Matched to zero setups.  The calibrations files are stale/wrong/corrupt.')
    elif len(mtch) > 1:
        msgs.warn('Matched to multiple setups.  This should not have happened, but we will take '
                  'the first.')

    # NOTE: Assumes ps_sci.par['calibrations']['calib_dir'] is the same as the
    # default, and the default was used when the calibrations were processed.
    return mtch[0], mtch[0].parent / ps_sci.par['calibrations']['calib_dir']


def get_setup_calib(calib_dir):
    """
    Using all of the files found in the provided directory, determine the setup
    and calibration group to use for the quicklook science frame(s).

    Args:
        calib_dir (:obj:`str`, `Path`_):
            Directory with the calibration files.

    Returns:
        :obj:`tuple`: The setup name and calibration group.
    """
    _calib_dir = Path(calib_dir).resolve()

    # Check there are files in the directory
    calib_files = sorted(_calib_dir.glob('*'))
    if len(calib_files) == 0:
        msgs.error(f'Calibrations directory is empty: {_calib_dir}')

    # For each file, try to parse the setup and calibration ID(s)
    setups = []
    calibs = []
    for f in calib_files:
        try:
            setup, calib, detname = CalibFrame.parse_calib_key(
                                        CalibFrame.parse_key_dir(str(f), from_filename=True)[0])
        except:
            continue
        setups += [setup]
        calibs += [calib]   # 'calib' is a comma-separated list of integers or 'all'

    # Find the unique setups
    setups = np.unique(setups)
    if len(setups) != 1:
        msgs.error(f'Calibration files for more than one setup found in {_calib_dir}, '
                    'according to their file names.  Calibration directory should only hold data '
                    'for *one* setup.')
    setup = setups[0]

    # Get the unique calibration groups identifiers
    unique_calibs = np.unique(np.concatenate([c.split(',') for c in calibs])).tolist()
    # Ignore any calibs assigned to all
    if 'all' in unique_calibs:
        unique_calibs.remove('all')
    if len(unique_calibs) == 0:
        # Everything assigned to all, so just use a default number
        return setup, '0'
    if len(unique_calibs) == 1:
        # There's only one calibration group, so use it.
        return setup, unique_calibs[0]

    # TODO: For now, we fault here.  For quick-look, this means there should
    # only be one calibration group, or all the calibrations are assigned to
    # 'all' calibration groups.  Otherwise, we could have the user supply the
    # calibration group.
    msgs.error(f'Calibrations in {_calib_dir} are part of multiple calibration groups.  Unclear '
               'how to proceed.')


class QL(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to produce quick-look PypeIt reductions',
                                    width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))

        parser.add_argument('--raw_files', type=str, nargs='+',
                            help='Either a PypeIt-formatted input file with the list of raw '
                                 'images to process and the relevant path, or a space-separated '
                                 'list of the filenames (e.g., "img1.fits img2.fits").  For the '
                                 'latter entry mode, the path containing the files is set using '
                                 '--raw_path.')
        parser.add_argument('--raw_path', type=str, default='current working directory',
                            help='Directory with the raw files to process.  Ignored if a '
                                 'PypeIt-formatted file is provided using the --rawfiles option.')
        parser.add_argument('--ext', type=str, default='.fits',
                            help='If raw file names are not provided directly using the '
                                 '--rawfiles option, this sets the extension used when searching '
                                 'for any files in the path defined by --raw_path.  All files '
                                 'found in the raw path with this extension will be processed.')

        parser.add_argument('--sci_files', type=str, nargs='+',
                            help='A space-separated list of raw file names that are science '
                                 'exposures.  These files must *also* be in the list of raw '
                                 'files.  Use of this option overrides the automated PypeIt '
                                 'frame typing.')

        parser.add_argument('--redux_path', type=str, default='current working directory',
                            help='Path for the QL reduction outputs.')

        parser.add_argument('--parent_calib_dir', type=str, 
                            help='Directory with/for calibrations for *all* instrument '
                                 'configurations/setups.  If provided, the data for your '
                                 'instrument configuration will be placed or pulled from a '
                                 'relevant sub-directory.  If None, the redux_path is used.')
        parser.add_argument('--setup_calib_dir', type=str, 
                            help='Directory with/for calibrations specific to your instrument '
                                 'configuration/setup.  Use of this option circumvents the '
                                 'automated naming system for the configuration/setup '
                                 'sub-directories.  If None, it is assumed that no calibrations '
                                 'exist and they must be created using the provided raw files.  '
                                 'The top-level directory is given by parent_calib_dir (or '
                                 'redux_path) and the sub-directories follow the normal PypeIt '
                                 'naming scheme.')

        # TODO: Add option to not overwrite existing reductions

        parser.add_argument('--calibs_only', default=False, action='store_true',
                            help='Reduce only the calibrations?')
        parser.add_argument('--overwrite_calibs', default=False, action='store_true',
                            help='Overwrite any existing calibration files?')
        parser.add_argument('--slitspatnum', type=str,
                            help='Reduce the slit(s) as specified by the slitspatnum value(s)')
        parser.add_argument('--maskID', type=int,
                            help='Reduce the slit(s) as specified by the maskID value(s)')
        parser.add_argument('--boxcar_radius', type=float,
                            help='Set the radius for the boxcar extraction in arcseconds')
        parser.add_argument('--det', type=str, nargs='+',
                            help='A space-separated set of detectors or detector mosaics to '
                                 'reduce.  By default, *all* detectors or default mosaics for '
                                 'this instrument will be reduced.  Detectors in a mosaic must '
                                 'be a mosaic "allowed" by PypeIt and should be provided as '
                                 'comma-separated integers (with no spaces).  For example, to '
                                 'separately reduce detectors 1 and 5 for Keck/DEIMOS, you would '
                                 'use --det 1 5; to reduce mosaics made up of detectors 1,5 and '
                                 '3,7, you would use --det 1,5 3,7')
        parser.add_argument('--no_stack', dest='stack', default=True, action='store_false',
                            help='Do *not* stack multiple science frames')
        parser.add_argument('--bkg_redux', default=False, action='store_true',
                            help='If set the script will perform difference imaging. Namely it '
                                 'will identify sequences of AB pairs based on the dither pattern '
                                 'and perform difference imaging sky subtraction and fit for '
                                 'residuals.  This functionality only works for instruments where '
                                 'PypeIt can automatically parse dither sequences from the file '
                                 'headers.')
        return parser


    @staticmethod
    def main(args):

        tstart = time.perf_counter()

        # Get the files to reduce
        files = None
        if args.raw_files is not None and len(args.raw_files) == 1:
            # Try to read the file names as if they're in a file_of_files
            try:
                files = io.grab_rawfiles(file_of_files=args.raw_files[0]) 
            except:
                # Failed, so assume its a raw file (only one calibration file?) and proceed
                pass
        if files is None: 
            try:
                files = io.grab_rawfiles(raw_paths=[args.raw_path], list_of_files=args.raw_files,
                                         extension=args.ext) 
            except PypeItError as e:
                msgs.error('Unable to parse provided input files.  Check --raw_files, --raw_path, '
                           'and/or --ext input.')
        if len(files) == 0:
            msgs.error('No files to read!  Check --raw_files, --raw_path, and/or --ext input.')

        # Include an option to save the ingested file list as a PypeIt RawFile
        # that can be edited?
        # from pypeit.inputfiles import RawFiles
        # tbl = Table()
        # tbl['filename'] = [Path(r).resolve().name for r in files]
        # RawFiles(file_paths=[args.raw_path], data_table=tbl).write('test.rawfiles')

        # Run PypeIt Setup on all the files
        ps = pypeitsetup.PypeItSetup.from_rawfiles(files, args.spectrograph)
        ps.run(setup_only=True)

        # Check for any untyped files
        unknown_types = [t is None for t in ps.fitstbl['frametype']]
        if any(unknown_types):
            # TODO: Remove them and keep going or fault?
            # TODO: Check these against ones that have been specified as science using 'sci_files'
            msgs.error('Could not determine frame types for the following files: ' +
                       ', '.join(ps.fitstbl['frametype'][unknown_types]))
        
        # Calibrate, if necessary
        parent_calib_dir = args.parent_calib_dir if args.parent_calib_dir is not None \
                                else args.redux_path
        if args.setup_calib_dir is None:
            # Generate PypeIt files (and folders)
            calib_pypeit_files = ps.generate_ql_calib_pypeit_files(
                parent_calib_dir, det=args.det, configs='all',
                overwrite=args.overwrite_calibs, bkg_redux=args.bkg_redux)

            # Process them
            for calib_pypeit_file in calib_pypeit_files: 
                # Path to PypeIt file
                redux_path = Path(calib_pypeit_file).resolve().parent
                # Check for existing masters
                calib_files = list((redux_path / ps.par['calibrations']['calib_dir']).glob('*'))
                if len(calib_files) > 0 and not args.overwrite_calibs:
                    msgs.info('Calibration files already exist.  Skipping calibration.')
                    continue

                # TODO: What happens for science frames that automatically get
                # typed as `arc,science,tilt` with separate calibration groups?

                # TODO: There may be a non-negligible overhead, but we should be
                # able to use Calibrations.get_association to get the list of
                # calibrations that *should* exist to figure out if we need to
                # do anything.  The above just checks for any calibration
                # files...

                # Run
                pypeIt = pypeit.PypeIt(calib_pypeit_file,
                                       redux_path=str(redux_path),
                                       calib_only=True)
                pypeIt.calib_all()

        if args.calibs_only:
            msgs.info('Calibrations only requested.  Exiting.')
            return

        # Science files
        # TODO: Include standards as well?  Add a command-line option?
        sci_idx = ps.fitstbl.find_frames('science') if args.sci_files is None \
                        else np.in1d(ps.fitstbl['filename'].data, args.sci_files)
        if np.sum(sci_idx) == 0:
            msgs.error('No science frames found in the provided files.  Add at least one or '
                       'specify using --sci_files.')

        # Dither pattern?
        if args.bkg_redux:
#            _det = None if args.det is None else [eval(d) for d in args.det]
#            detectors = pypeit.PypeIt.select_detectors(ps.spectrograph, _det, args.slitspatnum)
            # Binning
            binspectral, binspatial = parse_binning(ps.fitstbl['binning'][sci_idx][0])
            # Plate scale
            # NOTE: Assumes that the platescale does not change between
            # detectors or between observations!
            platescale = ps.spectrograph.get_detector_par(1)['platescale']*binspatial
            # Report
            print_offset_report(ps.fitstbl[sci_idx], platescale)

        # Generate science setup object
        # TODO: What happens for science frames that automatically get typed as
        # `arc,science,tilt` with separate calibration groups?
        sci_files = ps.fitstbl.frame_paths(sci_idx)
        ps_sci = pypeitsetup.PypeItSetup.from_rawfiles(sci_files, ps.spectrograph.name)
        ps_sci.run(setup_only=True)
        # Limit to a single science setup
        if len(ps_sci.fitstbl.configs.keys()) > 1:
            msgs.error('Your science files come from more than one setup.  They must be reduces '
                       'separately.')

        # Set the calibrations directory, either provided by the user or by
        # matching the available setups to that used for the science files
        setup_calib_dir = match_science_to_calibs(ps_sci, parent_calib_dir)[1] \
                            if args.setup_calib_dir is None \
                            else Path(args.setup_calib_dir).resolve()

        # Build the PypeIt file and link to Masters
        sci_pypeit_file = generate_sci_pypeitfile(
                    args.redux_path, 
                    setup_calib_dir, 
                    ps_sci, 
                    det=args.det,
                    slitspatnum=args.slitspatnum,
                    maskID=args.maskID, 
                    boxcar_radius=args.boxcar_radius,
                    bkg_redux=args.bkg_redux,
                    stack=args.stack)
        
        # Run it
        pypeIt = pypeit.PypeIt(sci_pypeit_file, reuse_calibs=True)
        pypeIt.reduce_all()
        pypeIt.build_qa()

        # COADD 2D GOES HERE
        exec_s = np.around(time.perf_counter()-tstart, decimals=1)
        msgs.info(f'Quicklook execution time: {datetime.timedelta(seconds=exec_s)}')


def print_offset_report(fitstbl:Table, platescale:float):
    """
    Print the dither pattern for a set of files to the screen.

    The files in the provided table must come from a single dither pattern!

    Args:
        fitstbl (`astropy.table.Table`_):
            Table with the file metadata.  This should be the table extracted
            from a :class:`~pypeit.metadata.PypeItMetaData` instance, selecting
            the relevant science frames.
        platescale (:obj:`float`):
            The relevant (binned) pixelscale in arcsec/pixel.
    """

    # Parse
    files = fitstbl['filename'].data
    offset_arcsec = fitstbl['dithoff'].data
    dither_pattern = fitstbl['dithpat'].data
    dither_id = fitstbl['dithpos'].data
    target = fitstbl['target'].data[0]

    # Proceed
    if len(np.unique(dither_pattern)) > 1:
        msgs.error('Script only supported for a single type of dither pattern.')

    # Print out a report on the offsets
    msg_string = msgs.newline() + '*******************************************************'
    msg_string += msgs.newline() + ' Summary of offsets for target {:s} with dither pattern:   {:s}'.format(target,
                                                                                                            dither_pattern[
                                                                                                                0])
    msg_string += msgs.newline() + '*******************************************************'
    msg_string += msgs.newline() + 'filename     Position         arcsec    pixels    '
    msg_string += msgs.newline() + '----------------------------------------------------'
    for iexp, file in enumerate(files):
        msg_string += msgs.newline() + '    {:s}    {:s}   {:6.2f}    {:6.2f}'.format(
            file, dither_id[iexp], offset_arcsec[iexp], offset_arcsec[iexp] / platescale)
    msg_string += msgs.newline() + '********************************************************'
    msgs.info(msg_string)

