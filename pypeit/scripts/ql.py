"""
Script for quick-look PypeIt reductions.

Use cases:

  1. DONE: User inputs N files: arc, flat, science(s)
  2. DONE: User inputs 1 or more science files for a fixed-format instrument (e.g. NIRES)
  3. DONE: User inputs 1 folder of files
  4. DONE: User inputs 1 folder of files including 1 new science frame
  5. DONE: User inputs an ASCII file of files
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

from pypeit.pypmsgs import PypeItError
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import io
from pypeit import pypeit
from pypeit.par.pypeitpar import PypeItPar
from pypeit.calibframe import CalibFrame
from pypeit.core.parse import parse_binning
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs
from pypeit.slittrace import SlitTraceSet 
from pypeit import inputfiles 

from pypeit.scripts.setup_coadd2d import SetupCoAdd2D
from pypeit.scripts.coadd_2dspec import CoAdd2DSpec

def get_files(raw_files, raw_path, ext):
    """
    Use the user-provided input to get the files to process.

    Args:
        raw_files (:obj:`list`):
            The list of strings parsed from the ``raw_files`` command line
            argument.  Can be None.
        raw_path (:obj:`str`):
            The path to the raw files parsed from the ``raw_path`` command line
            argument.
        ext (:obj:`str`):
            The file extension used to search for raw files parsed from the
            ``ext`` command line argument.

    Returns:
        :obj:`list`: List of strings providing the full path to each raw file.
    """
    # Get the files to reduce
    files = None
    if raw_files is not None and len(raw_files) == 1:
        # Try to read the file names as if they're in a file_of_files
        try:
            files = inputfiles.grab_rawfiles(file_of_files=raw_files[0]) 
        except:
            # Failed, so assume its a raw file (only one calibration file?) and proceed
            pass
    if files is None: 
        try:
            files = inputfiles.grab_rawfiles(raw_paths=[raw_path], list_of_files=raw_files,
                                             extension=ext) 
        except PypeItError as e:
            msgs.error('Unable to parse provided input files.  Check --raw_files, --raw_path, '
                        'and/or --ext input.')
    return files


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
                            clean:bool=False,
                            slitspatnum:str=None,
                            maskID:str=None,
                            boxcar_radius:float=None,
                            snr_thresh:float=None):
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
        clean (:obj:`bool`, optional):
            Remove the directory structure for these files; i.e., start a
            completely clean reduction.  If false, any existing directory
            structure will remain, but any existing science reduction products
            will still be overwritten.
        slitspatnum (:obj:`str`, optional):  
            Used to identify the slit that should be reduced; see
            :ref:`reduxpar`.  If None, all slits are reduced.
        maskID (:obj:`str`, optional):
            Slit identifier from the mask design, used to select a single slit
            to reduce.  If None, all slits are reduced.
        boxcar_radius (:obj:`float`, optional):
            Boxcar radius in arcsec used for extraction.  

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

    # Set the I/O directories
    _redux_path = Path(redux_path).resolve()
    sci_dir = _redux_path / folder_name_from_scifiles(ps_sci.fitstbl['filename'].data)
    calib_dir = sci_dir / ps_sci.par['calibrations']['calib_dir']

    # Science reduction folder
    if sci_dir.exists() and clean:
        shutil.rmtree(sci_dir)
    # NOTE: This is done here (as opposed to waiting for the pypeit.PypeIt call
    # to create it) so that the symlink to the calibrations directory can be
    # created.
    if not sci_dir.exists():
        sci_dir.mkdir(parents=True)

    # If the science directory wasn't removed and the new calibration directory
    # already exists, check that it points to the right directory.  If not,
    # raise an error.
    if calib_dir.exists() and calib_dir.is_symlink() and calib_dir.readlink() != ref_calib_dir:
        msgs.error(f'Symlink to calibrations directory ({calib_dir}) already exists and points '
                   f'to {calib_dir.readlink()} instead of {ref_calib_dir}.  Re-run quicklook '
                   f'forcing the existing reductions in {sci_dir} to be removed.')
    # Create the symlink if it doesn't already exist
    if not calib_dir.exists():
        calib_dir.symlink_to(ref_calib_dir, target_is_directory=True)

    # If a standard was provided, see if has been reduced already.  If so, use
    # the reduced standard and remove the standard from the setup object.
    std_spec1d = None
    is_std = ps_sci.fitstbl.find_frames('standard', index=True)
    if len(is_std) > 0 and not clean:
        for i in is_std:
            std_spec1d = pypeit.PypeIt.get_spec_file_name(
                            str(sci_dir / ps_sci.par['rdx']['scidir']),
                            ps_sci.fitstbl.construct_basename(i))
            if Path(std_spec1d).exists():
                break
            # File doesn't exist, so reset
            std_spec1d = None
        if std_spec1d is not None:
            # Found an existing reduction, so remove the standard frames.
            # NOTE: Should not need to regroup!
            msgs.warn(f'Found existing standard star reduction: {std_spec1d}.  This will be used '
                      'and the standards will not be re-reduced!  To force them to be '
                      're-reduced, use the --clean option.')
            ps_sci.remove_table_rows(is_std)

    # TODO: Enable code to find standard spec1d file in ref_calib_dir?  Something like:
#    elif find_archive_std:
#        std_files = sorted(calib_dir.glob('std_spec1d_*'))
#        if len(std_files) > 0:
#            std_spec1d = std_files[0]
    # TODO: Enable the user to provide the standard spec1d file directly?

    # TODO: Push this stuff into a function in par/pypeitpar.py

    # Get the quicklook parameters to use for this spectrograph
    spec_cfg = ps_sci.spectrograph.ql_par()

    # Configure
    cfg = {}
    cfg['rdx'] = {}
    cfg['rdx']['spectrograph'] = ps_sci.spectrograph.name
    cfg['rdx']['redux_path'] = str(sci_dir)
    cfg['rdx']['quicklook'] = True
    if det is not None:
        # TODO: Need to check that this works given how det is set by the
        # argument parser
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

    # Turn-off use of bias by default
    cfg['baseprocess'] = {}
    cfg['baseprocess']['use_biasimage'] = False

    # Add reduce dictionary?
    if any([k is not None for k in [snr_thresh, boxcar_radius, std_spec1d]]):
        cfg['reduce'] = {}
    if boxcar_radius is not None:
        cfg['reduce']['extraction'] = {'boxcar_radius': boxcar_radius}
    if snr_thresh is not None:
        cfg['reduce']['findobj'] = {'snr_thresh': snr_thresh}
    if std_spec1d is not None:
        cfg['reduce']['findobj'] = {'std_spec1d': std_spec1d}

    # Merge the spectrograph specific QL parameters with the ones above; the ones
    # above take precedence.
    cfg_lines = configobj.ConfigObj(spec_cfg)
    cfg_lines.merge(configobj.ConfigObj(cfg))
    cfg_lines = cfg_lines.write()

    # Write the pypeit file (note this returns the filename, not the list,
    # because of the [0] at the end of the call)
    return ps_sci.fitstbl.write_pypeit(output_path=sci_dir, cfg_lines=cfg_lines,
                                       write_bkg_pairs=True, configs=setup,
                                       config_subdir=False)[0]


def match_to_calibs(ps:pypeitsetup.PypeItSetup, calib_dir:str):
    """
    Match observations to a set of existing calibrations.

    The calibrations must exist within the provided parent directory.  The
    calibration directories must start with the PypeIt-specific name used for
    the relevant spectrograph (e.g., ``shane_kast_blue``), and the
    sub-directories must have a pypeit file that provides the instrument
    configuration (setup) to be matched against the provided observation
    metadata.
    
    Args:
        ps (:class:`pypeit.pypeitsetup.PypeItSetup`):
            Object providing metadata and parameters necessary to execute PypeIt
            data reduction.
        calib_dir (:obj:`str`):
            Parent directory with the calibrations; see above.  This directory
            must exist.

    Returns:
        `Path`_: Directory with existing calibrations that Provides the PypeIt file used to build the calibrations
        and the full path to the calibration directory.
    """
    # Check on one setup
    if len(ps.fitstbl.configs.keys()) > 1:
        msgs.error('Your observations use more than one setup.  They must be reduced separately.')

    # Get the setup identifier.  This should be 'A' because there is only one
    # setup, but this just ensures we're using the exact identifier.
    setup = list(ps.fitstbl.configs.keys())[0]

    # Check the calibration directory exists
    _calib_dir = Path(calib_dir).resolve()
    if not _calib_dir.exists():
        msgs.error(f'Calibration directory does not exist: {_calib_dir}')

    # Find the pypeit files
    pypeit_files = sorted(_calib_dir.glob('{0}_*/{0}*.pypeit'.format(ps.spectrograph.name)))

    if len(pypeit_files) == 0:
        msgs.error('Could not find any pypeit files!')

    matched_calib_dir = []
    for pypeit_file in pypeit_files:
        # Read
        pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)

        # Check for a match
        match = True
        for key in ps.spectrograph.configuration_keys():
            if ps.fitstbl.configs[setup][key] != pypeitFile.setup[key]:
                match = False
                break
        if not match:
            continue

        # Matched.  Use the pypeit file to set the directory with the
        # calibrations and only add it if the directory exists.
        par = PypeItPar.from_cfg_lines(pypeitFile.cfg_lines)
        _matched_dir = pypeit_file.parent / par['calibrations']['calib_dir']
        if _matched_dir.exists():
            matched_calib_dir += [_matched_dir]

    # Are we ok?
    if len(matched_calib_dir) == 0:
        msgs.error('Matched to zero setups.  The calibrations files are stale/wrong/corrupt.')
    elif len(matched_calib_dir) > 1:
        msgs.warn('Matched to multiple setups.  This should not have happened, but we will take '
                  'the first.')

    return matched_calib_dir[0]


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
        parser.add_argument('--clean', default=False, action='store_true',
                            help='Remove the existing output directories to force a fresh '
                                 'reduction.  If False, any existing directory structure will '
                                 'remain, but any existing science files will still be '
                                 'overwritten.')

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
        parser.add_argument('--ignore_std', default=False, action='store_true',
                            help='If standard star observations are automatically detected, '
                                 'ignore those frames.  Otherwise, they are included with the '
                                 'reduction of the science frames.')
        parser.add_argument('--snr_thresh', default=None, type=float,
                            help='Change the default S/N threshold used during source detection')

        parser.add_argument('--coadd', default=False, action='store_true',
                            help='Perform default 2D coadding.')
        # TODO: Add in relevant 2d coadding parameters, like detector, slits,
        # objects, etc.  2D coadding is basically just a default run of
        # `pypeit_setup_coadd2d` and `pypeit_coadd_2dspec`.

        return parser


    @staticmethod
    def main(args):

        tstart = time.perf_counter()

        #assert False

        # Parse the raw files
        files = get_files(args.raw_files, args.raw_path, args.ext)
        if len(files) == 0:
            msgs.error('No files to read!  Check --raw_files, --raw_path, and/or --ext input.')

        # TODO: Include an option to save the ingested file list as a PypeIt
        # RawFile that can be edited?
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

        # Find the raw science files
        sci_idx = ps.fitstbl.find_frames('science') if args.sci_files is None \
                        else np.in1d(ps.fitstbl['filename'].data, args.sci_files)
        if any(sci_idx) and not args.ignore_std:
            # Include standards
            sci_idx |= ps.fitstbl.find_frames('standard')

        # Set the directory with the calibrations
        setup_calib_dir = None if args.setup_calib_dir is None \
                            else Path(args.setup_calib_dir).resolve()

        # For any science files, independently prep their meta data and
        # associate them with the correct calibrations, if necessary.
        if any(sci_idx):
            # Generate science setup object
            sci_files = ps.fitstbl.frame_paths(sci_idx)
            ps_sci = pypeitsetup.PypeItSetup.from_rawfiles(sci_files, ps.spectrograph.name)
            # NOTE: By default, groupings=True in PypeItSetup.run, meaning that
            # the combination groups should be automatically set.  Setting the
            # combination groups is needed when standards are included and to
            # automatically identify dithered observations for some
            # spectrographs.
            ps_sci.run(setup_only=True)

            # Limit to a single setup
            if len(ps_sci.fitstbl.configs.keys()) > 1:
                msgs.error('Your science/standard files come from more than one setup.  Try '
                           'either ignoring the standard frames (if any are present and '
                           'auto-detected) and/or changing the list of science files.')

            # Setting the combination groups should also set the background IDs
            # for dithered observations.  That means we can automatically detect
            # whether or not the reductions use difference imaging
            bkg_redux = 'bkg_id' in ps_sci.fitstbl.keys() and any(ps_sci.fitstbl['bkg_id'] != -1)

            # If performing difference imaging, print the automatically detected
            # dither pattern
            if bkg_redux:
                sci_only_idx = ps.fitstbl.find_frames('science', index=True)
                # Binning
                binspectral, binspatial = parse_binning(ps.fitstbl['binning'][sci_only_idx[0]])
                # Plate scale
                # NOTE: Assumes that the platescale does not change between
                # detectors or between observations!
                platescale = ps.spectrograph.get_detector_par(1)['platescale']*binspatial
                # Report
                print_offset_report(ps.fitstbl[sci_only_idx], platescale)

            # Handle image stacking
            if args.stack:
                if bkg_redux:
                    msgs.warn('Dither pattern automatically detected for these observations.  '
                              'Ignoring request to stack all observations.  Image combination '
                              'and background subtraction sequences automatically set; confirm '
                              'the behavior is what you want by checking the auto-generated '
                              'pypeit file.')
                else:
                    # Stack all of the science and standard frames
                    ps_sci.fitstbl['comb_id'][ps_sci.fitstbl.find_frames('science')] = 1
                    ps_sci.fitstbl['comb_id'][ps_sci.fitstbl.find_frames('standard')] = 2

            # Check that all the frames are assigned to the same calibration group
            # NOTE: This is largely superfluous given the use of get_setup_calib
            # in generate_sci_pyepitfile, but it's useful to keep the warning
            # here.
            if any(ps_sci.fitstbl['calib'] != ps_sci.fitstbl['calib'][0]):
                msgs.warn('Automated configuration assigned multiple calibration groups to your '
                          'science frames.  Ignoring!  Assigning all frames to the same group.')
                ps_sci.fitstbl['calib'] = ps_sci.fitstbl['calib'][0]

            if args.parent_calib_dir is not None and setup_calib_dir is None:
                # The parent directory has been defined, so try to find a
                # relevant set of existing calibrations
                setup_calib_dir = match_to_calibs(ps_sci, args.parent_calib_dir)
                msgs.info(f'Attempting to use archived calibrations found in {setup_calib_dir}.')

        elif not args.calibs_only:
            msgs.warn('No science frames found among the files provided.  Will only process '
                      'calibration frames.  If you have provided science frames, you can specify '
                      'which ones they are using the --sci_files option.')

        # Calibrate, if necessary
        if setup_calib_dir is None:
            msgs.info('Building the processed calibration frames.')
            parent_calib_dir = args.redux_path if args.parent_calib_dir is None \
                                    else args.parent_calib_dir

            # Generate PypeIt files (and folders)
            calib_pypeit_files = ps.generate_ql_calib_pypeit_files(
                parent_calib_dir, det=args.det, configs='all',
                overwrite=args.overwrite_calibs)

            # Process them
            for calib_pypeit_file in calib_pypeit_files: 
                # Path to PypeIt file
                redux_path = Path(calib_pypeit_file).resolve().parent

                # Path for calibrations.
                # NOTE: When running with science frames, there will only be one
                # setup and one pypeit file.  So this sets the `setup_calib_dir`
                # used for the science frames.  If calibrations are present from
                # multiple setups, the code will have faulted above if there
                # were science frames included.  Otherwise, only the
                # calibrations will be processed and this object is only needed
                # locally, within this `for` loop.
                setup_calib_dir = redux_path / ps.par['calibrations']['calib_dir']

                # Check for existing calibrations
                # TODO: Although there may be a non-negligible overhead, we
                # should be able to use Calibrations.get_association to get the
                # list of calibrations that *should* exist to figure out if we
                # need to do anything.  Below just checks for *any* files in the
                # relevant directory.
                calib_files = list(setup_calib_dir.glob('*'))
                if len(calib_files) > 0 and not args.overwrite_calibs:
                    msgs.info('Calibration files already exist.  Skipping calibration.')
                    continue

                # Run
                pypeIt = pypeit.PypeIt(calib_pypeit_file,
                                       redux_path=str(redux_path),
                                       calib_only=True)
                pypeIt.calib_all()

        if args.calibs_only or not any(sci_idx):
            msgs.info('Only calibrations exist or request calibration processing only.  Exiting.')
            return

        # Build the PypeIt file for the science frames and link to the existing
        # calibrations.
        sci_pypeit_file = generate_sci_pypeitfile(
                    args.redux_path, 
                    setup_calib_dir, 
                    ps_sci, 
                    det=args.det,
                    clean=args.clean,
                    slitspatnum=args.slitspatnum,
                    maskID=args.maskID, 
                    boxcar_radius=args.boxcar_radius,
                    snr_thresh=args.snr_thresh)

#        from pypeit.inputfiles import PypeItFile
#        f = PypeItFile.from_file(sci_pypeit_file)
#        embed()
#        exit()

        # Run it
        pypeIt = pypeit.PypeIt(sci_pypeit_file, reuse_calibs=True)
        pypeIt.reduce_all()
        pypeIt.build_qa()

        # TODO: There are currently no tests for this!!
        # Perform coadding if requested
        if args.coadd:
            # Run the setup script to get the baseline coadd2d file
            # TODO: Add sensitivity functions, and other options
            SetupCoAdd2D.main(SetupCoAdd2D.parse_args([sci_pypeit_file]))

            # Find all the coadd2d scripts
            # TODO: Need to have SetupCoAdd2D.main return the names of the written files...
            coadd_files = sorted(Path(sci_pypeit_file).resolve().parent.glob('*.coadd2d'))
            
            # Run the coadding, only on those coadd files with more than one file
            for coadd_file in coadd_files:
                coadd2dFile = inputfiles.Coadd2DFile.from_file(coadd_file)
                if len(coadd2dFile.data) < 2:
                    msgs.warn(f'{coadd_file} only has one spec2d file.  Continuing...')
                    continue

                # TODO: Add options (e.g. spatial/spectral sampling...)
                CoAdd2DSpec.main(CoAdd2DSpec.parse_args([coadd_file]))

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

