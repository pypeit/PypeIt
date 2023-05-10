"""
Script for quick-look PypeIt reductions.

Use cases:

  #. User inputs N files: arc, flat, science(s)

  #. User inputs 1 or more science files for a fixed-format instrument (e.g.
     NIRES)

  #. User inputs 1 folder of files

  #. User inputs 1 folder of files including 1 new science frame

  #. User inputs an ASCII file of files

  #. User inputs 2 science files with A-B [and calibs or uses defaults]

  #. User inputs N science files with A and B (only), stacks all files at A and
     B independently, A-B, add pos+neg

  #. User inputs N science files with an arbitrary set of dither patterns that
     are encoded in the headers (e.g. MOSFIRE, currently this works for just one
     dither pattern, and that may be all we need). Total stack is computed

Notes with JFH:

  #. Label B images as "sky" for A-B redux

  #. Write spec2D A images to disk with a minus sign and call B

  #. Consider not writing out but return instead

.. include:: ../include/links.rst
"""
from pathlib import Path
import time
import datetime
import shutil
from copy import deepcopy

from IPython import embed

import numpy as np

import configobj

from astropy.table import Table

from pypeit.pypmsgs import PypeItError
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import metadata
from pypeit import io
from pypeit import inputfiles 
from pypeit import pypeit
from pypeit import coadd2d
from pypeit.par.pypeitpar import PypeItPar
from pypeit.calibframe import CalibFrame
from pypeit.core.parse import parse_binning
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs
from pypeit.slittrace import SlitTraceSet 

from pypeit.scripts.setup_coadd2d import SetupCoAdd2D
from pypeit.scripts.coadd_2dspec import CoAdd2DSpec
from pypeit.scripts.show_2dspec import Show2DSpec

def get_files(raw_files, raw_path):
    """
    Use the user-provided input to get the files to process.

    Args:
        raw_files (:obj:`list`):
            The list of strings parsed from the ``raw_files`` command line
            argument.  Can be None.
        raw_path (:obj:`str`):
            The path to the raw files parsed from the ``raw_path`` command line
            argument.

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
            files = inputfiles.grab_rawfiles(raw_paths=[raw_path], list_of_files=raw_files)
        except PypeItError as e:
            msgs.error('Unable to parse provided input files.  Check --raw_files and '
                       '--raw_path input.')
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


def quicklook_regroup(fitstbl):
    """
    Regroup frames for quick-look reductions.

    Restrictions/Alterations applied are:

        - Science and standard frames each must be observations of a single
          target (according to the ``'target'`` keyword).

        - For metadata tables where the background image *are not* set, all
          science and standard frames are each assigned to a single combination
          group, so that all observations are stacked.
        
        - For metadata tables where the background images *are* set (e.g., for
          difference-imaging dither sequences), all frames of a given offset are
          combined (the column ``dithoff`` *must* exist).  I.e., if there are
          multiple dither sequences peformed (e.g., ABBA, then ABBA), all with
          the same offset (e.g., +/- 5 arcsec), the image combination and
          background IDs are changed so that all the As are combined and all the
          Bs are combined before calculating the difference.

    **This function directly alters the input object!**

    Args:
        fitstbl (:class:~pypeit.metadata.PypeItMetaData`):
            Metadata table for frames to be processed.
    """
    comb_strt = 0
    for frametype in ['science', 'standard']:
        is_type = fitstbl.find_frames(frametype)
        if any(is_type):
            # All frames must be of the same target
            if 'target' in fitstbl.keys() \
                    and not all(fitstbl['target'][is_type] == fitstbl['target'][is_type][0]):
                msgs.error(f'All {frametype} frames must be of the same target.')

            # Regroup dithered observations so that all images at a unique
            # offset are combined.
            if 'bkg_id' in fitstbl.keys() and any(fitstbl['bkg_id'].data[is_type] != -1):
                if 'dithoff' not in fitstbl.keys():
                    msgs.error('CODING ERROR: Metadata does not include dithoff column!')
                # Group the unique dither positions
                dith, inv = np.unique(fitstbl['dithoff'].data[is_type], return_inverse=True)
                if len(dith) == 1:
                    msgs.warn('All exposures have the same offset!')
                    fitstbl['comb_id'][is_type] = comb_strt
                else:
                    # This creates comb+bkg pairs that match the absolute value of the offset
                    # - Find the unique absolute values
                    abs_dith, abs_inv = np.unique(np.absolute(dith), return_inverse=True)
                    # - Create the indices that will identify the background images
                    bkg_inv = np.arange(dith.size)
                    # - For each unique throw, swap the combination IDs for the background IDs
                    for i in range(abs_dith.size):
                        bkg_inv[abs_inv == i] = bkg_inv[abs_inv == i][::-1]
                    # - Get the background indices
                    bkg_inv = bkg_inv[inv]

                    if np.array_equal(bkg_inv, inv):
                        # All of the offsets are unique!  Do the "dumb" thing of
                        # just setting the background image for the first offset
                        # to the second one, and setting the background for all
                        # other offsets to be the first one.
                        fitstbl['bkg_id'][np.where(is_type)[0][inv == comb_strt]] = comb_strt+1
                        fitstbl['bkg_id'][np.where(is_type)[0][inv != comb_strt]] = comb_strt
                    else:
                        same_frame = np.equal(bkg_inv, inv)
                        if any(same_frame):
                            # One or more of the offsets has no matched
                            # absolute-value pair.  Just use one of the other
                            # background images.
                            indx = np.where(np.logical_not(same_frame))[0][-1]
                            bkg_inv[same_frame] = bkg_inv[indx]
                        # All the offsets have a matched absolute value pair
                        fitstbl['comb_id'][is_type] = inv + comb_strt
                        fitstbl['bkg_id'][is_type] = bkg_inv + comb_strt
            else:
                # Force all images to be combined
                fitstbl['comb_id'][is_type] = comb_strt
            # Prep for the next image type (only if there are standards!)
            comb_strt = np.amax(fitstbl['comb_id'].data[is_type])+1


def generate_sci_pypeitfile(redux_path:str, 
                            ref_calib_dir:Path,
                            ps_sci, 
                            det:str=None,
                            clear:bool=False,
                            slitspatnum:str=None,
                            maskID:str=None,
                            boxcar_radius:float=None,
                            snr_thresh:float=None):
    """
    Prepare to reduce the science frames:

        - Correct the setup and calibration group for the science frames to be
          the same as the associated calibration files.

        - Create the path for the science reductions, and including a symlink to
          the pre-processed (reference) calibration frames.

        - Write the pypeit file with the requested parameter adjustments.
    
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
        clear (:obj:`bool`, optional):
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
    if sci_dir.exists() and clear:
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
    if len(is_std) > 0 and not clear:
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
                      're-reduced, use the --clear_science option.')
            ps_sci.remove_table_rows(is_std)

    # TODO: Push this stuff into a function in par/pypeitpar.py?

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


def calib_manifest(calib_dir, spectrograph):
    """
    Collate the list of setups with processed calibrations.

    The calibrations must exist within the provided parent directory.  The
    calibration directories must start with the PypeIt-specific name used for
    the relevant spectrograph (e.g., ``shane_kast_blue``), and the
    sub-directories must have a pypeit file that provides the instrument
    configuration (setup).  Directories that contain a pypeit file but no
    calibration subdirectory are ignored.
    
    Args:
        calib_dir (:obj:`str`, `Path`_):
            Parent directory with the calibrations; see above.  This directory
            must exist.
        spectrograph (:obj:`str`):
            The PypeIt-specific name of the spectrograph.

    Returns:
        :obj:`dict`: Dictionary with the list of instruments setups with
        available calibrations.  If no calibrations are found, either because no
        pypeit files are found or none of the directories with pypeit files
        include a calibrations directory, None is returned.  The dictionary has
        one item for every setup found.  Each item includes all the
        configuration parameters used to define a unique setup for this
        instrument, and it also includes the item ``'calib_dir'`` giving the
        directory with the processed calibrations.
    """
    # Check the calibration directory exists
    _calib_dir = Path(calib_dir).resolve()
    if not _calib_dir.exists():
        return None

    # Find the pypeit files
    pypeit_files = sorted(_calib_dir.glob('{0}_*/{0}*.pypeit'.format(spectrograph)))
    if len(pypeit_files) == 0:
        return None

    # Build the dictionary
    setups = {}
    for pypeit_file in pypeit_files:
        # Read
        pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)
        # Copy the setup
        setups[pypeitFile.setup_name] = deepcopy(pypeitFile.setup)
        # Remove the 'Setup *' entry
        del setups[pypeitFile.setup_name][f'Setup {pypeitFile.setup_name}']
        # Add the calibrations directory
        par = PypeItPar.from_cfg_lines(pypeitFile.cfg_lines)
        setups[pypeitFile.setup_name]['calib_dir'] \
                = pypeit_file.parent / par['calibrations']['calib_dir']
        # If the calibrations directory doesn't exist, ignore the directory!
        if not setups[pypeitFile.setup_name]['calib_dir'].exists():
            del setups[pypeitFile.setup_name]
    return setups


def match_to_calibs(ps:pypeitsetup.PypeItSetup, calib_dir:str, calibrated_setups=None):
    """
    Match observations to existing calibrations.

    The function first builds the list of existing calibration setups; see
    :func:`calib_manifest`.  Then, each setup for data provided by the
    :class:`~pypeit.pypeitsetup.PypeItSetup` object is matched to the existing
    setups.

    Args:
        ps (:class:`~pypeit.pypeitsetup.PypeItSetup`):
            Object providing metadata and parameters necessary to execute PypeIt
            data reduction.  This must contain data from a *single*
            setup/configuration; an error is raised if not.
        calib_dir (:obj:`str`):
            Parent directory with the calibrations; see :func:`calib_manifest`.
            This directory must exist.

    Returns:
        :obj:`dict`: A nested dictionary with the matched setups.  There is one
        dictionary key for each setup in ``ps``.  Each setup dictionary item is
        itself a dictionary with two items, ``setup`` and ``calib_dir``, where
        ``setup`` is setup identifier of the existing calibrations (as opposed
        to the setups provided in ``ps``) and the `Path`_ to the calibrations.
        If a setup in ``ps`` is not matched, the associated dictionary item is
        None.
    """
    matched_configs = {}
    if calibrated_setups is None:
        # Try to build the calibration manifest
        calibrated_setups = calib_manifest(calib_dir, ps.spectrograph.name)
    if calibrated_setups is None:
        # If it still doesn't exist, return empty dictionary
        return matched_configs

    for setup, config in ps.fitstbl.configs.items():
        if calibrated_setups is None:
            matched_configs[setup] = None
            continue
        matched_configs[setup] = dict(setup=[], calib_dir=[])
        for cal_setup, cal_config in calibrated_setups.items():
            if not ps.spectrograph.same_configuration([config, cal_config]):
                continue

            matched_configs[setup]['setup'] += [cal_setup]
            # NOTE: calib_manifest() checks that the calibration directory
            # exists for *any* returned setup.
            matched_configs[setup]['calib_dir'] += [cal_config['calib_dir']]

        if len(matched_configs[setup]['setup']) == 0:
            matched_configs[setup] = None
            continue
        elif len(matched_configs[setup]['setup']) > 1:
            msgs.warn('Existing calibrations have degenerate configurations!  We recommend you '
                      'clean your calibrations parent directory.  For now, using the first match.')
        matched_configs[setup]['setup'] = matched_configs[setup]['setup'][0]
        matched_configs[setup]['calib_dir'] = matched_configs[setup]['calib_dir'][0]

    return matched_configs


def merge_setups(calibrated_setups, calib_match, frame_setup):
    """
    Merge the setup identifiers for new reductions to those matched to existing
    reductions.

    Args:
        calibrated_setups (array-like):
            The list of setup identifiers with existing calibrations.  This can
            be the list of keys in the top-level dictionary returned by 
            :func:`calib_manifest`.  This is used to determine the next
            identifier for new setups to be calibrated.  *All* existing setups
            and every setup in ``calib_match`` must be included in this list.
        calib_match (:obj:`dict`):
            A dictionary that provides, for each new setup to be reduced, the
            matched set of existing calibrations, if there are any.  See
            :func:`match_to_calibs`.  The identifier for any matched setup must
            be included in ``calibrated_setups``.
        frame_setup (array-like):
            The setups associated with each frame in the new set of reductions.

    Returns:
        :obj:`list`: The list of replacement setups for each frame that
        consolidates existing calibrations with new ones.
    """
    # These are the identifiers that have been assigned to the new setups
    new_setup = np.array(list(calib_match.keys()), dtype=object)
    # These are the existing setups.  Any value here that is None means we need
    # to get a new setup identifier.
    old_setup = np.array([None if calib_match[setup] is None else calib_match[setup]['setup'] 
                            for setup in new_setup], dtype=object)
    # Set any new setup identifiers so that they don't conflict with existing ones
    gen = metadata.PypeItMetaData.configuration_generator()
    while None in old_setup:
        setup = next(gen)
        if setup in calibrated_setups:
            continue
        old_setup[np.where(old_setup == None)[0][0]] = setup
    # Associate each frame with an existing setup or a new one
    _frame_setup = np.asarray(frame_setup, dtype=object)
    new_frame_setup = np.empty_like(_frame_setup)
    for i,setup in enumerate(new_setup):
        new_frame_setup[_frame_setup == setup] = old_setup[i]
    return new_frame_setup.tolist()


def get_setup_calib(calib_dir, calib_grp=None):
    """
    Using all of the files found in the provided directory, determine the setup
    and calibration group to use for the quicklook science frame(s).

    Args:
        calib_dir (:obj:`str`, `Path`_):
            Directory with the calibration files.
        calib_grp (:obj:`int`, optional):
            If the calibration directory contains results from more than one
            calibration group, this *must* be provided, specifying which
            calibration group should be used.

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
    if calib_grp is not None:
        if str(calib_grp) in unique_calibs:
            return setup, str(calib_grp)
        msgs.error(f'Selected calibration group {calib_grp} is not available in {_calib_dir}.  '
                   'Must select a valid group.  Directory currently contains the following '
                   f'calibration groups: {unique_calibs}')

    # Cannot determine which calibration group to use.
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
        # TODO: Allow user to specify target to be reduced?

        parser.add_argument('--sci_files', type=str, nargs='+',
                            help='A space-separated list of raw file names that are science '
                                 'exposures.  These files must *also* be in the list of raw '
                                 'files.  Use of this option overrides the automated PypeIt '
                                 'frame typing.  Should only be used of automatic frame typing '
                                 'fails or is undesirable.')

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
                                 'sub-directories.  If None, the code will try to find relevant '
                                 'calibrations in the parent_calib_dir.  If no calibrations exist '
                                 'in that directory that match the instrument setup/configuration '
                                 'of the provided data, the code will construct new calibrations '
                                 '(assuming relevant raw files are provided).')
        parser.add_argument('--clear_science', default=False, action='store_true',
                            help='Remove the existing output science directories to force a fresh '
                                 'reduction.  If False, any existing directory structure will '
                                 'remain, and any alterations to existing science files will '
                                 'follow the normal behavior of run_pypeit.')

        # TODO: Get rid of "calibs only"?  I.e., force user only provide
        # calibrations if they want to process only calibrations.  Maybe the
        # point is that the latter isn't clear for frames typed as
        # `arc,science,tilt`...
        parser.add_argument('--calibs_only', default=False, action='store_true',
                            help='Reduce only the calibrations?')

        parser.add_argument('--overwrite_calibs', default=False, action='store_true',
                            help='Re-process and overwrite any existing calibration files.')

        parser.add_argument('--det', type=str, nargs='+',
                            help='A space-separated set of detectors or detector mosaics to '
                                 'reduce.  By default, *all* detectors or default mosaics for '
                                 'this instrument will be reduced.  Detectors in a mosaic must '
                                 'be a mosaic "allowed" by PypeIt and should be provided as '
                                 'comma-separated integers (with no spaces).  For example, to '
                                 'separately reduce detectors 1 and 5 for Keck/DEIMOS, you would '
                                 'use --det 1 5; to reduce mosaics made up of detectors 1,5 and '
                                 '3,7, you would use --det 1,5 3,7')
        parser.add_argument('--slitspatnum', type=str,
                            help='Reduce the slit(s) as specified by the slitspatnum value(s)')
        parser.add_argument('--maskID', type=int,
                            help='Reduce the slit(s) as specified by the maskID value(s)')
        parser.add_argument('--boxcar_radius', type=float,
                            help='Set the radius for the boxcar extraction in arcseconds')
        parser.add_argument('--snr_thresh', default=None, type=float,
                            help='Change the default S/N threshold used during source detection')

        parser.add_argument('--ignore_std', default=False, action='store_true',
                            help='If standard star observations are automatically detected, '
                                 'ignore those frames.  Otherwise, they are included with the '
                                 'reduction of the science frames.')
        parser.add_argument('--skip_display', dest='show', default=True, action='store_false',
                            help='Run the quicklook without displaying any results.')

        # TODO: Add fluxing option?

        # Coadding options
        parser.add_argument('--coadd2d', default=False, action='store_true',
                            help='Perform default 2D coadding.')
        # TODO: Consolidate slitspatnum and only_slits!
        parser.add_argument('--only_slits', type=str, nargs='+',
                            help='If coadding, only coadd this space-separated set of slits.  If '
                                 'not provided, all slits are coadded.')
        parser.add_argument('--offsets', type=str, default=None,
                            help='If coadding, spatial offsets to apply to each image; see the '
                                 '[coadd2d][offsets] parameter.  Options are restricted here to '
                                 'either maskdef_offsets or auto.  If not specified, the '
                                 '(spectrograph-specific) default is used.')
        parser.add_argument('--weights', type=str, default=None,
                            help='If coadding, weights used to coadd images; see the '
                                 '[coadd2d][weights] parameter.  Options are restricted here to '
                                 'either uniform or auto.  If not specified, the '
                                 '(spectrograph-specific) default is used.')
        parser.add_argument('--spec_samp_fact', default=1.0, type=float,
                            help='If coadding, adjust the wavelength grid sampling by this '
                                 'factor.  For a finer grid, set value to <1.0; for coarser '
                                 'sampling, set value to >1.0).')
        parser.add_argument('--spat_samp_fact', default=1.0, type=float,
                            help='If coadding, adjust the spatial grid sampling by this '
                                 'factor.  For a finer grid, set value to <1.0; for coarser '
                                 'sampling, set value to >1.0).')

        return parser


    @staticmethod
    def main(args):

        tstart = time.perf_counter()

        # Parse the raw files
        files = get_files(args.raw_files, args.raw_path)
        if len(files) == 0:
            msgs.error('No files to read!  Check --raw_files and --raw_path input.')

        # TODO: Include an option to save the ingested file list as a PypeIt
        # RawFile that can be edited?
        # from pypeit.inputfiles import RawFiles
        # tbl = Table()
        # tbl['filename'] = [Path(r).resolve().name for r in files]
        # RawFiles(file_paths=[args.raw_path], data_table=tbl).write('test.rawfiles')

        # Run PypeIt Setup on all the files
        ps = pypeitsetup.PypeItSetup.from_rawfiles(files, args.spectrograph)
        ps.run(setup_only=True)

        # Find the raw science files
        sci_idx = ps.fitstbl.find_frames('science') if args.sci_files is None \
                        else np.in1d(ps.fitstbl['filename'].data, args.sci_files)
        # TODO: Allow for standard files to be identified?

        # Check for any untyped files (that have not been typed) as science
        # files.  NOTE: This works regardless of whether or not sci_files has
        # been defined.
        unknown_types = [t is None for t in ps.fitstbl['frametype']]
        if any(unknown_types & np.logical_not(sci_idx)):
            # TODO: Remove them and keep going instead?
            msgs.error('Could not determine frame types for the following files: ' +
                       ', '.join(ps.fitstbl['filename'][unknown_types & np.logical_not(sci_idx)]))

        # Include any standards? 
        # TODO: Standards are only reduced if there are also science frames?
        if any(sci_idx) and not args.ignore_std:
            std_idx = ps.fitstbl.find_frames('standard')
        else:
            std_idx = np.zeros(sci_idx.size, dtype=bool)

        # Set the directory with the calibrations
        setup_calib_dir = None if args.setup_calib_dir is None \
                            else Path(args.setup_calib_dir).resolve()

        # For any science files, independently prep their meta data and
        # associate them with the correct calibrations, if necessary.
        if any(sci_idx):
            # Generate science setup object
            sci_files = ps.fitstbl.frame_paths(sci_idx)
            use_type = 'arc,science,tilt' \
                            if 'arc' in ps.fitstbl['frametype'][sci_idx][0] else 'science'
            frametype = {Path(f).name : use_type for f in sci_files}
            if any(std_idx):
                std_files = ps.fitstbl.frame_paths(std_idx)
                use_type = 'arc,standard,tilt' \
                                if 'arc' in ps.fitstbl['frametype'][std_idx][0] else 'standard'
                frametype = {**frametype, **{Path(f).name : use_type for f in std_files}}
                sci_files += std_files
            ps_sci = pypeitsetup.PypeItSetup.from_rawfiles(sci_files, ps.spectrograph.name,
                                                           frametype=frametype)
            # NOTE: By default, groupings=True in PypeItSetup.run, meaning that
            # the combination groups should be automatically set.  Setting the
            # combination groups is needed when standards are included and to
            # automatically identify dithered observations for some
            # spectrographs.
            ps_sci.run(setup_only=True)

            # TODO: Ensure that the file types are correct

            # Limit to a single setup
            if len(ps_sci.fitstbl.configs.keys()) > 1:
                msgs.error('Your science/standard files come from more than one setup.  Try '
                           'either ignoring the standard frames (if any are present and '
                           'auto-detected) and/or changing the list of science files.')

            # Regroup the image combination so that all images are combined.
            # NOTE: This will fault if all of the images are not of the same
            # target (based on the target in the fitstbl).
            # TODO: Should we instead do this in PypeItSetup?  I.e., should
            # PypeItSetup have an option that allows *all* images taken with a
            # given offset to be combined, and an option to only perform, e.g.,
            # A-B reduction, instead of the A-B + B-A.
            quicklook_regroup(ps_sci.fitstbl)

            # Setting the combination groups should also set the background IDs
            # for dithered observations.  That means we can automatically detect
            # whether or not the reductions use difference imaging
            # TODO: This is now the only place bkg_redux is used...
            bkg_redux = 'bkg_id' in ps_sci.fitstbl.keys() and any(ps_sci.fitstbl['bkg_id'] != -1)
            if bkg_redux:
                msgs.warn('Dither pattern automatically detected for these observations.  Image '
                          'combination and background subtraction sequences automatically set; '
                          'confirm the behavior is what you want by checking the auto-generated '
                          'pypeit file.')
                # Print the automatically detected dither pattern
                sci_only_idx = ps.fitstbl.find_frames('science', index=True)
                # Binning
                binspectral, binspatial = parse_binning(ps.fitstbl['binning'][sci_only_idx[0]])
                # Plate scale
                # NOTE: Assumes that the platescale does not change between
                # detectors or between observations!
                platescale = ps.spectrograph.get_detector_par(1)['platescale']*binspatial
                # Report
                print_offset_report(ps.fitstbl[sci_only_idx], platescale)

            # Check that all the frames are assigned to the same calibration group
            # NOTE: This is largely superfluous given the use of get_setup_calib
            # in generate_sci_pypeitfile, but it's useful to keep the warning
            # here.
            if any(ps_sci.fitstbl['calib'] != ps_sci.fitstbl['calib'][0]):
                msgs.warn('Automated configuration assigned multiple calibration groups to your '
                          'science frames.  Ignoring!  Assigning all frames to the same group.')
                ps_sci.fitstbl['calib'] = ps_sci.fitstbl['calib'][0]

            if args.parent_calib_dir is not None and setup_calib_dir is None:
                # The parent directory has been defined, so try to find a
                # relevant set of existing calibrations
                setup_calib_dir = match_to_calibs(ps_sci, args.parent_calib_dir)
                if setup_calib_dir is None:
                    # TODO: Fault here, or keep going to the next step, which is
                    # to try to build the calibrations?
                    msgs.error('No calibrations exist or could not find appropriate setup match '
                               f'in provided parent directory: {args.parent_calib_dir}')
                # NOTE: Code above check that there is only one setup in ps_sci
                setup_calib_dir = setup_calib_dir[ps_sci.fitstbl['setup'][0]]['calib_dir']
                msgs.info(f'Attempting to use archived calibrations found in {setup_calib_dir}.')

        elif not args.calibs_only:
            msgs.warn('No science frames found among the files provided.  Will only process '
                      'calibration frames.  If you have provided science frames, you can specify '
                      'which ones they are using the --sci_files option.')

        # TODO: What happens if there are *only* science frames and no matching
        # calibrations?

        # Calibrate, if necessary
        if setup_calib_dir is None:
            msgs.info('Building the processed calibration frames.')
            # Set the parent directory
            parent_calib_dir = args.redux_path if args.parent_calib_dir is None \
                                    else args.parent_calib_dir

            # Find all of the existing calibrations for this setup in the parent
            # directory
            calibrated_setups = calib_manifest(parent_calib_dir, ps.spectrograph.name)
            if calibrated_setups is not None:
                # Match the setups of the frames to be reduced to the existing
                # setups
                calib_match = match_to_calibs(ps, parent_calib_dir,
                                              calibrated_setups=calibrated_setups)
                # Edit the setup column in the metadata table so that existing
                # setups are updated and new setups are reduced with a unique
                # identifier, following the nominal sequence
                ps.fitstbl['setup'] = merge_setups(list(calibrated_setups.keys()), calib_match,
                                                   ps.fitstbl['setup'])
                # Force an update of the configurations
                ps.fitstbl.unique_configurations(force=True)

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
                                       reuse_calibs=not args.overwrite_calibs,
                                       calib_only=True)
                pypeIt.calib_all()

        if args.calibs_only or not any(sci_idx):
            msgs.info('Only calibrations exist or requested calibration processing only.  Done.')
            return

        # Build the PypeIt file for the science frames and link to the existing
        # calibrations.
        sci_pypeit_file = generate_sci_pypeitfile(
                    args.redux_path, 
                    setup_calib_dir, 
                    ps_sci,
                    det=args.det,
                    clear=args.clear_science,
                    slitspatnum=args.slitspatnum,
                    maskID=args.maskID, 
                    boxcar_radius=args.boxcar_radius,
                    snr_thresh=args.snr_thresh)

        # Run it
        pypeIt = pypeit.PypeIt(sci_pypeit_file, reuse_calibs=True)
        pypeIt.reduce_all()
        pypeIt.build_qa()

        # TODO: There are currently no tests for this!!
        # Perform coadding if requested
        if args.coadd2d:
            # Run the setup script to get the baseline coadd2d file
            # NOTE: By running the setup script with ``--keep_par`` flag,
            # parameters that skip the sky subtraction and the optimal
            # extraction should be propagated to the 2D coadding...
            # TODO:
            #   - Add sensitivity function?
            command_line_args = ['-f', sci_pypeit_file, '--keep_par']
            if args.only_slits is not None:
                command_line_args += ['--only_slits'] + args.only_slits
            if args.offsets is not None:
                command_line_args += ['--offsets', args.offsets]
            if args.weights is not None:
                command_line_args += ['--weights', args.weights]
            SetupCoAdd2D.main(SetupCoAdd2D.parse_args(command_line_args))

            # Find all the coadd2d scripts
            # NOTE: This should only find *one* coadd2d file because quick-look
            # should be limited to performing reduction for one target at a
            # time.
            coadd_file = sorted(Path(sci_pypeit_file).resolve().parent.glob('*.coadd2d'))
            if len(coadd_file) != 1:
                msgs.error('There should be only one 2D coadd file.')
            coadd_file = coadd_file[0]
            
            # Run the coadding
            coadd2dFile = inputfiles.Coadd2DFile.from_file(coadd_file)
            CoAdd2DSpec.main(CoAdd2DSpec.parse_args([str(coadd_file),
                                                     '--spec_samp_fact', str(args.spec_samp_fact),
                                                     '--spat_samp_fact', str(args.spat_samp_fact)]))

            # Get the output file name
            spectrograph, par, _ = coadd2dFile.get_pypeitpar()
            spec2d_files = coadd2dFile.filenames
            coadd_scidir = Path(coadd2d.CoAdd2D.output_paths(spec2d_files, par)[0]).resolve()
            basename = coadd2d.CoAdd2D.default_basename(spec2d_files)
            spec2d_file = str(coadd_scidir / f'spec2d_{basename}.fits')
        else:
            # Grab the spec2d file (or at least the first one)
            frame = pypeIt.fitstbl.find_frames('science', index=True)[0]
            spec2d_file = pypeIt.spec_output_file(frame, twod=True)

        if args.show:
            # TODO: Need to parse detector here?
            Show2DSpec.main(Show2DSpec.parse_args([spec2d_file]))

        # TODO: 
        #   - Print a statement that allows users to copy-paste the correct
        #     pypeit_show_1dspec call?
        #   - Provide a summary of the objects (that's not buried in previous
        #     screen output)?

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

