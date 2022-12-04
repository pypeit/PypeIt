"""
Script for quick-look reductions for Multislit observations.

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

import os
import time
import glob
import shutil

import numpy as np

import configobj

from pypeit import utils
from pypeit import masterframe
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import io
from pypeit import pypeit
from pypeit.scripts import scriptbase
from pypeit.spectrographs import available_spectrographs
from pypeit.slittrace import SlitTraceSet 
from pypeit import inputfiles 

from IPython import embed


def folder_name_from_scifiles(sci_files:list):
    """ Folder name for output of QL on science file(s)

    If one file, use the filename minus any extensions (e.g. .fits.gz)
    If multiple files, use a - separated string of the first
    and last filenames

    Args:
        sci_files (list): List of science files

    Returns:
        str: Folder name
    """
    def strip_me_down(ifile):
        base = os.path.basename(ifile)
        ipos = base.find('.fits') # Will this fail for any spectrograph?
        return base[:ipos]

    first_file = strip_me_down(sci_files[0])
    if len(sci_files) == 1:
        return first_file
    else:
        last_file = strip_me_down(sci_files[-1])
        return f'{first_file}-{last_file}'

def generate_sci_pypeitfile(redux_path:str, 
                            sci_files:list, 
                            master_calib_dir:str,
                            master_setup_and_bit:list,
                            ps_sci, 
                            det:str=None,
                            input_cfg_dict:dict=None, 
                            remove_sci_dir:bool=True, 
                            slitspatnum:str=None,
                            maskID:str=None,
                            boxcar_radius:float=None,
                            stack:bool=True):
    """
    Generate the PypeIt file for the science frames
    from the calib PypeIt file.

    The primary steps are:
      - Genreate the science reduction folder based on the science filenames
      - Generate a soft-link to the Masters/ folder provided by master_calib_dir
      - Build the configuration lines for the PypeIt file
      - Write the PypeIt file to disk in the science reduction folder
    
    Args:
        redux_path (str): Path to the redux folder
        sci_files (list): List of science files (full path)
        master_calib_dir (str): Path to the master calib folder
        master_setup_and_bit (list): 
            Name of the master setup and bit (list of str)
            The latter is used to tie the science frames to the Masters
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
            Setup object for the science frame(s)
        input_cfg_dict (dict, optional): 
            Input configuration dictionary. Defaults to None.
        det (str, optional): Detector/mosaic. Defaults to None.
        remove_sci_dir (bool, optional): Remove the science directory if it exists. Defaults to True.
        slitspatnum (str, optional):  
            Value for slitspatnum, e.g. MSCO2:4244
        maskID (str, optional): Mask ID to isolate for QL.  Defaults to None.
        boxcar_radius (float, optional): Boxcar radius for extraction.  
            In units of arcsec.  Defaults to None.
        stack (bool, optional): Stack the science frames.  Defaults to True.

    Returns: 
        tuple: name of pypeit file (str), pypeitFile object (:class:`~pypeit.inputfiles.PypeItFile`)
    """

    # Parse science file info
    folder = folder_name_from_scifiles(sci_files)
    sci_dir = os.path.join(redux_path, folder)
    masters_dir = os.path.join(sci_dir, 'Masters')

    # Science reduction folder
    if os.path.isdir(sci_dir) and remove_sci_dir:
        shutil.rmtree(sci_dir)
    if not os.path.isdir(sci_dir):
        os.makedirs(sci_dir)
        
    # Link to Masters
    if not os.path.isdir(masters_dir):
        os.symlink(master_calib_dir, masters_dir)
        
    # Configure
    user_cfg = ['[rdx]', f'spectrograph = {ps_sci.spectrograph.name}']
    if det is not None:
        user_cfg += [f'detnum = {det}']
    user_cfg += ['quicklook = True']
    user_cfg += ['[baseprocess]', f'master_setup_and_bit = {master_setup_and_bit[0]}_{master_setup_and_bit[1]}']
    full_cfg = configobj.ConfigObj(user_cfg)

    # Add input configs
    if input_cfg_dict is not None:
        full_cfg.merge(configobj.ConfigObj(input_cfg_dict))

    # maskID specified?
    if maskID is not None:
        # Loop on SlitTrace files
        slittrace_files = glob.glob(os.path.join(
            masters_dir, 
            f'MasterSlits_{master_setup_and_bit[0]}_{master_setup_and_bit[1]}_*'))
        detname = None
        for sliittrace_file in slittrace_files:
            slitTrace = SlitTraceSet.from_file(sliittrace_file)
            if maskID in slitTrace.maskdef_id:
                detname = slitTrace.detname
                # Mosaic?
                mosaic = True if detname[0:3] == 'MSC' else False
                det_id = np.where(
                    ps_sci.spectrograph.list_detectors(
                        mosaic=mosaic) == detname)[0][0]
                # Set det
                if mosaic:
                    detnum = [ps_sci.spectrograph.allowed_mosaics[det_id]]
                else:
                    detnum = det_id+1 # 1-based indexing
                # Break
                break
        if detname is None:
            msgs.error('Could not find a SlitTrace file with maskID={}'.format(maskID))

        # Add to config
        maskID_dict = dict(rdx=dict(detnum=detnum,
                                    maskIDs=maskID))
        full_cfg.merge(configobj.ConfigObj(maskID_dict))
            
    # slitspatnum specified?
    if slitspatnum is not None: #'rdx' in full_cfg.keys() and 'slitspatnum' in full_cfg['rdx'].keys():
        ssn_dict = dict(rdx=dict(slitspatnum=slitspatnum))
        full_cfg.merge(configobj.ConfigObj(ssn_dict))

    # Boxcar radius?
    if boxcar_radius is not None:
        # Add to config
        boxcar_lines = ['[reduce]', '[[extraction]]', f'boxcar_radius = {boxcar_radius}']
        full_cfg.merge(configobj.ConfigObj(boxcar_lines))
        

    # Generate PypeIt file
    config_lines = full_cfg.write()

    # Setup, forcing name to match Masters 
    setup = ps_sci.fitstbl.configs.copy()
    key = list(setup.keys())[0]
    setup[f'Setup {master_setup_and_bit[0]}'] = setup[key].copy()
    setup.pop(key)
    # Odds and ends at the finish
    output_cols = ps_sci.fitstbl.set_pypeit_cols(write_bkg_pairs=True,
                                           write_manual=False)
    file_paths = np.unique([os.path.dirname(ff) for ff in sci_files]).tolist()

    # Stack?
    if len(sci_files) > 1 and stack:
        ps_sci.fitstbl['comb_id'] = 1

    # Generate
    pypeitFile = inputfiles.PypeItFile(
        config=config_lines, 
        file_paths=file_paths,
        data_table=ps_sci.fitstbl.table[output_cols],
        setup=setup)

    # Write
    pypeit_file = f'{ps_sci.spectrograph.name}_{master_setup_and_bit[0]}.pypeit' 
    science_pypeit_filename = os.path.join(sci_dir, pypeit_file)
    pypeitFile.write(science_pypeit_filename)

    # Return
    return science_pypeit_filename, pypeitFile



def match_science_to_calibs(ps_sci:pypeitsetup.PypeItSetup, 
                            calib_dir:str):
    """
    Match a given science frame to the set of pre-made calibrations
    in the specified reduction folder. If any exists 
    
    Args:
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
        calib_dir (str): Full path to the calibration directory

    Returns:
        tuple: str, str
            Name of PypeIt file for calibrations
            Full path to Masters

    """
    # Check on one setup
    if len(ps_sci.fitstbl.configs.keys()) > 1:
        msgs.error('Your science files come from more than one setup. Please reduce them separately.')
    # Check against existing calibration PypeIt files
    pypeit_files = glob.glob(os.path.join(
        calib_dir, f'{ps_sci.spectrograph.name}_*', 
        f'{ps_sci.spectrograph.name}_calib_*.pypeit'))
    mtch = []
    for pypeit_file in pypeit_files:
        # Read
        pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)

        # Check for a match
        match = True
        for key in ps_sci.spectrograph.configuration_keys():
            if ps_sci.fitstbl.configs['A'][key] != pypeitFile.setup[key]:
                match = False
        if match:
            mtch.append(pypeit_file)
    # Are we ok?
    if len(mtch) == 0:
        msgs.error("Matched to zero setups.  The calibrations files are stale/wrong/corrupt..")
    elif len(mtch) > 1:
        msgs.warn("Matched to multiple setups.  This should not have happened, but we will take the first.")

    # Master dir
    masters_dir = os.path.join(os.path.dirname(mtch[0]), 'Masters')

    return mtch[0], masters_dir

class QL(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Script to produce quick-look multislit PypeIt reductions', width=width)
        parser.add_argument('spectrograph', type=str,
                            help='A valid spectrograph identifier: {0}'.format(
                                 ', '.join(available_spectrographs)))
        parser.add_argument('--rawfile_list', type=str, 
                            help='File providing raw files to reduce including their path(s)')
        parser.add_argument('--full_rawpath', type=str, 
                            help='Full path to the raw files. Used with --rawfiles or --raw_extension')
        parser.add_argument('--raw_extension', type=str, default='.fits',
                            help='Extension for raw files in full_rawpath.  Only use if --rawfile_list and --rawfiles are not provided')
        parser.add_argument('--rawfiles', type=str, nargs='+',
                            help='space separated list of raw frames e.g. img1.fits img2.fits.  These must exist within --full_rawpath')
        parser.add_argument('--sci_files', type=str, nargs='+',
                            help='space separated list of raw frames to be specified as science exposures (over-rides PypeIt frame typing)')
        parser.add_argument("--redux_path", type=str, default=os.getcwd(),
                            help="Full path to where QL reduction should be run.")
        parser.add_argument("--calib_dir", type=str, 
                            help="Location folders of calibration reductions")
        parser.add_argument("--masters_dir", type=str, 
                            help="Location of PypeIt Master files used for the reduction.")
        parser.add_argument("--calibs_only", default=False, action="store_true",
                            help='Reduce only the calibrations?')
        parser.add_argument("--clobber_calibs", default=False, 
                            action="store_true",
                            help='Clobber existing calibration files?')
        parser.add_argument('--slitspatnum', type=str,
                            help='Reduce the slit(s) as specified by the slitspatnum value(s)')
        parser.add_argument('--maskID', type=int,
                            help='Reduce the slit(s) as specified by the maskID value(s)')
        parser.add_argument('--boxcar_radius', type=float,
                            help='Set the radius for the boxcar extraction in arcseconds')
        parser.add_argument('--det', type=str, help='Detector to reduce. Same format as detnum')
        parser.add_argument('--no_stack', dest='stack', default=True, 
                            action="store_false",
                            help='Do *not* stack multiple science frames')
        return parser


    @staticmethod
    def main(args):

        tstart = time.perf_counter()

        # Ingest Files 
        files = io.grab_rawfiles(
            raw_paths=[args.full_rawpath], 
            file_of_files=args.rawfile_list, 
            list_of_files=args.rawfiles) 

        # Run PypeIt Setup on all the files
        ps = pypeitsetup.PypeItSetup.from_rawfiles(files,
                                                args.spectrograph)
        ps.run(setup_only=True, no_write_sorted=True)

        # Calibrate, if necessary
        calib_dir = args.calib_dir if args.calib_dir is not None else args.redux_path
        if args.masters_dir is None:
            # Generate PypeIt files (and folders)
            calib_pypeit_files = ps.generate_ql_calib_pypeit_files(
                calib_dir, det=args.det, configs='all',
                clobber=args.clobber_calibs)
            # Process them
            for calib_pypeit_file in calib_pypeit_files: 
                redux_path = os.path.dirname(calib_pypeit_file)  # Path to PypeIt file
                # Check for existing masters
                master_files = glob.glob(os.path.join(
                    redux_path, 'Masters', 'Master*'))
                if len(master_files) > 0 and not args.clobber_calibs:
                    msgs.info('Master files already exist.  Skipping calibration.')
                    continue
                # Run
                pypeIt = pypeit.PypeIt(calib_pypeit_file,
                                       redux_path=redux_path, 
                                       calib_only=True)
                calib_dict = pypeIt.calib_all()
        if args.calibs_only:
            msgs.info("Calibrations only requested.  Exiting")
            return

        # Science files                                
        if args.sci_files is not None:
            sci_idx = np.in1d(ps.fitstbl['filename'], args.sci_files)
        else:
            sci_idx = ps.fitstbl.find_frames('science')

        if np.sum(sci_idx) == 0:
            msgs.error('No science frames found in the provided files.  Add at least one or specify using --sci_files.')

        # Generate science setup object
        full_scifiles = [os.path.join(dir_path, sci_file) for dir_path, sci_file in zip(
            ps.fitstbl['directory'][sci_idx], ps.fitstbl['filename'][sci_idx])]
        ps_sci = pypeitsetup.PypeItSetup.from_rawfiles(
                full_scifiles, ps.spectrograph.name)
        ps_sci.run(setup_only=True, no_write_sorted=True)
        # Limit to a single science setup
        if len(ps_sci.fitstbl.configs.keys()) > 1:
            msgs.error('Your science files come from more than one setup. Please reduce them separately.')

        # Masters dir and their setup
        if args.masters_dir is None:
            # Match to calibs
            _, masters_dir = match_science_to_calibs(
                ps_sci, calib_dir)
        else:
            masters_dir = args.masters_dir

        # Parse for setup
        master_files = glob.glob(os.path.join(
            masters_dir, 'Master*'))
        if len(master_files) == 0:
            msgs.error('No Master files found in {:s}'.format(masters_dir))
        master_key, _ = masterframe.grab_key_mdir(
            master_files[0], from_filename=True)
        masters_setup_and_bit =  master_key.split(masterframe.sep1)[0:2]

        # Build the PypeIt file and link to Masters
        sci_pypeit_file, _ = \
                generate_sci_pypeitfile(
                    args.redux_path, 
                    full_scifiles, 
                    masters_dir, 
                    masters_setup_and_bit, 
                    ps_sci, 
                    maskID=args.maskID, 
                    slitspatnum=args.slitspatnum,
                    det=args.det,
                    boxcar_radius=args.boxcar_radius,
                    stack=args.stack)
        
        # Run it
        redux_path = os.path.dirname(sci_pypeit_file)  # Path to PypeIt file
        pypeIt = pypeit.PypeIt(sci_pypeit_file, 
                               reuse_masters=True,
                               redux_path=redux_path) 
        pypeIt.reduce_all()
        pypeIt.build_qa()
        msgs.info(f'Quicklook completed in {utils.get_time_string(time.perf_counter()-tstart)} seconds')
