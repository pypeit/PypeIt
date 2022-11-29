""" Base routines for Quicklook scripts 

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
"""

import os
import shutil
import glob

import numpy as np

import configobj

from pypeit import inputfiles, msgs
from pypeit import pypeitsetup
from pypeit.slittrace import SlitTraceSet 

from IPython import embed

def is_on(config:configobj.ConfigObj):
    """ Check whether QL is set to "on"

    Args:
        config (configobj.ConfigObj): 
            parameters

    Returns:
        bool: True if QL is on
    """
    if 'rdx' in config.keys() and 'quicklook' in config['rdx'].keys()\
        and config['rdx']['quicklook']:
        return True
    else:
        False

def folder_name_from_scifiles(sci_files:list):
    """ Folder name for output of QL on science file(s)

    Currently, the code takes the name of the first file.
    This may evolve.. 

    Args:
        sci_files (list): List of science files

    Returns:
        str: Folder name
    """
    # For now, we return the first filename
    #  without .fits
    return os.path.splitext(os.path.basename(sci_files[0]))[0]

def generate_sci_pypeitfile(redux_path:str, 
                            sci_files:list, 
                            master_calib_dir:str,
                            master_setup_and_bit:list,
                            ps_sci, 
                            det:str=None,
                            input_cfg_dict:dict=None, 
                            remove_sci_dir:bool=True, 
                            maskID:str=None):
    """
    Generate the PypeIt file for the science frames
    from the calib PypeIt file
    
    Args:
        redux_path (str): Path to the redux folder
        sci_files (list): List of science files (full path)
        master_calib_dir (str): Path to the master calib folder
        master_setup_and_bit (list): 
            Name of the master setup and bit (list of str)
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
        input_cfg_dict (dict, optional): 
            Input configuration dictionary. Defaults to None.
        det (str, optional): Detector/mosaic. Defaults to None.
        remove_sci_dir (bool, optional): Remove the science directory if it exists. Defaults to True.
        maskID (str, optional): Mask ID to isolate for QL.  Defaults to None.

    Returns: 
        tuple: name of pypeit file (str), pypeitFile object (pypeit.inputfiles.PypeItFile)
    """

    # Parse science file info
    folder = folder_name_from_scifiles(sci_files)
    sci_dir = os.path.join(redux_path, folder)
    master_dir = os.path.join(sci_dir, 'Masters')

    # Science reduction folder
    if os.path.isdir(sci_dir) and remove_sci_dir:
        shutil.rmtree(sci_dir)
    if not os.path.isdir(sci_dir):
        os.makedirs(sci_dir)
        
    # Link to Masters
    if not os.path.isdir(master_dir):
        os.symlink(master_calib_dir, master_dir)
        
    # Configure
    user_cfg = ['[rdx]', 'spectrograph = {}'.format(ps_sci.spectrograph.name)]
    if det is not None:
        user_cfg += ['detnum = {}'.format(det)]
    user_cfg += ['quicklook = True']
    full_cfg = configobj.ConfigObj(user_cfg)

    # Add input configs
    if input_cfg_dict is not None:
        full_cfg.merge(configobj.ConfigObj(input_cfg_dict))

    # maskID specified?
    if maskID is not None:
        # Loop on SlitTrace files
        slittrace_files = glob.glob(os.path.join(
            master_dir, 
            f'MasterSlits_{master_setup_and_bit[0]}_{master_setup_and_bit[1]}_*'))
        detname = None
        for sliittrace_file in slittrace_files:
            slitTrace = SlitTraceSet.from_file(sliittrace_file)
            if maskID in slitTrace.maskdef_id:
                detname = slitTrace.detname
                mosaic_id = np.where(ps_sci.spectrograph.list_detectors(mosaic=True) == detname)[0][0]
                det_tuple = ps_sci.spectrograph.allowed_mosaics[mosaic_id]
                break
        if detname is None:
            msgs.error('Could not find a SlitTrace file with maskID={}'.format(maskID))

        # Add to config
        maskID_dict = dict(rdx=dict(detnum=[det_tuple],
                                    maskIDs=maskID))
        full_cfg.merge(configobj.ConfigObj(maskID_dict))
            
    # slitspatnum specified?
    if 'rdx' in full_cfg.keys() and 'slitspatnum' in full_cfg['rdx'].keys():
        # Remove detnum
        for kk, item in enumerate(config_lines):
            if 'detnum' in item:
                config_lines.pop(kk)

        # Add in name, slitspatnum
        ridx = config_lines.index('[rdx]')
        config_lines.insert(ridx+1, '    slitspatnum = {0}'.format(full_cfg['rdx']['slitspatnum']))

        # this is to avoid that the default detnum (which was introduced for mosaic)
        # will be passed to the reduction and crash it
        config_lines.insert(ridx+2, '    detnum = None')

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
    if len(mtch) != 1:
        msgs.error("Matched to zero or more than one setup.  Inconceivable!")

    # Master dir
    master_dir = os.path.join(os.path.dirname(mtch[0]), 'Masters')

    return mtch[0], master_dir