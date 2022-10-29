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

"""

from email import header
import os
import glob

import numpy as np

import configobj

from pypeit import inputfiles, msgs
from pypeit import pypeitsetup
from pypeit import utils
from pypeit.spectrographs.util import load_spectrograph
from pypeit.slittrace import SlitTraceSet 
from pypeit.scripts import run_pypeit

from IPython import embed

def default_par():

    cfg_default = {}
    cfg_default['rdx'] = dict(ignore_bad_headers=True)
    cfg_default['scienceframe'] = dict(process=dict(mask_cr=False))
    cfg_default['baseprocess'] = dict(use_biasimage=False)
    cfg_default['reduce'] = dict(extraction=dict(skip_optimal=True),
                                 findobj=dict(skip_second_find=True))
    return cfg_default

def generate_calib_pypeit_files(ps, output_path:str,
                   det:str=None,
                   configs:str='all'):
    # Grab setups
    setups, indx = ps.fitstbl.get_configuration_names(
        return_index=True)

    # Restrict on detector -- May remove this
    ps.user_cfg = ['[rdx]', 'spectrograph = {}'.format(ps.spectrograph.name)]
    if det is not None:
        ps.user_cfg += ['detnum = {}'.format(det)]
    # Avoid crash in flat fielding from saturated slits
    ps.user_cfg += ['[calibrations]', '[[flatfield]]', 'saturated_slits = mask']

    # TODO -- Remove the science files!  We want calibs only

    # Write the PypeIt files
    pypeit_files = ps.fitstbl.write_pypeit(output_path=output_path,
                                          cfg_lines=ps.user_cfg,
                                          configs=configs)

    # Rename calibs
    calib_pypeit_files = []
    for pypeit_file, setup in zip(pypeit_files, setups):

        # Rename with _calib
        calib_pypeit_file = pypeit_file.replace('_{}.pypeit'.format(setup),
                                                '_calib_{}.pypeit'.format(setup))
        os.rename(pypeit_file, calib_pypeit_file)
        calib_pypeit_files.append(calib_pypeit_file)

    return calib_pypeit_files

def process_calibs(calib_pypeit_files:list):

    # Loop on setups, rename + run calibs
    for calib_pypeit_file in calib_pypeit_files: 
        # Run me via the script
        redux_path = os.path.dirname(calib_pypeit_file)  # Path to PypeIt file
        run_pargs = run_pypeit.RunPypeIt.parse_args(
            [calib_pypeit_file, 
             '-r={}'.format(redux_path), '-c'])
        run_pypeit.RunPypeIt.main(run_pargs)


def folder_name_from_scifiles(sci_files:list):
    # For now, we return the first filename
    #  without .fits
    return os.path.splitext(os.path.basename(sci_files[0]))[0]

def generate_sci_pypeitfile(calib_pypeit_file:str, 
                            redux_path:str,
                            sci_files:list,
                            ps_sci_list:list,
                            input_cfg_dict:dict={},
                            remove_sci_dir:bool=True,
                            maskID:str=None):
    """
    Generate the PypeIt file for the science frames
    """

    # Parse science file info
    #science_file = os.path.join(pargs.full_rawpath, pargs.science)
    science_pypeit = calib_pypeit_file.replace('calib', 'science')

    folder = folder_name_from_scifiles(sci_files)
    sci_dir = os.path.join(redux_path, folder)
    master_dir = os.path.join(sci_dir, 'Masters')

    # Science reuction folder
    if os.path.isdir(sci_dir) and remove_sci_dir:
        os.system('rm -rf {}'.format(sci_dir))
    if not os.path.isdir(sci_dir):
        os.makedirs(sci_dir)
        
    # Link to Masters
    calib_dir = os.path.dirname(calib_pypeit_file)
    master_calib_dir = os.path.join(calib_dir, 'Masters')
    if not os.path.isdir(master_dir):
        os.symlink(master_calib_dir, master_dir)
        
    # Continuing..
    science_pypeit = os.path.join(sci_dir, os.path.basename(science_pypeit))
    calibPypeItFile = inputfiles.PypeItFile.from_file(calib_pypeit_file)

    # Add science file to data block?
    gd_files = calibPypeItFile.data['frametype'] != 'science'
    for ps_sci, science_file in zip(ps_sci_list, sci_files):
        if science_file not in calibPypeItFile.filenames:
            # NEED TO DEVELOP THIS
            embed(header='125 of ql')
            new_row = {}
            for key in calibPypeItFile.data.keys():
                new_row[key] = ps_sci.fitstbl[key][0]
            new_row['filename'] = science_file
        # Add to good files
        mt = calibPypeItFile.data['filename'] == os.path.basename(science_file)
        gd_files = gd_files | mt

    # Cut down
    cut_data = calibPypeItFile.data[gd_files]

    # Add to configs
    ql_cfg = configobj.ConfigObj(default_par())
    full_cfg = calibPypeItFile.config
    full_cfg.merge(ql_cfg)
    if len(input_cfg_dict) > 0:
        full_cfg.merge(configobj.ConfigObj(input_cfg_dict))


    if maskID is not None:
        # Loop on SlitTrace files
        slittrace_files = glob.glob(os.path.join(
            master_dir, 
            f'MasterSlits_{calibPypeItFile.setup_name}_1_*'))
        detname = None
        for sliittrace_file in slittrace_files:
            slitTrace = SlitTraceSet.from_file(sliittrace_file)
            if maskID in slitTrace.maskdef_id:
                detname = slitTrace.detname
                spectrograph = load_spectrograph(
                    calibPypeItFile.config['rdx']['spectrograph'])
                
                mosaic_id = np.where(spectrograph.list_detectors(mosaic=True) == detname)[0][0]
                det_tuple = spectrograph.allowed_mosaics[mosaic_id]
                break
        if detname is None:
            msgs.error('Could not find a SlitTrace file with maskID={}'.format(maskID))

        # Add to config
        maskID_dict = dict(rdx=dict(detnum=[det_tuple],
                                    maskIDs=maskID))
        full_cfg.merge(configobj.ConfigObj(maskID_dict))
            
    # A touch of voodoo for slitspat_num
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
    pypeitFile = inputfiles.PypeItFile(config=config_lines, 
                                       file_paths=calibPypeItFile.file_paths,
                                       data_table=cut_data,
                                       setup=calibPypeItFile.setup)
    pypeitFile.write(science_pypeit)

    # Return
    return science_pypeit, pypeitFile



def match_science_to_calibs(science_file:str,
                            ps_sci:pypeitsetup.PypeItSetup, 
                            spectrograph,
                            calib_dir:str):
    """
    Match a given science frame to the set of pre-made calibrations
    in the specified reduction folder. If any exists 
    
    Args:
        science_file (str): Full path to the science file
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`): 

    Returns:
        tuple: str, str
            Name of PypeIt file for calibrations
            PypeIt Setup class for the science file
            Name of setup key

    """
    # Check file exists
    if not os.path.isfile(science_file):
        msgs.error("Your science filename {} does not exist. Check your path".format(science_file))


    # Generate the setup dict and yamilfy (yes, this is necessary)
    setup_dict = {}
    for key in spectrograph.configuration_keys():
        setup_dict[key] = ps_sci.fitstbl[key][0]
    setup_dict = utils.yamlify(setup_dict)

    # Check against existing PypeIt files
    pypeit_files = glob.glob(os.path.join(
        calib_dir, f'{spectrograph.name}_*', 
        f'{spectrograph.name}_calib_*.pypeit'))
    mtch = []
    setup_key = None
    for pypeit_file in pypeit_files:
        # Read
        pypeitFile = inputfiles.PypeItFile.from_file(pypeit_file)

        # Check for a match
        match = True
        for key in spectrograph.configuration_keys():
            if setup_dict[key] != pypeitFile.setup[key]:
                match = False
        if match:
            mtch.append(pypeit_file)
            setup_key = pypeitFile.setup_name
    # Are we ok?
    if len(mtch) != 1:
        embed(header='224 of ql')
        msgs.error("Matched to zero or more than one setup.  Inconceivable!")

    return mtch[0], setup_key