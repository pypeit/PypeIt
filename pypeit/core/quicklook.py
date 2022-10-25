""" Base routines for Quicklook scripts """
import os
import glob

from pypeit import inputfiles, msgs
from pypeit import pypeitsetup
from pypeit import utils

from IPython import embed


def generate_calib_pypeit_files(ps, output_path:str,
                   det:str='1',
                   configs:str='all'):
    # Grab setups
    setups, indx = ps.fitstbl.get_configuration_names(
        return_index=True)

    # Restrict on detector -- May remove this
    ps.user_cfg = ['[rdx]', 'spectrograph = {}'.format(ps.spectrograph.name)]
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

def generate_science_pypeitfiles(calib_pypeit_file, ps_sci):
    """
    Process a science frame

    Args:
        pargs (argparse.ArgumentParser):
            Command line arguments
        calib_pypeit_file (str):
        ps_sci (:class:`pypeit.pypeitsetup.PypeItSetup`):
    """
    # Parse science file info
    science_file = os.path.join(pargs.full_rawpath, pargs.science)
    science_pypeit = calib_pypeit_file.replace('calib', 'science')

    calibPypeItFile = inputfiles.PypeItFile.from_file(calib_pypeit_file)

    # Add science file to data block?
    if science_file not in calibPypeItFile.filenames:
        new_row = {}
        for key in calibPypeItFile.data.keys():
            new_row[key] = ps_sci.fitstbl[key][0]
        new_row['filename'] = pargs.science

    # Generate data block
    # Remove any extraneous science files in the folder
    gd_files = (calibPypeItFile.data['filename'] == os.path.basename(science_file)) | (
        calibPypeItFile.data['frametype'] != 'science')
    cut_data = calibPypeItFile.data[gd_files]

    # Add to configs
    config_lines = calibPypeItFile.cfg_lines
    if pargs.slit_spat is not None:
        # Remove detnum
        for kk, item in enumerate(config_lines):
            if 'detnum' in item:
                config_lines.pop(kk)

        # Add in name, slitspatnum
        ridx = config_lines.index('[rdx]')
        config_lines.insert(ridx+1, '    slitspatnum = {0}'.format(pargs.slit_spat))

        # this is to avoid that the default detnum (which was introduced for mosaic)
        # will be passed to the reduction and crash it
        config_lines.insert(ridx+2, '    detnum = None')
    else:
        raise NotImplementedError('NOT READY:  118 of ql_deimos')

    # Generate PypeIt file
    pypeitFile = inputfiles.PypeItFile(config=config_lines, 
                                       file_paths=calibPypeItFile.file_paths,
                                       data_table=cut_data,
                                       setup=calibPypeItFile.setup)
    pypeitFile.write(science_pypeit)

    # Run me!
    redux_path = os.path.dirname(science_pypeit)  # Path to PypeIt file
    run_pargs = run_pypeit.RunPypeIt.parse_args([science_pypeit,
                                   '-r={}'.format(redux_path),
                                   ])
    run_pypeit.RunPypeIt.main(run_pargs)



def generate_sci_pypeitfile(calib_pypeit_file:str, 
                            sci_files:list,
                            ps_sci_list:list,
                            cfg_dict:dict=None):
    """
    Generate the PypeIt file for the science frames
    """

    # Parse science file info
    #science_file = os.path.join(pargs.full_rawpath, pargs.science)
    science_pypeit = calib_pypeit_file.replace('calib', 'science')

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
    config_lines = calibPypeItFile.cfg_lines

    # TODO -- GENERALIZE THIS
    if cfg_dict is not None:
        if pargs.slit_spat is not None:
            # Remove detnum
            for kk, item in enumerate(config_lines):
                if 'detnum' in item:
                    config_lines.pop(kk)

            # Add in name, slitspatnum
            ridx = config_lines.index('[rdx]')
            config_lines.insert(ridx+1, '    slitspatnum = {0}'.format(pargs.slit_spat))

            # this is to avoid that the default detnum (which was introduced for mosaic)
            # will be passed to the reduction and crash it
            config_lines.insert(ridx+2, '    detnum = None')
        else:
            raise NotImplementedError('NOT READY:  118 of ql_deimos')

    # Generate PypeIt file
    pypeitFile = inputfiles.PypeItFile(config=config_lines, 
                                       file_paths=calibPypeItFile.file_paths,
                                       data_table=cut_data,
                                       setup=calibPypeItFile.setup)
    pypeitFile.write(science_pypeit)

    # Return
    return science_pypeit, pypeitFile


def run_on_science(pargs, calib_pypeit_file, ps_sci):
    # Run me!
    redux_path = os.path.dirname(science_pypeit)  # Path to PypeIt file
    run_pargs = run_pypeit.RunPypeIt.parse_args([science_pypeit,
                                   '-r={}'.format(redux_path),
                                   ])
    run_pypeit.RunPypeIt.main(run_pargs)


def match_science_to_calibs(science_file:str, 
                            spectrograph,
                            redux_path:str):
    """
    Match a given science frame to the set of pre-made calibrations
    in the specified reduction folder. If any exists 
    
    Args:

    Returns:
        tuple: str, :class:`pypeit.pypeitsetup.PypeItSetup`, str
            Name of PypeIt file for calibrations
            PypeIt Setup class for the science file
            Name of setup key

    """
    # Check file exists
    if not os.path.isfile(science_file):
        msgs.error("Your science filename {} does not exist. Check your path".format(science_file))

    # Run setup on the single science file
    ps_sci = pypeitsetup.PypeItSetup.from_file_root(science_file, 
                                        spectrograph.name, 
                                        extension='')
    ps_sci.run(setup_only=True, no_write_sorted=True)

    # Generate the setup dict and yamilfy (yes, this is necessary)
    setup_dict = {}
    for key in spectrograph.configuration_keys():
        setup_dict[key] = ps_sci.fitstbl[key][0]
    setup_dict = utils.yamlify(setup_dict)

    # Check against existing PypeIt files
    pypeit_files = glob.glob(os.path.join(
        redux_path, f'{spectrograph.name}_*', 
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
        msgs.error("Matched to zero or more than one setup.  Inconceivable!")

    return mtch[0], ps_sci, setup_key