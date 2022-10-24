""" Base routines for Quicklook scripts """
import os

from pypeit import inputfiles


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