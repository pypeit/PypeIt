#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script examines a set of files and indicates which do and
which do not have sufficient calibs
"""
import argparse
from pypeit import defs

from IPython import embed

def parser(options=None):
    parser = argparse.ArgumentParser(description="Script to check for calibrations [v1]")
    parser.add_argument('root', type=str, default=None, help='File path+root, e.g. /data/Kast/b ')
    parser.add_argument('-s', '--spectrograph', default=None, type=str,
                        help='A valid spectrograph identifier: {0}'.format(
                                ', '.join(defs.pypeit_spectrographs)))
    parser.add_argument('-e', '--extension', default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """

    Args:
        args:

    Returns:
        astropy.table.Table:

    """

    import os
    import numpy as np

    from astropy import table

    from pypeit.pypeitsetup import PypeItSetup
    from pypeit import calibrations
    from pypeit import msgs
    from pypeit.par import PypeItPar

    # Check that the spectrograph is provided if using a file root
    if args.root is not None:
        if args.spectrograph is None:
            raise ValueError('Must provide spectrograph identifier with file root.')
        # Check that input spectrograph is supported
        instruments_served = defs.pypeit_spectrographs
        if args.spectrograph not in instruments_served:
            raise ValueError('Instrument \'{0}\' unknown to PypeIt.\n'.format(args.spectrograph)
                             + '\tOptions are: {0}\n'.format(', '.join(instruments_served))
                             + '\tSelect an available instrument or consult the documentation '
                             + 'on how to add a new instrument.')

    # Initialize PypeItSetup based on the arguments
    ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension)

    # Run the setup
    ps.run(setup_only=True)#, write_bkg_pairs=args.background)
    is_science = ps.fitstbl.find_frames('science')

    msgs.info('Loaded spectrograph {0}'.format(ps.spectrograph.spectrograph))

    # Unique configurations
    setups, indx = ps.fitstbl.get_configuration_names(return_index=True)

    answers = table.Table()
    answers['setups'] = setups
    passes, scifiles, cfgs = [], [], []

    for setup, i in zip(setups, indx):
        # Get the setup lines
        cfg = ps.fitstbl.get_setup(i, config_only=False)
        cfgs.append(cfg)
        if setup == 'None':
            print("There is a setup without science frames.  Skipping...")
            passes.append(False)
            scifiles.append(None)
            continue
        in_cfg = ps.fitstbl['setup'] == setup
        # TODO -- Make the snippet below, which is also in the init of PypeIt a method somewhere
        config_specific_file = None

        msgs.info('=======================================================================')
        msgs.info('Working on setup: {}'.format(setup))
        msgs.info(str(cfg))
        msgs.info('=======================================================================')

        # Grab a science/standard frame
        data_files = [os.path.join(row['directory'], row['filename']) for row in ps.fitstbl[in_cfg]]
        for idx, row in enumerate(ps.fitstbl[in_cfg]):
            if ('science' in row['frametype']) or ('standard' in row['frametype']):
                config_specific_file = data_files[idx]
        if config_specific_file is not None:
            msgs.info(
                'Setting configuration-specific parameters using {0}'.format(os.path.split(config_specific_file)[1]))
        else:
            msgs.warn('No science or standard frame.  Punting..')
            passes.append(False)
            scifiles.append(None)
            continue
        #
        spectrograph_cfg_lines = ps.spectrograph.config_specific_par(config_specific_file).to_config()

        #   - Build the full set, merging with any user-provided
        #     parameters
        par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines)
        # Print science frames
        if np.any(in_cfg & is_science):
            msgs.info("Your science frames are: {}".format(ps.fitstbl['filename'][in_cfg & is_science]))
            scifiles.append(','.join(ps.fitstbl['filename'][in_cfg & is_science]))
        else:
            msgs.warn("This setup has no science frames!")
            scifiles.append('')
        # Check!
        passed = calibrations.check_for_calibs(par, ps.fitstbl, raise_error=False,
                                               cut_cfg=in_cfg)
        if not passed:
            msgs.warn("Setup {} did not pass the calibration check!".format(setup))
        #
        passes.append(passed)

    print('= RESULTS ============================================')

    # Pass/fail
    answers['pass'] = passes

    # Parse the configs
    pcfg = dict(disperser=[], angle=[], dichroic=[], decker=[], slitwid=[], binning=[])
    for cfg in cfgs:
        # None?
        if len(cfg) == 0:
            for key in pcfg.keys():
                pcfg[key].append(None)
            continue
        # Load it up
        key0 = list(cfg.keys())[0]
        subd = cfg[key0]['--'] # for convenience
        pcfg['disperser'].append(subd['disperser']['name'])
        pcfg['angle'].append(subd['disperser']['angle'])
        pcfg['dichroic'].append(subd['dichroic'])
        pcfg['decker'].append(subd['slit']['decker'])
        pcfg['slitwid'].append(subd['slit']['slitwid'])
        pcfg['binning'].append(subd['binning'])

    # Add
    for key in pcfg.keys():
        answers[key] = pcfg[key]

    # Sci files [put this last as it can get large]
    answers['scifiles'] = scifiles

    # Print
    answers.pprint_all()
    print('= ===================================================================================')
    # Return
    return answers, ps
