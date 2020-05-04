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
    # TODO: Add argument that specifies the log file
    parser = argparse.ArgumentParser(description="Script to setup a PypeIt run [v3]")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-r', '--root', type=str, default=None,
                       help='File path+root, e.g. /data/Kast/b ')
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
    if args.root is not None:
        ps = PypeItSetup.from_file_root(args.root, args.spectrograph, extension=args.extension)
    else:
        # Should never reach here
        raise IOError('Need to set -r !!')

    # Run the setup
    ps.run(setup_only=True)#, write_bkg_pairs=args.background)
    is_science = ps.fitstbl.find_frames('science')

    msgs.info('Loaded spectrograph {0}'.format(ps.spectrograph.spectrograph))

    # Unique configurations
    setups, indx = ps.fitstbl.get_configuration_names(return_index=True)

    answers = table.Table()
    answers['setups'] = setups
    passes, scifiles = [], []

    for setup, i in zip(setups, indx):
        msgs.info('=======================================================================')
        msgs.info('Working on setup: {}'.format(setup))
        msgs.info('=======================================================================')
        # Get the setup lines
        cfg = ps.fitstbl.get_setup(i, config_only=False)
        in_cfg = ps.fitstbl['setup'] == setup
        # TODO -- Make the snippet below, which is also in the init of PypeIt a method somehwere
        config_specific_file = None

        data_files = [os.path.join(row['directory'], row['filename']) for row in ps.fitstbl[in_cfg]]
        for idx, row in enumerate(ps.fitstbl[in_cfg]):
            if ('science' in row['frametype']) or ('standard' in row['frametype']):
                config_specific_file = data_files[idx]
        # search for arcs, trace if no scistd was there
        if config_specific_file is None:
            for idx, row in enumerate(ps.fitstbl[in_cfg]):
                if ('arc' in row['frametype']) or ('trace' in row['frametype']):
                    config_specific_file = data_files[idx]
        if config_specific_file is not None:
            msgs.info(
                'Setting configuration-specific parameters using {0}'.format(os.path.split(config_specific_file)[1]))
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

    msgs.info('========================== ALL DONE =======================================')
    print('========================== RESULTS =======================================')
    print('=================================================================')
    #
    answers['passfail'] = passes
    answers['scifiles'] = scifiles
    # Print
    print(answers)
    print('=================================================================')
    # Return
    return answers
