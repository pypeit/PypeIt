
import os
import glob
import shutil

from IPython import embed

import numpy as np

import pytest

from pypeit.tests.tstutils import dev_suite_required
from pypeit.pypeitsetup import PypeItSetup
from pypeit.par.util import parse_pypeit_file

@dev_suite_required
def test_deimos():
    # Raw DEIMOS directory
    raw_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos')

    # Get the list of setup directories
    setups = glob.glob(os.path.join(raw_dir, '*'))

    # Set the output path and *remove if* if it already exists
    output_path = os.path.join(os.getcwd(), 'output')
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)

    # Iterate through the setups
    for setup in setups:
 
        # Find the relevant pypeit file constructed by hand.
        by_hand_pypeit = os.path.join(os.getenv('PYPEIT_DEV'), 'pypeit_files',
                                      'keck_deimos_{0}.pypeit'.format(
                                        os.path.split(setup)[1].lower()))

        if not os.path.isfile(by_hand_pypeit):
            # It doesn't exist, so assume there is no by-hand pypeit
            # file to compare to
            continue

        # Run pypeit_setup
        ps = PypeItSetup.from_file_root(setup, 'keck_deimos', output_path=output_path)
        ps.run(setup_only=True)
        # Write the automatically generated pypeit data
        pypeit_files = ps.fitstbl.write_pypeit(output_path, cfg_lines=ps.user_cfg)

        # Read the frame types from both the by-hand and automated
        # pypeit files
        _, _, by_hand_frametypes, _, _ = parse_pypeit_file(by_hand_pypeit, file_check=False)
        _, _, auto_frametypes, _, _ = parse_pypeit_file(pypeit_files[0], file_check=False)

        # For each file in the by-hand list, check that the frame types
        # in the automatically generated pypeit file are identical
        for f in by_hand_frametypes.keys():
            type_list = np.sort(by_hand_frametypes[f].split(','))
            if 'science' in type_list or 'standard' in type_list:
                # Only ensuring that calibrations are correctly typed
                continue
            assert f in auto_frametypes.keys(), \
                'Frame {0} not automatically parsed for setup {1}.'.format(f, setup)
            assert np.array_equal(type_list, np.sort(auto_frametypes[f].split(','))), \
                'Frame types differ for file {0} in setup {1}\n'.format(f, setup) \
                 + '    By-hand types: {0}'.format(by_hand_frametypes[f]) \
                 + '    Automated types: {0}'.format(auto_frametypes[f])

        # Clean up after every setup
        shutil.rmtree(output_path)

