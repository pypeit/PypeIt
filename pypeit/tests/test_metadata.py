from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import glob
import pytest

from pypeit.par.util import parse_pypeit_file
from pypeit.pypeitsetup import PypeItSetup
from pypeit.tests.tstutils import dev_suite_required

@dev_suite_required
def test_lris_red_multi_400():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                          'multi_400_8500_d560', '*.fits.gz'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red',
                 '[calibrations]',
                 '[[pixelflatframe]]',
                 'number = 3',
                 '[[standardframe]]',
                 'number = 0']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    ps.match_to_science(setup_only=True)


@dev_suite_required
def test_lris_red_multi():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red',
                 '[calibrations]',
                 '[[pixelflatframe]]',
                 'number = 3',
                 '[[standardframe]]',
                 'number = 0']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    ps.match_to_science(setup_only=True)


@dev_suite_required
def test_lris_red_multi_run():
    # Perform the setup
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red',
                 '[calibrations]',
                 '[[pixelflatframe]]',
                 'number = 3',
                 '[[standardframe]]',
                 'number = 0']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.run(setup_only=True)

    # Test
    assert len(ps.setup_dict) == 2, 'Should find two setups'
    assert len(ps.fitstbl) == 40, 'Should find 40 files'
    arcs = ps.fitstbl['filename'][ps.fitstbl.find_frames('arc')]
    assert len(arcs) == 2, 'Should find two arcs'
    assert 'r170320_2017.fits.gz' in arcs, \
            'Should have identified r170320_2017.fits.gz as an arc'
    assert 'r170816_0057.fits' in ps.fitstbl['filename'][ps.fitstbl.find_frames('science')], \
            'Should have identified r170816_0057.fits as a science frame'

    # Clean-up
    os.remove('keck_lris_red.lst')
    os.remove('keck_lris_red.setups')
    os.remove('keck_lris_red.sorted')

