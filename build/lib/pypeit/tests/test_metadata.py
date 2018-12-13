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
from pypeit.metadata import PypeItMetaData


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


@dev_suite_required
def test_lris_blue_pypeit_overwrite():

    # JFH The RAW_DATA is at PYPEIT_DEV, but the github DEV suite where
    # the pypeit files are at a different path. To fix this, I've just
    # made pypeit_files directory in Cooked and copied this file over.

    # KBW: I'd prefer that you put symlinks in the github dev suite:
    # e.g.:
    #   cd $PYPEIT_DEV
    #   ln -s /Volumes/GoogleDrive/Team\ Drives/PHYS-GP-Hennawi/PypeIt/PypeIt-development-suite/RAW_DATA  RAW_DATA

    # Read the dev suite pypeit file
    f = os.path.join(os.environ['PYPEIT_DEV'],
                     'Cooked/pypeit_files/keck_lris_blue_long_400_3400_d560.pypeit')
    if not os.path.isfile(f):
        f = os.path.join(os.environ['PYPEIT_DEV'],
                         'pypeit_files/keck_lris_blue_long_400_3400_d560.pypeit')
    assert os.path.isfile(f), 'Could not find pypeit file.'
        
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(f, file_check=False)

    # Change the dev path
    for i in range(len(data_files)):
        path_list = data_files[i].split('/')
        for j,p in enumerate(path_list):
            if p == 'RAW_DATA':
                break
        data_files[i] = os.path.join(os.environ['PYPEIT_DEV'], '/'.join(path_list[j:]))

    # Read the fits table with and without the user data
    fitstbl = PypeItMetaData('keck_lris_blue', file_list=data_files)
    fitstbl_usr = PypeItMetaData('keck_lris_blue', file_list=data_files, usrdata=usrdata)

    assert fitstbl['dispname'][0] == '600/7500', 'Grating name changed in file header'
    assert fitstbl_usr['dispname'][0] == '400/3400', 'Grating name changed in pypeit file'
    assert fitstbl['dispname'][0] != fitstbl_usr['dispname'][0], \
            'Fits header value and input pypeit file value expected to be different.'

