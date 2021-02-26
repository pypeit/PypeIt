import os
import glob
import shutil
import yaml

from IPython import embed

import pytest

import numpy as np

from pypeit.par.util import parse_pypeit_file
from pypeit.pypeitsetup import PypeItSetup
from pypeit.tests.tstutils import dev_suite_required, data_path
from pypeit.metadata import PypeItMetaData
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts import setup

def test_read_combid():

    # ------------------------------------------------------------------
    # In case of failed tests
    setup_dir = data_path('setup_files')
    if os.path.isdir(setup_dir):
        shutil.rmtree(setup_dir)
    config_dir = data_path('shane_kast_blue_A')
    if os.path.isdir(config_dir):
        shutil.rmtree(config_dir)
    # ------------------------------------------------------------------

    # Generate the pypeit file with the comb_id
    droot = data_path('b')
    pargs = setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c=all', '-b',
                             '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    setup.main(pargs)
    shutil.rmtree(setup_dir)

    pypeit_file = os.path.join(config_dir, 'shane_kast_blue_A.pypeit')
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(pypeit_file)

    # Get the spectrograph
    spectrograph = None
    for l in cfg_lines:
        if 'spectrograph' in l:
            spectrograph = load_spectrograph(l.split(' ')[-1])
            break
    assert spectrograph is not None, 'Did not appropriately read spectrograph'

    # Set the metadata
    pmd = PypeItMetaData(spectrograph, spectrograph.default_pypeit_par(), files=data_files,
                         usrdata=usrdata, strict=False)

    indx = pmd['filename'] == 'b27.fits.gz'
    assert pmd['comb_id'][indx] == [1], 'Incorrect combination group ID'
    assert pmd['comb_id'][~indx] == [-1], 'Incorrect combination group ID'

    shutil.rmtree(config_dir)

@dev_suite_required
def test_lris_red_multi_400():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi_400_8500_d560', '*.fits.gz'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])
    # Test
    assert np.all(ps.fitstbl['setup'] == 'A')


@dev_suite_required
def test_lris_red_multi():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])


@dev_suite_required
def test_lris_red_multi_calib():
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi_400_8500_d560', '*.fits.gz'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.build_fitstbl()
    ps.get_frame_types(flag_unknown=True)
    cfgs = ps.fitstbl.unique_configurations()
    ps.fitstbl.set_configurations(cfgs)
    ps.fitstbl.set_calibration_groups() #global_frames=['bias', 'dark'])

    cfile = data_path('test.calib') 
    ps.fitstbl.write_calib(cfile)
    with open(cfile, 'r') as f:
        calib = yaml.load(f, Loader=yaml.FullLoader)

    assert np.array_equal(list(calib['A'].keys()), ['--', 1]), \
            'Calibrations dictionary read incorrectly.'

    os.remove(cfile)


@dev_suite_required
def test_lris_red_multi_run():
    # Perform the setup
    file_list = glob.glob(os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'keck_lris_red',
                          'multi*', '*.fits*'))
    cfg_lines = ['[rdx]',
                 'spectrograph = keck_lris_red']
    ps = PypeItSetup(file_list, cfg_lines=cfg_lines)
    ps.run(setup_only=True)

    # Test
    #assert len(ps.setup_dict) == 2, 'Should find two setups'
    assert len(ps.fitstbl) >= 40, 'Should find 40+ files'
    arcs = ps.fitstbl['filename'][ps.fitstbl.find_frames('arc')]
    assert len(arcs) >= 2, 'Should find two or more arcs'
    assert 'r170320_2017.fits.gz' in arcs, \
            'Should have identified r170320_2017.fits.gz as an arc'
    assert 'r170816_0057.fits' in ps.fitstbl['filename'][ps.fitstbl.find_frames('science')], \
            'Should have identified r170816_0057.fits as a science frame'

    # Clean-up
    #os.remove('keck_lris_red.lst')
    #os.remove('keck_lris_red.setups')
    os.remove('keck_lris_red.sorted')


@dev_suite_required
def test_lris_blue_pypeit_overwrite():
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
    spectrograph = load_spectrograph('keck_lris_blue')
    par = spectrograph.default_pypeit_par()
    fitstbl = PypeItMetaData(spectrograph, par, files=data_files)
    fitstbl_usr = PypeItMetaData(spectrograph, par, files=data_files, usrdata=usrdata)

    assert fitstbl['target'][0] == 'unknown', 'Grating name changed in file header'
    assert fitstbl_usr['target'][0] == 'test', 'Grating name changed in pypeit file'
    assert fitstbl['target'][0] != fitstbl_usr['target'][0], \
            'Fits header value and input pypeit file value expected to be different.'

