"""
Module to run tests on PypeItSetup class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from astropy.table import Table


from pypeit import pypeitsetup
from pypeit.par import pypeitpar
from pypeit.metadata import PypeItMetaData

from pypeit.tests.tstutils import dev_suite_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def get_files():
    # Check for files
    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    return files


@dev_suite_required
def test_init():
    # Init
    files = get_files()
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    assert len(setupc.steps) == 0
    assert setupc.nfiles == 0


@dev_suite_required
def test_build_fitstbl():
    # Check for files
    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/shane_kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    #
    fitstbl = setupc.build_fitstbl(files)
    assert isinstance(fitstbl, Table)
    assert setupc.nfiles == 26

    # I/O
    setupc.write_metadata(ofile=data_path('fitstbl.fits'))
    tmp = setupc.load_metadata(data_path('fitstbl.fits'))
    assert len(tmp) == 26

    # Cleanup
    os.remove(data_path('fitstbl.fits'))


@dev_suite_required
def test_image_type():
    # Check for files
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    fitstbl = setupc.build_fitstbl(files)
    # Type
    setupc.get_frame_types(flag_unknown=True)
    assert 'framebit' in setupc.fitstbl.keys()
    assert np.sum(setupc.fitstbl.find_frames('arc')) == 1
    assert np.sum(setupc.fitstbl.find_frames('None')) == 0

    assert np.sum(setupc.fitstbl.find_frames('pixelflat')
                    & setupc.fitstbl.find_frames('trace')) == 12

@dev_suite_required
def test_type():
    # Check for files
    files = get_files()
    # Init
    cfg_lines = ['[rdx]',
                 'spectrograph = shane_kast_blue',
                 '[calibrations]',
                 '[[standardframe]]',
                 'exprng = None,60',
                 '[scienceframe]',
                 'exprng = 60,None']
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue',
                                     cfg_lines=cfg_lines, path=data_path(''))
    setupc.build_fitstbl(files)
    setupc.get_frame_types(flag_unknown=True)
    assert np.sum(setupc.fitstbl.find_frames('science')) == 2

#def test_match_ABBA():
#    if skip_test:
#        assert True
#        return
#    # Check for files
#    file_root = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA/Keck_NIRSPEC/NIRSPEC-1/NS')
#    files = glob.glob(file_root+'*')
#    assert len(files) > 0
#    # Settings
#    settings_argflag, settings_spect = settings_kludge('keck_nirspec')
#    settings_spect['bias']['number'] = 0
#    settings_spect['standard']['number'] = 0
#    # Init
#    setupc = pypeitsetup.PypeItSetup(settings_argflag, settings_spect)
#    # fitstbl
#    _ = setupc.build_fitstbl(files)
#
#    # Type
#    _ = setupc.type_data()
#
#    # Match to science
#    fitstbl = setupc.match_to_science()
#    fitstbl = setupc.match_ABBA()
#
#    assert fitstbl['AB_frame'][-1] == 'NS.20160414.55235.fits.gz'

@dev_suite_required
def test_run():
    # Check for files
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    # Run
    par, spectrograph, fitstbl = setupc.run(sort_dir=data_path(''))
    # Test
    assert isinstance(par, pypeitpar.PypeItPar)
    assert isinstance(fitstbl, PypeItMetaData)
    #assert isinstance(setup_dict, dict)

    # Cleanup
    os.remove(data_path('shane_kast_blue.calib'))

@dev_suite_required
def test_run_calcheck():
    # Check for files
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    # Run
    par, spectrograph, fitstbl = setupc.run(calibration_check=True, sort_dir=data_path(''))
    # Test
    assert isinstance(par, pypeitpar.PypeItPar)

    # Cleanup
    os.remove(data_path('shane_kast_blue.calib'))

@dev_suite_required
def test_run_setup():
    files = get_files()
    # Init
    setupc = pypeitsetup.PypeItSetup(files, spectrograph_name='shane_kast_blue', path=data_path(''))
    # Run
    par, spectrograph, fitstbl = setupc.run(setup_only=True, sort_dir=data_path(''))
    # Test
    assert par is None

    # Cleanup
    os.remove(data_path('shane_kast_blue.sorted'))
