# Module to run tests on PypitSetup class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from astropy.table import Table

from pypit import pypitsetup

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def settings_kludge():
    # We be replaced with the settings Refactor
    from pypit import arparse as settings
    settings.dummy_settings()
    settings.argflag['run']['spectrograph'] = 'shane_kast_blue'
    settings.argflag['reduce']['masters']['setup'] = 'C_01_aa'
    #
    # Load default spectrograph settings
    spect = settings.get_spect_class(('ARMLSD', 'shane_kast_blue', 'pypit'))  # '.'.join(redname.split('.')[:-1])))
    lines = spect.load_file(base=True)  # Base spectrograph settings
    spect.set_paramlist(lines)
    lines = spect.load_file()  # Instrument specific
    spect.set_paramlist(lines)
    return settings.argflag, settings.spect

def test_init():
    if skip_test:
        assert True
        return
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    assert len(setupc.steps) == 0
    assert setupc.nfiles == 0


def test_build_fitstbl():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # fitstlb
    fitstbl = setupc.build_fitstbl(files)
    assert isinstance(fitstbl, Table)
    assert setupc.nfiles == 27

    # I/O
    setupc.write_fitstbl(data_path('fitstbl.fits'))
    tmp = setupc.load_fitstbl(data_path('fitstbl.fits'))
    assert len(tmp) == 27


def test_image_type():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # fitstlb
    fitstbl = setupc.build_fitstbl(files)

    # Type
    filetypes = setupc.type_data()
    assert np.sum(filetypes['arc']) == 1
    assert 'arc' in setupc.fitstbl.keys()
    assert np.sum(filetypes['pixelflat'] & filetypes['trace']) == 12


def test_match():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # fitstlb
    _ = setupc.build_fitstbl(files)

    # Type
    _ = setupc.type_data()

    # Match to science
    fitstbl = setupc.match_to_science()
    assert fitstbl['sci_ID'][fitstbl['science']].tolist() == [1,2,4]


def test_run():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # Run
    code, fitstbl, setup_dict = setupc.run(files)
    # Test
    assert code == 'run'
    assert isinstance(fitstbl, Table)
    assert isinstance(setup_dict, dict)


def test_run_calcheck():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    settings_argflag['run']['calcheck'] = True

    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # Run
    code, fitstbl, setup_dict = setupc.run(files)
    # Test
    assert code == 'calcheck'


def test_run_setup():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.getenv('PYPIT_DEV') + 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b'
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge()
    settings_argflag['run']['setup'] = True

    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # Run
    code, fitstbl, setup_dict = setupc.run(files)
    # Test
    assert code == 'setup'
