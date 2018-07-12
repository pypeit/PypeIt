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
from pypit.par import pypitpar

from pypit.tests import tstutils

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def set_pars():
    reduce_par = pypitpar.ReducePar()
    run_par = pypitpar.RunPar()
    return run_par, reduce_par

def get_files():
    # Check for files
    file_root = os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    return files

# INIT
spectrograph = tstutils.load_kast_blue_masters(get_spectrograph=True)[0]

def test_init():
    if skip_test:
        assert True
        return
    # Init
    run_par, reduce_par = set_pars()
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
    assert len(setupc.steps) == 0
    assert setupc.nfiles == 0


def test_build_fitstbl():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA/Shane_Kast_blue/600_4310_d55/b')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Init
    run_par, reduce_par = set_pars()
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
    #
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
    files = get_files()
    # Init
    run_par, reduce_par = set_pars()
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
    fitstbl = setupc.build_fitstbl(files)
    # Type
    filetypes = setupc.type_data(flag_unknown=True)
    assert np.sum(filetypes['arc']) == 1
    assert np.sum(filetypes['unknown']) == 0
    assert 'arc' in setupc.fitstbl.keys()
    assert np.sum(filetypes['pixelflat'] & filetypes['trace']) == 12

def test_match():
    if skip_test:
        assert True
        return
    # Check for files
    files = get_files()
    # Init
    run_par, reduce_par = set_pars()
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
    _ = setupc.build_fitstbl(files)
    _ = setupc.type_data(flag_unknown=True)

    # Match to science
    fitstbl = setupc.match_to_science()
    assert fitstbl['sci_ID'][fitstbl['science']].tolist() == [1,2,4]

'''
def test_match_ABBA():
    if skip_test:
        assert True
        return
    # Check for files
    file_root = os.path.join(os.getenv('PYPIT_DEV'), 'RAW_DATA/Keck_NIRSPEC/NIRSPEC-1/NS')
    files = glob.glob(file_root+'*')
    assert len(files) > 0
    # Settings
    settings_argflag, settings_spect = settings_kludge('keck_nirspec')
    settings_spect['bias']['number'] = 0
    settings_spect['standard']['number'] = 0
    # Init
    setupc = pypitsetup.PypitSetup(settings_argflag, settings_spect)
    # fitstbl
    _ = setupc.build_fitstbl(files)

    # Type
    _ = setupc.type_data()

    # Match to science
    fitstbl = setupc.match_to_science()
    fitstbl = setupc.match_ABBA()

    assert fitstbl['AB_frame'][-1] == 'NS.20160414.55235.fits.gz'
'''

def test_run():
    if skip_test:
        assert True
        return
    # Check for files
    files = get_files()
    # Init
    run_par, reduce_par = set_pars()
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
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
    files = get_files()
    # Init
    run_par, reduce_par = set_pars()
    run_par['calcheck'] = True
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)
    # Run
    code, fitstbl, setup_dict = setupc.run(files)
    # Test
    assert code == 'calcheck'


def test_run_setup():
    if skip_test:
        assert True
        return
    files = get_files()
    # Init
    run_par, reduce_par = set_pars()
    run_par['setup'] = True
    setupc = pypitsetup.PypitSetup(spectrograph, run_par, reduce_par)

    # Run
    code, fitstbl, setup_dict = setupc.run(files)
    # Test
    assert code == 'setup'
