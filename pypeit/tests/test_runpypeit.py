"""
Module to run tests on scripts
"""
import os
import sys
import glob
import shutil

from configobj import ConfigObj

import pytest

import numpy as np

import matplotlib
matplotlib.use('agg')  # For Travis

from astropy.io import fits

from pypeit.scripts import setup
from pypeit.scripts import run_pypeit
from pypeit.tests.tstutils import dev_suite_required
from pypeit import specobjs


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@dev_suite_required
def test_run_pypeit_calib_only():
    # Get the directories
    rawdir = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'shane_kast_blue', '600_4310_d55')
    assert os.path.isdir(rawdir), 'Incorrect raw directory'

    # File list
    all_files = {
        'arcs': ['b1.fits.gz'],
        'flats': ['b11.fits.gz', 'b12.fits.gz', 'b13.fits.gz'],
        'bias': ['b21.fits.gz', 'b22.fits.gz', 'b23.fits.gz'],
    }
    all_masters = ['MasterArc_A_1_01.fits', 'MasterTiltimg_A_1_01.fits',
                   'MasterBias_A_1_01.fits',
                   'MasterTilts_A_1_01.fits', 'MasterEdges_A_1_01.fits.gz',
                   'MasterFlat_A_1_01.fits',
                   'MasterWaveCalib_A_1_01.fits']

    # Just get a few files
    for ss, sub_files, masters in zip(range(3),
            [['arcs', 'flats', 'bias'], ['arcs', 'bias'], ['flats', 'bias']],
            [all_masters, ['MasterArc_A_1_01.fits', 'MasterTiltimg_A_1_01.fits'],
             ['MasterEdges_A_1_01.fits.gz']]):
        # Grab the subset
        files = []
        for sub_file in sub_files:
            files += all_files[sub_file]
        #
        testrawdir = os.path.join(rawdir, 'TEST')
        if os.path.isdir(testrawdir):
            shutil.rmtree(testrawdir)
        os.makedirs(testrawdir)
        for f in files:
            shutil.copy(os.path.join(rawdir, f), os.path.join(testrawdir, f))

        outdir = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT_TEST')

        # For previously failed tests
        if os.path.isdir(outdir):
            shutil.rmtree(outdir)

        # Run the setup
        sargs = setup.parse_args(['-r', testrawdir, '-s', 'shane_kast_blue', '-c all', '-o',
                                  '--output_path', outdir])
        setup.main(sargs)

        # Change to the configuration directory and set the pypeit file
        configdir = os.path.join(outdir, 'shane_kast_blue_A')
        pyp_file = os.path.join(configdir, 'shane_kast_blue_A.pypeit')
        assert os.path.isfile(pyp_file), 'PypeIt file not written.'

        # Perform the calib-only reduction
        pargs = run_pypeit.parse_args([pyp_file, '-c', '-r', configdir])
        run_pypeit.main(pargs)

        # Test!
        for master_file in masters:
            assert os.path.isfile(os.path.join(configdir, 'Masters', master_file)
                                  ), 'Master File {:s} missing!'.format(master_file)

        # Clean-up
        shutil.rmtree(outdir)
        shutil.rmtree(testrawdir)


@dev_suite_required
def test_run_pypeit():
    # Get the directories
    rawdir = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'shane_kast_blue', '600_4310_d55')
    assert os.path.isdir(rawdir), 'Incorrect raw directory'

    # Just get a few files
    testrawdir = os.path.join(rawdir, 'TEST')
    if os.path.isdir(testrawdir):
        shutil.rmtree(testrawdir)
    os.makedirs(testrawdir)
    files = [ 'b21.fits.gz', 'b22.fits.gz', 'b23.fits.gz', 'b27.fits.gz', 'b1.fits.gz',
              'b11.fits.gz', 'b12.fits.gz', 'b13.fits.gz' ]
    for f in files:
        shutil.copy(os.path.join(rawdir, f), os.path.join(testrawdir, f))

    outdir = os.path.join(os.getenv('PYPEIT_DEV'), 'REDUX_OUT_TEST')

    # For previously failed tests
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Run the setup
    sargs = setup.parse_args(['-r', testrawdir, '-s', 'shane_kast_blue', '-c all', '-o',
                              '--output_path', outdir])
    setup.main(sargs)

    # Change to the configuration directory and set the pypeit file
    configdir = os.path.join(outdir, 'shane_kast_blue_A')
    pyp_file = os.path.join(configdir, 'shane_kast_blue_A.pypeit')
    assert os.path.isfile(pyp_file), 'PypeIt file not written.'

    # Try to run with -m and -o
    pargs = run_pypeit.parse_args([pyp_file, '-o', '-m', '-r', configdir])
    run_pypeit.main(pargs)

    # #########################################################
    # Test!!
    # Files exist
    assert os.path.isfile(os.path.join(configdir, 'Science', 'spec2d_b27-J1217p3905_KASTb_2015May20T045733.560.fits'))

    # spec1d
    spec1d_file = os.path.join(configdir, 'Science', 'spec1d_b27-J1217p3905_KASTb_2015May20T045733.560.fits')
    assert os.path.isfile(spec1d_file)
    specObjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)

    # Flexure
    assert abs(-0.03 - specObjs[0].FLEX_SHIFT_TOTAL) < 0.1  # difference must be less than 0.1 pixels

    # Helio
    assert abs(specObjs[0].VEL_CORR - 0.9999261685542624) < 1.0E-10

    # Now re-use those master files
    pargs = run_pypeit.parse_args([pyp_file, '-o', '-r', configdir])
    run_pypeit.main(pargs)

    # Clean-up
    shutil.rmtree(outdir)
    shutil.rmtree(testrawdir)

