"""
Module to do a full run of PypeIt on good ole Kastb
"""
import os
import sys
import glob
import shutil
from IPython.terminal.embed import embed

from configobj import ConfigObj

import pytest

import numpy as np

import matplotlib
matplotlib.use('agg')  

from pypeit.scripts.parse_calib_id import ParseCalibID
from pypeit.scripts.setup import Setup
from pypeit.scripts.run_pypeit import RunPypeIt
from pypeit.tests.tstutils import data_path
from pypeit import specobjs


def test_run_pypeit():

    # Just get a few files
    testrawdir = data_path('')

    outdir = data_path('REDUX_OUT_TEST')

    # For previously failed tests
    if os.path.isdir(outdir):
        shutil.rmtree(outdir)

    # Run the setup
    sargs = Setup.parse_args(['-r', testrawdir+'b', '-s', 
                              'shane_kast_blue', '-c all', '-o', 
                              '--output_path', outdir])
    Setup.main(sargs)

    # Change to the configuration directory and set the pypeit file
    configdir = os.path.join(outdir, 'shane_kast_blue_A')
    pyp_file = os.path.join(configdir, 'shane_kast_blue_A.pypeit')
    assert os.path.isfile(pyp_file), 'PypeIt file not written.'

    # Try to run with -m and -o
    pargs = RunPypeIt.parse_args([pyp_file, '-o', '-m', '-r', configdir])
    RunPypeIt.main(pargs)

    # TODO: Add some code that will try to open the QA HTML and check that it
    # has the correct PNGs in it.

    # #########################################################
    # Test!!
    # Files exist
    spec1d_file = os.path.join(configdir, 'Science',
                               'spec1d_b27-J1217p3905_KASTb_20150520T045733.560.fits')
    assert os.path.isfile(spec1d_file), 'spec1d file missing'

    # spec1d
    specObjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
    
    # Check RMS
    assert specObjs[0].WAVE_RMS < 0.02  # difference must be less than 0.02 pixels

    # Flexure
    assert abs(-0.03 - specObjs[0].FLEX_SHIFT_TOTAL) < 0.1  # difference must be less than 0.1 pixels

    # Helio
    assert abs(specObjs[0].VEL_CORR - 0.9999251762866389) < 1.0E-10

    # Now re-use those master files
    pargs = RunPypeIt.parse_args([pyp_file, '-o', '-r', configdir])
    RunPypeIt.main(pargs)

    # Clean-up
    shutil.rmtree(outdir)



