"""
Module to run tests on scripts
"""
import os
import glob
import shutil

from IPython import embed

import numpy as np

import pytest

from pypeit.scripts.setup import Setup
from pypeit.inputfiles import PypeItFile
from pypeit.tests.tstutils import data_path


def expected_file_extensions():
    return ['sorted']


def test_run_setup():
    """ Test the setup script
    """
    # Remove .setup if needed
    sfiles = glob.glob('*.setups')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c=all',
                              '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    Setup.main(pargs)

    #setup_file = glob.glob(data_path('setup_files/shane_kast_blue*.setups'))[0]
    ## Load
    #with open(setup_file, 'r') as infile:
    #    setup_dict = yaml.load(infile)
    ## Test
    #assert '01' in setup_dict['A'].keys()
    #assert setup_dict['A']['--']['disperser']['name'] == '600/4310'
    # Failures
    pargs2 = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blu', '-c=all',
                               '--extension=fits.gz', '--output_path={:s}'.format(data_path(''))])
    with pytest.raises(ValueError):
        Setup.main(pargs2)
    
    # Cleanup
    shutil.rmtree(data_path('setup_files'))


def test_setup_made_pypeit_file():
    """ Test the .pypeit file(s) made by pypeit_setup

    This test depends on the one above
    """
    pypeit_file = data_path('shane_kast_blue_A/shane_kast_blue_A.pypeit')
    pypeItFile = PypeItFile.from_file(pypeit_file)

    # Test
    assert len(pypeItFile.filenames) == 8
    assert sorted(pypeItFile.frametypes['b1.fits.gz'].split(',')) == ['arc', 'tilt']
    assert pypeItFile.setup_name == 'A'

    # Cleanup
    shutil.rmtree(data_path('shane_kast_blue_A'))
