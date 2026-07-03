"""
Module to run tests on scripts
"""
from pathlib import Path
import os
import shutil

from IPython import embed

import pytest

from astropy.table import Table

from pypeit.scripts.setup import Setup
from pypeit.inputfiles import PypeItFile
from pypeit.inputfiles import RawFiles
from pypeit.tests import tstutils


def expected_file_extensions():
    return ['sorted']


def test_read_list_rawfiles():
    """ Read in a file which is a 
    list of raw data files for setting up
    """
    tst_file = tstutils.data_output_path('test.rawfiles')
    if os.path.isfile(tst_file):
        os.remove(tst_file)

    # Download and move all the required b*fits.gz files into the local package
    # installation
    tstutils.install_shane_kast_blue_raw_data()

    # Bulid
    tbl = Table()
    tbl['filename'] = ['b11.fits.gz', 'b12.fits.gz']
    iRaw = RawFiles(file_paths=[tstutils.data_output_path('')],
                    data_table=tbl)

    # Write
    iRaw.write(tst_file)
    
    # Read
    tst = RawFiles.from_file(tst_file)

    # Test
    assert os.path.basename(tst.filenames[0]) == 'b11.fits.gz'

    # Clean up
    if os.path.isfile(tst_file):
        os.remove(tst_file)


def test_run_setup():
    """ Test the setup script
    """
    # Download and move all the required b*fits.gz files into the local package
    # installation
    tstutils.install_shane_kast_blue_raw_data()

    droot = tstutils.data_output_path('b')
    odir = Path(tstutils.data_output_path('')).absolute() / 'shane_kast_blue_A'
    pargs = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blue', '-c', 'all',
                              '--output_path', f'{odir.parent}'])
    Setup.main(pargs)

    # Fails because name of spectrograph is wrong
    pargs2 = Setup.parse_args(['-r', droot, '-s', 'shane_kast_blu', '-c', 'all',
                               '--output_path', f'{odir.parent}'])
    with pytest.raises(ValueError):
        Setup.main(pargs2)
    
    pypeit_file = tstutils.data_output_path('shane_kast_blue_A/shane_kast_blue_A.pypeit')
    pypeItFile = PypeItFile.from_file(pypeit_file)

    # Test
    assert len(pypeItFile.filenames) == 9
    assert sorted(pypeItFile.frametypes['b1.fits.gz'].split(',')) == ['arc', 'tilt']
    assert pypeItFile.setup_name == 'A'

    # Cleanup
    shutil.rmtree(tstutils.data_output_path('shane_kast_blue_A'))


