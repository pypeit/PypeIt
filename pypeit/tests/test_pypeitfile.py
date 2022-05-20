""" Tests on PypeItFile """

import os
import pytest

from astropy.table import Table

from pypeit import pypeitfile
from pypeit.tests.tstutils import data_path

import configobj

def test_read_pypeit_file():
    # Read the PypeIt file
    pypeItFile = pypeitfile.PypeItFile.from_file(
                data_path('example_pypeit_file.pypeit'))
    assert isinstance(pypeItFile.config, dict)

def test_write_pypeit_file():
    outfile = data_path('tmp_pypeit.file')
    if os.path.isfile(outfile):
        os.remove(outfile)

    data = Table()
    data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
    data['frametype'] = ['arc', 'science']
    data['exptime'] = [1., 10.]
    file_paths = [data_path('')]
    confdict = {'rdx': {'spectrograph': 'keck_hires'}}
    setup_dict = {'Setup A': ' '}


    # Instantiate
    pypeItFile = pypeitfile.PypeItFile(confdict, file_paths, 
                                       data, setup_dict)

    # Write
    pypeItFile.write(outfile)

    # Clean up
    os.remove(outfile)
