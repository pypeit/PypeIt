""" Tests on PypeItFile """

import os
import pytest

from astropy.table import Table

from pypeit import pypeitfile
from pypeit.tests.tstutils import data_path

# Bits needed to generate a PypeIt file
confdict = {'rdx': {'spectrograph': 'keck_hires'}}

data = Table()
data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
data['frametype'] = ['arc', 'science']
data['exptime'] = [1., 10.]

file_paths = [data_path('')]

setup_dict = {'Setup A': ' '}

# TESTS
def test_instantiate():
    # Test of instantiation
    pypeItFile = pypeitfile.PypeItFile(confdict, file_paths, 
                                       data, setup_dict)
    # Data files                                    
    data_files = pypeItFile.data_files
    assert 'b1.fits.gz' in data_files[0]

    # Frame types
    frame_type_dict = pypeItFile.frametypes
    assert frame_type_dict['b1.fits.gz'] == 'arc'

    # More tests
    assert pypeItFile.setup_name == 'A'


def test_read_pypeit_file():
    # Read the PypeIt file
    pypeItFile = pypeitfile.PypeItFile.from_file(
                data_path('example_pypeit_file.pypeit'))
    assert isinstance(pypeItFile.config, dict)

def test_write_pypeit_file():
    # Test writing a PypeIt file
    outfile = data_path('tmp_pypeit.file')
    if os.path.isfile(outfile):
        os.remove(outfile)

    # Instantiate
    pypeItFile = pypeitfile.PypeItFile(confdict, file_paths, 
                                       data, setup_dict)
    # Write
    pypeItFile.write(outfile)

    # Clean up
    os.remove(outfile)
