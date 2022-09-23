""" Tests on PypeItFile """

import os
import pytest

from astropy.table import Table

from pypeit.inputfiles import PypeItFile
from pypeit.tests.tstutils import data_path

# Bits needed to generate a PypeIt file

def pieces_of_pypeitfile():
    confdict = {'rdx': {'spectrograph': 'keck_hires'}}

    data = Table()
    data['filename'] = ['b1.fits.gz', 'b27.fits.gz']
    data['frametype'] = ['arc', 'science']
    data['exptime'] = [1., 10.]

    file_paths = [data_path('')]

    setup_dict = {'Setup A': ' '}

    # Return
    return confdict, data, file_paths, setup_dict

# TESTS
def test_instantiate():
    # Test of instantiation
    confdict, data, file_paths, setup_dict = pieces_of_pypeitfile()
    pypeItFile = PypeItFile(confdict, file_paths, 
                                       data, setup_dict)
    # Data files                                    
    filenames = pypeItFile.filenames
    assert 'b1.fits.gz' in filenames[0]

    # Frame types
    frame_type_dict = pypeItFile.frametypes
    assert frame_type_dict['b1.fits.gz'] == 'arc'

    # More tests
    assert pypeItFile.setup_name == 'A'

def test_read_pypeit_file():
    # Read the PypeIt file (backwards compatability)
    pypeItFile = PypeItFile.from_file(
                data_path('example_pypeit_file.pypeit'))
    assert isinstance(pypeItFile.config, dict)

def test_read_backwards_pypeit_file():
    # Read the PypeIt file (backwards compatability)
    pypeItFile = PypeItFile.from_file(
                data_path('example_pypeit_file_backwards.pypeit'))
    assert isinstance(pypeItFile.config, dict)

def test_write_pypeit_file():
    # Test writing a PypeIt file
    outfile = data_path('tmp_file.pypeit')
    if os.path.isfile(outfile):
        os.remove(outfile)

    # Instantiate
    confdict, data, file_paths, setup_dict = pieces_of_pypeitfile()
    pypeItFile = PypeItFile(confdict, file_paths, 
                                       data, setup_dict)
    # Write
    pypeItFile.write(outfile)

    # Let's read it too
    pypeItFile2 = PypeItFile.from_file(outfile)

    # Clean up
    os.remove(outfile)
