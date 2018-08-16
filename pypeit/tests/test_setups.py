# Module to run tests on scripts
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest
import glob

import yaml

from pypeit import msgs
from pypeit.par.util import parse_pypeit_file
from pypeit.scripts import setup

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_run_setup():
    """ Test the setup script
    """
    # Remove .setup if needed
    sfiles = glob.glob('*.setups')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = setup.parser([droot, 'shane_kast_blue', '-c',
                          '--extension=fits.gz', '--redux_path={:s}'.format(data_path(''))])
    setup.main(pargs)

    setup_file = glob.glob(data_path('setup_files/shane_kast_blue*.setups'))[0]
    # Load
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Test
    assert '01' in setup_dict['A'].keys()
    assert setup_dict['A']['--']['disperser']['name'] == '600/4310'
    # Failures
    pargs2 = setup.parser([droot, 'shane_kast_blu', '-d', '-c',
                              '--extension=fits.gz', '--redux_path={:s}'.format(data_path(''))])
    with pytest.raises(ValueError):
        setup.main(pargs2)

def test_setup_made_pypeit_file():
    """ Test the .pypeit file(s) made by pypeit_setup
    This test depends on the one above
    """
    pypeit_file = data_path('shane_kast_blue_setup_A/shane_kast_blue_setup_A.pypeit')
    cfg_lines, data_files, frametype, setups = parse_pypeit_file(pypeit_file)
    print(setups)
    # Test
    assert len(data_files) == 2
    assert frametype['b1.fits.gz'] == 'arc'
    assert setups[0] == 'A'

'''
def test_setup_pfile():
    """ This won't run because of msgs issue..
    """
    from pypit.scripts import setup
    # Try from .pypit file
    sfiles = glob.glob('*.setup')
    pypit_file = sfiles[0].replace('.setup', '.pypit')
    for sfile in sfiles:
        os.remove(sfile)
    pargs2 = setup.parser([pypit_file, 'shane_kast_blue', '--pypit_file', '-d'])
    setup.main(pargs2)
    pytest.set_trace()
'''
