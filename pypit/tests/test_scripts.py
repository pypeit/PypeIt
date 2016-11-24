# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest
import glob

from pypit import pyputils

msgs = pyputils.get_dummy_logger()

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


"""
def test_arcid_plot():
    json_file = data_path('LRISb_600_WaveCalib_01.json')
    pargs = arcid_plot.parser([json_file, 'LRISb', 'tmp.pdf'])
    # Run
    arcid_plot.main(pargs)
"""


'''
def test_show_1dspec():
    #from PyQt4 import QtGui
    #app = QtGui.QApplication(sys.argv)
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    pargs = show_1dspec.parser([spec_file, '--list'])
    # Run
    show_1dspec.main(pargs, unit_test=True)
'''


def test_setup():
    """ Test the setup script
    """
    from pypit.scripts import setup
    import yaml
    # Remove .setup if needed
    sfiles = glob.glob('*.setup')
    for sfile in sfiles:
        os.remove(sfile)
    #
    droot = data_path('b')
    pargs = setup.parser([droot, 'kast_blue', '--extension=fits.gz'])
    setup.main(pargs)
    setup_file = glob.glob('kast_blue*.setup')[0]
    # Load
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Test
    assert '01' in setup_dict.keys()
    assert setup_dict['01']['disperser']['name'] == '600/4310'

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
    pargs2 = setup.parser([pypit_file, 'kast_blue', '--pypit_file', '-d'])
    setup.main(pargs2)
    pytest.set_trace()

'''
def test_pypit_file():
    """ Generate .pypit files.
    Uses .pypit file from previous test!
    """
    from pypit.scripts import pypit_files
    pypit_file = glob.glob('kast_blue*.pypit')[0]
    pargs = pypit_files.parser([pypit_file, 'kast_blue'])
    pypit_files.main(pargs)
    assert os.path.isfile('./kast_blue_setup_01.pypit')

def test_view_fits():
    """ Only test the list option
    """
    from pypit.scripts import view_fits
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    pargs = view_fits.parser([spec_file, '--list'])
