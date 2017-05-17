# Module to run tests on scripts

import matplotlib
matplotlib.use('agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest
import glob

from pypit import pyputils
from pypit.scripts import coadd_1dspec

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


def test_view_fits():
    """ Only test the list option
    """
    from pypit.scripts import view_fits
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    pargs = view_fits.parser([spec_file, '--list'])

def test_coadd():
    coadd_file = data_path('coadd_UGC3672A_red.yaml')
    args = coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files = coadd_1dspec.main(args, unit_test=True)
    # Test
    assert len(gparam) == 0
    assert isinstance(gparam, dict)
    assert ex_value == 'opt'
    assert flux_value is True
    assert iobj == 'O210-S1467-D02-I0012'
    assert outfile == 'UGC3672A_r.fits'
    assert len(files) == 4

def test_coadd2():
    """ Test using a list of object names
    """
    coadd_file = data_path('coadd_UGC3672A_red_objlist.yaml')
    args = coadd_1dspec.parser([coadd_file])
    # Main
    gparam, ex_value, flux_value, iobj, outfile, files = coadd_1dspec.main(args, unit_test=True)
    # Test
    assert len(iobj) == len(files)
    # Crash it
    coadd_file = data_path('coadd_UGC3672A_red_badlist.yaml')
    args = coadd_1dspec.parser([coadd_file])
    with pytest.raises(IOError):
        gparam, ex_value, flux_value, iobj, outfile, files = coadd_1dspec.main(args, unit_test=True)

