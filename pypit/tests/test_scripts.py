# Module to run tests on scripts

import matplotlib
matplotlib.use('Agg')  # For Travis

# TEST_UNICODE_LITERALS

import sys, os
import pytest

skip = pytest.mark.skip

from pypit.scripts import arcid_plot, show_1dspec, view_fits
from pypit import pyputils
msgs = pyputils.get_dummy_logger()

from PyQt4 import QtGui
app = QtGui.QApplication(sys.argv)

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_arcid_plot():
    json_file = data_path('LRISb_600_WaveCalib_01.json')
    pargs = arcid_plot.parser([json_file, 'LRISb', 'tmp.pdf'])
    # Run
    arcid_plot.main(pargs)


"""
@skip
def test_show_1dspec():
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    pargs = show_1dspec.parser([spec_file])
    # Run
    show_1dspec.main(pargs, unit_test=True)
"""

def test_view_fits():
    """ Only test the list option
    """
    spec_file = data_path('spec1d_J0025-0312_KASTr_2015Jan23T025323.85.fits')
    pargs = view_fits.parser([spec_file, '--list'])
