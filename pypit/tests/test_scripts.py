# Module to run tests on scripts

import matplotlib
matplotlib.use('Agg')  # For Travis

# TEST_UNICODE_LITERALS

import os
import pytest

from pypit.scripts import arcid_plot

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_arcid_plot():
    options = data_path('LRISb_600_WaveCalib_01.json')
    pargs = arcid_plot.parser([options, 'LRISb', 'tmp.pdf'])
    # Run
    arcid_plot.main(pargs)

