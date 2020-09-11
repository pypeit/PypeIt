"""
Module to run tests on BPMImage class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import dev_suite_required
from pypeit.spectrographs import util
from pypeit.core import procimg


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@dev_suite_required
def test_keck_lris_red():
    # Spectrograph
    spectrograph = util.load_spectrograph('keck_lris_red')
    #
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_lris_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    # Get the shape
    det = 2
    # Simple
    bpm = spectrograph.bpm(example_file, det)
    assert np.sum(bpm) > 0


@dev_suite_required
def test_keck_deimos():
    spectrograph = util.load_spectrograph('keck_deimos')
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '830G_L_8400',
                                'd0914_0002.fits.gz')
    # Get the shape
    det = 4
    # Simple
    bpm = spectrograph.bpm(example_file, det)
    assert bpm[0,0] == 1

