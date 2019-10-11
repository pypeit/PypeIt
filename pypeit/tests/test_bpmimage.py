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



def test_dummy_image():
    # Simple
    shape=(2048,2048)
    spectrograph = util.load_spectrograph('shane_kast_blue')
    bpm = spectrograph.bpm(shape=shape)#, trim=False)
    assert isinstance(bpm, np.ndarray)
    assert bpm.shape == shape
    assert np.sum(bpm) == 0


@dev_suite_required
def test_keck_lris_red():
    # Spectrograph
    spectrograph = util.load_spectrograph('keck_lris_red')
    #
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    # Get the shape
    det = 2
    rdsec_img = spectrograph.get_rawdatasec_img(example_file, det=det)
    trim = procimg.trim_frame(rdsec_img, rdsec_img < 1)
    orient = spectrograph.orient_image(trim, det)
    shape = orient.shape
    # Simple
    bpm = spectrograph.bpm(shape=shape, filename=example_file, det=det)
    assert np.sum(bpm) > 0


@dev_suite_required
def test_keck_deimos():
    spectrograph = util.load_spectrograph('keck_deimos')
    example_file = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L_8400',
                                'd0914_0002.fits.gz')
    # Get the shape
    det = 2
    rdsec_img = spectrograph.get_rawdatasec_img(example_file, det=det)
    trim = procimg.trim_frame(rdsec_img, rdsec_img < 1)
    orient = spectrograph.orient_image(trim, det)
    shape = orient.shape
    # Simple
    bpm = spectrograph.bpm(shape=shape,det=4)
    assert bpm[0,0] == 1

