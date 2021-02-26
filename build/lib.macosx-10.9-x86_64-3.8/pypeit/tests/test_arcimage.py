"""
Module to run tests on ArcImage as a test
of master frame functionality
"""
import os

import pytest
import glob

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit.images import pypeitimage
from pypeit.images import buildimage
from pypeit import masterframe
from pypeit.tests import test_detector

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_init():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.fullmask = np.zeros((1000, 1000), dtype=np.int64)
    pypeitImage.detector = test_detector.detector_container.DetectorContainer(**test_detector.def_det)
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage)

def test_master_io():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.fullmask = np.zeros((1000, 1000), dtype=np.int64)
    pypeitImage.detector = test_detector.detector_container.DetectorContainer(**test_detector.def_det)
    pypeitImage.PYP_SPEC = 'shane_kast_blue'
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage)
    # Write
    master_filename = masterframe.construct_file_name(arcImage, 'A_01_22', master_dir=data_path(''))
    arcImage.to_master_file(master_filename)
    # Read
    _arcImage = buildimage.ArcImage.from_file(data_path('MasterArc_A_01_22.fits'))
    assert isinstance(_arcImage.detector, test_detector.detector_container.DetectorContainer)
    # Cleanup
    os.remove(master_filename)

