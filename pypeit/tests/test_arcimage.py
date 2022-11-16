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
from pypeit.images import detector_container
from pypeit.tests.test_detector import def_det
from pypeit.tests.tstutils import data_path

def test_init():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
    pypeitImage.detector = detector_container.DetectorContainer(**def_det)
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage)

def test_master_io():
    # Instantiate a simple pypeitImage
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
    pypeitImage.detector = detector_container.DetectorContainer(**def_det)
    pypeitImage.PYP_SPEC = 'shane_kast_blue'
    # Now the arcimage
    arcImage = buildimage.ArcImage.from_pypeitimage(pypeitImage)
    # Write
    master_filename = masterframe.construct_file_name(arcImage, 'A_01_22', master_dir=data_path(''))
    assert master_filename == data_path('MasterArc_A_01_22.fits'), 'Master filename changed'
    arcImage.to_master_file(master_filename)
    # Read
    _arcImage = buildimage.ArcImage.from_file(master_filename)
    assert isinstance(_arcImage.detector, detector_container.DetectorContainer), \
            'detector has wrong type'
    assert np.array_equal(_arcImage.fullmask.mask, arcImage.fullmask.mask), 'mask changed'
    # Cleanup
    os.remove(master_filename)


