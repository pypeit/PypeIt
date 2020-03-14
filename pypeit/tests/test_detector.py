"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.images import detector_container
from pypeit import io

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def_det = dict(
    dataext=0,
    specaxis=1,
    specflip=False,
    spatflip=False,
    platescale=0.43,
    saturation=65535.,
    mincounts=-1e10,
    nonlinear=0.76,
    numamplifiers=2,
    gain=np.asarray([1.2, 1.2]),
    ronoise=np.asarray([3.7, 3.7]),
    det=1,
    xgap=0.,
    ygap=0.,
    ysize=1.,
    darkcurr=0.0,
    binning='1,1',
    datasec=np.asarray(['[:, 1:1024]', '[:, 1025:2048]']),  # These are rows, columns on the raw frame, 1-indexed
    oscansec=np.asarray(['[:, 2050:2080]', '[:, 2081:2111]']))


def test_init():
    detector = detector_container.Detector(**def_det)

def test_bundle():
    detector = detector_container.Detector(**def_det)
    data = detector._bundle()
    assert len(data) == 1

def test_to_hdu():
    detector = detector_container.Detector(**def_det)
    hdul = detector.to_hdu()
    # Parse me
    dd, tst1, tst2 = detector._parse(hdul[0])

def test_write():
    detector = detector_container.Detector(**def_det)
    detector.to_file(data_path('tmp.fits'), overwrite=True)


def test_read():
    detector = detector_container.Detector(**def_det)
    new_detector = detector.from_file(data_path('tmp.fits'))

