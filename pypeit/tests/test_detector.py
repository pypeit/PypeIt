"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.tests.tstutils import data_path
from pypeit.images import detector_container


# Example (shane_kast_blue)
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
    detector = detector_container.DetectorContainer(**def_det)
    assert detector.specaxis == 1


def test_name():
    detector = detector_container.DetectorContainer(**def_det)
    assert detector.det == detector_container.DetectorContainer.parse_name(detector.name), \
            'name parsing mismatch'


def test_bundle():
    detector = detector_container.DetectorContainer(**def_det)
    data = detector._bundle()
    assert len(data) == 1


def test_io():
    detector = detector_container.DetectorContainer(**def_det)
    detector.to_file(data_path('tmp_detector.fits'), overwrite=True)

    _new_detector = detector.from_file(data_path('tmp_detector.fits'))

    # Check a few attributes are equal
    assert detector['dataext'] == _new_detector['dataext'], 'Bad read dataext'
    assert np.array_equal(detector['gain'], _new_detector['gain']), 'Bad read gain'
    assert detector['binning'] == _new_detector['binning'], 'Bad read binning'
    assert np.array_equal(detector['datasec'], _new_detector['datasec']), 'Bad read datasec'

    os.remove(data_path('tmp_detector.fits'))

