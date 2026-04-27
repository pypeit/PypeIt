"""
Module to run tests on SlitTraceSet
"""
import os

from IPython import embed

import numpy as np

from pypeit.tests import tstutils
from pypeit.images import detector_container


def test_init():
    detector = detector_container.DetectorContainer(**tstutils.default_detector())
    assert detector.specaxis == 1


def test_name():
    detector = detector_container.DetectorContainer(**tstutils.default_detector())
    assert detector.det == detector_container.DetectorContainer.parse_name(detector.name), \
            'name parsing mismatch'


def test_bundle():
    detector = detector_container.DetectorContainer(**tstutils.default_detector())
    data = detector._bundle()
    assert len(data) == 1


def test_io():
    detector = detector_container.DetectorContainer(**tstutils.default_detector())
    detector.to_file(tstutils.data_output_path('tmp_detector.fits'), overwrite=True)

    _new_detector = detector.from_file(tstutils.data_output_path('tmp_detector.fits'))

    # Check a few attributes are equal
    assert detector['dataext'] == _new_detector['dataext'], 'Bad read dataext'
    assert np.array_equal(detector['gain'], _new_detector['gain']), 'Bad read gain'
    assert detector['binning'] == _new_detector['binning'], 'Bad read binning'
    assert np.array_equal(detector['datasec'], _new_detector['datasec']), 'Bad read datasec'

    os.remove(tstutils.data_output_path('tmp_detector.fits'))


def test_copy():
    detector = detector_container.DetectorContainer(**tstutils.default_detector())
    detcopy = detector.copy()

    # Check that a couple of relevant attributes have different memory locations
    assert detector is not detcopy, 'Should not point to the same reference'
    assert detector.gain is not detcopy.gain, \
        'numpy array attributes should not point to the same reference'
    
