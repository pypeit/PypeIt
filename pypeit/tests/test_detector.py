"""
Module to run tests on SlitTraceSet
"""
import os
import pytest

import numpy as np

from pypeit.images import detector_container

def test_init():
    detector = detector_container.Detector(
        dataext=0,
        specaxis=1,
        specflip=False,
        spatflip=False,
        platescale=0.43,
        saturation=65535.,
        mincounts=-1e10,
        nonlinear=0.76,
        numamplifiers=2,
        gain=[1.2, 1.2],
        ronoise=[3.7, 3.7],
        det=1,
        xgap=0.,
        ygap=0.,
        ysize=1.,
        darkcurr=0.0,
        datasec=['[:, 1:1024]', '[:, 1025:2048]'],  # These are rows, columns on the raw frame, 1-indexed
        oscansec=['[:, 2050:2080]', '[:, 2081:2111]'])
        #suffix='_blue'

