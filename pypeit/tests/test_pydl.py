"""
Module to run tests on pyidl functions
"""

import numpy as np
from pypeit.core.pydl import bspline
import pytest

try:
    tsterror = FileExistsError
except NameError:
    FileExistsError = OSError


def test_bsplinetodict():
    """ Test for writing a bspline onto a dict
    (and also reading it out).
    """

    x = np.random.rand(500)

    # Create bspline
    init_bspline = bspline(x, bkspace=0.01*(np.max(x)-np.min(x)))
    # Write bspline to bspline_dict
    bspline_dict = init_bspline.to_dict()
    # Create bspline from bspline_dict
    bspline_fromdict = bspline(None, from_dict=bspline_dict)

    assert np.max(np.array(bspline_dict['breakpoints'])-bspline_fromdict.breakpoints) == 0.

