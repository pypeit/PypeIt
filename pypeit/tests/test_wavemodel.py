"""
Module to run tests on wavemodel.py
"""

import numpy as np
from pypeit import wavemodel
import pytest

try:
    tsterror = FileExistsError
except NameError:
    FileExistsError = OSError

# TODO -- Note that the transparency function is not used anywhere

def test_transparency():
    """ Test for creting the ski transmission model. It is basically
    testing if the skisim directory is reachable.
    """

    wave_test = np.arange(0.7,2.9,0.001)
    tran_test = wavemodel.transparency(wave_test, debug=False)

    assert np.max(tran_test) - np.min(tran_test) == 1.
