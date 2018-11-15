# Module to run tests on wavemodel.py

import numpy as np
from pypeit.wavemodel import transparency
import pytest

try:
    tsterror = FileExistsError
except NameError:
    FileExistsError = OSError

def test_transparency():
    """ Test for creting the ski transmission model
    """

    wave_test = np.arange(0.7,2.9,0.001)
    tran_test = transparency(wave_test)

    assert np.max(tran_test) - np.min(tran_test) == 1.
