"""
Module to run tests on scripts
"""
import os
import numpy as np
import pytest

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit.scripts import chk_calibs
from pypeit.tests.tstutils import dev_suite_required, cooked_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@dev_suite_required
def test_chk_calibs_not():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'
    pargs = chk_calibs.parser(['-r', droot, '-s', 'not_alfosc'])
    answers = chk_calibs.main(pargs)

    assert np.all(answers), 'One or more failures!'

