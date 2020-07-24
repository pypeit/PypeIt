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
    pargs = chk_calibs.parser([droot, '-s', 'not_alfosc'])
    answers, _ = chk_calibs.main(pargs)

    assert answers['pass'][0], 'One or more failures!'

@dev_suite_required
def test_chk_calibs_deimos():
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/*/')
    pargs = chk_calibs.parser([droot, '-s', 'keck_deimos'])
    answers, _ = chk_calibs.main(pargs)

    assert np.all(answers['pass'][:-1]), 'One or more failures!'

