"""
Module to run tests on scripts
"""
import os
import numpy as np
import pytest
import shutil

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit.scripts import chk_calibs
from pypeit.tests.tstutils import dev_suite_required, cooked_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@dev_suite_required
def test_chk_calibs_not():
    os.chdir(data_path(''))
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'
    pargs = chk_calibs.parser([droot, '-s', 'not_alfosc'])
    answers, _ = chk_calibs.main(pargs)

    assert answers['pass'][0], 'One or more failures!'

    # Cleanup
    shutil.rmtree(data_path('setup_files'))

@dev_suite_required
def test_chk_calibs_deimos():
    os.chdir(data_path(''))
    # 830G
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600/')
    pargs = chk_calibs.parser([droot, '-s', 'keck_deimos'])
    answers, _ = chk_calibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # 600ZD
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/600ZD_M_7500/')
    pargs = chk_calibs.parser([droot, '-s', 'keck_deimos'])
    answers, _ = chk_calibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # Cleanup
    shutil.rmtree(data_path('setup_files'))

