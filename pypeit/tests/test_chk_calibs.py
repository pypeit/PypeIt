"""
Module to run tests on scripts
"""
import os
import shutil

from IPython import embed
import pytest

import numpy as np

import matplotlib
matplotlib.use('agg')  # For Travis

from pypeit.scripts import chk_for_calibs
from pypeit.tests.tstutils import dev_suite_required, cooked_required, data_path

@dev_suite_required
def test_chk_calibs_not():
    os.chdir(data_path(''))
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/not_alfosc/grism4')
    droot += '/ALD'

    pargs = chk_for_calibs.parse_args([droot, '-s', 'not_alfosc'])
    answers, _ = chk_for_calibs.main(pargs)

    assert answers['pass'][0], 'One or more failures!'


@dev_suite_required
def test_chk_calibs_deimos():
    os.chdir(data_path(''))
    # 830G
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/830G_M_8600/')
    pargs = chk_for_calibs.parse_args([droot, '-s', 'keck_deimos'])
    answers, _ = chk_for_calibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # 600ZD
    droot = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_deimos/600ZD_M_7500/')
    pargs = chk_for_calibs.parse_args([droot, '-s', 'keck_deimos'])
    answers, _ = chk_for_calibs.main(pargs)
    assert np.all(answers['pass'][:-1]), 'One or more failures!'

    # Cleanup
    # NOTE: chk_for_calibs now removes the setup files by default
    #shutil.rmtree(data_path('setup_files'))


