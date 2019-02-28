"""
Module to run tests on arcoadd
"""
import os

import pytest
import numpy as np

from pypeit import msgs
from pypeit.tests.tstutils import dev_suite_required


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@dev_suite_required
def test_cooked_version():
    # Load up the version
    v_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'version')
    with open(v_file) as f:
        tmp = f.readlines()
    value = float(tmp[-1].strip())
    # Test
    assert value >= 1.0

