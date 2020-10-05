"""
Module to run tests on arcoadd
"""
import os

import pytest
import numpy as np

from pypeit import msgs
from pypeit.tests.tstutils import cooked_required
from IPython import embed


@cooked_required
def test_cooked_version():
    # TODO: Use pygit2 to both add the branch name when constructing
    # Cooked and here to check that the branch is correct?

    # Load up the version
    v_file = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'version')
    with open(v_file) as f:
        tmp = f.readlines()
    assert tmp[-1].strip() == '1.1.1'


