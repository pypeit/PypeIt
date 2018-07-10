# Module to run tests on TraceSlits class
#   Requires files in Development suite and an Environmental variable
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

import os

import pytest
import glob
import numpy as np

from pypit import traceslits

# These tests are not run on Travis
if os.getenv('PYPIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def chk_for_files(root):
    files = glob.glob(root+'*')
    if len(files) == 0:
        return False
    else:
        return True


def test_load_from_master_and_run():
    if skip_test:
        assert True
        return
    # Check for files
    mstrace_root = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_150420_402')
    assert chk_for_files(mstrace_root)
    # Load
    traceSlits = traceslits.TraceSlits.from_master_files(mstrace_root)
    assert isinstance(traceSlits.mstrace, np.ndarray)


def test_add_slit():
    if skip_test:
        assert True
        return
    # Check for files
    mstrace_root = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_150420_402')
    assert chk_for_files(mstrace_root)
    # Load
    traceSlits = traceslits.TraceSlits.from_master_files(mstrace_root)
    norig = traceSlits.nslit
    #  left edge, right edge, row on image
    add_user_slits = [[489, 563, 1024]]
    # run_to_finish resets things in a proper manner
    traceSlits.add_user_slits(add_user_slits, run_to_finish=True)
    assert traceSlits.nslit == (norig+1)


def test_remove_slit():
    if skip_test:
        assert True
        return
    # Check for files
    mstrace_root = os.path.join(os.getenv('PYPIT_DEV'), 'Cooked', 'Trace',
                                'MasterTrace_KeckLRISr_20160110_A')
    assert chk_for_files(mstrace_root)
    # Load
    traceSlits = traceslits.TraceSlits.from_master_files(mstrace_root)
    norig = traceSlits.nslit
    # Setup slit to remove --  xleft, yleft at yrow=nrow/2
    rm_slits = [[229, 380]]
    # Remove
    traceSlits.remove_slit(rm_slits)
    assert traceSlits.nslit == (norig-1)

