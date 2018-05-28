# Module to run tests on FlatField class
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

from astropy.table import Table

from pypit.tests import tstutils
from pypit import flatfield

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

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


def test_step_by_step():
    if skip_test:
        assert True
        return
    # Masters
    settings, TSlits, tilts, datasec_img = tstutils.load_kast_blue_masters(
        get_settings=True, tslits=True, tilts=True, datasec=True)
    # Instantiate
    ftField = flatfield.FlatField(spectrograph='shane_kast_blue', settings=settings, det=1,
                                  tilts=tilts, slits_dict=TSlits.slits_dict.copy())
    # Use mstrace
    ftField.mspixelflat = TSlits.mstrace.copy()
    # Gain
    ftField.apply_gain(datasec_img)
    # Normalize a slit
    slit=0
    modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = ftField.slit_profile(slit)
    assert np.isclose(iextrap_slit, 0.)
    # Apply
    word = np.where(ftField.slits_dict['slitpix'] == slit + 1)
    ftField.mspixelflatnrm = ftField.mspixelflat.copy()
    ftField.mspixelflatnrm[word] /= nrmvals
    assert np.isclose(np.median(ftField.mspixelflatnrm), 1.0011169)

def test_run():
    if skip_test:
        assert True
        return
    # Masters
    settings, TSlits, tilts, datasec_img = tstutils.load_kast_blue_masters(
        get_settings=True, tslits=True, tilts=True, datasec=True)
    # Instantiate
    ftField = flatfield.FlatField(spectrograph='shane_kast_blue', settings=settings, det=1,
                                  tilts=tilts, slits_dict=TSlits.slits_dict.copy())
    # Use mstrace
    ftField.mspixelflat = TSlits.mstrace.copy()
    mspixelflatnrm = ftField.run(datasec_img)
    assert np.isclose(np.median(mspixelflatnrm), 1.003463)

