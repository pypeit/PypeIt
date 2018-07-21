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

from pypeit.tests import tstutils
from pypeit import flatfield

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
    spectrograph, TSlits, tilts, datasec_img \
                = tstutils.load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True,
                                                  datasec=True)
    # Instantiate
    flatField = flatfield.FlatField(spectrograph=spectrograph, det=1, tilts=tilts,
                                    tslits_dict=TSlits.tslits_dict.copy())
    # Use mstrace
    flatField.mspixelflat = TSlits.mstrace.copy()
    # Normalize a slit
    slit=0
    flatField._prep_tck()
    modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = flatField.slit_profile(slit)
    assert np.isclose(iextrap_slit, 0.)
    # Apply
    word = np.where(flatField.tslits_dict['slitpix'] == slit + 1)
    flatField.mspixelflatnrm = flatField.mspixelflat.copy()
    flatField.mspixelflatnrm[word] /= nrmvals
    assert np.isclose(np.median(flatField.mspixelflatnrm), 1.0291458)

def test_run():
    if skip_test:
        assert True
        return
    # Masters
    spectrograph, TSlits, tilts, datasec_img \
                = tstutils.load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True,
                                                  datasec=True)
    # Instantiate
    flatField = flatfield.FlatField(spectrograph=spectrograph, det=1, tilts=tilts,
                                    tslits_dict=TSlits.tslits_dict.copy())
    # Use mstrace
    flatField.mspixelflat = TSlits.mstrace.copy()
    mspixelflatnrm = flatField.run()
    assert np.isclose(np.median(mspixelflatnrm), 1.0086422)

