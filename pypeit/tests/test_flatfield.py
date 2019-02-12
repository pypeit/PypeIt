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


from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters
from pypeit import flatfield
from pypeit.par import pypeitpar

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


# TODO: Bring this test back in some way?
#def test_step_by_step():
#    if skip_test:
#        assert True
#        return
#    # Masters
#    spectrograph, TSlits, tilts, datasec_img \
#                = load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True,
#                                         datasec=True)
#    # Instantiate
#    flatField = flatfield.FlatField(spectrograph, det=1, tilts=tilts,
#                                    tslits_dict=TSlits.tslits_dict.copy())
#    # Use mstrace
#    flatField.mspixelflat = TSlits.mstrace.copy()
#    # Normalize a slit
#    slit=0
#    flatField._prep_tck()
#    modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = flatField.slit_profile(slit)
#    assert np.isclose(iextrap_slit, 0.)
#    # Apply
#    word = np.where(flatField.tslits_dict['slitpix'] == slit + 1)
#    flatField.mspixelflatnrm = flatField.mspixelflat.copy()
#    flatField.mspixelflatnrm[word] /= nrmvals
#    assert np.isclose(np.median(flatField.mspixelflatnrm), 1.0267346)


@dev_suite_required
def test_run():
    # Masters
    spectrograph, tslits_dict, tilts_dict, datasec_img \
                = load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True,
                                         datasec=True)
    # Instantiate
    frametype = 'pixelflat'
    par = pypeitpar.FrameGroupPar(frametype)
    flatField = flatfield.FlatField(spectrograph, par, det=1, tilts_dict=tilts_dict,
                                    tslits_dict=tslits_dict.copy())

    # TODO mstrace is no longer stored in the master file so this needs to be run from an input file which must
    # be added to the cooked directory. Assigning this one to @profx.  Disabling until it is fixed.

    # Use mstrace
    #flatField.rawflatimg = tslits_dict['mstrace'].copy()
    #mspixelflatnrm, msillumflat = flatField.run()
    #assert np.isclose(np.median(mspixelflatnrm), 1.0)

