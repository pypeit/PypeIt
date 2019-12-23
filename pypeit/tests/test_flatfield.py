"""
Module to run tests on FlatField class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import dev_suite_required, load_kast_blue_masters, cooked_required
from pypeit import flatfield
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.images import pypeitimage

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

@cooked_required
def test_run():
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    edges, tilts_dict = load_kast_blue_masters(edges=True, tilts=True)
    # Instantiate
    frametype = 'pixelflat'
    par = pypeitpar.FrameGroupPar(frametype)
    flatField = flatfield.FlatField(spectrograph, par, det=1, tilts_dict=tilts_dict,
                                    tslits_dict=edges.convert_to_tslits_dict())

    # Use the trace image
    flatField.rawflatimg = pypeitimage.PypeItImage(edges.img.copy())
    mspixelflatnrm, msillumflat = flatField.run()
    assert np.isclose(np.median(mspixelflatnrm), 1.0)

