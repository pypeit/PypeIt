""" Tests for ManualExtraction object """
import pytest

import numpy as np

from pypeit.manual_extract import ManualExtractionObj
from pypeit.spectrographs.util import load_spectrograph

def test_instantiate():

    # Init
    mex = ManualExtractionObj(frame='tst.fits',
                              detname=np.array(['DET01']*10),
                              spat=np.arange(10).astype(float),
                              spec=np.arange(10).astype(float),
                              fwhm=np.arange(10).astype(float))

    # FITS table input
    spectrograph = load_spectrograph('keck_deimos')
    mex2 = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.', spectrograph)
    mex3 = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.;-1:1183.8:3820.6:3.', spectrograph)
    mex4 = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.:4.;2:1183.8:3820.6:3.', spectrograph)

    # Test
    assert np.all(mex3.detname == np.array(['DET01','DET01']))

    # Boxcar
    assert np.all(mex3.boxcar_rad == np.array([-1.]))
    assert np.all(mex4.boxcar_rad == np.array([4., -1.]))


def test_dict_for_obj():
    spectrograph = load_spectrograph('keck_deimos')
    mex3 = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.;-1:1183.8:3820.6:3.', spectrograph)

    dobj = mex3.dict_for_objfind('DET01', neg=False)
    assert dobj['detname'] == 'DET01'

    dobj2 = mex3.dict_for_objfind('DET01', neg=True)
    assert dobj2['detname'] == 'DET01'
    assert np.isclose(dobj2['spat'][0], 1183.8)
