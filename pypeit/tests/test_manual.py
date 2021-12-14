""" Tests for ManualExtraction object """
import pytest

import numpy as np

from pypeit.manual_extract import ManualExtractionObj

def test_instantiate():

    # Init
    mex = ManualExtractionObj(frame='tst.fits',
                              det = np.arange(10)+1,
                              spat = np.arange(10).astype(float),
                              spec = np.arange(10).astype(float),
                              fwhm = np.arange(10).astype(float))

    # FITS table input
    mex2 = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '1:1181.8:3820.6:3.')
    mex3 = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '1:1181.8:3820.6:3.,-1:1183.8:3820.6:3.')

    # Test
    assert np.all(mex3.det == np.array([1,-1]))
                                            
def test_dict_for_obj():
    mex3 = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '1:1181.8:3820.6:3.,-1:1183.8:3820.6:3.')

    dobj = mex3.dict_for_objfind(1, neg=False)
    assert dobj['det'] == 1

    dobj2 = mex3.dict_for_objfind(1, neg=True)
    assert dobj2['det'] == 1
    assert np.isclose(dobj2['spat'][0], 1183.8)