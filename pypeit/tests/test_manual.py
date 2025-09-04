""" Tests for ManualExtraction object """
from IPython import embed
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
    assert mex.neg is None, 'neg not defined'

    # FITS table input
    spectrograph = load_spectrograph('keck_deimos')

    mex = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.', spectrograph)
    assert not mex.neg[0], 'neg should be false'

    mex = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '1:1181.8:3820.6:3.;-1:1183.8:3820.6:3.',
                                               spectrograph)
    assert np.all(mex.detname == np.array(['DET01','DET01']))
    assert np.all(mex.boxcar_rad == np.array([-1.]))

    mex = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '1:1181.8:3820.6:3.:4.;2:1183.8:3820.6:3.',
                                               spectrograph)
    assert np.all(mex.boxcar_rad == np.array([4., -1.]))

    # Mosaic
    mex = ManualExtractionObj.by_fitstbl_input('tst.fits', '(1,5):1181.8:3820.6:3.', spectrograph)
    assert mex.detname[0] == 'MSC01'
    mex = ManualExtractionObj.by_fitstbl_input('tst.fits', '(2,6):1181.8:3820.6:3.', spectrograph)
    assert mex.detname[0] == 'MSC02'
    mex = ManualExtractionObj.by_fitstbl_input('tst.fits',
                                               '(1,5):1181.8:3820.6:3.;(2,6):1181.8:3820.6:3.',
                                               spectrograph)
    assert mex.detname[0] == 'MSC01'
    assert mex.detname[1] == 'MSC02'


def test_dict_for_obj():
    spectrograph = load_spectrograph('keck_deimos')
    mex3 = ManualExtractionObj.by_fitstbl_input('tst.fits', '1:1181.8:3820.6:3.;-1:1183.8:3820.6:3.', spectrograph)

    dobj = mex3.dict_for_objfind('DET01', neg=False)
    assert dobj['detname'] == 'DET01'

    dobj2 = mex3.dict_for_objfind('DET01', neg=True)
    assert dobj2['detname'] == 'DET01'
    assert np.isclose(dobj2['spat'][0], 1183.8)
