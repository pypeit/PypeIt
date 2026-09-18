""" Tests for ManualExtraction object """
from IPython import embed
import pytest

import numpy as np

from astropy.table import Table

from pypeit.manual_extract import ManualExtractionObj
from pypeit.spectrographs.util import load_spectrograph
from pypeit import PypeItError
from pypeit.pypeit_steps import get_manual_flexure

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


def test_manual_flexure():
    """
    Test the manual spatial-flexure sentinel (Issue #2180): a blank
    ``shift`` entry means "no manual value"; any numeric entry
    (including 0.) is a user-requested override.
    """
    # Mimic the fitstbl: blank default, a finite override, an explicit
    # 0. override (still "manual"), and a literal 'None'
    tbl = Table()
    tbl['filename'] = ['a.fits', 'b.fits', 'c.fits', 'd.fits']
    tbl['shift'] = ['', '2.5', '0.', 'None']

    # Blank default -> no manual flexure
    assert get_manual_flexure(tbl, 0) is None, \
        'Blank shift should mean no manual flexure'
    # Numeric value -> manual flexure
    assert get_manual_flexure(tbl, 1) == 2.5, \
        'Numeric shift should be returned as the manual flexure'
    # 0. is numeric and therefore a valid manual override
    assert get_manual_flexure(tbl, 2) == 0., \
        'A shift of 0. should still count as a manual flexure'
    # Literal None -> no manual flexure
    assert get_manual_flexure(tbl, 3) is None, \
        'A shift of "None" should mean no manual flexure'

    # Missing column -> no manual flexure
    tbl_noshift = Table()
    tbl_noshift['filename'] = ['a.fits']
    assert get_manual_flexure(tbl_noshift, 0) is None, \
        'Missing shift column should mean no manual flexure'

    # A non-numeric entry raises an error
    tbl_bad = Table()
    tbl_bad['filename'] = ['a.fits']
    tbl_bad['shift'] = ['bad']
    with pytest.raises(PypeItError):
        get_manual_flexure(tbl_bad, 0)
