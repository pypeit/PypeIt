""" Tests for ManualExtraction object """
from IPython import embed
import pytest

import numpy as np

from astropy.table import Table

from pypeit.manual_extract import ManualExtractionObj
from pypeit.spectrographs.util import load_spectrograph
from pypeit import PypeItError
from pypeit.pypeit_steps import get_manual_flexure
from pypeit.metadata import PypeItMetaData, parse_manual_spat_flexure

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


def test_parse_manual_spat_flexure():
    # Blank, None, nan, and masked entries mean "not set"
    for value in ['', '  ', 'None', 'none', 'nan', None, np.ma.masked]:
        assert parse_manual_spat_flexure(value) is None, \
            f'{value!r} should mean no manual flexure'
    # Numbers (including 0) are manual overrides
    assert parse_manual_spat_flexure(' -1.5 ') == -1.5, 'String number not parsed'
    assert parse_manual_spat_flexure(0.) == 0., 'A shift of 0. should count as a manual flexure'
    assert parse_manual_spat_flexure(3) == 3., 'Integer shift not parsed'
    # Anything else raises an error
    for value in ['bad', 'inf']:
        with pytest.raises(PypeItError):
            parse_manual_spat_flexure(value)


def test_get_shifts():
    # Manual spatial flexure for multiple frames, as used to build the
    # science and calibration images
    spec = load_spectrograph('shane_kast_blue')
    data = Table()
    data['filename'] = ['b1.fits', 'b27.fits', 'b28.fits', 'b29.fits']
    data['shift'] = ['', '2.5', '0', 'None']
    fitstbl = PypeItMetaData(spec, spec.default_pypeit_par(), data=data)

    shifts = fitstbl.get_shifts([0, 1, 2, 3])
    assert isinstance(shifts, np.ma.MaskedArray), 'Shifts should be a masked array'
    assert np.array_equal(shifts.mask, [True, False, False, True]), \
        'Only the frames with a numeric shift should be unmasked'
    assert np.array_equal(shifts.compressed(), [2.5, 0.]), 'Incorrect shifts'
    # Single index and boolean selection
    assert fitstbl.get_shifts(1)[0] == 2.5, 'Incorrect shift for a single frame'
    assert np.array_equal(fitstbl.get_shifts(np.array([False, True, True, False])), [2.5, 0.]), \
        'Incorrect shifts for a boolean selection'
    # Individual elements are masked when not set, as expected by
    # RawImage.spatial_flexure_shift
    assert np.ma.is_masked(fitstbl.get_shifts([0])[0]), 'Unset shift should be masked'

    # No shift column means no manual flexure for any frame
    fitstbl.table.remove_column('shift')
    assert np.all(fitstbl.get_shifts([0, 1]).mask), 'Missing shift column should mask all shifts'


@pytest.mark.remote_data
def test_process_manual_flexure_with_skip():
    # A manual spatial flexure is applied when processing a frame, even if the
    # spatial flexure method is 'skip'
    from pypeit import dataPaths
    from pypeit.images.rawimage import RawImage
    from pypeit.slittrace import SlitTraceSet
    from pypeit.spectrographs.util import load_spectrograph
    spectrograph = load_spectrograph('shane_kast_blue')
    par = spectrograph.default_pypeit_par()['scienceframe']['process']
    for key in ['use_biasimage', 'use_darkimage', 'use_pixelflat', 'use_illumflat',
                'use_specillum', 'subtract_scattlight']:
        par[key] = False
    par['spat_flexure_method'] = 'skip'
    raw_file = str(dataPaths.tests.get_file_path('b24.fits.gz'))
    # Without a manual value, no spatial flexure is determined
    image = RawImage(raw_file, spectrograph, 1).process(par)
    assert image.spat_flexure is None, 'No spatial flexure should be determined'
    # A single slit covering the middle of the detector
    nspec, nspat = image.shape
    left = np.full((nspec, 1), 0.3 * nspat)
    slits = SlitTraceSet(left_init=left, right_init=left + 0.4 * nspat, pypeline='MultiSlit',
                         nspat=nspat, PYP_SPEC=spectrograph.name)
    image = RawImage(raw_file, spectrograph, 1).process(par, slits=slits)
    assert image.spat_flexure is None, 'No spatial flexure should be determined'
    for manual in [1.5, 0.0]:
        image = RawImage(raw_file, spectrograph, 1).process(par, slits=slits,
                                                            manual_spat_flexure=manual)
        assert np.array_equal(image.spat_flexure, np.full((1, 2), manual)), \
            f'The manual spatial flexure ({manual}) should be applied with method=skip'
