"""
Module to run tests on sort and arsetup
"""
from IPython import embed

import pytest

import numpy as np

from pypeit.calibframe import CalibFrame
from pypeit.core import framematch
from pypeit.tests.tstutils import dummy_fitstbl
from pypeit.pypmsgs import PypeItError


@pytest.fixture
def fitstbl():
    return dummy_fitstbl()


def test_frame_selection(fitstbl):
    """ Test that the frame bits are successfully read
    """
    # Sort
    assert fitstbl.find_frames('bias')[0]
    assert fitstbl.find_frames('arc')[1]
    assert fitstbl.find_frames('trace')[2]
    assert fitstbl.find_frames('standard')[4]
    assert np.sum(fitstbl.find_frames('science')) == 5


def test_calibration_groups(fitstbl):
    """
    Test the frame selection specific to a provided calibration group
    """
    calib_ID = 0
    par = fitstbl.spectrograph.default_pypeit_par()
    assert fitstbl.find_frames('arc', calib_ID=calib_ID, index=True)[0] == 1
    assert fitstbl.find_frames('standard', calib_ID=calib_ID, index=True)[0] == 4
    assert fitstbl.find_frames('trace', calib_ID=calib_ID, index=True)[0] == 2


# TODO: This doesn't test anything
#def test_neg_match_science(fitstbl):
#    """ Test using negative number for calibs
#    """
#    par = fitstbl.spectrograph.default_pypeit_par()
#    # Use negative number
#    for ftype in ['arc', 'pixelflat', 'bias']:
#        par['calibrations']['{0}frame'.format(ftype)]['number'] = 1
#    par['calibrations']['traceframe']['number'] = -1
#    fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'])
#    assert np.sum(fitstbl.find_frames('trace')) == 2


# TODO: Need a function that checks the calibration groups have the
# correct number of calibration frames
#def test_match_science_errors(fitstbl):
#    par = fitstbl.spectrograph.default_pypeit_par()
#    par['calibrations']['traceframe']['number'] = 10
#    with pytest.raises(PypeItError):
#        fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'])


def test_instr_setup(fitstbl):
    """
    Test instrument setup naming convention
    """
    detname = fitstbl.spectrograph.get_det_name(1)
    calib_key = CalibFrame.construct_calib_key(fitstbl['setup'][0], fitstbl['calib'][0], detname)

    # Check the calibration  key
    assert calib_key == 'A_0_DET01'
    # Invalid detector
    with pytest.raises(PypeItError):
        # Shane kast blue doesn't have a second detector
        fitstbl.spectrograph.get_det_name(2)



# TODO: Need a test that adds a calibration group and checks the result
#    # New calib set
#    #  Turn exposure 9 into an arc
#    fitstbl.edit_frame_type(-1, 'arc')
#    fitstbl['sci_ID'][-1] = 2
#    # Turn off other arc
#    fitstbl['sci_ID'][1] = 1 + 4 + 8
#    # Run
#    setupID3, setup_dict = pypsetup.instr_setup(2, 1, fitstbl, setup_dict=setup_dict)
#    assert setupID3 == 'A_01_ab'
#    assert setup_dict['A']['ab']['arc'][0] == 'b009.fits.gz'


def test_exptime():
    exptime = np.array([0, 30, None, 900])
    assert np.array_equal(framematch.check_frame_exptime(exptime, [0,None]),
                          np.array([False, True, False, True]))
    assert np.array_equal(framematch.check_frame_exptime(exptime, [None,1000]),
                          np.array([True, True, False, True]))
    assert np.array_equal(framematch.check_frame_exptime(exptime, [None,None]),
                          np.array([True, True, False, True]))
    assert np.array_equal(framematch.check_frame_exptime(exptime, [None,500]),
                          np.array([True, True, False, False]))
    assert np.array_equal(framematch.check_frame_exptime(exptime, [10,20]),
                          np.array([False, False, False, False]))


def test_valid_frametype():

    assert framematch.valid_frametype('arc'), 'arc should be a valid frametype'
    assert framematch.valid_frametype('arc', quiet=True), 'arc should be a valid frametype'
    assert not framematch.valid_frametype('junk', quiet=True, raise_error=False), \
            'junk should not be a valid frametype'
    assert not framematch.valid_frametype('junk', quiet=False, raise_error=False), \
            'junk should not be a valid frametype'
    # These should raise errors
    with pytest.raises(PypeItError):
        framematch.valid_frametype('junk', quiet=False, raise_error=True)
    with pytest.raises(PypeItError):
        framematch.valid_frametype('junk', quiet=True, raise_error=True)


