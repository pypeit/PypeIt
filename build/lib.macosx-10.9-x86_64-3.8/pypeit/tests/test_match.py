"""
Module to run tests on sort and arsetup
"""
import pytest

import numpy as np

from pypeit.core import framematch
from pypeit.tests.tstutils import dummy_fitstbl


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
    """ Test instrument setup naming convention
    Tickles most of the arsetup methods
    """
    par = fitstbl.spectrograph.default_pypeit_par()

    # Check the master key
    assert fitstbl.master_key(0) == 'A_1_01'
    # Invalid detector
    with pytest.raises(IndexError):
        # Shane kast blue doesn't have a second detector
        fitstbl.master_key(0, det=2)

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


