"""
Module to run tests on sort and arsetup
"""
import pytest

import numpy as np

from pypeit.core import framematch
#from pypeit.pypmsgs import PypeItError

from pypeit.tests.tstutils import dummy_fitstbl


@pytest.fixture
def fitstbl():
    return dummy_fitstbl()


#@pytest.fixture
#def fitstblno():
#    return metadata.dummy_fitstbl(notype=True)


# TODO: These are out of date; removing them for now
#def test_chk_condition(fitstbl):
#    # Lamp (str)
#    cond = 'lampstat06=on'
#    ntmp = framematch.chk_condition(fitstbl, cond)
#    assert np.sum(ntmp) == 1
#    # exptime (float)
#    cond = 'exptime>30'
#    ntmp = framematch.chk_condition(fitstbl, cond)
#    assert np.sum(ntmp) == 6
#    cond = 'exptime<30'
#    ntmp = framematch.chk_condition(fitstbl, cond)
#    assert np.sum(ntmp) == 1
#    cond = 'exptime<=30'
#    ntmp = framematch.chk_condition(fitstbl, cond)
#    assert np.sum(ntmp) == 4
#    cond = 'exptime!=30'
#    ntmp = framematch.chk_condition(fitstbl, cond)
#    assert np.sum(ntmp) == 7


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
#    # Match and test
#    fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'],
#                             setup=True)
    assert fitstbl.find_frames('arc', calib_ID=calib_ID, index=True)[0] == 1
    assert fitstbl.find_frames('standard', calib_ID=calib_ID, index=True)[0] == 4
    assert fitstbl.find_frames('trace', calib_ID=calib_ID, index=True)[0] == 2
    #assert fitstbl['sci_ID'][0] == 31


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


