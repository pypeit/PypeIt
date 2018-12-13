# Module to run tests on sort and arsetup
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy as np

from astropy.table import hstack

from pypeit import metadata
from pypeit.core import framematch
from pypeit.core import pypsetup
from pypeit.pypmsgs import PypeItError


@pytest.fixture
def fitstbl():
    return metadata.dummy_fitstbl()


@pytest.fixture
def fitstblno():
    return metadata.dummy_fitstbl(notype=True)


# TODO: These are out of date (but still pass)
def test_chk_condition(fitstbl):
    # Lamp (str)
    cond = 'lampstat06=on'
    ntmp = framematch.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 1
    # exptime (float)
    cond = 'exptime>30'
    ntmp = framematch.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 6
    cond = 'exptime<30'
    ntmp = framematch.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 1
    cond = 'exptime<=30'
    ntmp = framematch.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 4
    cond = 'exptime!=30'
    ntmp = framematch.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 7


def test_sort_data(fitstbl):
    """ Test sort_data
    """
    # Sort
    assert fitstbl.find_frames('bias')[0]
    assert fitstbl.find_frames('arc')[1]
    assert fitstbl.find_frames('trace')[2]
    assert fitstbl.find_frames('standard')[4]
    assert np.sum(fitstbl.find_frames('science')) == 5


def test_match_science(fitstbl):
    """ Test match_science routine
    """
    calib_ID = 0
    par = fitstbl.spectrograph.default_pypeit_par()
    # Match and test
    fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'],
                             setup=True)
    assert fitstbl.find_frames('arc', calib_ID=calib_ID, index=True)[0] == 1
    assert fitstbl.find_frames('standard', calib_ID=calib_ID, index=True)[0] == 4
    assert fitstbl.find_frames('trace', calib_ID=calib_ID, index=True)[0] == 2
    #assert fitstbl['sci_ID'][0] == 31


def test_neg_match_science(fitstbl):
    """ Test using negative number for calibs
    """
    par = fitstbl.spectrograph.default_pypeit_par()
    # Use negative number
    for ftype in ['arc', 'pixelflat', 'bias']:
        par['calibrations']['{0}frame'.format(ftype)]['number'] = 1
    par['calibrations']['traceframe']['number'] = -1
    fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'])
    assert np.sum(fitstbl.find_frames('trace')) == 2


def test_match_science_errors(fitstbl):
    par = fitstbl.spectrograph.default_pypeit_par()
    par['calibrations']['traceframe']['number'] = 10
    with pytest.raises(PypeItError):
        fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'])


# TODO -- Need a new unit test

'''
def test_instr_setup(fitstbl):
    """ Test instrument setup naming convention
    Tickles most of the arsetup methods
    """
    par = fitstbl.spectrograph.default_pypeit_par()
    fitstbl.match_to_science(par['calibrations'], par['rdx']['calwin'], par['fluxcalib'],
                             setup=True)

    # Get an ID
    namp = 1
    setupID, setup_dict = pypsetup.instr_setup(1, 1, fitstbl)
    assert setupID == 'A_01_aa'
    # Should get same thing
    setupID, setup_dict = pypsetup.instr_setup(1, 1, fitstbl, setup_dict=setup_dict)
    assert setupID == 'A_01_aa'
    with pytest.raises(IndexError):
        # Shane kast blue doesn't have a second detector
        setupID2, setup_dict = pypsetup.instr_setup(1, 2, fitstbl, setup_dict=setup_dict)

    # New calib set
    #  Turn exposure 9 into an arc
    fitstbl.edit_frame_type(-1, 'arc')
    fitstbl['sci_ID'][-1] = 2
    # Turn off other arc
    fitstbl['sci_ID'][1] = 1 + 4 + 8
    # Run
    setupID3, setup_dict = pypsetup.instr_setup(2, 1, fitstbl, setup_dict=setup_dict)
    assert setupID3 == 'A_01_ab'
    assert setup_dict['A']['ab']['arc'][0] == 'b009.fits.gz'
'''


