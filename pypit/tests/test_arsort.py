# Module to run tests on arsort and arsetup
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

import numpy as np

from astropy.table import hstack

from pypit import arparse as settings
from pypit.core import arsort
from pypit.core import arsetup
from pypit.armsgs import PypitError


@pytest.fixture
def fitstbl():
    return arsort.dummy_fitstbl()

@pytest.fixture
def fitstblno():
    return arsort.dummy_fitstbl(notype=True)


def test_chk_condition(fitstbl):
    # Lamp (str)
    cond = 'lampstat06=on'
    ntmp = arsort.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 1
    # exptime (float)
    cond = 'exptime>30'
    ntmp = arsort.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 6
    cond = 'exptime<30'
    ntmp = arsort.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 1
    cond = 'exptime<=30'
    ntmp = arsort.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 4
    cond = 'exptime!=30'
    ntmp = arsort.chk_condition(fitstbl, cond)
    assert np.sum(ntmp) == 7


def test_sort_data(fitstbl):
    """ Test sort_data
    """
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Sort
    filesort = arsort.type_data(fitstbl, settings.spect, settings.argflag)
    assert filesort['bias'][0]
    assert filesort['arc'][1]
    assert filesort['trace'][2]
    assert filesort['standard'][4]
    assert np.sum(filesort['science']) == 5


def test_user_frametype(fitstbl):
    """ Test setting frametype manually
    """
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Modify settings -- WARNING: THIS IS GLOBAL!
    settings.spect['set'] = {}
    settings.spect['set']['standard'] = ['b009.fits.gz']
    filesort = arsort.type_data(fitstbl, settings.spect, settings.argflag)
    assert filesort['standard'][9]
    settings.spect['set'] = {}

def test_match_science(fitstblno):
    """ Test match_science routine
    """
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    settings.argflag['run']['setup'] = True  # Over-ride default numbers
    filesort = arsort.type_data(fitstblno, settings.spect, settings.argflag)
    # Match and test
    mtbl = hstack([fitstblno,filesort])
    fitstbl = arsort.match_to_science(mtbl, settings.spect, settings.argflag)
    assert arsort.ftype_indices(fitstbl, 'arc', 1)[0] == 1
    assert arsort.ftype_indices(fitstbl, 'standard', 4)[0] == 4
    assert arsort.ftype_indices(fitstbl, 'trace', 1)[0] == 2
    assert fitstbl['sci_ID'][0] == 31

def test_neg_match_science(fitstblno):
    """ Test using negative number for calibs
    """
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    filesort = arsort.type_data(fitstblno, settings.spect, settings.argflag)
    mtbl = hstack([fitstblno,filesort])
    # Use negative number
    for ftype in ['arc', 'pixelflat', 'trace', 'bias']:
        settings.spect[ftype]['number'] = 1
    settings.spect['trace']['number'] = -1
    fitstbl = arsort.match_to_science(mtbl, settings.spect, settings.argflag)
    assert np.sum(fitstbl['trace']) == 2


def test_match_science_errors(fitstblno):
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    filesort = arsort.type_data(fitstblno, settings.spect, settings.argflag)
    mtbl = hstack([fitstblno,filesort])
    settings.spect['trace']['number'] = 10
    with pytest.raises(PypitError):
        _ = arsort.match_to_science(mtbl, settings.spect, settings.argflag)


def test_instr_setup(fitstblno):
    """ Test instrument setup naming convention
    Tickles most of the arsetup methods
    """
    settings.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    settings.argflag['run']['setup'] = True # Over-ride default numbers
    filesort = arsort.type_data(fitstblno, settings.spect, settings.argflag)
    mtbl = hstack([fitstblno,filesort])
    # Match and test
    fitstbl = arsort.match_to_science(mtbl, settings.spect, settings.argflag)

    # Make a science frame (old)
    #sciexp = arsciexp.ScienceExposure(1, fitstbl, settings.argflag, settings.spect, do_qa=False)
    # Get an ID
    namp = 1
    setup_dict = {}
    setupID = arsetup.instr_setup(1, 1, fitstbl, setup_dict, namp)
    assert setupID == 'A_01_aa'
    # Should get same thing
    setupID = arsetup.instr_setup(1, 1, fitstbl, setup_dict, namp)
    assert setupID == 'A_01_aa'
    # New det (fake out kast_blue)
    settings.spect['det02'] = dict(numamplifiers=1)
    #setupID2 = arsetup.instr_setup(sciexp, 2, fitsdict, setup_dict)
    setupID2 = arsetup.instr_setup(1, 2, fitstbl, setup_dict, namp)
    assert setupID2 == 'A_02_aa'

    # New calib set
    #  Turn exposure 9 into an arc
    fitstbl['arc'][-1] = True
    fitstbl['science'][-1] = False
    fitstbl['sci_ID'][-1] = 2
    # Turn off other arc
    fitstbl['sci_ID'][1] = 1 + 4 + 8
    # Run
    setupID3 = arsetup.instr_setup(2, 1, fitstbl, setup_dict, namp)
    assert setupID3 == 'A_01_ab'
    assert setup_dict['A']['ab']['arc'][0] == 'b009.fits.gz'
