# Module to run tests on arsort

import pytest

import numpy as np
from pypit import pyputils
msgs = pyputils.get_dummy_logger(develop=True)
from pypit import arsort
from pypit import arutils
from pypit import arparse as settings

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


@pytest.fixture
def fitsdict():
    return arutils.dummy_fitsdict()


def test_chk_condition(fitsdict):
    # Lamp (str)
    cond = 'lampstat06=on'
    ntmp = arsort.chk_condition(fitsdict, cond)
    assert np.sum(ntmp) == 1
    # exptime (float)
    cond = 'exptime>30'
    ntmp = arsort.chk_condition(fitsdict, cond)
    assert np.sum(ntmp) == 6
    cond = 'exptime<30'
    ntmp = arsort.chk_condition(fitsdict, cond)
    assert np.sum(ntmp) == 1
    cond = 'exptime<=30'
    ntmp = arsort.chk_condition(fitsdict, cond)
    assert np.sum(ntmp) == 4
    cond = 'exptime!=30'
    ntmp = arsort.chk_condition(fitsdict, cond)
    assert np.sum(ntmp) == 7


def test_sort_data(fitsdict):
    """ Test sort_data
    """
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Sort
    filesort = arsort.sort_data(fitsdict)
    assert filesort['bias'][0] == 0
    assert filesort['arc'][0] == 1
    assert filesort['trace'][0] == 2
    assert filesort['standard'][0] == 4
    assert len(filesort['science']) == 5


def test_user_frametype(fitsdict):
    """ Test setting frametype manually
    """
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Modify settings -- WARNING: THIS IS GLOBAL!
    settings.spect['set'] = {}
    settings.spect['set']['standard'] = ['b009.fits']
    filesort = arsort.sort_data(fitsdict)
    assert 9 in filesort['standard']
    settings.spect['set'] = {}


def test_match_science(fitsdict):
    """ Test match_science routine
    """
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    settings.argflag['run']['setup'] = True  # Over-ride default numbers
    filesort = arsort.sort_data(fitsdict)
    # Match and test
    arsort.match_science(fitsdict, filesort)
    assert settings.spect['arc']['index'][1][0] == 1
    assert settings.spect['standard']['index'][1][0] == 4
    assert len(settings.spect['trace']['index'][0]) == 2


def test_neg_match_science(fitsdict):
    """ Test using negative number for calibs
    """
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    filesort = arsort.sort_data(fitsdict)
    # Use negative number
    for ftype in ['arc', 'pixelflat', 'trace', 'bias']:
        settings.spect[ftype]['number'] = 1
    settings.spect['trace']['number'] = -1
    arsort.match_science(fitsdict, filesort)
    assert len(settings.spect['trace']['index'][1]) == 2


def test_instr_setup(fitsdict):
    """ Test instrument setup naming convention
    """
    from pypit import arsciexp
    arutils.dummy_settings(spectrograph='shane_kast_blue', set_idx=False)
    # Load
    settings.argflag['run']['setup'] = True # Over-ride default numbers
    filesort = arsort.sort_data(fitsdict)
    # Match and test
    arsort.match_science(fitsdict, filesort)
    # Make a science frame
    sciexp = arsciexp.ScienceExposure(0, fitsdict, do_qa=False)
    # Get an ID
    setup_dict = {}
    setupID = arsort.instr_setup(sciexp, 1, fitsdict, setup_dict)
    assert setupID == 'A_01_aa'
    # Should get same thing
    setupID = arsort.instr_setup(sciexp, 1, fitsdict, setup_dict)
    assert setupID == 'A_01_aa'
    # New det (fake out kast_blue)
    settings.spect['det02'] = dict(numamplifiers=1)
    setupID2 = arsort.instr_setup(sciexp, 2, fitsdict, setup_dict)
    assert setupID2 == 'A_02_aa'
    # New calib set
    settings.spect['arc']['index'][1] = np.array([9])  # Not really an arc, but ok
    sciexp1 = arsciexp.ScienceExposure(1, fitsdict, do_qa=False)
    setupID3 = arsort.instr_setup(sciexp1, 1, fitsdict, setup_dict)
    assert setupID3 == 'A_01_ab'
    assert setup_dict['A']['ab']['arcs'][0] == 'b009.fits'
