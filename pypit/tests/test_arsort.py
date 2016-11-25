# Module to run tests on arsort

import pytest

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

def test_sort_data(fitsdict):
    """ Test sort_data
    """
    arutils.dummy_settings(spectrograph='kast_blue')
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
    arutils.dummy_settings(spectrograph='kast_blue')
    # Modify settings -- WARNING: THIS IS GLOBAL!
    settings.spect['set'] = {}
    settings.spect['set']['standard'] = ['b009.fits']
    filesort = arsort.sort_data(fitsdict)
    assert 9 in filesort['standard']
    settings.spect['set'] = {}


def test_match_science(fitsdict):
    """ Test match_science routine
    """
    arutils.dummy_settings(spectrograph='kast_blue')
    # Load
    settings.argflag['run']['setup'] = True # Over-ride default numbers
    filesort = arsort.sort_data(fitsdict)
    # Match and test
    arsort.match_science(fitsdict, filesort)
    assert settings.spect['arc']['index'][1][0] == 1
    assert settings.spect['standard']['index'][1][0] == 4
    assert len(settings.spect['trace']['index']) == 6


def test_neg_match_science(fitsdict):
    """ Test using negative number for calibs
    """
    arutils.dummy_settings(spectrograph='kast_blue')
    # Load
    filesort = arsort.sort_data(fitsdict)
    # Use negative number
    for ftype in ['arc', 'pixelflat', 'trace', 'slitflat', 'bias']:
        settings.spect[ftype]['number'] = 1
    settings.spect['trace']['number'] = -1
    arsort.match_science(fitsdict, filesort)
    assert len(settings.spect['trace']['index'][1]) == 2


