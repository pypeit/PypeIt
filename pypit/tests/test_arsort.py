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


def test_sort_data():
    """ Test sort_data
    """
    # Load settings
    arutils.dummy_settings(spectrograph='kast_blue')
    # Generate fitsdict
    fitsdict = arutils.dummy_fitsdict()
    # Sort
    filesort = arsort.sort_data(fitsdict)
    assert filesort['bias'][0] == 0
    assert filesort['arc'][0] == 1
    assert filesort['trace'][0] == 2
    assert filesort['standard'][0] == 4
    assert len(filesort['science']) == 5


def test_match_science():
    """ Test match_science routine
    """
    # Load
    arutils.dummy_settings(spectrograph='kast_blue')
    settings.argflag['run']['setup'] = True # Over-ride default numbers
    fitsdict = arutils.dummy_fitsdict()
    filesort = arsort.sort_data(fitsdict)
    # Match and test
    arsort.match_science(fitsdict, filesort)
    assert settings.spect['arc']['index'][1][0] == 1
    assert settings.spect['standard']['index'][1][0] == 4
    assert len(settings.spect['trace']['index']) == 6

def test_neg_match_science():
    """ Test using negative number for calibs
    """
    # Load
    arutils.dummy_settings(spectrograph='kast_blue')
    settings.argflag['run']['setup'] = True # Over-ride default numbers
    fitsdict = arutils.dummy_fitsdict()
    filesort = arsort.sort_data(fitsdict)
    # Use negative number
    for ftype in ['arc', 'pixelflat', 'trace', 'slitflat']:
        settings.spect[ftype]['number'] = 1
    settings.spect['trace']['number'] = -1
    arsort.match_science(fitsdict, filesort)
    assert len(settings.spect['trace']['index'][1]) == 2
