# Module to run tests on arparse

import pytest


from pypit import pyputils
msgs = pyputils.get_dummy_logger()
from pypit import arparse

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_armlsd():
    """ Test loading of ARMLSD
    """
    argf = arparse.get_argflag_class(('ARMLSD', ''))
    argf.init_param()
    argf.set_param('run pypitdir {0:s}'.format('pypit'))
    argf.set_param('run redname {0:s}'.format('kast_blue_setup_01'))
    # Test
    assert argf._argflag['run']['redname'] == 'kast_blue_setup_01'
    assert argf._argflag['reduce']['trim'] == True
    assert argf._argflag['reduce']['skysub']['perform'] == True
    # Save
    argf.save()


def test_load_spectrograph():
    """ Test loading of spectrograph
    Using kast_blue as a fiducial example
    """
    spect = arparse.get_spect_class(('ARMLSD', 'kast_blue', 'kast_blue_setup_01'))
    lines = spect.load_file()
    spect.set_paramlist(lines)
    # Test
    assert spect._spect['keyword']['dispname'] == '01.GRISM_N'
    assert spect._spect['keyword']['dichroic'] == '01.BSPLIT_N'
    spect.save()

def test_parse_binning():
    """ Test parse binning algorithm
    """
    bin1, bin2 = arparse.parse_binning('2,2')
    assert bin1 == 2
    assert bin2 == 2
    # Other input
    bin1, bin2 = arparse.parse_binning((2,2))
    assert bin1 == 1
    assert bin2 == 1
