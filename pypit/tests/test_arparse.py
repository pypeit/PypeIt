# Module to run tests on arparse
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

from pypit import pyputils
from pypit import arparse

msgs = pyputils.get_dummy_logger()

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_load_arms():
    """ Test loading of ARMS
    """
    argf = arparse.get_argflag_class(('ARMS', ''))
    argf.init_param()
    argf.set_param('run pypitdir {0:s}'.format('pypit'))
    argf.set_param('run redname {0:s}'.format('shane_kast_blue_setup_01'))
    # Test
    assert argf._argflag['run']['redname'] == 'shane_kast_blue_setup_01'
    assert argf._argflag['reduce']['trim'] == True
    assert argf._argflag['reduce']['skysub']['perform'] == True
    # Save
    argf.save()


def test_load_spectrograph():
    """ Test loading of spectrograph
    Using kast_blue as a fiducial example
    """
    spect = arparse.get_spect_class(('ARMS', 'shane_kast_blue', 'shane_kast_blue_setup_01'))
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
    bin1, bin2 = arparse.parse_binning((2,2))   # String output required so this returns 1,1 (the default)
    assert bin1 == 1
    assert bin2 == 1

def test_sec2slice():
    sub = ':10,10:'
    subslice = arparse.sec2slice(sub, require_dim=2)
    assert subslice[0].start is None

    subslice = arparse.sec2slice(sub, include_end=True, require_dim=2)
    assert subslice[0].stop == 11

    subslice = arparse.sec2slice(sub, one_indexed=True, require_dim=2)
    assert subslice[0].stop == 9

    
