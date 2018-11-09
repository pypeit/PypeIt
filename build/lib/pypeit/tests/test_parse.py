# Module to run tests on arparse
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import pytest

from pypeit.core import parse

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_parse_binning():
    """ Test parse binning algorithm
    """
    bin1, bin2 = parse.parse_binning('2,2')
    assert bin1 == 2
    assert bin2 == 2
    # Other input
    bin1, bin2 = parse.parse_binning((2,2))   # String output required so this returns 1,1 (the default)
    assert bin1 == 1
    assert bin2 == 1

def test_sec2slice():
    sub = ':10,10:'
    subslice = parse.sec2slice(sub, require_dim=2)
    assert subslice[0].start is None

    subslice = parse.sec2slice(sub, include_end=True, require_dim=2)
    assert subslice[0].stop == 11

    subslice = parse.sec2slice(sub, one_indexed=True, require_dim=2)
    assert subslice[0].stop == 9

    
