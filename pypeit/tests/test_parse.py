"""
Module to run tests on arparse
"""
from IPython import embed

import pytest

import numpy as np

from pypeit.pypmsgs import PypeItError
from pypeit.core import parse
from pypeit.spectrographs.util import load_spectrograph

def test_parse_binning():
    """ Test parse binning algorithm
    """
    bin1, bin2 = parse.parse_binning('2,2')
    assert bin1 == 2
    assert bin2 == 2
    # Other input
    bin1, bin2 = parse.parse_binning((2,2))
    assert bin1 == 2
    assert bin2 == 2


def test_sec2slice():
    sub = ':10,10:'
    subslice = parse.sec2slice(sub, require_dim=2)
    assert subslice[0].start is None

    subslice = parse.sec2slice(sub, include_end=True, require_dim=2)
    assert subslice[0].stop == 11

    subslice = parse.sec2slice(sub, one_indexed=True, require_dim=2)
    assert subslice[0].stop == 9


def test_str2list():
    """
    Test the conversion of the string to a list of integers
    """
    assert np.array_equal(parse.str2list('all', length=10), [0,1,2,3,4,5,6,7,8,9])
    assert np.array_equal(parse.str2list(':4', length=10), [0,1,2,3])
    assert np.array_equal(parse.str2list('3:5,8:', length=10), [3,4,8,9])
    assert np.array_equal(parse.str2list('3,1:5,6', length=10), [1,2,3,4,6])
    assert np.array_equal(parse.str2list('3,1:5,8:', length=10), [1,2,3,4,8,9])


def test_parse_slitspatnum():
    assert [x[0] for x in parse.parse_slitspatnum('DET01:224')] == ['DET01', 224], \
        'Bad parsing for single pair'
    assert [x[0] for x in parse.parse_slitspatnum(['DET01:224'])]  == ['DET01', 224], \
        'Bad parsing for single list pair'
    assert [x.tolist() for x in parse.parse_slitspatnum('DET01:224,DET02:331')] \
        == [['DET01', 'DET02'], [224, 331]], 'Bad parsing of comma-separated pairs'
    assert [x.tolist() for x in parse.parse_slitspatnum(['DET01:224,DET02:331'])] \
        == [['DET01', 'DET02'], [224, 331]], 'Bad parsing of comma-separated pairs list'
    assert [x.tolist() for x in parse.parse_slitspatnum(['DET01:224', 'DET02:331'])] \
        == [['DET01', 'DET02'], [224, 331]], 'Bad parsing of list of pairs'
    assert [x.tolist() for x in parse.parse_slitspatnum(['DET01:224,DET02:331', 'DET03:442'])] \
        == [['DET01', 'DET02', 'DET03'], [224, 331, 442]], 'Bad mixed parsing'


def test_parse_image_location():
    # keck_nires just used as an example of a single detector spectrograph
    spec = load_spectrograph('keck_nires')

    p = parse.parse_image_location('-1:34.5:400.1:4', spec)
    assert len(p) == 5, 'Incorrect number of returned objects'
    assert p[0], 'Should return that the detector integer is negative'
    assert p[1] == 'DET01', 'Wrong detector identifier'

    p = parse.parse_image_location('1:34.5:400.1:4', spec)
    assert len(p) == 5, 'Incorrect number of returned objects'
    assert not p[0], 'Should return that the detector integer is positive'
    assert p[1] == 'DET01', 'Wrong detector identifier'

    # This should fail because keck_nires does not have this mosaic
    with pytest.raises(PypeItError):
        p = parse.parse_image_location('(1,2,3):34.5:400.1:4', spec)

    spec = load_spectrograph('gemini_gmos_south_ham')
    p = parse.parse_image_location('(1,2,3):34.5:400.1:4', spec)
    assert len(p) == 5, 'Incorrect number of returned objects'
    assert not p[0], 'Should return that the detector integer is positive'
    assert p[1] == 'MSC01', 'Wrong mosaic identifier'

    p = parse.parse_image_location('(-1,-2,-3):34.5:400.1:4', spec)
    assert len(p) == 5, 'Incorrect number of returned objects'
    assert p[0], 'Should return that the detector integer is negative'
    assert p[1] == 'MSC01', 'Wrong mosaic identifier'

    p = parse.parse_image_location('2:34.5:400.1', spec)
    assert len(p) == 4, 'Incorrect number of returned objects'
    assert not p[0], 'Should return that the detector integer is positive'
    assert p[1] == 'DET02', 'Wrong mosaic identifier'


def test_join_image_location():
    # Fix without parentheses
    par = ['3:1500:331; 3:1500:635']
    list_fixed_par = parse.fix_config_par_image_location(par)

    par = '3:1500:331; 3:1500:635'
    fixed_par = parse.fix_config_par_image_location(par)
    assert fixed_par[0] == par.split(';')[0], 'Incorrect parsing of the 1st element'
    assert fixed_par[1] == par.split(';')[1].strip(), 'Incorrect parsing of the 2nd element'
    assert list_fixed_par == fixed_par, 'str and one-item list entry should yield identical result'

    # Fix with parentheses
    par = ['(1', '2', '3):1500:331; (1', '2', '3):1500:635']
    fixed_par = parse.fix_config_par_image_location(par)
    assert fixed_par[0] == (','.join(par[:3])).split(';')[0], '1st element fix is wrong'
    assert fixed_par[1] == (','.join(par[2:])).split(';')[1].strip(), '2nd element fix is wrong'

    # Mix mosaics and single detectors
    par = ['(1', '2', '3):1500:331; 3:1500:635']
    fixed_par = parse.fix_config_par_image_location(par)
    assert fixed_par[0] == (','.join(par[:3])).split(';')[0], '1st element fix is wrong'
    assert fixed_par[1] == (','.join(par[2:])).split(';')[1].strip(), '2nd element fix is wrong'

