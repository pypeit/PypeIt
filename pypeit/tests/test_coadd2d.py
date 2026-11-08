"""
Module to test functions in pypeit.coadd2d
"""
import numpy as np
import pytest

from pypeit import PypeItError
from pypeit.coadd2d import CoAdd2D


@pytest.fixture
def coadd():
    # check_input only depends on the number of exposures, so skip __init__, which
    # needs real spec2d files, and set nexp directly
    _coadd = CoAdd2D.__new__(CoAdd2D)
    _coadd.nexp = 2
    return _coadd


def test_check_input_list(coadd):
    weights = coadd.check_input([0.3, np.float32(0.7)], 'weights')
    assert isinstance(weights, list), f'weights should be returned as a list, got {type(weights)}'
    assert np.allclose(weights, [0.3, 0.7]), f'weights values changed: {weights}'
    offsets = coadd.check_input(np.array([0, -5]), 'offsets')
    assert isinstance(offsets, np.ndarray), f'offsets should be returned as an ndarray, got {type(offsets)}'
    assert np.array_equal(offsets, [0, -5]), f'offsets values changed: {offsets}'

    # Wrong number of elements
    with pytest.raises(PypeItError, match='same number of elements'):
        coadd.check_input([1.0, 2.0, 3.0], 'offsets')

    # Lists with the right length but non-numerical (or boolean) elements are rejected
    for bad in [['[1.0', '2.0'], [True, False], [1.0, None]]:
        with pytest.raises(PypeItError, match='Unrecognized format'):
            coadd.check_input(bad, 'weights')

    # Unknown type
    with pytest.raises(PypeItError, match='Unrecognized type'):
        coadd.check_input([1.0, 2.0], 'shifts')
