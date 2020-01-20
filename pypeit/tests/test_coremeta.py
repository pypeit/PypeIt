"""
Module to run tests on arcoadd
"""
import pytest

from pypeit.core import meta

def test_meta():
    d = meta.get_meta_data_model()
    assert 'mjd' in d.keys(), 'MJD not in keys'


