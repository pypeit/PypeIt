"""
Module to run tests on core/meta.py
"""
import pytest

from pypeit.core import meta

# TODO -- We could use *many* more tests on core/meta.py

def test_meta():
    # Constructr the simple dicg
    d = meta.get_meta_data_model()
    assert 'mjd' in d.keys(), 'MJD not in keys'
    assert isinstance(d, dict)


