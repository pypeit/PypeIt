"""
Module to run tests on arcoadd
"""
import os
import numpy

import pytest

from IPython import embed
from pkg_resources import iter_entry_points

from pypeit.display import plugins_available

def test_plugin():
    success, report = plugins_available(return_report=True)
    assert success, report
