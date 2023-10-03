"""
Module to run tests on arcoadd
"""
from IPython import embed

from pypeit.display import plugins_available

def test_plugin():
    success, report = plugins_available(return_report=True)
    assert success, report

