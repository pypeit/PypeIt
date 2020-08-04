"""
Module to run tests on arcoadd
"""
import os
import numpy

import pytest

from IPython import embed
from pkg_resources import iter_entry_points

from pypeit.display import required_plugins

def test_plugin():
    available_plugins = []
    for entry_point in iter_entry_points(group='ginga.rv.plugins', name=None):
        spec = entry_point.load()()
        available_plugins += [spec.get('name', spec.get('menu',
                                                        spec.get('klass', spec.get('module'))))]
    indx = numpy.isin(required_plugins, available_plugins)
    assert numpy.all(indx), \
            'Missing plugins: {0}'.format(required_plugins[numpy.logical_not(indx)])
