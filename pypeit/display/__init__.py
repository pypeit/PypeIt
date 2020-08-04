"""
Register the ginga global plugin(s).
"""
import os.path
from pkg_resources import iter_entry_points
import numpy

from ginga.misc.Bunch import Bunch

required_plugins = ['SlitWavelength']

def plugins_available(return_report=False):
    available_plugins = []
    for entry_point in iter_entry_points(group='ginga.rv.plugins', name=None):
        spec = entry_point.load()()
        available_plugins += [spec.get('name', spec.get('menu',
                                                        spec.get('klass', spec.get('module'))))]
    indx = numpy.isin(required_plugins, available_plugins)
    result = (numpy.all(indx),)
    if return_report:
        result += ('' if result[0] else \
                    'Missing plugins: {0}'.format(required_plugins[numpy.logical_not(indx)]),)
    return result

def setup_SlitWavelength():
    return Bunch(path=os.path.join(os.path.split(__file__)[0], 'ginga_plugins.py'),
                 module='ginga_plugins', klass='SlitWavelength',
                 ptype='global', workspace='right', start=False,
                 category='PypeIt', menu='SlitWavelength', tab='SlitWavelength')
