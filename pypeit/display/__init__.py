"""
Register the ginga global plugin(s).
"""
import os.path
from importlib import metadata
import numpy

from ginga.misc.Bunch import Bunch

required_plugins = ['Spec1dView', 'Spec2dView']

def plugins_available(return_report=False):
    available_plugins = []
    # WARNING:
    #   - metadata.entry_points(group='ginga.rv.plugins') doesn't work in python3.9
    #   - and metadata.entry_points()['ginga.rv.plugins'] doesn't work in python3.12
    for entry_point in metadata.entry_points(group='ginga.rv.plugins'):
        spec = entry_point.load()()
        available_plugins += [spec.get('name', spec.get('menu',
                                spec.get('klass', spec.get('module'))))]
    indx = numpy.isin(required_plugins, available_plugins)
    result = (numpy.all(indx),)
    if return_report:
        result += ('' if result[0] else \
                   f'Missing plugins: {numpy.asarray(required_plugins)[numpy.logical_not(indx)]}',)
    return result

def pypeit_plugin_specs():
    """
    Return the specifications of all PypeIt ginga plugins.

    These are added directly to the in-process viewer (see
    :mod:`pypeit.display.viewer`), so the viewer does not depend on the
    entry points registered when PypeIt is installed.

    Returns
    -------
    :obj:`list`
        List of ``ginga.misc.Bunch.Bunch`` plugin specifications.
    """
    return [setup_Spec1dView(), setup_Spec2dView()]

def setup_Spec1dView():
    return Bunch(path=os.path.join(os.path.split(__file__)[0], 'spec1dview.py'),
                 name='Spec1dView', module='spec1dview', klass='Spec1dView',
                 ptype='local', workspace='right', start=False,
                 category='PypeIt', menu='Spec1dView', tab='Spec1dView')

def setup_Spec2dView():
    return Bunch(path=os.path.join(os.path.split(__file__)[0], 'spec2dview.py'),
                 name='Spec2dView', module='spec2dview', klass='Spec2dView',
                 ptype='global', workspace='right', start=False,
                 category='PypeIt', menu='Spec2dView', tab='Spec2dView')
