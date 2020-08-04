"""
Register the ginga global plugin(s).
"""
import os.path

from ginga.misc.Bunch import Bunch

required_plugins = ['SlitWavelength']

def setup_SlitWavelength():
    return Bunch(path=os.path.join(os.path.split(__file__)[0], 'ginga_plugins.py'),
                 module='ginga_plugins', klass='SlitWavelength',
                 ptype='global', workspace='right', start=False,
                 category='PypeIt', menu='SlitWavelength', tab='SlitWavelength')
