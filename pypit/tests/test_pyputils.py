# Module to run tests on simple fitting routines for arrays
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# TEST_UNICODE_LITERALS

try:
    basestring
except NameError:
    basestring = str

from pypit import __version__, __last_updated__
from pypit import pyputils

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)

def test_version():
    # Dummy self
#    ver,upd = pyputils.get_version()
    ver,upd = __version__, __last_updated__
    assert isinstance(ver,basestring)


def test_dummy_logger():
    msgs = pyputils.get_dummy_logger()
    msgs.info('testing 123')

