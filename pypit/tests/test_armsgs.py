# Module to run tests on armsgs


import numpy as np
import pytest

from astropy import units as u

from pypit import pyputils

#def data_path(filename):
#    data_dir = os.path.join(os.path.dirname(__file__), 'files')
#    return os.path.join(data_dir, filename)


def test_log_write():
    from pypit import ardebug
    from pypit import armsgs as pyparm
    debug = ardebug.init()

    version, last_updated = pyputils.get_version()
    msgs = pyparm.Messages('tst.log', debug, last_updated, version, 1)
    # Insure scipy, numpy, astropy are being version
    pytest.set_trace()

def test_msgs():
    from pypit import ardebug
    from pypit import armsgs as pyparm
    debug = ardebug.init()

    version, last_updated = pyputils.get_version()
    msgs = pyparm.Messages(None, debug, last_updated, version, 1)
    pytest.set_trace()

