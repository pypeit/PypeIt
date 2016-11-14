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

    outfil = 'tst.log'
    msgs = pyparm.Messages(outfil, debug, 1)
    msgs.close()
    # Insure scipy, numpy, astropy are being version
    with open(outfil, 'r') as f:
        lines = f.readlines()
    pckgs = ['scipy', 'numpy', 'astropy']
    flgs = [False]*len(pckgs)
    for line in lines:
        for jj,pckg in enumerate(pckgs):
            if pckg in line:
                flgs[jj] = True
    for flg in flgs:
        assert flg is True


def test_msgs():
    from pypit import ardebug
    from pypit import armsgs as pyparm
    debug = ardebug.init()

    msgs = pyparm.Messages(None, debug, 1)
    msgs.info("test 123")
    msgs.warn("test 123")
    msgs.bug("test 123")
    msgs.work("test 123")
    msgs.close()

