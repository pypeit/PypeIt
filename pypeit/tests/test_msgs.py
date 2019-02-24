"""
Module to run tests on armsgs
"""
import numpy as np
import pytest

from pypeit import pypmsgs

def test_log_write():

    outfil = 'tst.log'
    msgs = pypmsgs.Messages(outfil, verbosity=1)
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
    msgs = pypmsgs.Messages(None, verbosity=1)
    msgs.info("test 123")
    msgs.warn("test 123")
    msgs.bug("test 123")
    msgs.work("test 123")
    msgs.close()

