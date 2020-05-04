"""
Module to run tests on skysub routines
"""
import pytest
import numpy as np

from pypeit.core import skysub
from pypeit.slittrace import SlitTraceSet


def test_userregions():
    """ Run the parameter setup script
    """
    regions = [":10",
               "30:60",
               "80:",
               ":10,35:65,80:",
               "10:20;",
               "10:20,"]
    result = [[[0, 100]],
              [[300, 599]],
              [[799, 1000]],
              [[0, 100], [350, 649], [799, 1000]],
              [],
              [[100, 200]]
              ]
    resstat = [0, 0, 0, 0, 1, 2]
    for rr, reg in enumerate(regions):
        status, regs = skysub.read_userregions(reg, resolution=1000)
        assert(regs == result[rr])
        assert(status == resstat[rr])


def test_generatemask():
    resolution = 1000
    reg = "80:"
    tstmsk = np.zeros((1000, 1000))
    tstmsk[:, 744:901] = 1
    status, regs = skysub.read_userregions(reg, resolution=resolution)
    slits = SlitTraceSet(left_init=np.full((1000, 1), 120, dtype=float),
                         right_init=np.full((1000, 1), 900, dtype=float), binspec=1, binspat=1,
                         pypeline='IFU', nspat=1000, PYP_SPEC='dummy')
    skymask = skysub.generate_mask("IFU", regs, slits, slits.left_init, slits.right_init, resolution=resolution)
    assert(np.array_equal(skymask, tstmsk))
