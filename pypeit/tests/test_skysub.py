"""
Module to run tests on skysub routines (mainly for IFU)
"""
from pathlib import Path
from IPython import embed

import numpy as np

from pypeit.core import skysub
from pypeit.images.buildimage import SkyRegions
from pypeit.slittrace import SlitTraceSet
from pypeit.tests.tstutils import data_path


def test_userregions():
    """ Run the parameter setup script
    """
    regions = [":10",
               "30:60",
               "80:",
               ":10,35:65,80:",
               "10:20;",
               "10:20,"]
    result = [[np.array([], dtype=int), np.array([99], dtype=int)],
              [np.array([299], dtype=int), np.array([598], dtype=int)],
              [np.array([798], dtype=int), np.array([], dtype=int)],
              [np.array([349, 798], dtype=int), np.array([99, 648], dtype=int)],
              [],
              [np.array([99], dtype=int), np.array([199], dtype=int)]
              ]
    resstat = [0, 0, 0, 0, 1, 2]
    nslits = 2
    maxsl = 100
    for rr, reg in enumerate(regions):
        status, regs = skysub.read_userregions(reg, nslits, maxslitlength=maxsl)
        if status != 1:
            assert len(regs) == nslits, 'Bad number of slits'
            assert np.array_equal(np.where(regs[0][1:] & ~regs[0][:-1])[0], result[rr][0]), \
                    'Bad result for slit 1'
            assert np.array_equal(np.where(~regs[0][1:] & regs[0][:-1])[0], result[rr][1]), \
                    'Bad result for slit 2'
        assert status == resstat[rr], 'Unexpected status'


def test_generatemask():
    maxsl = 1000
    nslits = 2
    reg = "80:"
    tstmsk = np.zeros((1000, 1000), dtype=bool)
    tstmsk[:, 744:901] = True
    status, regs = skysub.read_userregions(reg, nslits, maxslitlength=maxsl)
    assert status == 0, 'Bad generation of user regions'
    slits = SlitTraceSet(left_init=np.full((1000, 1), 120, dtype=float),
                         right_init=np.full((1000, 1), 900, dtype=float), binspec=1, binspat=1,
                         pypeline='IFU', nspat=1000, PYP_SPEC='dummy')
    skymask = skysub.generate_mask("IFU", regs, slits, slits.left_init, slits.right_init)
    assert np.array_equal(skymask, tstmsk), 'Sky mask mismatch'


def test_skyregions_io():
    tstmsk = np.zeros((1000, 1000), dtype=bool)
    tstmsk[:, 744:901] = True
    ofile = Path(data_path('test_skyregions.fits')).resolve()
    if ofile.exists():
        ofile.unlink()
    SkyRegions(image=tstmsk.astype(float), PYP_SPEC='dummy').to_file(file_path=ofile)
    assert ofile.exists(), 'SkyRegions file not written'
    reg = SkyRegions.from_file(ofile)
    assert np.array_equal(tstmsk, reg.image.astype(bool)), 'Bad read'

    # Clean-up
    ofile.unlink()


