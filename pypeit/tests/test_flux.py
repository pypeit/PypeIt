"""
Module to run tests on a few flux routines
"""

import numpy as np

from pypeit.core import flux_calib
from pypeit.par.pypeitpar import Coadd1DPar


def test_filter_scale():
    # Test scale_in_filter() method which is called in coadding
    wave = np.arange(3000.,10000.)
    flux = np.ones_like(wave)
    gdm = np.ones_like(wave, dtype=bool)
    #
    par = Coadd1DPar()
    par['filter'] = 'DECAM-R'
    par['filter_mag'] = 17.
    # Run
    scale = flux_calib.scale_in_filter(wave, flux, gdm, par)
    assert np.isclose(scale, 41.698475048180406, rtol=1e-3)
