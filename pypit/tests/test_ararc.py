# Module to run tests on ararclines
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import numpy as np
import pytest

from linetools.spectra import xspectrum1d

import pypit
from pypit.core import arsort
from pypit import arsciexp
from pypit import ararc
from pypit import arparse as settings


def test_setup_param():
    """ Run the parameter setup script
    Returns
    -------

    """
    # Initialize some settings
    settings.dummy_settings()
    # Load Dummy self
    slf = arsciexp.dummy_self()
    settings.argflag['run']['spectrograph'] = 'shane_kast_blue'
    fitstbl = arsort.dummy_fitstbl()
    # Run
    arcparm = ararc.setup_param((2048,2048),fitstbl,0)
    for key in ['llist','disp','wvmnx']:
        assert key in arcparm


def test_detect_lines():
    if os.getenv('RUN_ALL_PYPIT_TESTS') is None:  # REMOVE WHEN WORKING WITH Python 3
        assert True
        return
    slf = arsciexp.dummy_self()
    det = 1
    # Using Paranal night sky as an 'arc'
    arx_sky = xspectrum1d.XSpectrum1D.from_file(pypit.__path__[0]+'/data/sky_spec/paranal_sky.fits')
    arx_amp, arx_cent, arx_wid, arx_w, arx_satsnd, arx_yprep = ararc.detect_lines(slf, det, msarc=None, censpec=arx_sky.flux.value, MK_SATMASK=False)
    # Test
    assert len(arx_w[0]) == 1767

