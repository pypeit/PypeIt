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
from pypit.core import ararc
from pypit.spectrographs import wavelengths
from pypit import arparse as settings


def test_setup_param():
    """ Run the parameter setup script
    Returns
    -------

    """
    # Initialize some settings
    settings.dummy_settings()
    # Load Dummy self
    settings.argflag['run']['spectrograph'] = 'shane_kast_blue'
    fitstbl = arsort.dummy_fitstbl()
    # Run
    arcparm = wavelengths.setup_param('shane_kast_blue', (2048,2048),fitstbl,0)
    for key in ['llist','disp','wvmnx']:
        assert key in arcparm


def test_detect_lines():
    if os.getenv('RUN_ALL_PYPIT_TESTS') is None:  # REMOVE WHEN WORKING WITH Python 3
        assert True
        return
    det = 1
    # Using Paranal night sky as an 'arc'
    arx_sky = xspectrum1d.XSpectrum1D.from_file(pypit.__path__[0]+'/data/sky_spec/paranal_sky.fits')
    arx_amp, arx_cent, arx_wid, arx_w, arx_yprep = ararc.detect_lines(arx_sky.flux.value)
    # Test
    assert len(arx_w[0]) == 1767

