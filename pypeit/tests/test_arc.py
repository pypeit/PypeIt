# Module to run tests on ararclines
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import numpy as np
import pytest

from linetools.spectra import xspectrum1d

import pypeit
from pypeit.core import fsort
from pypeit.core import arc

from pypeit.spectrographs.util import load_spectrograph


def test_setup_param():
    """ Run the parameter setup script
    
    Returns
    -------
    """
    spectrograph = load_spectrograph(spectrograph='shane_kast_blue')
    fitstbl = fsort.dummy_fitstbl()
    # Run
    arcparm = arc.setup_param(spectrograph, (2048,2048), fitstbl, 0)
    for key in ['llist','disp','wvmnx']:
        assert key in arcparm


# TODO: Why is this not run for python 3?
'''
def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    arx_sky = xspectrum1d.XSpectrum1D.from_file(pypeit.__path__[0]+'/data/sky_spec/paranal_sky.fits')
    arx_amp, arx_cent, arx_wid, arx_w, arx_yprep, nsig  = arc.detect_lines(arx_sky.flux.value)
    # Test
    assert len(arx_w[0]) == 1767
'''

