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
from pypeit import metadata
from pypeit.core import arc

from pypeit.spectrographs.util import load_spectrograph

import pkg_resources

# JFH comementing out this test since the arcparm is now defunct.
#def test_setup_param():
#    """ Run the parameter setup script
#
#    Returns
#    -------
#    """
#    spectrograph = load_spectrograph('shane_kast_blue')
#    fitstbl = metadata.dummy_fitstbl()
#    # Run
#    arcparm = arc.setup_param(spectrograph, (2048,2048), fitstbl, 0)
#    for key in ['llist','disp','wvmnx']:
#        assert key in arcparm

def test_detect_lines():
    # Using Paranal night sky as an 'arc'
    sky_file = pkg_resources.resource_filename('pypeit', 'data/sky_spec/paranal_sky.fits')
    arx_sky = xspectrum1d.XSpectrum1D.from_file(sky_file)
    arx_amp, arx_cont, arx_cent, arx_wid, arx_centerr, arx_w, arx_yprep, _ \
            = arc.detect_lines(arx_sky.flux.value)
    assert len(arx_w[0]) == 4738

