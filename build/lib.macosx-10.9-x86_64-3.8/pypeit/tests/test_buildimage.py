"""
Module to run tests on ProcessImages class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images import buildimage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import procimg

kast_blue = load_spectrograph('shane_kast_blue')

@pytest.fixture
@dev_suite_required
def deimos_flat_files():
    # Longslit in dets 3,7
    deimos_flat_files = [os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos',
                                      '830G_L_8400', ifile) 
                            for ifile in ['d0914_0014.fits.gz', 'd0914_0015.fits.gz']]
    assert len(deimos_flat_files) == 2
    return deimos_flat_files

@pytest.fixture
@dev_suite_required
def kast_blue_bias_files():
    return glob.glob(os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue',
                                  '600_4310_d55', 'b1?.fits*'))

@dev_suite_required
def test_combine(deimos_flat_files):
    spectograph = load_spectrograph('keck_deimos')
    par = spectograph.default_pypeit_par()
    par['calibrations']['pixelflatframe']['process']['use_biasimage'] = False
    # DEIMOS
    deimos_flat = buildimage.buildimage_fromlist(spectograph, 3,
                                                  par['calibrations']['pixelflatframe'],
                                                  deimos_flat_files)
    # Process steps
    # Test
    assert deimos_flat.image.shape == (4096,2048)


