# Module to run loading images for different images served using the
# ProcessImages class
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import pytest
import glob
import numpy as np

from pypeit.processimages import ProcessImages
from pypeit.tests.tstutils import dev_suite_required

#@dev_suite_required
#def test_load_deimos():
#    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L_8400',
#                         'd0914_0014.fits.gz')
#    proc = processimages.ProcessImages('keck_deimos', None)
#    proc.load_images(files)
#    try:
#        # First amplifier
#        data_img = proc.raw_images[0][proc.datasec[0][0]]
#    except:
#        pytest.fail('DEIMOS test data section failed.')

@dev_suite_required
def test_load_vlt_xshooter_uvb():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER',
                         'UVB_1x1/XSHOO.2010-04-28T05:34:32.723.fits.gz')
    proc = ProcessImages('vlt_xshooter_uvb', None)
    proc.load_images(files)
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('VLT XSHOOTER UVB test data section failed: {0}'.format(files))


@dev_suite_required
def test_load_vlt_xshooter_vis():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER')
    files = [ os.path.join(root, 'VIS_1x1/XSHOO.2010-04-28T05:34:37.853.fits.gz'),
              os.path.join(root, 'VIS_1x2/XSHOO.2016-08-02T08:45:46.510.fits.gz'),
              os.path.join(root, 'VIS_2x2/XSHOO.2016-10-08T00:51:04.703.fits.gz') ]

    proc = ProcessImages('vlt_xshooter_vis', None)
    for f in files:
        proc.load_images(f)
        try:
            data_img = proc.raw_images[0][proc.datasec[0][0]]
        except:
            pytest.fail('VLT XSHOOTER VIS test data section failed: {0}'.format(f))


@dev_suite_required
def test_load_vlt_xshooter_nir():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER',
                         'NIR/XSHOO.2016-08-02T08:45:49.494.fits.gz')
    proc = ProcessImages('vlt_xshooter_nir', None)
    proc.load_images(files)
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('VLT XSHOOTER NIR test data section failed: {0}'.format(files))

