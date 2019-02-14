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
from pypeit.par import pypeitpar

par = pypeitpar.ProcessImagesPar()

@dev_suite_required
def test_load_deimos():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L_8400',
                         'd0914_0014.fits.gz')
    proc = ProcessImages('keck_deimos', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('DEIMOS test data section failed.')

@dev_suite_required
def test_load_lris():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_LRIS_blue',
                         'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    proc = ProcessImages('keck_lris_blue', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('LRIS test data section failed.')

@dev_suite_required
def test_load_nires():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRES', 'NIRES',
                         's180604_0004.fits.gz')
    proc = ProcessImages('keck_nires', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('NIRES test data section failed.')

@dev_suite_required
def test_load_nirspec():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRSPEC', 'NIRSPEC-1',
                         'NS.20160414.02604.fits.gz')
    proc = ProcessImages('keck_nirspec_low', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('NIRSPEC test data section failed.')

@dev_suite_required
def test_load_kast():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Shane_Kast_blue', '600_4310_d55',
                         'b1.fits.gz')
    proc = ProcessImages('shane_kast_blue', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Shane Kast test data section failed.')

@dev_suite_required
def test_load_isis():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'WHT_ISIS_blue', 'long_R300B_d5300',
                         'r2324566.fit.gz')
    proc = ProcessImages('wht_isis_blue', par, files)
    proc.load_images()
    try:
        # First amplifier
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('WHT ISIS test data section failed.')

@dev_suite_required
def test_load_vlt_xshooter_uvb():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER',
                         'UVB_1x1/XSHOO.2010-04-28T05:34:32.723.fits.gz')
    proc = ProcessImages('vlt_xshooter_uvb', par, files)
    proc.load_images()
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

    proc = ProcessImages('vlt_xshooter_vis', par, None)
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
    proc = ProcessImages('vlt_xshooter_nir', par, files)
    proc.load_images()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('VLT XSHOOTER NIR test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_gnirs():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Gemini_GNIRS/GNIRS/',
                         'cN20170331S0206.fits')
    proc = ProcessImages('gemini_gnirs', par, files)
    proc.load_images()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Gemini GNIRS test data section failed: {0}'.format(files))

'''
@dev_suite_required
def test_load_fire():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Magellan_FIRE/FIRE',
                         'fire_0029.fits.gz')
    proc = ProcessImages('magellan_fire', par, files)
    proc.load_images()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Magellan FIRE test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_mage():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Magellan_MAGE/MAGE',
                         'mage1002.fits.gz')
    proc = ProcessImages('magellan_mage', par, files)
    proc.load_images()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Magellan MAGE test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_hires():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_HIRES/RED',
                         'hires0009.fits.gz')
    proc = ProcessImages('keck_hires_red', par, files)
    proc.load_images()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Keck HIRES test data section failed: {0}'.format(files))
'''

