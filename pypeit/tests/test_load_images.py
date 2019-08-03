"""
Module to run loading images for different images served using the
ProcessImages class
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images.calibrationimage import CalibrationImage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core import procimg

par = pypeitpar.ProcessImagesPar()

# Dumb wrapper because I am too lazy to replace the old approach
def ProcessImages(specstr, par, files, det=1):
    spec = load_spectrograph(specstr)
    calibImage = CalibrationImage(spec, det, par, files=files)
    return calibImage

def grab_img(proc, filename):
    img, hdu, exptime, rawdatasec_img, oscansec_img = proc.spectrograph.get_rawimage(filename, proc.det)
    data_img, _ = procimg.rect_slice_with_mask(img, rawdatasec_img)
    return data_img


@dev_suite_required
def test_load_deimos():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_DEIMOS', '830G_L_8400',
                         'd0914_0014.fits.gz')
    proc = ProcessImages('keck_deimos', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc, files)
    except:
        pytest.fail('DEIMOS test data section failed.')

@dev_suite_required
def test_load_lris():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_LRIS_blue',
                         'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    proc = ProcessImages('keck_lris_blue', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc, files)
    except:
        pytest.fail('LRIS test data section failed.')

@dev_suite_required
def test_load_nires():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRES', 'NIRES',
                         's180604_0004.fits.gz')
    proc = ProcessImages('keck_nires', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc, files)
    except:
        pytest.fail('NIRES test data section failed.')

@dev_suite_required
def test_load_nirspec():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Keck_NIRSPEC', 'NIRSPEC-1',
                         'NS.20160414.02604.fits.gz')
    proc = ProcessImages('keck_nirspec_low', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc, files)
    except:
        pytest.fail('NIRSPEC test data section failed.')

@dev_suite_required
def test_load_kast():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'Shane_Kast_blue', '600_4310_d55',
                         'b1.fits.gz')
    proc = ProcessImages('shane_kast_blue', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc, files)
    except:
        pytest.fail('Shane Kast test data section failed.')


@dev_suite_required
def test_load_vlt_xshooter_uvb():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER',
                         'UVB_1x1/XSHOO.2010-04-28T05:34:32.723.fits.gz')
    proc = ProcessImages('vlt_xshooter_uvb', par, files)
    proc.build_image()
    try:
        data_img = grab_img(proc, files)
    except:
        pytest.fail('VLT XSHOOTER UVB test data section failed: {0}'.format(files))


@dev_suite_required
def test_load_vlt_xshooter_vis():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER')
    files = [ os.path.join(root, 'VIS_1x1/XSHOO.2010-04-28T05:34:37.853.fits.gz'),
              os.path.join(root, 'VIS_2x1/XSHOO.2016-08-02T08:45:46.510.fits.gz'),
              os.path.join(root, 'VIS_2x2/XSHOO.2016-10-08T00:51:04.703.fits.gz') ]

    for f in files:
        proc = ProcessImages('vlt_xshooter_vis', par, f)
        proc.build_image()
        try:
            data_img = grab_img(proc, f)
        except:
            pytest.fail('VLT XSHOOTER VIS test data section failed: {0}'.format(f))

@dev_suite_required
def test_load_vlt_xshooter_nir():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/VLT_XSHOOTER',
                         'NIR/XSHOO.2016-08-02T08:45:49.494.fits.gz')
    proc = ProcessImages('vlt_xshooter_nir', par, files)
    proc.build_image()
    try:
        data_img = grab_img(proc, files)
    except:
        pytest.fail('VLT XSHOOTER NIR test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_gnirs():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Gemini_GNIRS/GNIRS/',
                         'cN20170331S0206.fits')
    proc = ProcessImages('gemini_gnirs', par, files)
    proc.build_image()
    try:
        data_img = grab_img(proc, files)
    except:
        pytest.fail('Gemini GNIRS test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_mage():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Magellan_MAGE/1x1',
                         'mage0050.fits')
    proc = ProcessImages('magellan_mage', par, files)
    proc.build_image()
    try:
        data_img = grab_img(proc, files)
    except:
        pytest.fail('Magellan MAGE test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_gmos():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Gemini_GMOS/GS_HAM_R400_700',
                         'S20181005S0086.fits.gz')
    proc = ProcessImages('gemini_gmos_south_ham', par, files)
    proc.build_image()
    try:
        data_img = grab_img(proc, files)
    except:
        pytest.fail('Gemini GMOS test data section failed: {0}'.format(files))

'''
@dev_suite_required
def test_load_fire():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Magellan_FIRE/FIRE',
                         'fire_0029.fits.gz')
    proc = ProcessImages('magellan_fire', par, files)
    proc.build_image()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Magellan FIRE test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_hires():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/Keck_HIRES/RED',
                         'hires0009.fits.gz')
    proc = ProcessImages('keck_hires_red', par, files)
    proc.build_image()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Keck HIRES test data section failed: {0}'.format(files))
        
@dev_suite_required
def test_load_isis():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'WHT_ISIS_blue', 'long_R300B_d5300',
                         'r2324566.fit.gz')
    proc = ProcessImages('wht_isis_blue', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc)
    except:
        pytest.fail('WHT ISIS test data section failed.')
'''

