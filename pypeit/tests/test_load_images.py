"""
Module to run loading images for different images served using the
RawImage class
"""
import os

import pytest
import glob
import numpy as np

from pypeit.images.rawimage import RawImage
from pypeit.tests.tstutils import dev_suite_required
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph

par = pypeitpar.ProcessImagesPar()

'''
# Dumb wrapper because I am too lazy to replace the old approach
def load_RawImage(specstr, par, files, det=1):
    rawImage = RawImage(files[0], spec, det)
    #CalibrationImage(spec, det, par, files=files)
    #return calibImage
    return rawImage

def grab_img(proc, filename):
    spec = load_spectrograph(specstr)
    img, hdu, exptime, rawdatasec_img, oscansec_img = proc.spectrograph.get_rawimage(filename, proc.det)
    data_img, _ = procimg.rect_slice_with_mask(img, rawdatasec_img)
    return data_img
'''

def grab_img(specstr, rawfile, det=1):
    spec = load_spectrograph(specstr)
    rawImage = RawImage(rawfile, spec, det)
    return rawImage


@dev_suite_required
def test_load_deimos():
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_deimos', '830G_L_8400',
                         'd0914_0014.fits.gz')
    try:
        # First amplifier
        data_img = grab_img('keck_deimos', ifile)
    except:
        pytest.fail('DEIMOS test data section failed.')

@dev_suite_required
def test_load_lris():
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_lris_blue',
                         'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    try:
        # First amplifier
        data_img = grab_img('keck_lris_blue', ifile)
    except:
        pytest.fail('LRIS test data section failed.')

@dev_suite_required
def test_load_nires():
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_nires', 'NIRES',
                         's180604_0004.fits.gz')
    try:
        # First amplifier
        data_img = grab_img('keck_nires', ifile)
    except:
        pytest.fail('NIRES test data section failed.')

@dev_suite_required
def test_load_nirspec():
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'keck_nirspec', 'LOW_NIRSPEC-1',
                         'NS.20160414.02604.fits.gz')
    try:
        # First amplifier
        data_img = grab_img('keck_nirspec_low', ifile)
    except:
        pytest.fail('NIRSPEC test data section failed.')

@dev_suite_required
def test_load_kast():
    ifile = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'shane_kast_blue', '600_4310_d55',
                         'b1.fits.gz')
    try:
        # First amplifier
        data_img = grab_img('shane_kast_blue', ifile)
    except:
        pytest.fail('Shane Kast test data section failed.')


@dev_suite_required
def test_load_vlt_xshooter_uvb():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter',
                         'UVB_1x1/XSHOO.2010-04-28T05:34:32.723.fits.gz')
    try:
        data_img = grab_img('vlt_xshooter_uvb', ifile)
    except:
        pytest.fail('VLT XSHOOTER UVB test data section failed: {0}'.format(ifile))


@dev_suite_required
def test_load_vlt_xshooter_vis():

    root = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter')
    files = [ os.path.join(root, 'VIS_1x1/XSHOO.2010-04-28T05:34:37.853.fits.gz'),
              os.path.join(root, 'VIS_2x1/XSHOO.2016-08-02T08:45:46.510.fits.gz'),
              os.path.join(root, 'VIS_2x2/XSHOO.2016-10-08T00:51:04.703.fits.gz') ]

    for f in files:
        try:
            data_img = grab_img('vlt_xshooter_vis', f)
        except:
            pytest.fail('VLT XSHOOTER VIS test data section failed: {0}'.format(f))

@dev_suite_required
def test_load_vlt_xshooter_nir():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/vlt_xshooter',
                         'NIR/XSHOO.2016-08-02T08:45:49.494.fits.gz')
    try:
        data_img = grab_img('vlt_xshooter_nir', ifile)
    except:
        pytest.fail('VLT XSHOOTER NIR test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_gnirs():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/gemini_gnirs/32_SB_SXD/',
                         'cN20170331S0206.fits')
    try:
        data_img = grab_img('gemini_gnirs', ifile)
    except:
        pytest.fail('Gemini GNIRS test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_mage():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/magellan_mage/1x1',
                         'mage0050.fits')
    try:
        data_img = grab_img('magellan_mage', ifile)
    except:
        pytest.fail('Magellan MAGE test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_gmos():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/gemini_gmos/GS_HAM_R400_700',
                         'S20181005S0086.fits.gz')
    try:
        data_img = grab_img('gemini_gmos_south_ham', ifile)
    except:
        pytest.fail('Gemini GMOS test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_osiris():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/gtc_osiris/R2500R',
                         '0002851159-20210217-OSIRIS-OsirisBias.fits')
    try:
        data_img = grab_img('gtc_osiris', ifile)
    except:
        pytest.fail('GTC OSIRIS test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_bok():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/bok_bc/600',
                         'g0005.fits')
    try:
        data_img = grab_img('bok_bc', ifile)
    except:
        pytest.fail('Bok BC test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_efosc2():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/ntt_efosc2/gr6',
                         'EFOSC.2020-02-12T02:03:38.359.fits')
    try:
        data_img = grab_img('ntt_efosc2', ifile)
    except:
        pytest.fail('NTT/EFOSC2 test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_goodman():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/soar_goodman_red/M2',
                         '0320_FRB210320_host_05-04-2021.fits.fz')
    try:
        data_img = grab_img('soar_goodman_red', ifile)
    except:
        pytest.fail('Bok BC test data section failed: {0}'.format(ifile))

@dev_suite_required
def test_load_deveny():
    ifile = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/ldt_deveny/DV2',
                         '20210522.0001.fits')
    try:
        data_img = grab_img('ldt_deveny', ifile)
    except:
        pytest.fail('LDT DeVeny test data section failed: {0}'.format(ifile))

'''
@dev_suite_required
def test_load_fire():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/magellan_fire/FIRE',
                         'fire_0029.fits.gz')
    proc = ProcessImages('magellan_fire', par, files)
    proc.build_image()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Magellan FIRE test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_hires():
    files = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA/keck_hires/RED',
                         'hires0009.fits.gz')
    proc = ProcessImages('keck_hires_red', par, files)
    proc.build_image()
    try:
        data_img = proc.raw_images[0][proc.datasec[0][0]]
    except:
        pytest.fail('Keck HIRES test data section failed: {0}'.format(files))

@dev_suite_required
def test_load_isis():
    files = os.path.join(os.getenv('PYPEIT_DEV'), 'RAW_DATA', 'wht_isis_blue', 'long_R300B_d5300',
                         'r2324566.fit.gz')
    proc = ProcessImages('wht_isis_blue', par, files)
    proc.build_image()
    try:
        # First amplifier
        data_img = grab_img(proc)
    except:
        pytest.fail('WHT ISIS test data section failed.')
'''
