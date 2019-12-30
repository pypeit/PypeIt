"""
Module to test spectrograph read functions
"""
import os

import pytest
import glob

from pkg_resources import resource_filename

from pypeit import spectrographs
from pypeit.core import procimg

from pypeit.tests.tstutils import dev_suite_required

# TODO: Add a test for Gemini GNIRS

@dev_suite_required
def test_keckdeimos():
    s = spectrographs.keck_deimos.KeckDEIMOSSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_DEIMOS',
                                '830G_L_8400', 'd0914_0002.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck DEIMOS read.'
    det = 2
    #data, _ = s.load_raw_frame(example_file, det=det)
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det) #shape=shape) # filename=example_file)
    assert data.shape == (4096,2128)
    assert bpm.shape == (4096,2048)


@dev_suite_required
def test_kecklrisblue():
    s = spectrographs.keck_lris.KeckLRISBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_blue',
                                'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS blue read.'
    det = 2
    #data, _ = s.load_raw_frame(example_file, det=det)
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det) #shape=shape) # filename=example_file)
    assert data.shape == (2048,1154)
    assert bpm.shape == (2048,1024)


@dev_suite_required
def test_kecklrisred():
    s = spectrographs.keck_lris.KeckLRISRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS red read.'
    det = 1
    #data, _ = s.load_raw_frame(example_file, det=det)
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    #
    bpm = s.bpm(example_file, det)#, debug=True) #shape=shape) # filename=example_file)
    assert data.shape == (2068,1110)
    assert bpm.shape == (2048,1024)


@dev_suite_required
def test_kecknires():
    s = spectrographs.keck_nires.KeckNIRESSpectrograph()
    # TODO: Any Keck NIRES files to read?


@dev_suite_required
def test_kecknirspec():
    s = spectrographs.keck_nirspec.KeckNIRSPECSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_NIRSPEC',
                                'LOW_NIRSPEC-1', 'NS.20160414.02637.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck NIRSPEC read.'
    #data, _ = s.load_raw_frame(example_file)
    det=1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = os.path.join(resource_filename('pypeit', 'tests'), 'files',
                                'b1.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast blue read.'
    det=1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (350, 2112)
    assert bpm.shape == (2048,350)


@dev_suite_required
def test_shanekastredret():
    s = spectrographs.shane_kast.ShaneKastRedRetSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Shane_Kast_red',
                                '600_7500_d55_ret', 'r112.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast red read.'
    det = 1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (250, 1232)
    assert bpm.shape == (1200, 250)


def test_shanekastred():
    spectrographs.shane_kast.ShaneKastRedSpectrograph()
    # TODO: Any Shane Kast Red files to read?


def test_tngdolores():
    s = spectrographs.tng_dolores.TNGDoloresSpectrograph()
    # TODO: Any TNG Dolores files to read?


@dev_suite_required
def test_vltxshooteruvb():
    s = spectrographs.vlt_xshooter.VLTXShooterUVBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'UVB_1x1', 'XSHOO.2010-04-28T05:34:32.723.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter UVB read.'
    det = 1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (3000, 2144)
    assert bpm.shape == (3000, 2048)


@dev_suite_required
def test_vltxshootervis():
    s = spectrographs.vlt_xshooter.VLTXShooterVISSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'VIS_1x1', 'XSHOO.2010-04-28T05:34:37.853.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter VIS read.'
    det = 1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (4000, 2106)
    assert bpm.shape == (4000, 2048)


@dev_suite_required
def test_vltxshooternir():
    s = spectrographs.vlt_xshooter.VLTXShooterNIRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'NIR', 'XSHOO.2016-08-02T08:45:49.494.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter NIR read.'
    det = 1
    data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (1100,2048)
    assert bpm.shape == (2045, 1097)



