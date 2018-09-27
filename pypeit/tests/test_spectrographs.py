# Module to test spectrograph read functions
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

import pytest
import glob

from pypeit.par.util import pypeit_root_directory
from pypeit import spectrographs

from pypeit.tests.tstutils import dev_suite_required

@dev_suite_required
def test_keckdeimos():
    s = spectrographs.keck_deimos.KeckDEIMOSSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_DEIMOS',
                                '830G_L_8400', 'd0914_0002.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck DEIMOS read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_kecklrisblue():
    s = spectrographs.keck_lris.KeckLRISBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_blue',
                                'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS blue read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_kecklrisred():
    s = spectrographs.keck_lris.KeckLRISRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck LRIS red read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_kecknires():
    s = spectrographs.keck_nires.KeckNIRESSpectrograph()
    # TODO: Any Keck NIRES files to read?


@dev_suite_required
def test_kecknirspec():
    s = spectrographs.keck_nirspec.KeckNIRSPECSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_NIRSPEC',
                                'NIRSPEC-1', 'NS.20160414.02171.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Keck NIRSPEC read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = os.path.join(pypeit_root_directory(), 'pypeit', 'tests', 'files',
                                'b1.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast blue read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_shanekastred():
    s = spectrographs.shane_kast.ShaneKastRedSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Shane_Kast_red',
                                '600_7500_d55', 'r112.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast red read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastredret():
    spectrographs.shane_kast.ShaneKastRedRetSpectrograph()
    # TODO: Any Shane Kast Red Ret files to read?


def test_tngdolores():
    s = spectrographs.tng_dolores.TNGDoloresSpectrograph()
    # TODO: Any TNG Dolores files to read?


@dev_suite_required
def test_vltxshooteruvb():
    s = spectrographs.vlt_xshooter.VLTXShooterUVBSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'UVB_1x1', 'XSHOO.2010-04-28T05_34_32.723.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter UVB read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_vltxshootervis():
    s = spectrographs.vlt_xshooter.VLTXShooterVISSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'VIS_1x1', 'XSHOO.2010-04-28T05_34_37.853.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter VIS read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_vltxshooternir():
    s = spectrographs.vlt_xshooter.VLTXShooterNIRSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'VLT_XSHOOTER',
                                'NIR', 'XSHOO.2016-08-02T08_45_49.494.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for VLT Xshooter NIR read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


@dev_suite_required
def test_whtisisblue():
    s = spectrographs.wht_isis.WhtIsisBlueSpectrograph()
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'WHT_ISIS_blue',
                                'long_R300B_d5300', 'r2324566.fit.gz')
    assert os.path.isfile(example_file), 'Could not find example file for WHT ISIS blue read.'
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


