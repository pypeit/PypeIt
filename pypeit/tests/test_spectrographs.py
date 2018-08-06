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

# These tests are not run on Travis
if os.getenv('PYPEIT_DEV') is None:
    skip_test=True
else:
    skip_test=False

def test_keckdeimos():
    s = spectrographs.keck_deimos.KeckDEIMOSSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_DEIMOS',
                                '830G_L', 'd0914_0002.fits')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Keck DEIMOS read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_kecklrisblue():
    s = spectrographs.keck_lris.KeckLRISBSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_blue',
                                'long_400_3400_d560', 'LB.20160109.14149.fits.gz')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Keck LRIS Blue read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_kecklrisred():
    s = spectrographs.keck_lris.KeckLRISRSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_LRIS_red',
                                'long_600_7500_d560', 'LR.20160216.05529.fits')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Keck LRIS Red read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(filename=example_file)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_kecknirspec():
    s = spectrographs.keck_nirspec.KeckNIRSPECSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Keck_NIRSPEC',
                                'NIRSPEC-1', 'NS.20160414.02171.fits.gz')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Keck LRIS Red read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = os.path.join(pypeit_root_directory(), 'pypeit', 'tests', 'files',
                                'b1.fits.gz')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Shane Kast Blue read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastred():
    s = spectrographs.shane_kast.ShaneKastRedSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'Shane_Kast_red',
                                '600_7500_d55', 'r112.fits.gz')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for Shane Kast Red read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'


def test_shanekastredret():
    spectrographs.shane_kast.ShaneKastRedRetSpectrograph()
    # TODO: Any Shane Kast Red Ret files to read?


def test_tngdolores():
    s = spectrographs.tng_dolores.TngDoloresSpectrograph()
    # TODO: Any TNG Dolores files to read?


def test_whtisisblue():
    s = spectrographs.wht_isis.WhtIsisBlueSpectrograph()
    if skip_test:
        return
    example_file = os.path.join(os.environ['PYPEIT_DEV'], 'RAW_DATA', 'WHT_ISIS_blue',
                                'long_R300B_d5300', 'r2324566.fit.gz')
    if not os.path.isfile(example_file):
        raise FileNotFoundError('Could not find example file for WHT ISIS Blue read.')
    data, _ = s.load_raw_frame(example_file)
    bpm = s.bpm(shape=data.shape)
    assert data.shape == bpm.shape, 'Image and BPM have different shapes!'

