"""
Module to test spectrograph read functions
"""
import os
from copy import deepcopy

from IPython import embed

import pytest

from pypeit.pypmsgs import PypeItError
from pypeit import spectrographs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import pypeitsetup
from pypeit.tests import tstutils
from pypeit.tests.tstutils import data_path



def test_shanekastblue():
    s = spectrographs.shane_kast.ShaneKastBlueSpectrograph()
    example_file = data_path('b1.fits.gz')
    assert os.path.isfile(example_file), 'Could not find example file for Shane Kast blue read.'
    det=1
    _, data, hdu, exptime, rawdatasec_img, oscansec_img = s.get_rawimage(example_file, det)
    bpm = s.bpm(example_file, det)
    assert data.shape == (350, 2112)
    assert bpm.shape == (2048,350)


def test_select_detectors_pypeit_file():
    # Generate a PypeIt file
    pypeItFile = tstutils.make_shane_kast_blue_pypeitfile()
    pypeit_file = data_path('test.pypeit')
    pypeItFile.write(pypeit_file)

    # Perform the setup
    setup = pypeitsetup.PypeItSetup.from_pypeit_file(pypeit_file)
    par, spectrograph, fitstbl = setup.run()

    assert spectrograph.select_detectors(subset=par['rdx']['detnum']) == [1], \
            'Incorrect detectors selected.'

    # Clean-up
    os.remove(data_path('test.pypeit'))


def test_select_detectors_mosaic():

    spec = load_spectrograph('gemini_gmos_north_ham')

    # Invalid detector
    with pytest.raises(PypeItError):
        spec.select_detectors(subset=4)
    # Invalid mosaic
    with pytest.raises(PypeItError):
        spec.select_detectors(subset=(2,3))

    # Valid
    assert spec.select_detectors() == [1,2,3], 'Bad detector selection'
    # Valid
    assert spec.select_detectors(subset=[3,(1,2,3)]) == [3,(1,2,3)], 'Bad detector selection'


def test_list_detectors_deimos():
    deimos = load_spectrograph('keck_deimos')
    dets = deimos.list_detectors()
    assert dets.ndim == 2, 'DEIMOS has a 2D array of detectors'
    assert dets.size == 8, 'DEIMOS has 8 detectors'
    mosaics = deimos.list_detectors(mosaic=True)
    assert mosaics.ndim == 1, 'Mosaics are listed as 1D arrays'
    assert mosaics.size == 4, 'DEIMOS has 4 predefined mosaics'


def test_list_detectors_mosfire():
    mosfire = load_spectrograph('keck_mosfire')
    dets = mosfire.list_detectors()
    assert dets.ndim == 1, 'MOSFIRE has a 1D array of detectors'
    assert dets.size == 1, 'MOSFIRE has 1 detector'
    with pytest.raises(PypeItError):
        mosaics = mosfire.list_detectors(mosaic=True)


def test_list_detectors_mods():
    mods = load_spectrograph('lbt_mods1r')
    dets = mods.list_detectors()
    assert dets.ndim == 1, 'MODS1R has a 1D array of detectors'
    assert dets.size == 1, 'MODS1R has 1 detector'
    with pytest.raises(PypeItError):
        mosaics = mods.list_detectors(mosaic=True)


def test_list_detectors_hires():
    hires = load_spectrograph('keck_hires')
    dets = hires.list_detectors()
    assert dets.ndim == 1, 'HIRES has a 1D array of detectors'
    assert dets.size == 3, 'HIRES has 3 detectors'
    mosaics = hires.list_detectors(mosaic=True)
    assert mosaics.ndim == 1, 'Mosaics are listed as 1D arrays'
    assert mosaics.size == 1, 'HIRES has 1 predefined mosaic'


def test_configs():

    spec = load_spectrograph('keck_deimos')
    cfg1 = dict(amp='"SINGLE:B"',
                binning='1,1',
                decker='LongMirr',
                dispangle=8099.98291016,
                dispname='830G',
                filter1='OG550')
    cfg2 = dict(amp='"SINGLE:B"',
                binning='1,1',
                decker='LongMirr',
                dispangle=8399.93554688,
                dispname='830G',
                filter1='OG550')

    assert spec.same_configuration([cfg1,cfg1]), 'Configurations should be the same'
    assert not spec.same_configuration([cfg1,cfg2]), 'Configurations should be different'

    cfg3 = deepcopy(cfg1)
    cfg3['dispangle'] *= (1.+spec.meta['dispangle']['rtol']/2)

    assert spec.same_configuration([cfg1,cfg3]), \
        'Configurations should be the same within tolerance'

    cfg3 = deepcopy(cfg1)
    cfg3['dispangle'] *= (1.+2*spec.meta['dispangle']['rtol'])

    assert not spec.same_configuration([cfg1,cfg3]), \
        'Configurations should not be the same within tolerance'


