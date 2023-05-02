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


