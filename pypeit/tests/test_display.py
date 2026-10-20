"""
Module to run tests on arcoadd
"""
from IPython import embed

import numpy as np
from astropy.table import Table

from pypeit.display import plugins_available
from pypeit.slittrace import SlitTraceSet
from pypeit.scripts.show_2dspec import wave_rms_image

def test_plugin():
    success, report = plugins_available(return_report=True)
    assert success, report



def test_wave_rms_image():
    # Three 10-pixel wide slits separated by 10-pixel gaps
    nspec, nspat = 50, 70
    left = np.tile(np.array([10., 30., 50.]), (nspec, 1))
    slits = SlitTraceSet(left_init=left, right_init=left + 10., pypeline='MultiSlit',
                         nspat=nspat, PYP_SPEC='dummy')

    # No diagnostics table means no RMS image
    assert wave_rms_image(slits, None) is None, 'Should not build an image without a table'

    # Table rows are matched by ID, not order; the last slit has no solution
    wavesol = Table()
    wavesol['SpatOrderID'] = slits.spat_id[[1, 0]]
    wavesol['RMS'] = [0.2, 0.1]
    rms_img = wave_rms_image(slits, wavesol)

    assert rms_img.shape == (nspec, nspat), 'Bad image shape'
    spec = nspec // 2
    assert rms_img[spec, 15] == 0.1, 'Bad RMS for the first slit'
    assert rms_img[spec, 35] == 0.2, 'Bad RMS for the second slit'
    assert rms_img[spec, 55] == 0., 'Slit without a solution should be 0'
    assert rms_img[spec, 25] == 0., 'Pixels between slits should be 0'

    # Masked slits are not included
    slits.mask[0] = slits.bitmask.turn_on(slits.mask[0], 'BADWVCALIB')
    assert wave_rms_image(slits, wavesol)[spec, 15] == 0., 'Masked slit should be 0'

    # The image follows the slits when they are shifted for flexure
    assert wave_rms_image(slits, wavesol, flexure=5.)[spec, 43] == 0.2, 'Bad flexure shift'
