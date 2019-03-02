"""
Module to run tests on WaveImage class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import load_kast_blue_masters, cooked_required
from pypeit import waveimage


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@cooked_required
def test_build_me():
    # Masters
    spectrograph, tslits_dict, tilts_dict, wv_calib \
            = load_kast_blue_masters(get_spectrograph=True, tslits=True, tilts=True, wvcalib=True)
    # Instantiate
    master_key = 'A_01_aa'
    # TODO: This should always point to Cooked
    root_path = data_path('MF') if os.getenv('PYPEIT_DEV') is None \
                    else os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'MF')
    master_dir = root_path+'_'+spectrograph.spectrograph
    nslits = tslits_dict['slit_left'].shape[1]
    maskslits = np.zeros(nslits, dtype=bool)
    wvImg = waveimage.WaveImage(tslits_dict, tilts_dict['tilts'], wv_calib, spectrograph, maskslits,
                                master_key=master_key, master_dir=master_dir, reuse_masters=True)
    # Build
    wave = wvImg.build_wave()
    assert int(np.max(wave)) > 5510

