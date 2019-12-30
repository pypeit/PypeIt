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
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)


@cooked_required
def test_build_me():
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    edges, tilts_dict, wv_calib = load_kast_blue_masters(edges=True, tilts=True, wvcalib=True)
    tslits_dict = edges.convert_to_tslits_dict()
    # Instantiate
    master_key = 'A_01_aa'
    master_dir = os.path.join(os.getenv('PYPEIT_DEV'), 'Cooked', 'Shane_Kast_blue')
    nslits = tslits_dict['nslits']
    maskslits = np.zeros(nslits, dtype=bool)
    det = 1
    wvImg = waveimage.WaveImage(tslits_dict, tilts_dict['tilts'], wv_calib,
                                spectrograph, det, maskslits,
                                master_key=master_key, master_dir=master_dir,
                                reuse_masters=True)
    # Build
    wave = wvImg.build_wave()
    assert int(np.max(wave)) > 5510

