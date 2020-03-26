"""
Module to run tests on WaveImage class
Requires files in Development suite and an Environmental variable
"""
import os

import pytest
import glob
import numpy as np

from pypeit.tests.tstutils import load_kast_blue_masters, cooked_required
from pypeit.spectrographs.util import load_spectrograph


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

@cooked_required
def test_build_me():
    # Masters
    spectrograph = load_spectrograph('shane_kast_blue')
    edges, tilts_dict, wv_calib = load_kast_blue_masters(edges=True, tilts=True, wvcalib=True)
    slits = edges.get_slits()
    # Instantiate
    wvImg = wv_calib.BuildWaveImage(spectrograph, tilts_dict['tilts'], slits, wv_calib)
    # Build
    assert int(np.max(wvImg)) > 5510

