"""
Module to run tests on armasters
"""
import os
import numpy as np
import pytest

from astropy.io import fits

from pypeit import masterframe
from pypeit.images import buildimage

def data_root():
    return os.path.join(os.path.dirname(__file__), 'files')

def test_masterframe_methods():
    master_key = 'A_1_01'
    master_dir = data_root()

    # Filename
    filename = masterframe.construct_file_name(buildimage.ArcImage, master_key, master_dir=master_dir)
    assert isinstance(filename, str)

    # Header
    hdr = masterframe.build_master_header(buildimage.ArcImage, master_key, master_dir, 'shane_kast_blue')
    assert isinstance(hdr, fits.Header)

    # Get those keys!
    _master_key, _master_dir = masterframe.grab_key_mdir(hdr)
    assert _master_key == master_key
