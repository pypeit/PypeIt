"""
Module to run tests on armasters
"""
import os

from IPython import embed

import numpy as np

from astropy.io import fits

from pypeit import masterframe
from pypeit.images import buildimage
from pypeit.tests.tstutils import data_path


def test_masterframe_methods():
    master_key = 'A_1_DET01'
    master_dir = data_path('')

    # Filename
    Aimg = buildimage.ArcImage(None)
    Aimg.PYP_SPEC = 'shane_kast_blue'
    filename = masterframe.construct_file_name(Aimg, master_key, master_dir=master_dir)
    assert isinstance(filename, str)
    assert filename == data_path(f'Master{Aimg.master_type}_{master_key}.{Aimg.master_file_format}')

    # Header
    hdr = masterframe.build_master_header(Aimg, master_key, master_dir)
    assert isinstance(hdr, fits.Header)

    # Get those keys!
    _master_key, _master_dir = masterframe.grab_key_mdir(hdr)
    assert _master_key == master_key

    _master_key2, _master_dir2 = masterframe.grab_key_mdir(filename, from_filename=True)
    assert _master_key2 == master_key


