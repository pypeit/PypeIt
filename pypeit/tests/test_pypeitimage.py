"""
Module to run tests on PypeItImage class
"""
import os

import pytest
import numpy as np

from pypeit.images import pypeitimage

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_full():
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    # Mask
    pypeitImage.fullmask = np.zeros((1000, 1000), dtype=np.int64)
    # Full datamodel
    full_datamodel = pypeitImage.full_datamodel()
    assert 'gain' in full_datamodel.keys()
    assert 'detector' in pypeitImage.keys()

    # I/O
    outfile = data_path('tst_pypeitimage.fits')
    pypeitImage.to_file(outfile, overwrite=True)
    _pypeitImage = pypeitimage.PypeItImage.from_file(outfile)

    # Cleanup
    os.remove(outfile)

    # Test
    assert isinstance(_pypeitImage.image, np.ndarray)
    assert _pypeitImage.ivar is None

