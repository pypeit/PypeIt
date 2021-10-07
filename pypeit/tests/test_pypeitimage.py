"""
Module to run tests on PypeItImage class
"""
import os

from IPython import embed

import numpy as np

from pypeit.images import pypeitimage
from pypeit.images import imagebitmask

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_full():
    pypeitImage = pypeitimage.PypeItImage(np.ones((1000, 1000)))
    pypeitImage.reinit_mask()
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


def test_bitmask():
    bm = imagebitmask.ImageBitMask()

    # Check that a few of the bit numbers are correct (i.e., that the order
    # didn't get messed up)
    assert bm.bits['BPM'] == 0, 'BPM bit number changed'
    assert bm.bits['CR'] == 1, 'CR bit number changed'
    assert bm.bits['OFFSLITS'] == 4, 'OFFSLITS bit number changed'
    assert bm.bits['EXTRACT'] == 8, 'EXTRACT bit number changed'



