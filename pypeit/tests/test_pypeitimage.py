"""
Module to run tests on PypeItImage class
"""
import os

import pytest
import glob
import numpy as np

from astropy.io import fits

from pypeit.images import pypeitimage
from pypeit.tests.tstutils import dev_suite_required

def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

def test_dumb_instantiate():
    pypeitImage = pypeitimage.PypeItImage(None)
    assert pypeitImage.image is None


def test_image():
    # Just for a dummy image
    one_file = data_path('b1.fits.gz')
    data = fits.open(one_file)[0].data.astype(float)
    #
    pypeitImage = pypeitimage.PypeItImage(data)
    #
    assert pypeitImage.image.shape == (350, 2112)


def test_write():
    # Just for a dummy image
    data = np.ones((1000,1000))
    ivar = np.ones_like(data)
    mask = np.ones_like(data).astype(int)
    #
    pypeitImage = pypeitimage.PypeItImage(data, ivar=ivar, mask=mask)
    #
    outfile = data_path('tst.fits')
    pypeitImage.write(outfile)
    # Test
    hdul = fits.open(outfile)
    assert len(hdul) == 4
    assert hdul[2].name == 'IVAR'


def test_load():
    # This depends on the save method above
    tst_file = data_path('tst.fits')
    pypeitImage = pypeitimage.PypeItImage.from_file(tst_file)
    # Test
    assert pypeitImage.image is not None
    assert pypeitImage.ivar is not None
    assert pypeitImage.mask is not None
    assert pypeitImage.rn2img is None
    assert pypeitImage.head0 is not None
