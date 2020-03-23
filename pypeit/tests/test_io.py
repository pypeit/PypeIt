"""
Module to run tests on I/O core routines
"""
import os

import pytest

import numpy as np

from astropy.io import fits
from pypeit import io


def data_path(filename):
    data_dir = os.path.join(os.path.dirname(__file__), 'files')
    return os.path.join(data_dir, filename)

from IPython import embed


def test_init_hdus():
    # Dummy file

    ofile = data_path('tst_io.fits')
    if os.path.isfile(ofile):
        os.remove(ofile)

    # Write a file with several HDUs
    hdu0 = fits.PrimaryHDU()
    hdu0.header['EXT0001'] = 'DET01'
    hdu0.header['EXT0002'] = 'DET02'

    hdu1 = fits.ImageHDU()
    hdu1.name = 'DUMMY-DET01'
    hdu2 = fits.ImageHDU()
    hdu2.name = 'DUMMY-DET02'

    hdul = fits.HDUList([hdu0, hdu1, hdu2])
    hdul.writeto(ofile)

    # Now
    hdus, pri_hdr = io.init_hdus(1, ofile)
    assert 'EXT0001' not in hdus[0].header
    assert 'EXT0002' in hdus[0].header



