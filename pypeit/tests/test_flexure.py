"""
Module to run tests on simple fitting routines for arrays
"""

import numpy as np

from linetools.spectra.io import readspec

from pypeit.core import flexure
from pypeit import dataPaths
from pypeit.core.wavecal import autoid

from IPython import embed


def test_flex_shift():
    # Dummy slf
    # Read spectra
    spec_file = dataPaths.tests.get_file_path('obj_lrisb_600_sky.fits', to_pkg='symlink')
    obj_spec = readspec(str(spec_file))
    arx_file = dataPaths.sky_spec.get_file_path('sky_LRISb_600.fits', to_pkg='symlink')

    # Call
    flex_dict = flexure.spec_flex_shift(obj_spec, sky_file=arx_file, mxshft=60)

    assert np.abs(flex_dict['shift'] - 43.5) < 0.1


def test_flex_image():
    """ Test the image alignment """
    # Generate some fake data
    sz = 100
    xx, yy = np.meshgrid(np.arange(sz),np.arange(sz))
    img = np.exp(-0.5*((sz/2-xx)**2 - (sz/2-yy)**2)/16)
    # Check odd/even image sizes to make sure the offset is always zero for identical input images
    xshft, yshft = flexure.calculate_image_offset(img, img)
    assert(np.abs(xshft) < 1.0e-6 and np.abs(yshft) < 1.0e-6)
    xshft, yshft = flexure.calculate_image_offset(img[:, :-1], img[:, :-1])
    assert(np.abs(xshft) < 1.0e-6 and np.abs(yshft) < 1.0e-6)
    xshft, yshft = flexure.calculate_image_offset(img[:-1, :], img[:-1, :])
    assert(np.abs(xshft) < 1.0e-6 and np.abs(yshft) < 1.0e-6)
    xshft, yshft = flexure.calculate_image_offset(img[:-1, :-1], img[:-1, :-1])
    assert(np.abs(xshft) < 1.0e-6 and np.abs(yshft) < 1.0e-6)

