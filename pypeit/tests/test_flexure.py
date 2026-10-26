"""
Module to run tests on simple fitting routines for arrays
"""

import numpy as np

import pytest

from pypeit.core import flexure
from pypeit import dataPaths
from pypeit.core.wavecal import autoid
from pypeit import onespec

from IPython import embed


@pytest.mark.remote_data
def test_flex_shift():
    # Dummy slf
    # Read spectra
    spec_file = dataPaths.tests.get_file_path('obj_lrisb_600_sky.fits', to_pkg='symlink')
    #obj_spec = readspec(str(spec_file))
    obj_spec = onespec.OneSpec.from_xspec_file(str(spec_file))
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



def test_archive_fwhm_cache(monkeypatch):
    """
    The FWHM of an archived sky spectrum should be measured only once per sky
    file, no matter how many times get_archive_spectrum is called.
    """
    ncalls = []
    _measure_fwhm = autoid.measure_fwhm

    def counting_measure_fwhm(*args, **kwargs):
        ncalls.append(1)
        return _measure_fwhm(*args, **kwargs)

    monkeypatch.setattr(flexure.autoid, 'measure_fwhm', counting_measure_fwhm)
    flexure._archive_sky_fwhm.cache_clear()
    try:
        spec1, fwhm1 = flexure.get_archive_spectrum('paranal_sky.fits')
        assert len(ncalls) == 1, 'FWHM should be measured on the first call'
        spec2, fwhm2 = flexure.get_archive_spectrum('paranal_sky.fits')
        assert len(ncalls) == 1, 'FWHM should not be re-measured for the same sky file'
        assert fwhm1 == fwhm2
        # A new spectrum object must be returned each time
        assert spec1 is not spec2
        assert np.array_equal(spec1.flux, spec2.flux)
        # The cached value matches a direct measurement
        assert fwhm1 == _measure_fwhm(spec1.flux, sigdetect=4., fwhm=4.)

        # A different sky file is measured separately
        _, fwhm3 = flexure.get_archive_spectrum('sky_kastb_600.fits')
        assert len(ncalls) == 2, 'A different sky file should be measured'

        # Clearing the cache forces a re-measurement
        flexure._archive_sky_fwhm.cache_clear()
        _, fwhm4 = flexure.get_archive_spectrum('paranal_sky.fits')
        assert len(ncalls) == 3, 'FWHM should be re-measured after clearing the cache'
        assert fwhm4 == fwhm1
    finally:
        flexure._archive_sky_fwhm.cache_clear()
