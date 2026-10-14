"""
Module to run tests on simple fitting routines for arrays
"""

import numpy as np

import pytest
from scipy.special import erf

from pypeit import PypeItError
from pypeit.core import flexure
from pypeit.slittrace import SlitTraceSet
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
    flex_dict = flexure.spec_flexure_shift(obj_spec, sky_file=arx_file, mxshft=60)

    assert np.abs(flex_dict['spec_flexure'] - 43.5) < 0.1


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



def _synthetic_slit_image(dl, dr, nslit=2, nspec=200, sig=1.2, drift=3.0):
    """
    Build a synthetic (noisy) slit image with nslit slits, whose left and right
    edges are offset from the nominal slit edges by dl and dr pixels.

    The slit edges drift by a few pixels along the spectral direction, as is
    typical of real data, which provides the sub-pixel sampling required to
    refine the edge positions.
    """
    shift = np.linspace(0, drift, nspec)
    left = 30. + 60. * np.arange(nslit)[None, :] + shift[:, None]
    right = left + 30.
    nspat = int(left[0, -1]) + 70
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit', nspat=nspat,
                         PYP_SPEC='dummy')
    xx = np.arange(nspat, dtype=float)[None, :]
    img = np.zeros((nspec, nspat))
    for slit in range(nslit):
        img += 0.5 * (erf((xx - left[:, slit, None] - dl) / (np.sqrt(2) * sig))
                      - erf((xx - right[:, slit, None] - dr) / (np.sqrt(2) * sig)))
    rng = np.random.default_rng(42)
    return 1000. * img + rng.normal(0, 1, img.shape), slits


@pytest.mark.parametrize('method', ['detector', 'slit', 'edge', 'leftedge', 'rightedge'])
def test_spat_flexure_shift_methods(method):
    # Different spatial flexure for the left and right edges of the slits
    dl, dr = 3.4, 2.6
    sciimg, slits = _synthetic_slit_image(dl, dr)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method=method,
                                              bpm=np.zeros(sciimg.shape, dtype=bool),
                                              maxlag=10, sigdetect=5.)
    assert spat_flexure.shape == (slits.nslits, 2), 'Spatial flexure has the wrong shape'
    expected = {'detector': None,
                'slit': [0.5*(dl+dr), 0.5*(dl+dr)],
                'edge': [dl, dr],
                'leftedge': [dl, dl],
                'rightedge': [dr, dr]}[method]
    if expected is None:
        # A single shift for all slits and edges, between the two edge shifts
        assert np.all(spat_flexure == spat_flexure[0, 0]), 'Detector flexure should be constant'
        assert dr - 0.5 < spat_flexure[0, 0] < dl + 0.5, 'Detector flexure is incorrect'
    else:
        assert np.allclose(spat_flexure, np.array(expected)[None, :], atol=0.15), \
            f'Incorrect spatial flexure for method={method}'


def test_spat_flexure_shift_bad_method():
    sciimg, slits = _synthetic_slit_image(1., 1.)
    with pytest.raises(PypeItError):
        flexure.spat_flexure_shift(sciimg, slits, method='bad')


def test_spat_flexure_shift_multislit():
    # With many slits, the cross-correlation of the slit edges has many
    # (genuine) secondary peaks; the detector-level flexure must still be found.
    sciimg, slits = _synthetic_slit_image(3., 3., nslit=10)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method='detector',
                                              bpm=np.zeros(sciimg.shape, dtype=bool),
                                              maxlag=10, sigdetect=5.)
    assert np.allclose(spat_flexure, 3., atol=0.1), 'Incorrect detector flexure for many slits'


@pytest.mark.parametrize('shift', [-12., 15.])
def test_spat_flexure_shift_beyond_maxlag(shift):
    # A flexure beyond the search window (+/- maxlag) cannot be measured, and
    # must not be replaced by a spurious secondary peak of the cross-correlation.
    sciimg, slits = _synthetic_slit_image(shift, shift)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method='detector',
                                              bpm=np.zeros(sciimg.shape, dtype=bool),
                                              maxlag=10, sigdetect=5.)
    assert np.all(spat_flexure == 0.), 'Flexure beyond maxlag should be rejected'
