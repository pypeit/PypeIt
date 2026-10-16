"""
Module to run tests on simple fitting routines for arrays
"""

import numpy as np

import pytest
from scipy.special import erf

from pypeit import PypeItError
from pypeit.core import flexure
from pypeit.slittrace import SlitTraceSet
from pypeit.spectrographs.util import load_spectrograph
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


@pytest.mark.parametrize('method', ['slit', 'edge', 'leftedge'])
def test_spat_flexure_shift_edge_near_detector(method):
    # The first slit starts close to the detector boundary, and the detector-level
    # flexure moves its left edge to within 2 pixels of the boundary for part of
    # the spectral range.  The edge refinement must not fail.
    # NOTE: An integer shift is used, because the cross-correlation of the (one
    # pixel wide) synthetic edge images is then at least 3 samples wide, which is
    # required to fit the detector-level peak.
    nspec, nspat, sig, dx = 200, 100, 1.2, -3.0
    drift = np.linspace(0., 3., nspec)
    left = np.stack([3. + drift, 50. + drift], axis=1)
    right = left + 30.
    slits = SlitTraceSet(left_init=left, right_init=right, pypeline='MultiSlit', nspat=nspat,
                         PYP_SPEC='dummy')
    xx = np.arange(nspat, dtype=float)[None, :]
    img = np.zeros((nspec, nspat))
    for slit in range(2):
        img += 0.5 * (erf((xx - left[:, slit, None] - dx) / (np.sqrt(2) * sig))
                      - erf((xx - right[:, slit, None] - dx) / (np.sqrt(2) * sig)))
    sciimg = 1000. * img + np.random.default_rng(42).normal(0, 1, img.shape)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method=method,
                                              bpm=np.zeros(sciimg.shape, dtype=bool),
                                              maxlag=10, sigdetect=5.)
    assert np.all(np.isfinite(spat_flexure)), 'Spatial flexure should be finite'
    assert np.allclose(spat_flexure, dx, atol=0.2), \
        f'Incorrect spatial flexure for method={method} with an edge near the detector boundary'


@pytest.mark.parametrize('method', ['detector', 'edge'])
def test_spat_flexure_shift_ignores_tweaked_edges(method):
    # The spatial flexure is measured relative to the initial (traced) slit
    # edges, which are the edges of the slits in the image.  Tweaked edges
    # (e.g., inset by a few pixels, based on the slit illumination profile)
    # must not change the result.
    dx = 3.
    sciimg, slits = _synthetic_slit_image(dx, dx)
    bpm = np.zeros(sciimg.shape, dtype=bool)
    reference = flexure.spat_flexure_shift(sciimg, slits, method=method, bpm=bpm, maxlag=10,
                                           sigdetect=5.)
    slits.left_tweak = slits.left_init + 3.2
    slits.right_tweak = slits.right_init - 1.0
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method=method, bpm=bpm, maxlag=10,
                                              sigdetect=5.)
    assert np.allclose(spat_flexure, dx, atol=0.15), \
        f'Incorrect spatial flexure for method={method} with tweaked slit edges'
    assert np.array_equal(spat_flexure, reference), 'Tweaked slit edges changed the spatial flexure'


def test_spat_flexure_shift_non_integer():
    # A non-integer detector-level shift must be measured (the cross-correlation
    # of the single-pixel edge images is then only 2 samples wide)
    dx = -2.6
    sciimg, slits = _synthetic_slit_image(dx, dx)
    bpm = np.zeros(sciimg.shape, dtype=bool)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method='detector', bpm=bpm, maxlag=10,
                                              sigdetect=5.)
    assert np.allclose(spat_flexure, dx, atol=0.1), 'Incorrect non-integer detector flexure'


def test_spat_flexure_shift_vertical_edges():
    # The refined flexure of each edge must be accurate, even if the slit edges
    # are aligned with the detector columns
    dl, dr = 3.4, 2.6
    sciimg, slits = _synthetic_slit_image(dl, dr, drift=0.)
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method='edge',
                                              bpm=np.zeros(sciimg.shape, dtype=bool),
                                              maxlag=10, sigdetect=5.)
    assert np.allclose(spat_flexure, np.array([dl, dr])[None, :], atol=0.02), \
        'Incorrect flexure for slit edges aligned with the detector columns'


def test_spat_flexure_shift_masked_slit():
    # A slit excluded from the flexure calculation is not refined; it is
    # assigned the detector level flexure, and flagged in the returned mask
    dx = 3.
    sciimg, slits = _synthetic_slit_image(dx + 0.3, dx + 0.3, nslit=3)
    slits.mask[1] = slits.bitmask.turn_on(slits.mask[1], 'BADREDUCE')
    bpm = np.zeros(sciimg.shape, dtype=bool)
    spat_flexure, gpm = flexure.spat_flexure_shift(sciimg, slits, method='edge', bpm=bpm,
                                                   maxlag=10, sigdetect=5., return_gpm=True)
    assert np.array_equal(gpm, np.array([[True, True], [False, False], [True, True]])), \
        'Incorrect good-pixel mask for the slit edges'
    assert np.allclose(spat_flexure[[0, 2]], dx + 0.3, atol=0.02), 'Incorrect refined flexure'
    det_flexure = flexure.spat_flexure_shift(sciimg, slits, method='detector', bpm=bpm,
                                             maxlag=10, sigdetect=5.)
    assert np.array_equal(spat_flexure[1], det_flexure[1]), \
        'Excluded slit should be assigned the detector level flexure'



@pytest.mark.parametrize('method', ['detector', 'edge'])
def test_spat_flexure_shift_padded_slits(method):
    # The spatial flexure is measured relative to the traced slit edges, even
    # if the slit edges are padded (see the ``pad`` parameter of EdgeTracePar)
    dl, dr = 3.4, 2.6
    sciimg, slits = _synthetic_slit_image(dl, dr)
    bpm = np.zeros(sciimg.shape, dtype=bool)
    reference = flexure.spat_flexure_shift(sciimg, slits, method=method, bpm=bpm, maxlag=10,
                                           sigdetect=5.)
    slits.pad = np.array([5., 3.])
    spat_flexure = flexure.spat_flexure_shift(sciimg, slits, method=method, bpm=bpm, maxlag=10,
                                              sigdetect=5.)
    assert np.array_equal(spat_flexure, reference), 'Padding changed the spatial flexure'


def test_spat_flexure_edge_match():
    # Traced edges (left = -1, right = 1) in all rows; the image shows the
    # edges shifted by 3 pixels, but only in the first 30 of 100 rows
    nspec, nspat = 100, 60
    slits_edges = np.zeros((nspec, nspat), dtype=int)
    slits_edges[:, 10] = -1
    slits_edges[:, 40] = 1
    sci_edges = np.zeros_like(slits_edges)
    sci_edges[:30, 13] = -1
    sci_edges[:30, 43] = 1
    assert np.isclose(flexure.spat_flexure_edge_match(sci_edges, slits_edges, 3.), 0.3)
    # Matched within the tolerance of 2 pixels
    assert np.isclose(flexure.spat_flexure_edge_match(sci_edges, slits_edges, 1.6), 0.3)
    assert flexure.spat_flexure_edge_match(sci_edges, slits_edges, 0.) == 0.
    # Edges of the wrong type (left vs. right) are not matched
    assert flexure.spat_flexure_edge_match(-sci_edges, slits_edges, 3.) == 0.
    # Traced edges shifted off the detector are not matched
    assert flexure.spat_flexure_edge_match(sci_edges, slits_edges, 100.) == 0.
    # No traced edges
    assert flexure.spat_flexure_edge_match(sci_edges, np.zeros_like(slits_edges), 3.) == 0.


@pytest.mark.parametrize('method', ['detector', 'edge'])
def test_spat_flexure_shift_undetected_edges(method):
    # The slits are only visible in the first nvis spectral rows (e.g. a faint
    # exposure), so few of the traced slit-edge pixels are matched
    dx = 3.
    sciimg, slits = _synthetic_slit_image(dx, dx, nslit=3)
    bpm = np.zeros(sciimg.shape, dtype=bool)
    rng = np.random.default_rng(1)
    for nvis, detected in [(10, False), (40, True)]:
        _sciimg = sciimg.copy()
        _sciimg[nvis:] = rng.normal(0, 1, _sciimg[nvis:].shape)
        spat_flexure, gpm = flexure.spat_flexure_shift(_sciimg, slits, method=method, bpm=bpm,
                                                       maxlag=10, sigdetect=5., minfrac=0.1,
                                                       return_gpm=True)
        if detected:
            # NOTE: The tolerance is larger than for the other tests, because
            # the slit edges drift by less than a pixel in the visible rows
            assert np.allclose(spat_flexure, dx, atol=0.3), 'Incorrect spatial flexure'
        else:
            assert np.all(spat_flexure == 0.), 'Undetected slit edges should give no flexure'
            assert not np.any(gpm), 'Undetected slit edges should be flagged'
    # The check can be disabled
    spat_flexure = flexure.spat_flexure_shift(_sciimg, slits, method='detector', bpm=bpm,
                                              maxlag=10, sigdetect=5., minfrac=0.)
    assert np.allclose(spat_flexure, dx, atol=0.3), 'Incorrect spatial flexure with minfrac=0'



def _fractional_slit_image(slits, spat_flexure, offslit=1e6):
    """
    Image with a value that depends linearly on the fractional position across
    each slit (and on the slit index), for the given spatial flexure.
    Off-slit pixels are set to ``offslit``.
    """
    left, right, _ = slits.select_edges(initial=True, spat_flexure=spat_flexure, pad=0)
    slitidx_img = slits.slit_img(pad=0, initial=True, spat_flexure=spat_flexure,
                                 use_spatial=False)
    img = np.full((slits.nspec, slits.nspat), offslit)
    for slit_idx in range(slits.nslits):
        spec, spat = np.where(slitidx_img == slit_idx)
        frac = (spat - left[spec, slit_idx]) / (right[spec, slit_idx] - left[spec, slit_idx])
        img[spec, spat] = 10. * slit_idx + 1. + frac
    return img


def test_spat_flexure_remap():
    _, slits = _synthetic_slit_image(0., 0., nslit=2)
    src_flexure = np.array([[1.3, 0.8], [1.3, 0.8]])
    flexure_ = np.array([[-0.6, -0.2], [2.1, 2.4]])
    src_img = _fractional_slit_image(slits, src_flexure)
    expected = _fractional_slit_image(slits, flexure_, offslit=0.)
    remap = flexure.spat_flexure_remap(src_img, slits, spat_flexure=flexure_,
                                       src_spat_flexure=src_flexure)
    onslit = slits.slit_img(pad=0, initial=True, spat_flexure=flexure_) != -1
    # Off-slit pixels of the source image are never used for on-slit pixels
    assert np.all(remap[onslit] < 100.), 'Off-slit source pixels used in the remapped image'
    # Each pixel is mapped to the same fractional position across the slit,
    # except for the pixels at the slit edges, where the source image is
    # limited to the pixels on the slit (i.e., the nearest source pixel on the
    # slit is used, which is less than one pixel, or ~0.035 in fractional
    # position, away)
    err = np.abs(remap - expected)
    exact = onslit & (err < 1e-10)
    assert np.sum(exact) > 0.85 * np.sum(onslit), 'Pixels should be remapped by fractional position'
    assert np.all(err[onslit] < 0.05), 'Incorrect remapped values at the slit edges'
    left, right, _ = slits.select_edges(initial=True, spat_flexure=flexure_, pad=0)
    spec, spat = np.where(onslit & ~exact)
    dist = np.min(np.abs(np.concatenate([left[spec], right[spec]], axis=1) - spat[:, None]), axis=1)
    assert np.all(dist < 2.), 'Only the pixels at the slit edges should be approximate'
    # Off-slit pixels keep the input values
    assert np.array_equal(remap[~onslit], src_img[~onslit]), 'Off-slit pixels should not change'
    # Identical slits and flexure: nothing to remap
    assert np.array_equal(flexure.spat_flexure_remap(src_img, slits, spat_flexure=src_flexure,
                                                     src_spat_flexure=src_flexure), src_img)
    # The flexure can also be included in the slit edges of the source frame
    src_slits = SlitTraceSet(left_init=slits.left_init + src_flexure[None, :, 0],
                             right_init=slits.right_init + src_flexure[None, :, 1],
                             pypeline='MultiSlit', nspat=slits.nspat, PYP_SPEC='dummy')
    assert np.allclose(flexure.spat_flexure_remap(src_img, slits, spat_flexure=flexure_,
                                                  src_slits=src_slits), remap)


def test_spat_flexure_remap_bad_input():
    sciimg, slits = _synthetic_slit_image(0., 0., nslit=2)
    _, slits3 = _synthetic_slit_image(0., 0., nslit=3)
    with pytest.raises(PypeItError):
        flexure.spat_flexure_remap(sciimg, slits, src_slits=slits3)
    with pytest.raises(PypeItError):
        flexure.spat_flexure_remap(sciimg, slits, spat_flexure=np.zeros((3, 2)))
    with pytest.raises(PypeItError):
        flexure.spat_flexure_remap(sciimg[:, :-1], slits)


def test_fit_spat_flexure():
    # A linear trend of the flexure with slit position, with one outlier and
    # one edge without a measured flexure
    slit_cen = np.array([10., 50., 90., 130., 170., 210.])
    true = np.stack([1. + 0.01 * slit_cen, 2. - 0.005 * slit_cen], axis=1)
    spat_flexure = true.copy()
    spat_flexure[3, 0] += 5.    # Outlier
    spat_flexure[4, 1] = 0.     # Not measured
    gpm = np.ones(spat_flexure.shape, dtype=bool)
    gpm[4, 1] = False
    fit = flexure.fit_spat_flexure(spat_flexure, slit_cen, 1, gpm=gpm)
    assert np.allclose(fit, true), 'Incorrect polynomial fit to the spatial flexure'
    # The spectrograph method uses the same fit
    spec = load_spectrograph('shane_kast_blue')
    assert np.allclose(spec.spatial_flexure(spat_flexure, polyord=1, slit_cen=slit_cen, gpm=gpm), true)
    # Too few edges to fit: the measured values are returned
    gpm[:, 0] = False
    fit = flexure.fit_spat_flexure(spat_flexure, slit_cen, 1, gpm=gpm)
    assert np.array_equal(fit[:, 0], spat_flexure[:, 0]), 'Unfit edges should not be changed'
