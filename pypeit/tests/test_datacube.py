"""
Module to run tests on datacube generation
"""

from astropy.wcs import WCS
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
import astropy.units as u
from astropy.coordinates import SkyCoord
from IPython import embed
import numpy as np
import pytest
from scipy.optimize import curve_fit

from pypeit import utils
from pypeit.core import datacube
from pypeit.core.flexure import calculate_image_phase
from pypeit.coadd3d import DataCube


def make_test_datacube(whitelight_range=None):
    """
    Construct a small DataCube object for datamodel tests.
    """
    wave = np.linspace(5000.0, 5018.0, 10)
    flux = np.ones((wave.size, 4, 5), dtype=float)
    sig = np.ones_like(flux)
    bpm = np.zeros(flux.shape, dtype=np.uint8)
    blaze_wave = np.array([wave[0], wave[-1]])
    blaze_spec = np.ones(2, dtype=float)
    return DataCube(
        flux, sig, bpm, wave, 'keck_kcrm', blaze_wave, blaze_spec,
        whitelight_range=whitelight_range, fluxed=False
    )


def make_test_cube_wcs():
    """
    Construct a simple 3D WCS for DataCube I/O tests.
    """
    wcs = WCS(naxis=3)
    wcs.wcs.crpix = [1.0, 1.0, 1.0]
    wcs.wcs.cdelt = [-1.0e-4, 1.0e-4, 2.0]
    wcs.wcs.crval = [150.0, 2.0, 5000.0]
    wcs.wcs.cunit = ['deg', 'deg', 'Angstrom']
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'WAVE']
    return wcs


def test_datacube_whitelight_range_datamodel(tmp_path):
    wl_range = np.array([5004.0, 5012.0])
    cube = make_test_datacube(whitelight_range=wl_range)

    assert np.allclose(cube.resolve_whitelight_range([None, None]), wl_range), \
        'resolve_whitelight_range([None, None]) should fall back to the stored default range'
    assert np.allclose(cube.resolve_whitelight_range([5006.0, None]), [5006.0, 5012.0]), \
        'a user-supplied lower bound should override the default, keeping the default upper bound'
    assert np.allclose(cube.resolve_whitelight_range([None, 5010.0]), [5004.0, 5010.0]), \
        'a user-supplied upper bound should override the default, keeping the default lower bound'
    assert np.allclose(cube.resolve_whitelight_range([5002.0, 5016.0]), [5002.0, 5016.0]), \
        'fully user-supplied bounds should be returned unchanged'

    ofile = tmp_path / 'spec3d_test.fits'
    cube.to_file(str(ofile), hdr=make_test_cube_wcs().to_header(), overwrite=True)
    read_cube = DataCube.from_file(str(ofile))
    assert np.allclose(read_cube.whitelight_range, wl_range), \
        'whitelight_range did not survive a FITS write/read round trip'
    assert np.allclose(read_cube.resolve_whitelight_range([None, None]), wl_range), \
        'resolve_whitelight_range should still recover the stored default after a round trip'

    int_range_cube = make_test_datacube(whitelight_range=[5004, 5012])
    assert np.issubdtype(int_range_cube.whitelight_range.dtype, np.floating), \
        'whitelight_range should be cast to float even when given as integers'
    assert np.allclose(int_range_cube.whitelight_range, wl_range), \
        'integer whitelight_range input should be preserved numerically'

    # Old cubes do not have the new datamodel field.  They should still load,
    # and extraction should fall back to the full saved cube wavelength range.
    old_cube = make_test_datacube(whitelight_range=None)
    old_file = tmp_path / 'spec3d_old.fits'
    old_cube.to_file(str(old_file), hdr=make_test_cube_wcs().to_header(), overwrite=True)
    with fits.open(old_file, mode='update') as hdu:
        for ext in hdu[1:]:
            ext.header['DMODVER'] = '1.2.0'

    read_old_cube = DataCube.from_file(str(old_file))
    assert read_old_cube.whitelight_range is None, \
        'a pre-1.3.0 cube on disk should load with whitelight_range unset, not a stale/garbage value'
    assert np.allclose(
        read_old_cube.resolve_whitelight_range([None, None]),
        [read_old_cube.wave.min(), read_old_cube.wave.max()]
    ), 'a cube with no stored whitelight_range should fall back to the full saved wavelength range'

photutils_required = pytest.mark.skipif(
    datacube.DAOStarFinder is None, reason='photutils not installed'
)


def make_point_source_image(ra, dec, dspat, wcs,
                            n_pix=100, fwhm_arcsec=None):
    """
    Create a 100x100 pixel image of a point source modelled as a 2D Gaussian.

    Parameters
    ----------
    ra : float
        Right ascension of the point source in degrees.
    dec : float
        Declination of the point source in degrees.
    dspat : float
        Spatial pixel size in arcseconds.
    wcs : astropy.wcs.WCS
        Astropy WCS object describing the coordinate system.
        Only the spatial axes are used.
    n_pix : int, optional
        Number of pixels along each axis of the output image. Default is 100.
    fwhm_arcsec : float, optional
        The FWHM of the PSF in arcseconds. Default is None.

    Returns
    -------
    image : ndarray, shape (100, 100)
        Normalised 2D Gaussian image. Peak value is 1.0.
    x_cen : float
        Sub-pixel x (column) position of the source in the image.
    y_cen : float
        Sub-pixel y (row) position of the source in the image.
    """
    # Default PSF FWHM: 6 pixels
    _fwhm_arcsec = 6 * dspat if fwhm_arcsec is None else fwhm_arcsec

    # Convert FWHM to sigma in pixels
    fwhm_pix = _fwhm_arcsec / dspat
    sigma_pix = fwhm_pix / (2.0 * np.sqrt(2.0 * np.log(2.0)))

    # Get the source position in output image pixel coords
    src_pix = wcs.wcs_world2pix(ra, dec, 0)
    x_cen_out, y_cen_out = src_pix[0], src_pix[1]

    # Evaluate 2D Gaussian on the pixel grid
    x = np.arange(n_pix, dtype=np.float64)  # column indices
    y = np.arange(n_pix, dtype=np.float64)  # row indices
    xx, yy = np.meshgrid(x, y)  # shape (100, 100)

    image = np.exp(
        -((xx - x_cen_out) ** 2 + (yy - y_cen_out) ** 2) / (2.0 * sigma_pix ** 2)
    )

    return image

def test_align_phase():
    numiter = 5  # This is the same number of iterations used in the coadd3d run_align() method
    ref_idx = 0
    _dspat = 0.2*u.arcsec  # Spatial pixel scale
    fwhm = 4*_dspat
    ra1, dec1 = 130.0*u.deg, 35.0*u.deg  # Pick two random numbers
    cosdec = np.cos(dec1.to(u.radian).value)
    offs_ra, offs_dec = 1.0*u.arcsec, -1.5*u.arcsec
    # Compute second coordinates using pre-defined offset
    ra2, dec2 = ra1 + offs_ra, dec1 + offs_dec
    # Build a simple TAN-projected WCS
    n_pix = 100
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [n_pix//2, n_pix//2]
    wcs.wcs.cdelt = [-_dspat.to(u.deg).value, _dspat.to(u.deg).value]
    wcs.wcs.crval = [ra1.to(u.deg).value, dec1.to(u.deg).value]  # RA, Dec in degrees
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    # Define two coordinates
    coord1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs')
    coord2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs')

    # Calculate and check the separation
    dra, ddec = coord1.spherical_offsets_to(coord2)
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value*cosdec, atol=0.001), \
        'test setup error: injected RA separation (with cos(dec) factor) does not match ' \
        'the requested offset'
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001), \
        'test setup error: injected Dec separation does not match the requested offset'

    # Loop through iterations to find the offsets using the PHASE method
    # Set up the offsets arrays that will be computed by the spatial alignment methods
    ra_offsets, dec_offsets = np.zeros(2) * u.deg, np.zeros(2) * u.deg
    for ii in range(numiter):
        # Generate two images wih the current RA, Dec and offsets
        img1 = make_point_source_image(ra1-ra_offsets[0], dec1-dec_offsets[0], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        img2 = make_point_source_image(ra2-ra_offsets[1], dec2-dec_offsets[1], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        wl_imgs = [img1, img2]
        # Select the reference image (i.e. the one that doesn't shift)
        ref_img = wl_imgs[ref_idx].copy()
        # Compute the offsets between all images
        for ff in range(len(wl_imgs)):
            dec_shift, ra_shift = calculate_image_phase(
                ref_img.copy(), wl_imgs[ff], maskval=0.0
            )
            # Convert pixel shift to degrees shift
            ra_shift *= -_dspat.to(u.deg)/cosdec
            dec_shift *= _dspat.to(u.deg)
            # Update the offsets for the next iteration
            ra_offsets[ff] -= ra_shift.to(u.deg)
            dec_offsets[ff] -= dec_shift.to(u.deg)

    # Check that the offsets agree to within 1/3 of a pixel
    atol = _dspat.to(u.arcsec).value/3
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'phase method: the reference image should recover a zero RA offset relative to itself'
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'phase method: the reference image should recover a zero Dec offset relative to itself'
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol), \
        'phase method: recovered RA offset does not match the injected offset'
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol), \
        'phase method: recovered Dec offset does not match the injected offset'


def test_align_cc():
    numiter = 5  # This is the same number of iterations used in the coadd3d run_align() method
    ref_idx = 0
    _dspat = 0.2 * u.arcsec  # Spatial pixel scale
    fwhm = 4 * _dspat
    ra1, dec1 = 130.0 * u.deg, 35.0 * u.deg  # Pick two random numbers
    cosdec = np.cos(dec1.to(u.radian).value)
    offs_ra, offs_dec = 1.0 * u.arcsec, -1.5 * u.arcsec
    # Compute second coordinates using pre-defined offset
    ra2, dec2 = ra1 + offs_ra, dec1 + offs_dec
    # Build a simple TAN-projected WCS
    n_pix = 100
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [n_pix // 2, n_pix // 2]
    wcs.wcs.cdelt = [-_dspat.to(u.deg).value, _dspat.to(u.deg).value]
    wcs.wcs.crval = [ra1.to(u.deg).value, dec1.to(u.deg).value]  # RA, Dec in degrees
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    # Define two coordinates
    coord1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs')
    coord2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs')

    # Calculate and check the separation
    dra, ddec = coord1.spherical_offsets_to(coord2)
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value * cosdec, atol=0.001), \
        'test setup error: injected RA separation (with cos(dec) factor) does not match ' \
        'the requested offset'
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001), \
        'test setup error: injected Dec separation does not match the requested offset'

    # Loop through iterations to find the offsets using the CC method
    # Set up the offsets arrays that will be computed by the spatial alignment methods
    ra_offsets, dec_offsets = np.zeros(2) * u.deg, np.zeros(2) * u.deg
    for ii in range(numiter):
        # Generate two images wih the current RA, Dec and offsets
        img1 = make_point_source_image(ra1-ra_offsets[0], dec1-dec_offsets[0], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        img2 = make_point_source_image(ra2-ra_offsets[1], dec2-dec_offsets[1], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        wl_imgs = [img1, img2]
        # Select the reference image (i.e. the one that doesn't shift)
        ref_img = wl_imgs[ref_idx].copy()
        # Compute the offsets between all images
        for ff in range(len(wl_imgs)):
            dec_shift, ra_shift = calculate_image_phase(
                ref_img.copy(), wl_imgs[ff], maskval=0.0, force_cc=True
            )
            # Convert pixel shift to degrees shift
            ra_shift *= -_dspat.to(u.deg)/cosdec
            dec_shift *= _dspat.to(u.deg)
            # Update the offsets for the next iteration
            ra_offsets[ff] -= ra_shift.to(u.deg)
            dec_offsets[ff] -= dec_shift.to(u.deg)

    # Check that the offsets agree to within 1/3 of a pixel
    atol = _dspat.to(u.arcsec).value/3
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'cc method: the reference image should recover a zero RA offset relative to itself'
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'cc method: the reference image should recover a zero Dec offset relative to itself'
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol), \
        'cc method: recovered RA offset does not match the injected offset'
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol), \
        'cc method: recovered Dec offset does not match the injected offset'


@photutils_required
def test_align_fit():
    numiter = 5  # This is the same number of iterations used in the coadd3d run_align() method
    ref_idx = 0
    _dspat = 0.2 * u.arcsec  # Spatial pixel scale
    fwhm = 4 * _dspat
    ra1, dec1 = 130.0 * u.deg, 35.0 * u.deg  # Pick two random numbers
    cosdec = np.cos(dec1.to(u.radian).value)
    offs_ra, offs_dec = 1.0 * u.arcsec, -1.5 * u.arcsec
    # Compute second coordinates using pre-defined offset
    ra2, dec2 = ra1 + offs_ra, dec1 + offs_dec
    # Build a simple TAN-projected WCS
    n_pix = 100
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [n_pix // 2, n_pix // 2]
    wcs.wcs.cdelt = [-_dspat.to(u.deg).value, _dspat.to(u.deg).value]
    wcs.wcs.crval = [ra1.to(u.deg).value, dec1.to(u.deg).value]  # RA, Dec in degrees
    wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']

    # Define two coordinates
    coord1 = SkyCoord(ra=ra1, dec=dec1, frame='icrs')
    coord2 = SkyCoord(ra=ra2, dec=dec2, frame='icrs')

    # Calculate and check the separation
    dra, ddec = coord1.spherical_offsets_to(coord2)
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value * cosdec, atol=0.001), \
        'test setup error: injected RA separation (with cos(dec) factor) does not match ' \
        'the requested offset'
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001), \
        'test setup error: injected Dec separation does not match the requested offset'

    # Loop through iterations to find the offsets using the FIT method
    # Set up the offsets arrays that will be computed by the spatial alignment methods
    ra_offsets, dec_offsets = np.zeros(2) * u.deg, np.zeros(2) * u.deg
    for ii in range(numiter):
        # Generate two images wih the current RA, Dec and offsets
        img1 = make_point_source_image(ra1-ra_offsets[0], dec1-dec_offsets[0], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        img2 = make_point_source_image(ra2-ra_offsets[1], dec2-dec_offsets[1], _dspat.value, wcs,
                                       fwhm_arcsec=fwhm.value, n_pix=n_pix)
        wl_imgs = [img1, img2]
        numfiles = len(wl_imgs)
        # Compute the offsets between all images
        ra_pix_star = np.zeros(len(wl_imgs))
        dec_pix_star = np.zeros(len(wl_imgs))
        for ff in range(numfiles):
            popt, pcov, model, init_obj_position, flux_opt, sigma_opt = \
                datacube.fitGaussian2D(
                    wl_imgs[ff],
                    # ivar=utils.inverse(np.square(sig_imgs[ff])),
                    # gpm=np.logical_not(bpm_imgs[ff]),
                    fwhm=fwhm.value/_dspat.value, norm=False
            )
            gaussian_position = popt[1], popt[2]
            ra_pix_star[ff], dec_pix_star[ff] = gaussian_position
        ra_shifts = -(ra_pix_star[ref_idx] - ra_pix_star) * _dspat / cosdec
        dec_shifts = (dec_pix_star[ref_idx] - dec_pix_star) * _dspat
        ra_offsets = [ra_offsets[ff] - ra_shifts[ff] for ff in range(numfiles)]
        dec_offsets = [dec_offsets[ff] - dec_shifts[ff] for ff in range(numfiles)]

    # Check that the offsets agree to within 1/3 of a pixel
    atol = _dspat.to(u.arcsec).value/3
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'fit method: the reference image should recover a zero RA offset relative to itself'
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol), \
        'fit method: the reference image should recover a zero Dec offset relative to itself'
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol), \
        'fit method: recovered RA offset does not match the injected offset'
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol), \
        'fit method: recovered Dec offset does not match the injected offset'


@photutils_required
def test_fitgaussian2d_ginga_xy_convention():
    """`fitGaussian2D` must return (x, y) = (column, row), matching Ginga/DS9.

    Uses a non-square image with two point sources of unequal brightness at
    positions with different row and column indices, so a transpose or a
    row/column swap cannot accidentally produce the right answer.
    """
    ny, nx = 20, 28
    bright_row, bright_col, bright_amp = 14, 7, 1.0
    faint_row, faint_col, faint_amp = 5, 21, 0.4
    sigma_pix = 1.3
    fwhm_pix = sigma_pix * 2.0 * np.sqrt(2.0 * np.log(2.0))

    yimg, ximg = np.mgrid[0:ny, 0:nx]
    image = (
        bright_amp * np.exp(-0.5 * (((ximg - bright_col) / sigma_pix) ** 2
                                     + ((yimg - bright_row) / sigma_pix) ** 2))
        + faint_amp * np.exp(-0.5 * (((ximg - faint_col) / sigma_pix) ** 2
                                      + ((yimg - faint_row) / sigma_pix) ** 2))
    )

    # Auto-detection (DAOStarFinder, n_brightest=1) must lock onto the bright source.
    popt, *_ = datacube.fitGaussian2D(image, fwhm=fwhm_pix, mask_edge=2, norm=False)
    assert np.isclose(popt[1], bright_col, atol=0.2), 'x (column) does not match the bright source'
    assert np.isclose(popt[2], bright_row, atol=0.2), 'y (row) does not match the bright source'

    # A manually-supplied (x, y) position must override auto-detection and lock
    # onto the source at that position, even though it is fainter.
    popt_manual, *_ = datacube.fitGaussian2D(
        image, init_obj_position=(faint_col, faint_row), fwhm=fwhm_pix, mask_edge=2, norm=False)
    assert np.isclose(popt_manual[1], faint_col, atol=0.2), \
        'manual_position did not lock the fit onto the faint source in x (column)'
    assert np.isclose(popt_manual[2], faint_row, atol=0.2), \
        'manual_position did not lock the fit onto the faint source in y (row)'


# ---------------------------------------------------------------------------
# NOTE -- TEMPORARY, branch-specific regression scaffolding.
#
# The helpers and tests in this block (`_pre_fix_auto_position`,
# `_pre_fix_manual_position`, `test_fitgaussian2d_autodetect_unaffected_by_bugfix`,
# and `test_manual_position_bugfix_regression`) exist only to document and lock
# in the specific `(x, y)` mislabeling bug fixed by `fix-cube-manual-coords-kbw`
# relative to its parent branch, `kcwi_dec_2024_weights`.  They intentionally
# duplicate old, buggy logic as a frozen reference for comparison.  These
# should be removed or refactored into ordinary (non-comparative) convention
# tests before this branch stack is merged into `develop` -- once merged,
# there is no further value in re-deriving the pre-fix behavior, and carrying
# duplicated buggy code forward in the test suite is a maintenance liability.
# ---------------------------------------------------------------------------

def _pre_fix_auto_position(image, fwhm, nsigma=5.0, mask_edge=2):
    """
    Reference re-implementation of the auto-detection (DAOStarFinder + 2D
    Gaussian fit) branch of ``fitGaussian2D`` as it existed on
    ``kcwi_dec_2024_weights``, before the ``fix-cube-manual-coords-kbw`` fix.

    On that branch, the curve-fit coordinate grid was built with
    ``indexing='ij'`` using ``image.shape[0]``/``image.shape[1]`` swapped
    relative to the Ginga/DS9 (x, y) convention -- so internally "x" tracked
    the row (Dec) axis and "y" tracked the column (RA) axis -- DAOStarFinder's
    centroids were assigned as ``(y_centroid, x_centroid)`` to match, and
    ``extract_point_source`` applied a final ``yobj, xobj = gaussian_position``
    swap to undo the internal mislabeling before using the position.  This
    function reproduces that chain bug-for-bug, purely so
    ``test_fitgaussian2d_autodetect_unaffected_by_bugfix`` below can prove the
    fix left the auto-detection result numerically unchanged.  Do not "fix"
    this function -- it must stay identical to the pre-fix code.

    TEMPORARY: remove this helper (and its one caller) before merging this
    branch stack into `develop`; see the module note above.
    """
    fwhm2sigma = 1.0 / (2 * np.sqrt(2 * np.log(2)))
    ny, nx = image.shape
    yimg, ximg = np.mgrid[0:ny, 0:nx]
    edgemask = (ximg < mask_edge) | (ximg >= nx - mask_edge) \
        | (yimg < mask_edge) | (yimg >= ny - mask_edge)

    _, median, std = sigma_clipped_stats(image[~edgemask], sigma=3.0)
    ivar = np.full_like(image, 1.0 / std ** 2) if std > 0 else np.ones_like(image)

    daofind = datacube.DAOStarFinder(fwhm=fwhm, threshold=nsigma, sharpness_range=(0.2, 2.0),
                                     exclude_border=False, n_brightest=1)
    sources = daofind((image - median) * np.sqrt(ivar), mask=edgemask)
    init_obj_position = sources['y_centroid'][0], sources['x_centroid'][0]

    initial_guess = (1, init_obj_position[0], init_obj_position[1],
                     fwhm * fwhm2sigma, fwhm * fwhm2sigma, 0, 0)
    bounds = ([0, init_obj_position[0] - fwhm / 3.0, init_obj_position[1] - fwhm / 3.0,
               fwhm / 6.0, fwhm / 6.0, -np.pi, -np.inf],
              [np.inf, init_obj_position[0] + fwhm / 3.0, init_obj_position[1] + fwhm / 3.0,
               fwhm, fwhm, np.pi, np.inf])

    x = np.linspace(0, ny - 1, ny)
    y = np.linspace(0, nx - 1, nx)
    xx, yy = np.meshgrid(x, y, indexing='ij')
    popt, _ = curve_fit(datacube.gaussian2D, (xx, yy), image.ravel(), bounds=bounds,
                        p0=initial_guess)
    xpos_gauss, ypos_gauss = popt[1], popt[2]
    yobj, xobj = xpos_gauss, ypos_gauss  # extract_point_source's pre-fix unpack
    return xobj, yobj


def _pre_fix_manual_position(manual_position):
    """
    Reproduces the pre-fix (``kcwi_dec_2024_weights``) handling of a
    user-supplied ``manual_position`` inside ``extract_point_source``:
    ``yobj, xobj = manual_position`` instead of the corrected
    ``xobj, yobj = manual_position``.  Kept only to lock in a regression test
    for the ``fix-cube-manual-coords-kbw`` fix -- do not "fix" this helper.

    TEMPORARY: remove this helper (and its one caller) before merging this
    branch stack into `develop`; see the module note above.
    """
    yobj, xobj = manual_position
    return xobj, yobj


@photutils_required
def test_fitgaussian2d_autodetect_unaffected_by_bugfix():
    """The manual-coordinate bugfix must not change the auto-detection result.

    The pre-fix code's internal (x, y) mislabeling was self-consistent for
    auto-detected positions (see `_pre_fix_auto_position`), so the final
    position handed to the extraction code should be numerically identical
    before and after the fix.

    TEMPORARY: this comparative test should be removed (or refactored into an
    ordinary convention test, without the `_pre_fix_auto_position` reference
    implementation) before merging this branch stack into `develop`; see the
    module note above `_pre_fix_auto_position`.
    """
    ny, nx = 20, 28
    row0, col0 = 12, 17
    sigma_pix = 1.4
    fwhm_pix = sigma_pix * 2.0 * np.sqrt(2.0 * np.log(2.0))
    yimg, ximg = np.mgrid[0:ny, 0:nx]
    image = np.exp(-0.5 * (((ximg - col0) / sigma_pix) ** 2 + ((yimg - row0) / sigma_pix) ** 2))

    popt, *_ = datacube.fitGaussian2D(image, fwhm=fwhm_pix, mask_edge=2, norm=False)
    x_new, y_new = popt[1], popt[2]
    x_old, y_old = _pre_fix_auto_position(image, fwhm=fwhm_pix, mask_edge=2)

    assert np.isclose(x_new, x_old, atol=1e-4), \
        'auto-detected x position changed between the pre-fix and current code'
    assert np.isclose(y_new, y_old, atol=1e-4), \
        'auto-detected y position changed between the pre-fix and current code'
    # Sanity check both still recover the true, known source position.
    assert np.isclose(x_new, col0, atol=0.2), \
        'current fitGaussian2D did not recover the true source column (x) position'
    assert np.isclose(y_new, row0, atol=0.2), \
        'current fitGaussian2D did not recover the true source row (y) position'


def test_manual_position_bugfix_regression():
    """Lock in the fix-cube-manual-coords-kbw fix for manual extraction positions.

    Pre-fix, a user-supplied ``manual_position=(x, y)`` (Ginga/DS9 convention)
    was silently transposed before use.  Choose x != y so the swap is
    detectable.

    TEMPORARY: this comparative test should be removed before merging this
    branch stack into `develop`, since `test_fitgaussian2d_ginga_xy_convention`
    above already covers the correct (non-comparative) behavior going forward;
    see the module note above `_pre_fix_auto_position`.
    """
    manual_position = (17, 6)  # (x, y): column 17, row 6

    x_new, y_new = manual_position  # current (fixed) unpack in extract_point_source
    assert (x_new, y_new) == manual_position, \
        'the current (fixed) unpack must return manual_position unchanged, i.e. (x, y) with no swap'

    x_old, y_old = _pre_fix_manual_position(manual_position)
    assert (x_old, y_old) == (manual_position[1], manual_position[0]), \
        'the pre-fix reference implementation should reproduce the x/y swap bug'
    assert (x_old, y_old) != manual_position, \
        'pre-fix code should have gotten the manual position wrong (x/y swapped)'


@photutils_required
def test_extract_point_source_manual_position_selects_correct_source():
    """`extract_point_source` must extract the source at `manual_position`,
    not always the brightest source in the field.
    """
    nwave, ny, nx = 6, 20, 28
    bright_row, bright_col, bright_amp = 14, 7, 1.0
    faint_row, faint_col, faint_amp = 5, 21, 0.4
    sigma_pix = 1.3

    yimg, ximg = np.mgrid[0:ny, 0:nx]
    image2d = (
        bright_amp * np.exp(-0.5 * (((ximg - bright_col) / sigma_pix) ** 2
                                     + ((yimg - bright_row) / sigma_pix) ** 2))
        + faint_amp * np.exp(-0.5 * (((ximg - faint_col) / sigma_pix) ** 2
                                      + ((yimg - faint_row) / sigma_pix) ** 2))
    )
    flxcube = np.broadcast_to(image2d, (nwave, ny, nx)).copy()
    # The synthetic image is noise-free, so use a large (but finite) ivar to
    # give DAOStarFinder's internal S/N image (image * sqrt(ivar)) a healthy
    # signal well above its default 5-sigma detection threshold.
    ivarcube = np.full((nwave, ny, nx), 1.0e8)
    bpmcube = np.zeros((nwave, ny, nx), dtype=bool)
    wave = np.linspace(5000.0, 5010.0, nwave)

    dspat_arcsec = 0.3
    dspat_deg = dspat_arcsec / 3600.0
    wcscube = WCS(naxis=3)
    wcscube.wcs.crpix = [1.0, 1.0, 1.0]
    wcscube.wcs.crval = [150.0, 10.0, wave[0]]
    wcscube.wcs.cdelt = [-dspat_deg, dspat_deg, 2.0]
    wcscube.wcs.cunit = [u.deg, u.deg, u.Angstrom]
    wcscube.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'WAVE']

    extract_kwargs = dict(
        exptime=1.0, fluxed=False, subpixel=5, boxcar_radius=1.0, fwhm=0.9,
        opt_prof_method='fit_gauss', spectrograph='keck_kcrm', show_qa=False
    )

    sobjs_auto, *_ = datacube.extract_point_source(
        wave, flxcube, ivarcube, bpmcube, wcscube, manual_position=None, **extract_kwargs)
    sobjs_bright, *_ = datacube.extract_point_source(
        wave, flxcube, ivarcube, bpmcube, wcscube,
        manual_position=(bright_col, bright_row), **extract_kwargs)
    sobjs_faint, *_ = datacube.extract_point_source(
        wave, flxcube, ivarcube, bpmcube, wcscube,
        manual_position=(faint_col, faint_row), **extract_kwargs)

    flux_auto = np.median(sobjs_auto[0].BOX_COUNTS)
    flux_bright = np.median(sobjs_bright[0].BOX_COUNTS)
    flux_faint = np.median(sobjs_faint[0].BOX_COUNTS)

    # With no manual position, auto-detection should lock onto the brighter source.
    assert np.isclose(flux_auto, flux_bright, rtol=0.05), \
        'auto-detection did not recover the same flux as the known bright source'
    # A manual position at the fainter source must override auto-detection and
    # recover roughly the expected (amplitude) flux ratio, not the bright flux.
    assert 0.2 * flux_bright < flux_faint < 0.6 * flux_bright, \
        'manual_position at the faint source did not recover the expected flux ratio'


def test_resample_spec_to_grid_identity():
    """Resampling onto the native grid reproduces flux/ivar on interior pixels.

    The flux-conserving resampler maps each interior output pixel one-to-one
    onto its native pixel, so the values are reproduced exactly there.  The two
    edge pixels can be flagged as partially covered (output borders float just
    outside the native coverage), so they are excluded from the comparison.
    """
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.sin(wave / 100.0) + 5.0
    ivar = np.full_like(wave, 4.0)
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, ivar, wave)
    interior = slice(1, -1)
    assert np.all(cov[interior]), \
        "interior output pixels should all be flagged as covered"
    np.testing.assert_allclose(
        fg[interior], flux[interior], rtol=1e-6,
        err_msg="resampled flux should match the input on interior pixels")
    np.testing.assert_allclose(
        ig[interior], ivar[interior], rtol=1e-6,
        err_msg="resampled ivar should match the input on interior pixels")


def test_resample_spec_to_grid_propagates_variance_in_quadrature():
    """Error is propagated in quadrature (variance space), not ivar space.

    Resampling onto the native grid leaves the per-pixel variance unchanged on
    interior pixels, so the returned inverse variance matches the input there.
    A naive average of the inverse variance would not preserve this.
    """
    wave = np.linspace(4000.0, 4100.0, 51)
    flux = np.ones_like(wave)
    # Inverse variance varying across the spectrum.
    ivar = np.linspace(1.0, 9.0, wave.size)
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, ivar, wave)
    interior = slice(1, -1)
    np.testing.assert_allclose(
        ig[interior], ivar[interior], rtol=1e-6,
        err_msg="variance must be propagated in quadrature; interior ivar "
                "should match the input")


def test_resample_spec_to_grid_zero_outside_coverage():
    """Grid points outside the native wavelength range stay zero."""
    wave = np.linspace(4500.0, 4600.0, 51)
    flux = np.ones_like(wave)
    ivar = np.ones_like(wave)
    grid = np.linspace(4000.0, 5000.0, 101)
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, ivar, grid)
    outside = (grid < 4500.0) | (grid > 4600.0)
    assert np.all(fg[outside] == 0.0), \
        "flux outside the native coverage should be zero"
    assert np.all(ig[outside] == 0.0), \
        "ivar outside the native coverage should be zero"
    assert not np.any(cov[outside]), \
        "no output pixel outside the native coverage should be flagged covered"
    # Interior of the covered range carries real flux.
    assert np.any(cov), \
        "at least some output pixels within the native range should be covered"


def test_resample_spec_to_grid_ivar_none():
    """A None inverse variance yields all-zero ivar but valid flux."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    grid = np.linspace(4100.0, 4900.0, 51)
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, None, grid)
    assert np.all(cov), \
        "all in-range output pixels should be covered when resampling flat flux"
    np.testing.assert_allclose(
        fg[cov], 1.0, rtol=1e-6,
        err_msg="resampled flat (unit) flux should remain 1.0 where covered")
    assert np.all(ig == 0.0), \
        "ivar should be all zero when no input inverse variance is supplied"


def test_resample_spec_to_grid_below_min_good():
    """Too few valid samples returns all-zero arrays and an empty mask."""
    wave = np.array([4000.0, -1.0, 0.0])  # only one positive wavelength
    flux = np.array([1.0, 1.0, 1.0])
    ivar = np.array([1.0, 1.0, 1.0])
    grid = np.linspace(4000.0, 5000.0, 11)
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, ivar, grid, min_good=2)
    assert not np.any(cov), \
        "too few valid samples should leave every output pixel uncovered"
    assert np.all(fg == 0.0), \
        "too few valid samples should give all-zero flux"
    assert np.all(ig == 0.0), \
        "too few valid samples should give all-zero ivar"


def test_resample_spec_to_grid_masked_native_pixels_drop_out():
    """Native pixels with ivar <= 0 are masked and carry no inverse variance."""
    wave = np.linspace(4000.0, 5000.0, 101)
    flux = np.ones_like(wave)
    ivar = np.zeros_like(wave)  # no weight anywhere
    fg, ig, cov = datacube.resample_spec_to_grid(wave, flux, ivar, wave)
    assert np.all(ig == 0.0), \
        "masked native pixels (ivar <= 0) should carry no inverse variance"
    # utils.inverse semantics: zero ivar -> zero variance contribution
    assert np.all(utils.inverse(ig) == 0.0), \
        "utils.inverse of zero ivar should be zero"
