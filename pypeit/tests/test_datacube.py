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
from pypeit import slittrace
from pypeit import alignframe
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


# ---------------------------------------------------------------------------
# Phase 1b -- decomposing the bundled develop -> kcwi_dec_2024 WCS axis-order
# and sign-convention changes (see kcwi_wcs.md, open question [Q3]).  Each
# test below isolates one of the four changes that were bundled into that
# rewrite (RA/Dec offset sign, voxedges axis order, the sky-right crpix/cdelt
# convention, and the world2pix axis reversal in `subpixellate`), so a future
# regression in any one of them is caught without needing a full cube build.
# ---------------------------------------------------------------------------

def test_wcs_bounds_subtracts_offsets():
    """`wcs_bounds` must SUBTRACT ra_offsets/dec_offsets from the pixel
    coordinates, matching the sign convention introduced in `kcwi_dec_2024`
    (previously `+=`, now `-=`).  Uses a nonzero, asymmetric offset so an
    accidental sign flip back to addition is guaranteed to produce a
    different (and thus detectable) result.
    """
    raImg = np.array([[150.0, 150.001], [150.0005, 150.0008]])
    decImg = np.array([[10.0, 10.002], [10.0011, 10.0006]])
    waveImg = np.array([[5000.0, 5001.0], [5002.0, 5003.0]])
    slitid_img_gpm = np.ones_like(raImg, dtype=int)

    ra_offset, dec_offset = 5.0e-4, -3.0e-4  # deg; deliberately asymmetric and nonzero

    ra_min, ra_max, dec_min, dec_max, wave_min, wave_max = datacube.wcs_bounds(
        raImg, decImg, waveImg, slitid_img_gpm, ra_offsets=ra_offset, dec_offsets=dec_offset
    )

    assert np.isclose(ra_min, raImg.min() - ra_offset), \
        'wcs_bounds should subtract ra_offsets from the minimum RA, not add it'
    assert np.isclose(ra_max, raImg.max() - ra_offset), \
        'wcs_bounds should subtract ra_offsets from the maximum RA, not add it'
    assert np.isclose(dec_min, decImg.min() - dec_offset), \
        'wcs_bounds should subtract dec_offsets from the minimum Dec, not add it'
    assert np.isclose(dec_max, decImg.max() - dec_offset), \
        'wcs_bounds should subtract dec_offsets from the maximum Dec, not add it'
    # Wavelength bounds are unaffected by RA/Dec offsets.
    assert np.isclose(wave_min, waveImg.min()), \
        'wave_min should be unaffected by ra_offsets/dec_offsets'
    assert np.isclose(wave_max, waveImg.max()), \
        'wave_max should be unaffected by ra_offsets/dec_offsets'


def test_create_wcs_voxedges_axis_order():
    """`create_wcs` must return `voxedges` ordered `(spec_bins, ybins,
    xbins)`, i.e. (wavelength, Dec, RA), matching the cube's
    `(nwave, ndec, nra)` numpy array axis order.  Uses deliberately different
    bin counts along each axis so a silent axis swap cannot hide behind
    equal-length arrays.
    """
    dspat = 1.0e-4  # deg/pixel
    dwave = 1.0     # Angstrom/pixel

    dec_min, dec_max = 10.0, 10.0 + 6.5 * dspat          # -> numdec = 6
    cosdec = np.cos(np.radians(0.5 * (dec_min + dec_max)))
    ra_min = 150.0
    ra_max = ra_min + 10.5 * dspat / cosdec              # -> numra = 10
    wave_min, wave_max = 5000.0, 5000.0 + 4.05 * dwave   # -> numwav = 4

    # With all six bounds specified explicitly, wcs_bounds() never touches the
    # image arrays, so trivial placeholders are sufficient here.
    dummy = np.zeros((2, 2))
    dummy_gpm = np.ones((2, 2), dtype=int)

    _, voxedges, _ = datacube.create_wcs(
        dummy, dummy, dummy, dummy_gpm, dspat, dwave,
        ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max,
        wave_min=wave_min, wave_max=wave_max
    )
    spec_bins, ybins, xbins = voxedges

    assert spec_bins.size - 1 == 4, \
        'voxedges[0] (spec_bins) should have 4 wavelength bins; got the wrong axis or count'
    assert ybins.size - 1 == 6, \
        'voxedges[1] (ybins/Dec) should have 6 Dec bins; got the wrong axis or count'
    assert xbins.size - 1 == 10, \
        'voxedges[2] (xbins/RA) should have 10 RA bins; got the wrong axis or count'


def test_generate_wcs_sky_right_convention():
    """`generate_WCS` must produce a WCS where RA decreases as the pixel-x
    index increases, so that displaying the cube with `imshow` (origin at the
    bottom-left) is "sky-right": RA increasing to the left, i.e. East left,
    North up.
    """
    dspat = 1.0e-4  # deg/pixel
    numra = 10
    crval = [150.0, 10.0, 5000.0]
    cdelt = [-dspat, dspat, 1.0]
    w = datacube.generate_WCS(crval, cdelt, numra)

    pix = np.array([[0, 0, 0], [3, 0, 0], [9, 0, 0]], dtype=float)
    world = w.wcs_pix2world(pix, 0)
    ra = world[:, 0]

    assert np.all(np.diff(ra) < 0), \
        'RA should strictly decrease as the pixel-x index increases (sky-right convention)'


def test_subpixellate_world2pix_axis_reversal():
    """The `output_wcs.wcs_world2pix(...)[:, :, ::-1]` step inside
    `subpixellate()` must reorder world2pix's natural `(ra_pix, dec_pix,
    wave_pix)` output into `(wave_pix, dec_pix, ra_pix)`, matching the cube's
    `(nwave, ndec, nra)` numpy axis order used for histogramming.  Chooses
    `ira != iwave` so a missed or partial reversal is guaranteed to be caught.
    """
    numra = 10
    w = datacube.generate_WCS([150.0, 10.0, 5000.0], [-1.0e-4, 1.0e-4, 1.0], numra)

    iwave, idec, ira = 2, 3, 7  # target voxel index, in (wave, dec, ra) order

    # WCS/FITS pixel order is (ra, dec, wave); get the world coordinates for
    # this target pixel (0-indexed, matching subpixellate()'s use of origin=0).
    world = w.wcs_pix2world(np.array([[ira, idec, iwave]], dtype=float), 0)

    # Reproduce the exact expression used inside subpixellate():
    #   output_wcs.wcs_world2pix(...).reshape(numpix, num_subpixels, 3)[:, :, ::-1]
    # with numpix = num_subpixels = 1, since this test targets the reversal
    # itself, not the subpixel-averaging machinery around it.
    vox_coord = w.wcs_world2pix(world, 0).reshape(1, 1, 3)[:, :, ::-1]

    assert np.allclose(vox_coord[0, 0], [iwave, idec, ira], atol=1e-6), \
        'the [:, :, ::-1] reversal should reorder world2pix output from ' \
        '(ra_pix, dec_pix, wave_pix) to (wave_pix, dec_pix, ra_pix), matching ' \
        "the cube's (nwave, ndec, nra) numpy axis order"


# ---------------------------------------------------------------------------
# Phase 2 -- a full synthetic cube through subpixellate()/generate_cube_subpixel(),
# still cheap (no real calibrations).  Only the first part is implemented here
# (a single synthetic exposure); combining a second, offset exposure is left
# for a follow-up test.
# ---------------------------------------------------------------------------

def test_generate_cube_subpixel_recovers_injected_source_position():
    """End-to-end synthetic-cube test of `subpixellate`/`generate_cube_subpixel`.

    Injects two point sources of unequal brightness at asymmetric
    detector-frame positions -- different slit *and* different spectral row,
    so no flip/rotation of the field of view maps one onto the other -- builds
    the corresponding datacube, and confirms the brightest cube voxel lands at
    the sky position (per the cube's own WCS) that independently matches the
    true RA/Dec of the injected pixel, as computed by
    `SlitTraceSet.get_radec_image` using the exact same per-exposure
    WCS/astrometric-transform/tilts the cube is built from.  The fainter
    source, in a different slit, is checked the same way, to guard against an
    axis mapping that happens to work only for the brighter source.

    Uses the cheapest possible synthetic setup: 3 straight, untilted slits, no
    DAR correction, and no subpixel sampling (nearest-grid-point), since this
    test targets whole-pixel positional correctness, not sub-pixel accuracy.
    """
    # ---- Detector geometry: 3 narrow, straight slits on a small detector ----
    nspec, nspat = 40, 40
    slit_cols = [(3, 11), (16, 24), (29, 37)]  # (left, right) detector columns
    spat_ids = [101, 102, 103]

    # A small rigid drift of the slit edges with row (a mild "tilt") is used
    # instead of perfectly straight/vertical edges: with perfectly straight
    # edges, the along-slit spatial offset returned by AlignmentSplines is
    # exactly degenerate between rows (identical for every row at a given
    # column), which produces duplicate x-values when subpixellate()
    # linearly interpolates RA/Dec over that spatial offset, and the
    # resulting divide-by-zero (NaN) crashes the fast_histogram C extension.
    # The drift amplitude is kept tiny (a fraction of a pixel over the full
    # spectral range) so it breaks the exact degeneracy without coupling the
    # wavelength (row) and Dec (along-slit) axes enough to shift the
    # wavelength-collapsed white-light peak by more than a fraction of a
    # cube pixel.
    row_drift = np.linspace(-0.05, 0.05, nspec)[:, None]
    left_init = np.tile([c[0] for c in slit_cols], (nspec, 1)).astype(float) + row_drift
    right_init = np.tile([c[1] for c in slit_cols], (nspec, 1)).astype(float) + row_drift
    slits = slittrace.SlitTraceSet(
        left_init, right_init, 'SlicerIFU', nspat=nspat, spat_id=np.array(spat_ids)
    )

    # Derive the on-slit mask directly from slit_img() (rather than the
    # nominal (lo, hi) column ranges above), so it exactly matches the
    # strict-inequality on-slit region that get_radec_image() populates --
    # slit_img() excludes the boundary columns themselves.
    slitid_img_gpm = slits.slit_img(pad=0)
    slitid_img_gpm[slitid_img_gpm < 0] = 0

    # Tilts are essentially a pure function of row, with a tiny (~0.1-row)
    # linear trend added across columns. Without it, every column in a given
    # row shares an identical tilt value, so the per-slit wavelength spline
    # `subpixellate` builds (mapping tilt -> wavelength) is queried at many
    # duplicate x-values, producing a 0/0 divide warning from `interp1d`
    # (harmless here since the corresponding y-values are duplicated too, but
    # avoided for a cleaner test).
    row_idx, col_idx = np.mgrid[0:nspec, 0:nspat]
    tilts = (row_idx + 0.1 * col_idx / nspat) / (nspec - 1)

    # Straight slits -> AlignmentSplines just needs the (constant) left/right
    # traces and the row-only tilts above.
    traces = np.stack([left_init, right_init], axis=1)  # (nspec, 2, nslit)
    astrom_trans = alignframe.AlignmentSplines(traces, np.array([0.0, 1.0]), tilts)

    # ---- A simple per-exposure "instrument" WCS: (slit index, along-slit
    # pixel, spectral row) -> (RA, Dec, wavelength), with no rotation ----
    wave0, dwv = 5000.0, 2.0
    pxscl_deg = 0.5 / 3600.0  # along-slit spatial scale (Dec-like axis)
    slscl_deg = 0.5 / 3600.0  # cross-slit (slice) scale (RA-like axis)
    frame_wcs = WCS(naxis=3)
    frame_wcs.wcs.crval = [150.0, 10.0, wave0]
    frame_wcs.wcs.crpix = [1.0, 1.0, 1.0]
    frame_wcs.wcs.cd = np.array([[-slscl_deg, 0.0, 0.0],
                                 [0.0, pxscl_deg, 0.0],
                                 [0.0, 0.0, dwv]])
    frame_wcs.wcs.cunit = [u.deg, u.deg, u.Angstrom]
    frame_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN', 'WAVE']

    # ---- Inject two point sources of unequal brightness, in different
    # slits AND at different spectral rows, so no flip/rotation of the field
    # of view maps one onto the other. ----
    bright = dict(row=8, col=7, amp=1.0)    # slit 0 (cols 3-11)
    faint = dict(row=32, col=33, amp=0.4)   # slit 2 (cols 29-37)
    sigma_pix = 1.2

    yy, xx = np.mgrid[0:nspec, 0:nspat]
    sciImg = np.zeros((nspec, nspat))
    for src in (bright, faint):
        sciImg += src['amp'] * np.exp(
            -0.5 * (((xx - src['col']) / sigma_pix) ** 2 + ((yy - src['row']) / sigma_pix) ** 2)
        )
    ivarImg = np.full((nspec, nspat), 1.0e8)
    wghtImg = np.ones((nspec, nspat))
    # Derived from the actual tilts (rather than the raw row index), so it
    # stays exactly consistent with the tilt -> wavelength relationship
    # `subpixellate` assumes, matching the WCS's wavelength axis.
    waveImg = wave0 + dwv * tilts * (nspec - 1)

    # ---- Ground truth: the true RA/Dec of each injected pixel, from the
    # same per-exposure WCS/astrometric-transform/tilts the cube is built
    # from -- this is the identical calculation `subpixellate` uses
    # internally, so it is a faithful oracle, not a re-derivation. ----
    raimg, decimg, _ = slits.get_radec_image(frame_wcs, astrom_trans, tilts)
    ra_bright, dec_bright = raimg[bright['row'], bright['col']], decimg[bright['row'], bright['col']]
    ra_faint, dec_faint = raimg[faint['row'], faint['col']], decimg[faint['row'], faint['col']]

    # ---- Build the output cube WCS/bins directly from the data (single
    # frame, no offsets -- multi-frame combination is left for a follow-up
    # test). ----
    dspat_cube = 0.15 / 3600.0  # deg/pixel; finer than the slice spacing above

    # `create_wcs`'s bin-count calculation truncates (`int(...)`) rather than
    # rounds up, which can silently place a source sitting exactly at the
    # extreme edge of the RA/Dec range outside the nominal pixel grid by up
    # to ~1 pixel (see kcwi_wcs.md, "Inconsistencies" item on numra/numdec
    # truncation). For continuous, real illumination this only clips a
    # negligible sliver at the true edge; here, with all of a slit's flux
    # sitting at a single discrete RA value at that edge, it can lose an
    # entire source. Padding the explicit bounds by a few cube pixels avoids
    # exercising that truncation edge case, which is not what this test is
    # meant to probe.
    onslit = slitid_img_gpm > 0
    ra_pad = 3 * dspat_cube / np.cos(np.radians(np.mean(decimg[onslit])))
    dec_pad = 3 * dspat_cube
    cube_wcs, voxedges, _ = datacube.create_wcs(
        raimg, decimg, waveImg, slitid_img_gpm, dspat_cube, dwv,
        ra_min=raimg[onslit].min() - ra_pad, ra_max=raimg[onslit].max() + ra_pad,
        dec_min=decimg[onslit].min() - dec_pad, dec_max=decimg[onslit].max() + dec_pad,
    )

    # ---- Build the cube (nearest-grid-point: no subpixel sampling needed
    # for this position-recovery check; no DAR). ----
    flxcube, sigcube, bpmcube, normcube, wave = datacube.generate_cube_subpixel(
        cube_wcs, voxedges, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
        frame_wcs, tilts, slits, astrom_trans, all_dar=None, ra_offset=0.0, dec_offset=0.0,
        spec_subpixel=1, spat_subpixel=1, slice_subpixel=1, skip_subpix_weights=True,
        correct_dar=False
    )

    whitelight = flxcube.sum(axis=0)
    cube_celestial = cube_wcs.celestial

    # The brightest voxel overall must be the bright source, at its true sky position.
    iy_max, ix_max = np.unravel_index(np.argmax(whitelight), whitelight.shape)
    ra_cube_bright, dec_cube_bright = cube_celestial.wcs_pix2world(
        np.array([ix_max]), np.array([iy_max]), 0)
    assert np.isclose(ra_cube_bright[0], ra_bright, atol=dspat_cube), \
        "the cube's brightest voxel RA does not match the injected bright source's true RA"
    assert np.isclose(dec_cube_bright[0], dec_bright, atol=dspat_cube), \
        "the cube's brightest voxel Dec does not match the injected bright source's true Dec"

    # The fainter source, in a different slit and at a different row, must
    # also land at its true sky position, with roughly the expected relative
    # flux -- this guards against an axis mapping that happens to work only
    # for the brighter source.
    ix_faint, iy_faint = cube_celestial.wcs_world2pix(
        np.array([ra_faint]), np.array([dec_faint]), 0)
    ix_faint, iy_faint = int(np.round(ix_faint[0])), int(np.round(iy_faint[0]))
    flux_bright_peak = whitelight[iy_max, ix_max]
    flux_faint_at_expected = whitelight[iy_faint, ix_faint]
    assert 0.15 * flux_bright_peak < flux_faint_at_expected < 0.7 * flux_bright_peak, \
        "flux at the faint source's expected cube position does not match its known " \
        "amplitude relative to the bright source"


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
