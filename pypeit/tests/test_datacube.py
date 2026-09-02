"""
Module to run tests on datacube generation
"""
from pathlib import Path

from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord
from IPython import embed
import numpy as np
import pytest

from pypeit import utils
from pypeit.core import datacube
from pypeit.core.flexure import calculate_image_phase

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
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value*cosdec, atol=0.001)
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001)

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
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol)
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol)


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
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value * cosdec, atol=0.001)
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001)

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
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol)
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol)


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
    assert np.isclose(dra.to(u.arcsec).value, offs_ra.to(u.arcsec).value * cosdec, atol=0.001)
    assert np.isclose(ddec.to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=0.001)

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
            dec_pix_star[ff], ra_pix_star[ff] = gaussian_position
        ra_shifts = -(ra_pix_star[ref_idx] - ra_pix_star) * _dspat / cosdec
        dec_shifts = (dec_pix_star[ref_idx] - dec_pix_star) * _dspat
        ra_offsets = [ra_offsets[ff] - ra_shifts[ff] for ff in range(numfiles)]
        dec_offsets = [dec_offsets[ff] - dec_shifts[ff] for ff in range(numfiles)]

    # Check that the offsets agree to within 1/3 of a pixel
    atol = _dspat.to(u.arcsec).value/3
    assert np.isclose(ra_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(dec_offsets[0].to(u.arcsec).value, 0.0, atol=atol)
    assert np.isclose(ra_offsets[1].to(u.arcsec).value, offs_ra.to(u.arcsec).value, atol=atol)
    assert np.isclose(dec_offsets[1].to(u.arcsec).value, offs_dec.to(u.arcsec).value, atol=atol)


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
