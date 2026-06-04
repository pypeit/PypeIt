"""
Module to run tests on datacube generation
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.wcs import WCS
import astropy.units as u
from astropy.coordinates import SkyCoord

from pypeit.core import datacube
from pypeit.core.flexure import calculate_image_phase


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

def test_align():
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
