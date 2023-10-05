"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit.core import coadd, flux_calib

from IPython import embed


def gaussian2D(tup, intflux, xo, yo, sigma_x, sigma_y, theta, offset):
    """
    Fit a 2D Gaussian function to an image.

    Args:
        tup (:obj:`tuple`):
            A two element tuple containing the x and y coordinates of each pixel
            in the image
        intflux (float):
            The Integrated flux of the 2D Gaussian
        xo (float):
            The centre of the Gaussian along the x-coordinate when z=0
        yo (float):
            The centre of the Gaussian along the y-coordinate when z=0
        sigma_x (float):
            The standard deviation in the x-direction
        sigma_y (float):
            The standard deviation in the y-direction
        theta (float):
            The orientation angle of the 2D Gaussian
        offset (float):
            Constant offset

    Returns:
        `numpy.ndarray`_: The 2D Gaussian evaluated at the coordinate (x, y)
    """
    # Extract the (x, y, z) coordinates of each pixel from the tuple
    (x, y) = tup
    # Ensure these are floating point
    xo = float(xo)
    yo = float(yo)
    # Account for a rotated 2D Gaussian
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    # Normalise so that the integrated flux is a parameter, instead of the amplitude
    norm = 1/(2*np.pi*np.sqrt(a*c-b*b))
    gtwod = offset + norm*intflux*np.exp(-(a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return gtwod.ravel()


def fitGaussian2D(image, norm=False):
    """
    Fit a 2D Gaussian to an input image. It is recommended that the input image
    is scaled to a maximum value that is ~1, so that all fit parameters are of
    the same order of magnitude. Set norm=True if you do not care about the
    amplitude or integrated flux. Otherwise, make sure you scale the image by
    a known value prior to passing it into this function.

    Parameters
    ----------
    image : `numpy.ndarray`_
        A 2D input image
    norm : bool, optional
        If True, the input image will be normalised to the maximum value
        of the input image.

    Returns
    -------
    popt : `numpy.ndarray`_
       The optimum parameters of the Gaussian in the following order: Integrated
       flux, x center, y center, sigma_x, sigma_y, theta, offset. See
       :func:`~pypeit.core.datacube.gaussian2D` for a more detailed description
       of the model.
    pcov : `numpy.ndarray`_
        Corresponding covariance matrix
    """
    # Normalise if requested
    wlscl = np.max(image) if norm else 1
    # Setup the coordinates
    x = np.linspace(0, image.shape[0] - 1, image.shape[0])
    y = np.linspace(0, image.shape[1] - 1, image.shape[1])
    xx, yy = np.meshgrid(x, y, indexing='ij')
    # Setup the fitting params
    idx_max = [image.shape[0]/2, image.shape[1]/2]  # Just use the centre of the image as the best guess
    #idx_max = np.unravel_index(np.argmax(image), image.shape)
    initial_guess = (1, idx_max[0], idx_max[1], 2, 2, 0, 0)
    bounds = ([0, 0, 0, 0.5, 0.5, -np.pi, -np.inf],
              [np.inf, image.shape[0], image.shape[1], image.shape[0], image.shape[1], np.pi, np.inf])
    # Perform the fit
    popt, pcov = opt.curve_fit(gaussian2D, (xx, yy), image.ravel() / wlscl, bounds=bounds, p0=initial_guess)
    # Return the fitting results
    return popt, pcov


def dar_fitfunc(radec, coord_ra, coord_dec, datfit, wave, obstime, location, pressure,
                temperature, rel_humidity):
    """
    Generates a fitting function to calculate the offset due to differential
    atmospheric refraction

    Args:
        radec (tuple):
            A tuple containing two floats representing the shift in ra and dec
            due to DAR.
        coord_ra (float):
            RA in degrees
        coord_dec (float):
            Dec in degrees
        datfit (`numpy.ndarray`_):
            The RA and DEC that the model needs to match
        wave (float):
            Wavelength to calculate the DAR
        location (`astropy.coordinates.EarthLocation`_):
            observatory location
        pressure (float):
            Outside pressure at `location`
        temperature (float):
            Outside ambient air temperature at `location`
        rel_humidity (float):
            Outside relative humidity at `location`. This should be between 0 to 1.

    Returns:
        float: chi-squared difference between datfit and model
    """
    (diff_ra, diff_dec) = radec
    # Generate the coordinate with atmospheric conditions
    coord_atmo = SkyCoord(coord_ra + diff_ra, coord_dec + diff_dec, unit=(units.deg, units.deg))
    coord_altaz = coord_atmo.transform_to(AltAz(obstime=obstime, location=location, obswl=wave,
                                          pressure=pressure, temperature=temperature,
                                          relative_humidity=rel_humidity))
    # Return chi-squared value
    return np.sum((np.array([coord_altaz.alt.value, coord_altaz.az.value])-datfit)**2)


def correct_grating_shift(wave_eval, wave_curr, spl_curr, wave_ref, spl_ref, order=2):
    """
    Using spline representations of the blaze profile, calculate the grating
    correction that should be applied to the current spectrum (suffix ``curr``)
    relative to the reference spectrum (suffix ``ref``). The grating correction
    is then evaluated at the wavelength array given by ``wave_eval``.

    Args:
        wave_eval (`numpy.ndarray`_):
            Wavelength array to evaluate the grating correction
        wave_curr (`numpy.ndarray`_):
            Wavelength array used to construct spl_curr
        spl_curr (`scipy.interpolate.interp1d`_):
            Spline representation of the current blaze function (based on the illumflat).
        wave_ref (`numpy.ndarray`_):
            Wavelength array used to construct spl_ref
        spl_ref (`scipy.interpolate.interp1d`_):
            Spline representation of the reference blaze function (based on the illumflat).
        order (int):
            Polynomial order used to fit the grating correction.

    Returns:
        `numpy.ndarray`_: The grating correction to apply
    """
    msgs.info("Calculating the grating correction")
    # Calculate the grating correction
    grat_corr_tmp = spl_curr(wave_eval) / spl_ref(wave_eval)
    # Determine the useful overlapping wavelength range
    minw, maxw = max(np.min(wave_curr), np.min(wave_ref)), max(np.min(wave_curr), np.max(wave_ref))
    # Perform a low-order polynomial fit to the grating correction (should be close to linear)
    wave_corr = (wave_eval - minw) / (maxw - minw)  # Scale wavelengths to be of order 0-1
    wblz = np.where((wave_corr > 0.1) & (wave_corr < 0.9))  # Remove the pixels that are within 10% of the edges
    coeff_gratcorr = np.polyfit(wave_corr[wblz], grat_corr_tmp[wblz], order)
    grat_corr = np.polyval(coeff_gratcorr, wave_corr)
    # Return the estimates grating correction
    return grat_corr


def extract_standard_spec(stdcube, subpixel=20):
    """
    Extract a spectrum of a standard star from a datacube

    Parameters
    ----------
    std_cube : `astropy.io.fits.HDUList`_
        An HDU list of fits files
    subpixel : int
        Number of pixels to subpixelate spectrum when creating mask

    Returns
    -------
    wave : `numpy.ndarray`_
        Wavelength of the star.
    Nlam_star : `numpy.ndarray`_
        counts/second/Angstrom
    Nlam_ivar_star : `numpy.ndarray`_
        inverse variance of Nlam_star
    gpm_star : `numpy.ndarray`_
        good pixel mask for Nlam_star
    """
    # Extract some information from the HDU list
    flxcube = stdcube['FLUX'].data.T.copy()
    varcube = stdcube['SIG'].data.T.copy()**2
    bpmcube = stdcube['BPM'].data.T.copy()
    numwave = flxcube.shape[2]

    # Setup the WCS
    stdwcs = wcs.WCS(stdcube['FLUX'].header)

    wcs_scale = (1.0 * stdwcs.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
    wave = wcs_scale * stdwcs.spectral.wcs_pix2world(np.arange(numwave), 0)[0]

    # Generate a whitelight image, and fit a 2D Gaussian to estimate centroid and width
    wl_img = make_whitelight_fromcube(flxcube)
    popt, pcov = fitGaussian2D(wl_img, norm=True)
    wid = max(popt[3], popt[4])

    # Setup the coordinates of the mask
    x = np.linspace(0, flxcube.shape[0] - 1, flxcube.shape[0] * subpixel)
    y = np.linspace(0, flxcube.shape[1] - 1, flxcube.shape[1] * subpixel)
    xx, yy = np.meshgrid(x, y, indexing='ij')

    # Generate a mask
    newshape = (flxcube.shape[0] * subpixel, flxcube.shape[1] * subpixel)
    mask = np.zeros(newshape)
    nsig = 4  # 4 sigma should be far enough... Note: percentage enclosed for 2D Gaussian = 1-np.exp(-0.5 * nsig**2)
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < nsig * wid))
    mask[ww] = 1
    mask = utils.rebinND(mask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)

    # Generate a sky mask
    newshape = (flxcube.shape[0] * subpixel, flxcube.shape[1] * subpixel)
    smask = np.zeros(newshape)
    nsig = 8  # 8 sigma should be far enough
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < nsig * wid))
    smask[ww] = 1
    smask = utils.rebinND(smask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)
    smask -= mask

    # Subtract the residual sky
    skymask = np.logical_not(bpmcube) * smask
    skycube = flxcube * skymask
    skyspec = skycube.sum(0).sum(0)
    nrmsky = skymask.sum(0).sum(0)
    skyspec *= utils.inverse(nrmsky)
    flxcube -= skyspec.reshape((1, 1, numwave))

    # Subtract the residual sky from the whitelight image
    sky_val = np.sum(wl_img[:, :, np.newaxis] * smask) / np.sum(smask)
    wl_img -= sky_val

    msgs.info("Extracting a boxcar spectrum of datacube")
    # Construct an image that contains the fraction of flux included in the
    # boxcar extraction at each wavelength interval
    norm_flux = wl_img[:,:,np.newaxis] * mask
    norm_flux /= np.sum(norm_flux)
    # Extract boxcar
    cntmask = np.logical_not(bpmcube) * mask  # Good pixels within the masked region around the standard star
    flxscl = (norm_flux * cntmask).sum(0).sum(0)  # This accounts for the flux that is missing due to masked pixels
    scimask = flxcube * cntmask
    varmask = varcube * cntmask**2
    nrmcnt = utils.inverse(flxscl)
    box_flux = scimask.sum(0).sum(0) * nrmcnt
    box_var = varmask.sum(0).sum(0) * nrmcnt**2
    box_gpm = flxscl > 1/3  # Good pixels are those where at least one-third of the standard star flux is measured
    # Setup the return values
    ret_flux, ret_var, ret_gpm = box_flux, box_var, box_gpm

    # Convert from counts/s/Ang/arcsec**2 to counts/s/Ang
    arcsecSQ = 3600.0*3600.0*(stdwcs.wcs.cdelt[0]*stdwcs.wcs.cdelt[1])
    ret_flux *= arcsecSQ
    ret_var *= arcsecSQ**2
    # Return the box extraction results
    return wave, ret_flux, utils.inverse(ret_var), ret_gpm


def make_sensfunc(ss_file, senspar, blaze_wave=None, blaze_spline=None, grating_corr=False):
    """
    Generate the sensitivity function from a standard star DataCube.

    Args:
        ss_file (:obj:`str`_):
            The relative path and filename of the standard star datacube. It should be fits format, and
            for full functionality, should ideally of the form `pypeit.coadd3d.DataCube`_
        senspar (:class:`~pypeit.par.pypeitpar.SensFuncPar`):
            The parameters required for the sensitivity function computation.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array used to construct blaze_spline
        blaze_spline (`scipy.interpolate.interp1d`_, optional):
            Spline representation of the reference blaze function (based on the illumflat).
        grating_corr (:obj:`bool`_, optional):
            If a grating correction should be performed, set this variable to True.

    Returns:
        `numpy.ndarray`_: A mask of the good sky pixels (True = good)
    """
    # Check if the standard star datacube exists
    if not os.path.exists(ss_file):
        msgs.error("Standard cube does not exist:" + msgs.newline() + ss_file)
    msgs.info(f"Loading standard star cube: {ss_file:s}")
    # Load the standard star cube and retrieve its RA + DEC
    stdcube = fits.open(ss_file)
    star_ra, star_dec = stdcube[1].header['CRVAL1'], stdcube[1].header['CRVAL2']

    # Extract a spectrum of the standard star
    wave, Nlam_star, Nlam_ivar_star, gpm_star = extract_standard_spec(stdcube)

    # Extract the information about the blaze
    if grating_corr:
        blaze_wave_curr, blaze_spec_curr = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
        blaze_spline_curr = interp1d(blaze_wave_curr, blaze_spec_curr,
                                     kind='linear', bounds_error=False, fill_value="extrapolate")
        # Perform a grating correction
        grat_corr = correct_grating_shift(wave, blaze_wave_curr, blaze_spline_curr, blaze_wave, blaze_spline)
        # Apply the grating correction to the standard star spectrum
        Nlam_star /= grat_corr
        Nlam_ivar_star *= grat_corr ** 2

    # Read in some information above the standard star
    std_dict = flux_calib.get_standard_spectrum(star_type=senspar['star_type'],
                                                star_mag=senspar['star_mag'],
                                                ra=star_ra, dec=star_dec)
    # Calculate the sensitivity curve
    # TODO :: This needs to be addressed... unify flux calibration into the main PypeIt routines.
    msgs.warn("Datacubes are currently flux-calibrated using the UVIS algorithm... this will be deprecated soon")
    zeropoint_data, zeropoint_data_gpm, zeropoint_fit, zeropoint_fit_gpm = \
        flux_calib.fit_zeropoint(wave, Nlam_star, Nlam_ivar_star, gpm_star, std_dict,
                                 mask_hydrogen_lines=senspar['mask_hydrogen_lines'],
                                 mask_helium_lines=senspar['mask_helium_lines'],
                                 hydrogen_mask_wid=senspar['hydrogen_mask_wid'],
                                 nresln=senspar['UVIS']['nresln'],
                                 resolution=senspar['UVIS']['resolution'],
                                 trans_thresh=senspar['UVIS']['trans_thresh'],
                                 polyorder=senspar['polyorder'],
                                 polycorrect=senspar['UVIS']['polycorrect'],
                                 polyfunc=senspar['UVIS']['polyfunc'])
    wgd = np.where(zeropoint_fit_gpm)
    sens = np.power(10.0, -0.4 * (zeropoint_fit[wgd] - flux_calib.ZP_UNIT_CONST)) / np.square(wave[wgd])
    return interp1d(wave[wgd], sens, kind='linear', bounds_error=False, fill_value="extrapolate")


def make_good_skymask(slitimg, tilts):
    """
    Mask the spectral edges of each slit (i.e. the pixels near the ends of the
    detector in the spectral direction). Some extreme values of the tilts are
    only sampled with a small fraction of the pixels of the slit width. This
    leads to a bad extrapolation/determination of the sky model.

    Args:
        slitimg (`numpy.ndarray`_):
            An image of the slit indicating which slit each pixel belongs to
        tilts (`numpy.ndarray`_):
            Spectral tilts.

    Returns:
        `numpy.ndarray`_: A mask of the good sky pixels (True = good)
    """
    msgs.info("Masking edge pixels where the sky model is poor")
    # Initialise the GPM
    gpm = np.zeros(slitimg.shape, dtype=bool)
    # Find unique slits
    unq = np.unique(slitimg[slitimg>0])
    for uu in range(unq.size):
        # Find the x,y pixels in this slit
        ww = np.where((slitimg == unq[uu]) & (tilts != 0.0))
        # Mask the bottom pixels first
        wb = np.where(ww[0] == 0)[0]
        wt = np.where(ww[0] == np.max(ww[0]))[0]
        # Calculate the maximum tilt from the bottom row, and the miminum tilt from the top row
        maxtlt = np.max(tilts[0,  ww[1][wb]])
        mintlt = np.min(tilts[-1, ww[1][wt]])
        # Mask all values below this maximum
        gpm[ww] = (tilts[ww] >= maxtlt) & (tilts[ww] <= mintlt)  # The signs are correct here.
    return gpm


def get_output_filename(fil, par_outfile, combine, idx=1):
    """
    Get the output filename of a datacube, given the input

    Args:
        fil (str):
            The spec2d filename.
        par_outfile (str):
            The user-specified output filename (see cubepar['output_filename'])
        combine (bool):
            Should the input frames be combined into a single datacube?
        idx (int, optional):
            Index of filename to be saved. Required if combine=False.

    Returns:
        str: The output filename to use.
    """
    if combine:
        if par_outfile == '':
            par_outfile = 'datacube.fits'
        # Check if we needs to append an extension
        return par_outfile if '.fits' in par_outfile else f'{par_outfile}.fits'
    if par_outfile == '':
        return fil.replace('spec2d_', 'spec3d_')
    # Finally, if nothing else, use the output filename as a prefix, and a numerical suffic
    return os.path.splitext(par_outfile)[0] + f'_{idx:03}.fits'


def get_output_whitelight_filename(outfile):
    """
    Given the output filename of a datacube, create an appropriate whitelight
    fits file name

    Args:
        outfile (str):
            The output filename used for the datacube.

    Returns:
        A string containing the output filename to use for the whitelight image.
    """
    return os.path.splitext(outfile)[0] + "_whitelight.fits"


def get_whitelight_pixels(all_wave, min_wl, max_wl):
    """
    Determine which pixels are included within the specified wavelength range

    Args:
        all_wave (`numpy.ndarray`_):
            The wavelength of each individual pixel
        min_wl (float):
            Minimum wavelength to consider
        max_wl (float):
            Maximum wavelength to consider

    Returns:
        :obj:`tuple`: A `numpy.ndarray`_ object with the indices of all_wave
        that contain pixels within the requested wavelength range, and a float
        with the wavelength range (i.e. maximum wavelength - minimum wavelength)
    """
    wavediff = np.max(all_wave) - np.min(all_wave)
    if min_wl < max_wl:
        ww = np.where((all_wave > min_wl) & (all_wave < max_wl))
        wavediff = max_wl - min_wl
    else:
        msgs.warn("Datacubes do not completely overlap in wavelength. Offsets may be unreliable...")
        ww = (np.arange(all_wave.size),)
    return ww, wavediff


def get_whitelight_range(wavemin, wavemax, wl_range):
    """
    Get the wavelength range to use for the white light images

    Parameters
    ----------
    wavemin : float
        Automatically determined minimum wavelength to use for making the white
        light image.
    wavemax : float
        Automatically determined maximum wavelength to use for making the white
        light image.
    wl_range : list
        Two element list containing the user-specified values to manually
        override the automated values determined by PypeIt.

    Returns
    -------
    wlrng : list
        A two element list containing the minimum and maximum wavelength to use
        for the white light images
    """
    wlrng = [wavemin, wavemax]
    if wl_range[0] is not None:
        if wl_range[0] < wavemin:
            msgs.warn("The user-specified minimum wavelength ({0:.2f}) to use for the white light".format(wl_range[0]) +
                      msgs.newline() + "images is lower than the recommended value ({0:.2f}),".format(wavemin) +
                      msgs.newline() + "which ensures that all spaxels cover the same wavelength range.")
        wlrng[0] = wl_range[0]
    if wl_range[1] is not None:
        if wl_range[1] > wavemax:
            msgs.warn("The user-specified maximum wavelength ({0:.2f}) to use for the white light".format(wl_range[1]) +
                      msgs.newline() + "images is greater than the recommended value ({0:.2f}),".format(wavemax) +
                      msgs.newline() + "which ensures that all spaxels cover the same wavelength range.")
        wlrng[1] = wl_range[1]
    msgs.info("The white light images will cover the wavelength range: {0:.2f}A - {1:.2f}A".format(wlrng[0], wlrng[1]))
    return wlrng


def make_whitelight_fromcube(cube, wave=None, wavemin=None, wavemax=None):
    """
    Generate a white light image using an input cube.

    Args:
        cube (`numpy.ndarray`_):
            3D datacube (the final element contains the wavelength dimension)
        wave (`numpy.ndarray`_, optional):
            1D wavelength array. Only required if wavemin or wavemax are not
            None.
        wavemin (float, optional):
            Minimum wavelength (same units as wave) to be included in the
            whitelight image.  You must provide wave as well if you want to
            reduce the wavelength range.
        wavemax (float, optional):
            Maximum wavelength (same units as wave) to be included in the
            whitelight image.  You must provide wave as well if you want to
            reduce the wavelength range.

    Returns:
        A whitelight image of the input cube (of type `numpy.ndarray`_).
    """
    # Make a wavelength cut, if requested
    cutcube = cube.copy()
    if wavemin is not None or wavemax is not None:
        # Make some checks on the input
        if wave is None:
            msgs.error("wave variable must be supplied to create white light image with wavelength cuts")
        else:
            if wave.size != cube.shape[2]:
                msgs.error("wave variable should have the same length as the third axis of cube.")
        # assign wavemin & wavemax if one is not provided
        if wavemin is None:
            wavemin = np.min(wave)
        if wavemax is None:
            wavemax = np.max(wave)
        ww = np.where((wave >= wavemin) & (wave <= wavemax))[0]
        wmin, wmax = ww[0], ww[-1]+1
        cutcube = cube[:, :, wmin:wmax]
    # Now sum along the wavelength axis
    nrmval = np.sum(cutcube != 0.0, axis=2)
    nrmval[nrmval == 0.0] = 1.0
    wl_img = np.sum(cutcube, axis=2) / nrmval
    return wl_img


def load_imageWCS(filename, ext=0):
    """
    Load an image and return the image and the associated WCS.

    Args:
        filename (str):
            A fits filename of an image to be used when generating white light
            images. Note, the fits file must have a valid 3D WCS.
        ext (bool, optional):
            The extension that contains the image and WCS

    Returns:
        :obj:`tuple`: An `numpy.ndarray`_ with the 2D image data and a
        `astropy.wcs.WCS`_ with the image WCS.
    """
    imghdu = fits.open(filename)
    image = imghdu[ext].data.T
    imgwcs = wcs.WCS(imghdu[ext].header)
    # Return required info
    return image, imgwcs


def align_user_offsets(all_ra, all_dec, all_idx, ifu_ra, ifu_dec, ra_offset, dec_offset):
    """
    Align the RA and DEC of all input frames, and then
    manually shift the cubes based on user-provided offsets.
    The offsets should be specified in arcseconds, and the
    ra_offset should include the cos(dec) factor.

    Args:
        all_ra (`numpy.ndarray`_):
            A 1D array containing the RA values of each detector pixel of every frame.
        all_dec (`numpy.ndarray`_):
            A 1D array containing the Dec values of each detector pixel of every frame.
            Same size as all_ra.
        all_idx (`numpy.ndarray`_):
            A 1D array containing an ID value for each detector frame (0-indexed).
            Same size as all_ra.
        ifu_ra (`numpy.ndarray`_):
            A list of RA values of the IFU (one value per frame)
        ifu_dec (`numpy.ndarray`_):
            A list of Dec values of the IFU (one value per frame)
        ra_offset (`numpy.ndarray`_):
            A list of RA offsets to be applied to the input pixel values (one value per frame).
            Note, the ra_offset MUST contain the cos(dec) factor. This is the number of arcseconds
            on the sky that represents the telescope offset.
        dec_offset (`numpy.ndarray`_):
            A list of Dec offsets to be applied to the input pixel values (one value per frame).

    Returns:
        A tuple containing a new set of RA and Dec values that have been aligned. Both arrays
        are of type `numpy.ndarray`_.
    """
    # First, translate all coordinates to the coordinates of the first frame
    # Note: You do not need cos(dec) here, this just overrides the IFU coordinate centre of each frame
    #       The cos(dec) factor should be input by the user, and should be included in the self.opts['ra_offset']
    ref_shift_ra = ifu_ra[0] - ifu_ra
    ref_shift_dec = ifu_dec[0] - ifu_dec
    numfiles = ra_offset.size
    for ff in range(numfiles):
        # Apply the shift
        all_ra[all_idx == ff] += ref_shift_ra[ff] + ra_offset[ff] / 3600.0
        all_dec[all_idx == ff] += ref_shift_dec[ff] + dec_offset[ff] / 3600.0
        msgs.info("Spatial shift of cube #{0:d}:" + msgs.newline() +
                  "RA, DEC (arcsec) = {1:+0.3f} E, {2:+0.3f} N".format(ff + 1, ra_offset[ff], dec_offset[ff]))
    return all_ra, all_dec


def create_wcs(all_ra, all_dec, all_wave, dspat, dwave,
               ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
               reference=None, collapse=False, equinox=2000.0, specname="PYP_SPEC"):
    """
    Create a WCS and the expected edges of the voxels, based on user-specified
    parameters or the extremities of the data.

    Parameters
    ----------
    all_ra : `numpy.ndarray`_
        1D flattened array containing the RA values of each pixel from all
        spec2d files
    all_dec : `numpy.ndarray`_
        1D flattened array containing the DEC values of each pixel from all
        spec2d files
    all_wave : `numpy.ndarray`_
        1D flattened array containing the wavelength values of each pixel from
        all spec2d files
    dspat : float
        Spatial size of each square voxel (in arcsec). The default is to use the
        values in cubepar.
    dwave : float
        Linear wavelength step of each voxel (in Angstroms)
    ra_min : float, optional
        Minimum RA of the WCS (degrees)
    ra_max : float, optional
        Maximum RA of the WCS (degrees)
    dec_min : float, optional
        Minimum Dec of the WCS (degrees)
    dec_max : float, optional
        Maximum Dec of the WCS (degrees)
    wave_min : float, optional
        Minimum wavelength of the WCS (degrees)
    wave_max : float, optional
        Maximum wavelength of the WCS (degrees)
    reference : str, optional
        Filename of a fits file that contains a WCS in the Primary HDU.
    collapse : bool, optional
        If True, the spectral dimension will be collapsed to a single channel
        (primarily for white light images)
    equinox : float, optional
        Equinox of the WCS
    specname : str, optional
        Name of the spectrograph

    Returns
    -------
    cubewcs : `astropy.wcs.WCS`_
        astropy WCS to be used for the combined cube
    voxedges : tuple
        A three element tuple containing the bin edges in the x, y (spatial) and
        z (wavelength) dimensions
    reference_image : `numpy.ndarray`_
        The reference image to be used for the cross-correlation. Can be None.
    """
    # Grab cos(dec) for convenience
    cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)

    # Setup the cube ranges
    _ra_min = ra_min if ra_min is not None else np.min(all_ra)
    _ra_max = ra_max if ra_max is not None else np.max(all_ra)
    _dec_min = dec_min if dec_min is not None else np.min(all_dec)
    _dec_max = dec_max if dec_max is not None else np.max(all_dec)
    _wav_min = wave_min if wave_min is not None else np.min(all_wave)
    _wav_max = wave_max if wave_max is not None else np.max(all_wave)
    # dwave = self.cubepar['wave_delta'] if self.cubepar['wave_delta'] is not None else dwv

    # Number of voxels in each dimension
    numra = int((_ra_max - _ra_min) * cosdec / dspat)
    numdec = int((_dec_max - _dec_min) / dspat)
    numwav = int(np.round((_wav_max - _wav_min) / dwave))

    # If a white light WCS is being generated, make sure there's only 1 wavelength bin
    if collapse:
        _wav_min = np.min(all_wave)
        _wav_max = np.max(all_wave)
        dwave = _wav_max - _wav_min
        numwav = 1

    # Generate a master WCS to register all frames
    coord_min = [_ra_min, _dec_min, _wav_min]
    coord_dlt = [dspat, dspat, dwave]

    # If a reference image is being used and a white light image is requested (collapse=True) update the celestial parts
    reference_image = None
    if reference is not None:
        # Load the requested reference image
        reference_image, imgwcs = load_imageWCS(reference)
        # Update the celestial WCS
        coord_min[:2] = imgwcs.wcs.crval
        coord_dlt[:2] = imgwcs.wcs.cdelt
        numra, numdec = reference_image.shape

    cubewcs = generate_WCS(coord_min, coord_dlt, equinox=equinox, name=specname)
    msgs.info(msgs.newline() + "-" * 40 +
              msgs.newline() + "Parameters of the WCS:" +
              msgs.newline() + "RA   min = {0:f}".format(coord_min[0]) +
              msgs.newline() + "DEC  min = {0:f}".format(coord_min[1]) +
              msgs.newline() + "WAVE min, max = {0:f}, {1:f}".format(_wav_min, _wav_max) +
              msgs.newline() + "Spaxel size = {0:f} arcsec".format(3600.0 * dspat) +
              msgs.newline() + "Wavelength step = {0:f} A".format(dwave) +
              msgs.newline() + "-" * 40)

    # Generate the output binning
    xbins = np.arange(1 + numra) - 0.5
    ybins = np.arange(1 + numdec) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    voxedges = (xbins, ybins, spec_bins)
    return cubewcs, voxedges, reference_image


def generate_WCS(crval, cdelt, equinox=2000.0, name="PYP_SPEC"):
    """
    Generate a WCS that will cover all input spec2D files

    Args:
        crval (list):
            3 element list containing the [RA, DEC, WAVELENGTH] of
            the reference pixel
        cdelt (list):
            3 element list containing the delta values of the [RA,
            DEC, WAVELENGTH]
        equinox (float, optional):
            Equinox of the WCS

    Returns:
        `astropy.wcs.WCS`_ : astropy WCS to be used for the combined cube
    """
    # Create a new WCS object.
    msgs.info("Generating WCS")
    w = wcs.WCS(naxis=3)
    w.wcs.equinox = equinox
    w.wcs.name = name
    w.wcs.radesys = 'FK5'
    # Insert the coordinate frame
    w.wcs.cname = ['RA', 'DEC', 'Wavelength']
    w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
    w.wcs.crval = crval  # RA, DEC, and wavelength zeropoints
    w.wcs.crpix = [0, 0, 0]  # RA, DEC, and wavelength reference pixels
    #w.wcs.cd = np.array([[cdval[0], 0.0, 0.0], [0.0, cdval[1], 0.0], [0.0, 0.0, cdval[2]]])
    w.wcs.cdelt = cdelt
    w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
    w.wcs.latpole = 0.0  # Native latitude of the Celestial pole
    return w


def compute_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_idx, whitelight_img, dspat, dwv,
                    sn_smooth_npix=None, relative_weights=False):
    r"""
    Calculate wavelength dependent optimal weights. The weighting is currently
    based on a relative :math:`(S/N)^2` at each wavelength

    Args:
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the RA values of each pixel from all
            spec2d files
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the DEC values of each pixel from all
            spec2d files
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength values of each pixel
            from all spec2d files
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all
            spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel
            from all spec2d files
        all_idx (`numpy.ndarray`_):
            1D flattened array containing an integer identifier indicating which
            spec2d file each pixel originates from. For example, a 0 would
            indicate that a pixel originates from the first spec2d frame listed
            in the input file. a 1 would indicate that this pixel originates
            from the second spec2d file, and so forth.
        whitelight_img (`numpy.ndarray`_):
            A 2D array containing a whitelight image, that was created with the
            input ``all_`` arrays.
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        dwv (float):
            The size of each wavelength pixel (in Angstroms)
        sn_smooth_npix (float, optional):
            Number of pixels used for determining smoothly varying S/N ratio
            weights.  This is currently not required, since a relative weighting
            scheme with a polynomial fit is used to calculate the S/N weights.
        relative_weights (bool, optional):
            Calculate weights by fitting to the ratio of spectra?

    Returns:
        `numpy.ndarray`_ : a 1D array the same size as all_sci, containing
        relative wavelength dependent weights of each input pixel.
    """
    msgs.info("Calculating the optimal weights of each pixel")
    # Determine number of files
    numfiles = np.unique(all_idx).size

    # Find the location of the object with the highest S/N in the combined white light image
    idx_max = np.unravel_index(np.argmax(whitelight_img), whitelight_img.shape)
    msgs.info("Highest S/N object located at spaxel (x, y) = {0:d}, {1:d}".format(idx_max[0], idx_max[1]))

    # Generate a 2D WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, dwv]
    whitelightWCS = generate_WCS(coord_min, coord_dlt)
    # Make the bin edges to be at +/- 1 pixels around the maximum (i.e. summing 9 pixels total)
    numwav = int((np.max(all_wave) - np.min(all_wave)) / dwv)
    xbins = np.array([idx_max[0]-1, idx_max[0]+2]) - 0.5
    ybins = np.array([idx_max[1]-1, idx_max[1]+2]) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    bins = (xbins, ybins, spec_bins)

    # Extract the spectrum of the highest S/N object
    flux_stack = np.zeros((numwav, numfiles))
    ivar_stack = np.zeros((numwav, numfiles))
    for ff in range(numfiles):
        msgs.info("Extracting spectrum of highest S/N detection from frame {0:d}/{1:d}".format(ff + 1, numfiles))
        ww = (all_idx == ff)
        # Extract the spectrum
        pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
        spec, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci[ww])
        var, edges = np.histogramdd(pix_coord, bins=bins, weights=1/all_ivar[ww])
        norm, edges = np.histogramdd(pix_coord, bins=bins)
        normspec = (norm > 0) / (norm + (norm == 0))
        var_spec = var[0, 0, :]
        ivar_spec = (var_spec > 0) / (var_spec + (var_spec == 0))
        # Calculate the S/N in a given spectral bin
        flux_stack[:, ff] = spec[0, 0, :] * np.sqrt(normspec)  # Note: sqrt(nrmspec), is because we want the S/N in a _single_ pixel (i.e. not spectral bin)
        ivar_stack[:, ff] = ivar_spec

    mask_stack = (flux_stack != 0.0) & (ivar_stack != 0.0)
    # Obtain a wavelength of each pixel
    wcs_res = whitelightWCS.wcs_pix2world(np.vstack((np.zeros(numwav), np.zeros(numwav), np.arange(numwav))).T, 0)
    wcs_scale = (1.0 * whitelightWCS.wcs.cunit[2]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
    wave_spec = wcs_scale * wcs_res[:, 2]
    # Compute the smoothing scale to use
    if sn_smooth_npix is None:
        sn_smooth_npix = int(np.round(0.1 * wave_spec.size))
    rms_sn, weights = coadd.sn_weights(utils.array_to_explist(flux_stack), utils.array_to_explist(ivar_stack), utils.array_to_explist(mask_stack),
                                       sn_smooth_npix=sn_smooth_npix, relative_weights=relative_weights)

    # Because we pass back a weights array, we need to interpolate to assign each detector pixel a weight
    all_wghts = np.ones(all_idx.size)
    for ff in range(numfiles):
        ww = (all_idx == ff)
        all_wghts[ww] = interp1d(wave_spec, weights[ff], kind='cubic',
                                 bounds_error=False, fill_value="extrapolate")(all_wave[ww])
    msgs.info("Optimal weighting complete")
    return all_wghts
