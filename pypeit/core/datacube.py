"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy import signal
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit.core import coadd, flux_calib

# Use a fast histogram for speed!
from fast_histogram import histogramdd

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
        ss_file (:obj:`str`):
            The relative path and filename of the standard star datacube. It
            should be fits format, and for full functionality, should ideally of
            the form :class:`~pypeit.coadd3d.DataCube`.
        senspar (:class:`~pypeit.par.pypeitpar.SensFuncPar`):
            The parameters required for the sensitivity function computation.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array used to construct blaze_spline
        blaze_spline (`scipy.interpolate.interp1d`_, optional):
            Spline representation of the reference blaze function (based on the illumflat).
        grating_corr (:obj:`bool`, optional):
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


def get_whitelight_pixels(all_wave, all_slitid, min_wl, max_wl):
    """
    Determine which pixels are included within the specified wavelength range

    Args:
        all_wave (`numpy.ndarray`_, list):
            List of `numpy.ndarray`_ wavelength images. The length of the list is the number of spec2d frames.
            Each element of the list contains a wavelength image that provides the wavelength at each pixel on
            the detector, with shape is (nspec, nspat).
        all_slitid (`numpy.ndarray`_, list):
            List of `numpy.ndarray`_ slitid images. The length of the list is the number of spec2d frames.
            Each element of the list contains a slitid image that provides the slit number at each pixel on
            the detector, with shape (nspec, nspat).
        min_wl (float):
            Minimum wavelength to consider
        max_wl (float):
            Maximum wavelength to consider

    Returns:
        :obj:`tuple`: The first element of the tuple is a list of `numpy.ndarray`_ slitid images
        (or a single `numpy.ndarray`_ slitid image if only one spec2d frame is provided),
        shape is (nspec, nspat), where a zero value corresponds to an excluded pixel
        (either outside the desired wavelength range, a bad pixel, a pixel not on the slit).
        All other pixels have a value equal to the slit number. The second element of the tuple
        is the wavelength difference between the maximum and minimum wavelength in the desired
        wavelength range.
    """
    # Check if lists or numpy arrays are input
    list_inputs = [all_wave, all_slitid]
    if all([isinstance(l, list) for l in list_inputs]):
        numframes = len(all_wave)
        if not all([len(l) == numframes for l in list_inputs]):
            msgs.error("All input lists must have the same length")
        # Store in the following variables
        _all_wave, _all_slitid = all_wave, all_slitid
    elif all([not isinstance(l, list) for l in list_inputs]):
        _all_wave, _all_slitid = [all_wave], [all_slitid]
        numframes = 1
    else:
        msgs.error("The input lists must either all be lists (of the same length) or all be numpy arrays")
    if max_wl < min_wl:
        msgs.error("The maximum wavelength must be greater than the minimum wavelength")
    # Initialise the output
    out_slitid = [np.zeros(_all_slitid[0].shape, dtype=int) for _ in range(numframes)]
    # Loop over all frames and find the pixels that are within the wavelength range
    if min_wl < max_wl:
        # Loop over files and determine which pixels are within the wavelength range
        for ff in range(numframes):
            ww = np.where((_all_wave[ff] > min_wl) & (_all_wave[ff] < max_wl))
            out_slitid[ff][ww] = _all_slitid[ff][ww]
    else:
        msgs.warn("Datacubes do not completely overlap in wavelength.")
        out_slitid = _all_slitid
        min_wl, max_wl = None, None
        for ff in range(numframes):
            this_wave = _all_wave[ff][_all_slitid[ff] > 0]
            tmp_min = np.min(this_wave)
            tmp_max = np.max(this_wave)
            if min_wl is None or tmp_min < min_wl:
                min_wl = tmp_min
            if max_wl is None or tmp_max > max_wl:
                max_wl = tmp_max
    # Determine the wavelength range
    wavediff = max_wl - min_wl
    # Need to return a single slitid image if only one frame, otherwise return a list of slitid images.
    # Also return the wavelength difference
    return out_slitid[0] if numframes == 1 else out_slitid, wavediff


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


def align_user_offsets(ifu_ra, ifu_dec, ra_offset, dec_offset):
    """
    Align the RA and DEC of all input frames, and then
    manually shift the cubes based on user-provided offsets.
    The offsets should be specified in arcseconds, and the
    ra_offset should include the cos(dec) factor.

    Args:
        ifu_ra (`numpy.ndarray`_):
            A list of RA values of the IFU (one value per frame)
        ifu_dec (`numpy.ndarray`_):
            A list of Dec values of the IFU (one value per frame)
        ra_offset (`numpy.ndarray`_):
            A list of RA offsets to be applied to the input pixel values (one value per frame).
            Note, the ra_offset MUST contain the cos(dec) factor. This is the number of degrees
            on the sky that represents the telescope offset.
        dec_offset (`numpy.ndarray`_):
            A list of Dec offsets to be applied to the input pixel values (one value per frame).
            This is the number of degrees on the sky that represents the telescope offset.

    Returns:
        A tuple containing a new set of RA and Dec offsets for each frame.
        Both arrays are of type `numpy.ndarray`_, and are in units of degrees.
    """
    # First, translate all coordinates to the coordinates of the first frame
    # Note: You do not need cos(dec) here, this just overrides the IFU coordinate centre of each frame
    #       The cos(dec) factor should be input by the user, and should be included in the self.opts['ra_offset']
    ref_shift_ra = ifu_ra[0] - ifu_ra
    ref_shift_dec = ifu_dec[0] - ifu_dec
    numfiles = len(ra_offset)
    out_ra_offsets = [0.0 for _ in range(numfiles)]
    out_dec_offsets = [0.0 for _ in range(numfiles)]
    for ff in range(numfiles):
        # Apply the shift
        out_ra_offsets[ff] = ref_shift_ra[ff] + ra_offset[ff]
        out_dec_offsets[ff] = ref_shift_dec[ff] + dec_offset[ff]
        msgs.info("Spatial shift of cube #{0:d}:".format(ff + 1) + msgs.newline() +
                  "RA, DEC (arcsec) = {0:+0.3f} E, {1:+0.3f} N".format(ra_offset[ff]*3600.0, dec_offset[ff]*3600.0))
    return out_ra_offsets, out_dec_offsets


def set_voxel_sampling(spatscale, specscale, dspat=None, dwv=None):
    """
    This function checks if the spatial and spectral scales of all frames are consistent.
    If the user has not specified either the spatial or spectral scales, they will be set here.

    Parameters
    ----------
    spatscale : `numpy.ndarray`_
        2D array, shape is (N, 2), listing the native spatial scales of N spec2d frames.
        spatscale[:,0] refers to the spatial pixel scale of each frame
        spatscale[:,1] refers to the slicer scale of each frame
        Each element of the array must be in degrees
    specscale : `numpy.ndarray`_
        1D array listing the native spectral scales of multiple frames. The length of this array should be equal
        to the number of frames you are using. Each element of the array must be in Angstrom
    dspat: :obj:`float`, optional
        Spatial scale to use as the voxel spatial sampling. If None, a new value will be derived based on the inputs
    dwv: :obj:`float`, optional
        Spectral scale to use as the voxel spectral sampling. If None, a new value will be derived based on the inputs

    Returns
    -------
    _dspat : :obj:`float`
        Spatial sampling
    _dwv : :obj:`float`
        Wavelength sampling
    """
    # Make sure all frames have consistent pixel scales
    ratio = (spatscale[:, 0] - spatscale[0, 0]) / spatscale[0, 0]
    if np.any(np.abs(ratio) > 1E-4):
        msgs.warn("The pixel scales of all input frames are not the same!")
        spatstr = ", ".join(["{0:.6f}".format(ss) for ss in spatscale[:,0]*3600.0])
        msgs.info("Pixel scales of all input frames:" + msgs.newline() + spatstr + "arcseconds")
    # Make sure all frames have consistent slicer scales
    ratio = (spatscale[:, 1] - spatscale[0, 1]) / spatscale[0, 1]
    if np.any(np.abs(ratio) > 1E-4):
        msgs.warn("The slicer scales of all input frames are not the same!")
        spatstr = ", ".join(["{0:.6f}".format(ss) for ss in spatscale[:,1]*3600.0])
        msgs.info("Slicer scales of all input frames:" + msgs.newline() + spatstr + "arcseconds")
    # Make sure all frames have consistent wavelength sampling
    ratio = (specscale - specscale[0]) / specscale[0]
    if np.any(np.abs(ratio) > 1E-2):
        msgs.warn("The wavelength samplings of the input frames are not the same!")
        specstr = ", ".join(["{0:.6f}".format(ss) for ss in specscale])
        msgs.info("Wavelength samplings of all input frames:" + msgs.newline() + specstr + "Angstrom")

    # If the user has not specified the spatial scale, then set it appropriately now to the largest spatial scale
    _dspat = np.max(spatscale) if dspat is None else dspat
    msgs.info("Adopting a square pixel spatial scale of {0:f} arcsec".format(3600.0 * _dspat))
    # If the user has not specified the spectral sampling, then set it now to the largest value
    _dwv = np.max(specscale) if dwv is None else dwv
    msgs.info("Adopting a wavelength sampling of {0:f} Angstrom".format(_dwv))
    return _dspat, _dwv


def check_inputs(list_inputs):
    """
    This function checks the inputs to several of the cube building routines, and makes sure they are all consistent.
    Often, this is to make check if all inputs are lists of the same length, or if all inputs are 2D `numpy.ndarray`_.
    The goal of the routine is to return a consistent set of lists of the input.

    Parameters
    ----------
    list_inputs : :obj:`list`
        A list of inputs to check.

    Returns
    -------
    list_inputs : :obj:`list`
        A list of inputs that have been checked for consistency.
    """
    if all([isinstance(l, list) for l in list_inputs]):
        # Several frames are being combined. Check the lists have the same length
        numframes = len(list_inputs[0])
        if not all([len(l) == numframes for l in list_inputs]):
            msgs.error("All input lists must have the same length")
        # The inputs are good, return as is
        return tuple(list_inputs)
    elif all([not isinstance(l, list) for l in list_inputs]):
        # Just a single frame - store as single element lists
        ret_list = ()
        for l in list_inputs:
            ret_list += ([l],)
        return ret_list
    else:
        msgs.error("The input arguments should all be of type 'list', or all not be of type 'list':")


def wcs_bounds(raImg, decImg, waveImg, slitid_img_gpm, ra_offsets=None, dec_offsets=None,
               ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None):
    """
    Calculate the bounds of the WCS and the expected edges of the voxels, based
    on user-specified parameters or the extremities of the data. This is a
    convenience function that calls the core function in
    :mod:`~pypeit.core.datacube`.

    Parameters
    ----------
    raImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg (`numpy.ndarray`_, list):
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    slitid_img_gpm : (`numpy.ndarray`_, list):
        A list of 2D array containing the spat ID of each pixel, with shape (nspec, nspat).
        A value of 0 indicates that the pixel is not on a slit. All other values indicate the
        slit spatial ID.
    ra_offsets : list, optional
        A list of the RA offsets for each frame
    dec_offsets : list, optional
        A list of the Dec offsets for each frame
    ra_min : :obj:`float`, optional
        Minimum RA of the WCS
    ra_max : :obj:`float`, optional
        Maximum RA of the WCS
    dec_min : :obj:`float`, optional
        Minimum Dec of the WCS
    dec_max : :obj:`float`, optional
        Maximum Dec of the WCS
    wave_min : :obj:`float`, optional
        Minimum wavelength of the WCS
    wave_max : :obj:`float`, optional
        Maximum wavelength of the WCS

    Returns
    -------
    _ra_min : :obj:`float`
        Minimum RA of the WCS
    _ra_max : :obj:`float`
        Maximum RA of the WCS
    _dec_min : :obj:`float`
        Minimum Dec of the WCS
    _dec_max : :obj:`float`
        Maximum Dec of the WCS
    _wave_min : :obj:`float`
        Minimum wavelength of the WCS
    _wave_max : :obj:`float`
        Maximum wavelength of the WCS
    """
    # Check if the ra_offsets and dec_offsets are specified
    if ra_offsets is None or dec_offsets is None:
        if isinstance(raImg, list):
            ra_offsets = [0.0]*len(raImg)
            dec_offsets = [0.0]*len(raImg)
        else:
            ra_offsets = 0.0
            dec_offsets = 0.0
    # Check the inputs
    _raImg, _decImg, _waveImg, _slitid_img_gpm, _ra_offsets, _dec_offsets = \
        check_inputs([raImg, decImg, waveImg, slitid_img_gpm, ra_offsets, dec_offsets])
    numframes = len(_raImg)

    # Loop over the frames and get the bounds - start by setting the default values
    _ra_min, _ra_max = ra_min, ra_max
    _dec_min, _dec_max = dec_min, dec_max
    _wave_min, _wave_max = wave_min, wave_max
    for fr in range(numframes):
        # Only do calculations if the min/max inputs are not specified
        # Get the RA, Dec, and wavelength of the pixels on the slit
        if ra_min is None or ra_max is None:
            this_ra = _raImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min, tmp_max = np.min(this_ra)+_ra_offsets[fr], np.max(this_ra)+_ra_offsets[fr]
            if fr == 0 or tmp_min < _ra_min:
                _ra_min = tmp_min
            if fr == 0 or tmp_max > _ra_max:
                _ra_max = tmp_max
        if dec_min is None or dec_max is None:
            this_dec = _decImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min, tmp_max = np.min(this_dec)+_dec_offsets[fr], np.max(this_dec)+_dec_offsets[fr]
            if fr == 0 or tmp_min < _dec_min:
                _dec_min = tmp_min
            if fr == 0 or tmp_max > _dec_max:
                _dec_max = tmp_max
        if wave_min is None or wave_max is None:
            this_wave = _waveImg[fr][_slitid_img_gpm[fr] > 0]
            tmp_min, tmp_max = np.min(this_wave), np.max(this_wave)
            if fr == 0 or tmp_min < _wave_min:
                _wave_min = tmp_min
            if fr == 0 or tmp_max > _wave_max:
                _wave_max = tmp_max
    # Return the bounds
    return _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max


def create_wcs(raImg, decImg, waveImg, slitid_img_gpm, dspat, dwave,
               ra_offsets=None, dec_offsets=None,
               ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
               reference=None, collapse=False, equinox=2000.0, specname="PYP_SPEC"):
    """
    Create a WCS and the expected edges of the voxels, based on user-specified
    parameters or the extremities of the data.

    Parameters
    ----------
    raImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : (`numpy.ndarray`_, list):
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg (`numpy.ndarray`_, list):
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    slitid_img_gpm : (`numpy.ndarray`_, list):
        A list of 2D array containing the spat ID of each pixel, with shape (nspec, nspat).
        A value of 0 indicates that the pixel is not on a slit. All other values indicate the
        slit spatial ID.
    dspat : float
        Spatial size of each square voxel (in arcsec). The default is to use the
        values in cubepar.
    dwave : float
        Linear wavelength step of each voxel (in Angstroms)
    ra_offsets : list, optional
        List of RA offsets for each frame (degrees)
    dec_offsets : list, optional
        List of Dec offsets for each frame (degrees)
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
    # Setup the cube ranges
    _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max = \
        wcs_bounds(raImg, decImg, waveImg, slitid_img_gpm,
                   ra_offsets=ra_offsets, dec_offsets=dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max)

    # Grab cos(dec) for convenience. Use the average of the min and max dec
    cosdec = np.cos(0.5*(_dec_min+_dec_max) * np.pi / 180.0)

    # Number of voxels in each dimension
    numra = int((_ra_max - _ra_min) * cosdec / dspat)
    numdec = int((_dec_max - _dec_min) / dspat)
    numwav = int(np.round((_wave_max - _wave_min) / dwave))

    # If a white light WCS is being generated, make sure there's only 1 wavelength bin
    if collapse:
        dwave = _wave_max - _wave_min
        numwav = 1

    # Generate a master WCS to register all frames
    coord_min = [_ra_min, _dec_min, _wave_min]
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
              msgs.newline() + "WAVE min, max = {0:f}, {1:f}".format(_wave_min, _wave_max) +
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


def compute_weights_frompix(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg, dspat, dwv, mnmx_wv, wghtsImg,
                            all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                            ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
                            sn_smooth_npix=None, weight_method='auto', reference_image=None, whitelight_range=None,
                            correct_dar=True, specname="PYPSPEC"):
    r"""
    Calculate wavelength dependent optimal weights. The weighting is currently
    based on a relative :math:`(S/N)^2` at each wavelength. Note, this function
    first prepares a whitelight image, and then calls compute_weights() to
    determine the appropriate weights of each pixel.

    Parameters
    ----------

    raImg : `numpy.ndarray`_, list
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : `numpy.ndarray`_, list
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg : `numpy.ndarray`_, list
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    sciImg : `numpy.ndarray`_, list
        A list of 2D array containing the science image of each pixel, with shape (nspec, nspat)
    ivarImg : `numpy.ndarray`_, list
        A list of 2D array containing the inverse variance image of each pixel, with shape (nspec, nspat)
    slitidImg : `numpy.ndarray`_, list
        A list of 2D array containing the slit ID of each pixel, with shape (nspec, nspat)
    dspat : float
        The size of each spaxel on the sky (in degrees)
    dwv : float
        The size of each wavelength pixel (in Angstroms)
    mnmx_wv : `numpy.ndarray`_
        The minimum and maximum wavelengths of every slit and frame. The shape is (Nframes, Nslits, 2),
        The minimum and maximum wavelengths are stored in the [:,:,0] and [:,:,1] indices, respectively.
    wghtsImg : `numpy.ndarray`_, list
        A list of 2D array containing the weights of each pixel, with shape (nspec, nspat)
    all_wcs : `astropy.wcs.WCS`_, list
        A list of WCS objects, one for each frame.
    all_tilts : `numpy.ndarray`_, list
        2D wavelength tilts frame, or a list of tilt frames
    all_slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        Information stored about the slits, or a list of SlitTraceSet objects
    all_align : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel
        coordinates and WCS pixel coordinates, or a list of Alignment
        Splines.
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of DARcorrection
        classes. If a list, it must be the same length as astrom_trans.
    ra_offsets : float, list
        RA offsets for each frame in units of degrees
    dec_offsets : float, list
        Dec offsets for each frame in units of degrees
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
    sn_smooth_npix : float, optional
        Number of pixels used for determining smoothly varying S/N ratio
        weights.  This is currently not required, since a relative weighting
        scheme with a polynomial fit is used to calculate the S/N weights.
    weight_method : `str`, optional
        Weight method to be used in :func:`~pypeit.coadd.sn_weights`.
        Options are ``'auto'``, ``'constant'``, ``'uniform'``, ``'wave_dependent'``, ``'relative'``, or
        ``'ivar'``. The default is ``'auto'``.  Behavior is as follows:

            - ``'auto'``: Use constant weights if rms_sn < 3.0, otherwise
                use wavelength dependent.

            - ``'constant'``: Constant weights based on rms_sn**2

            - ``'uniform'``: Uniform weighting.

            - ``'wave_dependent'``: Wavelength dependent weights will be
                used irrespective of the rms_sn ratio. This option will not
                work well at low S/N ratio although it is useful for objects
                where only a small fraction of the spectral coverage has high
                S/N ratio (like high-z quasars).

            - ``'relative'``: Calculate weights by fitting to the ratio of
                spectra? Note, relative weighting will only work well when
                there is at least one spectrum with a reasonable S/N, and a
                continuum.  RJC note - This argument may only be better when
                the object being used has a strong continuum + emission lines.
                The reference spectrum is assigned a value of 1 for all
                wavelengths, and the weights of all other spectra will be
                determined relative to the reference spectrum. This is
                particularly useful if you are dealing with highly variable
                spectra (e.g. emission lines) and require a precision better
                than ~1 per cent.

            - ``'ivar'``: Use inverse variance weighting. This is not well
                tested and should probably be deprecated.

    reference_image : `numpy.ndarray`_
        Reference image to use for the determination of the highest S/N spaxel in the image.
    correct_dar : bool, optional
        Correct for the differential atmospheric refraction.  Default is False.
    specname : str
        Name of the spectrograph

    Returns
    -------
    weights : `numpy.ndarray`_ 
        a 1D array the same size as all_sci, containing relative wavelength
        dependent weights of each input pixel.
    """
    # Find the wavelength range where all frames overlap
    min_wl, max_wl = get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                          np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                          whitelight_range)  # The user-specified values (if any)
    # Get the good white light pixels
    slitid_img_gpm, wavediff = get_whitelight_pixels(waveImg, slitidImg, min_wl, max_wl)

    # Generate the WCS
    image_wcs, voxedge, reference_image = \
        create_wcs(raImg, decImg, waveImg, slitid_img_gpm, dspat, wavediff,
                   ra_offsets=ra_offsets, dec_offsets=dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max,
                   reference=reference_image, collapse=True, equinox=2000.0, specname=specname)

    # Generate the white light image
    # NOTE: hard-coding subpixel=1 in both directions for speed, and combining into a single image
    wl_full = generate_image_subpixel(image_wcs, voxedge, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtsImg,
                                      all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                                      spec_subpixel=1, spat_subpixel=1, slice_subpixel=1, combine=True,
                                      correct_dar=correct_dar)

    # Compute the weights
    return compute_weights(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                           all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                           wl_full[:, :, 0], dspat, dwv,
                           ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min,
                           sn_smooth_npix=sn_smooth_npix, weight_method=weight_method, correct_dar=correct_dar)


def compute_weights(raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                    all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets,
                    whitelight_img, dspat, dwv,
                    ra_min=None, ra_max=None, dec_min=None, dec_max=None, wave_min=None, wave_max=None,
                    sn_smooth_npix=None, weight_method='auto', correct_dar=True):
    r"""
    Calculate wavelength dependent optimal weights. The weighting is currently
    based on a relative :math:`(S/N)^2` at each wavelength

    Parameters
    ----------

    raImg : `numpy.ndarray`_, list
        A list of 2D array containing the RA of each pixel, with shape (nspec, nspat)
    decImg : `numpy.ndarray`_, list
        A list of 2D array containing the Dec of each pixel, with shape (nspec, nspat)
    waveImg : `numpy.ndarray`_, list
        A list of 2D array containing the wavelength of each pixel, with shape (nspec, nspat)
    sciImg : `numpy.ndarray`_, list
        A list of 2D array containing the science image of each pixel, with shape (nspec, nspat)
    ivarImg : `numpy.ndarray`_, list
        A list of 2D array containing the inverse variance image of each pixel, with shape (nspec, nspat)
    slitidImg : `numpy.ndarray`_, list
        A list of 2D array containing the slit ID of each pixel, with shape (nspec, nspat)
    all_wcs : `astropy.wcs.WCS`_, list
        A list of WCS objects, one for each frame.
    all_tilts : `numpy.ndarray`_, list
        2D wavelength tilts frame, or a list of tilt frames
    all_slits : :class:`~pypeit.slittrace.SlitTraceSet`, list
        Information stored about the slits, or a list of SlitTraceSet objects
    all_align : :class:`~pypeit.alignframe.AlignmentSplines`, list
        A Class containing the transformation between detector pixel
        coordinates and WCS pixel coordinates, or a list of Alignment
        Splines.
    all_dar : :class:`~pypeit.coadd3d.DARcorrection`, list
        A Class containing the DAR correction information, or a list of DARcorrection
        classes. If a list, it must be the same length as astrom_trans.
    ra_offsets : float, list
        RA offsets for each frame in units of degrees
    dec_offsets : float, list
        Dec offsets for each frame in units of degrees
    whitelight_img : `numpy.ndarray`_
        A 2D array containing a white light image, that was created with the
        input ``all`` arrays.
    dspat : float
        The size of each spaxel on the sky (in degrees)
    dwv : float
        The size of each wavelength pixel (in Angstroms)
    sn_smooth_npix : float, optional
        Number of pixels used for determining smoothly varying S/N ratio
        weights.  This is currently not required, since a relative weighting
        scheme with a polynomial fit is used to calculate the S/N weights.
    correct_dar : bool, optional
        Apply the DAR correction to the input data.  The default is True.
    weight_method : `str`, optional
        Weight method to be used in :func:`~pypeit.coadd.sn_weights`.
        Options are ``'auto'``, ``'constant'``, ``'uniform'``, ``'wave_dependent'``, ``'relative'``, or
        ``'ivar'``. The default is ``'auto'``.  Behavior is as follows:

            - ``'auto'``: Use constant weights if rms_sn < 3.0, otherwise
                use wavelength dependent.

            - ``'constant'``: Constant weights based on rms_sn**2

            - ``'uniform'``: Uniform weighting.

            - ``'wave_dependent'``: Wavelength dependent weights will be
                used irrespective of the rms_sn ratio. This option will not
                work well at low S/N ratio although it is useful for objects
                where only a small fraction of the spectral coverage has high
                S/N ratio (like high-z quasars).

            - ``'relative'``: Calculate weights by fitting to the ratio of
                spectra? Note, relative weighting will only work well when
                there is at least one spectrum with a reasonable S/N, and a
                continuum.  RJC note - This argument may only be better when
                the object being used has a strong continuum + emission lines.
                The reference spectrum is assigned a value of 1 for all
                wavelengths, and the weights of all other spectra will be
                determined relative to the reference spectrum. This is
                particularly useful if you are dealing with highly variable
                spectra (e.g. emission lines) and require a precision better
                than ~1 per cent.

            - ``'ivar'``: Use inverse variance weighting. This is not well
                tested and should probably be deprecated.

    Returns
    -------
    all_wghts: list
        Either a 2D `numpy.ndarray`_ or a list of 2D `numpy.ndarray`_ arrays
        containing the optimal weights of each pixel for all frames, with shape
        (nspec, nspat).
    """
    msgs.info("Calculating the optimal weights of each pixel")
    # Check the inputs for combinations of lists or not, and then determine the number of frames
    _raImg, _decImg, _waveImg, _sciImg, _ivarImg, _slitidImg, \
        _all_wcs, _all_tilts, _all_slits, _all_align, _all_dar, _ra_offsets, _dec_offsets = \
            check_inputs([raImg, decImg, waveImg, sciImg, ivarImg, slitidImg,
                          all_wcs, all_tilts, all_slits, all_align, all_dar, ra_offsets, dec_offsets])
    numframes = len(_sciImg)

    # If there's only one frame, use uniform weighting
    if numframes == 1:
        msgs.warn("Only one frame provided.  Using uniform weighting.")
        return np.ones_like(sciImg)

    # Check the WCS bounds
    _ra_min, _ra_max, _dec_min, _dec_max, _wave_min, _wave_max = \
        wcs_bounds(_raImg, _decImg, _waveImg, _slitidImg, ra_offsets=_ra_offsets, dec_offsets=_dec_offsets,
                   ra_min=ra_min, ra_max=ra_max, dec_min=dec_min, dec_max=dec_max, wave_min=wave_min, wave_max=wave_max)

    # Find the location of the object with the highest S/N in the combined white light image
    med_filt_whitelight = signal.medfilt2d(whitelight_img, kernel_size=3)
    idx_max = np.unravel_index(np.argmax(med_filt_whitelight), med_filt_whitelight.shape)
    # TODO: Taking the maximum pixel of the whitelight image is extremely brittle to the case where
    #  their are hot pixels in the white light image, which there are plenty of since the edges of the slits are very
    #  poorly behaved.
    #idx_max = np.unravel_index(np.argmax(whitelight_img), whitelight_img.shape)
    msgs.info("Highest S/N object located at spaxel (x, y) = {0:d}, {1:d}".format(idx_max[0], idx_max[1]))

    # Generate a 2D WCS to register all frames
    coord_min = [_ra_min, _dec_min, _wave_min]
    coord_dlt = [dspat, dspat, dwv]
    whitelightWCS = generate_WCS(coord_min, coord_dlt)
    wcs_scale = (1.0 * whitelightWCS.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
    # Make the bin edges to be at +/- 1 pixels around the maximum (i.e. summing 9 pixels total)
    numwav = int((_wave_max - _wave_min) / dwv)
    xbins = np.array([idx_max[0]-1, idx_max[0]+2]) - 0.5
    ybins = np.array([idx_max[1]-1, idx_max[1]+2]) - 0.5
    spec_bins = np.arange(1 + numwav) - 0.5
    bins = (xbins, ybins, spec_bins)

    # Extract the spectrum of the highest S/N object
    flux_stack = np.zeros((numwav, numframes))
    ivar_stack = np.zeros((numwav, numframes))
    for ff in range(numframes):
        msgs.info("Extracting spectrum of highest S/N detection from frame {0:d}/{1:d}".format(ff + 1, numframes))
        flxcube, sigcube, bpmcube, wave = \
            generate_cube_subpixel(whitelightWCS, bins, _sciImg[ff], _ivarImg[ff], _waveImg[ff],
                                   _slitidImg[ff], np.ones(_sciImg[ff].shape), _all_wcs[ff],
                                   _all_tilts[ff], _all_slits[ff], _all_align[ff], _all_dar[ff],
                                   _ra_offsets[ff], _dec_offsets[ff],
                                   spec_subpixel=1, spat_subpixel=1, slice_subpixel=1,
                                   skip_subpix_weights=True, correct_dar=correct_dar)
        # Store the flux and ivar spectra of the highest S/N object.
        # TODO :: This is the flux per spectral pixel, and not per detector pixel.  Is this correct?
        flux_stack[:, ff] = flxcube[:, 0, 0]
        ivar_stack[:, ff] = utils.inverse(sigcube[:, 0, 0])**2

    # Mask out any pixels that are zero in the flux or ivar stack
    mask_stack = (flux_stack != 0.0) & (ivar_stack != 0.0)
    # Obtain a wavelength of each pixel
    wcs_res = whitelightWCS.wcs_pix2world(np.vstack((np.zeros(numwav), np.zeros(numwav), np.arange(numwav))).T, 0)
    wcs_scale = (1.0 * whitelightWCS.wcs.cunit[2]).to_value(units.Angstrom)  # Ensures the WCS is in Angstroms
    wave_spec = wcs_scale * wcs_res[:, 2]
    # Compute the smoothing scale to use
    if sn_smooth_npix is None:
        sn_smooth_npix = int(np.round(0.1 * wave_spec.size))
    rms_sn, weights = coadd.sn_weights(utils.array_to_explist(flux_stack), utils.array_to_explist(ivar_stack), utils.array_to_explist(mask_stack),
                                       sn_smooth_npix=sn_smooth_npix, weight_method=weight_method)

    # Because we pass back a weights array, we need to interpolate to assign each detector pixel a weight
    all_wghts = [np.ones(_sciImg[0].shape) for _ in range(numframes)]
    for ff in range(numframes):
        ww = (slitidImg[ff] > 0)
        all_wghts[ff][ww] = interp1d(wave_spec, weights[ff], kind='cubic',
                                 bounds_error=False, fill_value="extrapolate")(waveImg[ff][ww])
    msgs.info("Optimal weighting complete")
    return all_wghts


def generate_image_subpixel(image_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
                            all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset,
                            spec_subpixel=5, spat_subpixel=5, slice_subpixel=5, combine=False, correct_dar=True):
    """
    Generate a white light image from the input pixels

    Args:
        image_wcs (`astropy.wcs.WCS`_):
            World coordinate system to use for the white light images.
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        sciImg (`numpy.ndarray`_, list):
            A list of 2D science images, or a single 2D image containing the
            science data.
        ivarImg (`numpy.ndarray`_, list):
            A list of 2D inverse variance images, or a single 2D image
            containing the inverse variance data.
        waveImg (`numpy.ndarray`_, list):
            A list of 2D wavelength images, or a single 2D image containing the
            wavelength data.
        slitid_img_gpm (`numpy.ndarray`_, list):
            A list of 2D slit ID images, or a single 2D image containing the
            slit ID data.
        wghtImg (`numpy.ndarray`_, list):
            A list of 2D weight images, or a single 2D image containing the
            weight data.
        all_wcs (`astropy.wcs.WCS`_, list):
            A list of WCS objects, or a single WCS object containing the WCS
            information of each image.
        tilts (`numpy.ndarray`_, list):
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`~pypeit.slittrace.SlitTraceSet`, list):
            Information stored about the slits, or a list of SlitTraceSet (see
            all_idx)
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        all_dar (:class:`~pypeit.coadd3d.DARcorrection`, list):
            A Class containing the DAR correction information, or a list of DARcorrection
            classes. If a list, it must be the same length as astrom_trans.
        ra_offset (:obj:`float`, list):
            The RA offset to apply to each image, or a list of RA offsets.
        dec_offset (:obj:`float`, list):
            The DEC offset to apply to each image, or a list of DEC offsets.
        spec_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spectral direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
            default value is 5, which divides each detector pixel into 5 subpixels
            in the spectral direction.
        spat_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spatial direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
            default value is 5, which divides each detector pixel into 5 subpixels
            in the spatial direction.
        slice_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the slice direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel * slice_subpixel``). The
            default value is 5, which divides each IFU slice into 5 subpixels
            in the slice direction.
        combine (:obj:`bool`, optional):
            If True, all of the input frames will be combined into a single
            output. Otherwise, individual images will be generated.
        correct_dar (:obj:`bool`, optional):
            If True, the DAR correction will be applied to the input images
            before generating the white light images. If False, the DAR
            correction will not be applied.

    Returns:
        `numpy.ndarray`_: The white light images for all frames
    """
    # Perform some checks on the input -- note, more complete checks are performed in subpixellate()
    _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg, _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset = \
        check_inputs([sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset])
    numframes = len(_sciImg)

    # Prepare the array of white light images to be stored
    numra = bins[0].size-1
    numdec = bins[1].size-1
    all_wl_imgs = np.zeros((numra, numdec, numframes))

    # Loop through all frames and generate white light images
    for fr in range(numframes):
        msgs.info(f"Creating image {fr+1}/{numframes}")
        if combine:
            # Subpixellate
            img, _, _ = subpixellate(image_wcs, bins, _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg,
                                     _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset,
                                     spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
                                     skip_subpix_weights=True, correct_dar=correct_dar)
        else:
            # Subpixellate
            img, _, _ = subpixellate(image_wcs, bins, _sciImg[fr], _ivarImg[fr], _waveImg[fr], _slitid_img_gpm[fr], _wghtImg[fr],
                                     _all_wcs[fr], _tilts[fr], _slits[fr], _astrom_trans[fr], _all_dar[fr], _ra_offset[fr], _dec_offset[fr],
                                     spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
                                     skip_subpix_weights=True, correct_dar=correct_dar)
        all_wl_imgs[:, :, fr] = img[:, :, 0]
    # Return the constructed white light images
    return all_wl_imgs


def generate_cube_subpixel(output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
                           all_wcs, tilts, slits, astrom_trans, all_dar,
                           ra_offset, dec_offset,
                           spec_subpixel=5, spat_subpixel=5, slice_subpixel=5, skip_subpix_weights=False,
                           overwrite=False, outfile=None, whitelight_range=None, correct_dar=True):
    """
    Save a datacube using the subpixel algorithm. Refer to the subpixellate()
    docstring for further details about this algorithm

    Args:
        output_wcs (`astropy.wcs.WCS`_):
            Output world coordinate system.
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        sciImg (`numpy.ndarray`_, list):
            A list of 2D array containing the counts of each pixel. If a list,
            the shape of each numpy array is (nspec, nspat).
        ivarImg (`numpy.ndarray`_, list):
            A list of 2D array containing the inverse variance of each pixel. If a list,
            the shape of each numpy array is (nspec, nspat).
        waveImg (`numpy.ndarray`_, list):
            A list of 2D array containing the wavelength of each pixel. If a list,
            the shape of each numpy array is (nspec, nspat).
        slitid_img_gpm (`numpy.ndarray`_, list):
            A list of 2D array containing the slitmask of each pixel. If a list,
            the shape of each numpy array is (nspec, nspat).
            A zero value indicates that a pixel is either not on a slit or it is a bad pixel.
            All other values are the slit spatial ID number.
        wghtImg (`numpy.ndarray`_, list):
            A list of 2D array containing the weights of each pixel to be used in the
            combination. If a list, the shape of each numpy array is (nspec, nspat).
        all_wcs (`astropy.wcs.WCS`_, list):
            A list of `astropy.wcs.WCS`_ objects, one for each spec2d file.
        tilts (list):
            A list of `numpy.ndarray`_ objects, one for each spec2d file,
            containing the tilts of each pixel. The shape of each numpy array
            is (nspec, nspat).
        slits (:class:`pypeit.slittrace.SlitTraceSet`, list):
            A list of :class:`pypeit.slittrace.SlitTraceSet` objects, one for each
            spec2d file, containing the properties of the slit for each spec2d file
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        all_dar (:class:`~pypeit.coadd3d.DARcorrection`, list):
            A Class containing the DAR correction information, or a list of DARcorrection
            classes. If a list, it must be the same length as astrom_trans.
        ra_offset (float, list):
            A float or list of floats containing the RA offset of each spec2d file
        dec_offset (float, list):
            A float or list of floats containing the DEC offset of each spec2d file
        spec_subpixel (int, optional):
            What is the subpixellation factor in the spectral direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral
            direction.
        spat_subpixel (int, optional):
            What is the subpixellation factor in the spatial direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial
            direction.
        slice_subpixel (int, optional):
            What is the subpixellation factor in the slice direction. Higher
            values give more reliable results, but note that the time required
            goes as (``slice_subpixel``). The default value is 5, which divides
            each IFU slice into 5 subslices in the slice direction.
        skip_subpix_weights (bool, optional):
            If True, the computationally expensive step to calculate the
            subpixellation weights will be skipped. If set the True, note that
            the variance cubes returned will not be accurate. However, if you
            are not interested in the variance cubes, this can save a lot of
            time, and this is an example where you might consider setting this
            variable to True. The flux datacube is unaffected by this variable.
            The default is False.
        overwrite (bool, optional):
            If True, the output cube will be overwritten.
        outfile (str, optional):
            Filename to be used to save the datacube
        whitelight_range (None, list, optional):
            A two element list that specifies the minimum and maximum
            wavelengths (in Angstroms) to use when constructing the white light
            image (format is: [min_wave, max_wave]). If None, the cube will be
            collapsed over the full wavelength range. If a list is provided an
            either element of the list is None, then the minimum/maximum
            wavelength range of that element will be set by the minimum/maximum
            wavelength of all_wave.
        correct_dar (bool, optional):
            If True, the DAR correction will be applied to the datacube. If the
            DAR correction is not available, the datacube will not be corrected.

    Returns:
        :obj:`tuple`: Four `numpy.ndarray`_ objects containing
        (1) the datacube generated from the subpixellated inputs. The shape of
        the datacube is (nwave, nspat1, nspat2).
        (2) the corresponding error cube (standard deviation). The shape of the
        error cube is (nwave, nspat1, nspat2).
        (3) the corresponding bad pixel mask cube. The shape of the bad pixel
        mask cube is (nwave, nspat1, nspat2).
        (4) a 1D array containing the wavelength at each spectral coordinate of the datacube. The
        shape of the wavelength array is (nwave,).
    """
    # Check the inputs
    if whitelight_range is not None and outfile is None:
            msgs.error("Must provide an outfile name if whitelight_range is set")

    # Subpixellate
    flxcube, varcube, bpmcube = subpixellate(output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
                                             all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset,
                                             spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel,
                                             slice_subpixel=slice_subpixel, skip_subpix_weights=skip_subpix_weights,
                                             correct_dar=correct_dar)

    # Get wavelength of each pixel
    nspec = flxcube.shape[2]
    wcs_scale = (1.0*output_wcs.spectral.wcs.cunit[0]).to(units.Angstrom).value  # Ensures the WCS is in Angstroms
    wave = wcs_scale * output_wcs.spectral.wcs_pix2world(np.arange(nspec), 0)[0]

    # Check if the user requested a white light image
    if whitelight_range is not None:
        # Grab the WCS of the white light image
        whitelight_wcs = output_wcs.celestial
        # Determine the wavelength range of the whitelight image
        if whitelight_range[0] is None:
            whitelight_range[0] = wave[0]
        if whitelight_range[1] is None:
            whitelight_range[1] = wave[-1]
        msgs.info("White light image covers the wavelength range {0:.2f} A - {1:.2f} A".format(
            whitelight_range[0], whitelight_range[1]))
        # Get the output filename for the white light image
        out_whitelight = get_output_whitelight_filename(outfile)
        whitelight_img = make_whitelight_fromcube(flxcube, wave=wave, wavemin=whitelight_range[0], wavemax=whitelight_range[1])
        msgs.info("Saving white light image as: {0:s}".format(out_whitelight))
        img_hdu = fits.PrimaryHDU(whitelight_img.T, header=whitelight_wcs.to_header())
        img_hdu.writeto(out_whitelight, overwrite=overwrite)

    # TODO :: Avoid transposing these large cubes
    return flxcube.T, np.sqrt(varcube.T), bpmcube.T, wave


def subpixellate(output_wcs, bins, sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg,
                 all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset,
                 spec_subpixel=5, spat_subpixel=5, slice_subpixel=5, skip_subpix_weights=False,
                 correct_dar=True):
    r"""
    Subpixellate the input data into a datacube. This algorithm splits each
    detector pixel into multiple subpixels and each IFU slice into multiple subslices.
    Then, the algorithm assigns each subdivided detector pixel to a
    voxel. For example, if ``spec_subpixel = spat_subpixel = slice_subpixel = 5``, then each
    detector pixel is divided into :math:`5^3=125` subpixels. Alternatively,
    when spec_subpixel = spat_subpixel = slice_subpixel = 1, this corresponds to the nearest grid
    point (NGP) algorithm.

    Important Note: If spec_subpixel > 1 or spat_subpixel > 1 or slice_subpixel > 1,
    the errors will be correlated, and the covariance is not being tracked, so the
    errors will not be (quite) right. There is a tradeoff one has to make between
    sampling and better looking cubes, versus no sampling and better behaved errors.

    Args:
        output_wcs (`astropy.wcs.WCS`_):
            Output world coordinate system.
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        sciImg (`numpy.ndarray`_, list):
            A list of 2D array containing the counts of each pixel. The shape of
            each 2D array is (nspec, nspat).
        ivarImg (`numpy.ndarray`_, list):
            A list of 2D array containing the inverse variance of each pixel. The shape of
            each 2D array is (nspec, nspat).
        waveImg (`numpy.ndarray`_, list):
            A list of 2D array containing the wavelength of each pixel. The shape of
            each 2D array is (nspec, nspat).
        slitid_img_gpm (`numpy.ndarray`_, list):
            A list of 2D array containing the slitmask of each pixel. The shape of
            each 2D array is (nspec, nspat).
            A zero value indicates that a pixel is either not on a slit or it is a bad pixel.
            All other values are the slit spatial ID number.
        wghtImg (`numpy.ndarray`_, list):
            A list of 2D array containing the weights of each pixel to be used in the
            combination. The shape of each 2D array is (nspec, nspat).
        all_wcs (`astropy.wcs.WCS`_, list):
            A list of `astropy.wcs.WCS`_ objects, one for each spec2d file
        tilts (list):
            A list of `numpy.ndarray`_ objects, one for each spec2d file,
            containing the tilts of each pixel. The shape of each 2D array is
            (nspec, nspat).
        slits (:class:`pypeit.slittrace.SlitTraceSet`, list):
            A list of :class:`pypeit.slittrace.SlitTraceSet` objects, one for each
            spec2d file, containing the properties of the slit for each spec2d file
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        all_dar (:class:`~pypeit.coadd3d.DARcorrection`, list):
            A Class containing the DAR correction information, or a list of DARcorrection
            classes. If a list, it must be the same length as astrom_trans.
        ra_offset (float, list):
            A float or list of floats containing the RA offset of each spec2d file
            relative to the first spec2d file
        dec_offset (float, list):
            A float or list of floats containing the DEC offset of each spec2d file
            relative to the first spec2d file
        spec_subpixel (int, optional):
            What is the subpixellation factor in the spectral direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral
            direction.
        spat_subpixel (int, optional):
            What is the subpixellation factor in the spatial direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial
            direction.
        slice_subpixel (int, optional):
            What is the subpixellation factor in the slice direction. Higher
            values give more reliable results, but note that the time required
            goes as (``slice_subpixel``). The default value is 5, which divides
            each IFU slice into 5 subslices in the slice direction.
        skip_subpix_weights (bool, optional):
            If True, the computationally expensive step to calculate the
            subpixellation weights will be skipped. If set the True, note that
            the variance cubes returned will not be accurate. However, if you
            are not interested in the variance cubes, this can save a lot of
            time, and this is an example where you might consider setting this
            variable to True. The flux datacube is unaffected by this variable.
            The default is False.
        correct_dar (bool, optional):
            If True, the DAR correction will be applied to the datacube. The
            default is True.

    Returns:
        :obj:`tuple`: Three or four `numpy.ndarray`_ objects containing (1) the
        datacube generated from the subpixellated inputs, (2) the corresponding
        variance cube, (3) the corresponding bad pixel mask cube, and (4) the
        residual cube.  The latter is only returned if debug is True.
    """
    # Check the inputs for combinations of lists or not
    _sciImg, _ivarImg, _waveImg, _gpmImg, _wghtImg, _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset = \
        check_inputs([sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset])
    numframes = len(_sciImg)

    # Prepare the output arrays
    outshape = (bins[0].size-1, bins[1].size-1, bins[2].size-1)
    binrng = np.array([[bins[0][0], bins[0][-1]], [bins[1][0], bins[1][-1]], [bins[2][0], bins[2][-1]]])
    flxcube, varcube, normcube = np.zeros(outshape), np.zeros(outshape), np.zeros(outshape)
    # Divide each pixel into subpixels
    spec_offs = np.arange(0.5/spec_subpixel, 1, 1/spec_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    spat_offs = np.arange(0.5/spat_subpixel, 1, 1/spat_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    slice_offs = np.arange(0.5/slice_subpixel, 1, 1/slice_subpixel) - 0.5  # -0.5 is to offset from the centre of each slice.
    spat_x, spec_y = np.meshgrid(spat_offs, spec_offs)
    num_subpixels = spec_subpixel * spat_subpixel  # Number of subpixels (spat & spec) per detector pixel
    num_all_subpixels = num_subpixels * slice_subpixel  # Number of subpixels, including slice subpixels
    # Loop through all exposures
    for fr in range(numframes):
        onslit_gpm = _gpmImg[fr]
        this_onslit_gpm = onslit_gpm > 0
        this_specpos, this_spatpos = np.where(this_onslit_gpm)
        this_spatid = onslit_gpm[this_onslit_gpm]

        # Extract tilts and slits for convenience
        this_tilts = _tilts[fr]
        this_slits = _slits[fr]
        this_wcs = _all_wcs[fr]
        this_astrom_trans = _astrom_trans[fr]
        this_wght_subpix = _wghtImg[fr][this_onslit_gpm]
        this_sci = _sciImg[fr][this_onslit_gpm]
        this_var = utils.inverse(_ivarImg[fr][this_onslit_gpm])
        this_wav = _waveImg[fr][this_onslit_gpm]
        # Loop through all slits
        for sl, spatid in enumerate(this_slits.spat_id):
            if numframes == 1:
                msgs.info(f"Resampling slit {sl + 1}/{this_slits.nslits}")
            else:
                msgs.info(f"Resampling slit {sl + 1}/{this_slits.nslits} of frame {fr + 1}/{numframes}")
            # Find the pixels on this slit
            this_sl = np.where(this_spatid == spatid)
            wpix = (this_specpos[this_sl], this_spatpos[this_sl])
            # Create an array to index each subpixel
            numpix = wpix[0].size
            # Generate a spline between spectral pixel position and wavelength
            yspl = this_tilts[wpix] * (this_slits.nspec - 1)
            tiltpos = np.add.outer(yspl, spec_y).flatten()
            wspl = this_wav[this_sl]
            asrt = np.argsort(yspl)
            wave_spl = interp1d(yspl[asrt], wspl[asrt], kind='linear', bounds_error=False, fill_value='extrapolate')
            # Calculate the wavelength at each subpixel
            this_wave_subpix = wave_spl(tiltpos)
            # Calculate the DAR correction at each sub pixel
            ra_corr, dec_corr = 0.0, 0.0
            if correct_dar:
                # NOTE :: This routine needs the wavelengths to be expressed in Angstroms
                ra_corr, dec_corr = _all_dar[fr].correction( this_wave_subpix)
            # Calculate spatial and spectral positions of the subpixels
            spat_xx = np.add.outer(wpix[1], spat_x.flatten()).flatten()
            spec_yy = np.add.outer(wpix[0], spec_y.flatten()).flatten()
            # Transform this to spatial location
            spatpos_subpix = _astrom_trans[fr].transform(sl, spat_xx, spec_yy)
            spatpos = _astrom_trans[fr].transform(sl, wpix[1], wpix[0])
            ssrt = np.argsort(spatpos)
            # Initialize the voxel coordinates for each spec2D pixel
            vox_coord = np.full((numpix, num_all_subpixels, 3), -1, dtype=float)
            # Loop over the subslices
            for ss in range(slice_subpixel):
                if slice_subpixel > 1:
                    # Only print this if there are multiple subslices
                    msgs.info(f"Resampling subslice {ss+1}/{slice_subpixel}")
                # Generate an RA/Dec image for this subslice
                raimg, decimg, minmax = this_slits.get_radec_image(this_wcs, this_astrom_trans, this_tilts,
                                                                   slit_compute=sl, slice_offset=slice_offs[ss])
                this_ra = raimg[this_onslit_gpm]
                this_dec = decimg[this_onslit_gpm]
                # Interpolate the RA/Dec over the subpixel spatial positions
                tmp_ra = this_ra[this_sl]
                tmp_dec = this_dec[this_sl]
                ra_spl = interp1d(spatpos[ssrt], tmp_ra[ssrt], kind='linear', bounds_error=False, fill_value='extrapolate')
                dec_spl = interp1d(spatpos[ssrt], tmp_dec[ssrt], kind='linear', bounds_error=False, fill_value='extrapolate')
                # Evaluate the RA/Dec at the subpixel spatial positions
                this_ra_int = ra_spl(spatpos_subpix)
                this_dec_int = dec_spl(spatpos_subpix)
                # Now apply the DAR correction and any user-supplied offsets
                this_ra_int += ra_corr + _ra_offset[fr]
                this_dec_int += dec_corr + _dec_offset[fr]
                # Convert world coordinates to voxel coordinates, then histogram
                sslo = ss * num_subpixels
                sshi = (ss + 1) * num_subpixels
                vox_coord[:,sslo:sshi,:] = output_wcs.wcs_world2pix(np.vstack((this_ra_int, this_dec_int, this_wave_subpix * 1.0E-10)).T, 0).reshape(numpix, num_subpixels, 3)
            # Convert the voxel coordinates to a bin index
            if num_all_subpixels == 1 or skip_subpix_weights:
                subpix_wght = 1.0
            else:
                msgs.info("Preparing subpixel weights")
                vox_index = np.floor(outshape * (vox_coord - binrng[:,0].reshape((1, 1, 3))) /
                                                (binrng[:,1] - binrng[:,0]).reshape((1, 1, 3))).astype(int)
                # Convert to a unique index
                vox_index = np.dot(vox_index, np.array([1, outshape[0], outshape[0]*outshape[1]]))
                # Calculate the number of repeated indices for each subpixel - this is the subpixel weights
                subpix_wght = np.apply_along_axis(utils.occurrences, 1, vox_index).flatten()
            # Reshape the voxel coordinates
            vox_coord = vox_coord.reshape(numpix * num_all_subpixels, 3)
            # Use the "fast histogram" algorithm, that assumes regular bin spacing
            flxcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_sci[this_sl] * this_wght_subpix[this_sl], num_all_subpixels) * subpix_wght)
            varcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_var[this_sl] * this_wght_subpix[this_sl]**2, num_all_subpixels) * subpix_wght**3)
            normcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(this_wght_subpix[this_sl], num_all_subpixels) * subpix_wght)

    # Normalise the datacube and variance cube
    nc_inverse = utils.inverse(normcube)
    flxcube *= nc_inverse
    varcube *= nc_inverse**2
    bpmcube = (normcube == 0).astype(np.uint8)
<<<<<<< HEAD
=======
    if debug:
        residcube *= nc_inverse
        return datacube, varcube, bpmcube, residcube
    return datacube, varcube, bpmcube


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
        if par_outfile == "":
            par_outfile = "datacube.fits"
        # Check the output files don't exist
        outfile = par_outfile if ".fits" in par_outfile else par_outfile + ".fits"
    else:
        if par_outfile == "":
            outfile = fil.replace("spec2d_", "spec3d_")
        else:
            # Use the output filename as a prefix
            outfile = os.path.splitext(par_outfile)[0] + "_{0:03d}.fits".format(idx)
    # Return the outfile
    return outfile


def get_output_whitelight_filename(outfile):
    """
    Given the output filename of a datacube, create an appropriate whitelight
    fits file name

    Args:
        outfile (str):
            The output filename used for the datacube.

    Returns:
        str: The output filename to use for the whitelight image.
    """
    out_wl_filename = os.path.splitext(outfile)[0] + "_whitelight.fits"
    return out_wl_filename


def coadd_cube(files, opts, spectrograph=None, parset=None, overwrite=False):
    """ Main routine to coadd spec2D files into a 3D datacube

    Args:
        files (:obj:`list`):
            List of all spec2D files
        opts (:obj:`dict`):
            coadd2d options associated with each spec2d file
        spectrograph (:obj:`str`, :class:`~pypeit.spectrographs.spectrograph.Spectrograph`, optional):
            The name or instance of the spectrograph used to obtain the data.
            If None, this is pulled from the file header.
        parset (:class:`~pypeit.par.pypeitpar.PypeItPar`, optional):
            An instance of the parameter set.  If None, assumes that detector 1
            is the one reduced and uses the default reduction parameters for the
            spectrograph (see
            :func:`~pypeit.spectrographs.spectrograph.Spectrograph.default_pypeit_par`
            for the relevant spectrograph class).
        overwrite (:obj:`bool`, optional):
            Overwrite the output file, if it exists?
    """
    if spectrograph is None:
        with fits.open(files[0]) as hdu:
            spectrograph = hdu[0].header['PYP_SPEC']

    if isinstance(spectrograph, str):
        spec = load_spectrograph(spectrograph)
        specname = spectrograph
    else:
        # Assume it's a Spectrograph instance
        spec = spectrograph
        specname = spectrograph.name

    # Grab the parset, if not provided
    if parset is None:
        # TODO :: Use config_specific_par instead?
        parset = spec.default_pypeit_par()
    cubepar = parset['reduce']['cube']
    flatpar = parset['calibrations']['flatfield']
    senspar = parset['sensfunc']

    # prep
    numfiles = len(files)
    method = cubepar['method'].lower()
    combine = cubepar['combine']
    align = cubepar['align']
    # If there is only one frame being "combined" AND there's no reference image, then don't compute the translation.
    if numfiles == 1 and cubepar["reference_image"] is None:
        if not align:
            msgs.warn("Parameter 'align' should be False when there is only one frame and no reference image")
            msgs.info("Setting 'align' to False")
        align = False
    if opts['ra_offset'] is not None:
        if not align:
            msgs.warn("When 'ra_offset' and 'dec_offset' are set, 'align' must be True.")
            msgs.info("Setting 'align' to True")
        align = True
    # TODO :: The default behaviour (combine=False, align=False) produces a datacube that uses the instrument WCS
    #  It should be possible (and perhaps desirable) to do a spatial alignment (i.e. align=True), apply this to the
    #  RA,Dec values of each pixel, and then use the instrument WCS to save the output (or, just adjust the crval).
    #  At the moment, if the user wishes to spatially align the frames, a different WCS is generated.
    if histogramdd is None:
        msgs.warn("Generating a datacube is faster if you install fast-histogram:"+msgs.newline()+
                  "https://pypi.org/project/fast-histogram/")
        if method != 'ngp':
            msgs.warn("Forcing NGP algorithm, because fast-histogram is not installed")
            method = 'ngp'

    # Determine what method is requested
    spec_subpixel, spat_subpixel = 1, 1
    if method == "subpixel":
        msgs.info("Adopting the subpixel algorithm to generate the datacube.")
        spec_subpixel, spat_subpixel = cubepar['spec_subpixel'], cubepar['spat_subpixel']
    elif method == "ngp":
        msgs.info("Adopting the nearest grid point (NGP) algorithm to generate the datacube.")
    else:
        msgs.error(f"The following datacube method is not allowed: {method}")

    # Get the detector number and string representation
    det = 1 if parset['rdx']['detnum'] is None else parset['rdx']['detnum']
    detname = spec.get_det_name(det)

    # Check if the output file exists
    if combine:
        outfile = get_output_filename("", cubepar['output_filename'], combine)
        out_whitelight = get_output_whitelight_filename(outfile)
        if os.path.exists(outfile) and not overwrite:
            msgs.error("Output filename already exists:"+msgs.newline()+outfile)
        if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
            msgs.error("Output filename already exists:"+msgs.newline()+out_whitelight)
    else:
        # Finally, if there's just one file, check if the output filename is given
        if numfiles == 1 and cubepar['output_filename'] != "":
            outfile = get_output_filename("", cubepar['output_filename'], True, -1)
            out_whitelight = get_output_whitelight_filename(outfile)
            if os.path.exists(outfile) and not overwrite:
                msgs.error("Output filename already exists:" + msgs.newline() + outfile)
            if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
                msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)
        else:
            for ff in range(numfiles):
                outfile = get_output_filename(files[ff], cubepar['output_filename'], combine, ff+1)
                out_whitelight = get_output_whitelight_filename(outfile)
                if os.path.exists(outfile) and not overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + outfile)
                if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)

    # Check the reference cube and image exist, if requested
    fluxcal = False
    blaze_wave, blaze_spec = None, None
    blaze_spline, flux_spline = None, None
    if cubepar['standard_cube'] is not None:
        fluxcal = True
        ss_file = cubepar['standard_cube']
        if not os.path.exists(ss_file):
            msgs.error("Standard cube does not exist:" + msgs.newline() + ss_file)
        msgs.info(f"Loading standard star cube: {ss_file:s}")
        # Load the standard star cube and retrieve its RA + DEC
        stdcube = fits.open(ss_file)
        star_ra, star_dec = stdcube[1].header['CRVAL1'], stdcube[1].header['CRVAL2']

        # Extract a spectrum of the standard star
        wave, Nlam_star, Nlam_ivar_star, gpm_star = extract_standard_spec(stdcube)

        # Extract the information about the blaze
        if cubepar['grating_corr']:
            blaze_wave_curr, blaze_spec_curr = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
            blaze_spline_curr = interp1d(blaze_wave_curr, blaze_spec_curr,
                                         kind='linear', bounds_error=False, fill_value="extrapolate")
            # The first standard star cube is used as the reference blaze spline
            if blaze_spline is None:
                blaze_wave, blaze_spec = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
                blaze_spline = interp1d(blaze_wave, blaze_spec,
                                        kind='linear', bounds_error=False, fill_value="extrapolate")
            # Perform a grating correction
            grat_corr = correct_grating_shift(wave.value, blaze_wave_curr, blaze_spline_curr, blaze_wave, blaze_spline)
            # Apply the grating correction to the standard star spectrum
            Nlam_star /= grat_corr
            Nlam_ivar_star *= grat_corr**2

        # Read in some information above the standard star
        std_dict = flux_calib.get_standard_spectrum(star_type=senspar['star_type'],
                                         star_mag=senspar['star_mag'],
                                         ra=star_ra, dec=star_dec)
        # Calculate the sensitivity curve
        # TODO :: This needs to be addressed... unify flux calibration into the main PypeIt routines.
        msgs.warn("Datacubes are currently flux-calibrated using the UVIS algorithm... this will be deprecated soon")
        zeropoint_data, zeropoint_data_gpm, zeropoint_fit, zeropoint_fit_gpm =\
            flux_calib.fit_zeropoint(wave.value, Nlam_star, Nlam_ivar_star, gpm_star, std_dict,
                          mask_hydrogen_lines=senspar['mask_hydrogen_lines'],
                          mask_helium_lines=senspar['mask_helium_lines'],
                          hydrogen_mask_wid=senspar['hydrogen_mask_wid'],
                          nresln=senspar['UVIS']['nresln'], resolution=senspar['UVIS']['resolution'],
                          trans_thresh=senspar['UVIS']['trans_thresh'], polyorder=senspar['polyorder'],
                          polycorrect=senspar['UVIS']['polycorrect'], polyfunc=senspar['UVIS']['polyfunc'])
        wgd = np.where(zeropoint_fit_gpm)
        sens = np.power(10.0, -0.4 * (zeropoint_fit[wgd] - flux_calib.ZP_UNIT_CONST)) / np.square(wave[wgd])
        flux_spline = interp1d(wave[wgd], sens, kind='linear', bounds_error=False, fill_value="extrapolate")

    # If a reference image has been set, check that it exists
    if cubepar['reference_image'] is not None:
        if not os.path.exists(cubepar['reference_image']):
            msgs.error("Reference image does not exist:" + msgs.newline() + cubepar['reference_image'])

    # Initialise arrays for storage
    ifu_ra, ifu_dec = np.array([]), np.array([])  # The RA and Dec at the centre of the IFU, as stored in the header
    all_ra, all_dec, all_wave = np.array([]), np.array([]), np.array([])
    all_sci, all_ivar, all_idx, all_wghts = np.array([]), np.array([]), np.array([]), np.array([])
    all_spatpos, all_specpos, all_spatid = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=int)
    all_tilts, all_slits, all_align = [], [], []
    all_wcs = []
    dspat = None if cubepar['spatial_delta'] is None else cubepar['spatial_delta']/3600.0  # binning size on the sky (/3600 to convert to degrees)
    dwv = cubepar['wave_delta']       # binning size in wavelength direction (in Angstroms)
    wave_ref = None
    mnmx_wv = None  # Will be used to store the minimum and maximum wavelengths of every slit and frame.
    weights = np.ones(numfiles)  # Weights to use when combining cubes
    flat_splines = dict()   # A dictionary containing the splines of the flatfield
    # Load the default scaleimg frame for the scale correction
    scalecorr_default = "none"
    relScaleImgDef = np.array([1])
    if cubepar['scale_corr'] is not None:
        if cubepar['scale_corr'] == "image":
            msgs.info("The default relative spectral illumination correction will use the science image")
            scalecorr_default = "image"
        else:
            msgs.info("Loading default scale image for relative spectral illumination correction:" +
                      msgs.newline() + cubepar['scale_corr'])
            try:
                spec2DObj = spec2dobj.Spec2DObj.from_file(cubepar['scale_corr'], detname)
                relScaleImgDef = spec2DObj.scaleimg
                scalecorr_default = cubepar['scale_corr']
            except:
                msgs.warn("Could not load scaleimg from spec2d file:" + msgs.newline() +
                          cubepar['scale_corr'] + msgs.newline() +
                          "scale correction will not be performed unless you have specified the correct" + msgs.newline() +
                          "scale_corr file in the spec2d block")
                cubepar['scale_corr'] = None
                scalecorr_default = "none"

    # Load the default sky frame to be used for sky subtraction
    skysub_default = "image"
    skyImgDef, skySclDef = None, None  # This is the default behaviour (i.e. to use the "image" for the sky subtraction)
    if cubepar['skysub_frame'] in [None, 'none', '', 'None']:
        skysub_default = "none"
        skyImgDef = np.array([0.0])  # Do not perform sky subtraction
        skySclDef = np.array([0.0])  # Do not perform sky subtraction
    elif cubepar['skysub_frame'].lower() == "image":
        msgs.info("The sky model in the spec2d science frames will be used for sky subtraction" +msgs.newline() +
                  "(unless specific skysub frames have been specified)")
        skysub_default = "image"
    else:
        msgs.info("Loading default image for sky subtraction:" +
                  msgs.newline() + cubepar['skysub_frame'])
        try:
            spec2DObj = spec2dobj.Spec2DObj.from_file(cubepar['skysub_frame'], detname)
            skysub_exptime = fits.open(cubepar['skysub_frame'])[0].header['EXPTIME']
        except:
            msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + cubepar['skysub_frame'])
        skysub_default = cubepar['skysub_frame']
        skyImgDef = spec2DObj.skymodel/skysub_exptime  # Sky counts/second
        skySclDef = spec2DObj.scaleimg

    # Load all spec2d files and prepare the data for making a datacube
    for ff, fil in enumerate(files):
        # Load it up
        msgs.info("Loading PypeIt spec2d frame:" + msgs.newline() + fil)
        spec2DObj = spec2dobj.Spec2DObj.from_file(fil, detname)
        detector = spec2DObj.detector
        spat_flexure = None  #spec2DObj.sci_spat_flexure

        # Load the header
        hdr = spec2DObj.head0
        ifu_ra = np.append(ifu_ra, spec.compound_meta([hdr], 'ra'))
        ifu_dec = np.append(ifu_dec, spec.compound_meta([hdr], 'dec'))

        # Get the exposure time
        exptime = hdr['EXPTIME']

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        # TODO :: Consider loading all calibrations into a single variable.

        # Initialise the slit edges
        msgs.info("Constructing slit image")
        slits = spec2DObj.slits
        slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=spat_flexure)
        slits_left, slits_right, _ = slits.select_edges(initial=True, flexure=spat_flexure)

        # The order of operations below proceeds as follows:
        #  (1) Get science image
        #  (2) Subtract sky (note, if a joint fit has been performed, the relative scale correction is applied in the reduction!)
        #  (3) Apply relative scale correction to both science and ivar

        # Set the default behaviour if a global skysub frame has been specified
        this_skysub = skysub_default
        if skysub_default == "image":
            skyImg = spec2DObj.skymodel
            skyScl = spec2DObj.scaleimg
        else:
            skyImg = skyImgDef.copy() * exptime
            skyScl = skySclDef.copy()
        # See if there's any changes from the default behaviour
        if opts['skysub_frame'][ff] is not None:
            if opts['skysub_frame'][ff].lower() == 'default':
                if skysub_default == "image":
                    skyImg = spec2DObj.skymodel
                    skyScl = spec2DObj.scaleimg
                    this_skysub = "image"  # Use the current spec2d for sky subtraction
                else:
                    skyImg = skyImgDef.copy() * exptime
                    skyScl = skySclDef.copy() * exptime
                    this_skysub = skysub_default  # Use the global value for sky subtraction
            elif opts['skysub_frame'][ff].lower() == 'image':
                skyImg = spec2DObj.skymodel
                skyScl = spec2DObj.scaleimg
                this_skysub = "image"  # Use the current spec2d for sky subtraction
            elif opts['skysub_frame'][ff].lower() == 'none':
                skyImg = np.array([0.0])
                skyScl = np.array([1.0])
                this_skysub = "none"  # Don't do sky subtraction
            else:
                # Load a user specified frame for sky subtraction
                msgs.info("Loading skysub frame:" + msgs.newline() + opts['skysub_frame'][ff])
                try:
                    spec2DObj_sky = spec2dobj.Spec2DObj.from_file(opts['skysub_frame'][ff], detname)
                    skysub_exptime = fits.open(opts['skysub_frame'][ff])[0].header['EXPTIME']
                except:
                    msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + opts['skysub_frame'][ff])
                # TODO :: Consider allowing the actual frame (instead of the skymodel) to be used as the skysub image - make sure the BPM is carried over.
                #      :: Allow sky data fitting (i.e. scale the flux of a skysub frame to the science frame data)
                skyImg = spec2DObj_sky.skymodel * exptime / skysub_exptime  # Sky counts
                skyScl = spec2DObj_sky.scaleimg
                this_skysub = opts['skysub_frame'][ff]  # User specified spec2d for sky subtraction
        if this_skysub == "none":
            msgs.info("Sky subtraction will not be performed.")
        else:
            msgs.info("Using the following frame for sky subtraction:"+msgs.newline()+this_skysub)

        # Load the relative scale image, if something other than the default has been provided
        this_scalecorr = scalecorr_default
        relScaleImg = relScaleImgDef.copy()
        if opts['scale_corr'][ff] is not None:
            if opts['scale_corr'][ff].lower() == 'default':
                if scalecorr_default == "image":
                    relScaleImg = spec2DObj.scaleimg
                    this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
                else:
                    this_scalecorr = scalecorr_default  # Use the default value for the scale correction
            elif opts['scale_corr'][ff].lower() == 'image':
                relScaleImg = spec2DObj.scaleimg
                this_scalecorr = "image"  # Use the current spec2d for the relative spectral illumination scaling
            elif opts['scale_corr'][ff].lower() == 'none':
                relScaleImg = np.array([1])
                this_scalecorr = "none"  # Don't do relative spectral illumination scaling
            else:
                # Load a user specified frame for sky subtraction
                msgs.info("Loading the following frame for the relative spectral illumination correction:" +
                          msgs.newline() + opts['scale_corr'][ff])
                try:
                    spec2DObj_scl = spec2dobj.Spec2DObj.from_file(opts['scale_corr'][ff], detname)
                except:
                    msgs.error("Could not load skysub image from spec2d file:" + msgs.newline() + opts['skysub_frame'][ff])
                relScaleImg = spec2DObj_scl.scaleimg
                this_scalecorr = opts['scale_corr'][ff]
        if this_scalecorr == "none":
            msgs.info("Relative spectral illumination correction will not be performed.")
        else:
            msgs.info("Using the following frame for the relative spectral illumination correction:" +
                      msgs.newline()+this_scalecorr)

        # Prepare the relative scaling factors
        relSclSky = skyScl/spec2DObj.scaleimg  # This factor ensures the sky has the same relative scaling as the science frame
        relScale = spec2DObj.scaleimg/relScaleImg  # This factor is applied to the sky subtracted science frame

        # Extract the relevant information from the spec2d file
        sciImg = (spec2DObj.sciimg - skyImg*relSclSky)*relScale  # Subtract sky and apply relative illumination
        ivar = spec2DObj.ivarraw / relScale**2
        waveimg = spec2DObj.waveimg
        bpmmask = spec2DObj.bpmmask

        # TODO :: Really need to write some detailed information in the docs about all of the various corrections that can optionally be applied

        # TODO :: Include a flexure correction from the sky frame? Note, you cannot use the waveimg from a sky frame,
        #  since the heliocentric correction may have been applied to the sky frame. Need to recalculate waveimg using
        #  the slitshifts from a skyimage, and then apply the vel_corr from the science image.

        wnonzero = (waveimg != 0.0)
        if not np.any(wnonzero):
            msgs.error("The wavelength image contains only zeros - You need to check the data reduction.")
        wave0 = waveimg[wnonzero].min()
        # Calculate the delta wave in every pixel on the slit
        waveimp = np.roll(waveimg, 1, axis=0)
        waveimn = np.roll(waveimg, -1, axis=0)
        dwaveimg = np.zeros_like(waveimg)
        # All good pixels
        wnz = np.where((waveimg!=0) & (waveimp!=0))
        dwaveimg[wnz] = np.abs(waveimg[wnz]-waveimp[wnz])
        # All bad pixels
        wnz = np.where((waveimg!=0) & (waveimp==0))
        dwaveimg[wnz] = np.abs(waveimg[wnz]-waveimn[wnz])
        # All endpoint pixels
        dwaveimg[0, :] = np.abs(waveimg[0, :] - waveimn[0, :])
        dwaveimg[-1, :] = np.abs(waveimg[-1, :] - waveimp[-1, :])
        dwv = np.median(dwaveimg[dwaveimg != 0.0]) if cubepar['wave_delta'] is None else cubepar['wave_delta']

        msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

        # Obtain the minimum and maximum wavelength of all slits
        if mnmx_wv is None:
            mnmx_wv = np.zeros((len(files), slits.nslits, 2))
        for slit_idx, slit_spat in enumerate(slits.spat_id):
            onslit_init = (slitid_img_init == slit_spat)
            mnmx_wv[ff, slit_idx, 0] = np.min(waveimg[onslit_init])
            mnmx_wv[ff, slit_idx, 1] = np.max(waveimg[onslit_init])

        # Remove edges of the spectrum where the sky model is bad
        sky_is_good = make_good_skymask(slitid_img_init, spec2DObj.tilts)

        # Construct a good pixel mask
        # TODO: This should use the mask function to figure out which elements are masked.
        onslit_gpm = (slitid_img_init > 0) & (bpmmask.mask == 0) & sky_is_good

        # Grab the WCS of this frame
        frame_wcs = spec.get_wcs(spec2DObj.head0, slits, detector.platescale, wave0, dwv)
        all_wcs.append(copy.deepcopy(frame_wcs))

        # Find the largest spatial scale of all images being combined
        # TODO :: probably need to put this in the DetectorContainer
        pxscl = detector.platescale * parse.parse_binning(detector.binning)[1] / 3600.0  # This should be degrees/pixel
        slscl = spec.get_meta_value([spec2DObj.head0], 'slitwid')
        if dspat is None:
            dspat = max(pxscl, slscl)
        if pxscl > dspat:
            msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the pixel scale ({1:f} arcsec)".format(3600.0*dspat, 3600.0*pxscl))
        if slscl > dspat:
            msgs.warn("Spatial scale requested ({0:f} arcsec) is less than the slicer scale ({1:f} arcsec)".format(3600.0*dspat, 3600.0*slscl))

        # Loading the alignments frame for these data
        alignments = None
        if cubepar['astrometric']:
            key = alignframe.Alignments.calib_type.upper()
            if key in spec2DObj.calibs:
                alignfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
                if os.path.exists(alignfile) and cubepar['astrometric']:
                    msgs.info("Loading alignments")
                    alignments = alignframe.Alignments.from_file(alignfile)
            else:
                msgs.warn(f'Processed alignment frame not recorded or not found!')
                msgs.info("Using slit edges for astrometric transform")
        else:
            msgs.info("Using slit edges for astrometric transform")
        # If nothing better was provided, use the slit edges
        if alignments is None:
            left, right, _ = slits.select_edges(initial=True, flexure=spat_flexure)
            locations = [0.0, 1.0]
            traces = np.append(left[:,None,:], right[:,None,:], axis=1)
        else:
            locations = parset['calibrations']['alignment']['locations']
            traces = alignments.traces
        # Generate an RA/DEC image
        msgs.info("Generating RA/DEC image")
        alignSplines = alignframe.AlignmentSplines(traces, locations, spec2DObj.tilts)
        raimg, decimg, minmax = slits.get_radec_image(frame_wcs, alignSplines, spec2DObj.tilts,
                                                      initial=True, flexure=spat_flexure)
        # Perform the DAR correction
        if wave_ref is None:
            wave_ref = 0.5*(np.min(waveimg[onslit_gpm]) + np.max(waveimg[onslit_gpm]))
        # Get DAR parameters
        raval = spec.get_meta_value([spec2DObj.head0], 'ra')
        decval = spec.get_meta_value([spec2DObj.head0], 'dec')
        obstime = spec.get_meta_value([spec2DObj.head0], 'obstime')
        pressure = spec.get_meta_value([spec2DObj.head0], 'pressure')
        temperature = spec.get_meta_value([spec2DObj.head0], 'temperature')
        rel_humidity = spec.get_meta_value([spec2DObj.head0], 'humidity')
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))
        location = spec.location  # TODO :: spec.location should probably end up in the TelescopePar (spec.telescope.location)
        if pressure == 0.0:
            msgs.warn("Pressure is set to zero - DAR correction will not be performed")
        else:
            msgs.info("DAR correction parameters:"+msgs.newline() +
                      "   Pressure = {0:f} bar".format(pressure) + msgs.newline() +
                      "   Temperature = {0:f} deg C".format(temperature) + msgs.newline() +
                      "   Humidity = {0:f}".format(rel_humidity))
            ra_corr, dec_corr = correct_dar(waveimg[onslit_gpm], coord, obstime, location,
                                            pressure * units.bar, temperature * units.deg_C, rel_humidity, wave_ref=wave_ref)
            raimg[onslit_gpm] += ra_corr*np.cos(np.mean(decimg[onslit_gpm]) * np.pi / 180.0)
            decimg[onslit_gpm] += dec_corr

        # Get copies of arrays to be saved
        wave_ext = waveimg[onslit_gpm].copy()
        flux_ext = sciImg[onslit_gpm].copy()
        ivar_ext = ivar[onslit_gpm].copy()
        dwav_ext = dwaveimg[onslit_gpm].copy()

        # Correct for sensitivity as a function of grating angle
        # (this assumes the spectrum of the flatfield lamp has the same shape for all setups)
        key = flatfield.FlatImages.calib_type.upper()
        if key not in spec2DObj.calibs:
            msgs.error('Processed flat calibration file not recorded by spec2d file!')
        flatfile = os.path.join(spec2DObj.calibs['DIR'], spec2DObj.calibs[key])
        if cubepar['grating_corr'] and flatfile not in flat_splines.keys():
            msgs.info("Calculating relative sensitivity for grating correction")
            # Check if the Flat file exists
            if not os.path.exists(flatfile):
                msgs.error("Grating correction requested, but the following file does not exist:" +
                           msgs.newline() + flatfile)
            # Load the Flat file
            flatimages = flatfield.FlatImages.from_file(flatfile)
            total_illum = flatimages.fit2illumflat(slits, finecorr=False, frametype='illum', initial=True, spat_flexure=spat_flexure) * \
                          flatimages.fit2illumflat(slits, finecorr=True, frametype='illum', initial=True, spat_flexure=spat_flexure)
            flatframe = flatimages.pixelflat_raw / total_illum
            if flatimages.pixelflat_spec_illum is None:
                # Calculate the relative scale
                scale_model = flatfield.illum_profile_spectral(flatframe, waveimg, slits,
                                                               slit_illum_ref_idx=flatpar['slit_illum_ref_idx'], model=None,
                                                               skymask=None, trim=flatpar['slit_trim'], flexure=spat_flexure,
                                                               smooth_npix=flatpar['slit_illum_smooth_npix'])
            else:
                msgs.info("Using relative spectral illumination from FlatImages")
                scale_model = flatimages.pixelflat_spec_illum
            # Apply the relative scale and generate a 1D "spectrum"
            onslit = waveimg != 0
            wavebins = np.linspace(np.min(waveimg[onslit]), np.max(waveimg[onslit]), slits.nspec)
            hist, edge = np.histogram(waveimg[onslit], bins=wavebins, weights=flatframe[onslit]/scale_model[onslit])
            cntr, edge = np.histogram(waveimg[onslit], bins=wavebins)
            cntr = cntr.astype(float)
            norm = (cntr != 0) / (cntr + (cntr == 0))
            spec_spl = hist * norm
            wave_spl = 0.5 * (wavebins[1:] + wavebins[:-1])
            flat_splines[flatfile] = interp1d(wave_spl, spec_spl, kind='linear',
                                              bounds_error=False, fill_value="extrapolate")
            flat_splines[flatfile+"_wave"] = wave_spl.copy()
            # Check if a reference blaze spline exists (either from a standard star if fluxing or from a previous
            # exposure in this for loop)
            if blaze_spline is None:
                blaze_wave, blaze_spec = wave_spl, spec_spl
                blaze_spline = interp1d(wave_spl, spec_spl, kind='linear',
                                        bounds_error=False, fill_value="extrapolate")

        # Perform extinction correction
        msgs.info("Applying extinction correction")
        longitude = spec.telescope['longitude']
        latitude = spec.telescope['latitude']
        airmass = spec2DObj.head0[spec.meta['airmass']['card']]
        extinct = flux_calib.load_extinction_data(longitude, latitude, senspar['UVIS']['extinct_file'])
        # extinction_correction requires the wavelength is sorted
        wvsrt = np.argsort(wave_ext)
        ext_corr = flux_calib.extinction_correction(wave_ext[wvsrt] * units.AA, airmass, extinct)
        # Grating correction
        grat_corr = 1.0
        if cubepar['grating_corr']:
            grat_corr = correct_grating_shift(wave_ext[wvsrt], flat_splines[flatfile + "_wave"], flat_splines[flatfile],
                                              blaze_wave, blaze_spline)
        # Sensitivity function
        sens_func = 1.0
        if fluxcal:
            msgs.info("Calculating the sensitivity function")
            sens_func = flux_spline(wave_ext[wvsrt])
        # Convert the flux_sav to counts/s,  correct for the relative sensitivity of different setups
        ext_corr *= sens_func / (exptime * grat_corr)
        # Correct for extinction
        flux_sav = flux_ext[wvsrt] * ext_corr
        ivar_sav = ivar_ext[wvsrt] / ext_corr ** 2

        # Convert units to Counts/s/Ang/arcsec2
        # Slicer sampling * spatial pixel sampling
        sl_deg = np.sqrt(frame_wcs.wcs.cd[0, 0] ** 2 + frame_wcs.wcs.cd[1, 0] ** 2)
        px_deg = np.sqrt(frame_wcs.wcs.cd[1, 1] ** 2 + frame_wcs.wcs.cd[0, 1] ** 2)
        scl_units = dwav_ext[wvsrt] * (3600.0 * sl_deg) * (3600.0 * px_deg)
        flux_sav /= scl_units
        ivar_sav *= scl_units ** 2

        # sort back to the original ordering
        resrt = np.argsort(wvsrt)
        numpix = raimg[onslit_gpm].size

        # Calculate the weights relative to the zeroth cube
        weights[ff] = 1.0#exptime  #np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2

        # Get the slit image and then unset pixels in the slit image that are bad
        this_specpos, this_spatpos = np.where(onslit_gpm)
        this_spatid = slitid_img_init[onslit_gpm]

        # If individual frames are to be output without aligning them,
        # there's no need to store information, just make the cubes now
        if not combine and not align:
            # Get the output filename
            if numfiles == 1 and cubepar['output_filename'] != "":
                outfile = get_output_filename("", cubepar['output_filename'], True, -1)
            else:
                outfile = get_output_filename(fil, cubepar['output_filename'], combine, ff+1)
            # Get the coordinate bounds
            slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
            numwav = int((np.max(waveimg) - wave0) / dwv)
            bins = spec.get_datacube_bins(slitlength, minmax, numwav)
            # Generate the output WCS for the datacube
            crval_wv = cubepar['wave_min'] if cubepar['wave_min'] is not None else 1.0E10 * frame_wcs.wcs.crval[2]
            cd_wv = cubepar['wave_delta'] if cubepar['wave_delta'] is not None else 1.0E10 * frame_wcs.wcs.cd[2, 2]
            output_wcs = spec.get_wcs(spec2DObj.head0, slits, detector.platescale, crval_wv, cd_wv)
            # Set the wavelength range of the white light image.
            wl_wvrng = None
            if cubepar['save_whitelight']:
                wl_wvrng = get_whitelight_range(np.max(mnmx_wv[ff, :, 0]),
                                                np.min(mnmx_wv[ff, :, 1]),
                                                cubepar['whitelight_range'])
            # Make the datacube
            if method in ['subpixel', 'ngp']:
                # Generate the datacube
                generate_cube_subpixel(outfile, output_wcs, raimg[onslit_gpm], decimg[onslit_gpm], wave_ext,
                                       flux_sav[resrt], ivar_sav[resrt], np.ones(numpix),
                                       this_spatpos, this_specpos, this_spatid,
                                       spec2DObj.tilts, slits, alignSplines, bins,
                                       all_idx=None, overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                                       fluxcal=fluxcal, specname=specname, whitelight_range=wl_wvrng,
                                       spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
            continue

        # Store the information if we are combining multiple frames
        all_ra = np.append(all_ra, raimg[onslit_gpm].copy())
        all_dec = np.append(all_dec, decimg[onslit_gpm].copy())
        all_wave = np.append(all_wave, wave_ext.copy())
        all_sci = np.append(all_sci, flux_sav[resrt].copy())
        all_ivar = np.append(all_ivar, ivar_sav[resrt].copy())
        all_idx = np.append(all_idx, ff*np.ones(numpix))
        all_wghts = np.append(all_wghts, weights[ff]*np.ones(numpix)/weights[0])
        all_spatpos = np.append(all_spatpos, this_spatpos)
        all_specpos = np.append(all_specpos, this_specpos)
        all_spatid = np.append(all_spatid, this_spatid)
        all_tilts.append(spec2DObj.tilts)
        all_slits.append(slits)
        all_align.append(alignSplines)

    # No need to continue if we are not combining nor aligning frames
    if not combine and not align:
        return

    # Grab cos(dec) for convenience
    cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)

    # Register spatial offsets between all frames
    if align:
        if opts['ra_offset'] is not None:
            # First, translate all coordinates to the coordinates of the first frame
            # Note :: Don't need cosdec here, this just overrides the IFU coordinate centre of each frame
            ref_shift_ra = ifu_ra[0] - ifu_ra
            ref_shift_dec = ifu_dec[0] - ifu_dec
            for ff in range(numfiles):
                # Apply the shift
                all_ra[all_idx == ff] += ref_shift_ra[ff] + opts['ra_offset'][ff]/3600.0
                all_dec[all_idx == ff] += ref_shift_dec[ff] + opts['dec_offset'][ff]/3600.0
                msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f} E, {2:+0.3f} N".format(ff + 1, opts['ra_offset'][ff], opts['dec_offset'][ff]))
        else:
            # Find the wavelength range where all frames overlap
            min_wl, max_wl = get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                                  np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                                  cubepar['whitelight_range'])  # The user-specified values (if any)
            # Get the good whitelight pixels
            ww, wavediff = get_whitelight_pixels(all_wave, min_wl, max_wl)
            # Iterate over white light image generation and spatial shifting
            numiter = 2
            for dd in range(numiter):
                msgs.info(f"Iterating on spatial translation - ITERATION #{dd+1}/{numiter}")
                # Setup the WCS to use for all white light images
                ref_idx = None  # Don't use an index - This is the default behaviour when a reference image is supplied
                image_wcs, voxedge, reference_image = create_wcs(cubepar, all_ra[ww], all_dec[ww], all_wave[ww],
                                                                 dspat, wavediff, collapse=True)
                if voxedge[2].size != 2:
                    msgs.error("Spectral range for WCS is incorrect for white light image")

                wl_imgs = generate_image_subpixel(image_wcs, all_ra[ww], all_dec[ww], all_wave[ww],
                                                  all_sci[ww], all_ivar[ww], all_wghts[ww],
                                                  all_spatpos[ww], all_specpos[ww], all_spatid[ww],
                                                  all_tilts, all_slits, all_align, voxedge, all_idx=all_idx[ww],
                                                  spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
                if reference_image is None:
                    # ref_idx will be the index of the cube with the highest S/N
                    ref_idx = np.argmax(weights)
                    reference_image = wl_imgs[:, :, ref_idx].copy()
                    msgs.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
                else:
                    msgs.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")

                # Calculate the image offsets relative to the reference image
                for ff in range(numfiles):
                    # Calculate the shift
                    ra_shift, dec_shift = calculate_image_phase(reference_image.copy(), wl_imgs[:, :, ff], maskval=0.0)
                    # Convert pixel shift to degrees shift
                    ra_shift *= dspat/cosdec
                    dec_shift *= dspat
                    msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f} E, {2:+0.3f} N".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
                    # Apply the shift
                    all_ra[all_idx == ff] += ra_shift
                    all_dec[all_idx == ff] += dec_shift

    # Calculate the relative spectral weights of all pixels
    if numfiles == 1:
        # No need to calculate weights if there's just one frame
        all_wghts = np.ones_like(all_sci)
    else:
        # Find the wavelength range where all frames overlap
        min_wl, max_wl = get_whitelight_range(np.max(mnmx_wv[:, :, 0]),  # The max blue wavelength
                                              np.min(mnmx_wv[:, :, 1]),  # The min red wavelength
                                              cubepar['whitelight_range'])  # The user-specified values (if any)
        # Get the good white light pixels
        ww, wavediff = get_whitelight_pixels(all_wave, min_wl, max_wl)
        # Get a suitable WCS
        image_wcs, voxedge, reference_image = create_wcs(cubepar, all_ra, all_dec, all_wave, dspat, wavediff, collapse=True)
        # Generate the white light image (note: hard-coding subpixel=1 in both directions, and combining into a single image)
        wl_full = generate_image_subpixel(image_wcs, all_ra, all_dec, all_wave,
                                          all_sci, all_ivar, all_wghts,
                                          all_spatpos, all_specpos, all_spatid,
                                          all_tilts, all_slits, all_align, voxedge, all_idx=all_idx,
                                          spec_subpixel=1, spat_subpixel=1, combine=True)
        # Compute the weights
        all_wghts = compute_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_idx, wl_full[:, :, 0],
                                    dspat, dwv, relative_weights=cubepar['relative_weights'])

    # Generate the WCS, and the voxel edges
    cube_wcs, vox_edges, _ = create_wcs(cubepar, all_ra, all_dec, all_wave, dspat, dwv)

    sensfunc = None
    if flux_spline is not None:
        # Get wavelength of each pixel, and note that the WCS gives this in m, so convert to Angstroms (x 1E10)
        numwav = vox_edges[2].size-1
        senswave = cube_wcs.spectral.wcs_pix2world(np.arange(numwav), 0)[0] * 1.0E10
        sensfunc = flux_spline(senswave)

    # Generate a datacube
    outfile = get_output_filename("", cubepar['output_filename'], True, -1)
    if method in ['subpixel', 'ngp']:
        # Generate the datacube
        wl_wvrng = None
        if cubepar['save_whitelight']:
            wl_wvrng = get_whitelight_range(np.max(mnmx_wv[:, :, 0]),
                                            np.min(mnmx_wv[:, :, 1]),
                                            cubepar['whitelight_range'])
        if combine:
            generate_cube_subpixel(outfile, cube_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar,
                                   np.ones(all_wghts.size),  # all_wghts,
                                   all_spatpos, all_specpos, all_spatid, all_tilts, all_slits, all_align, vox_edges,
                                   all_idx=all_idx, overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                                   fluxcal=fluxcal, sensfunc=sensfunc, specname=specname, whitelight_range=wl_wvrng,
                                   spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
        else:
            for ff in range(numfiles):
                outfile = get_output_filename("", cubepar['output_filename'], False, ff)
                ww = np.where(all_idx == ff)
                generate_cube_subpixel(outfile, cube_wcs, all_ra[ww], all_dec[ww], all_wave[ww], all_sci[ww], all_ivar[ww], np.ones(all_wghts[ww].size),
                                       all_spatpos[ww], all_specpos[ww], all_spatid[ww], all_tilts[ff], all_slits[ff], all_align[ff], vox_edges,
                                       all_idx=all_idx[ww], overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                                       fluxcal=fluxcal, sensfunc=sensfunc, specname=specname, whitelight_range=wl_wvrng,
                                       spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)

>>>>>>> 3d081acc5 (Revert "Merge branch 'nirspec' into APF_Levy")

    # Return the datacube, variance cube and bad pixel cube
    return flxcube, varcube, bpmcube
