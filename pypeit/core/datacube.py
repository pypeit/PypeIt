"""
Module containing routines used by 3D datacubes.

.. include:: ../include/links.rst
"""

import os
import copy
import inspect

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import alignframe, datamodel, flatfield, io, specobj, spec2dobj, utils
from pypeit.core.flexure import calculate_image_phase
from pypeit.core import coadd, extract, findobj_skymask, flux_calib, parse, skysub
from pypeit.core.procimg import grow_mask
from pypeit.spectrographs.util import load_spectrograph

# Use a fast histogram for speed!
try:
    from fast_histogram import histogramdd
except ImportError:
    histogramdd = None

from IPython import embed


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


def correct_dar(wave_arr, coord, obstime, location, pressure, temperature, rel_humidity,
                wave_ref=None, numgrid=10):
    """
    Apply a differental atmospheric refraction correction to the
    input ra/dec.

    This implementation is based on ERFA, which is called through
    astropy.

    .. todo::
        There's probably going to be issues when the RA angle is
        either side of RA=0.

    Parameters
    ----------
    wave_arr : `numpy.ndarray`_
        wavelengths to obtain ra and dec offsets
    coord : `astropy.coordinates.SkyCoord`_
        ra, dec positions at the centre of the field
    obstime : `astropy.time.Time`_
        time at the midpoint of observation
    location : `astropy.coordinates.EarthLocation`_
        observatory location
    pressure : :obj:`float`
        Outside pressure at `location`
    temperature : :obj:`float`
        Outside ambient air temperature at `location`
    rel_humidity : :obj:`float`
        Outside relative humidity at `location`. This should be between 0 to 1.
    wave_ref : :obj:`float`
        Reference wavelength (The DAR correction will be performed relative to this wavelength)
    numgrid : :obj:`int`
        Number of grid points to evaluate the DAR correction.

    Returns
    -------
    ra_diff : `numpy.ndarray`_
        Relative RA shift at each wavelength given by ``wave_arr``
    dec_diff : `numpy.ndarray`_
        Relative DEC shift at each wavelength given by ``wave_arr``
    """
    msgs.info("Performing differential atmospheric refraction correction")
    if wave_ref is None:
        wave_ref = 0.5*(wave_arr.min() + wave_arr.max())

    # First create the reference frame and wavelength grid
    coord_altaz = coord.transform_to(AltAz(obstime=obstime, location=location))
    wave_grid = np.linspace(wave_arr.min(), wave_arr.max(), numgrid) * units.AA
    # Prepare the fit
    ra_grid, dec_grid = np.zeros(numgrid), np.zeros(numgrid)
    datfit = np.array([coord_altaz.alt.value, coord_altaz.az.value])
    # Loop through all wavelengths
    for ww in range(numgrid):
        # Fit the differential
        args = (coord.ra.value, coord.dec.value, datfit, wave_grid[ww], obstime, location, pressure, temperature, rel_humidity)
        #b_popt, b_pcov = opt.curve_fit(dar_fitfunc, tmp, datfit, p0=(0.0, 0.0))
        res_lsq = opt.least_squares(dar_fitfunc, [0.0, 0.0], args=args, xtol=1.0e-10, ftol=None, gtol=None)
        if not res_lsq.success:
            msgs.warn("DAR correction failed")
        # Store the result
        ra_grid[ww] = res_lsq.x[0]
        dec_grid[ww] = res_lsq.x[1]

    # Generate spline of differentials
    spl_ra = interp1d(wave_grid, ra_grid, kind='cubic')
    spl_dec = interp1d(wave_grid, dec_grid, kind='cubic')

    # Evaluate the differentials at the input wave_arr
    ra_diff = spl_ra(wave_arr) - spl_ra(wave_ref)
    dec_diff = spl_dec(wave_arr) - spl_dec(wave_ref)

    return ra_diff, dec_diff


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


def gaussian2D_cube(tup, intflux, xo, yo, dxdz, dydz, sigma_x, sigma_y, theta, offset):
    """
    Fit a 2D Gaussian function to a datacube. This function assumes that each
    wavelength slice of the datacube is well-fit by a 2D Gaussian. The centre of
    the Gaussian is allowed to vary linearly as a function of wavelength.

    .. note::

        The integrated flux does not vary with wavelength.

    Args:
        tup (:obj:`tuple`):
            A three element tuple containing the x, y, and z locations of each
            pixel in the cube
        intflux (float):
            The Integrated flux of the Gaussian
        xo (float):
            The centre of the Gaussian along the x-coordinate when z=0
        yo (float):
            The centre of the Gaussian along the y-coordinate when z=0
        dxdz (float):
            The change of xo with increasing z
        dydz (float):
            The change of yo with increasing z
        sigma_x (float):
            The standard deviation in the x-direction
        sigma_y (float):
            The standard deviation in the y-direction
        theta (float):
            The orientation angle of the 2D Gaussian
        offset (float):
            Constant offset

    Returns:
        `numpy.ndarray`_: The 2D Gaussian evaluated at the coordinate (x, y, z)
    """
    # Extract the (x, y, z) coordinates of each pixel from the tuple
    (x, y, z) = tup
    # Calculate the centre of the Gaussian for each z coordinate
    xo = float(xo) + z*dxdz
    yo = float(yo) + z*dydz
    # Account for a rotated 2D Gaussian
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    # Normalise so that the integrated flux is a parameter, instead of the amplitude
    norm = 1/(2*np.pi*np.sqrt(a*c-b*b))
    gtwod = offset + norm*intflux*np.exp(-(a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return gtwod.ravel()


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


def extract_standard_spec(stdcube, subpixel=20, method='boxcar'):
    """
    Extract a spectrum of a standard star from a datacube

    Parameters
    ----------
    std_cube : `astropy.io.fits.HDUList`_
        An HDU list of fits files
    subpixel : int
        Number of pixels to subpixelate spectrum when creating mask
    method : str
        Method used to extract standard star spectrum. Currently, only 'boxcar'
        is supported

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
    wcs_wav = stdwcs.wcs_pix2world(np.vstack((np.zeros(numwave), np.zeros(numwave), np.arange(numwave))).T, 0)
    wave = wcs_wav[:, 2] * 1.0E10 * units.AA

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
    sky_val = np.sum(wl_img[:,:,np.newaxis] * smask) / np.sum(smask)
    wl_img -= sky_val

    if method == 'boxcar':
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
    elif method == 'gauss2d':
        msgs.error("Use method=boxcar... this method has not been thoroughly tested")
        # Generate a mask
        fitmask = np.logical_not(bpmcube) * mask
        # Setup the coordinates
        x = np.linspace(0, flxcube.shape[0] - 1, flxcube.shape[0])
        y = np.linspace(0, flxcube.shape[1] - 1, flxcube.shape[1])
        z = np.linspace(0, flxcube.shape[2] - 1, flxcube.shape[2])
        xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
        # Normalise the flux in each wavelength channel
        scispec = (flxcube * fitmask).sum(0).sum(0).reshape((1, 1, flxcube.shape[2]))
        cntspec = fitmask.sum(0).sum(0).reshape((1, 1, flxcube.shape[2]))
        # These operations are all inverted, because we need to divide flxcube by scispec
        cntspec *= utils.inverse(scispec)
        cubefit = flxcube * cntspec
        cubesigfit = np.sqrt(varcube) * cntspec
        # Setup the fit params
        ww = np.where(fitmask)
        initial_guess = (1, idx_max[0], idx_max[1], 0.0, 0.0, 2, 2, 0, 0)
        bounds = ([-np.inf, 0, 0, -np.inf, -np.inf, 0.5, 0.5, -np.pi, -np.inf],
                  [np.inf,wl_img.shape[0],wl_img.shape[1],np.inf, np.inf, wl_img.shape[0],wl_img.shape[0],np.pi,np.inf])
        msgs.info("Fitting a 2D Gaussian to the datacube")
        popt, pcov = opt.curve_fit(gaussian2D_cube, (xx[ww], yy[ww], zz[ww]), cubefit[ww],
                                   sigma=cubesigfit[ww], bounds=bounds, p0=initial_guess)
        # Subtract off the best-fitting continuum
        popt[-1] = 0
        # Generate the best-fitting model to be used as an optimal profile
        model = gaussian2D_cube((xx, yy, zz), *popt).reshape(flxcube.shape)
        numim = flxcube.shape[0]*flxcube.shape[1]

        # Optimally extract
        msgs.info("Optimally extracting...")
        sciimg = (flxcube*mask).reshape((numim, numwave)).T
        ivar = utils.inverse((varcube*mask**2).reshape((numim, numwave)).T)
        optmask = fitmask.reshape((numim, numwave)).T
        waveimg = np.ones((numwave, numim))  # Just a dummy array - not needed
        skyimg = np.zeros((numwave, numim))  # Just a dummy array - not needed
        thismask = np.ones((numwave, numim))  # Just a dummy array - not needed
        oprof = model.reshape((numim, numwave)).T
        sobj = specobj.SpecObj('SlicerIFU', 'DET01', SLITID=0)
        extract.extract_optimal(sciimg, ivar, optmask, waveimg, skyimg, thismask, oprof, sobj)
        opt_flux, opt_var, opt_gpm = sobj.OPT_COUNTS, sobj.OPT_COUNTS_SIG**2, sobj.OPT_MASK
        # Setup the return values
        ret_flux, ret_var, ret_gpm = opt_flux, opt_var, opt_gpm
    elif method == 'optimal':
        msgs.error("Use method=boxcar... this method has not been thoroughly tested")
        # First do a boxcar along one dimension
        msgs.info("Collapsing datacube to a 2D image")
        omask = mask+smask
        idx_sum = 0
        cntmask = np.logical_not(bpmcube) * omask
        scimask = flxcube * cntmask
        varmask = varcube * cntmask**2
        cnt_spec = cntmask.sum(idx_sum) * utils.inverse(omask.sum(idx_sum))
        nrmcnt = utils.inverse(cnt_spec)
        box_sciimg = scimask.sum(idx_sum) * nrmcnt
        box_scivar = varmask.sum(idx_sum) * nrmcnt**2
        box_sciivar = utils.inverse(box_scivar)
        # Transpose for optimal
        box_sciimg = box_sciimg.T
        box_sciivar = box_sciivar.T

        # Prepare for optimal
        msgs.info("Starting optimal extraction")
        thismask = np.ones(box_sciimg.shape, dtype=bool)
        nspec, nspat = thismask.shape[0], thismask.shape[1]
        slit_left = np.zeros(nspec)
        slit_right = np.ones(nspec)*(nspat-1)
        tilts = np.outer(np.linspace(0.0,1.0,nspec), np.ones(nspat))
        waveimg = np.outer(wave.value, np.ones(nspat))
        global_sky = np.zeros_like(box_sciimg)
        # Find objects and then extract
        sobj = findobj_skymask.objs_in_slit(box_sciimg, thismask, slit_left, slit_right)
        skysub.local_skysub_extract(box_sciimg, box_sciivar, tilts, waveimg, global_sky, thismask, slit_left,
                             slit_right, sobj, model_noise=False)
        opt_flux, opt_var, opt_gpm = sobj.OPT_COUNTS[0,:], sobj.OPT_COUNTS_SIG[0,:]**2, sobj.OPT_MASK[0,:]
        # Setup the return values
        ret_flux, ret_var, ret_gpm = opt_flux, opt_var, opt_gpm
    else:
        msgs.error("Unknown extraction method: ", method)

    # Convert from counts/s/Ang/arcsec**2 to counts/s/Ang
    arcsecSQ = 3600.0*3600.0*(stdwcs.wcs.cdelt[0]*stdwcs.wcs.cdelt[1])
    ret_flux *= arcsecSQ
    ret_var *= arcsecSQ**2
    # Return the box extraction results
    return wave, ret_flux, utils.inverse(ret_var), ret_gpm


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
        `numpy.ndarray`_: Whitelight image of the input cube.
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


def make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                               all_ivar=None, whitelightWCS=None, numra=None, numdec=None, trim=1):
    """
    Generate a whitelight image using the individual pixels of every input frame

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
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights attributed to each pixel
            from all spec2d files
        all_idx (`numpy.ndarray`_):
            1D flattened array containing an integer identifier indicating which
            spec2d file each pixel originates from. For example, a 0 would
            indicate that a pixel originates from the first spec2d frame listed
            in the input file. a 1 would indicate that this pixel originates
            from the second spec2d file, and so forth.
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        all_ivar (`numpy.ndarray`_, optional):
            1D flattened array containing of the inverse variance of each pixel
            from all spec2d files.  If provided, inverse variance images will be
            calculated and returned for each white light image.
        whitelightWCS (`astropy.wcs.WCS`_, optional):
            The WCS of a reference white light image. If supplied, you must also
            supply numra and numdec.
        numra (int, optional):
            Number of RA spaxels in the reference white light image
        numdec (int, optional):
            Number of DEC spaxels in the reference white light image
        trim (int, optional):
            Number of pixels to grow around a masked region

    Returns:
        tuple: two 3D arrays will be returned, each of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light
        images.  The first array is a white light image, and the second array is
        the corresponding inverse variance image. If all_ivar is None, this will
        be an empty array.
    """
    # Determine number of files
    numfiles = np.unique(all_idx).size

    if whitelightWCS is None:
        # Generate a 2D WCS to register all frames
        coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
        coord_dlt = [dspat, dspat, np.max(all_wave) - np.min(all_wave)]
        whitelightWCS = generate_WCS(coord_min, coord_dlt)

        # Generate coordinates
        cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)
        numra = 1+int((np.max(all_ra) - np.min(all_ra)) * cosdec / dspat)
        numdec = 1+int((np.max(all_dec) - np.min(all_dec)) / dspat)
    else:
        # If a WCS is supplied, the numra and numdec must be specified
        if (numra is None) or (numdec is None):
            msgs.error("A WCS has been supplied to make_whitelight." + msgs.newline() +
                       "numra and numdec must also be specified")
    xbins = np.arange(1 + numra) - 1
    ybins = np.arange(1 + numdec) - 1
    spec_bins = np.arange(2) - 1
    bins = (xbins, ybins, spec_bins)

    whitelight_Imgs = np.zeros((numra, numdec, numfiles))
    whitelight_ivar = np.zeros((numra, numdec, numfiles))
    for ff in range(numfiles):
        msgs.info("Generating white light image of frame {0:d}/{1:d}".format(ff + 1, numfiles))
        ww = (all_idx == ff)
        # Make the cube
        pix_coord = whitelightWCS.wcs_world2pix(np.vstack((all_ra[ww], all_dec[ww], all_wave[ww] * 1.0E-10)).T, 0)
        wlcube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci[ww] * all_wghts[ww])
        norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts[ww])
        nrmCube = (norm > 0) / (norm + (norm == 0))
        whtlght = (wlcube * nrmCube)[:, :, 0]
        # Create a mask of good pixels (trim the edges)
        gpm = grow_mask(whtlght == 0, trim) == 0  # A good pixel = 1
        whtlght *= gpm
        # Set the masked regions to the minimum value
        minval = np.min(whtlght[gpm == 1])
        whtlght[gpm == 0] = minval
        # Store the white light image
        whitelight_Imgs[:, :, ff] = whtlght.copy()
        # Now operate on the inverse variance image
        if all_ivar is not None:
            ivar_img, _ = np.histogramdd(pix_coord, bins=bins, weights=all_ivar[ww])
            ivar_img = ivar_img[:, :, 0]
            ivar_img *= gpm
            minval = np.min(ivar_img[gpm == 1])
            ivar_img[gpm == 0] = minval
            whitelight_ivar[:, :, ff] = ivar_img.copy()
    return whitelight_Imgs, whitelight_ivar, whitelightWCS


def create_wcs(cubepar, all_ra, all_dec, all_wave, dspat, dwv, collapse=False, equinox=2000.0,
               specname="PYP_SPEC"):
    """
    Create a WCS and the expected edges of the voxels, based on user-specified
    parameters or the extremities of the data.

    Parameters
    ----------
    cubepar : :class:`~pypeit.par.pypeitpar.CubePar`
        An instance of the CubePar parameter set, contained parameters of the
        datacube reduction
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
    dwv : float
        Linear wavelength step of each voxel (in Angstroms)
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
    reference_image = None  # The default behaviour is that the reference image is not used
    ra_min = cubepar['ra_min'] if cubepar['ra_min'] is not None else np.min(all_ra)
    ra_max = cubepar['ra_max'] if cubepar['ra_max'] is not None else np.max(all_ra)
    dec_min = cubepar['dec_min'] if cubepar['dec_min'] is not None else np.min(all_dec)
    dec_max = cubepar['dec_max'] if cubepar['dec_max'] is not None else np.max(all_dec)
    wav_min = cubepar['wave_min'] if cubepar['wave_min'] is not None else np.min(all_wave)
    wav_max = cubepar['wave_max'] if cubepar['wave_max'] is not None else np.max(all_wave)
    dwave = cubepar['wave_delta'] if cubepar['wave_delta'] is not None else dwv

    # Number of voxels in each dimension
    numra = int((ra_max-ra_min) * cosdec / dspat)
    numdec = int((dec_max-dec_min)/dspat)
    numwav = int(np.round((wav_max-wav_min)/dwave))

    # If a white light WCS is being generated, make sure there's only 1 wavelength bin
    if collapse:
        wav_min = np.min(all_wave)
        wav_max = np.max(all_wave)
        dwave = wav_max - wav_min
        numwav = 1

    # Generate a master WCS to register all frames
    coord_min = [ra_min, dec_min, wav_min]
    coord_dlt = [dspat, dspat, dwave]

    # If a reference image is being used and a white light image is requested (collapse=True) update the celestial parts
    if cubepar["reference_image"] is not None:
        # Load the requested reference image
        reference_image, imgwcs = load_imageWCS(cubepar["reference_image"])
        # Update the celestial WCS
        coord_min[:2] = imgwcs.wcs.crval
        coord_dlt[:2] = imgwcs.wcs.cdelt
        numra, numdec = reference_image.shape

    cubewcs = generate_WCS(coord_min, coord_dlt, equinox=equinox, name=specname)
    msgs.info(msgs.newline() + "-" * 40 +
              msgs.newline() + "Parameters of the WCS:" +
              msgs.newline() + "RA   min = {0:f}".format(coord_min[0]) +
              msgs.newline() + "DEC  min = {0:f}".format(coord_min[1]) +
              msgs.newline() + "WAVE min, max = {0:f}, {1:f}".format(wav_min, wav_max) +
              msgs.newline() + "Spaxel size = {0:f} arcsec".format(3600.0*dspat) +
              msgs.newline() + "Wavelength step = {0:f} A".format(dwave) +
              msgs.newline() + "-" * 40)

    # Generate the output binning
    xbins = np.arange(1+numra)-0.5
    ybins = np.arange(1+numdec)-0.5
    spec_bins = np.arange(1+numwav)-0.5
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
    wave_spec = wcs_res[:, 2] * 1.0E10
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


def generate_image_subpixel(image_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_spatpos, all_specpos,
                            all_spatid, tilts, slits, astrom_trans, bins, all_idx=None,
                            spec_subpixel=10, spat_subpixel=10, combine=False):
    """
    Generate a white light image from the input pixels

    Args:
        image_wcs (`astropy.wcs.WCS`_):
            World coordinate system to use for the white light images.
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the right ascension of each pixel
            (units = degrees)
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the declination of each pixel (units =
            degrees)
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength of each pixel (units =
            Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all
            spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel
            from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used
            in the combination
        all_spatpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spatial direction
        all_specpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spectral direction
        all_spatid (`numpy.ndarray`_):
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list):
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`~pypeit.slittrace.SlitTraceSet`, list):
            Information stored about the slits, or a list of SlitTraceSet (see
            all_idx)
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional):
            If tilts, slits, and astrom_trans are lists, this should contain a
            1D flattened array, of the same length as all_sci, containing the
            index the tilts, slits, and astrom_trans lists that corresponds to
            each pixel. Note that, in this case all of these lists need to be
            the same length.
        spec_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spectral direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral
            direction.
        spat_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spatial direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial
            direction.
        combine (:obj:`bool`, optional):
            If True, all of the input frames will be combined into a single
            output. Otherwise, individual images will be generated.

    Returns:
        `numpy.ndarray`_: The white light images for all frames
    """
    # Perform some checks on the input -- note, more complete checks are performed in subpixellate()
    _all_idx = np.zeros(all_sci.size) if all_idx is None else all_idx
    if combine:
        numfr = 1
    else:
        numfr = np.unique(_all_idx).size
        if len(tilts) != numfr or len(slits) != numfr or len(astrom_trans) != numfr:
            msgs.error("The following arguments must be the same length as the expected number of frames to be combined:"
                       + msgs.newline() + "tilts, slits, astrom_trans")
    # Prepare the array of white light images to be stored
    numra = bins[0].size-1
    numdec = bins[1].size-1
    all_wl_imgs = np.zeros((numra, numdec, numfr))

    # Loop through all frames and generate white light images
    for fr in range(numfr):
        msgs.info(f"Creating image {fr+1}/{numfr}")
        if combine:
            # Subpixellate
            img, _, _ = subpixellate(image_wcs, all_ra, all_dec, all_wave,
                                     all_sci, all_ivar, all_wghts, all_spatpos,
                                     all_specpos, all_spatid, tilts, slits, astrom_trans, bins,
                                     spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, all_idx=_all_idx)
        else:
            ww = np.where(_all_idx == fr)
            # Subpixellate
            img, _, _ = subpixellate(image_wcs, all_ra[ww], all_dec[ww], all_wave[ww],
                                     all_sci[ww], all_ivar[ww], all_wghts[ww], all_spatpos[ww],
                                     all_specpos[ww], all_spatid[ww], tilts[fr], slits[fr], astrom_trans[fr], bins,
                                     spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
        all_wl_imgs[:, :, fr] = img[:, :, 0]
    # Return the constructed white light images
    return all_wl_imgs


def generate_cube_subpixel(outfile, output_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts,
                           all_spatpos, all_specpos, all_spatid, tilts, slits, astrom_trans, bins,
                           all_idx=None, spec_subpixel=10, spat_subpixel=10, overwrite=False, blaze_wave=None,
                           blaze_spec=None, fluxcal=False, sensfunc=None, whitelight_range=None,
                           specname="PYP_SPEC", debug=False):
    r"""
    Save a datacube using the subpixel algorithm. Refer to the subpixellate()
    docstring for further details about this algorithm

    Args:
        outfile (str):
            Filename to be used to save the datacube
        output_wcs (`astropy.wcs.WCS`_):
            Output world coordinate system.
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the right ascension of each pixel
            (units = degrees)
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the declination of each pixel (units =
            degrees)
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength of each pixel (units =
            Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all
            spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel
            from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used
            in the combination
        all_spatpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spatial direction
        all_specpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spectral direction
        all_spatid (`numpy.ndarray`_):
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list):
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`~pypeit.slittrace.SlitTraceSet`, list):
            Information stored about the slits, or a list of SlitTraceSet (see
            all_idx)
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional):
            If tilts, slits, and astrom_trans are lists, this should contain a
            1D flattened array, of the same length as all_sci, containing the
            index the tilts, slits, and astrom_trans lists that corresponds to
            each pixel. Note that, in this case all of these lists need to be
            the same length.
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
        overwrite (bool, optional):
            If True, the output cube will be overwritten.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_, optional):
            Spectral blaze function
        fluxcal (bool, optional):
            Are the data flux calibrated? If True, the units are: :math:`{\rm
            erg/s/cm}^2{\rm /Angstrom/arcsec}^2` multiplied by the
            PYPEIT_FLUX_SCALE.  Otherwise, the units are: :math:`{\rm
            counts/s/Angstrom/arcsec}^2`.
        sensfunc (`numpy.ndarray`_, None, optional):
            Sensitivity function that has been applied to the datacube
        whitelight_range (None, list, optional):
            A two element list that specifies the minimum and maximum
            wavelengths (in Angstroms) to use when constructing the white light
            image (format is: [min_wave, max_wave]). If None, the cube will be
            collapsed over the full wavelength range. If a list is provided an
            either element of the list is None, then the minimum/maximum
            wavelength range of that element will be set by the minimum/maximum
            wavelength of all_wave.
        specname (str, optional):
            Name of the spectrograph
        debug (bool, optional):
            If True, a residuals cube will be output. If the datacube generation
            is correct, the distribution of pixels in the residual cube with no
            flux should have mean=0 and std=1.
    """
    # Prepare the header, and add the unit of flux to the header
    hdr = output_wcs.to_header()
    if fluxcal:
        hdr['FLUXUNIT'] = (flux_calib.PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
    else:
        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")

    # Subpixellate
    subpix = subpixellate(output_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_spatpos, all_specpos,
                          all_spatid, tilts, slits, astrom_trans, bins, all_idx=all_idx,
                          spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, debug=debug)
    # Extract the variables that we need
    if debug:
        datacube, varcube, bpmcube, residcube = subpix
        # Save a residuals cube
        outfile_resid = outfile.replace(".fits", "_resid.fits")
        msgs.info("Saving residuals datacube as: {0:s}".format(outfile_resid))
        hdu = fits.PrimaryHDU(residcube.T, header=hdr)
        hdu.writeto(outfile_resid, overwrite=overwrite)
    else:
        datacube, varcube, bpmcube = subpix

    # Check if the user requested a white light image
    if whitelight_range is not None:
        # Grab the WCS of the white light image
        whitelight_wcs = output_wcs.celestial
        # Determine the wavelength range of the whitelight image
        if whitelight_range[0] is None:
            whitelight_range[0] = np.min(all_wave)
        if whitelight_range[1] is None:
            whitelight_range[1] = np.max(all_wave)
        msgs.info("White light image covers the wavelength range {0:.2f} A - {1:.2f} A".format(
            whitelight_range[0], whitelight_range[1]))
        # Get the output filename for the white light image
        out_whitelight = get_output_whitelight_filename(outfile)
        nspec = datacube.shape[2]
        # Get wavelength of each pixel, and note that the WCS gives this in m, so convert to Angstroms (x 1E10)
        wave = 1.0E10 * output_wcs.spectral.wcs_pix2world(np.arange(nspec), 0)[0]
        whitelight_img = make_whitelight_fromcube(datacube, wave=wave, wavemin=whitelight_range[0], wavemax=whitelight_range[1])
        msgs.info("Saving white light image as: {0:s}".format(out_whitelight))
        img_hdu = fits.PrimaryHDU(whitelight_img.T, header=whitelight_wcs.to_header())
        img_hdu.writeto(out_whitelight, overwrite=overwrite)

    # Write out the datacube
    msgs.info("Saving datacube as: {0:s}".format(outfile))
    final_cube = DataCube(datacube.T, np.sqrt(varcube.T), bpmcube.T, specname, blaze_wave, blaze_spec,
                          sensfunc=sensfunc, fluxed=fluxcal)
    final_cube.to_file(outfile, hdr=hdr, overwrite=overwrite)


def subpixellate(output_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_spatpos, all_specpos,
                 all_spatid, tilts, slits, astrom_trans, bins, all_idx=None,
                 spec_subpixel=10, spat_subpixel=10, debug=False):
    r"""
    Subpixellate the input data into a datacube. This algorithm splits each
    detector pixel into multiple subpixels, and then assigns each subpixel to a
    voxel. For example, if ``spec_subpixel = spat_subpixel = 10``, then each
    detector pixel is divided into :math:`10^2=100` subpixels. Alternatively,
    when spec_subpixel = spat_subpixel = 1, this corresponds to the nearest grid
    point (NGP) algorithm.

    Important Note: If spec_subpixel > 1 or spat_subpixel > 1, the errors will
    be correlated, and the covariance is not being tracked, so the errors will
    not be (quite) right. There is a tradeoff one has to make between sampling
    and better looking cubes, versus no sampling and better behaved errors.

    Args:
        output_wcs (`astropy.wcs.WCS`_):
            Output world coordinate system.
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the right ascension of each pixel
            (units = degrees)
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the declination of each pixel (units =
            degrees)
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength of each pixel (units =
            Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all
            spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel
            from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used
            in the combination
        all_spatpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spatial direction
        all_specpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the
            spectral direction
        all_spatid (`numpy.ndarray`_):
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list):
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`~pypeit.slittrace.SlitTraceSet`, list):
            Information stored about the slits, or a list of SlitTraceSet (see
            all_idx)
        astrom_trans (:class:`~pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel
            coordinates and WCS pixel coordinates, or a list of Alignment
            Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial
            and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional):
            If tilts, slits, and astrom_trans are lists, this should contain a
            1D flattened array, of the same length as all_sci, containing the
            index the tilts, slits, and astrom_trans lists that corresponds to
            each pixel. Note that, in this case all of these lists need to be
            the same length.
        spec_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spectral direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral
            direction.
        spat_subpixel (:obj:`int`, optional):
            What is the subpixellation factor in the spatial direction. Higher
            values give more reliable results, but note that the time required
            goes as (``spec_subpixel * spat_subpixel``). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial
            direction.
        debug (bool):
            If True, a residuals cube will be output. If the datacube generation
            is correct, the distribution of pixels in the residual cube with no
            flux should have mean=0 and std=1.

    Returns:
        :obj:`tuple`: Three or four `numpy.ndarray`_ objects containing (1) the
        datacube generated from the subpixellated inputs, (2) the corresponding
        variance cube, (3) the corresponding bad pixel mask cube, and (4) the
        residual cube.  The latter is only returned if debug is True. 
    """
    # Check for combinations of lists or not
    if type(tilts) is list and type(slits) is list and type(astrom_trans) is list:
        # Several frames are being combined. Check the lists have the same length
        numframes = len(tilts)
        if len(slits) != numframes or len(astrom_trans) != numframes:
            msgs.error("The following lists must have the same length:" + msgs.newline() +
                       "tilts, slits, astrom_trans")
        # Check all_idx has been set
        if all_idx is None:
            if numframes != 1:
                msgs.error("Missing required argument for combining frames: all_idx")
            else:
                all_idx = np.zeros(all_sci.size)
        else:
            tmp = np.unique(all_idx).size
            if tmp != numframes:
                msgs.warn("Indices in argument 'all_idx' does not match the number of frames expected.")
        # Store in the following variables
        _tilts, _slits, _astrom_trans = tilts, slits, astrom_trans
    elif type(tilts) is not list and type(slits) is not list and \
            type(astrom_trans) is not list:
        # Just a single frame - store as lists for this code
        _tilts, _slits, _astrom_trans = [tilts], [slits], [astrom_trans],
        all_idx = np.zeros(all_sci.size)
        numframes = 1
    else:
        msgs.error("The following input arguments should all be of type 'list', or all not be type 'list':" +
                   msgs.newline() + "tilts, slits, astrom_trans")
    # Prepare the output arrays
    outshape = (bins[0].size-1, bins[1].size-1, bins[2].size-1)
    binrng = [[bins[0][0], bins[0][-1]], [bins[1][0], bins[1][-1]], [bins[2][0], bins[2][-1]]]
    datacube, varcube, normcube = np.zeros(outshape), np.zeros(outshape), np.zeros(outshape)
    if debug:
        residcube = np.zeros(outshape)
    # Divide each pixel into subpixels
    spec_offs = np.arange(0.5/spec_subpixel, 1, 1/spec_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    spat_offs = np.arange(0.5/spat_subpixel, 1, 1/spat_subpixel) - 0.5  # -0.5 is to offset from the centre of each pixel.
    spat_x, spec_y = np.meshgrid(spat_offs, spec_offs)
    num_subpixels = spec_subpixel * spat_subpixel
    area = 1 / num_subpixels
    all_wght_subpix = all_wghts * area
    all_var = utils.inverse(all_ivar)
    # Loop through all exposures
    for fr in range(numframes):
        # Extract tilts and slits for convenience
        this_tilts = _tilts[fr]
        this_slits = _slits[fr]
        # Loop through all slits
        for sl, spatid in enumerate(this_slits.spat_id):
            if numframes == 1:
                msgs.info(f"Resampling slit {sl+1}/{this_slits.nslits}")
            else:
                msgs.info(f"Resampling slit {sl+1}/{this_slits.nslits} of frame {fr+1}/{numframes}")
            this_sl = np.where((all_spatid == spatid) & (all_idx == fr))
            wpix = (all_specpos[this_sl], all_spatpos[this_sl])
            # Generate a spline between spectral pixel position and wavelength
            yspl = this_tilts[wpix]*(this_slits.nspec - 1)
            tiltpos = np.add.outer(yspl, spec_y).flatten()
            wspl = all_wave[this_sl]
            asrt = np.argsort(yspl)
            wave_spl = interp1d(yspl[asrt], wspl[asrt], kind='linear', bounds_error=False, fill_value='extrapolate')
            # Calculate spatial and spectral positions of the subpixels
            spat_xx = np.add.outer(wpix[1], spat_x.flatten()).flatten()
            spec_yy = np.add.outer(wpix[0], spec_y.flatten()).flatten()
            # Transform this to spatial location
            spatpos_subpix = _astrom_trans[fr].transform(sl, spat_xx, spec_yy)
            spatpos = _astrom_trans[fr].transform(sl, all_spatpos[this_sl], all_specpos[this_sl])
            ra_coeff = np.polyfit(spatpos, all_ra[this_sl], 1)
            dec_coeff = np.polyfit(spatpos, all_dec[this_sl], 1)
            this_ra = np.polyval(ra_coeff, spatpos_subpix)#ra_spl(spatpos_subpix)
            this_dec = np.polyval(dec_coeff, spatpos_subpix)#dec_spl(spatpos_subpix)
            # ssrt = np.argsort(spatpos)
            # ra_spl = interp1d(spatpos[ssrt], all_ra[this_sl][ssrt], kind='linear', bounds_error=False, fill_value='extrapolate')
            # dec_spl = interp1d(spatpos[ssrt], all_dec[this_sl][ssrt], kind='linear', bounds_error=False, fill_value='extrapolate')
            # this_ra = ra_spl(spatpos_subpix)
            # this_dec = dec_spl(spatpos_subpix)
            this_wave = wave_spl(tiltpos)
            # Convert world coordinates to voxel coordinates, then histogram
            vox_coord = output_wcs.wcs_world2pix(np.vstack((this_ra, this_dec, this_wave * 1.0E-10)).T, 0)
            if histogramdd is not None:
                # use the "fast histogram" algorithm, that assumes regular bin spacing
                datacube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(all_sci[this_sl] * all_wght_subpix[this_sl], num_subpixels))
                varcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(all_var[this_sl] * all_wght_subpix[this_sl]**2, num_subpixels))
                normcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(all_wght_subpix[this_sl], num_subpixels))
                if debug:
                    residcube += histogramdd(vox_coord, bins=outshape, range=binrng, weights=np.repeat(all_sci[this_sl] * np.sqrt(all_ivar[this_sl]), num_subpixels))
            else:
                datacube += np.histogramdd(vox_coord, bins=outshape, weights=np.repeat(all_sci[this_sl] * all_wght_subpix[this_sl], num_subpixels))[0]
                varcube += np.histogramdd(vox_coord, bins=outshape, weights=np.repeat(all_var[this_sl] * all_wght_subpix[this_sl]**2, num_subpixels))[0]
                normcube += np.histogramdd(vox_coord, bins=outshape, weights=np.repeat(all_wght_subpix[this_sl], num_subpixels))[0]
                if debug:
                    residcube += np.histogramdd(vox_coord, bins=outshape, weights=np.repeat(all_sci[this_sl] * np.sqrt(all_ivar[this_sl]), num_subpixels))[0]
    # Normalise the datacube and variance cube
    nc_inverse = utils.inverse(normcube)
    datacube *= nc_inverse
    varcube *= nc_inverse**2
    bpmcube = (normcube == 0).astype(np.uint8)
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

