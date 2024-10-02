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

from pypeit import msgs, utils, specobj, specobjs
from pypeit.core import coadd, extract, flux_calib

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
            The centre of the Gaussian along the x-coordinate when z=0 (units of pixels)
        yo (float):
            The centre of the Gaussian along the y-coordinate when z=0 (units of pixels)
        sigma_x (float):
            The standard deviation in the x-direction (units of pixels)
        sigma_y (float):
            The standard deviation in the y-direction (units of pixels)
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
    # Setup the fitting params - Estimate a starting point for the fit using a median filter
    med_filt_image = signal.medfilt2d(image, kernel_size=3)
    idx_max = np.unravel_index(np.argmax(med_filt_image), image.shape)
    initial_guess = (1, idx_max[0], idx_max[1], 2, 2, 0, 0)
    bounds = ([0, 0, 0, 0.5, 0.5, -np.pi, -np.inf],
              [np.inf, image.shape[0], image.shape[1], image.shape[0], image.shape[1], np.pi, np.inf])
    # Perform the fit
    # TODO :: May want to generate the image on a finer pixel scale first
    popt, pcov = opt.curve_fit(gaussian2D, (xx, yy), image.ravel() / wlscl, bounds=bounds, p0=initial_guess)
    # Generate a best fit model
    model = gaussian2D((xx, yy), *popt).reshape(image.shape) * wlscl
    # Return the fitting results
    return popt, pcov, model


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


def extract_point_source(wave, flxcube, ivarcube, bpmcube, wcscube, exptime,
                         subpixel=20, boxcar_radius=None, optfwhm=None, whitelight_range=None,
                         pypeline="SlicerIFU", fluxed=False):
    """
    Extract a spectrum of a standard star from a datacube

    Parameters
    ----------
    wave : `numpy.ndarray`_
        Wavelength array for the datacube
    flxcube : `numpy.ndarray`_
        Datacube of the flux
    ivarcube : `numpy.ndarray`_
        Datacube of the inverse variance
    bpmcube : `numpy.ndarray`_
        Datacube of the bad pixel mask
    wcscube : `astropy.wcs.WCS`_
        WCS of the datacube
    exptime : float
        Exposure time listed in the header of the datacube
    subpixel : int, optional
        Number of pixels to subpixelate spectrum when creating mask
    boxcar_radius : float, optional
        Radius of the circular boxcar (in arcseconds) to use for the extraction
    optfwhm : float, optional
        FWHM of the PSF in pixels that is used to generate a Gaussian profile
        for the optimal extraction.
    pypeline : str, optional
        PypeIt pipeline used to reduce the datacube
    fluxed : bool, optional
        Is the datacube fluxed?

    Returns
    -------
    sobjs : :class:`~pypeit.specobjs.SpecObjs`
        SpecObjs object containing the extracted spectrum
    """
    if whitelight_range is None:
        whitelight_range = [np.min(wave), np.max(wave)]

    # Generate a spec1d object to hold the extracted spectrum
    msgs.info("Initialising a PypeIt SpecObj spec1d file")
    sobj = specobj.SpecObj(pypeline, "DET01", SLITID=0)
    sobj.RA = wcscube.wcs.crval[0]
    sobj.DEC = wcscube.wcs.crval[1]
    sobj.SLITID = 0

    # Convert from counts/s/Ang/arcsec**2 to counts. The sensitivity function expects counts as input
    numxx, numyy, numwave = flxcube.shape
    arcsecSQ = (wcscube.wcs.cdelt[0] * wcscube.wcs.cunit[0].to(units.arcsec)) * \
               (wcscube.wcs.cdelt[1] * wcscube.wcs.cunit[1].to(units.arcsec))
    if fluxed:
        # The datacube is flux calibrated, in units of 10^-17 erg/s/cm**2/Ang/arcsec**2
        # Scale the flux and ivar cubes to be in units of erg/s/cm**2/Ang
        unitscale = arcsecSQ
    else:
        # Scale the flux and ivar cubes to be in units of counts. pypeit_sensfunc expects counts as input
        deltawave = wcscube.wcs.cdelt[2]*wcscube.wcs.cunit[2].to(units.Angstrom)
        unitscale = exptime * deltawave * arcsecSQ

    # Apply the relevant scaling
    _flxcube = flxcube * unitscale
    _ivarcube = ivarcube / unitscale**2

    # Calculate the variance cube
    _varcube = utils.inverse(_ivarcube)

    # Generate a whitelight image, and fit a 2D Gaussian to estimate centroid and width
    msgs.info("Making white light image")
    wl_img = make_whitelight_fromcube(_flxcube, bpmcube, wave=wave, wavemin=whitelight_range[0], wavemax=whitelight_range[1])
    popt, pcov, model = fitGaussian2D(wl_img, norm=True)
    if boxcar_radius is None:
        nsig = 4  # 4 sigma should be far enough... Note: percentage enclosed for 2D Gaussian = 1-np.exp(-0.5 * nsig**2)
        wid = nsig * max(popt[3], popt[4])
    else:
        # Set the user-defined radius
        wid = boxcar_radius / np.sqrt(arcsecSQ)
    # Set the radius of the extraction boxcar for the sky determination
    msgs.info("Using a boxcar radius of {:0.2f} arcsec".format(wid*np.sqrt(arcsecSQ)))
    widsky = 2 * wid

    # Setup the coordinates of the mask
    x = np.linspace(0, numxx - 1, numxx * subpixel)
    y = np.linspace(0, numyy - 1, numyy * subpixel)
    xx, yy = np.meshgrid(x, y, indexing='ij')

    # Generate a mask
    msgs.info("Generating an object mask")
    newshape = (numxx * subpixel, numyy * subpixel)
    mask = np.zeros(newshape)
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < wid))
    mask[ww] = 1
    mask = utils.rebinND(mask, (numxx, numyy)).reshape(numxx, numyy, 1)

    # Generate a sky mask
    msgs.info("Generating a sky mask")
    newshape = (numxx * subpixel, numyy * subpixel)
    smask = np.zeros(newshape)
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < widsky))
    smask[ww] = 1
    smask = utils.rebinND(smask, (numxx, numyy)).reshape(numxx, numyy, 1)
    # Subtract off the object mask region, so that we just have an annulus around the object
    smask -= mask

    msgs.info("Subtracting the residual sky")
    # Subtract the residual sky from the datacube
    skymask = np.logical_not(bpmcube) * smask
    skycube = _flxcube * skymask
    skyspec = skycube.sum(axis=(0,1))
    nrmsky = skymask.sum(axis=(0,1))
    skyspec *= utils.inverse(nrmsky)
    _flxcube -= skyspec.reshape((1, 1, numwave))
    # Now subtract the residual sky from the white light image
    sky_val = np.sum(wl_img[:, :, np.newaxis] * smask) / np.sum(smask)
    wl_img -= sky_val

    msgs.info("Extracting a boxcar spectrum of datacube")
    # Construct an image that contains the fraction of flux included in the
    # boxcar extraction at each wavelength interval
    norm_flux = wl_img[:,:,np.newaxis] * mask
    norm_flux /= np.sum(norm_flux)
    # Extract boxcar
    cntmask = np.logical_not(bpmcube) * mask  # Good pixels within the masked region around the standard star
    flxscl = (norm_flux * cntmask).sum(axis=(0,1))  # This accounts for the flux that is missing due to masked pixels
    scimask = _flxcube * cntmask
    varmask = _varcube * cntmask**2
    nrmcnt = utils.inverse(flxscl)
    box_flux = scimask.sum(axis=(0,1)) * nrmcnt
    box_var = varmask.sum(axis=(0,1)) * nrmcnt**2
    box_gpm = flxscl > 1/3  # Good pixels are those where at least one-third of the standard star flux is measured

    # Store the BOXCAR extraction information
    sobj.BOX_RADIUS = wid  # Size of boxcar radius in pixels
    sobj.BOX_WAVE = wave.astype(float)
    if fluxed:
        sobj.BOX_FLAM = box_flux
        sobj.BOX_FLAM_SIG = np.sqrt(box_var)
        sobj.BOX_FLAM_IVAR = utils.inverse(box_var)
    else:
        sobj.BOX_COUNTS = box_flux
        sobj.BOX_COUNTS_SIG = np.sqrt(box_var)
        sobj.BOX_COUNTS_IVAR = utils.inverse(box_var)
        sobj.BOX_COUNTS_SKY = skyspec  # This is not the real sky, it is the residual sky. The datacube is presumed to be sky subtracted
    sobj.BOX_MASK = box_gpm
    sobj.S2N = np.median(box_flux * np.sqrt(utils.inverse(box_var)))

    # Now do the OPTIMAL extraction
    msgs.info("Extracting an optimal spectrum of datacube")
    # First, we need to rearrange the datacube and inverse variance cube into a 2D array.
    # The 3D -> 2D conversion is done so that there is a spectral and spatial dimension,
    # and the brightest white light pixel is transformed to be at the centre column of the 2D
    # array. Then, the second brightest white light pixel is transformed to be next to the centre
    # column of the 2D array, and so on. This is done so that the optimal extraction algorithm
    # can be applied.
    optkern = wl_img
    if optfwhm is not None:
        msgs.info("Generating a 2D Gaussian kernel for the optimal extraction, with FWHM = {:.2f} pixels".format(optfwhm))
        x = np.linspace(0, wl_img.shape[0] - 1, wl_img.shape[0])
        y = np.linspace(0, wl_img.shape[1] - 1, wl_img.shape[1])
        xx, yy = np.meshgrid(x, y, indexing='ij')
        # Generate a Gaussian kernel
        intflux = 1
        xo, yo = popt[1], popt[2]
        fwhm2sigma = 1.0 / (2 * np.sqrt(2 * np.log(2)))
        sigma_x, sigma_y = optfwhm*fwhm2sigma, optfwhm*fwhm2sigma
        theta, offset, = 0.0, 0.0
        optkern = gaussian2D((xx, yy), intflux, xo, yo, sigma_x, sigma_y, theta, offset).reshape(wl_img.shape)
        # Normalise the kernel
        optkern /= np.sum(optkern)

    optkern_masked = optkern * mask[:,:,0]
    # Normalise the white light image
    optkern_masked /= np.sum(optkern_masked)
    asrt = np.argsort(optkern_masked, axis=None)
    # Need to ensure that the number of pixels in the object profile is even
    if asrt.size % 2 != 0:
        # Remove the pixel with the lowest kernel weight.
        # It should be a zero value (given the mask), so it doesn't matter if we remove it
        asrt = asrt[1:]
    # Now sort the indices of the pixels in the object profile
    tmp = asrt.reshape((asrt.size//2, 2))
    objprof_idx = np.append(tmp[:,0], tmp[::-1,1])
    objprof = optkern_masked[np.unravel_index(objprof_idx, optkern.shape)]

    # Now slice the datacube and inverse variance cube into a 2D array
    spat, spec = np.meshgrid(objprof_idx, np.arange(numwave), indexing='ij')
    spatspl = np.apply_along_axis(np.unravel_index, 1, spat, optkern.shape)
    # Now slice the datacube and corresponding cubes/vectors into a series of 2D arrays
    numspat = objprof_idx.size
    flxslice = (spatspl[:,0,:], spatspl[:,1,:], spec)
    flxcube2d = _flxcube[flxslice].T
    ivarcube2d = _ivarcube[flxslice].T
    gpmcube2d = np.logical_not(bpmcube[flxslice].T)
    waveimg = wave.reshape((numwave,1)).repeat(numspat, axis=1)
    skyimg = np.zeros((numwave, numspat))  # Note, the residual sky has already been subtracted off _flxcube
    oprof = objprof.reshape((1,numspat)).repeat(numwave, axis=0)
    thismask = np.ones_like(flxcube2d, dtype=bool)

    # Now do the optimal extraction
    extract.extract_optimal(flxcube2d, ivarcube2d, gpmcube2d, waveimg, skyimg, thismask, oprof,
                            sobj, min_frac_use=0.05, fwhmimg=None, base_var=None, count_scale=None, noise_floor=None)

    # TODO :: The optimal extraction may suffer from residual DAR correction issues. This is because the
    #      :: object profile assumes that the white light image represents the true spatial profile of the
    #      :: object. One possibility is to fit a (linear?) model to the ratio of box/optimal extraction
    #      :: and then apply this model to the optimal extraction. This is a bit of a fudge.
    # Note that extract.extract_optimal() stores the optimal extraction in the
    # sobj.OPT_COUNTS, sobj.OPT_COUNTS_SIG, and sobj.OPT_COUNTS_IVAR attributes.
    # We need to store the fluxed extraction into the FLAM attributes (a slight fudge).
    if fluxed:
        sobj.OPT_FLAM = sobj.OPT_COUNTS
        sobj.OPT_FLAM_SIG = sobj.OPT_COUNTS_SIG
        sobj.OPT_FLAM_IVAR = sobj.OPT_COUNTS_IVAR

    # Make a specobjs object
    sobjs = specobjs.SpecObjs()
    sobjs.add_sobj(sobj)
    # Return the specobj object
    return sobjs


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


def make_whitelight_fromcube(cube, bpmcube, wave=None, wavemin=None, wavemax=None):
    """
    Generate a white light image using an input cube.

    Args:
        cube (`numpy.ndarray`_):
            3D datacube (the final element contains the wavelength dimension)
        bpmcube (`numpy.ndarray`_):
            3D bad pixel mask cube (the final element contains the wavelength dimension).
            A value of 1 indicates a bad pixel.
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
        # Cut the bad pixel mask and convert it to a good pixel mask
        cutgpmcube = np.logical_not(bpmcube[:, :, wmin:wmax])
    else:
        cutcube = cube.copy()
        cutgpmcube = np.logical_not(bpmcube)
    # Now sum along the wavelength axis
    nrmval = np.sum(cutgpmcube, axis=2)
    nrmval[nrmval == 0] = 1.0
    wl_img = np.sum(cutcube*cutgpmcube, axis=2) / nrmval
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
                           wl_full, dspat, dwv,
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
        `numpy.ndarray`_: The white light images for all frames. If combine=True,
        this will be a single 2D image. Otherwise, it will be a 3D array with
        dimensions (numra, numdec, numframes).
    """
    # Perform some checks on the input -- note, more complete checks are performed in subpixellate()
    _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg, _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset = \
        check_inputs([sciImg, ivarImg, waveImg, slitid_img_gpm, wghtImg, all_wcs, tilts, slits, astrom_trans, all_dar, ra_offset, dec_offset])

    # Generate the white light images
    if combine:
        # Subpixellate
        img, _, _ = subpixellate(image_wcs, bins, _sciImg, _ivarImg, _waveImg, _slitid_img_gpm, _wghtImg,
                                 _all_wcs, _tilts, _slits, _astrom_trans, _all_dar, _ra_offset, _dec_offset,
                                 spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel, slice_subpixel=slice_subpixel,
                                 skip_subpix_weights=True, correct_dar=correct_dar)
        return img[:, :, 0]
    else:
        # Prepare the array of white light images to be stored
        numframes = len(_sciImg)
        numra = bins[0].size - 1
        numdec = bins[1].size - 1
        all_wl_imgs = np.zeros((numra, numdec, numframes))
        # Loop through all frames and generate white light images
        for fr in range(numframes):
            msgs.info(f"Creating image {fr + 1}/{numframes}")
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
        whitelight_img = make_whitelight_fromcube(flxcube, bpmcube, wave=wave, wavemin=whitelight_range[0], wavemax=whitelight_range[1])
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
        variance cube, and (3) the corresponding bad pixel mask cube.
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
            asrt = np.argsort(yspl, kind='stable')
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
            ssrt = np.argsort(spatpos, kind='stable')
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

    # Return the datacube, variance cube and bad pixel cube
    return flxcube, varcube, bpmcube
