"""
Module containing routines used by 3D datacubes.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import os
import copy
import inspect

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy.interpolate import interp1d, RegularGridInterpolator, RBFInterpolator
import numpy as np

from pypeit import msgs
from pypeit import alignframe, datamodel, flatfield, io, specobj, spec2dobj, utils, wavecalib
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


class DataCube(datamodel.DataContainer):
    """
    DataContainer to hold the products of a datacube

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_datacube.rst

    Args:
        flux (`numpy.ndarray`_):
            The science datacube (nwave, nspaxel_y, nspaxel_x)
        sig (`numpy.ndarray`_):
            The error datacube (nwave, nspaxel_y, nspaxel_x)
        bpm (`numpy.ndarray`_):
            The bad pixel mask of the datacube (nwave, nspaxel_y, nspaxel_x).
            True values indicate a bad pixel
        blaze_wave (`numpy.ndarray`_):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_):
            The spectral blaze function
        sensfunc (`numpy.ndarray`_, None):
            Sensitivity function (nwave,). Only saved if the data are fluxed.
        PYP_SPEC (str):
            Name of the PypeIt Spectrograph
        fluxed (bool):
            If the cube has been flux calibrated, this will be set to "True"

    Attributes:
        head0 (`astropy.io.fits.Header`):
            Primary header
        filename (str):
            Filename to use when loading from file
        spect_meta (:obj:`dict`):
            Parsed meta from the header
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Build from PYP_SPEC

    """
    version = '1.1.0'

    datamodel = {'flux': dict(otype=np.ndarray, atype=np.floating, descr='Flux datacube in units of counts/s/Ang/arcsec^2'
                                                                         'or 10^-17 erg/s/cm^2/Ang/arcsec^2'),
                 'sig': dict(otype=np.ndarray, atype=np.floating, descr='Error datacube (matches units of flux)'),
                 'bpm': dict(otype=np.ndarray, atype=np.uint8, descr='Bad pixel mask of the datacube (0=good, 1=bad)'),
                 'blaze_wave': dict(otype=np.ndarray, atype=np.floating, descr='Wavelength array of the spectral blaze function'),
                 'blaze_spec': dict(otype=np.ndarray, atype=np.floating, descr='The spectral blaze function'),
                 'sensfunc': dict(otype=np.ndarray, atype=np.floating, descr='Sensitivity function 10^-17 erg/(counts/cm^2)'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt: Spectrograph name'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the datacube is fluxed.')}

    internals = ['head0',
                 'filename',
                 'spectrograph',
                 'spect_meta'
                ]

    def __init__(self, flux, sig, bpm, PYP_SPEC, blaze_wave, blaze_spec, sensfunc=None, fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k, values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _bundle(self):
        """
        Over-write default _bundle() method to separate the DetectorContainer
        into its own HDU

        Returns:
            :obj:`list`: A list of dictionaries, each list element is
            written to its own fits extension. See the description
            above.
        """
        d = []
        # Rest of the datamodel
        for key in self.keys():
            # Skip Nones
            if self[key] is None:
                continue
            # Array?
            if self.datamodel[key]['otype'] == np.ndarray:
                tmp = {}
                if self.datamodel[key]['atype'] == np.floating:
                    tmp[key] = self[key].astype(np.float32)
                else:
                    tmp[key] = self[key]
                d.append(tmp)
            else:
                # Add to header of the primary image
                d[0][key] = self[key]
        # Return
        return d

    def to_file(self, ofile, primary_hdr=None, hdr=None, **kwargs):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.to_file`
        to deal with the header

        Args:
            ofile (:obj:`str`): Filename
            primary_hdr (`astropy.io.fits.Header`_, optional):
            wcs (`astropy.io.fits.Header`_, optional):
                The World Coordinate System, represented by a fits header
            **kwargs:  Passed to super.to_file()

        """
        if primary_hdr is None:
            primary_hdr = io.initialize_header()
        # Build the header
        if self.head0 is not None and self.PYP_SPEC is not None:
            spectrograph = load_spectrograph(self.PYP_SPEC)
            subheader = spectrograph.subheader_for_spec(self.head0, self.head0)
        else:
            subheader = {}
        # Add em in
        for key in subheader:
            primary_hdr[key] = subheader[key]
        # Do it
        super(DataCube, self).to_file(ofile, primary_hdr=primary_hdr, hdr=hdr, **kwargs)

    @classmethod
    def from_file(cls, ifile):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.from_file`
        to deal with the header

        Args:
            ifile (str):  Filename holding the object
        """
        with io.fits_open(ifile) as hdu:
            # Read using the base class
            self = super().from_hdu(hdu)
            # Internals
            self.filename = ifile
            self.head0 = hdu[1].header  # Actually use the first extension here, since it contains the WCS
            # Meta
            self.spectrograph = load_spectrograph(self.PYP_SPEC)
            self.spect_meta = self.spectrograph.parse_spec_header(hdu[0].header)
        return self

    @property
    def ivar(self):
        return utils.inverse(self.sig**2)

    @property
    def wcs(self):
        return wcs.WCS(self.head0)


def dar_fitfunc(radec, coord_ra, coord_dec, datfit, wave, obstime, location, pressure, temperature, rel_humidity):
    """ Generates a fitting function to calculate the offset due to differential atmospheric refraction

    Args:
        radec (tuple):
            A tuple containing two floats representing the shift in ra and dec due to DAR.
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
        chisq (float):
            chi-squared difference between datfit and model
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
        Relative RA shift at each wavelength given by `wave_arr`
    dec_diff : `numpy.ndarray`_
        Relative DEC shift at each wavelength given by `wave_arr`
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
    """ Using spline representations of the blaze profile, calculate the grating correction
    that should be applied to the current spectrum (suffix 'curr') relative to the reference
    spectrum (suffix 'ref'). The grating correction is then evaluated at the wavelength
    array given by 'wave_eval'.

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
        grat_corr (`numpy.ndarray`_): The grating correction to apply
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
    """ Fit a 2D Gaussian function to a datacube. This function assumes that each wavelength
    slice of the datacube is well-fit by a 2D Gaussian. The centre of the Gaussian is allowed
    to vary linearly as a function of wavelength.

    NOTE : the integrated flux does not vary with wavelength.

    Args:
        tup (:obj:`tuple`):
            A three element tuple containing the x, y, and z locations of each pixel in the cube
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
        gtwod (`numpy.ndarray`_): The 2D Gaussian evaluated at the coordinate (x, y, z)
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
    """ Fit a 2D Gaussian function to an image.

    Args:
        tup (:obj:`tuple`):
            A two element tuple containing the x and y coordinates of each pixel in the image
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
        gtwod (`numpy.ndarray`_): The 2D Gaussian evaluated at the coordinate (x, y)
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


def rebinND(img, shape):
    """
    Rebin a 2D image to a smaller shape. For example, if img.shape=(100,100),
    then shape=(10,10) would take the mean of the first 10x10 pixels into a
    single output pixel, then the mean of the next 10x10 pixels will be output
    into the next pixel

    Args:
        img (`numpy.ndarray`_):
            A 2D input image
        shape (:obj:`tuple`):
            The desired shape to be returned. The elements of img.shape
            should be an integer multiple of the elements of shape.

    Returns:
        img_out (`numpy.ndarray`_): The input image rebinned to shape
    """
    # Convert input 2D image into a 4D array to make the rebinning easier
    sh = shape[0], img.shape[0]//shape[0], shape[1], img.shape[1]//shape[1]
    # Rebin to the 4D array and then average over the second and last elements.
    img_out = img.reshape(sh).mean(-1).mean(1)
    return img_out


def extract_standard_spec(stdcube, subpixel=20, method='boxcar'):
    """ Extract a spectrum of a standard star from a datacube

    Args:
        std_cube (`astropy.io.fits.HDUList`_):
            An HDU list of fits files
        subpixel (int):
            Number of pixels to subpixelate spectrum when creating mask
        method (str):
            Method used to extract standard star spectrum. Currently, only 'boxcar' is supported

    Returns:
        wave (`numpy.ndarray`_): Wavelength of the star.
        Nlam_star (`numpy.ndarray`_): counts/second/Angstrom
        Nlam_ivar_star (`numpy.ndarray`_): inverse variance of Nlam_star
        gpm_star (`numpy.ndarray`_): good pixel mask for Nlam_star
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
    mask = rebinND(mask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)

    # Generate a sky mask
    newshape = (flxcube.shape[0] * subpixel, flxcube.shape[1] * subpixel)
    smask = np.zeros(newshape)
    nsig = 8  # 8 sigma should be far enough
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < nsig * wid))
    smask[ww] = 1
    smask = rebinND(smask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)
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
        sobj = specobj.SpecObj('IFU', 'DET01', SLITID=0)
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
    """ Mask the spectral edges of each slit (i.e. the pixels near the ends of
    the detector in the spectral direction). Some extreme values of the tilts are
    only sampled with a small fraction of the pixels of the slit width. This leads
    to a bad extrapolation/determination of the sky model.

    Args:
        slitimg (`numpy.ndarray`_):
            An image of the slit indicating which slit each pixel belongs to
        tilts (`numpy.ndarray`_):
            Spectral tilts.

    Returns:
        gpm (`numpy.ndarray`_): A mask of the good sky pixels (True = good)
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
    """ Determine which pixels are included within the specified wavelength range
    Args:
        all_wave (`numpy.ndarray`_):
            The wavelength of each individual pixel
        min_wl (float):
            Minimum wavelength to consider
        max_wl (float):
            Maximum wavelength to consider

    Returns:
        ww (`numpy.ndarray`_): The indices of all_wave that contain pixels within the requested wavelength range
        wavediff (float): The wavelength range (i.e. maximum wavelength - minimum wavelength)
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
    """ Get the wavelength range to use for the white light images
    Args:
        wavemin (float):
            Automatically determined minimum wavelength to use for making the white light image.
        wavemax (float):
            Automatically determined maximum wavelength to use for making the white light image.
        wl_range (list):
            Two element list containing the user-specified values to manually override the automated
            values determined by PypeIt.

    Returns:
        wlrng (list): A two element list containing the minimum and maximum wavelength
        to use for the white light images
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
    """ Generate a white light image using an input cube.

    Args:
        cube (`numpy.ndarray`_):
            3D datacube (the final element contains the wavelength dimension)
        wave (`numpy.ndarray`_, optional):
            1D wavelength array. Only required if wavemin or wavemax are not None.
        wavemin (float, None):
            Minimum wavelength (same units as wave) to be included in the whitelight image.
            You must provide wave as well if you want to reduce the wavelength range.
        wavemax (float, None):
            Maximum wavelength (same units as wave) to be included in the whitelight image.
            You must provide wave as well if you want to reduce the wavelength range.

    Returns:
        wl_img (`numpy.ndarray`_): Whitelight image of the input cube.
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
    """ Load an image and return the image and the associated WCS.

    Args:
        filename (str):
            A fits filename of an image to be used when generating white light
            images. Note, the fits file must have a valid 3D WCS.
        ext (bool, optional):
            The extension that contains the image and WCS

    Returns:
        image (`numpy.ndarray`_): 2D image data
        imgwcs (`astropy.wcs.wcs.WCS`_): The WCS of the image
    """
    imghdu = fits.open(filename)
    image = imghdu[ext].data.T
    imgwcs = wcs.WCS(imghdu[ext].header)
    # Return required info
    return image, imgwcs


def make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                               all_ivar=None, whitelightWCS=None, numra=None, numdec=None, trim=1):
    """ Generate a whitelight image using the individual pixels of every input frame

    Args:
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the RA values of each pixel from all spec2d files
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the DEC values of each pixel from all spec2d files
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength values of each pixel from all spec2d files
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights attributed to each pixel from all spec2d files
        all_idx (`numpy.ndarray`_):
            1D flattened array containing an integer identifier indicating which spec2d file
            each pixel originates from. For example, a 0 would indicate that a pixel originates
            from the first spec2d frame listed in the input file. a 1 would indicate that this
            pixel originates from the second spec2d file, and so forth.
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        all_ivar (`numpy.ndarray`_, optional):
            1D flattened array containing of the inverse variance of each pixel from all spec2d files.
            If provided, inverse variance images will be calculated and returned for each white light image.
        whitelightWCS (`astropy.wcs.wcs.WCS`_, optional):
            The WCS of a reference white light image. If supplied, you must also
            supply numra and numdec.
        numra (int, optional):
            Number of RA spaxels in the reference white light image
        numdec (int, optional):
            Number of DEC spaxels in the reference white light image
        trim (int, optional):
            Number of pixels to grow around a masked region

    Returns:
        tuple : two 3D arrays will be returned, each of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light images.
        The first array is a white light image, and the second array is the corresponding
        inverse variance image. If all_ivar is None, this will be an empty array.
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


def create_wcs(cubepar, all_ra, all_dec, all_wave, dspat, dwv, collapse=False, equinox=2000.0, specname="PYP_SPEC"):
    """
    Create a WCS and the expected edges of the voxels, based on user-specified parameters
    or the extremities of the data.


    Args:
        cubepar (:class:`~pypeit.par.pypeitpar.CubePar`):
            An instance of the CubePar parameter set, contained parameters of the datacube reduction
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the RA values of each pixel from all spec2d files
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the DEC values of each pixel from all spec2d files
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength values of each pixel from all spec2d files
        dspat (float):
            Spatial size of each square voxel (in arcsec). The default is to use the values in cubepar.
        dwv (float):
            Linear wavelength step of each voxel (in Angstroms)
        collapse (bool, optional):
            If True, the spectral dimension will be collapsed to a single channel (primarily for white light images)
        equinox (float, optional):
            Equinox of the WCS
        specname (str, optional):
            Name of the spectrograph

    Returns:
        cubewcs (`astropy.wcs.wcs.WCS`_):
            astropy WCS to be used for the combined cube
        voxedges (tuple) :
            A three element tuple containing the bin edges in the x, y (spatial) and z (wavelength) dimensions
        reference_image (`numpy.ndarray`_, None):
            The reference image to be used for the cross-correlation
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
        `astropy.wcs.wcs.WCS`_ : astropy WCS to be used for the combined cube
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
    """ Calculate wavelength dependent optimal weights. The weighting
        is currently based on a relative (S/N)^2 at each wavelength

    Args:
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the RA values of each pixel from all spec2d files
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the DEC values of each pixel from all spec2d files
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength values of each pixel from all spec2d files
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_idx (`numpy.ndarray`_):
            1D flattened array containing an integer identifier indicating which spec2d file
            each pixel originates from. For example, a 0 would indicate that a pixel originates
            from the first spec2d frame listed in the input file. a 1 would indicate that this
            pixel originates from the second spec2d file, and so forth.
        whitelight_img (`numpy.ndarray`_):
            A 2D array containing a whitelight image, that was created with the input all_* arrays.
        dspat (float):
            The size of each spaxel on the sky (in degrees)
        dwv (float):
            The size of each wavelength pixel (in Angstroms)
        sn_smooth_npix (float, optional):
            Number of pixels used for determining smoothly varying S/N ratio weights.
            This is currently not required, since a relative weighting scheme with a
            polynomial fit is used to calculate the S/N weights.
        relative_weights (bool, optional):
            Calculate weights by fitting to the ratio of spectra?
    Returns:
        `numpy.ndarray`_ : a 1D array the same size as all_sci, containing relative wavelength
                           dependent weights of each input pixel.
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
    rms_sn, weights = coadd.sn_weights(wave_spec, flux_stack, ivar_stack, mask_stack, sn_smooth_npix,
                                       relative_weights=relative_weights)

    # Because we pass back a weights array, we need to interpolate to assign each detector pixel a weight
    all_wghts = np.ones(all_idx.size)
    for ff in range(numfiles):
        ww = (all_idx == ff)
        all_wghts[ww] = interp1d(wave_spec, weights[:, ff], kind='cubic',
                                 bounds_error=False, fill_value="extrapolate")(all_wave[ww])
    msgs.info("Optimal weighting complete")
    return all_wghts


def generate_image_subpixel(image_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar, all_wghts, all_spatpos, all_specpos,
                            all_spatid, tilts, slits, astrom_trans, bins, all_idx=None,
                            spec_subpixel=10, spat_subpixel=10, combine=False):
    """
    Generate a white light image from the input pixels

    Args:
        image_wcs (`astropy.wcs.wcs.WCS`_):
            World coordinate system to use for the white light images.
        all_ra (`numpy.ndarray`_)
            1D flattened array containing the right ascension of each pixel (units = degrees)
        all_dec (`numpy.ndarray`_)
            1D flattened array containing the declination of each pixel (units = degrees)
        all_wave (`numpy.ndarray`_)
            1D flattened array containing the wavelength of each pixel (units = Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used in the combination
        all_spatpos (`numpy.ndarray`_)
            1D flattened array containing the detector pixel location in the spatial direction
        all_specpos (`numpy.ndarray`_)
            1D flattened array containing the detector pixel location in the spectral direction
        all_spatid (`numpy.ndarray`_)
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list)
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`pypeit.slittrace.SlitTraceSet`_, list)
            Information stored about the slits, or a list of SlitTraceSet (see all_idx)
        astrom_trans (:class:`pypeit.alignframe.AlignmentSplines`_, list):
            A Class containing the transformation between detector pixel coordinates
            and WCS pixel coordinates, or a list of Alignment Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional)
            If tilts, slits, and astrom_trans are lists, this should contain a 1D flattened array, of
            the same length as all_sci, containing the index the tilts, slits, and astrom_trans lists
            that corresponds to each pixel. Note that, in this case all of these lists need to be the same length.
        spec_subpixel (`int`, optional):
            What is the subpixellation factor in the spectral direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral direction.
        spat_subpixel (`int`, optional):
            What is the subpixellation factor in the spatial direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial direction.
        combine (`bool`, optional):
            If True, all of the input frames will be combined into a single output. Otherwise, individual images
            will be generated.
    Returns:
        all_wl_imgs (`numpy.ndarray`_): The white light images for all frames
    """
    # Perform some checks on the input -- note, more complete checks are performed in subpixellate()
    _all_idx = np.zeros(all_sci.size) if all_idx is None else all_idx
    numfr = 1 if combine else np.unique(_all_idx).size
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
            ww = np.where(_all_idx >= 0)
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
    """
    Save a datacube using the subpixel algorithm. Refer to the subpixellate() docstring
    for further details about this algorithm

    Args:
        outfile (str):
            Filename to be used to save the datacube
        output_wcs (`astropy.wcs.wcs.WCS`_):
            Output world coordinate system.
        all_ra (`numpy.ndarray`_):
            1D flattened array containing the right ascension of each pixel (units = degrees)
        all_dec (`numpy.ndarray`_):
            1D flattened array containing the declination of each pixel (units = degrees)
        all_wave (`numpy.ndarray`_):
            1D flattened array containing the wavelength of each pixel (units = Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used in the combination
        all_spatpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the spatial direction
        all_specpos (`numpy.ndarray`_):
            1D flattened array containing the detector pixel location in the spectral direction
        all_spatid (`numpy.ndarray`_):
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list):
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`pypeit.slittrace.SlitTraceSet`, list):
            Information stored about the slits, or a list of SlitTraceSet (see all_idx)
        astrom_trans (:class:`pypeit.alignframe.AlignmentSplines`, list):
            A Class containing the transformation between detector pixel coordinates
            and WCS pixel coordinates, or a list of Alignment Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional)
            If tilts, slits, and astrom_trans are lists, this should contain a 1D flattened array, of
            the same length as all_sci, containing the index the tilts, slits, and astrom_trans lists
            that corresponds to each pixel. Note that, in this case all of these lists need to be the same length.
        spec_subpixel (int, optional):
            What is the subpixellation factor in the spectral direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral direction.
        spat_subpixel (int, optional):
            What is the subpixellation factor in the spatial direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial direction.
        overwrite (bool, optional):
            If True, the output cube will be overwritten.
        blaze_wave (`numpy.ndarray`_, optional):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_, optional):
            Spectral blaze function
        fluxcal (bool, optional):
            Are the data flux calibrated? If True, the units are: erg/s/cm^2/Angstrom/arcsec^2
            multiplied by the PYPEIT_FLUX_SCALE. Otherwise, the units are: counts/s/Angstrom/arcsec^2")
        sensfunc (`numpy.ndarray`_, None, optional):
            Sensitivity function that has been applied to the datacube
        whitelight_range (None, list, optional):
            A two element list that specifies the minimum and maximum wavelengths (in Angstroms) to use
            when constructing the white light image (format is: [min_wave, max_wave]). If None, the cube
            will be collapsed over the full wavelength range. If a list is provided an either element of
            the list is None, then the minimum/maximum wavelength range of that element will be set by
            the minimum/maximum wavelength of all_wave.
        specname (str, optional):
            Name of the spectrograph
        debug (bool, optional):
            If True, a residuals cube will be output. If the datacube generation is correct, the
            distribution of pixels in the residual cube with no flux should have mean=0 and std=1.
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
    """
    Subpixellate the input data into a datacube. This algorithm splits
    each detector pixel into multiple subpixels, and then assigns each
    subpixel to a voxel. For example, if spec_subpixel = spat_subpixel = 10,
    then each detector pixel is divided into 10^2=100 subpixels. Alternatively,
    when spec_subpixel = spat_subpixel = 1, this corresponds to the nearest
    grid point (NGP) algorithm.

    Important Note: If spec_subpixel > 1 or spat_subpixel > 1, the errors will
    be correlated, and the covariance is not being tracked, so the errors will
    not be (quite) right. There is a tradeoff one has to make between sampling
    and better looking cubes, versus no sampling and better behaved errors.

    Args:
        output_wcs (`astropy.wcs.wcs.WCS`_):
            Output world coordinate system.
        all_ra (`numpy.ndarray`_)
            1D flattened array containing the right ascension of each pixel (units = degrees)
        all_dec (`numpy.ndarray`_)
            1D flattened array containing the declination of each pixel (units = degrees)
        all_wave (`numpy.ndarray`_)
            1D flattened array containing the wavelength of each pixel (units = Angstroms)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used in the combination
        all_spatpos (`numpy.ndarray`_)
            1D flattened array containing the detector pixel location in the spatial direction
        all_specpos (`numpy.ndarray`_)
            1D flattened array containing the detector pixel location in the spectral direction
        all_spatid (`numpy.ndarray`_)
            1D flattened array containing the spatid of each pixel
        tilts (`numpy.ndarray`_, list)
            2D wavelength tilts frame, or a list of tilt frames (see all_idx)
        slits (:class:`pypeit.slittrace.SlitTraceSet`_, list)
            Information stored about the slits, or a list of SlitTraceSet (see all_idx)
        astrom_trans (:class:`pypeit.alignframe.AlignmentSplines`_, list):
            A Class containing the transformation between detector pixel coordinates
            and WCS pixel coordinates, or a list of Alignment Splines (see all_idx)
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial and z wavelength coordinates
        all_idx (`numpy.ndarray`_, optional)
            If tilts, slits, and astrom_trans are lists, this should contain a 1D flattened array, of
            the same length as all_sci, containing the index the tilts, slits, and astrom_trans lists
            that corresponds to each pixel. Note that, in this case all of these lists need to be the same length.
        spec_subpixel (`int`, optional):
            What is the subpixellation factor in the spectral direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spectral direction.
        spat_subpixel (`int`, optional):
            What is the subpixellation factor in the spatial direction. Higher values give more reliable results,
            but note that the time required goes as (spec_subpixel * spat_subpixel). The default value is 5,
            which divides each detector pixel into 5 subpixels in the spatial direction.
        debug (bool):
            If True, a residuals cube will be output. If the datacube generation is correct, the
            distribution of pixels in the residual cube with no flux should have mean=0 and std=1.

    Returns:
        datacube (`numpy.ndarray`_): The datacube generated from the subpixellated inputs
        varcube (`numpy.ndarray`_): The corresponding variance cube
        bpmcube (`numpy.ndarray`_): The corresponding bad pixel mask cube
        residcube (`numpy.ndarray`_): If debug=True, the resid cube is also returned.
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
        outfile (str): The output filename to use.
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
    """ Given the output filename of a datacube, create an appropriate whitelight fits file name

    Args:
        outfile (str):
            The output filename used for the datacube.

    Returns:
        out_wl_filename (str): The output filename to use for the whitelight image.
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
        align = False
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
        # TODO: Check that the file exists?
        if cubepar['grating_corr'] and flatfile not in flat_splines.keys():
            msgs.info("Calculating relative sensitivity for grating correction")
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
                msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f}, {2:+0.3f}".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
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
        # TODO :: THIS IS JUST TEMPORARY (probably)... the first bit outputs separate files, the second bit combines all frames.
        #  Might want to have an option where if reference_image is provided, but not combine, the first option is done.
        if combine:
            generate_cube_subpixel(outfile, cube_wcs, all_ra, all_dec, all_wave, all_sci, all_ivar,
                                   np.ones(all_wghts.size),  # all_wghts,
                                   all_spatpos, all_specpos, all_spatid, all_tilts, all_slits, all_align, vox_edges,
                                   all_idx=all_idx, overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                                   fluxcal=fluxcal, sensfunc=sensfunc, specname=specname, whitelight_range=wl_wvrng,
                                   spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
        else:
            # embed()
            # assert(False)
            for ff in range(numfiles):
                outfile = get_output_filename("", cubepar['output_filename'], False, ff)
                ww = np.where(all_idx == ff)
                generate_cube_subpixel(outfile, cube_wcs, all_ra[ww], all_dec[ww], all_wave[ww], all_sci[ww], all_ivar[ww], np.ones(all_wghts[ww].size),
                                       all_spatpos[ww], all_specpos[ww], all_spatid[ww], all_tilts[ff], all_slits[ff], all_align[ff], vox_edges,
                                       all_idx=all_idx[ww], overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                                       fluxcal=fluxcal, sensfunc=sensfunc, specname=specname, whitelight_range=wl_wvrng,
                                       spec_subpixel=spec_subpixel, spat_subpixel=spat_subpixel)
