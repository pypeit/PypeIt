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
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit import alignframe, datamodel, flatfield, io, masterframe, specobj, spec2dobj, utils
from pypeit.core.flux_calib import load_extinction_data, extinction_correction, fit_zeropoint, get_standard_spectrum, ZP_UNIT_CONST, PYPEIT_FLUX_SCALE
from pypeit.core.flexure import calculate_image_offset
from pypeit.core import coadd, extract, findobj_skymask, parse, skysub
from pypeit.core.procimg import grow_mask
from pypeit.spectrographs.util import load_spectrograph

from IPython import embed


class DataCube(datamodel.DataContainer):
    """
    DataContainer to hold the products of a datacube

    See the datamodel for argument descriptions

    Args:
        flux (`numpy.ndarray`_):
            The science datacube (nwave, nspaxel_y, nspaxel_x)
        variance (`numpy.ndarray`_):
            The variance datacube (nwave, nspaxel_y, nspaxel_x)
        blaze_wave (`numpy.ndarray`_):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_):
            The spectral blaze function
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
    version = '1.0.2'

    datamodel = {'flux': dict(otype=np.ndarray, atype=np.floating, descr='Flux array in units of counts/s/Ang/arcsec^2'
                                                                         'or 10^-17 erg/s/cm^2/Ang/arcsec^2'),
                 'variance': dict(otype=np.ndarray, atype=np.floating, descr='Variance array (matches units of flux)'),
                 'blaze_wave': dict(otype=np.ndarray, atype=np.floating, descr='Wavelength array of the spectral blaze function'),
                 'blaze_spec': dict(otype=np.ndarray, atype=np.floating, descr='The spectral blaze function'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt: Spectrograph name'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the datacube is fluxed.')}

    @classmethod
    def from_file(cls, ifile):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.from_file`
        to deal with the header

        Args:
            ifile (str):  Filename holding the object
        """
        hdul = fits.open(ifile)
        slf = super(DataCube, cls).from_hdu(hdul)

        # Internals
        slf.filename = ifile
        slf.head0 = hdul[0].header
        # Meta
        slf.spectrograph = load_spectrograph(slf.PYP_SPEC)
        slf.spect_meta = slf.spectrograph.parse_spec_header(slf.head0)
        return slf

    def __init__(self, flux, variance, PYP_SPEC, blaze_wave, blaze_spec, fluxed=None):

        args, _, _, values = inspect.getargvalues(inspect.currentframe())
        _d = dict([(k,values[k]) for k in args[1:]])
        # Setup the DataContainer
        datamodel.DataContainer.__init__(self, d=_d)

    def _init_internals(self):
        self.head0 = None
        self.filename = None
        self.spectrograph = None
        self.spect_meta = None

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
            primary_hdr = io.initialize_header(primary=True)
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
    # Generate the coordinate with atmopheric conditions
    coord_atmo = SkyCoord(coord_ra + diff_ra, coord_dec + diff_dec, unit=(units.deg, units.deg))
    coord_altaz = coord_atmo.transform_to(AltAz(obstime=obstime, location=location, obswl=wave,
                                          pressure=pressure, temperature=temperature,
                                          relative_humidity=rel_humidity))
    # Return chi-squared value
    return np.sum((np.array([coord_altaz.alt.value, coord_altaz.az.value])-datfit)**2)


def dar_correction(wave_arr, coord, obstime, location, pressure, temperature, rel_humidity,
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
    idx_max = np.unravel_index(np.argmax(image), image.shape)
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


def extract_standard_spec(stdcube, subsample=20, method='boxcar'):
    """ Extract a spectrum of a standard star from a datacube

    Args:
        std_cube (`astropy.io.fits.HDUList`_):
            An HDU list of fits files
        subsample (int):
            Number of pixels to subpixellate spectrum when creating mask
        method (str):
            Method used to extract standard star spectrum

    Returns:
        wave (`numpy.ndarray`_): Wavelength of the star.
        Nlam_star (`numpy.ndarray`_): counts/second/Angstrom
        Nlam_ivar_star (`numpy.ndarray`_): inverse variance of Nlam_star
        gpm_star (`numpy.ndarray`_): good pixel mask for Nlam_star
    """
    # Extract some information from the HDU list
    flxcube = stdcube['FLUX'].data.T.copy()
    varcube = stdcube['VARIANCE'].data.T.copy()
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
    x = np.linspace(0, flxcube.shape[0] - 1, flxcube.shape[0] * subsample)
    y = np.linspace(0, flxcube.shape[1] - 1, flxcube.shape[1] * subsample)
    xx, yy = np.meshgrid(x, y, indexing='ij')

    # Generate a mask
    newshape = (flxcube.shape[0] * subsample, flxcube.shape[1] * subsample)
    mask = np.zeros(newshape)
    nsig = 4  # 4 sigma should be far enough... Note: percentage enclosed for 2D Gaussian = 1-np.exp(-0.5 * nsig**2)
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < nsig * wid))
    mask[ww] = 1
    mask = rebinND(mask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)

    # Generate a sky mask
    newshape = (flxcube.shape[0] * subsample, flxcube.shape[1] * subsample)
    smask = np.zeros(newshape)
    nsig = 8  # 8 sigma should be far enough
    ww = np.where((np.sqrt((xx - popt[1]) ** 2 + (yy - popt[2]) ** 2) < nsig * wid))
    smask[ww] = 1
    smask = rebinND(smask, (flxcube.shape[0], flxcube.shape[1])).reshape(flxcube.shape[0], flxcube.shape[1], 1)
    smask -= mask

    # Subtract the residual sky
    skymask = (varcube > 0.0) * smask
    skycube = flxcube * skymask
    skyspec = skycube.sum(0).sum(0)
    nrmsky = skymask.sum(0).sum(0)
    skyspec *= utils.inverse(nrmsky)
    flxcube -= skyspec.reshape((1, 1, numwave))

    if method == 'boxcar':
        msgs.info("Extracting a boxcar spectrum of datacube")
        # Extract boxcar
        cntmask = (varcube > 0.0) * mask
        scimask = flxcube * cntmask
        varmask = varcube * cntmask**2
        cnt_spec = cntmask.sum(0).sum(0) * utils.inverse(mask.sum())
        nrmcnt = utils.inverse(cnt_spec)
        box_flux = scimask.sum(0).sum(0) * nrmcnt
        box_var = varmask.sum(0).sum(0) * nrmcnt**2
        box_gpm = np.ones(box_flux.size, dtype=np.bool)
        # Setup the return values
        ret_flux, ret_var, ret_gpm = box_flux, box_var, box_gpm
    elif method == 'gauss2d':
        # Generate a mask
        fitmask = (varcube > 0.0) * mask
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
        # First do a boxcar along one dimension
        msgs.info("Collapsing datacube to a 2D image")
        omask = mask+smask
        idx_sum = 0
        cntmask = (varcube > 0.0) * omask
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
        thismask = np.ones(box_sciimg.shape, dtype=np.bool)
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


def make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat, ref_filename):
    """ Generate a whitelight image using the individual pixels of every
    input frame, based on a reference image. Note the, the reference
    image must have a well-defined WCS.

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
        ref_filename (str):
            A fits filename of a reference image to be used when generating white light
            images. Note, the fits file must have a valid 3D WCS.

    Returns:
        tuple : two `numpy.ndarray`_ and one WCS will be returned. The first is a 2D reference image
        loaded from ref_filename. The second element is a 3D array of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light images. The third is
        the WCS of the white light image.
    """
    refhdu = fits.open(ref_filename)
    reference_image = refhdu[0].data.T[:, :, 0]
    refwcs = wcs.WCS(refhdu[0].header)
    numra, numdec = reference_image.shape
    # Generate coordinate system (i.e. update wavelength range to include all values)
    coord_min = refwcs.wcs.crval
    coord_dlt = refwcs.wcs.cdelt
    coord_min[2] = np.min(all_wave)
    coord_dlt[2] = np.max(all_wave) - np.min(all_wave)  # For white light, we want to bin all wavelength pixels
    wlwcs = generate_masterWCS(coord_min, coord_dlt)

    # Generate white light images
    whitelight_imgs, _, _ = make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                                                       whitelightWCS=wlwcs, numra=numra, numdec=numdec)
    # Return required info
    return reference_image, whitelight_imgs, wlwcs


def make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                               all_ivar=None, whitelightWCS=None, numra=None, numdec=None):
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
            Inverse variance of each pixel from all spec2d files. If provided,
            inverse variance images will be calculated and return for each white light image.
        whitelightWCS (`astropy.wcs.wcs.WCS`_, optional):
            The WCS of a reference white light image. If supplied, you must also
            supply numra and numdec.
        numra (int, optional):
            Number of RA spaxels in the reference white light image
        numdec (int, optional):
            Number of DEC spaxels in the reference white light image

    Returns:
        tuple : two 3D arrays will be returned, each of shape [N, M, numfiles],
        where N and M are the spatial dimensions of the combined white light images.
        The first array is a white light image, and the second array is the corresponding
        inverse variance image. If all_ivar is None, this will be an empty array.
    """
    # Determine number of files
    numfiles = np.unique(all_idx).size

    if whitelightWCS is None:
        # Generate a master 2D WCS to register all frames
        coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
        coord_dlt = [dspat, dspat, np.max(all_wave) - np.min(all_wave)]
        whitelightWCS = generate_masterWCS(coord_min, coord_dlt)

        # Generate coordinates
        cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)
        numra = int((np.max(all_ra) - np.min(all_ra)) * cosdec / dspat)
        numdec = int((np.max(all_dec) - np.min(all_dec)) / dspat)
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
    trim = 3
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
#        gpm = grow_masked(whtlght == 0, trim, 1) == 0  # A good pixel = 1
        # TODO: NEED TO CHECK THIS IS OKAY!!
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


def generate_masterWCS(crval, cdelt, equinox=2000.0, name="Instrument Unknown"):
    """
    Generate a WCS that will cover all input spec2D files

    Args:
        crval (list):
            3 element list containing the [RA, DEC, WAVELENGTH] of
            the reference pixel
        cdelt (list):
            3 element list containing the delta values of the [RA,
            DEC, WAVELENGTH]
        equinox (float):
            Equinox of the WCS

    Returns:
        `astropy.wcs.wcs.WCS`_ : astropy WCS to be used for the combined cube
    """
    # Create a new WCS object.
    msgs.info("Generating Master WCS")
    w = wcs.WCS(naxis=3)
    w.wcs.equinox = equinox
    w.wcs.name = name
    w.wcs.radesys = 'FK5'
    # Insert the coordinate frame
    w.wcs.cname = ['RA', 'DEC', 'Wavelength']
    w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
    w.wcs.ctype = ["RA---TAN", "DEC--TAN", "AWAV"]
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

    # Generate a master 2D WCS to register all frames
    coord_min = [np.min(all_ra), np.min(all_dec), np.min(all_wave)]
    coord_dlt = [dspat, dspat, dwv]
    whitelightWCS = generate_masterWCS(coord_min, coord_dlt)
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


def generate_cube_ngp(outfile, hdr, all_sci, all_ivar, all_wghts, pix_coord, bins,
                      overwrite=False, blaze_wave=None, blaze_spec=None, fluxcal=False,
                      specname="PYP_SPEC", debug=False):
    """
    Save a datacube using the Nearest Grid Point (NGP) algorithm.

    Args:
        outfile (`str`):
            Filename to be used to save the datacube
        hdr (`astropy.io.fits.header_`):
            Header of the output datacube (must contain WCS)
        all_sci (`numpy.ndarray`_):
            1D flattened array containing the counts of each pixel from all spec2d files
        all_ivar (`numpy.ndarray`_):
            1D flattened array containing the inverse variance of each pixel from all spec2d files
        all_wghts (`numpy.ndarray`_):
            1D flattened array containing the weights of each pixel to be used in the combination
        pix_coord (`numpy.ndarray`_):
            The NGP pixel coordinates corresponding to the RA,DEC,WAVELENGTH of each individual
            pixel in the processed spec2d frames. After setting up an astropy WCS, pix_coord is
            returned by the function: `astropy.wcs.WCS.wcs_world2pix_`
        bins (tuple):
            A 3-tuple (x,y,z) containing the histogram bin edges in x,y spatial and z wavelength coordinates    :param overwrite:
        blaze_wave (`numpy.ndarray`_):
            Wavelength array of the spectral blaze function
        blaze_spec (`numpy.ndarray`_):
            Spectral blaze function
        fluxcal (bool):
            Are the data flux calibrated?
        specname (str):
            Name of the spectrograph
        debug (bool):
            Debug the code by writing out a residuals cube?
    """
    # Add the unit of flux to the header
    if fluxcal:
        hdr['FLUXUNIT'] = (PYPEIT_FLUX_SCALE, "Flux units -- erg/s/cm^2/Angstrom/arcsec^2")
    else:
        hdr['FLUXUNIT'] = (1, "Flux units -- counts/s/Angstrom/arcsec^2")

    # Use NGP to generate the cube - this ensures errors between neighbouring voxels are not correlated
    datacube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci * all_wghts)
    norm, edges = np.histogramdd(pix_coord, bins=bins, weights=all_wghts)
    norm_cube = utils.inverse(norm)
    datacube *= norm_cube
    # Create the variance cube, including weights
    msgs.info("Generating variance cube")
    all_var = utils.inverse(all_ivar)
    var_cube, edges = np.histogramdd(pix_coord, bins=bins, weights=all_var * all_wghts**2)
    var_cube *= norm_cube**2

    # Save the datacube
    if debug:
        datacube_resid, edges = np.histogramdd(pix_coord, bins=bins, weights=all_sci*np.sqrt(all_ivar))
        norm, edges = np.histogramdd(pix_coord, bins=bins)
        norm_cube = utils.inverse(norm)
        outfile_resid = "datacube_resid.fits"
        msgs.info("Saving datacube as: {0:s}".format(outfile_resid))
        hdu = fits.PrimaryHDU((datacube_resid*norm_cube).T, header=hdr)
        hdu.writeto(outfile_resid, overwrite=overwrite)

    msgs.info("Saving datacube as: {0:s}".format(outfile))
    final_cube = DataCube(datacube.T, var_cube.T, specname, blaze_wave, blaze_spec, fluxed=fluxcal)
    final_cube.to_file(outfile, hdr=hdr, overwrite=overwrite)


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


def coadd_cube(files, spectrograph=None, parset=None, overwrite=False):
    """ Main routine to coadd spec2D files into a 3D datacube

    Args:
        files (:obj:`list`):
            List of all spec2D files
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
        parset = spec.default_pypeit_par()
    cubepar = parset['reduce']['cube']
    flatpar = parset['calibrations']['flatfield']

    # Get the detector number and string representation
    det = 1 if parset['rdx']['detnum'] is None else parset['rdx']['detnum']
    detname = spec.get_det_name(det)

    # prep
    numfiles = len(files)
    combine = cubepar['combine']

    # Check if the output file exists
    if combine:
        outfile = get_output_filename("", cubepar['output_filename'], combine)
        out_whitelight = os.path.splitext(outfile)[0] + "_whitelight.fits"
        if os.path.exists(outfile) and not overwrite:
            msgs.error("Output filename already exists:"+msgs.newline()+outfile)
        if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
            msgs.error("Output filename already exists:"+msgs.newline()+out_whitelight)
    else:
        # Finally, if there's just one file, check if the output filename is given
        if numfiles == 1 and cubepar['output_filename'] != "":
            outfile = get_output_filename("", cubepar['output_filename'], True, -1)
            out_whitelight = os.path.splitext(outfile)[0] + "_whitelight.fits"
            if os.path.exists(outfile) and not overwrite:
                msgs.error("Output filename already exists:" + msgs.newline() + outfile)
            if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
                msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)
        else:
            for ff in range(numfiles):
                outfile = get_output_filename(files[ff], cubepar['output_filename'], combine, ff+1)
                out_whitelight = os.path.splitext(outfile)[0] + "_whitelight.fits"
                if os.path.exists(outfile) and not overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + outfile)
                if os.path.exists(out_whitelight) and cubepar['save_whitelight'] and not overwrite:
                    msgs.error("Output filename already exists:" + msgs.newline() + out_whitelight)

    # Check the reference cube and image exist, if requested
    fluxcal = False
    blaze_wave, blaze_spec = None, None
    blaze_spline, flux_spline = None, None
    if cubepar['standard_cube'] is not None:
        if not os.path.exists(cubepar['standard_cube']):
            msgs.error("Standard cube does not exist:" + msgs.newline() + cubepar['reference_cube'])
        fluxcal = True
        senspar = parset['sensfunc']
        # Load the standard star cube and retrieve its RA + DEC
        stdcube = fits.open(cubepar['standard_cube'])
        star_ra, star_dec = stdcube[1].header['CRVAL1'], stdcube[1].header['CRVAL2']
        # Extract the information about the blaze
        if cubepar['grating_corr']:
            blaze_wave, blaze_spec = stdcube['BLAZE_WAVE'].data, stdcube['BLAZE_SPEC'].data
            blaze_spline = interp1d(blaze_wave, blaze_spec, kind='linear', bounds_error=False, fill_value="extrapolate")
        # Extract a spectrum of the standard star
        wave, Nlam_star, Nlam_ivar_star, gpm_star = extract_standard_spec(stdcube)
        # Read in some information above the standard star
        std_dict = get_standard_spectrum(star_type=senspar['star_type'],
                                         star_mag=senspar['star_mag'],
                                         ra=star_ra, dec=star_dec)
        # Calculate the sensitivity curve
        zeropoint_data, zeropoint_data_gpm, zeropoint_fit, zeropoint_fit_gpm =\
            fit_zeropoint(wave.value, Nlam_star, Nlam_ivar_star, gpm_star, std_dict,
                          mask_abs_lines=senspar['mask_abs_lines'], balm_mask_wid=senspar['UVIS']['balm_mask_wid'],
                          nresln=senspar['UVIS']['nresln'], resolution=senspar['UVIS']['resolution'],
                          trans_thresh=senspar['UVIS']['trans_thresh'], polyorder=senspar['polyorder'],
                          polycorrect=senspar['UVIS']['polycorrect'], polyfunc=senspar['UVIS']['polyfunc'])
        wgd = np.where(zeropoint_fit_gpm)
        sens = np.power(10.0, -0.4 * (zeropoint_fit[wgd] - ZP_UNIT_CONST)) / np.square(wave[wgd])
        flux_spline = interp1d(wave[wgd], sens, kind='linear', bounds_error=False, fill_value="extrapolate")
    if cubepar['reference_image'] is not None:
        if not os.path.exists(cubepar['reference_image']):
            msgs.error("Reference cube does not exist:" + msgs.newline() + cubepar['reference_image'])

    # Initialise arrays for storage
    all_ra, all_dec, all_wave = np.array([]), np.array([]), np.array([])
    all_sci, all_ivar, all_idx, all_wghts = np.array([]), np.array([]), np.array([]), np.array([])
    all_wcs = []
    dspat = None if cubepar['spatial_delta'] is None else cubepar['spatial_delta']/3600.0  # binning size on the sky (/3600 to convert to degrees)
    dwv = cubepar['wave_delta']       # binning size in wavelength direction (in Angstroms)
    wave_ref = None
    weights = np.ones(numfiles)  # Weights to use when combining cubes
    locations = parset['calibrations']['alignment']['locations']
    flat_splines = dict()   # A dictionary containing the splines of the flatfield
    # Load the scaleimg frame for the scale correction
    relScaleImg = np.array([1])
    if cubepar['scale_corr'] is not None:
        msgs.info("Loading scale image for relative spectral illumination correction:" +
                  msgs.newline() + cubepar['scale_corr'])
        try:
            spec2DObj = spec2dobj.Spec2DObj.from_file(cubepar['scale_corr'], detname)
            relScaleImg = spec2DObj.scaleimg
        except:
            msgs.warn("Could not load scaleimg from spec2d file:" + msgs.newline() + cubepar['scale_corr'] +
                      "scale correction will not be performed")
            cubepar['scale_corr'] = None
    # Load all spec2d files and prepare the data for making a datacube
    for ff, fil in enumerate(files):
        # Load it up
        spec2DObj = spec2dobj.Spec2DObj.from_file(fil, detname)
        detector = spec2DObj.detector
        flexure = None  #spec2DObj.sci_spat_flexure

        # Load the header
        hdr = fits.open(fil)[0].header

        # Get the exposure time
        exptime = hdr['EXPTIME']

        # Setup for PypeIt imports
        msgs.reset(verbosity=2)

        # Extract the information
        relscl = spec2DObj.scaleimg/relScaleImg
        sciimg = (spec2DObj.sciimg-spec2DObj.skymodel)*relscl  # Subtract sky
        ivar = spec2DObj.ivarraw / relscl**2
        waveimg = spec2DObj.waveimg
        bpmmask = spec2DObj.bpmmask

        # Grab the slit edges
        slits = spec2DObj.slits

        wave0 = waveimg[waveimg != 0.0].min()
        # Calculate the delta wave in every pixel on the slit
        waveimp = np.roll(waveimg, 1, axis=0)
        waveimn = np.roll(waveimg, -1, axis=0)
        dwaveimg = np.zeros_like(waveimg)
        # All good pixels
        wnz = np.where((waveimg!=0)&(waveimp!=0))
        dwaveimg[wnz] = np.abs(waveimg[wnz]-waveimp[wnz])
        # All bad pixels
        wnz = np.where((waveimg!=0)&(waveimp==0))
        dwaveimg[wnz] = np.abs(waveimg[wnz]-waveimn[wnz])
        # All endpoint pixels
        dwaveimg[0, :] = np.abs(waveimg[0, :] - waveimn[0, :])
        dwaveimg[-1, :] = np.abs(waveimg[-1, :] - waveimp[-1, :])
        dwv = np.median(dwaveimg[dwaveimg != 0.0]) if cubepar['wave_delta'] is None else cubepar['wave_delta']

        msgs.info("Using wavelength solution: wave0={0:.3f}, dispersion={1:.3f} Angstrom/pixel".format(wave0, dwv))

        msgs.info("Constructing slit image")
        slitid_img_init = slits.slit_img(pad=0, initial=True, flexure=flexure)
        onslit_gpm = (slitid_img_init > 0) & (bpmmask == 0)

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
            msgs.warn("Spatial scale requested ({0:f}'') is less than the pixel scale ({1:f}'')".format(dspat,pxscl))
        if pxscl > dspat:
            msgs.warn("Spatial scale requested ({0:f}'') is less than the slicer scale ({1:f}'')".format(dspat, slscl))

        # Loading the alignments frame for these data
        astrometric = cubepar['astrometric']
        msgs.info("Loading alignments")
        alignfile = masterframe.construct_file_name(alignframe.Alignments, hdr['TRACMKEY'], master_dir=hdr['PYPMFDIR'])
        alignments = None
        if os.path.exists(alignfile) and cubepar['astrometric']:
            alignments = alignframe.Alignments.from_file(alignfile)
        else:
            msgs.warn("Could not find Master Alignment frame:"+msgs.newline()+alignfile)
            msgs.warn("Astrometric correction will not be performed")
            astrometric = False

        # Generate an RA/DEC image
        msgs.info("Generating RA/DEC image")
        raimg, decimg, minmax = slits.get_radec_image(frame_wcs, alignments, spec2DObj.tilts, locations,
                                                      astrometric=astrometric, initial=True, flexure=flexure)

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
            ra_corr, dec_corr = dar_correction(waveimg[onslit_gpm], coord, obstime, location,
                                               pressure*units.bar, temperature*units.deg_C, rel_humidity, wave_ref=wave_ref)
            raimg[onslit_gpm] += ra_corr*np.cos(np.mean(decimg[onslit_gpm]) * np.pi / 180.0)
            decimg[onslit_gpm] += dec_corr

        # Get copies of arrays to be saved
        wave_ext = waveimg[onslit_gpm].copy()
        flux_ext = sciimg[onslit_gpm].copy()
        ivar_ext = ivar[onslit_gpm].copy()
        dwav_ext = dwaveimg[onslit_gpm].copy()

        # Correct for sensitivity as a function of grating angle
        # (this assumes the spectrum of the flatfield lamp has the same shape for all setups)
        flatfile = masterframe.construct_file_name(flatfield.FlatImages, hdr['FLATMKEY'], master_dir=hdr['PYPMFDIR'])
        if cubepar['grating_corr'] and flatfile not in flat_splines.keys():
            msgs.info("Calculating relative sensitivity for grating correction")
            flatimages = flatfield.FlatImages.from_file(flatfile)
            flatframe = flatimages.pixelflat_model
            flatframe /= flatimages.fit2illumflat(slits, frametype='pixel', initial=True, flexure_shift=flexure)
            # Calculate the relative scale
            scale_model = flatfield.illum_profile_spectral(flatframe, waveimg, slits,
                                                           slit_illum_ref_idx=flatpar['slit_illum_ref_idx'], model=None,
                                                           skymask=None, trim=flatpar['slit_trim'], flexure=flexure,
                                                           smooth_npix=flatpar['slit_illum_smooth_npix'])
            # Apply the relative scale and generate a 1D "spectrum"
            onslit = waveimg != 0
            wavebins = np.linspace(np.min(waveimg[onslit]), np.max(waveimg[onslit]), slits.nspec)
            hist, edge = np.histogram(waveimg[onslit], bins=wavebins, weights=flatframe[onslit]/scale_model[onslit])
            cntr, edge = np.histogram(waveimg[onslit], bins=wavebins)
            cntr = cntr.astype(np.float)
            norm = (cntr != 0) / (cntr + (cntr == 0))
            spec_spl = hist * norm
            wave_spl = 0.5 * (wavebins[1:] + wavebins[:-1])
            flat_splines[flatfile] = interp1d(wave_spl, spec_spl, kind='linear',
                                              bounds_error=False, fill_value="extrapolate")
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
        extinct = load_extinction_data(longitude, latitude)
        # extinction_correction requires the wavelength is sorted
        wvsrt = np.argsort(wave_ext)
        ext_corr = extinction_correction(wave_ext[wvsrt] * units.AA, airmass, extinct)
        # Grating correction
        grat_corr = 1.0
        if cubepar['grating_corr']:
            msgs.info("Calculating the grating correction")
            grat_corr_tmp = flat_splines[flatfile](wave_ext[wvsrt]) / blaze_spline(wave_ext[wvsrt])
            # Fit a low order polynomial to this correction
            minw, maxw = max(np.min(wave_spl), np.min(blaze_wave)), max(np.min(wave_spl), np.max(blaze_wave))
            wblz = np.where((wave_ext[wvsrt] > minw) & (wave_ext[wvsrt] < maxw))
            wave_corr = (wave_ext[wvsrt] - minw) / (maxw - minw)
            coeff_gratcorr = np.polyfit(wave_corr[wblz], grat_corr_tmp[wblz], 2)
            grat_corr = np.polyval(coeff_gratcorr, wave_corr)
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
        scl_units = dwav_ext[wvsrt] * 3600.0 * 3600.0 * (frame_wcs.wcs.cdelt[0] * frame_wcs.wcs.cdelt[1])
        flux_sav /= scl_units
        ivar_sav *= scl_units ** 2

        # sort back to the original ordering
        resrt = np.argsort(wvsrt)
        numpix = raimg[onslit_gpm].size

        # Calculate the weights relative to the zeroth cube
        # TODO :: This is dodgy when there is mostly sky...
        weights[ff] = np.median(flux_sav[resrt]*np.sqrt(ivar_sav[resrt]))**2

        # If individual frames are to be output, there's no need to store information, just make the cubes now
        if not combine:
            # Get the output filename
            if numfiles == 1 and cubepar['output_filename'] != "":
                outfile = get_output_filename("", cubepar['output_filename'], True, -1)
            else:
                outfile = get_output_filename(fil, cubepar['output_filename'], combine, ff+1)
            # Generate individual whitelight images of each spec2d file
            if cubepar['save_whitelight']:
                out_whitelight = os.path.splitext(outfile)[0] + "_whitelight.fits"
                whitelight_img, _, wlwcs = make_whitelight_frompixels(raimg[onslit_gpm], decimg[onslit_gpm], wave_ext,
                                                                      flux_sav[resrt], np.ones(numpix), np.zeros(numpix), dspat)
                msgs.info("Saving white light image as: {0:s}".format(out_whitelight))
                img_hdu = fits.PrimaryHDU(whitelight_img.T, header=wlwcs.to_header())
                img_hdu.writeto(out_whitelight, overwrite=overwrite)
            # Make the datacube
            slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))
            numwav = int((np.max(waveimg) - wave0) / dwv)
            bins = spec.get_datacube_bins(slitlength, minmax, numwav)
            msgs.info("Generating pixel coordinates")
            pix_coord = frame_wcs.wcs_world2pix(np.vstack((raimg[onslit_gpm], decimg[onslit_gpm], wave_ext * 1.0E-10)).T, 0)
            hdr = frame_wcs.to_header()
            generate_cube_ngp(outfile, hdr, flux_sav[resrt], ivar_sav[resrt], np.ones(numpix), pix_coord, bins,
                              overwrite=overwrite, blaze_wave=blaze_wave, blaze_spec=blaze_spec,
                              fluxcal=fluxcal, specname=specname)
            continue

        # Store the information
        all_ra = np.append(all_ra, raimg[onslit_gpm].copy())
        all_dec = np.append(all_dec, decimg[onslit_gpm].copy())
        all_wave = np.append(all_wave, wave_ext.copy())
        all_sci = np.append(all_sci, flux_sav[resrt].copy())
        all_ivar = np.append(all_ivar, ivar_sav[resrt].copy())
        all_idx = np.append(all_idx, ff*np.ones(numpix))
        all_wghts = np.append(all_wghts, weights[ff]*np.ones(numpix)/weights[0])

    # No need to continue if we are not combining frames
    if not combine:
        return

    # Grab cos(dec) for convenience
    cosdec = np.cos(np.mean(all_dec) * np.pi / 180.0)

    # Register spatial offsets between all frames
    # Check if a reference whitelight image should be used to register the offsets
    if cubepar["reference_image"] is None:
        # Generate white light images
        whitelight_imgs, _, _ = make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat)
        # ref_idx will be the index of the cube with the highest S/N
        ref_idx = np.argmax(weights)
        reference_image = whitelight_imgs[:, :, ref_idx].copy()
        msgs.info("Calculating spatial translation of each cube relative to cube #{0:d})".format(ref_idx+1))
    else:
        ref_idx = -1  # Don't use an index
        # Load reference information
        reference_image, whitelight_imgs, wlwcs = \
            make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                                    cubepar['reference_image'])
        msgs.info("Calculating the spatial translation of each cube relative to user-defined 'reference_image'")

    # Calculate the image offsets - check the reference is a zero shift
    ra_shift_ref, dec_shift_ref = calculate_image_offset(reference_image.copy(), reference_image.copy())
    for ff in range(numfiles):
        # Don't correlate the reference image with itself
        if ff == ref_idx:
            continue
        # Calculate the shift
        ra_shift, dec_shift = calculate_image_offset(whitelight_imgs[:, :, ff], reference_image.copy())
        # Convert to reference
        ra_shift -= ra_shift_ref
        dec_shift -= dec_shift_ref
        # Convert pixel shift to degress shift
        ra_shift *= dspat/cosdec
        dec_shift *= dspat
        msgs.info("Spatial shift of cube #{0:d}: RA, DEC (arcsec) = {1:+0.3f}, {2:+0.3f}".format(ff+1, ra_shift*3600.0, dec_shift*3600.0))
        # Apply the shift
        all_ra[all_idx == ff] += ra_shift
        all_dec[all_idx == ff] += dec_shift

    # Generate a white light image of *all* data
    msgs.info("Generating global white light image")
    if cubepar["reference_image"] is None:
        whitelight_img, _, wlwcs = make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                                              np.zeros(all_ra.size), dspat)
    else:
        _, whitelight_img, wlwcs = \
            make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, np.zeros(all_ra.size),
                                    dspat, cubepar['reference_image'])

    # Calculate the relative spectral weights of all pixels
    all_wghts = compute_weights(all_ra, all_dec, all_wave, all_sci, all_ivar, all_idx,
                                whitelight_img[:, :, 0], dspat, dwv,
                                relative_weights=cubepar['relative_weights'])

    # Check if a whitelight image should be saved
    if cubepar['save_whitelight']:
        # Check if the white light image still needs to be generated - if so, generate it now
        if whitelight_img is None:
            msgs.info("Generating global white light image")
            if cubepar["reference_image"] is None:
                whitelight_img, _, wlwcs = make_whitelight_frompixels(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                                                      np.zeros(all_ra.size), dspat)
            else:
                _, whitelight_img, wlwcs = \
                    make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts,
                                            np.zeros(all_ra.size),
                                            dspat, cubepar['reference_image'])
        # Prepare and save the fits file
        msgs.info("Saving white light image as: {0:s}".format(out_whitelight))
        img_hdu = fits.PrimaryHDU(whitelight_img.T, header=wlwcs.to_header())
        img_hdu.writeto(out_whitelight, overwrite=overwrite)

    # Setup the cube ranges
    ra_min = cubepar['ra_min'] if cubepar['ra_min'] is not None else np.min(all_ra)
    ra_max = cubepar['ra_max'] if cubepar['ra_max'] is not None else np.max(all_ra)
    dec_min = cubepar['dec_min'] if cubepar['dec_min'] is not None else np.min(all_dec)
    dec_max = cubepar['dec_max'] if cubepar['dec_max'] is not None else np.max(all_dec)
    wav_min = cubepar['wave_min'] if cubepar['wave_min'] is not None else np.min(all_wave)
    wav_max = cubepar['wave_max'] if cubepar['wave_max'] is not None else np.max(all_wave)
    if cubepar['wave_delta'] is not None: dwv = cubepar['wave_delta']

    # Generate a master WCS to register all frames
    coord_min = [ra_min, dec_min, wav_min]
    coord_dlt = [dspat, dspat, dwv]
    masterwcs = generate_masterWCS(coord_min, coord_dlt, name=specname)
    msgs.info(msgs.newline()+"-"*40 +
              msgs.newline() + "Parameters of the WCS:" +
              msgs.newline() + "RA   min, max = {0:f}, {1:f}".format(ra_min, ra_max) +
              msgs.newline() + "DEC  min, max = {0:f}, {1:f}".format(dec_min, dec_max) +
              msgs.newline() + "WAVE min, max = {0:f}, {1:f}".format(wav_min, wav_max) +
              msgs.newline() + "Spaxel size = {0:f}''".format(3600.0*dspat) +
              msgs.newline() + "Wavelength step = {0:f} A".format(dwv) +
              msgs.newline() + "-" * 40)

    # Generate the output binning
    numra = int((ra_max-ra_min) * cosdec / dspat)
    numdec = int((dec_max-dec_min)/dspat)
    numwav = int((wav_max-wav_min)/dwv)
    xbins = np.arange(1+numra)-0.5
    ybins = np.arange(1+numdec)-0.5
    spec_bins = np.arange(1+numwav)-0.5
    bins = (xbins, ybins, spec_bins)

    # Make the cube
    msgs.info("Generating pixel coordinates")
    pix_coord = masterwcs.wcs_world2pix(all_ra, all_dec, all_wave * 1.0E-10, 0)
    hdr = masterwcs.to_header()

    # Generate a datacube using nearest grid point (NGP)
    msgs.info("Generating data cube")
    generate_cube_ngp(outfile, hdr, all_sci, all_ivar, all_wghts, pix_coord, bins, overwrite=overwrite,
                      blaze_wave=blaze_wave, blaze_spec=blaze_spec, fluxcal=fluxcal, specname=specname)
