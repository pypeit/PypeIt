"""
Module containing routines used by 3D datacubes.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import inspect

from astropy import wcs, units
from astropy.coordinates import AltAz, SkyCoord
from astropy.io import fits
import scipy.optimize as opt
from scipy.interpolate import interp1d
import numpy as np

from pypeit import msgs
from pypeit.core.procimg import grow_masked
from pypeit.core import coadd
from pypeit.spectrographs.util import load_spectrograph
from pypeit import datamodel
from pypeit import io


class DataCube(datamodel.DataContainer):
    """
    DataContainer to hold the products of a datacube

    See the datamodel for argument descriptions

    Args:
        flux (`numpy.ndarray`_):
            The science datacube (nwave, nspaxel_y, nspaxel_x)
        variance (`numpy.ndarray`_):
            The variance datacube (nwave, nspaxel_y, nspaxel_x)
        refscale (`numpy.ndarray`_):
            An image containing the relative scale of each pixel on the detector (nspec, nspat)
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
    version = '1.0.0'

    datamodel = {'flux': dict(otype=np.ndarray, atype=np.floating, descr='Flux array in units of counts/s or 10^-17 erg/s/cm^2/Ang'),
                 'variance': dict(otype=np.ndarray, atype=np.floating, descr='Variance array (matches units of flux)'),
                 'refscale': dict(otype=np.ndarray, atype=np.floating, descr='Reference scaling used for each slit'),
                 'PYP_SPEC': dict(otype=str, descr='PypeIt: Spectrograph name'),
                 'fluxed': dict(otype=bool, descr='Boolean indicating if the datacube is fluxed.')}

    @classmethod
    def from_file(cls, ifile):
        """
        Over-load :func:`pypeit.datamodel.DataContainer.from_file`
        to deal with the header

        Args:
            ifile (str):  Filename holding the object

        Returns:
            :class:`OneSpec`:

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

    def __init__(self, flux, variance, PYP_SPEC, refscale=None, fluxed=None):

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
        res_lsq = opt.least_squares(dar_fitfunc, [0.0, 0.0], args=args, xtol=1.0e-6, ftol=None, gtol=None)
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


def make_whitelight_fromref(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat, ref_filename):
    """ Generate a whitelight image of every input frame,
    based on a reference image. Note the, the reference
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
    whitelight_imgs, _, _ = make_whitelight(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                                            whitelightWCS=wlwcs, numra=numra, numdec=numdec)
    # Return required info
    return reference_image, whitelight_imgs, wlwcs


def make_whitelight(all_ra, all_dec, all_wave, all_sci, all_wghts, all_idx, dspat,
                    all_ivar=None, whitelightWCS=None, numra=None, numdec=None):
    """ Generate a whitelight image of every input frame

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
    #from IPython import embed
    #embed()
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
        gpm = grow_masked(whtlght == 0, trim, 1) == 0  # A good pixel = 1
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

