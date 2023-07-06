"""
Module for GTC OSIRIS specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from astropy import wcs, units
import astropy.io.fits as fits
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation


class GTCOSIRISPlusSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle GTC/OSIRIS specific code
    """
    ndet = 1
    name = 'gtc_osiris_plus'
    telescope = telescopes.GTCTelescopePar()
    camera = 'OSIRIS'
    url = 'http://www.gtc.iac.es/instruments/osiris/'
    header_name = 'OSIRIS'
    supported = True
    comment = 'See :doc:`gtc_osiris`'

    def __init__(self):
        super().__init__()
        self.location = EarthLocation.of_site('lapalma')

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.gtc.iac.es/instruments/osiris/>`__.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')
        gain = 1.90 if hdu is None else self.get_headarr(hdu)[0]['GAIN']
        ronoise = 4.3 if hdu is None else self.get_headarr(hdu)[0]['RDNOISE']

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.125,  # arcsec per pixel
            darkcurr        = 5.0,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = 0,
            numamplifiers   = 1,
            gain            = np.atleast_1d([gain]),
            ronoise         = np.atleast_1d([ronoise]),
            datasec         = np.atleast_1d('[180:4112,50:4096]'),
            oscansec        = np.atleast_1d('[180:4112,8:46]')  # Trim down the oscansec - looks like some bad pixels
            )

        detectors = [detector_dict1]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI,ArI']

        # Set the default exposure time ranges for the frame typing
        par['scienceframe']['exprng'] = [90, None]
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['standardframe']['exprng'] = [None, 180]
        # Multiple arcs with different lamps, so can't median combine nor clip, also need to remove continuum
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Increase the wave tilts order, given the longish slit
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5

        # Only extract one object per standard frame
        par['reduce']['findobj']['maxnumber_std'] = 1

        # Turn off the 2D fit - this seems to be giving bad results for OSIRIS
        par['reduce']['skysub']['no_poly'] = True
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='object')
        self.meta['idname'] = dict(ext=0, card='obsmode')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['detector'] = dict(ext=0, card='detector')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['datasec'] = dict(ext=0, card='DETSIZE')
        self.meta['dichroic'] = dict(ext=0, card='FILTER1')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['slitwid'] = dict(card=None, compound=True)

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)[::-1]
            return binning
        elif meta_key == 'pressure':
            try:
                return headarr[0]['PRESSURE'] * 0.001  # Must be in astropy.units.bar
            except KeyError:
                msgs.warn("Pressure is not in header")
                return 0.0
        elif meta_key == 'temperature':
            try:
                return headarr[0]['TAMBIENT']  # Must be in astropy.units.deg_C
            except KeyError:
                msgs.warn("Temperature is not in header")
                return 0.0
        elif meta_key == 'humidity':
            try:
                return headarr[0]['HUMIDITY']
            except KeyError:
                msgs.warn("Humidity is not in header")
                return 0.0
        elif meta_key == 'obstime':
            return Time(headarr[0]['DATE-END'])
        elif meta_key == 'gain':
            return headarr[0]['GAIN']
        elif meta_key == 'slitwid':
            if self.name == "gtc_maat":
                msgs.warn("HACK FOR MAAT SIMS --- NEED TO GET SLICER SCALE FROM HEADER, IDEALLY")
                return 0.305 / 3600.0
            elif self.name == "gtc_osiris_plus":
                return headarr[0]['SLITW']/3600.0   # Convert slit width from arcseconds to degrees
            else:
                msgs.error("Could not determine slit width from header information")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'decker', 'binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRISM', 'MASKNAME', 'CCDSUM', 'obsmode']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (np.logical_not(np.char.startswith(np.char.lower(fitstbl['target']), 'arclamp'))) & \
                   (np.char.lower(fitstbl['target']) != 'spectralflat') & \
                   (np.char.lower(fitstbl['target']) != 'bias')
        if ftype in ['arc', 'tilt']:
            return good_exp & (np.char.startswith(np.char.lower(fitstbl['target']), 'arclamp'))
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (np.char.lower(fitstbl['target']) == 'spectralflat')
        if ftype == 'bias':
            return good_exp & (np.char.lower(fitstbl['target']) == 'bias')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are considered independent of a configuration.
        Standards are assigned to the correct configuration frame group by
        grism (i.e. ignoring that they are taken with a wider slit).
        See :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'standard': 'dispname', 'bias': None, 'dark': None}

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'idname') == 'OsirisMOS':
            par['reduce']['findobj']['find_trim_edge'] = [1,1]
            par['calibrations']['slitedges']['sync_predict'] = 'pca'
            par['calibrations']['slitedges']['det_buffer'] = 1
        elif self.get_meta_value(scifile, 'idname') == 'OsirisLongSlitSpectroscopy':
            # Do not tweak the slit edges for longslit
            par['calibrations']['flatfield']['tweak_slits'] = False

        # Wavelength calibration and setup-dependent parameters
        if self.get_meta_value(scifile, 'dispname') == 'R300B':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300B.fits'
            par['reduce']['findobj']['find_min_max'] = [750, 2051]
            par['calibrations']['slitedges']['det_min_spec_length'] = 0.25
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.25
            par['calibrations']['slitedges']['smash_range'] = [0.38, 0.62]
            par['calibrations']['flatfield']['slit_illum_finecorr'] = False
            par['reduce']['cube']['wave_min'] = 3600.0
            par['reduce']['cube']['wave_max'] = 7200.0
        elif self.get_meta_value(scifile, 'dispname') == 'R300R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300R.fits'
            par['reduce']['findobj']['find_min_max'] = [750, 2051]
            par['calibrations']['slitedges']['det_min_spec_length'] = 0.25
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.25
            par['calibrations']['slitedges']['smash_range'] = [0.38, 0.62]
            par['calibrations']['flatfield']['slit_illum_finecorr'] = False
            par['reduce']['cube']['wave_min'] = 4800.0
            par['reduce']['cube']['wave_max'] = 10000.0
        elif self.get_meta_value(scifile, 'dispname') == 'R500B':
            par['calibrations']['wavelengths']['lamps'] = ['HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500B.fits'
            par['reduce']['findobj']['find_min_max'] = [500, 2051]
            par['reduce']['cube']['wave_min'] = 3600.0
            par['reduce']['cube']['wave_max'] = 7200.0
        elif self.get_meta_value(scifile, 'dispname') == 'R500R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500R.fits'
            par['reduce']['findobj']['find_min_max'] = [450, 2051]
            par['reduce']['cube']['wave_min'] = 4800.0
            par['reduce']['cube']['wave_max'] = 10000.0
        elif self.get_meta_value(scifile, 'dispname') == 'R1000B':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R1000R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2000B':
            par['calibrations']['wavelengths']['fwhm'] = 15.0
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500U':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500U.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500V':
            par['calibrations']['wavelengths']['lamps'] = ['HgI','NeI','XeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500V.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500R':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500I':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,XeI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500I.fits'
            par['sensfunc']['algorithm'] = 'IR'
            par['sensfunc']['IR']['telgridfile'] = "TelFit_MaunaKea_3100_26100_R20000.fits"
        else:
            msgs.warn('gtc_osiris.py: template arc missing for this grism! Trying holy-grail...')
            par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Return
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels. **This is
                ignored for KCWI.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm; msbias is always set to None.
        bpm_img = super().bpm(filename, det, shape=shape, msbias=None)

        # Extract some header info
        head0 = fits.getheader(filename, ext=0)
        binning = self.get_meta_value([head0], 'binning')

        msgs.warn("Bad pixel mask is not available for det={0:d} binning={1:s}".format(det, binning))
        # Construct a list of the bad columns
        bc = []
        # TODO :: Add BPM

        # Apply these bad columns to the mask
        for bb in range(len(bc)):
            bpm_img[bc[bb][2]:bc[bb][3] + 1, bc[bb][0]:bc[bb][1] + 1] = 1

        return bpm_img


class GTCMAATSpectrograph(GTCOSIRISPlusSpectrograph):
    pypeline = 'IFU'
    name = 'gtc_maat'

    def init_meta(self):
        super().init_meta()
        self.meta['obstime'] = dict(card=None, compound=True, required=False)
        self.meta['pressure'] = dict(card=None, compound=True, required=False)
        self.meta['temperature'] = dict(card=None, compound=True, required=False)
        self.meta['humidity'] = dict(card=None, compound=True, required=False)

    @classmethod
    def default_pypeit_par(cls):
        par = super().default_pypeit_par()

        # LACosmics parameters
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5
        par['scienceframe']['process']['use_illumflat'] = False  # illumflat is applied when building the relative scale image in reduce.py, so should be applied to scienceframe too.
        par['scienceframe']['process']['use_specillum'] = False  # apply relative spectral illumination
        par['scienceframe']['process']['spat_flexure_correct'] = False  # don't correct for spatial flexure - varying spatial illumination profile could throw this correction off. Also, there's no way to do astrometric correction if we can't correct for spatial flexure of the contbars frames
        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_darkimage'] = False
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        # Don't do 1D extraction for 3D data - it's meaningless because the DAR correction must be performed on the 3D data.
        par['reduce']['extraction']['skip_extraction'] = True  # Because extraction occurs before the DAR correction, don't extract

        # Decrease the wave tilts order, given the shorter slits of the IFU
        par['calibrations']['tilts']['spat_order'] = 1
        par['calibrations']['tilts']['spec_order'] = 1

        # Make sure that this is reduced as a slit (as opposed to fiber) spectrograph
        par['reduce']['cube']['slit_spec'] = True
        par['reduce']['cube']['combine'] = False  # Make separate spec3d files from the input spec2d files

        # Sky subtraction parameters
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = False
        par['reduce']['findobj']['skip_skysub'] = True
        par['reduce']['findobj']['skip_final_global'] = True

        # Don't correct flexure by default, but you should use slitcen,
        # because this is a slit-based IFU where no objects are extracted.
        par['flexure']['spec_method'] = 'skip'
        par['flexure']['spec_maxshift'] = 3  # Just in case someone switches on spectral flexure, this needs to be minimal

        # Flux calibration parameters
        par['sensfunc']['UVIS']['extinct_correct'] = False  # This must be False - the extinction correction is performed when making the datacube

        return par

    def get_wcs(self, hdr, slits, platescale, wave0, dwv, spatial_scale=None):
        """
        Construct/Read a World-Coordinate System for a frame.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame. The information in this
                header will be extracted and returned as a WCS.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces.
            platescale (:obj:`float`):
                The platescale of an unbinned pixel in arcsec/pixel (e.g.
                detector.platescale).
            wave0 (:obj:`float`):
                The wavelength zeropoint.
            dwv (:obj:`float`):
                Change in wavelength per spectral pixel.

        Returns:
            `astropy.wcs.wcs.WCS`_: The world-coordinate system.
        """
        msgs.info("Calculating the WCS")
        # Get the x and y binning factors, and the typical slit length
        binspec, binspat = parse.parse_binning(self.get_meta_value([hdr], 'binning'))

        # Get the pixel and slice scales
        pxscl = platescale * binspat / 3600.0  # Need to convert arcsec to degrees
        slscl = self.get_meta_value([hdr], 'slitwid')
        if spatial_scale is not None:
            if pxscl > spatial_scale / 3600.0:
                msgs.warn("Spatial scale requested ({0:f}'') is less than the pixel scale ({1:f}'')".format(spatial_scale, pxscl*3600.0))
            # Update the pixel scale
            pxscl = spatial_scale / 3600.0  # 3600 is to convert arcsec to degrees

        # Get the typical slit length (this changes by ~0.3% over all slits, so a constant is fine for now)
        slitlength = int(np.round(np.median(slits.get_slitlengths(initial=True, median=True))))

        # Get RA/DEC
        raval = self.get_meta_value([hdr], 'ra')
        decval = self.get_meta_value([hdr], 'dec')

        # Create a coordinate
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))

        # Get rotator position
        msgs.warn("HACK FOR MAAT SIMS --- NEED TO FIGURE OUT RPOS and RREF FOR MAAT FROM HEADER INFO")
        if 'ROTPOSN' in hdr:
            rpos = hdr['ROTPOSN']
        else:
            rpos = 0.
        if 'ROTREFAN' in hdr:
            rref = hdr['ROTREFAN']
        else:
            rref = 0.
        # Get the offset and PA
        rotoff = 0.0  # IFU-SKYPA offset (degrees)
        skypa = rpos + rref  # IFU position angle (degrees)
        crota = np.radians(-(skypa + rotoff))

        # Calculate the fits coordinates
        cdelt1 = -slscl
        cdelt2 = pxscl
        if coord is None:
            ra = 0.
            dec = 0.
            crota = 1
        else:
            ra = coord.ra.degree
            dec = coord.dec.degree
        # Calculate the CD Matrix
        cd11 = cdelt1 * np.cos(crota)                          # RA degrees per column
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)   # RA degrees per row
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)  # DEC degress per column
        cd22 = cdelt2 * np.cos(crota)                          # DEC degrees per row
        # Get reference pixels (set these to the middle of the FOV)
        crpix1 = 11   # i.e. see get_datacube_bins (11 is used as the reference point - somewhere in the middle of the FOV)
        crpix2 = slitlength / 2.
        crpix3 = 1.
        # Get the offset
        msgs.warn("HACK FOR MAAT SIMS --- Need to obtain offset from header?")
        off1 = 0.
        off2 = 0.
        off1 /= binspec
        off2 /= binspat
        crpix1 += off1
        crpix2 += off2

        # Create a new WCS object.
        msgs.info("Generating MAAT WCS")
        w = wcs.WCS(naxis=3)
        w.wcs.equinox = hdr['EQUINOX']
        w.wcs.name = 'MAAT'
        w.wcs.radesys = 'FK5'
        # Insert the coordinate frame
        w.wcs.cname = ['MAAT RA', 'MAAT DEC', 'MAAT Wavelength']
        w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
        w.wcs.crval = [ra, dec, wave0]  # RA, DEC, and wavelength zeropoints
        w.wcs.crpix = [crpix1, crpix2, crpix3]  # RA, DEC, and wavelength reference pixels
        w.wcs.cd = np.array([[cd11, cd12, 0.0], [cd21, cd22, 0.0], [0.0, 0.0, dwv]])
        w.wcs.lonpole = 180.0  # Native longitude of the Celestial pole
        w.wcs.latpole = 0.0  # Native latitude of the Celestial pole

        return w

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        r"""
        Calculate the bin edges to be used when making a datacube.

        Args:
            slitlength (:obj:`int`):
                Length of the slit in pixels
            minmax (`numpy.ndarray`_):
                An array with the minimum and maximum pixel locations on each
                slit relative to the reference location (usually the centre
                of the slit). Shape must be :math:`(N_{\rm slits},2)`, and is
                typically the array returned by
                :func:`~pypeit.slittrace.SlitTraceSet.get_radec_image`.
            num_wave (:obj:`int`):
                Number of wavelength steps.  Given by::
                    int(round((wavemax-wavemin)/delta_wave))

        Args:
            :obj:`tuple`: Three 1D `numpy.ndarray`_ providing the bins to use
            when constructing a histogram of the spec2d files. The elements
            are :math:`(x,y,\lambda)`.
        """
        xbins = np.arange(1 + 23) - 11.0 - 0.5  # 23 is for 23 slices, and 11 is the reference slit
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]), 1+slitlength) - 0.5
        spec_bins = np.arange(1+num_wave) - 0.5
        return xbins, ybins, spec_bins


class GTCOSIRISSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle GTC/OSIRIS specific code (old detector: MAT-44-82)
    """
    ndet = 2
    name = 'gtc_osiris'
    telescope = telescopes.GTCTelescopePar()
    camera = 'OSIRIS'
    url = 'http://www.gtc.iac.es/instruments/osiris/'
    header_name = 'OSIRIS'
    supported = True
    comment = 'See :doc:`gtc_osiris`'

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Detector data from `here
        <http://www.gtc.iac.es/instruments/osiris/>`__.

        Args:
            det (:obj:`int`):
                1-indexed detector number.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')
        #detwin2 = '[1:4102,300:600]' if self.par['rdx']['quicklook'] else '[1:4102,52:1920]'

        # Detector 1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 1,  # Not sure this is used
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.127,  # arcsec per pixel
            darkcurr        = 0.0,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = 0,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,280:2048]'),
            oscansec        = np.atleast_1d('[1:4102,6:44]')
            )
        # Detector 2
        detector_dict2 = dict(
            binning         = binning,
            det             = 2,
            dataext         = 2,  # Not sure this is used
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.127,
            darkcurr        = 0.0,
            saturation      = 65535., # ADU
            nonlinear       = 0.95,
            mincounts       = 0,
            numamplifiers   = 1,
            gain            = np.atleast_1d([0.95]),
            ronoise         = np.atleast_1d([4.5]),
            datasec         = np.atleast_1d('[1:4102,52:1920]'),
            oscansec        = np.atleast_1d('[1:4102,6:40]')
            )

        detectors = [detector_dict1, detector_dict2]
        # Return
        return detector_container.DetectorContainer(**detectors[det-1])


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI,ArI']

        # Set the default exposure time ranges for the frame typing
        par['scienceframe']['exprng'] = [90, None]
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['standardframe']['exprng'] = [None, 180]
        # Multiple arcs with different lamps, so can't median combine nor clip, also need to remove continuum
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True

        # Increase the wave tilts order, given the longish slit
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5

        #Only extract one object per standard frame
        par['reduce']['findobj']['maxnumber_std']=1

        # Turn off the 2D fit - this seems to be giving bad results for OSIRIS
        par['reduce']['skysub']['no_poly'] = True
        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='object')
        self.meta['idname'] = dict(ext=0, card='obsmode')
        self.meta['decker'] = dict(ext=0, card='MASKNAME')
        self.meta['binning'] = dict(card=None, compound=True)  # Uses CCDSUM
        self.meta['detector']=dict(ext=0,card='detector')
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        self.meta['dispname'] = dict(ext=0, card='GRISM')
        self.meta['datasec'] = dict(ext=1, card='DATASEC')
        self.meta['dichroic'] = dict(ext=0, card='FILTER1')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')

    def compound_meta(self, headarr, meta_key):
        """
        Methods to generate metadata requiring interpretation of the header
        data, instead of simply reading the value of a header card.

        Args:
            headarr (:obj:`list`):
                List of `astropy.io.fits.Header`_ objects.
            meta_key (:obj:`str`):
                Metadata keyword to construct.

        Returns:
            object: Metadata value read from the header(s).
        """
        if meta_key == 'binning':
            binspatial, binspec = parse.parse_binning(headarr[0]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'obstime':
            return Time(headarr[0]['DATE-END'])

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to constuct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['dispname', 'decker', 'binning']

    def raw_header_cards(self):
        """
        Return additional raw header cards to be propagated in
        downstream output files for configuration identification.

        The list of raw data FITS keywords should be those used to populate
        the :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.configuration_keys`
        or are used in :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.config_specific_par`
        for a particular spectrograph, if different from the name of the
        PypeIt metadata keyword.

        This list is used by :meth:`~pypeit.spectrographs.spectrograph.Spectrograph.subheader_for_spec`
        to include additional FITS keywords in downstream output files.

        Returns:
            :obj:`list`: List of keywords from the raw data files that should
            be propagated in output files.
        """
        return ['GRISM', 'MASKNAME', 'CCDSUM', 'obsmode']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Args:
            ftype (:obj:`str`):
                Type of frame to check. Must be a valid frame type; see
                frame-type :ref:`frame_type_defs`.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science','standard']:
            return good_exp & (fitstbl['target'] != 'ArcLamp_Xe') \
            & (fitstbl['target'] != 'ArcLamp_HgAr') \
            & (fitstbl['target'] != 'ArcLamp_Ne') \
            & (fitstbl['target'] != 'SpectralFlat') \
                    & (fitstbl['target'] != 'BIAS')
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'ArcLamp_Xe') \
            | (fitstbl['target'] == 'ArcLamp_HgAr') \
            | (fitstbl['target'] == 'ArcLamp_Ne'))
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['target'] == 'SpectralFlat')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are considered independent of a configuration.
        Standards are assigned to the correct configuration frame group by
        grism (i.e. ignoring that they are taken with a wider slit).
        See :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'standard': 'dispname','bias': None, 'dark': None}

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            scifile (:obj:`str`):
                File to use when determining the configuration and how
                to adjust the input parameters.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument wide
        par = super().config_specific_par(scifile, inp_par=inp_par)

        if self.get_meta_value(scifile, 'idname') == 'OsirisMOS':
            par['reduce']['findobj']['find_trim_edge'] = [1,1]
            par['calibrations']['slitedges']['sync_predict'] = 'pca'
            par['calibrations']['slitedges']['det_buffer'] = 1
        elif self.get_meta_value(scifile, 'idname') == 'OsirisLongSlitSpectroscopy':
            # Do not tweak the slit edges for longslit
            par['calibrations']['flatfield']['tweak_slits'] = False

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == 'R300B':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300B.fits'
            par['reduce']['findobj']['find_min_max']=[750,2051]
        elif self.get_meta_value(scifile, 'dispname') == 'R300R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R300R.fits'
            par['reduce']['findobj']['find_min_max']=[750,2051]
        elif self.get_meta_value(scifile, 'dispname') == 'R500B':
            par['calibrations']['wavelengths']['lamps'] = ['HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500B.fits'
            par['reduce']['findobj']['find_min_max']=[500,2051]
        elif self.get_meta_value(scifile, 'dispname') == 'R500R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R500R.fits'
            par['reduce']['findobj']['find_min_max']=[450,2051]
        elif self.get_meta_value(scifile, 'dispname') == 'R1000B':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R1000R':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R1000R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2000B':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2000B.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500U':
            par['calibrations']['wavelengths']['lamps'] = ['XeI,HgI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500U.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500V':
            par['calibrations']['wavelengths']['lamps'] = ['HgI','NeI','XeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500V.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500R':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,HgI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500R.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'R2500I':
            par['calibrations']['wavelengths']['lamps'] = ['ArI,XeI,NeI']
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gtc_osiris_R2500I.fits'
            par['sensfunc']['algorithm'] = 'IR'
            par['sensfunc']['IR']['telgridfile'] = "TelFit_MaunaKea_3100_26100_R20000.fits"
        else:
            msgs.warn('gtc_osiris.py: template arc missing for this grism! Trying holy-grail...')
            par['calibrations']['wavelengths']['method'] = 'holy-grail'

        # Return
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Even though they are both optional, either the precise shape for
        the image (``shape``) or an example file that can be read to get
        the shape (``filename`` using :func:`get_image_shape`) *must* be
        provided.

        Args:
            filename (:obj:`str` or None):
                An example file to use to get the image shape.
            det (:obj:`int`):
                1-indexed detector number to use when getting the image
                shape from the example file.
            shape (tuple, optional):
                Processed image shape
                Required if filename is None
                Ignored if filename is not None
            msbias (`numpy.ndarray`_, optional):
                Processed bias frame used to identify bad pixels. **This is
                ignored for KCWI.**

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm; msbias is always set to None.
        bpm_img = super().bpm(filename, det, shape=shape, msbias=None)

        # Extract some header info
        head0 = fits.getheader(filename, ext=0)
        binning = head0['CCDSUM']

        # Construct a list of the bad columns
        bc = []
        if det == 1:
            # No bad pixel columns on detector 1
            pass
        elif det == 2:
            if binning == '1 1':
                # The BPM is based on 2x2 binning data, so the 2x2 numbers are just multiplied by two
                msgs.warn("BPM is likely over-estimated for 1x1 binning")
                bc = [[220, 222, 3892, 4100],
                      [952, 954, 2304, 4100]]
            elif binning == '2 2':
                bc = [[110, 111, 1946, 2050],
                      [476, 477, 1154, 2050]]
        else:
            msgs.warn("Bad pixel mask is not available for det={0:d} binning={1:s}".format(det, binning))
            bc = []

        # Apply these bad columns to the mask
        for bb in range(len(bc)):
            bpm_img[bc[bb][2]:bc[bb][3] + 1, bc[bb][0]:bc[bb][1] + 1] = 1

        return bpm_img
