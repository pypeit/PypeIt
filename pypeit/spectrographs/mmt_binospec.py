"""
Module for MMT/BINOSPEC specific methods.

.. include:: ../include/links.rst
"""
from itertools import chain
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.time import Time
from astropy import wcs
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import uniform_filter1d
from scipy.signal import correlate
from scipy.special import erf

from pypeit import dataPaths
from pypeit.pkg.pypeitdata import PypeItDataPath
from pypeit import io
from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit.core import framematch
from pypeit.core import parse
from pypeit.core import procimg
from pypeit.core.moment import moment1d
from pypeit.images import detector_container
from pypeit.par import parset
from pypeit.spectrographs import spectrograph
from pypeit.spectrographs.slitmask import SlitMask


class MMTBINOSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MMT/BINOSPEC specific code
    """
    ndet = 2
    name = 'mmt_binospec'

    # Per-amplifier nonlinearity correction coefficients from IDL calibration
    # file (scicam_bino_sep2017.fits, measured 2017-09). Polynomial form:
    # C_corr = c[0] + c[1]*C + c[2]*C^2 + c[3]*C^3 + c[4]*C^4
    # Compatible with np.polynomial.polynomial.polyval.
    # Shape: (8 amplifiers, 5 coefficients for degree-4 polynomial)
    nonlinearity_coeffs = np.array([
        [0.00000000e+00, 1.00400089e+00, -1.39235362e-06, 8.31711824e-12, -1.20653479e-17],
        [0.00000000e+00, 1.00361458e+00, -1.29223833e-06, 6.93723177e-12, -9.67406255e-18],
        [0.00000000e+00, 1.00269542e+00, -9.29361806e-07, 5.97902827e-12, -2.30257302e-17],
        [0.00000000e+00, 1.00339616e+00, -8.47134521e-07, 7.92441693e-12, -4.46542834e-17],
        [0.00000000e+00, 1.00727205e+00, -1.69093388e-06, 2.07225055e-11, -1.62655178e-16],
        [0.00000000e+00, 1.00858745e+00, -2.35668901e-06, 2.40641019e-11, -1.50286358e-16],
        [0.00000000e+00, 1.00728526e+00, -1.80779473e-06, 1.73427719e-11, -1.01685780e-16],
        [0.00000000e+00, 1.00845168e+00, -2.02050567e-06, 2.97587091e-11, -2.65508521e-16],
    ])
    telescope = telescopes.MMTTelescopePar()
    camera = 'BINOSPEC'
    url = 'https://lweb.cfa.harvard.edu/mmti/binospec.html'
    header_name = 'Binospec'
    supported = True

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

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
        # Binning
        binning = '1,1' if hdu is None else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict1 = dict(
                            binning         = binning,
                            det             = 1,
                            dataext         = 1,
                            specaxis        = 0,
                            specflip        = False,
                            spatflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.24,
                            darkcurr        = 3.6,  #e-/pixel/hour  (=0.001 e-/pixel/s)  --  pulled from the ETC
                            saturation      = 65535.,
                            nonlinear       = 0.95,  #ToDO: To Be update
                            mincounts       = -1e10,
                            numamplifiers   = 4,
                            gain            = np.atleast_1d([1.085,1.046,1.042,0.975]),
                            ronoise         = np.atleast_1d([3.2,3.2,3.2,3.2]),
                            )
        # Detector 2
        detector_dict2 = detector_dict1.copy()
        detector_dict2.update(dict(
            det=2,
            dataext=2,
            gain=np.atleast_1d([1.028,1.115,1.047,1.045]), #ToDo: FW measures 1.115 for amp2 but 1.163 in IDL pipeline
            ronoise=np.atleast_1d([3.6,3.6,3.6,3.6])
        ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer(**detector_dicts[det-1])

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=1, card='RA')
        self.meta['dec'] = dict(ext=1, card='DEC')
        self.meta['target'] = dict(ext=1, card='OBJECT')
        self.meta['decker'] = dict(ext=1, card='MASK')

        self.meta['dichroic'] = dict(ext=1, card=None, default='default')
        self.meta['binning'] = dict(ext=1, card='CCDSUM', compound=True)

        self.meta['mjd'] = dict(ext=1, card='MJD')
        self.meta['exptime'] = dict(ext=1, card='EXPTIME')
        self.meta['airmass'] = dict(ext=1, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=1, card='DISPERS1')
        self.meta['idname'] = dict(ext=1, card='IMAGETYP')

        # used for arclamp
        self.meta['lampstat01'] = dict(ext=1, card='HENEAR')
        # used for flatlamp, SCRN is actually telescope status
        self.meta['lampstat02'] = dict(ext=1, card='SCRN')
        # incandescent (flat-field) lamp; required to be 'on' for flat/trace frames
        self.meta['lampstat03'] = dict(ext=1, card='INCAN')
        self.meta['instrument'] = dict(ext=1, card='INSTRUME')

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
            binspatial, binspec = parse.parse_binning(headarr[1]['CCDSUM'])
            binning = parse.binning2string(binspec, binspatial)
            return binning

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
        return ['dispname', 'decker']

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
        return ['DISPERS1', 'MASK']

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.125
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm']= 4.0
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'NeI', 'ArI', 'ArII']

        # Tilt and slit parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        par['calibrations']['tilts']['spat_order'] = 6
        par['calibrations']['tilts']['spec_order'] = 6
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Processing steps
        turn_off = dict(use_biasimage=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0
        ## Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['global_sky_std']  = False

        par['flexure']['spec_method'] = 'boxcar'

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Cosmic-ray cleanup parameters for arc and tilt frames.  Setting
        # mask_cr=True together with cr_median_width > 1 triggers the
        # residual-CR path inside PypeItImage.build_crmask
        # (pypeit.core.procimg.lacosmic_spatial_median_residual): a
        # row-local median is subtracted and LA Cosmic is run on the
        # residual, catching extended trails that the standard 3x3
        # Laplacian misses.  CR-flagged pixels are masked (CR bit set in
        # fullmask); image values are unmodified, and downstream tilt and
        # arc fitting honor the CR flag.
        for _frame in ('arcframe', 'tiltframe'):
            par['calibrations'][_frame]['process']['mask_cr'] = True
            par['calibrations'][_frame]['process']['cr_median_width'] = 31
            par['calibrations'][_frame]['process']['sigclip'] = 10.0
            par['calibrations'][_frame]['process']['sigfrac'] = 0.3
            par['calibrations'][_frame]['process']['objlim'] = 0.0
            par['calibrations'][_frame]['process']['lamaxiter'] = 2
            par['calibrations'][_frame]['process']['grow'] = 2.0
            par['calibrations'][_frame]['process']['rmcompact'] = False

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 100]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp (:obj:`str`, :obj:`list`, `Path`_, `astropy.io.fits.Header`_, `astropy.table.Table`_):
                Input filename, an `astropy.io.fits.Header`_ object, or a list
                of `astropy.io.fits.Header`_ objects.  Or a row from the
                metadata table.
            inp_par (:class:`~pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: The PypeIt parameter set
            adjusted for configuration specific parameter values.
        """
        # Start with instrument-wide parameters
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Adjust parameters based on instrument configuration
        grating = self.get_meta_value(inp, 'dispname')
        decker = self.get_meta_value(inp, 'decker')

        # wavelengths
        match grating:
            case 'x270':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_binospec_270.fits'
            case 'x600':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_binospec_600.fits'
            case 'x1000':
                par['calibrations']['wavelengths']['reid_arxiv'] = 'mmt_binospec_1000.fits'

        if 'Longslit' in decker:
            # Observations use a longslit so we skip the parameters primarily
            # used for multislit data
            return par

        # Turn on the use of mask design
        par['calibrations']['slitedges']['use_maskdesign'] = True
        # Since we use the slitmask info to find the alignment boxes, I don't need `minimum_slit_length_sci`
        par['calibrations']['slitedges']['minimum_slit_length_sci'] = None
        # Sometime the added missing slits at the edge of the detector are to small to be useful.
        par['calibrations']['slitedges']['minimum_slit_length'] = 3.
        # Since we use the slitmask info to add and remove traces, 'minimum_slit_gap' may undo the matching effort.
        par['calibrations']['slitedges']['minimum_slit_gap'] = 0.
        # Lower edge_thresh works better
        par['calibrations']['slitedges']['edge_thresh'] = 10.
        # Assign RA, DEC, OBJNAME to detected objects
        par['reduce']['slitmask']['assign_obj'] = True
        # force extraction of undetected objects
        par['reduce']['slitmask']['extract_missing_objs'] = True
        # Adjust sky subtraction parameters

        # lower tilts spat_order and higher spec_order for multislits (i.e., generally not very long slits)
        par['calibrations']['tilts']['spat_order'] = 2  # Default: 3
        par['calibrations']['tilts']['spec_order'] = 5  # Default: 4
        # pca
        par['calibrations']['slitedges']['sync_predict'] = 'auto'

        par['coadd2d']['offsets'] = 'maskdef_offsets'

        return par

    def update_edgetracepar(self, par):
        """
        This method is used in :func:`pypeit.edgetrace.EdgeTraceSet.maskdesign_matching`
        to update EdgeTraceSet parameters when the slitmask design matching is not feasible
        because too few slits are present in the detector.

        Args:
            par (:class:`pypeit.par.pypeitpar.EdgeTracePar`):
                The parameters used to guide slit tracing.

        Returns:
            :class:`pypeit.par.pypeitpar.EdgeTracePar`
            The modified parameters used to guide slit tracing.
        """

        par['minimum_slit_gap'] = 0.25
        par['minimum_slit_length_sci'] = 4.5
        return par

    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Generate a default bad-pixel mask.

        Loads a pre-built static BPM from the IDL pipeline calibration data
        (``badpix_binospec.fits`` + hard-coded bad columns and detector trap
        regions from ``bino_mosaic.pro``).

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
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        # Load and apply the static BPM from IDL pipeline calibration
        bpm_file = dataPaths.static_calibs.get_file_path(
            f'mmt_binospec/bpm_binospec_det{det}.fits.gz')
        static_bpm = fits.getdata(bpm_file)
        bpm_img |= static_bpm

        return bpm_img

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
        if ftype == 'science':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'stowed') & (fitstbl['exptime'] > 100.0)
        if ftype == 'standard':
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'stowed') & (fitstbl['exptime'] <= 100.0)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] == 'on')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Require the incandescent (flat) lamp to be on; INCAN=off frames
            # are dark/noisy exposures and must not be combined into flats.
            return good_exp & (fitstbl['lampstat01'] == 'off') \
                & (fitstbl['lampstat02'] == 'deployed') & (fitstbl['lampstat03'] == 'on')

        log.debug('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_slitmask(self, filename:str, det:int=1):
        """
        Parse the slitmask data from a raw file into :attr:`slitmask`, a
        :class:`~pypeit.spectrographs.slitmask.SlitMask` object.

        Parameters
        ----------
        filename : :obj:`str`
            Name of the file to read.
        det : :obj:`int`, optional
            1-indexed detector number to read the slitmask for.  Must be either
            1 or 2 for MMT/Binospec.

        Returns
        -------
        :class:`~pypeit.spectrographs.slitmask.SlitMask`
            The slitmask data read from the file. The returned object is the
            same as :attr:`slitmask`.

        Notes
        -----
        - Target-slit alignment is characterized via distances from slit edges.
        - Slit corners and on-sky positions are stored for each target.
        """
        slit_id, slit_width, slit_x, slit_y, poly_x, poly_y, \
            obj_id, obj_name, obj_ra, obj_dec, obj_mag, \
            mm_arcsec, rac, decc, posx_pa \
            = self._parse_slitmask_data(filename, det)
        
        # Number of slits
        numslits = slit_id.size

        # The polygon coordinates have a shape that is (Nslits,4), their order
        # is: [0]: bottom left, [1]: top left, [2]: top right, [3]: bottom right
        # TODO: I assume the focal plane is flipped wrt the coordinates
        # provided, which is why the code is as given below.

        # Compute projected distances from target to slit edges in arcseconds
        # left
        topdist = (slit_y - poly_y[:,0]) / mm_arcsec
        # right
        botdist = (poly_y[:,1] - slit_y) / mm_arcsec
        if det == 2:
            # flip for detector 2
            topdist, botdist = botdist, topdist
        slit_length_arcsec = topdist + botdist

        # Assemble object array: [slit_id, id, ra, dec, name, mag, mag_band, top, bot]
        # TODO: I don't know why we need to use round
        objects = np.array([
            slit_id,
            obj_id,
            obj_ra,
            obj_dec,
            obj_name,
            obj_mag,
            ['None'] * numslits,
            np.round(topdist,2),
            np.round(botdist,2)
        ], dtype=object).T

        # Compute slit centers offsets from object positions
        xcen_slit = (poly_x[:,0] + poly_x[:,3]) / 2.
        ycen_slit = (poly_y[:,0] + poly_y[:,1]) / 2.
        slit_xoff = (xcen_slit - slit_x) / mm_arcsec  # in arcseconds
        slit_yoff = (ycen_slit - slit_y) / mm_arcsec  # in arcseconds
        slit_offset = np.sqrt(slit_xoff ** 2 + slit_yoff ** 2)

        # Compute slit center RA/Dec via spherical offset from target position
        obj_coord = SkyCoord(ra=obj_ra, dec=obj_dec, unit='deg')

        # Slit position angles and sign of offset (accounting for up/down location)
        slit_pas = np.full(numslits, posx_pa, dtype=float)
        off_signs = np.ones_like(slit_pas)
        negy = slit_yoff < 0.
        off_signs[negy] = -1.

        # Compute slit center RA/Dec
        slit_ra = np.empty(numslits, dtype=float)
        slit_dec = np.empty(numslits, dtype=float)
        for i, (slit_off, obj_coo, slit_pa, off_sign) in enumerate(zip(
            slit_offset, obj_coord, slit_pas, off_signs
        )):
            slit_coord = obj_coo.directional_offset_by(
                slit_pa * units.deg, off_sign * slit_off * units.arcsec
            )
            slit_ra[i] = slit_coord.ra.deg
            slit_dec[i] = slit_coord.dec.deg

        # TODO: This ordering doesn't seem to match the documentation for
        # SlitMask, but there may be a reflection involved.
        corners = np.stack((poly_x, poly_y), axis=-1)
        corners = corners[:, [3, 0, 1, 2], :]

        # Slitmask pointing coordinates (mask center RA/Dec)
        mask_coord = SkyCoord(rac, decc, unit=('hourangle', 'deg'))

        # Construct and return the slitmask object
        self.slitmask = SlitMask(
            corners,
            slitid=slit_id,
            onsky=np.asarray([
                slit_ra, slit_dec, np.round(slit_length_arcsec, 2),
                slit_width, slit_pas
            ]).T,
            objects=objects,
            mask_radec=(mask_coord.ra.deg, mask_coord.dec.deg),
            posx_pa=posx_pa
        )

        return self.slitmask

    @staticmethod
    def _parse_slitmask_data(filename, det):

        # Open the FITS file
        hdu = io.fits_open(filename)

        # Position angle corresponding to detector +x axis (spatial direction)
        posx_pa = float(hdu[1].header['POSANG']) - 180
        if posx_pa < 0:
            posx_pa += 360.

        # Select appropriate extension for detector 1 or 2
        match det:
            case 1:
                mask_hdu = hdu[9].data[0]
            case 2:
                mask_hdu = hdu[10].data[0]
            case _:
                raise PypeItError(f'Detector number must be 1 or 2 for MMT/Binospec, not {det}.')

        targ = mask_hdu['TARGET_TYPE'] == 'TARGET'
        numslits = mask_hdu['NTARGETS']

        if np.sum(targ) != numslits:
            raise PypeItError(
                f'Expected {numslits} TARGET slits but found {np.sum(targ)} in mask design file.'
            )
        
        # NOTE: The use of np.atleast_* here is to handle the case when there is
        # only one target.

        # Slit properties
        slit_id = np.atleast_1d(mask_hdu['SLIT_ID'])[targ]          # ID number
        slit_width = np.atleast_1d(mask_hdu['SLIT_WIDTH'])[targ]    # in arcsec
        # Target positions in mm
        slit_x = np.atleast_1d(mask_hdu['SLITX'])[targ]
        slit_y = np.atleast_1d(mask_hdu['SLITY'])[targ]
        # Slit polygon coordinates in mm
        poly_x = np.atleast_2d(mask_hdu['POLY_X']).T[targ]
        poly_y = np.atleast_2d(mask_hdu['POLY_Y']).T[targ]

        # Target properties
        obj_id = np.atleast_1d(mask_hdu['TARGET_ID'])[targ]
        obj_name = np.atleast_1d(mask_hdu['TARGET_NAME'])[targ]
        obj_ra = np.atleast_1d(mask_hdu['RA'])[targ]
        obj_dec = np.atleast_1d(mask_hdu['DEC'])[targ]
        obj_mag = np.atleast_1d(mask_hdu['MAG'])[targ]

        # Scalars with the platescale and center coordinates of the slit mask
        mm_arcsec = mask_hdu['MM_PER_ARCSEC']
        rac = mask_hdu['CENTERRA']      # in hours
        decc = mask_hdu['CENTERDEC']    # in degrees

        hdu.close()

        return (
            slit_id, slit_width, slit_x, slit_y, poly_x, poly_y,
            obj_id, obj_name, obj_ra, obj_dec, obj_mag,
            mm_arcsec, rac, decc, posx_pa
        )

    def get_maskdef_slitedges(self, filename:str=None, det:int=1, debug:bool=None, 
                              binning:str=None, trc_path:str=None):
        """
        Provides the slit edges positions predicted by the slitmask design.

        This method is not defined for all spectrographs. This base-class
        method raises an exception. This may be because ``use_maskdesign``
        has been set to True for a spectrograph that does not support it.

        Parameters
        ---------- 
        filename : :obj:`str`, :obj:`list`, optional:
            Name of the file holding the mask design info or the maskfile and
            wcs_file in that order
        det : :obj:`int`, optional
            Detector number
        debug : :obj:`bool`, optional
            Flag to run in debugging mode
        trc_path : str, optional
            Path to the first trace file used to generate the trace flat
        binning : str, optional
            String with the comma-separated number of pixels binned in each
            dimension of the flat-field image.  Order must be spectral then
            spatial.

        Returns
        -------
        top_edges : :class:`numpy.ndarray`
            Predicted locations of the top edges of the slits in spatial pixel
            coordinates.
        bot_edges : :class:`numpy.ndarray`
            Predicted locations of the bottom edges of the slits in spatial pixel
            coordinates.
        sortindx : :class:`numpy.ndarray`
            Indices of the slits in the provided ``slitmask`` object that orders
            the slits from left to right, in the PypeIt orientation.
        slitmask : :class:`~pypeit.spectrographs.slitmask.SlitMask`
            Slit mask metadata read from the provided input file(s).

        Notes
        -----
        - Edges are sorted by bottom edge y-coordinate to order slits spatially.
        """
        if det is None:
            raise ValueError("A valid detector number must be provided.")

        if filename is None:
            raise ValueError("A valid slitmask filename must be provided.")

        # get the full path to the mask design file
        _maskfile = str(Path(trc_path) / filename) if not Path(filename).exists() else filename

        # check if the mask design file exists
        if not Path(_maskfile).exists():
            raise PypeItError(f'The mask design file {_maskfile} does not exist.')

        # Load slitmask information if a file is provided
        self.get_slitmask(_maskfile, det=det)

        if self.slitmask is None:
            raise ValueError("Unable to read slitmask design info. Provide a file.")

        # Open FITS file and read mask data for the correct detector
        hdu = io.fits_open(filename)
        mask_fits = hdu[9].data[0] if det == 1 else hdu[10].data[0]
        # keep only the TARGET slits
        targ = mask_fits['TARGET_TYPE'] == 'TARGET'

        # Define det buffer and mm/pixel scale factor
        # NOTE: these are hard-coded and not sure if there is a more robust way to determine them
        # slitmask offset from the detector edge in pixels
        mask_edge_off = 200
        # scale factor to convert mm to pixel. The value should be equal to
        # 1/mask_fits['MM_PER_ARCSEC']/(platescale * bin_spat), but for some reason it's not,
        # and it's also different for the two detectors
        mm_pixel = 24.555832 if det == 1 else 24.548194

        left_edges = (mask_fits['POLY_Y'][0][targ] - mask_fits['MASK_CORNERS'][1])*mm_pixel + mask_edge_off
        right_edges = (mask_fits['POLY_Y'][1][targ] - mask_fits['MASK_CORNERS'][1])*mm_pixel + mask_edge_off
        if det == 2:
            # flip and reverse for detector 2
            Nx = self.get_rawimage(filename, det)[1].shape[1]
            left_edges, right_edges = Nx - right_edges, Nx - left_edges

        # Sort slits by their bottom edge position in ascending y-coordinate
        sortindx = np.argsort(left_edges)

        # Return the slit edges, sorted indices, and slitmask object
        return left_edges.astype(float), right_edges.astype(float), sortindx, self.slitmask

    def get_rawimage(self, raw_file, det):
        """
        Read raw images and generate a few other bits and pieces
        that are key for image processing.

        Parameters
        ----------
        raw_file : :obj:`str`
            File to read
        det : :obj:`int`
            1-indexed detector to read

        Returns
        -------
        detector_par : :class:`pypeit.images.detector_container.DetectorContainer`
            Detector metadata parameters.
        raw_img : `numpy.ndarray`_
            Raw image for this detector.
        hdu : `astropy.io.fits.HDUList`_
            Opened fits file
        exptime : :obj:`float`
            Exposure time read from the file header
        rawdatasec_img : `numpy.ndarray`_
            Data (Science) section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        oscansec_img : `numpy.ndarray`_
            Overscan section of the detector as provided by setting the
            (1-indexed) number of the amplifier used to read each detector
            pixel. Pixels unassociated with any amplifier are set to 0.
        """
        fil = utils.find_single_file(f'{raw_file}*', required=True)

        # Read
        log.info(f'Reading BINOSPEC file: {fil}')
        hdu = io.fits_open(fil)
        head1 = hdu[1].header

        # TOdO Store these parameters in the DetectorPar.
        # Number of amplifiers
        detector_par = self.get_detector_par(det if det is not None else 1, hdu=hdu)
        numamp = detector_par['numamplifiers']

        # get the x and y binning factors...
        binning = head1['CCDSUM']
        xbin, ybin = [int(ibin) for ibin in binning.split(' ')]

        # First read over the header info to determine the size of the output array...
        datasec = head1['DATASEC']
        x1, x2, y1, y2 = chain.from_iterable(parse.load_sections(datasec, fmt_iraf=False))
        nxb = x1 - 1

        # determine the output array size...
        nx = (x2 - x1 + 1) * int(numamp/2) + nxb * int(numamp/2)
        ny = (y2 - y1 + 1) * int(numamp/2)

        # allocate output array...
        array = np.zeros((nx, ny))
        rawdatasec_img = np.zeros_like(array, dtype=int)
        oscansec_img = np.zeros_like(array, dtype=int)

        if det == 1:  # A DETECTOR
            order = range(1, 5, 1)
        elif det == 2:  # B DETECTOR
            order = range(5, 9, 1)

        # insert extensions into calibration image...
        for kk, jj in enumerate(order):
            # grab complete extension...
            data, overscan, datasec, biassec = binospec_read_amp(hdu, jj)

            # insert components into output array...
            inx = data.shape[0]
            iny = data.shape[1]

            b1, b2, b3, b4 = chain.from_iterable(parse.load_sections(biassec, fmt_iraf=False))

            if kk == 0:
                array[b2:inx+b2,:iny] = data #*1.028
                rawdatasec_img[b2:inx+b2,:iny] = kk + 1
                array[:b2,:iny] = overscan
                oscansec_img[2:b2,:iny] = kk + 1
            elif kk == 1:
                array[b2+inx:2*inx+b2,:iny] = np.flipud(data) #* 1.115
                rawdatasec_img[b2+inx:2*inx+b2:,:iny] = kk + 1
                array[2*inx+b2:,:iny] = overscan
                oscansec_img[2*inx+b2:,:iny] = kk + 1
            elif kk == 2:
                array[b2+inx:2*inx+b2,iny:] = np.fliplr(np.flipud(data)) #* 1.047
                rawdatasec_img[b2+inx:2*inx+b2,iny:] = kk + 1
                array[2*inx+b2:, iny:] = overscan
                oscansec_img[2*inx+b2:, iny:] = kk + 1
            elif kk == 3:
                array[b2:inx+b2,iny:] = np.fliplr(data) #* 1.045
                rawdatasec_img[b2:inx+b2,iny:] = kk + 1
                array[:b2,iny:] = overscan
                oscansec_img[2:b2,iny:] = kk + 1

        # Need the exposure time
        exptime = hdu[self.meta['exptime']['ext']].header[self.meta['exptime']['card']]
        # Return, transposing array back to orient the overscan properly
        return detector_par, np.fliplr(np.flipud(array)), hdu, exptime, np.fliplr(np.flipud(rawdatasec_img)), \
               np.fliplr(np.flipud(oscansec_img))

    def bino_get_slit_region(self, filename, det=None, Nx=4096, Ny=4112, pady=0):
        """
        Compute the pixel-space rectangular regions for each slit in a Binospec mask.

        This function reads the slitmask design from a FITS file (or an already-loaded
        `SlitMask` object), converts slit and object positions from mask coordinates to
        pixel coordinates, and determines the x/y pixel boundaries for each slit on the
        detector. It returns these boundaries along with the updated slitmask object.

        Parameters
        ----------
        filename : :obj:`str`
            Path to the slitmask FITS file. Must be provided unless the slitmask
            is already loaded via `self.get_slitmask`.
        det : :obj:`int`, optional
            Detector number (1 or 2). Must be specified.
        Nx : :obj:`int`, optional
            Detector size in the x-direction (default: 4096 pixels).
        Ny : :obj:`int`, optional
            Detector size in the y-direction (default: 4112 pixels).
        pady : :obj:`float`, optional
            Additional padding (in pixels) applied to the slit boundaries (default: 0).

        Returns
        -------
        region : :obj:`list`
            A list containing:
            - slit_x_range : array of x-boundaries for each slit [Nslits, 2]
            - slit_y_range : array of y-boundaries for each slit [Nslits, 2]
            - x_slitobj_pix : array of x pixel positions for slit objects
            - y_slitobj_pix : array of y pixel positions for slit objects
        slitmask : :class:`SlitMask`
            The updated `SlitMask` object containing slit geometry and metadata.

        Notes
        -----
        - Converts mask coordinates to pixel coordinates using the appropriate scale factor.
        - Handles detector 2 by reversing slit order and applying a vertical flip.
        - Slit boundaries are clipped to remain within detector dimensions.
        """

        if det is None:
            raise ValueError("A valid detector number must be provided.")

        # Load slitmask information if a file is provided
        if filename is None:
            raise ValueError("The name of a science file should be provided")
        self.get_slitmask(filename, det=det)

        if self.slitmask is None:
            raise ValueError("Unable to read slitmask design info. Provide a file.")

        # Open FITS file and read mask data for the correct detector
        hdu = io.fits_open(filename)
        mask_fits = hdu[9].data[0] if det == 1 else hdu[10].data[0]
        numslits = len(self.slitmask.slitid)

        # Initialize arrays to hold slit x/y boundaries
        res_x = np.zeros((2, numslits))
        res_y = np.zeros((2, numslits))

        # Extract target distances from slit edges and slit widths
        topdist = np.asarray(self.slitmask.objects[:, 7], dtype=float)
        botdist = np.asarray(self.slitmask.objects[:, 8], dtype=float)
        width = np.asarray(self.slitmask.width)

        # Extract slit center positions in mask coordinates
        x_slits = np.asarray(self.slitmask.center[:, 0])
        x_obj = x_slits
        y_slits = -np.asarray(self.slitmask.center[:, 1])

        # Compute object y-position relative to slit center
        y_obj = y_slits + (topdist - botdist) / 2

        # Extract slit lengths and widths (in mask coordinates)
        dx_slits = self.slitmask.length
        dy_slits = width

        # Define scale factor and detector offsets
        dy0 = -200.0
        y_scl = 24.555832 if det == 1 else 24.548194

        # Extract mask corner reference point
        mask_corners = np.asarray(mask_fits['MASK_CORNERS'])
        corner_x = mask_corners[0]
        corner_y = mask_corners[1]

        # Convert slit center positions to pixel coordinates
        x_slits_pix = (x_slits - corner_x) * y_scl + Nx / 2.0
        x_slitobj_pix = (x_obj - corner_x) * y_scl + Nx / 2.0
        y_slits_pix = Ny - 1 - ((y_slits - corner_y) * y_scl) + dy0
        y_slitobj_pix = Ny - 1 - ((y_obj - corner_y) * y_scl) + dy0

        # Convert slit lengths and widths to pixel units
        dx_slits_pix = dx_slits * y_scl
        dy_slits_pix = dy_slits * y_scl

        # Loop through slits to compute pixel-space rectangular boundaries
        for i in range(numslits):
            xmin = round(x_slits_pix[i] - dx_slits_pix[i] / 2.0 - pady)
            xmax = round(x_slits_pix[i] + dx_slits_pix[i] / 2.0 - 1 + pady)
            res_x[0, i] = max(0, xmin)
            res_x[1, i] = min(Ny - 1, xmax)

            ymin = round(y_slits_pix[i] - dy_slits_pix[i] / 2.0 - pady)
            ymax = round(y_slits_pix[i] + dy_slits_pix[i] / 2.0 - 1 + pady)
            res_y[0, i] = max(0, ymin)
            res_y[1, i] = min(Ny - 1, ymax)

        # Handle detector 2: reverse slit order and flip vertically
        if det == 2:
            res_y = res_y[:, ::-1]

            # Apply vertical flip relative to detector height (Ny) and offset
            res_y_flipped = np.zeros_like(res_y)
            res_y_flipped[0, :] = -1 * (res_y[1, :] - Ny - 14)
            res_y_flipped[1, :] = -1 * (res_y[0, :] - Ny - 14)
            res_y = res_y_flipped

        # Package results and return
        slit_x_range, slit_y_range = res_x.T, res_y.T
        region = [slit_x_range, slit_y_range, x_slitobj_pix, y_slitobj_pix]

        return region, self.slitmask


    def plot_mask(self, filename, det=None, save_dir=None):
        """
        Plot the slit mask layout and target positions for one or both detectors.

        This function retrieves slit region data for a given Binospec mask and
        plots the rectangular slit outlines and target positions for detector 1,
        detector 2, or both. It is useful for visually validating mask design and
        target alignment.

        Parameters
        ----------
        filename : :obj:`str`
            Path to the mask design file (e.g., a JSON file containing slit definitions).
        det : :obj:`int` or :obj:`str`
            Specifies which detector(s) to plot. Accepts 1, 2, or 'both'.
        save_dir : :obj:`str`, optional
            If provided, the plot will be saved as a PNG in the given directory.

        Returns
        -------
        region_1 : :obj:`tuple`, optional
            Slit region and target position data for detector 1, if requested.
        region_2 : :obj:`tuple`, optional
            Slit region and target position data for detector 2, if requested.
        """

        if det is None:
            raise ValueError("A valid detector number must be provided: 1, 2, or 'both'")

        if filename is None:
            raise ValueError("A valid filename must be provided.")

        # Build save filename from FITS header
        hdu = io.fits_open(filename)
        basename = Path(filename).name
        save_filename = Path(f"plot_mask_{hdu[1].header['MASK']}_{basename}").with_suffix('.png')

        plt.rcParams.update({"font.size": 20})

        # Load slit regions depending on the selected detector(s)
        if det == 'both':
            fig, (axA, axB) = plt.subplots(ncols=2, figsize=(16, 16))
            region_1 = self.bino_get_slit_region(filename, det=1)[0]
            region_2 = self.bino_get_slit_region(filename, det=2)[0]

        elif det == 1:
            fig, axA = plt.subplots(figsize=(8, 8))
            region_1 = self.bino_get_slit_region(filename, det=1)[0]

        elif det == 2:
            fig, axB = plt.subplots(figsize=(8, 8))
            region_2 = self.bino_get_slit_region(filename, det=2)[0]

        else:
            raise ValueError("det must be 1, 2, or 'both'.")

        # Plot based on detector selection
        if det == 'both':
            _plot_region(axA, region_1, color="red", side_label="1")
            _plot_region(axB, region_2, color="green", side_label="2")
        elif det == 1:
            _plot_region(axA, region_1, color="red", side_label="1")
        elif det == 2:
            _plot_region(axB, region_2, color="green", side_label="2")

        # Save to file if directory provided
        if save_dir is not None:
            _save_dir = Path(save_dir).absolute()
            _save_dir.mkdir(parents=True, exist_ok=True)
            plt.tight_layout()
            plt.savefig(_save_dir / save_filename)
            plt.close(fig)
        else:
            plt.tight_layout()
            plt.show()

        # Return the plotted region data
        if det == 'both':
            return region_1, region_2
        elif det == 1:
            return region_1
        elif det == 2:
            return region_2


# Internal helper to draw slits and targets on a given axis
def _plot_region(ax, region, color, side_label):
    num_targets = len(region[0])
    label = f" N = {num_targets}"

    for i in range(len(region[0])):
        slit_x_range = region[0][i]
        slit_y_range = region[1][i]
        width = slit_x_range[1] - slit_x_range[0]
        height = slit_y_range[1] - slit_y_range[0]

        rect = patches.Rectangle(
            (slit_x_range[0], slit_y_range[0]),
            width,
            height,
            linewidth=1,
            edgecolor="blue",
            facecolor="none"
        )
        ax.add_patch(rect)

    ax.scatter(region[2], region[3], s=10, color=color, label=label)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_title(f"Detector {side_label}")
    ax.set_aspect("equal")
    ax.grid(True)
    ax.legend()


def binospec_read_amp(inp, ext):
    """
    Read one amplifier of an MMT BINOSPEC multi-extension FITS image

    Parameters
    ----------
    inp : str, :class:`astropy.io.fits.HDUList`
        The input FITS file name or already opened HDU list.
    ext : :obj:`int`
        FITS extension to read

    Returns
    -------
    data : :class:`numpy.ndarray`
        Array with data from the data section of the image.    
    overscan : :class:`numpy.ndarray`
        Array with the overscan section of the image.
    datasec : :obj:`str`
        String with the data section in IRAF format, e.g. '[x1:x2,y1:y2]'.
    biassec : :obj:`str`
        String with the bias section in IRAF format, e.g. '[x1:x2,y1:y2]'.
    """
    # Parse input
    hdu = io.fits_open(inp) if isinstance(inp, str) else inp

    # get entire extension...
    temp = hdu[ext].data.transpose()
    nxt = temp.shape[0]
    nyt = temp.shape[1]

    # parse the DETSEC keyword to determine the size of the array.
    header = hdu[ext].header

    # parse the DATASEC keyword to determine the size of the science region (unbinned)
    datasec = header['DATASEC']

    x1, x2, y1, y2 = chain.from_iterable(parse.load_sections(datasec, fmt_iraf=False))
    datasec = f'[{x1-1}:{x2},{y1-1}:{y2}]'

    # Overscan subtraction following IDL pipeline (bino_mosaic.pro):
    # Y-axis first, then X-axis. Uses sigma-clipped mean (resistant_mean)
    # with outlier cleaning, matching IDL defaults (clean_w=9, clean_nsig=1.0).

    # Y-axis overscan: postscan rows after datasec
    if y2 < nyt:
        overscan_y = temp[:, y2:nyt]
        overscan_vec, _, _ = sigma_clipped_stats(overscan_y, sigma=3.0, axis=1)
        overscan_vec = procimg.clean_overscan_vector(overscan_vec, w=9, nsig=1.0)
        temp = temp - overscan_vec[:, None]

    # X-axis overscan: prescan + postscan columns
    overscan_x_regions = []
    if x1 > 1:
        overscan_x_regions.append(temp[0:x1-1, :])
    if x2 < nxt:
        overscan_x_regions.append(temp[x2:nxt, :])
    if len(overscan_x_regions) > 0:
        overscan_x = np.concatenate(overscan_x_regions, axis=0)
        overscan_x_vec, _, _ = sigma_clipped_stats(overscan_x, sigma=3.0, axis=0)
        overscan_x_vec = procimg.clean_overscan_vector(overscan_x_vec, w=9, nsig=1.0)
        temp = temp - overscan_x_vec[None, :]

    # Crop to datasec
    data = temp[x1-1:x2, y1-1:y2]

    # Apply per-amplifier nonlinearity correction (IDL: poly(im_cur, c_poly))
    data = np.polynomial.polynomial.polyval(
        data, MMTBINOSPECSpectrograph.nonlinearity_coeffs[ext - 1])

    # Fake overscan for PypeIt's general pipeline (effectively a no-op)
    biassec = f'[0:{x1-1},{y1-1}:{y2}]'
    xos1, xos2, yos1, yos2 = chain.from_iterable(parse.load_sections(biassec, fmt_iraf=False))
    overscan = np.zeros_like(temp[xos1:xos2, yos1:yos2])

    return data, overscan, datasec, biassec


class MMTBINOSPECIFUSpectrograph(MMTBINOSPECSpectrograph):
    """
    Child to handle MMT/BINOSPEC IFU specific code.

    The Binospec IFU is a fiber-fed integral field unit with a hexagonal
    lenslet array feeding ~360 fibers per side into the spectrograph.
    Each side has 40 dedicated sky fibers at the outermost ring of each
    sub-bundle (indices [0-7, 88-95, 176-183, 264-271, 352-359]).
    """
    name = 'mmt_binospec_ifu'
    pypeline = 'Fiber'
    supported = True

    # IFU fiber geometry constants
    # On-sky fiber pitch in arcsec (hexagonal lenslet array)
    ifu_fiber_pitch = 0.6
    # Number of fibers per side
    nfibers_a = 360
    nfibers_b = 356

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        Adds 'decker' to the parent keys so that IFU frames are not
        grouped with MOS frames in the same configuration.

        Returns:
            :obj:`list`: List of configuration keys.
        """
        return super().configuration_keys() + ['decker']

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        Extends the parent class metadata with IFU-specific fields
        required by the Fiber pipeline (atmospheric parameters
        for DAR correction).
        """
        super().init_meta()
        # IFU-specific metadata for Fiber pipeline
        self.meta['slitwid'] = dict(card=None, compound=True)
        self.meta['obstime'] = dict(card=None, compound=True, required=False)
        self.meta['pressure'] = dict(card=None, compound=True, required=False)
        self.meta['temperature'] = dict(card=None, compound=True, required=False)
        self.meta['humidity'] = dict(card=None, compound=True, required=False)
        self.meta['parangle'] = dict(card=None, compound=True, required=False)

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
        if meta_key == 'slitwid':
            # IFU fiber pitch on sky, converted to degrees for WCS
            return self.ifu_fiber_pitch / 3600.0
        elif meta_key == 'obstime':
            try:
                return Time(headarr[1]['DATE-OBS'])
            except KeyError:
                log.warning("Time of observation not in header")
                return None
        elif meta_key == 'pressure':
            # MMT at ~2600m elevation, typical pressure ~730 mbar
            try:
                return headarr[1]['PRESSURE']
            except KeyError:
                log.warning("Pressure not in header - using default for MMT "
                            "elevation (730 mbar)")
                return 730.0
        elif meta_key == 'temperature':
            try:
                return headarr[1]['TEMP']
            except KeyError:
                log.warning("Temperature not in header - using default (5 deg C)")
                return 5.0
        elif meta_key == 'humidity':
            try:
                return headarr[1]['HUMID']
            except KeyError:
                log.warning("Humidity not in header - using default (20%%)")
                return 20.0
        elif meta_key == 'parangle':
            try:
                # PypeIt stores parangle in radians (see pypeit.core.meta).
                return np.deg2rad(headarr[1]['PA'])
            except KeyError:
                log.warning("Parallactic angle not in header - using default (0)")
                return 0.0
        else:
            return super().compound_meta(headarr, meta_key)

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.

        Overrides the parent to ensure only IFU frames (MASK == 'IFU')
        are selected for this spectrograph.

        Args:
            ftype (:obj:`str`):
                Type of frame to check.
            fitstbl (`astropy.table.Table`_):
                The table with the metadata for one or more frames to check.
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type ``ftype``.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        # Use parent frame typing logic, then restrict to IFU frames only
        is_type = super().check_frame_type(ftype, fitstbl, exprng=exprng)
        is_ifu = np.array([d.strip().upper() == 'IFU' for d in fitstbl['decker']])
        return is_type & is_ifu

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # IFU science frame processing
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 1.5
        par['scienceframe']['process']['use_illumflat'] = False
        par['scienceframe']['process']['use_specillum'] = False
        par['scienceframe']['process']['spat_flexure_correct'] = False
        par['scienceframe']['process']['use_biasimage'] = False
        par['scienceframe']['process']['use_darkimage'] = False

        # FiberFindObjects creates one SpecObj per fiber from slit edges
        # (no peak detection needed) and handles sky subtraction in its
        # own run() method, so skip the second find and final global.
        par['reduce']['findobj']['skip_second_find'] = True
        par['reduce']['findobj']['skip_final_global'] = True

        # Sky subtraction: use joint fit across all fibers
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['joint_fit'] = True
        # Avoid trimming edges of narrow IFU fibers (~5-6 pixels wide)
        par['reduce']['trim_edge'] = [0, 0]

        # Slit edge parameters: Sobel edge detection is bypassed for this
        # spectrograph via get_block_slit_edges(), which defines block-slit
        # edges directly from the reference fiber profile. Only pad is used
        # (passed through to SlitTraceSet).
        par['calibrations']['slitedges']['pad'] = 0

        # Scattered light correction for ALL frames (IDL pipeline approach).
        # Inter-block gaps contain ~3000+ counts of scattered light that must
        # be removed from both flats and science frames for accurate
        # fiber throughput ratios (Fabricant et al. 2025, Section 6).
        par['calibrations']['scattlight_pad'] = 5
        par['calibrations']['pixelflatframe']['process']['subtract_scattlight'] = True
        par['calibrations']['pixelflatframe']['process']['scattlight']['method'] = 'gaps'
        par['scienceframe']['process']['subtract_scattlight'] = True
        par['scienceframe']['process']['scattlight']['method'] = 'gaps'

        # Flat field: no edge tweaking for fiber-fed IFU (fixed positions)
        par['calibrations']['flatfield']['tweak_slits'] = False
        par['calibrations']['flatfield']['slit_trim'] = 0
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False

        # Flexure: Binospec has active flexure control, so spectral
        # flexure correction is not needed for IFU mode
        par['flexure']['spec_method'] = 'skip'

        # Flux calibration: extinction correction is done during datacube
        # construction, not during 1D extraction
        par['sensfunc']['UVIS']['extinct_correct'] = False

        return par

    def config_specific_par(
            self,
            inp:str|list|Path|fits.Header|Table,
            inp_par:parset.ParSet|None=None
        ) -> parset.ParSet:
        """
        Modify the PypeIt parameters to hard-wired values used for
        specific instrument configurations.

        Args:
            inp: Input filename, header, or metadata table row.
            inp_par: Parameter set. If None, use default.

        Returns:
            :class:`~pypeit.par.parset.ParSet`: Adjusted parameters.
        """
        par = super().config_specific_par(inp, inp_par=inp_par)

        # Grating-dependent bspline spacing for sky subtraction
        # (adopted from IDL pipeline knot spacings)
        grating = self.get_meta_value(inp, 'dispname')
        match grating:
            case 'x270':
                par['reduce']['skysub']['bspline_spacing'] = 1.05
            case 'x600':
                par['reduce']['skysub']['bspline_spacing'] = 0.5
            case 'x1000':
                par['reduce']['skysub']['bspline_spacing'] = 0.35

        # Override MOS-specific settings from parent's config_specific_par
        # that are inappropriate for IFU fibers
        par['calibrations']['slitedges']['use_maskdesign'] = False
        par['calibrations']['slitedges']['minimum_slit_length'] = 0.5
        par['calibrations']['slitedges']['edge_thresh'] = 5.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['reduce']['slitmask']['assign_obj'] = False
        par['reduce']['slitmask']['extract_missing_objs'] = False

        # Restore the IFU tilt model after the MOS parent overrides it. The IFU fibers
        # are lined up much more like the long-slit case. The fiber blocks share roughly
        # the same wavelength solution, but with higher order effects due to optical distortion
        # and mechanical placement of the fiber heads.
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev2d'] = 0.5
        par['calibrations']['wavelengths']['n_final'] = 5

        return par

    def get_wcs(self, hdr, slits, platescale, wave0, dwv, spatial_scale=None):
        """
        Construct a World-Coordinate System for the IFU datacube.

        Args:
            hdr (`astropy.io.fits.Header`_):
                The header of the raw frame.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces.
            platescale (:obj:`float`):
                The platescale of an unbinned pixel in arcsec/pixel.
            wave0 (:obj:`float`):
                The wavelength zeropoint.
            dwv (:obj:`float`):
                Change in wavelength per spectral pixel.
            spatial_scale (:obj:`float`, optional):
                User-specified spatial scale in arcsec.

        Returns:
            `astropy.wcs.WCS`_: The world-coordinate system.
        """
        log.info("Calculating the WCS for Binospec IFU")

        # Get binning
        binspec, binspat = parse.parse_binning(self.get_meta_value([hdr], 'binning'))

        # Spatial scales
        pxscl = platescale * binspat / 3600.0  # arcsec -> degrees
        slscl = self.get_meta_value([hdr], 'slitwid')  # already in degrees

        if spatial_scale is not None:
            pxscl = spatial_scale / 3600.0

        # Typical slit length
        slitlength = int(np.round(np.median(slits.get_slitlengths(median=True))))

        # Pointing coordinates
        raval = self.get_meta_value([hdr], 'ra')
        decval = self.get_meta_value([hdr], 'dec')
        coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))

        # Position angle from POSANG header keyword
        posang = hdr.get('POSANG', 0.0)
        crota = np.radians(-posang)

        # CD matrix
        cdelt1 = -slscl
        cdelt2 = pxscl
        cd11 = cdelt1 * np.cos(crota)
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)
        cd22 = cdelt2 * np.cos(crota)

        # Reference pixels (center of FOV)
        nslits = slits.nslits
        crpix1 = nslits / 2.0
        crpix2 = slitlength / 2.0
        crpix3 = 1.0

        # Create WCS
        log.info("Generating Binospec IFU WCS")
        w = wcs.WCS(naxis=3)
        w.wcs.equinox = hdr.get('EQUINOX', 2000.0)
        w.wcs.name = 'Binospec IFU'
        w.wcs.radesys = 'ICRS'
        w.wcs.cname = ['RA', 'DEC', 'Wavelength']
        w.wcs.cunit = [units.degree, units.degree, units.Angstrom]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN", "WAVE"]
        w.wcs.crval = [coord.ra.degree, coord.dec.degree, wave0]
        w.wcs.crpix = [crpix1, crpix2, crpix3]
        w.wcs.cd = np.array([[cd11, cd12, 0.0],
                             [cd21, cd22, 0.0],
                             [0.0, 0.0, dwv]])
        w.wcs.lonpole = 180.0
        w.wcs.latpole = 0.0

        return w

    def get_datacube_bins(self, slitlength, minmax, num_wave):
        r"""
        Calculate the bin edges to be used when making a datacube.

        Args:
            slitlength (:obj:`int`):
                Length of the slit in pixels.
            minmax (`numpy.ndarray`_):
                An array with the minimum and maximum pixel locations on
                each slit relative to the reference location. Shape must
                be :math:`(N_{\rm slits},2)`.
            num_wave (:obj:`int`):
                Number of wavelength steps.

        Returns:
            :obj:`tuple`: Three 1D `numpy.ndarray`_ providing the
            :math:`(x,y,\lambda)` bins for datacube construction.
        """
        # Number of fiber traces (slits) is determined from minmax
        nslits = minmax.shape[0]
        ref_slit = nslits // 2
        xbins = np.arange(1 + nslits) - ref_slit - 0.5
        ybins = np.linspace(np.min(minmax[:, 0]), np.max(minmax[:, 1]),
                            1 + slitlength) - 0.5
        spec_bins = np.arange(1 + num_wave) - 0.5
        return xbins, ybins, spec_bins

    @staticmethod
    def _ifu_calib_path() -> PypeItDataPath:
        """Return the path to the IFU calibration data directory."""
        return dataPaths.spectrographs / 'mmt_binospec'

    def load_fiber_ref_profile(self, det: int) -> fits.FITS_rec:
        """
        Load the reference fiber trace profile for fiber identification.

        The reference profile contains the expected pixel positions and
        Gaussian-Hermite profile parameters for each fiber, obtained
        from a high-quality flat field observation. This is used to
        cross-match detected fiber traces against known fiber IDs.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).

        Returns:
            `astropy.io.fits.FITS_rec`_: Table with columns:
                FIB_ID, X, Y, SIDE, FIB_NAME, FIB_TYPE, FIB_BLOCK,
                FIB_DEAD_FLAG, TR_A0, TR_PIX, TR_SIGMA, TR_BGR,
                TR_H3, TR_H4, TR_H5, TR_H6.
        """
        if det not in (1, 2):
            raise PypeItError(f'Detector number must be 1 or 2 for MMT/Binospec, not {det}.')
        ref_file = self._ifu_calib_path() / 'fiber_ref_profile.fits'
        # ext 1 = IFUTRACES_A (side A, det 1), ext 2 = IFUTRACES_B (side B, det 2)
        return fits.getdata(ref_file, ext=det)

    def load_sky_layout(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Load the IFU fiber-to-sky position mapping.

        Returns the on-sky x,y positions (in arcsec) for all 640 fibers
        in the hexagonal IFU field of view.

        Returns:
            :obj:`tuple`:
                - targetx_asec (`numpy.ndarray`_): x positions in arcsec (640,)
                - targety_asec (`numpy.ndarray`_): y positions in arcsec (640,)
        """
        sky_file = self._ifu_calib_path() / 'bino_IFU_sky_layout.fits'
        with fits.open(sky_file) as hdu:
            data = hdu[1].data[0]
            targetx = data['TARGETX_ASEC'].copy()
            targety = data['TARGETY_ASEC'].copy()
        return targetx, targety

    def load_fiber_illumination(self, det: int) -> np.ndarray:
        """
        Load the fiber-to-fiber illumination correction (throughput map).

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).

        Returns:
            `numpy.ndarray`_: Relative illumination correction per fiber (nfibers,).
        """
        illum_file = self._ifu_calib_path() / 'fiber_illumination.fits'
        # Row 0 = side A, row 1 = side B
        row = 0 if det == 1 else 1
        with fits.open(illum_file) as hdu:
            f_illum = hdu[1].data['F_ILLUM'][row].copy()
        # The IDL flat-fielding code applies this vector directly to the
        # extracted side-B fiber axis.  PypeIt's DET02 detector image is
        # mirror-flipped relative to that axis, so map the vector back onto
        # the reference-profile row order used by MASKDEF_ID lookup.
        return f_illum[::-1] if det == 2 else f_illum

    def get_ifu_datacube_meta(self, raw_hdr):
        """
        Return datacube header/WCS metadata for the Binospec fiber IFU.

        Binospec has a single fiber-IFU mode, so the returned label is fixed.
        See
        :meth:`pypeit.spectrographs.spectrograph.Spectrograph.get_ifu_datacube_meta`.

        Parameters
        ----------
        raw_hdr : `astropy.io.fits.Header`_
            Primary header of the input file (unused; Binospec has one mode).

        Returns
        -------
        :obj:`dict`
            ``{'name': 'BINOSPEC IFU', 'mode': 'FIBER'}``.
        """
        return {'name': f'{self.camera} IFU', 'mode': 'FIBER'}

    @staticmethod
    def ifu_sky_wcs(raw_hdr, scale_arcsec):
        """
        Build the celestial reference coordinate and CD matrix for the IFU.

        Encapsulates the single TAN/POSANG sign convention shared by the
        datacube builder and the 1D fiber extractor so the two stay in sync.
        The returned 2x2 CD matrix maps a (+x east, +y north) instrument
        offset of ``scale_arcsec`` arcsec per unit step to RA/Dec degrees.

        Parameters
        ----------
        raw_hdr : `astropy.io.fits.Header`_
            Primary header providing ``RA``, ``DEC`` and (optionally)
            ``POSANG``.
        scale_arcsec : :obj:`float`
            Spatial scale in arcsec per unit step along the WCS axes.

        Returns
        -------
        coord : `astropy.coordinates.SkyCoord`_
            Reference pointing (the WCS ``crval``).
        cd : `numpy.ndarray`_
            2x2 CD matrix ``[[cd11, cd12], [cd21, cd22]]`` in degrees.
        """
        raval = raw_hdr.get('RA', 0.0)
        decval = raw_hdr.get('DEC', 0.0)
        try:
            coord = SkyCoord(raval, decval, unit=(units.hourangle, units.deg))
        except Exception:
            coord = SkyCoord(raval, decval, unit=(units.deg, units.deg))

        crota = np.radians(-raw_hdr.get('POSANG', 0.0))
        cdelt1 = -scale_arcsec / 3600.0  # RA decreases with +x (east is +RA)
        cdelt2 = scale_arcsec / 3600.0   # Dec increases with +y

        cd11 = cdelt1 * np.cos(crota)
        cd12 = abs(cdelt2) * np.sign(cdelt1) * np.sin(crota)
        cd21 = -abs(cdelt1) * np.sign(cdelt2) * np.sin(crota)
        cd22 = cdelt2 * np.cos(crota)
        return coord, np.array([[cd11, cd12], [cd21, cd22]])

    def get_block_slit_edges(self, traceimg, det):
        """
        Define block-slit edges from the reference fiber profile.

        Instead of using Sobel edge detection (which fails due to scattered
        light in inter-block gaps), this method defines slit edges at the
        midpoints between adjacent fiber blocks. A bulk pixel shift is
        determined by cross-correlating the trace image against the expected
        fiber pattern.

        Args:
            traceimg (`numpy.ndarray`_):
                Trace image (raw flat), shape ``(nspec, nspat)``.
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :obj:`tuple`: ``(left_edges, right_edges)`` arrays of shape
                ``(nspec, nblocks)`` with constant slit edge positions.
        """
        blocks = self.get_fiber_blocks(det)
        nspec, nspat = traceimg.shape

        # Build expected spatial profile from reference fiber positions as a
        # sum of pixelated Gaussians: each fiber contributes the integral of
        # its Gaussian over the pixel (via erf), not the Gaussian sampled at
        # the pixel center.  This is only a cross-correlation template for the
        # bulk shift below, so the absolute normalization is irrelevant.
        all_positions = np.concatenate([b['fiber_positions'] for b in blocks])
        sigma = 2.5  # typical fiber sigma in pixels
        x = np.arange(nspat)
        # Pixel i spans [i-0.5, i+0.5]; integral of a unit-area Gaussian is
        # 0.5*(erf((hi-pos)/(sigma*sqrt2)) - erf((lo-pos)/(sigma*sqrt2))).
        edge = (x[None, :] - all_positions[:, None]) / (sigma * np.sqrt(2.0))
        ref_profile = 0.5 * (erf(edge + 0.5 / (sigma * np.sqrt(2.0)))
                             - erf(edge - 0.5 / (sigma * np.sqrt(2.0)))).sum(axis=0)

        # Collapse trace image to spatial profile (median of central half)
        obs_profile = np.median(traceimg[nspec // 4:3 * nspec // 4, :], axis=0)

        # Cross-correlate to find bulk shift (search within +/-20 pixels)
        max_shift = 20
        cc = correlate(obs_profile, ref_profile, mode='full')
        # The peak of cc is at index (len(obs) - 1 + shift)
        mid = len(obs_profile) - 1
        search_range = cc[mid - max_shift:mid + max_shift + 1]
        best_offset = np.argmax(search_range) - max_shift

        log.info(f"DET{det:02d}: bulk fiber shift = {best_offset} pixels "
                 f"(cross-correlation)")

        # Compute block slit edges at midpoints between blocks
        nblocks = len(blocks)
        left_edges = np.zeros((nspec, nblocks))
        right_edges = np.zeros((nspec, nblocks))

        for i, block in enumerate(blocks):
            shifted_positions = block['fiber_positions'] + best_offset
            block_min = shifted_positions.min()
            block_max = shifted_positions.max()

            # Typical fiber spacing within block
            if len(shifted_positions) > 1:
                fiber_spacing = np.median(np.diff(shifted_positions))
            else:
                fiber_spacing = 6.6  # fallback

            # Left edge: midpoint to previous block, or extend by fiber_spacing
            if i > 0:
                prev_max = blocks[i - 1]['fiber_positions'].max() + best_offset
                left_edges[:, i] = (prev_max + block_min) / 2.0
            else:
                left_edges[:, i] = block_min - fiber_spacing

            # Right edge: midpoint to next block, or extend by fiber_spacing
            if i < nblocks - 1:
                next_min = blocks[i + 1]['fiber_positions'].min() + best_offset
                right_edges[:, i] = (block_max + next_min) / 2.0
            else:
                right_edges[:, i] = block_max + fiber_spacing

            # Clip to detector bounds
            left_edges[:, i] = np.clip(left_edges[:, i], 0, nspat - 1)
            right_edges[:, i] = np.clip(right_edges[:, i], 0, nspat - 1)

        return left_edges, right_edges

    def get_arc_extract_center(self, slitcen, slits, det):
        """
        Snap arc extraction center to the nearest fiber in each block.

        The default ``slitcen`` is the midpoint of the block-slit edges,
        which may fall in an inter-fiber gap.  A gap-centered extraction
        with the default 3-pixel boxcar yields a noisy arc spectrum that
        degrades the wavelength solution for all fibers in the block.

        This method shifts each block's extraction center to the reference
        fiber position closest to the geometric center.

        Parameters
        ----------
        slitcen : `numpy.ndarray`_
            Slit center traces, shape ``(nspec, nslits)``.
        slits : :class:`~pypeit.slittrace.SlitTraceSet`
            Slit traces.
        det : :obj:`int`
            1-indexed detector number.

        Returns
        -------
        `numpy.ndarray`_
            Adjusted slit center traces, same shape as ``slitcen``.
        """
        nslits = slitcen.shape[1]
        blocks = self._require_block_slit_match(det, nslits)

        # Bulk shift between reference profile and detected slit positions
        mid = slitcen.shape[0] // 2
        slit_centers_mid = slitcen[mid, :]
        ref_centers = np.array([0.5 * (b['min_pix'] + b['max_pix'])
                                for b in blocks])
        shift = np.median(slit_centers_mid - ref_centers)

        adjusted = slitcen.copy()
        for i, block in enumerate(blocks):
            fiber_pos = block['fiber_positions'] + shift
            center = slit_centers_mid[i]
            nearest_idx = np.argmin(np.abs(fiber_pos - center))
            offset = fiber_pos[nearest_idx] - center
            adjusted[:, i] += offset
            log.info(f"Block {block['block_id']}: arc extract center "
                     f"shifted {offset:+.1f} px to nearest fiber")

        return adjusted

    def get_fiber_position_shift(self, slits, det):
        """
        Measure the bulk shift between reference fiber positions and the
        observed slit traces.

        The Binospec IFU slit definitions are built from the static reference
        fiber profile, shifted to match the observed trace image.  Any later
        fiber extraction must apply the same shift to the reference fiber
        centers; otherwise apertures are centered on the unshifted reference
        positions.

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Observed block-slit traces.
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).

        Returns:
            :obj:`float`: Bulk spatial shift in detector pixels.
        """
        blocks = self._require_block_slit_match(det, slits.nslits)

        mid_row = slits.nspec // 2
        left, right, _ = slits.select_edges()
        slit_centers = 0.5 * (left[mid_row] + right[mid_row])
        ref_centers = np.array([0.5 * (b['min_pix'] + b['max_pix'])
                                for b in blocks])
        shift = float(np.median(slit_centers - ref_centers))
        log.info(f"DET{det:02d}: fiber position shift from slits = "
                 f"{shift:.2f} px")
        return shift

    def adjust_slit_edges_to_fibers(self, slits, det):
        """
        Shrink slit edges to tightly wrap fiber positions from the reference
        profile, exposing inter-block gaps for scattered light modeling.

        The edge detection places slit boundaries at the midpoints of inter-block
        gaps (~70 px wide), consuming all off-slit pixels.  This method moves
        the edges inward to the outermost fiber positions + a small margin,
        leaving ~55-60 px gaps between blocks for the scattered light model.

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit traces to modify **in place**.
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).
        """
        margin = 7  # pixels beyond outermost fiber center

        blocks = self._require_block_slit_match(det, slits.nslits)

        # Determine the bulk shift between reference and detected positions.
        mid_row = slits.nspec // 2
        left, right, _ = slits.select_edges()
        shift = self.get_fiber_position_shift(slits, det)
        log.info(f"Slit edge adjustment: ref→detector shift = {shift:.1f} px")

        # Build new edge arrays (constant across spectral rows)
        new_left = np.zeros(slits.nslits)
        new_right = np.zeros(slits.nslits)
        for i, block in enumerate(blocks):
            new_left[i] = block['min_pix'] + shift - margin
            new_right[i] = block['max_pix'] + shift + margin

        # Ensure edges don't extend beyond detector
        new_left = np.clip(new_left, 0, None)

        # Log the change
        old_gaps = np.array([left[mid_row, i+1] - right[mid_row, i]
                             for i in range(slits.nslits - 1)])
        new_gaps = np.array([new_left[i+1] - new_right[i]
                             for i in range(slits.nslits - 1)])
        log.info(f"Inter-block gaps: old median={np.median(old_gaps):.0f} px, "
                 f"new median={np.median(new_gaps):.0f} px")

        # Apply to all spectral rows (edges are constant for fiber IFU)
        for row in range(slits.nspec):
            slits.left_init[row, :] = new_left
            slits.right_init[row, :] = new_right
        # Update tweaks if they exist, otherwise they'll be None
        # and select_edges() will fall back to left_init/right_init
        if slits.left_tweak is not None:
            slits.left_tweak[:] = slits.left_init
        if slits.right_tweak is not None:
            slits.right_tweak[:] = slits.right_init

    def subtract_scattered_light_gaps(self, image, offslitmask):
        """
        Subtract scattered light by measuring signal in inter-block gaps
        and interpolating across fiber blocks.

        For each spectral bin, measures the median signal in each off-slit
        gap region and linearly interpolates a smooth scattered light model
        across the spatial direction.

        Args:
            image (`numpy.ndarray`_):
                2D image (nspec, nspat) to measure scattered light from.
            offslitmask (`numpy.ndarray`_):
                Boolean mask, True for off-slit (gap) pixels.

        Returns:
            `numpy.ndarray`_: 2D scattered light model, same shape as image.
        """
        nspec, nspat = image.shape
        scatt_img = np.zeros_like(image)

        # Find contiguous off-slit gap regions from the mask at midpoint row
        mid_row = nspec // 2
        offslit = offslitmask[mid_row, :]

        # Identify gap boundaries
        edges = np.diff(offslit.astype(np.int8))
        gap_starts = np.where(edges == 1)[0] + 1
        gap_ends = np.where(edges == -1)[0] + 1

        # Handle cases where the image starts/ends in a gap
        if offslit[0]:
            gap_starts = np.concatenate(([0], gap_starts))
        if offslit[-1]:
            gap_ends = np.concatenate((gap_ends, [nspat]))

        n_gaps = min(len(gap_starts), len(gap_ends))
        if n_gaps < 2:
            log.warning("Fewer than 2 inter-block gaps found; "
                        "cannot model scattered light from gaps")
            return scatt_img

        gap_starts = gap_starts[:n_gaps]
        gap_ends = gap_ends[:n_gaps]
        gap_mids = 0.5 * (gap_starts + gap_ends)

        log.info(f"Scattered light from gaps: {n_gaps} gaps, "
                 f"median width {np.median(gap_ends - gap_starts):.0f} px")

        # Process in spectral bins for efficiency
        bin_size = 32
        n_bins = (nspec + bin_size - 1) // bin_size

        for ibin in range(n_bins):
            row_lo = ibin * bin_size
            row_hi = min(row_lo + bin_size, nspec)

            # Median across the spectral bin
            strip = np.median(image[row_lo:row_hi, :], axis=0)

            # Measure median signal in each gap (excluding edge pixels
            # that may contain fiber wing flux)
            gap_vals = np.zeros(n_gaps)
            for ig in range(n_gaps):
                gs, ge = int(gap_starts[ig]), int(gap_ends[ig])
                width = ge - gs
                margin = max(3, width // 5)
                inner = strip[gs + margin:ge - margin]
                if len(inner) > 0:
                    gap_vals[ig] = np.median(inner)
                else:
                    gap_vals[ig] = np.median(strip[gs:ge])

            # Interpolate between gap midpoints (linear, constant beyond)
            scatt_profile = np.interp(np.arange(nspat), gap_mids, gap_vals)

            # Apply to all rows in the bin
            scatt_img[row_lo:row_hi, :] = scatt_profile[np.newaxis, :]

        # Smooth along spectral direction to remove binning steps
        scatt_img = uniform_filter1d(scatt_img, size=bin_size, axis=0,
                                     mode='nearest')

        log.info(f"Scattered light model: "
                 f"median={np.median(scatt_img):.0f}, "
                 f"range=[{np.min(scatt_img):.0f}, {np.max(scatt_img):.0f}]")

        return scatt_img

    def get_fiber_blocks(self, det):
        """
        Return the fiber block structure from the reference profile.

        Each block is a group of fibers that will become a single "slit"
        in the block-slit extraction approach. Blocks are defined by the
        FIB_BLOCK column in the reference profile.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).

        Returns:
            :obj:`list` of :obj:`dict`: One dict per block with keys:
                - 'block_id': int, block number from reference profile
                - 'nfibers': int, number of fibers in block
                - 'type': str, 'sky' or 'science'
                - 'fiber_positions': ndarray, reference pixel positions (TR_PIX)
                - 'fiber_names': list of str, fiber names
                - 'fiber_ids': ndarray, fiber IDs
                - 'min_pix': float, minimum pixel position in block
                - 'max_pix': float, maximum pixel position in block
        """
        ref = self.load_fiber_ref_profile(det)
        # Exclude dead fibers assigned to block -1 (no valid trace position)
        live = ref['FIB_BLOCK'] >= 0
        ref = ref[live]
        block_ids = np.unique(ref['FIB_BLOCK'])
        blocks = []
        for bid in block_ids:
            mask = ref['FIB_BLOCK'] == bid
            names = [n.strip() for n in ref['FIB_NAME'][mask]]
            n_sky = sum(1 for n in names if n.startswith('SKY'))
            positions = ref['TR_PIX'][mask]
            blocks.append({
                'block_id': int(bid),
                'nfibers': int(np.sum(mask)),
                'type': 'sky' if n_sky > 0 else 'science',
                'fiber_positions': positions,
                'fiber_names': names,
                'fiber_ids': ref['FIB_ID'][mask],
                'min_pix': float(positions.min()),
                'max_pix': float(positions.max()),
            })
        return blocks

    def _require_block_slit_match(self, det, nslits):
        """
        Return the reference fiber blocks for ``det``, requiring that their
        count matches the number of traced block-slits.

        The block-slits are defined directly from the static reference fiber
        profile (one block-slit per ``FIB_BLOCK``), so the traced slit count
        should always equal the reference block count.  Every fiber operation
        (arc-center snapping, edge adjustment, throughput, sky identification)
        pairs the i-th reference block with the i-th block-slit positionally.
        A mismatch means a block-slit was dropped or merged during edge
        tracing/QA, which breaks that correspondence; continuing would
        silently misassign fibers, throughputs, and wavelength solutions, so
        the reduction is faulted instead.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).
            nslits (:obj:`int`):
                Number of traced block-slits.

        Returns:
            :obj:`list`: The fiber blocks from :func:`get_fiber_blocks`.

        Raises:
            :class:`~pypeit.PypeItError`:
                If the reference block count does not equal ``nslits``.
        """
        blocks = self.get_fiber_blocks(det)
        if len(blocks) != nslits:
            raise PypeItError(
                f"DET{det:02d}: reference fiber-block count ({len(blocks)}) "
                f"does not match the traced block-slit count ({nslits}).  A "
                f"block-slit was likely dropped or merged during edge "
                f"tracing, breaking the block-to-slit correspondence required "
                f"for fiber assignment.  Check the slit-edge tracing for this "
                f"detector.")
        return blocks

    def identify_fibers_in_block(self, det, block_idx, detected_positions):
        """
        Identify fibers within a block-slit by matching detected peak positions
        to reference fiber positions.

        Args:
            det (:obj:`int`): 1-indexed detector number.
            block_idx (:obj:`int`): 0-based block index.
            detected_positions (`numpy.ndarray`_): Detected fiber peak pixel
                positions within the block-slit, sorted by position.

        Returns:
            :obj:`dict`: Keys 'fiber_id', 'fiber_name', 'fiber_type' — arrays
                aligned with detected_positions. Unmatched fibers get
                fiber_id=-1, fiber_name='UNKNOWN', fiber_type='unknown'.
        """
        blocks = self.get_fiber_blocks(det)
        block = blocks[block_idx]
        ref_positions = block['fiber_positions']
        ref_ids = block['fiber_ids']
        ref_names = block['fiber_names']

        n_det = len(detected_positions)
        n_ref = len(ref_positions)

        fiber_id = np.full(n_det, -1, dtype=int)
        fiber_name = np.array(['UNKNOWN'] * n_det, dtype='U20')
        fiber_type = np.array(['unknown'] * n_det, dtype='U10')

        # Compute offset: median shift between detected and reference
        if n_det == n_ref:
            offset = np.median(detected_positions - ref_positions)
        else:
            offset = 0.0

        shifted_ref = ref_positions + offset

        # Match by nearest neighbor with threshold
        threshold = 3.5  # pixels
        for i, dpos in enumerate(detected_positions):
            dists = np.abs(shifted_ref - dpos)
            best = np.argmin(dists)
            if dists[best] < threshold:
                fiber_id[i] = int(ref_ids[best])
                fiber_name[i] = ref_names[best]
                fiber_type[i] = 'sky' if ref_names[best].startswith('SKY') else 'science'

        return {'fiber_id': fiber_id, 'fiber_name': fiber_name, 'fiber_type': fiber_type}

    def modify_pixelflat(self, flatimages, slits, det):
        """
        Override to skip baking fiber illumination into the pixel flat.

        In the block-slit extraction approach, throughput corrections
        are applied to extracted 1D spectra via ``FiberFlatImages``, not
        to the pixel flat.
        """
        pass

    def measure_fiber_flat_flux(self, flatimg, slits, det):
        """
        Measure integrated flat field flux for each fiber within block-slits.

        Used to compute the bulk throughput ratio between sky fibers (bare)
        and science fibers (lenslet-fed).

        Args:
            flatimg (`numpy.ndarray`_):
                Flat field image, shape ``(nspec, nspat)``.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Block-slit traces (21 per detector).
            det (:obj:`int`):
                1-indexed detector number.

        Returns:
            :obj:`dict`: Dictionary with keys:
                - 'fiber_flux': per-fiber integrated flat flux (nfibers,)
                - 'fiber_type': per-fiber type ('sky' or 'science')
                - 'sky_avg': mean flux of sky fibers
                - 'sci_avg': mean flux of science fibers
                - 'bulk_scale': sci_avg / sky_avg (scalar)
        """
        blocks = self._require_block_slit_match(det, slits.nslits)
        slitmask = slits.slit_img(pad=0)

        all_flux = []
        all_type = []

        for block_idx, block in enumerate(blocks):
            slit_spat_id = slits.spat_id[block_idx]
            thismask = slitmask == slit_spat_id
            if not np.any(thismask):
                continue

            for j, fpos in enumerate(block['fiber_positions']):
                # Boxcar integrate the flat flux around each fiber
                trace = np.full(flatimg.shape[0], fpos)
                # Half-spacing from neighbors
                if block['nfibers'] > 1:
                    positions = block['fiber_positions']
                    spacings = np.diff(positions)
                    if j == 0:
                        half_sp = spacings[0] / 2.0
                    elif j == len(positions) - 1:
                        half_sp = spacings[-1] / 2.0
                    else:
                        half_sp = min(spacings[j-1], spacings[j]) / 2.0
                else:
                    half_sp = 3.3  # single fiber fallback
                box_flux = moment1d(flatimg * thismask, trace, 2 * half_sp,
                                    row=np.arange(flatimg.shape[0]))[0]
                med_flux = np.median(box_flux[box_flux > 0]) if np.any(box_flux > 0) else 0.0
                all_flux.append(med_flux)
                all_type.append(block['type'])

        all_flux = np.array(all_flux)
        all_type = np.array(all_type)
        sky_mask = all_type == 'sky'
        sci_mask = all_type == 'science'

        sky_avg = np.median(all_flux[sky_mask]) if np.any(sky_mask) else 1.0
        sci_avg = np.median(all_flux[sci_mask]) if np.any(sci_mask) else 1.0
        bulk_scale = sci_avg / sky_avg if sky_avg > 0 else 1.0

        return {
            'fiber_flux': all_flux,
            'fiber_type': all_type,
            'sky_avg': sky_avg,
            'sci_avg': sci_avg,
            'bulk_scale': bulk_scale,
        }

    def get_science_fiber_layout_indices(self, det: int,
                                         fiber_ids: np.ndarray,
                                         fiber_types: np.ndarray) -> np.ndarray:
        """
        Map detected fibers to layout file indices using fiber IDs.

        Uses fiber IDs from :func:`get_fiber_metadata` to look up each
        fiber's name in the reference profile, then matches that name to
        the layout file entry.  This works correctly even when fibers are
        missing from the input data.

        The layout file (``bino_IFU_sky_layout.fits``) contains 640 entries
        (indices 0-319 for side A, 320-639 for side B).  Live science
        fibers from the reference profile are sorted by detector position
        and paired with live layout entries: in forward order for side A,
        and in reverse order for side B (because the two detectors produce
        mirror-image spectra).  Dead fibers (``_DEAD`` suffix) are excluded
        from both lists so they do not disrupt the pairing.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).
            fiber_ids (`numpy.ndarray`_):
                Physical fiber IDs for each detected fiber, from
                ``fiber_meta['fiber_id']``.  Unmatched fibers have ID < 0.
            fiber_types (`numpy.ndarray`_):
                Fiber type strings (``'SCI'``, ``'SKY'``, or ``'UNKNOWN'``),
                from ``fiber_meta['fiber_type']``.

        Returns:
            `numpy.ndarray`_: Array of shape ``(nfibers,)`` with layout file
            indices (0-639) for each fiber. Sky, dead, and unmatched fibers
            are assigned -1.
        """
        ref = self.load_fiber_ref_profile(det)
        ref_ids = ref['FIB_ID']
        ref_names = [n.strip() for n in ref['FIB_NAME']]
        ref_pix = ref['TR_PIX']

        # Identify live science fibers in the reference profile and sort
        # by detector position (TR_PIX).
        ref_is_sci = np.array([not n.startswith('SKY') for n in ref_names])
        ref_is_live = np.array(['_DEAD' not in n for n in ref_names])
        ref_live_sci = ref_is_sci & ref_is_live
        live_sci_ids = ref_ids[ref_live_sci]
        live_sci_pix = ref_pix[ref_live_sci]
        sort_idx = np.argsort(live_sci_pix)
        live_sci_ids_sorted = live_sci_ids[sort_idx]

        # Get live layout entries for this side in layout-file order.
        sky_file = self._ifu_calib_path() / 'bino_IFU_sky_layout.fits'
        with fits.open(sky_file) as hdu:
            layout_names = [n.strip() for n in hdu[1].data[0]['TARGET_NAME']]
        start, end = (0, 320) if det == 1 else (320, 640)
        live_layout_indices = [i for i in range(start, end)
                               if '_DEAD' not in layout_names[i]]

        # Pair reference fibers (by detector position) with layout entries.
        # Side A: forward (leftmost on detector = first layout entry).
        # Side B: reversed (leftmost on detector = last layout entry)
        # because the two detectors produce mirror-image spectra.
        n_map = min(len(live_sci_ids_sorted), len(live_layout_indices))
        id_to_layout = {}
        for rank in range(n_map):
            fid = int(live_sci_ids_sorted[rank])
            if det == 1:
                id_to_layout[fid] = live_layout_indices[rank]
            else:
                id_to_layout[fid] = live_layout_indices[n_map - 1 - rank]

        nfibers = len(fiber_ids)
        layout_indices = np.full(nfibers, -1, dtype=int)

        for i in range(nfibers):
            if fiber_types[i] == 'SKY' or fiber_ids[i] < 0:
                continue
            layout_idx = id_to_layout.get(int(fiber_ids[i]), -1)
            if layout_idx >= 0:
                layout_indices[i] = layout_idx

        n_mapped = int(np.sum(layout_indices >= 0))
        log.info(f"Mapped {n_mapped} science fibers to layout indices "
                 f"(det={det})")
        return layout_indices

    def _load_ref_spatial_profile(self, det):
        """
        Load the reference spatial profile for cross-correlation.

        The reference profile (PROF_REF) is a 1D array representing the
        summed Gaussian-Hermite profiles of all fibers, stored in
        extension 3 (IFUPROF) of the reference fiber profile FITS file.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).

        Returns:
            `numpy.ndarray`_: 1D reference spatial profile.
        """
        ref_file = self._ifu_calib_path() / 'fiber_ref_profile.fits'
        side_idx = 0 if det == 1 else 1
        with fits.open(ref_file) as hdu:
            return hdu[3].data['PROF_REF'][side_idx].copy()

    def match_fibers_to_reference(self, det, detected_positions):
        """
        Cross-match detected fiber trace positions against the reference
        profile to assign physical fiber IDs.

        Follows the algorithm from the IDL pipeline
        (bino_ifu_fiber_id.pro):

        1. Build a synthetic spatial profile from detected positions.
        2. Cross-correlate against the reference profile in 5 segments
           with cosine apodization.
        3. Fit a linear polynomial to the segment offsets to capture
           position-dependent shifts (flexure, scale, distortion).
        4. Match individual fibers using a distance threshold after
           applying the per-trace polynomial shift.

        Args:
            det (:obj:`int`):
                1-indexed detector number (1=side A, 2=side B).
            detected_positions (`numpy.ndarray`_):
                Detected fiber center positions in pixels at a reference
                column (e.g., center of detector).

        Returns:
            :obj:`tuple`:
                - fiber_ids (`numpy.ndarray`_): Physical fiber IDs for each
                  detected trace, or -1 if unmatched.
                - is_sky (`numpy.ndarray`_): Boolean array, True for sky fibers.
                - is_dead (`numpy.ndarray`_): Boolean array, True for dead fibers
                  in the reference that were not detected.
        """
        ref = self.load_fiber_ref_profile(det)
        ref_positions = ref['TR_PIX']
        ref_ids = ref['FIB_ID']
        ref_dead = ref['FIB_DEAD_FLAG'].astype(bool)
        ref_names = [n.strip() for n in ref['FIB_NAME']]

        # Load the reference spatial profile (from ext 3 of the FITS file)
        ref_profile = self._load_ref_spatial_profile(det)
        nx = len(ref_profile)

        # Build a synthetic spatial profile from detected positions
        # using Gaussians with typical fiber width (sigma ~ 1.25 px)
        sigma = 1.25
        xvec = np.arange(nx, dtype=float)
        det_profile = np.zeros(nx)
        for pos in detected_positions:
            ipos = int(np.round(pos))
            # Only compute Gaussian within +/- 5 sigma for efficiency
            hw = int(5 * sigma) + 1
            lo = max(0, ipos - hw)
            hi = min(nx, ipos + hw + 1)
            if lo < hi:
                det_profile[lo:hi] += np.exp(
                    -0.5 * ((xvec[lo:hi] - pos) / sigma) ** 2)

        # 5-segment cross-correlation (following IDL bino_ifu_fiber_id.pro)
        M = 80          # maximum lag (+/- pixels)
        n_seg = 5
        overlap = 0.1
        skippix1 = 150  # skip edge pixels
        skippix2 = 150

        x_seg = np.zeros(n_seg)
        c_seg = np.full(n_seg, np.nan)

        for iseg in range(n_seg):
            # Segment boundaries with overlap (matching IDL logic)
            seg_width = (nx - skippix1 - skippix2) / n_seg
            if iseg == 0:
                nmin = skippix1
            else:
                nmin = int(skippix1 + (iseg - overlap) * seg_width)
            nmin = max(nmin, 0)

            if iseg == n_seg - 1:
                nmax = nx - skippix2 - 1
            else:
                nmax = int(skippix1 + (iseg + 1.0 + overlap) * seg_width)
            nmax = min(nmax, nx - 1)

            x_seg[iseg] = (nmax + nmin) / 2.0
            c_seg[iseg] = self._segment_xcorr_offset(
                det_profile, ref_profile, nmin, nmax, M)

        # Fit linear polynomial to segment offsets
        valid = np.isfinite(c_seg)
        n_valid = np.sum(valid)
        if n_valid >= 2:
            poly_coeff = np.polyfit(x_seg[valid], c_seg[valid], 1)
        elif n_valid == 1:
            poly_coeff = np.array([0.0, c_seg[valid][0]])
        else:
            log.warning("All cross-correlation segments failed; "
                        "falling back to zero offset")
            poly_coeff = np.array([0.0, 0.0])

        # Two-pass iterative matching (following IDL bino_ifu_fiber_id.pro):
        # Pass 1: match with 1.7 px threshold using cross-correlation polynomial
        # Pass 2: refit polynomial from matched pairs, then match remaining
        #         fibers with relaxed 3.5 px threshold to catch physically
        #         displaced fibers (broken/swapped bundles on DET02)
        dx_thr_tight = 1.7
        dx_thr_relaxed = 3.5
        fiber_ids = np.full(len(detected_positions), -1, dtype=int)
        matched = np.zeros(len(ref_positions), dtype=bool)

        for pass_num, dx_thr in enumerate([dx_thr_tight, dx_thr_relaxed]):
            # Evaluate per-trace shift from the polynomial
            dx_all = np.polyval(poly_coeff, detected_positions)

            for i, dpos in enumerate(detected_positions):
                if fiber_ids[i] >= 0:
                    continue    # already matched
                # Distance = |detected_pos - shift - ref_pos|
                dists = np.abs(dpos - dx_all[i] - ref_positions)
                dists[matched] = np.inf
                dists[ref_dead] = np.inf    # exclude dead fibers
                best = np.argmin(dists)
                if dists[best] < dx_thr:
                    fiber_ids[i] = ref_ids[best]
                    matched[best] = True

            # After pass 1, refit polynomial using matched pairs for
            # a more accurate alignment before the relaxed pass
            if pass_num == 0:
                matched_det = detected_positions[fiber_ids >= 0]
                matched_ref_pos = np.array([
                    ref_positions[np.where(ref_ids == fid)[0][0]]
                    for fid in fiber_ids if fid >= 0])
                offsets = matched_det - matched_ref_pos
                if len(offsets) >= 2:
                    poly_coeff = np.polyfit(matched_det, offsets, 1)

        # Determine sky fiber status from reference fiber names
        # (the FIB_TYPE field is unreliable; use FIB_NAME instead)
        is_sky = np.zeros(len(detected_positions), dtype=bool)
        for i in range(len(detected_positions)):
            if fiber_ids[i] >= 0:
                idx = np.where(ref_ids == fiber_ids[i])[0]
                if len(idx) > 0:
                    is_sky[i] = ref_names[idx[0]].startswith('SKY')

        # Dead fibers: reference fibers that were not matched
        is_dead = np.logical_not(matched) & np.logical_not(ref_dead)

        n_matched = np.sum(fiber_ids >= 0)
        log.info(f"Fiber matching: {n_matched}/{len(detected_positions)} "
                 f"detected traces matched to reference "
                 f"({np.sum(is_sky)} sky fibers, "
                 f"{np.sum(is_dead)} dead fibers)")

        return fiber_ids, is_sky, is_dead

    @staticmethod
    def _segment_xcorr_offset(det_profile, ref_profile, nmin, nmax, max_lag,
                              ntaper=10, min_corr=0.3, dx_cpeak=4):
        """
        Cross-correlate one spatial segment and return its sub-pixel offset.

        Helper for :meth:`match_fibers_to_reference`, invoked once per
        segment.  The segment ``[nmin, nmax]`` of the detected and reference
        profiles is cosine-apodized at both ends, the two are cross-correlated
        over ``+/- max_lag`` pixels and normalized, and the correlation peak
        is refined by parabolic interpolation.  Mirrors the per-segment logic
        of the IDL pipeline (``bino_ifu_fiber_id.pro``).

        Args:
            det_profile (`numpy.ndarray`_):
                Synthetic spatial profile built from the detected fiber
                positions.
            ref_profile (`numpy.ndarray`_):
                Reference spatial profile.
            nmin, nmax (:obj:`int`):
                Inclusive segment boundary columns.
            max_lag (:obj:`int`):
                Maximum cross-correlation lag in pixels.
            ntaper (:obj:`int`, optional):
                Width of the cosine taper at each segment end.
            min_corr (:obj:`float`, optional):
                Minimum normalized correlation peak required to accept the
                segment.
            dx_cpeak (:obj:`int`, optional):
                Required margin (in lag bins) between the peak and the
                lag-window edge for sub-pixel interpolation.

        Returns:
            :obj:`float`: Sub-pixel lag offset of the correlation peak, or
            ``np.nan`` if the segment is rejected (weak correlation or peak
            too close to the lag-window edge).
        """
        # Cosine apodization window
        seg_len = nmax - nmin + 1
        c_ap = np.ones(seg_len)
        taper = 0.5 * (1 - np.cos(np.pi * np.arange(ntaper) / ntaper))
        c_ap[:ntaper] = taper
        c_ap[-ntaper:] = taper[::-1]

        # Apply apodization to absolute values (like IDL)
        seg_det = np.abs(det_profile[nmin:nmax + 1]) * c_ap
        seg_ref = np.abs(ref_profile[nmin:nmax + 1]) * c_ap

        # Cross-correlate and extract +/- max_lag lag range
        corr = np.correlate(seg_det, seg_ref, mode='full')
        n = len(seg_det)
        lag = np.arange(len(corr)) - (n - 1)
        lag_mask = np.abs(lag) <= max_lag
        corr_sub = corr[lag_mask]
        lag_sub = lag[lag_mask]

        # Normalize
        norm = np.sqrt(np.sum(seg_det ** 2) * np.sum(seg_ref ** 2))
        if norm > 0:
            corr_sub = corr_sub / norm

        max_idx = np.argmax(corr_sub)
        max_val = corr_sub[max_idx]

        # Quality check and sub-pixel peak fit (like IDL)
        if not (max_val >= min_corr
                and max_idx > dx_cpeak
                and max_idx < len(corr_sub) - dx_cpeak):
            return np.nan

        # Parabolic sub-pixel interpolation
        left = corr_sub[max_idx - 1]
        center = corr_sub[max_idx]
        right = corr_sub[max_idx + 1]
        denom = left - 2 * center + right
        delta = 0.5 * (left - right) / denom if denom != 0 else 0.0
        return lag_sub[max_idx] + delta

    def get_fiber_metadata(self, det, slit_spat_ids, slit_centers=None):
        """
        Map detected fiber traces to Binospec IFU fiber identifiers.

        Uses cross-correlation of detected trace positions against the
        reference fiber profile to assign physical fiber IDs and names.

        See base class for parameter and return value documentation.

        Parameters
        ----------
        slit_centers : `numpy.ndarray`_, optional
            Float-valued slit center positions at the spectral midpoint.
            If provided, these are used instead of the integer
            ``slit_spat_ids`` for more accurate fiber matching.
        """
        ref = self.load_fiber_ref_profile(det)

        # Use float centers when available for sub-pixel accuracy
        detected_positions = slit_centers if slit_centers is not None \
            else slit_spat_ids.astype(float)

        fiber_ids, is_sky, _ = self.match_fibers_to_reference(det, detected_positions)

        # Build name and type arrays from the reference profile
        nslits = len(slit_spat_ids)
        fiber_names = np.full(nslits, 'UNKNOWN', dtype='U20')
        fiber_types = np.full(nslits, 'UNKNOWN', dtype='U10')

        # Map matched fiber IDs back to reference table for names,
        # and derive type from FIB_NAME (SKY* = sky fiber, else science).
        # The FIB_TYPE field in the reference file is unreliable (all SKY).
        ref_ids = ref['FIB_ID']
        ref_names = [n.strip() for n in ref['FIB_NAME']]

        for i in range(nslits):
            if fiber_ids[i] >= 0:
                idx = np.where(ref_ids == fiber_ids[i])[0]
                if len(idx) > 0:
                    fiber_names[i] = ref_names[idx[0]]
                    fiber_types[i] = 'SKY' if ref_names[idx[0]].startswith('SKY') \
                        else 'SCI'

        n_matched = np.sum(fiber_ids >= 0)
        n_sky = np.sum(fiber_types == 'SKY')
        log.info(f"Fiber metadata: {n_matched}/{nslits} fibers identified "
                 f"({n_sky} sky, {n_matched - n_sky} science)")

        return {'fiber_id': fiber_ids,
                'fiber_name': fiber_names,
                'fiber_type': fiber_types}
