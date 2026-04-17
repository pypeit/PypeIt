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
from IPython import embed
import matplotlib.pyplot as plt
from matplotlib import patches
import numpy as np

from pypeit import dataPaths
from pypeit import io
from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit import utils
from pypeit.core import framematch
from pypeit.core import parse
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
        return ['dispname']

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
        return ['DISPERS1']

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
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['lampstat02'] == 'deployed')

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
            xs = inx * kk
            xe = xs + inx

            iny = data.shape[1]
            ys = iny * kk
            yn = ys + iny

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

        # Extract slit corner y-coordinates (for top/bottom edges)
        y_slitsh = -np.asarray(self.slitmask.corners[:, 0, 1])
        y_slitsl = -np.asarray(self.slitmask.corners[:, 2, 1])

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
        y_slitsl_pix = Ny - 1 - ((y_slitsl - corner_y) * y_scl) + dy0
        y_slitsh_pix = Ny - 1 - ((y_slitsh - corner_y) * y_scl) + dy0

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

    # NOTE: Since pypeit can only subtract overscan along one axis, I'm subtract
    # the overscan here using median method.
    # Overscan X-axis
    if x1 > 1:
        overscanx = temp[2:x1-1, :]
        overscanx_vec = np.median(overscanx, axis=0)
        temp = temp - overscanx_vec[None,:]
    data = temp[x1-1:x2, y1-1:y2]

    ## Overscan Y-axis
    if y2 < nyt:
        os1, os2 = y2+1, nyt-1
        overscany = temp[x1 - 1:x2, y2:os2]
        overscany_vec = np.median(overscany, axis=1)
        data = data -  overscany_vec[:,None]

    # Overscan
    biassec = f'[0:{x1-1},{y1-1}:{y2}]'
    xos1, xos2, yos1, yos2 = chain.from_iterable(parse.load_sections(biassec, fmt_iraf=False))
    overscan = np.zeros_like(temp[xos1:xos2, yos1:yos2]) # Give a zero fake overscan at the edge of each amplifiers

    return data, overscan, datasec, biassec

