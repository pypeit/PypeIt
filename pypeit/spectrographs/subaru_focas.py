"""
Module for Subaru FOCAS

.. include:: ../include/links.rst
"""
import numpy as np
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.core import meta
from pypeit.spectrographs import spectrograph
from pypeit import io
from pypeit.images import detector_container
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io import fits
from IPython import embed

class SubaruFOCASSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle Subaru/FOCAS specific code
    """
    ndet = 1  # Because each detector is written to a separate FITS file
    telescope = telescopes.SubaruTelescopePar()
    url = 'https://www.naoj.org/Observing/Instruments/FOCAS/index.html'

    name = 'subaru_focas'
    camera = 'FOCAS'
    header_name = 'FOCAS'
    supported = True
    comment = 'Supported grisms: 300B, 300R, VPH520, VPH850; see :doc:`subaru_focas`'


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Always correct for flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'

        # Median overscan
        #   IF YOU CHANGE THIS, YOU WILL NEED TO DEAL WITH THE OVERSCAN GOING ALONG ROWS
        turn_off_on = dict(overscan_method='median')
        par.reset_all_processimages_par(**turn_off_on)

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        # NOTE: after removing overscan regions, finding edges does not
        # work well and translate into slits--the error message suggested
        # to add this
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] = 25.0
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 4

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        #par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.07
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0  # Good for 2x binning
        par['calibrations']['wavelengths']['n_final'] = 4

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5
        #par['sensfunc']['IR']['telgridfile'] = 'TelFit_MaunaKea_3100_26100_R20000.fits'
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_10500_R120000.fits'

        # Frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        # Think we needed to change lower bound on this to recognize some arcs
        #par['calibrations']['arcframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [1, None]
        # Think we needed to change the upper bound on this for some standards
        #par['calibrations']['standardframe']['exprng'] = [1, None]
        par['calibrations']['standardframe']['exprng'] = [1, 61]
        #par['scienceframe']['exprng'] = [1, None]
        par['scienceframe']['exprng'] = [61, None]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])  # Need to convert to : separated
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True)

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        #
        self.meta['decker'] = dict(ext=0, card='SLIT')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERSR', required_ftypes=['science', 'standard'])
        # 
        self.meta['dispangle'] = dict(ext=0, card='BZERO', rtol=2.0)#, required_ftypes=['science', 'standard'])
        self.meta['idname'] = dict(ext=0, card='DATA-TYP')
        self.meta['detector'] = dict(ext=0, card='DET-ID')
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
            binspatial = headarr[0]['BIN-FCT1'] # X
            binspec = headarr[0]['BIN-FCT2'] # Y
            return parse.binning2string(binspec, binspatial)
        else:
            msgs.error("Not ready for this compound meta")

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
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return good_exp & (fitstbl['idname'] == 'DOMEFLAT')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'COMPARISON')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        Args:
            det (:obj:`int`):
                1-indexed detector number.  ESO writes each of the two detectors
                to separate files.  When ``hdu`` is provided, this is ignored
                and instead the chip is determined by the header parameter
                "EXTNAME".  If ``hdu`` is None (for automatically generated
                documentation only), this can be used to set the chip (1 or 2)
                that is returned.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent parameters are set to a default.

        Returns:
            :class:`~pypeit.images.detector_container.DetectorContainer`:
            Object with the detector metadata.
        """
        if hdu is None:
            binning = '1,1'
            chip = '1' if det == 1 else '2'
        else:
            # Binning
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            chip = self.get_meta_value(self.get_headarr(hdu), 'detector')

        '''
        # From Kentaro Aoki on Slack on 2024-06-11
        Each chip has four amplifiers.

        The gains of chip 1 are
        gain1=2.081
        gain2=2.047
        gain3=2.111
        gain4=2.087.
        Those of chip 2 are
        gain1=2.105
        gain2=1.968
        gain3=1.999
        gain4=1.918.

        Readout noise, chip 1:
        noise1=4.2
        noise2=3.8
        noise3=3.6
        noise4=4.0
        Readout noise, chip 2:
        noise1=4.3
        noise2=3.7
        noise3=3.4
        noise4=3.6
        '''
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            # FOCAS puts all data in the primary hdu
            dataext         = 0,
            # FOCAS spectral axis is Y (0 in python/PypeIt indexing)
            specaxis        = 0,
            # FOCAS spectral axis increases in the opposite direction
            # PypeIt expects
            specflip        = True,
            spatflip        = False,
            # arcsec per pixel in the spatial dimension for an unbinned pixel
            platescale      = 0.1038,
            # Value from 2010, probably should be remeasured
            darkcurr        = 0.8,  # e-/pixel/hour
            saturation      = 65535.,
            # need to check on this, 40000 ADU count was provided by Aoki
            nonlinear       = 40000 / 65535.,
            # FIX! not a FOCAS value I could find
            mincounts       = -1e10,
            numamplifiers   = 4,
            # [gain1, gain2, gain3, gain4
            gain            = np.atleast_1d([2.081, 2.047, 2.111, 2.087]),
            # [noise1, noise2, noise3, noise4
            ronoise         = np.atleast_1d([4.2, 3.8, 3.6, 4.0]),
            # NOTE: PypeIt binning variable is "<spectral>,<spatial>"
            # but FOCAS only removes overscan in the spatial
        )

        # CHIP2
        detector_dict2 = dict(
            binning         = binning,
            det             = 1,  # ESO writes these to separate FITS images!!
            # FOCAS puts all data in the primary hdu
            dataext         = 0,
            # FOCAS spectral axis is Y (0 in python/PypeIt indexing)
            specaxis        = 0,
            # FOCAS spectral axis increases in the opposite direction
            # PypeIt expects
            specflip        = True,
            spatflip        = False,
            # arcsec per pixel in the spatial dimension for an unbinned pixel
            platescale      = 0.1038,
            # Value from 2010, probably should be remeasured
            darkcurr        = 0.7,  # e-/pixel/hour
            saturation      = 65535.,
            # need to check on this, 40000 ADU count was provided by Aoki
            nonlinear       = 40000 / 65535.,
            # FIX! not a FOCAS value I could find
            mincounts       = -1e10,
            numamplifiers   = 4,
            # [gain1, gain2, gain3, gain4
            gain            = np.atleast_1d([2.105, 1.968, 1.999, 1.918]),
            # [noise1, noise2, noise3, noise4
            ronoise         = np.atleast_1d([4.3, 3.7, 3.4, 3.6]),
            # NOTE: PypeIt binning variable is "<spectral>,<spatial>"
            # but FOCAS only removes overscan in the spatial
        )
        # Finish
        if chip == '1':
            return detector_container.DetectorContainer(**detector_dict1)
        elif chip == '2':
            return detector_container.DetectorContainer(**detector_dict2)
        else:
            msgs.error(f'Unknown chip: {chip}!')

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

        grism_ID = self.get_meta_value(scifile, 'dispname')

        if grism_ID == 'SCFCGRMB01':
            # 300B grism
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRMB01.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        elif grism_ID == 'SCFCGRMR01':
            # 300R grism
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRMR01.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGREL01':
        #     # 75/mm grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGREL01.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHD45':
        #     # VPH450 grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD45.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        elif grism_ID == 'SCFCGRHD52':
            # VPH520 grism
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD52.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHD65':
        #     # VPH650 grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD65.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHD68':
        #     # VPH680 grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD68.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHD80':
        #      # VPH800 grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD80.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        elif grism_ID == 'SCFCGRHD85':
             # VPH850 grism
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD85.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHD95':
        #     # VPH950 grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHD95.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRLD01':
        #     # 150/mm grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRLD01.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # elif grism_ID == 'SCFCGRHDEC':
        #     # Echelle grism
        #     par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRHDEC.fits'
        #     par['calibrations']['wavelengths']['method'] = 'full_template'
        #     par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'

        # ---- NOTE: from Debora ----
        # The pypeit_sensfunc script uses the config_specific_par() method
        # with a reduced spec1d file to pull out some info, although the method
        # is meant for raw frames. In this case, the spec1d file does not have
        # the dispname meta value in the header and PypeIt tries to run those
        # 2 lines of code (just a message to the terminal) and crashes.
        # So, removing them should fix the problem.
        # ----------------------------
        # else:
        #     msgs.error(f'Not ready for this grism {self.get_meta_value(scifile, "dispname")}')

        return par

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        This method returns a dictionary where the keys of the dictionary are
        the list of configuration-independent frame types. The value of each
        dictionary element can be set to one or more metadata keys that can
        be used to assign each frame type to a given configuration group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations` and how it
        interprets the dictionary values, which can be None.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration-independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': 'detector', 'dark': 'detector'}

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
        # TODO -- Consider dispangle
        # NOTE: FOCAS sometimes substitutes a different slit width for say,
        #       science than for standard, so we don't want to distingish
        #       a unique configuration for 'decker'
        #return ['dispname', 'detector']
        return ['dispname', 'decker', 'detector']

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
        # Definitions for the over scan regions in the DS9 image coordinate.
        # Format
        # 1: start of left overscan region
        # 2: end of left overscan region
        # 3: start of image region
        # 4: end of image region
        # 5: start of right overscan region
        # 6: end of right overscan region
        #
        overscan = {}
        # binning == 1
        overscan[1] = np.asarray([
            # For right (DET_ID == 1) image
            [2, 8, 9, 520, 521, 536],
            [537, 552, 553, 1064, 1065, 1071],
            [1074, 1080, 1081, 1592, 1593, 1608],
            [1610, 1625, 1626, 2137, 2138, 2143],
            #
            # For left (DET_ID == 2) image
            [2, 8, 9, 520, 521, 536],
            [537, 552, 553, 1064, 1065, 1071],
            [1074, 1080, 1081, 1592, 1593, 1608],
            [1609, 1624, 1625, 2136, 2137, 2142],
            ])

        # binning == 2
        overscan[2] = np.asarray([
            # For right (DET-ID == 1) image
            [1, 4, 5, 260, 261, 276],
            [277, 292, 293, 548, 549, 551],
            [553, 556, 557, 812, 813, 828],
            [829, 845, 846, 1101, 1102, 1104],
            #
            # For left (DET-ID == 2) image
            [2, 4, 5, 260, 261, 276],
            [277, 292, 293, 548, 549, 551],
            [553, 556, 557, 812, 813, 828],
            [829, 844, 845, 1100, 1101, 1104],
            ])

        # binning == 4
        overscan[4] = np.asarray([
            # For right (DET-ID == 1) image
            [1, 2, 3, 130, 131, 146],
            [147, 162, 163, 290, 291, 292],
            [293, 294, 295, 422, 423, 438],
            [439, 455, 456, 583, 584, 584],
            #
            # For left (DET-ID == 2) image
            [1, 2, 3, 130, 131, 146],
            [147, 162, 163, 290, 291, 292],
            [293, 294, 295, 422, 423, 438],
            [439, 454, 455, 582, 583, 584],
            ])

        # Read image
        msgs.info(f'Attempting to read FOCAS file: {raw_file}, det={det}')

        # NOTE: io.fits_open checks that the file exists
        hdu_l = io.fits_open(raw_file)
        head = hdu_l[0].header
        data = hdu_l[0].data

        # Grab the detector or mosaic parameters
        # TODO
        # nimg = 1
        # mosaic = None
        detector = self.get_detector_par(det, hdu=hdu_l)

        # Need the exposure time
        exptime = float(head['EXPTIME'])
        bin_x = int(head['BIN-FCT1'])

        arr_shape = data.shape
        # allocate output arrays
        rawdata = data.astype(float)
        rawdatasec_img = np.zeros(arr_shape, dtype=int)
        oscansec_img = np.zeros(arr_shape, dtype=int)

        # collect correct data & overscan for binning in spatial (X) axis
        oscan_arr = overscan[bin_x][(det - 1)*4:][:4]
        if len(oscan_arr) != 4:
            msgs.error(f'FOCAS detector {det} has an unexpected number of overscan regions: {len(oscan_arr)}. '
                       f'Expected 4 (2 for each chip). Please check the overscan definitions in the code.')

        # fill in rawdatasec_img and oscansec_img arrays according to
        # specification above
        for i, oscan in enumerate(oscan_arr):
            n_amp = i + 1

            if n_amp in [2, 4]:
                lft_oscan_start, lft_oscan_end = oscan[0] - 1, oscan[1]
                oscansec_img[:, lft_oscan_start:lft_oscan_end] = n_amp

            data_start, data_end = oscan[2] - 1, oscan[3]
            rawdatasec_img[:, data_start:data_end] = n_amp

            if n_amp in [1, 3]:
                rgt_oscan_start, rgt_oscan_end = oscan[4] - 1, oscan[5]
                oscansec_img[:, rgt_oscan_start:rgt_oscan_end] = n_amp

        return (detector, rawdata, hdu_l, exptime, rawdatasec_img, oscansec_img)

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
        return ['DISPERSR', 'FILTER02', 'SLIT', 'SLT-LEN', 'SLT-PA',
                'SLT-WID']  #, 'SLTCPIX1', 'SLTCPIX2', 'SLTC-RA', 'SLTC-DEC']
