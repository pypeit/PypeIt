"""
Module for Shane/Kast specific methods.

.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from astropy.time import Time
from astropy.coordinates import SkyCoord
from astropy.coordinates import EarthLocation
from astropy.coordinates import AltAz

from pypeit import log
from pypeit import PypeItError
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core import meta
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.core import parse


class ShaneHamspecSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Hamspec specific code
    """
    ndet = 1
    telescope = telescopes.ShaneTelescopePar()
    url = 'https://mthamilton.ucolick.org/techdocs/instruments/hamspec/index.html/'
    ql_supported = False

    name = 'shane_hamspec'
    camera = 'Hamspec'
    supported = True
    header_name = 'hamspec'
    pypeline = 'Echelle'
    ech_fixed_format = False

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Copied from Keck/HIRES

        # No bias or dark frames are taken with Hamspec, so subtract the
        # overscan (the raw frames have a 64-column overscan, two amps) and
        # disable bias subtraction for every processed image type.
        turn_off_on = dict(use_biasimage=False, use_overscan=True,
                           overscan_method='median')
        par.reset_all_processimages_par(**turn_off_on)

        # Slit tracing
        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 5  # 8 was too high
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 0.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.25

        # The real edge tracing already finds ~the true number of orders
        # (~100).  Leaving add_missed_orders on inserts ~140 spurious orders in
        # the dark detector regions (Hamspec orders do not fill the detector),
        # which corrupts the echelle order-number identification.  Disable it.
        par['calibrations']['slitedges']['add_missed_orders'] = False
        par['calibrations']['slitedges']['order_width_poly'] = 2
        par['calibrations']['slitedges']['order_gap_poly'] = 3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5  # [5, 5, 5] + 12*[7] # + [5]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.1
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 4

        par['calibrations']['wavelengths']['match_toler'] = 1.5
        # Reidentification parameters.  We follow the Keck/NIRES recipe: a
        # single composite-arc template (built from the legacy ThAr solution,
        # one ThAr spectrum + wavelength solution per order) is reidentified
        # against each observed order by cross-correlation.  This avoids the
        # HIRES-style 'echelle' angle_fits model, which would need ThAr at
        # many cross-disperser settings we do not have (see Report 04).
        par['calibrations']['wavelengths']['method'] = 'echelle'
        par['calibrations']['wavelengths']['cc_shift_range'] = (-200.,200.)
        # NOTE (Q36a): lowering cc_thresh to 0.5 and enabling reid_cont_sub were
        # tried to pick up the bluer orders but did not help (no extra orders,
        # higher RMS), so keep the original values.  The blue-order coverage
        # limit is in the observed-data matching, not these thresholds.
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.25
        par['calibrations']['wavelengths']['reid_cont_sub'] = False

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 3
        par['calibrations']['wavelengths']['ech_sigrej'] = 2.0
        # Hamspec is a single detector (unlike the HIRES mosaic this was copied
        # from), so do not separate the 2D wavelength fit by detector — there is
        # no det_img and ech_separate_2d=True crashes echelle_2dfit.
        par['calibrations']['wavelengths']['ech_separate_2d'] = False
        par['calibrations']['wavelengths']['bad_orders_maxfrac'] = 0.5

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['global_sky_std'] = False
        # local sky subtraction operates on entire slit
        par['reduce']['extraction']['model_full_slit'] = True
        # The echelle orders are very short in the spatial direction (~7 px on
        # the Loral CCD), so trim only 1 edge pixel each side — the default
        # (5,5), or even (3,3), leaves too few pixels for object finding and no
        # objects are detected on the standard/science.
        par['reduce']['findobj']['find_trim_edge'] = [1, 1]
        # Continnum order for determining thresholds

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5 #[9, 11, 11, 9, 9, 8, 8, 7, 7, 7, 7, 7, 7, 7, 7]
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_10500_R120000.fits'
        par['sensfunc']['IR']['pix_shift_bounds'] = (-40.0,40.0)
        
        # Telluric parameters
        # HIRES is usually oversampled, so the helio shift can be large
        par['telluric']['pix_shift_bounds'] = (-40.0,40.0)
        # Similarly, the resolution guess is higher than it should be
        par['telluric']['resln_frac_bounds'] = (0.25,1.25)

        # Coadding
        par['coadd1d']['wave_method'] = 'log10'        


        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, 61]
        par['calibrations']['standardframe']['exprng'] = [1, 301]
        #
        par['scienceframe']['exprng'] = [301, None]
        #par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_26000_R10000.fits'
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
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(ext=0, card='PLATENAM')  # Aperture plate
        # Follow the XIDL hamspec convention (Q&A 43): the echelle angle is the
        # dewar height (DHEITRAW) and the cross-disperser angle is the grating
        # tilt (GTILTRAW).  Both are raw stepper counts present in every era's
        # headers (DHEITVAL is absent in the older Loral data).
        self.meta['echangle'] = dict(ext=0, card='DHEITRAW', rtol=1e-2)
        self.meta['xdangle'] = dict(ext=0, card='GTILTRAW', rtol=1e-2)

        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['dispname'] = dict(ext=0, card=None, default='Hamspec')

        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card=None, compound=True)
        # Additional ones, generally for configuration determination or time
        self.meta['filter1'] = dict(ext=0, card='DFILTNAM')
        self.meta['lampstat01'] = dict(ext=0, card='LAMPPOS')

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
        if meta_key == 'mjd':
            return Time(self._obs_isot(headarr[0]), format='isot').mjd
        elif meta_key == 'airmass':
            # The raw headers have no AIRMASS card, so compute it from the
            # pointing (RA/DEC), the readout time, and the Shane/Lick location
            # via the reusable pypeit.core.meta.airmass helper.
            ra, dec = meta.convert_radec(headarr[0]['RA'], headarr[0]['DEC'])
            obstime = Time(self._obs_isot(headarr[0]), format='isot')
            return meta.airmass(ra, dec, obstime,
                                self.telescope['longitude'],
                                self.telescope['latitude'],
                                self.telescope['elevation'])
        raise PypeItError("Not ready for this compound meta")

    @staticmethod
    def _obs_isot(hdr):
        """
        Return the UT observation time as an ISOT string, handling both
        header eras: the 4kx4k data use the ``DATE`` card, the older Loral
        data use ``DATE-OBS``.

        Generated by JXP and Claude.

        Args:
            hdr (`astropy.io.fits.Header`_): Primary header of the raw frame.

        Returns:
            :obj:`str`: ISOT UT date/time string.
        """
        if 'DATE' in hdr and hdr['DATE'] not in (None, ''):
            return hdr['DATE']
        return hdr['DATE-OBS']

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
        # decker is not included because arcs are often taken with a 0.5" slit.
        # The configuration is set by the grating tilt (echangle <- GTILTRAW,
        # the along-dispersion wavelength) and the cross-disperser / dewar
        # height (xdangle <- DHEITRAW, which orders land where); both must
        # match for frames to share a calibration (see Report 04, Q&A 21).
        #return ['decker', 'filter1', 'xdangle', 'binning']
        return ['echangle', 'xdangle', 'binning']

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
            return good_exp & self.lamps(fitstbl, 'Off')
        if ftype == 'bias':
            return good_exp # & (fitstbl['target'] == 'Bias')
        # The aperture plate (PLATENAM) is stored under the PypeIt 'decker'
        # meta key.  The wide plate floods the detector (pixelflat); the narrow
        # science plate gives the order edges (trace) and slit illumination
        # (illumflat).  The plate names differ by era: 4kx4k (2014) uses
        # 800:6.0 wide / 640:5.0 narrow; the Loral (2010) uses 800:5.0 wide /
        # 800:2.5 narrow.
        wide_plates = ['800:6.0', '800:5.0']           # -> pixelflat
        narrow_plates = ['640:5.0', '800:2.5']         # -> trace/illumflat
        if ftype in ['pixelflat']:
            return good_exp & self.lamps(fitstbl, 'flat') \
                & np.isin(fitstbl['decker'], wide_plates)
        if ftype in ['trace', 'illumflat']:
            return good_exp & self.lamps(fitstbl, 'flat') \
                & np.isin(fitstbl['decker'], narrow_plates)
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & self.lamps(fitstbl, 'arcs')

        log.warning('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (`astropy.table.Table`_):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check. Can be ``'off'``, ``'arcs'``, or
                ``'dome'``.

        Returns:
            `numpy.ndarray`_: A boolean array selecting fits files that meet
            the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status.lower() == 'off':
            # Check if all are off.  The LAMPPOS card reads 'Off' (capital
            # O), so compare case-insensitively against 'off'/'none'.
            return np.all(np.array(
                [np.isin(np.char.lower(np.asarray(fitstbl[k], dtype=str)),
                         ['off', 'none'])
                 for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        elif status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,2) ]
            return np.any(np.array([ fitstbl[k] == 'Thorium-Argon' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        elif status == 'flat':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,2) ]
            return np.any(np.array([ fitstbl[k] == 'PolarQuartz' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        else:
            raise ValueError('No implementation for status = {0}'.format(status))



    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.
        https://mthamilton.ucolick.org/techdocs/instruments/hamspec/CCDs/

        Hamspec has used two CCDs over the years, which differ in size, number
        of amplifiers, gain, and read noise.  We pick the right one from the
        raw-frame column count (the e2v has 4160 columns, the Loral 1485x2080).
        The gain/read-noise/amp values come from the XIDL hamspec pipeline
        (``hamspec_strct.pro``).

        Generated by JXP and Claude.

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
        binning = '1,1' if hdu is None \
            else self.get_meta_value(self.get_headarr(hdu), 'binning')
        # Decide which CCD this is from the raw column count (e2v=4160).
        # Default (no hdu) is the e2v 4kx4k.
        ncol = None if hdu is None else hdu[0].header.get('NAXIS1', None)
        is_loral = ncol is not None and ncol < 3000

        if is_loral:
            # Loral 2Kx2K (the older detector, used e.g. for the 2010 Cooke
            # data): single amplifier, 2048 data + 32 overscan columns.
            # gain/readno from XIDL hamspec_strct.pro (Loral2Kx2K).
            detector_dict = dict(
                binning=binning, det=1, dataext=0,
                specaxis=1, specflip=False, spatflip=True,
                platescale=0.5,   # arcsec; approximate, verify if needed
                saturation=65535., mincounts=-1e10, nonlinear=0.86,
                numamplifiers=1,
                gain=np.asarray([1.2]),
                ronoise=np.asarray([18.0]),
                xgap=0., ygap=0., ysize=1., darkcurr=0.0,
                # rows, columns on the raw frame, 1-indexed
                datasec=np.asarray(['[:, 1:2048]']),
                oscansec=np.asarray(['[:, 2049:2080]']),
            )
            return detector_container.DetectorContainer(**detector_dict)

        # e2v CCD203-82 4kx4k thin: two amplifiers, split in columns.
        detector_dict = dict(
            binning=binning,
            det=1,
            dataext=0,
            specaxis=1,
            specflip=False,
            spatflip=True,
            platescale=0.5, # arcsec
            saturation=65535.,
            mincounts=-1e10,
            nonlinear=0.86,
            numamplifiers=2,
            gain=np.asarray([0.92,0.92]),
            ronoise=np.asarray([3.2,3.2]),
            xgap=0.,
            ygap=0.,
            ysize=1.,
            darkcurr=0.0,  # e-/pixel/hour
            # These are rows, columns on the raw frame, 1-indexed
            datasec=np.asarray(['[:, 1:2048]', '[:, 2049:4096]']),
            oscansec=np.asarray(['[:, 4098:4125]', '[:, 4129:4157]']),
        )
        return detector_container.DetectorContainer(**detector_dict)

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
        par = super().config_specific_par(scifile, inp_par=inp_par)
        
        #headarr = self.get_headarr(scifile)

        # wavelength
        #par['calibrations']['wavelengths']['fwhm'] = 8.0/bin_spec

        # Return
        return par



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
        #return ['GRISM_N', 'BSPLIT_N']
        return ['DFILTNAM', 'GTILTRAW', 'PLATENAM']



    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        # These archive files do not exist yet; they will be built as part of
        # the wavelength-calibration development.  Named to match the
        # spectrograph (shane_hamspec), not the legacy lick_hamspec prefix.
        angle_fits_file = 'shane_hamspec_angle_fits.fits'
        composite_arc_file = 'shane_hamspec_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]
        

    def order_platescale(self, order_vec, binning=None):
        """
        Return the platescale for each echelle order.

        This routine is only defined for echelle spectrographs, and it is
        undefined in the base class.

        Args:
            order_vec (`numpy.ndarray`_):
                The vector providing the order numbers.
            binning (:obj:`str`, optional):
                The string defining the spectral and spatial binning.

        Returns:
            `numpy.ndarray`_: An array with the platescale for each order
            provided by ``order``.
        """
        det = self.get_detector_par(1)
        binspectral, binspatial = parse.parse_binning(binning)

        # Assume no significant variation (which is likely true)
        return np.ones_like(order_vec)*det.platescale*binspatial