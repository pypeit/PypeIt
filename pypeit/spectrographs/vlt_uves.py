"""
Module for VLT/UVES

.. include:: ../include/links.rst
"""
from pathlib import Path

from IPython import embed

import numpy as np

from astropy.io import fits
from astropy.table import Table

from pypeit import log
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from pypeit.par import parset
from pypeit.images.mosaic import Mosaic
from pypeit.core.mosaic import build_image_mosaic_transform


class UVESMosaicLookUp:
    """
    Provides the geometry required to mosaic VLT/UVES data.
    Similar to :class:`~pypeit.spectrographs.gemini_gmos.GeminiGMOSMosaicLookUp`

    """
    # This is onle the red mosaic. The blue arm just has a single detector.
    geometry = {
        # red -- Note: Dekker et al. (2000) say that the red arm has a mosaic of two 2k x 4k CCDs, with a gap of 0.96mm,
        # corresponding to 64 pixels (15 um pixels), however, RJC trimmed each detector by 3 pixels on each side in the
        # cross-dispersion direction. Therefore, instead of 2048 it's 2042, and instead of 64 pixel gap, it's 70 pixels
        # (3 pixels on each side of the two detectors, plus the 64 pixel gap).
        # Using the fit_mosaic_parameters (on the 564, 580, 760, 860 setups) in the dev-suite, the best fit for the gap is 104 pixels
        'MSC01': {'default_shape': (2042 * 2 + 104.0, 4096),
                  'det1': {'shift': (0., 0.), 'rotation': 0.},
                  'det2': {'shift': (2042.0 + 104.0, 0.0), 'rotation': 0.0}},
    }

class VLTUVESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/UVES specific code.

    This spectrograph is not yet supported.
    """

    telescope = telescopes.VLTTelescopePar()
    url = 'https://www.eso.org/sci/facilities/paranal/instruments/uves.html'
    header_name = 'UVES'
    pypeline = 'Echelle'
    ech_fixed_format = False
    supported = True

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA',
            required_ftypes=['science', 'standard'])  # Need to convert to : separated
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['dispname'] = dict(card=None, compound=True)
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR TYPE')
        self.meta['arm'] = dict(card=None, compound=True)
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        self.meta['echangle'] = dict(card=None, default=0.0)  # There is no header card for this, but it is required
        self.meta['xdangle'] = dict(card=None, compound=True, rtol=0.01)  # There is no tolerance, really, because it's the central wavelength.

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        #par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        #par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames on UVES
        par['calibrations']['pixelflatframe']['exprng'] = [None, 120]
        par['calibrations']['traceframe']['exprng'] = [None, 120]
        par['calibrations']['illumflatframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [1, 600]
        par['scienceframe']['exprng'] = [30, None]

        # Slit tracing
        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3
        par['calibrations']['slitedges']['max_nudge'] = 0.
        par['calibrations']['slitedges']['overlap'] = True
        par['calibrations']['slitedges']['dlength_range'] = 0.25
        par['calibrations']['slitedges']['mask_off_detector'] = False

        par['calibrations']['slitedges']['add_missed_orders'] = True
        par['calibrations']['slitedges']['order_gap_poly'] = 3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 0.1
        par['calibrations']['wavelengths']['sigdetect'] = 4.
        par['calibrations']['wavelengths']['n_first'] = 3
        par['calibrations']['wavelengths']['n_final'] = 4

        par['calibrations']['wavelengths']['match_toler'] = 1.5
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'echelle'
        par['calibrations']['wavelengths']['cc_shift_range'] = (-80.,80.)
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.25
        par['calibrations']['wavelengths']['reid_cont_sub'] = False

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 6
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 2.0
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
        # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['reduce']['findobj']['find_trim_edge'] = [3, 3]
        # number of objects
        par['reduce']['findobj']['maxnumber_sci'] = 2  # Assume that there is max two object in each order.
        par['reduce']['findobj']['maxnumber_std'] = 1  # Assume that there is only one object in each order.

        # Coadding
        par['coadd1d']['wave_method'] = 'log10'

        return par

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
            try:
                binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            except KeyError:
                log.warning("Cannot determine spatial binning from the header. Setting to 1")
                binspatial = 1
            try:
                binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            except KeyError:
                log.warning("Cannot determine spectral binning from the header. Setting to 1")
                binspec = 1
            # Parse the binning information into a string
            return parse.binning2string(binspec, binspatial)
        elif meta_key == 'arm' or meta_key == 'dispname':
            if 'HIERARCH ESO TPL NAME' in headarr[0]:
                tplid = headarr[0]['HIERARCH ESO TPL NAME'].lower()
                if 'blue' in tplid:
                    arm = 'BLUE'
                elif 'red' in tplid:
                    arm = ('RED')
                elif 'HIERARCH ESO INS PATH' in headarr[0]:
                    arm = headarr[0]['HIERARCH ESO INS PATH']
                else:
                    arm = 'None'
            return arm
        elif meta_key == 'xdangle':
            if 'HIERARCH ESO INS GRAT1 WLEN' in headarr[0]:
                cwlen = headarr[0]['HIERARCH ESO INS GRAT1 WLEN']
            elif 'HIERARCH ESO INS GRAT2 WLEN' in headarr[0]:
                cwlen = headarr[0]['HIERARCH ESO INS GRAT2 WLEN']
            else:
                cwlen = 'None'
            return cwlen
        else:
            log.error("Not ready for this compound meta")

    def configuration_keys(self):
        """
        Return the metadata keys that define a unique instrument
        configuration.

        This list is used by :class:`~pypeit.metadata.PypeItMetaData` to
        identify the unique configurations among the list of frames read
        for a given reduction.

        Returns:
            :obj:`list`: List of keywords of data pulled from file headers
            and used to construct the :class:`~pypeit.metadata.PypeItMetaData`
            object.
        """
        return ['xdangle', 'arm', 'binning']

    def config_independent_frames(self):
        """
        Define frame types that are independent of the fully defined
        instrument configuration.

        Bias and dark frames are considered independent of a configuration,
        but the DATE-OBS keyword is used to assign each to the most-relevant
        configuration frame group. See
        :func:`~pypeit.metadata.PypeItMetaData.set_configurations`.

        Returns:
            :obj:`dict`: Dictionary where the keys are the frame types that
            are configuration independent and the values are the metadata
            keywords that can be used to assign the frames to a configuration
            group.
        """
        return {'bias': ['binning', 'arm'], 'dark': ['binning', 'arm']}

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
        return ['HIERARCH ESO SEQ ARM']

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
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'OBJECT')
                | (fitstbl['idname'] == 'OBJECT,POINT')
                | (fitstbl['idname'] == 'SCIENCE')
                | (fitstbl['idname'] == 'STD,TELLURIC')
                | (fitstbl['idname'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'DARK')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'LAMP,FLAT')
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'LAMP,WAVE')

        log.warning('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

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


class VLTUVESBlueSpectrograph(VLTUVESSpectrograph):
    
    name = 'vlt_uves_blue'
    camera = 'VLT_UVES_blue'
    ndet = 1

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS SLIT2 WID')

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Adjustments to parameters for Keck HIRES
        turn_off_on = dict(use_biasimage=False, use_overscan=True, overscan_method='median')
        par.reset_all_processimages_par(**turn_off_on)
        # Right now we are using the overscan and not biases becuase the
        # standards are read with a different read mode and we don't yet have
        # the option to use different sets of biases for different standards,
        # or use the overscan for standards but not for science frames

        # Slit tracing
        par['calibrations']['slitedges']['order_width_poly'] = 2

        # Extraction
        par['reduce']['skysub']['sky_sigrej'] = 4.0

        return par
        
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
        binfact = int(binning.split(',')[0]) * int(binning.split(',')[1])
        ronoise = 3.8/np.sqrt(binfact) if hdu is None else hdu[0].header['HIERARCH ESO DET OUT1 RON']
        gain = 0.54 if hdu is None else hdu[0].header['HIERARCH ESO DET OUT1 GAIN']

        # Detector
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.22,
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.7, # Website says 0.6, but we'll push it a bit
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d([gain]),
            ronoise         = np.atleast_1d([ronoise]),
            datasec         = np.atleast_1d('[:,51:2098]'), # Any changes to this line should also be made in the mosaic geometry at the top of this file.
            oscansec        = np.atleast_1d('[:,4:50]'),
            )

        return detector_container.DetectorContainer(**detector_dict)

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
        par = super().config_specific_par(inp, inp_par=inp_par)

        bin_spec, bin_spat = parse.parse_binning(self.get_meta_value(inp, 'binning'))

        # slit edges
        # NOTE: With add_missed_orders set to True and order_spat_range set to the
        # default (None), the code will try to add missing orders over the full
        # range of the detector mosaic!
        par['calibrations']['slitedges']['order_spat_range'] = [-50., (2048.0+50.0)/bin_spat]

        # wavelength
        par['calibrations']['wavelengths']['fwhm'] = 8.0/bin_spec

        # Return
        return par

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        angle_fits_file = 'vlt_uves_blue_angle_fits.fits'
        composite_arc_file = 'vlt_uves_blue_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]


class VLTUVESRedSpectrograph(VLTUVESSpectrograph):

    name = 'vlt_uves_red'
    camera = 'VLT_UVES_red'
    ndet = 2

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        super().init_meta()
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS SLIT3 WID')

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.

        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['rdx']['detnum'] = [(1,2)]

        # Set default processing for slitless_pixflat
        par['calibrations']['slitless_pixflatframe']['process']['scale_to_mean'] = True

        # Slit tracing
        par['calibrations']['slitedges']['mask_off_detector'] = True
        par['calibrations']['slitedges']['order_width_poly'] = 4

        # Wavelength calibration should be done separately for a mosaic
        par['calibrations']['wavelengths']['ech_separate_2d'] = True  # Doesn't seem like there's an offset+rotation that works for VLT/UVES red (note that blue doesn't have a mosaic)

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5
        par['sensfunc']['IR']['telgridfile'] = 'TellPCA_3000_10500_R120000.fits'
        par['sensfunc']['IR']['pix_shift_bounds'] = (-40.0,40.0)

        # Telluric parameters
        # Allow for a large helio shift with UVES
        par['telluric']['pix_shift_bounds'] = (-40.0,40.0)
        # Similarly, the resolution guess is higher than it should be
        par['telluric']['resln_frac_bounds'] = (0.25,1.25)

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
        par = super().config_specific_par(inp, inp_par=inp_par)

        bin_spec, bin_spat = parse.parse_binning(self.get_meta_value(inp, 'binning'))

        # slit edges
        # NOTE: With add_missed_orders set to True and order_spat_range set to the
        # default (None), the code will try to add missing orders over the full
        # range of the detector mosaic!
        # par['calibrations']['slitedges']['order_spat_range'] = [-50/bin_spat, (2042+50)/bin_spat]
        offset = 50.0  # Extra number of pixels to add to the end of the mosaic to allow for missed orders.
        par['calibrations']['slitedges']['order_spat_range'] = [-offset/bin_spat, (2042 * 2 + 104.0 + offset)/bin_spat]

        # wavelength
        par['calibrations']['wavelengths']['fwhm'] = 8.0/bin_spec

        # Return
        return par

    @property
    def allowed_mosaics(self):
        """
        Return the list of allowed detector mosaics.

        Only red arm on VLT/UVES requires mosaicing.

        Returns:
            :obj:`list`: List of tuples, where each tuple provides the 1-indexed
            detector numbers that can be combined into a mosaic and processed by
            PypeIt.
        """
        return [(1,2)]

    @property
    def default_mosaic(self):
        return self.allowed_mosaics[0]

    def get_mosaic_par(self, mosaic, hdu=None, msc_ord=0):
        """
        Return the hard-coded parameters needed to construct detector mosaics
        from unbinned images.

        The parameters expect the images to be trimmed and oriented to follow
        the PypeIt shape convention of ``(nspec,nspat)``.  For returned
        lists, the length of the list is the same as the number of detectors in
        the mosaic, and they are ordered by the detector number.

        Args:
            mosaic (:obj:`tuple`):
                Tuple of detector numbers used to construct the mosaic.  Must be
                one among the list of possible mosaics as hard-coded by the
                :func:`allowed_mosaics` function.
            hdu (`astropy.io.fits.HDUList`_, optional):
                The open fits file with the raw image of interest.  If not
                provided, frame-dependent detector parameters are set to a
                default.  BEWARE: If ``hdu`` is not provided, the binning is
                assumed to be `1,1`, which will cause faults if applied to
                binned images!
            msc_ord (:obj:`int`, optional):
                Order of the interpolation used to construct the mosaic.

        Returns:
            :class:`~pypeit.images.mosaic.Mosaic`: Object with the mosaic *and*
            detector parameters.
        """

        # Validate the entered (list of) detector(s)
        nimg, _ = self.validate_det(mosaic)

        # Index of mosaic in list of allowed detector combinations
        mosaic_id = self.allowed_mosaics.index(mosaic)+1
        detid = f'MSC0{mosaic_id}'

        # Get the detectors
        detectors = np.array([self.get_detector_par(det, hdu=hdu) for det in mosaic])
        # Binning *must* be consistent for all detectors
        if any(d.binning != detectors[0].binning for d in detectors[1:]):
            log.error('Binning is somehow inconsistent between detectors in the mosaic!')

        # Collect the offsets and rotations for *all unbinned* detectors in the
        # full instrument, ordered by the number of the detector.  Detector
        # numbers must be sequential and 1-indexed.
        # See the mosaic documentation.
        msc_geometry = UVESMosaicLookUp.geometry
        expected_shape = msc_geometry[detid]['default_shape']
        shift = np.array([(msc_geometry[detid]['det1']['shift'][0], msc_geometry[detid]['det1']['shift'][1]),
                          (msc_geometry[detid]['det2']['shift'][0], msc_geometry[detid]['det2']['shift'][1])])

        rotation = np.array([msc_geometry[detid]['det1']['rotation'], msc_geometry[detid]['det2']['rotation']])

        # The binning and process image shape must be the same for all images in
        # the mosaic
        binning = tuple(int(b) for b in detectors[0].binning.split(','))
        shape = tuple(n // b for n, b in zip(expected_shape, binning))

        msc_sft = [None]*nimg
        msc_rot = [None]*nimg
        msc_tfm = [None]*nimg

        for ii in range(nimg):
            msc_sft[ii] = shift[ii]
            msc_rot[ii] = rotation[ii]
            # binning is here in the PypeIt convention of (binspec, binspat), but the mosaic transformations
            # occur in the raw data frame, which flips spectral and spatial
            msc_tfm[ii] = build_image_mosaic_transform(shape, msc_sft[ii], msc_rot[ii], tuple(reversed(binning)))

        return Mosaic(mosaic_id, detectors, shape, np.array(msc_sft), np.array(msc_rot),
                      np.array(msc_tfm), msc_ord)

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

        # Detector base parameters
        detector_base = dict(
            binning         = binning,
            specaxis        = 0,
            specflip        = True,
            spatflip        = True,
            platescale      = 0.135,
            darkcurr        = 0.0,  # e-/pixel/hour
            saturation      = 65535.,
            nonlinear       = 0.7, # Website says 0.6, but we'll push it a bit
            mincounts       = -1e10,
            numamplifiers   = 1,
            # Placeholders, will be updated for each detector -- must all be None, so that  the code will break
            # if it tries to obtain information without passing in the hdu argument.
            det=None,
            dataext=None,
            gain=None,
            ronoise=None,
            datasec=None,
            oscansec=None,
        )
        # Now, depending on the HDU format (which changed at some point into multi-extension fits files),
        # we need to extract information from different HDUs for the two detectors.
        detector_dict1 = detector_base.copy()
        detector_dict2 = detector_base.copy()
        if hdu is not None:
            if len(hdu) == 1:
                # This is the old format, where the two detectors are stored in a single HDU
                # Detector 1
                detector_dict1.update(dict(
                    det=1,
                    dataext=0,
                    gain=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT4 GAIN']]),
                    ronoise=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT4 RON']]),
                    datasec=np.atleast_1d('[:,2201:4242]'), # Any changes to this line requires changes to the mosaic at the top of this file.
                    oscansec=np.atleast_1d('[:,4246:4292]'),
                ))

                # Detector 2
                detector_dict2.update(dict(
                    det=2,
                    dataext=0,
                    gain=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[0].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'),
                    oscansec=np.atleast_1d('[:,2102:2148]'),
                ))
            else:
                # This is the new format, where the two detectors are stored in separate HDUs.
                # Detector 1.
                detector_dict1.update(dict(
                    det=1,
                    dataext=2,
                    gain=np.atleast_1d([hdu[2].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[2].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'), # Any changes to this line requires changes to the mosaic at the top of this file.
                    oscansec=np.atleast_1d('[:,2102:]'),
                ))

                # Detector 2.
                detector_dict2.update(dict(
                    det=2,
                    dataext=1,
                    gain=np.atleast_1d([hdu[1].header['HIERARCH ESO DET OUT1 GAIN']]),
                    ronoise=np.atleast_1d([hdu[1].header['HIERARCH ESO DET OUT1 RON']]),
                    datasec=np.atleast_1d('[:,57:2098]'),
                    oscansec=np.atleast_1d('[:,2102:]'),
                ))

        # Instantiate
        detector_dicts = [detector_dict1, detector_dict2]
        return detector_container.DetectorContainer( **detector_dicts[det-1])

    def get_echelle_angle_files(self):
        """ Pass back the files required
        to run the echelle method of wavecalib

        Returns:
            list: List of files
        """
        angle_fits_file = 'vlt_uves_red_angle_fits.fits'
        composite_arc_file = 'vlt_uves_red_composite_arc.fits'

        return [angle_fits_file, composite_arc_file]
