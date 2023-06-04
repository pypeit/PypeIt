"""
Module for VLT FORS (1 and 2)

.. include:: ../include/links.rst
"""
import numpy as np
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.core import meta
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container
from astropy.coordinates import SkyCoord
from astropy import units
from astropy.io import fits
from IPython import embed

class VLTFORSSpectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle VLT/FORS specific code
    Parent for FORS1 and FORS2
    """
    ndet = 1  # Because each detector is written to a separate FITS file
    telescope = telescopes.VLTTelescopePar()
    url = 'https://www.eso.org/sci/facilities/paranal/instruments/fors.html'

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
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan_method'] = 'median'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 3
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] = 25.0
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 4

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'ArI']  # Grating dependent
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0  # Good for 2x binning
        par['calibrations']['wavelengths']['n_final'] = 4

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 5
        par['sensfunc']['IR']['telgridfile'] = 'TelFit_Paranal_VIS_9800_25000_R25000.fits'



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

        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        #
        self.meta['decker'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='HIERARCH ESO INS GRIS1 NAME', required_ftypes=['science', 'standard'])
        self.meta['dispangle'] = dict(ext=0, card='HIERARCH ESO INS GRIS1 WLEN', rtol=2.0, required_ftypes=['science', 'standard'])
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
        self.meta['detector'] = dict(ext=0, card='EXTNAME')
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
            binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key == 'decker':
            try:  # Science
                decker = headarr[0]['HIERARCH ESO INS SLIT NAME']
            except KeyError:  # Standard!
                try:
                    decker = headarr[0]['HIERARCH ESO SEQ SPEC TARG']
                except KeyError:
                    return None
            return decker
        else:
            msgs.error("Not ready for this compound meta")

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
        return []

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
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & ((fitstbl['target'] == 'STD,FLUX')
                               | (fitstbl['target'] == 'STD'))
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'FLAT,LAMP')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'LAMP,WAVE')
                               | (fitstbl['target'] == 'WAVE,LAMP'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class VLTFORS2Spectrograph(VLTFORSSpectrograph):
    """
    Child to handle VLT/FORS2 specific code
    """

    name = 'vlt_fors2'
    camera = 'FORS2'
    header_name = 'FORS2'
    supported = True
    comment = '300I, 300V gratings'

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
            chip = 'CHIP1' if det == 1 else 'CHIP2'
        else:
            # Binning
            # TODO: Could this be detector dependent??
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            chip = self.get_meta_value(self.get_headarr(hdu), 'detector')

        # These numbers are from the ESO FORS2 user manual at: 0
        # http://www.eso.org/sci/facilities/paranal/instruments/fors/doc/VLT-MAN-ESO-13100-1543_P01.1.pdf
        # They are for the MIT CCD (which is the default detector) for the high-gain, 100 khZ readout mode used for
        # spectroscpy. The other readout modes are not yet implemented. The E2V detector is not yet supported!!

        # CHIP1
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 2.1,
            saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.70),
            ronoise         = np.atleast_1d(2.9), # High gain
            datasec         = np.atleast_1d('[11:2059,:]'),  # For 1x binning, I think
            #oscansec=np.atleast_1d('[2062:,:]'),
            oscansec=np.atleast_1d('[1:10,:]'), # Overscan has artifacts so use pre-scan
        )
        # CHIP2
        detector_dict2 = dict(
            binning         = binning,
            det             = 1,  # ESO writes these to separate FITS images!!
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 1.4,
            saturation      = 2.0e5,  # I think saturation may never be a problem here since there are many DITs
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.70),
            ronoise         = np.atleast_1d(3.15),  # High gain
            datasec=np.atleast_1d('[11:2059,:]'),
            oscansec=np.atleast_1d('[2062:,:]'), # Pre-scan has artifacts, so use overscan
            #datasec=np.atleast_1d('[20:,0:2048]'),
            #oscansec=np.atleast_1d('[4:20,4:2044]'),
        )
        # Finish
        if chip == 'CHIP1':
            return detector_container.DetectorContainer(**detector_dict1)
        elif chip == 'CHIP2':
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
        # TODO: Should we allow the user to override these?

        #detector = self.get_meta_value(scifile, 'detector')
        #self.set_detector(detector)
        # Wavelengths
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        if self.get_meta_value(scifile, 'dispname') == 'GRIS_300I':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_fors2_300I.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == 'GRIS_300V':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_fors2_300V.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
        elif self.get_meta_value(scifile, 'dispname') == 'GRIS_600z':
            par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
            par['calibrations']['wavelengths']['method'] = 'holy-grail'
            # Since we are using the sky to fit the wavelengths don't correct for flexure
            par['flexure']['spec_method'] = 'skip'
            #par['reduce']['skysub']['bspline_spacing'] = 0.6

        if 'lSlit' in self.get_meta_value(scifile, 'decker') or 'LSS' in self.get_meta_value(scifile, 'decker'):
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'


        return par

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
        return ['dispname', 'dispangle', 'decker', 'detector']

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
        return ['HIERARCH ESO INS GRIS1 NAME', 'HIERARCH ESO INS GRIS1 WLEN',
                'HIERARCH ESO INS SLIT NAME', 'HIERARCH ESO SEQ SPEC TARG']

    # TODO -- Convert this into get_comb_group()
    def parse_dither_pattern(self, file_list, ext=None):
        """
        Parse headers from a file list to determine the dither pattern.

        Parameters
        ----------
        file_list (list of strings):
            List of files for which dither pattern is desired
        ext (int, optional):
            Extension containing the relevant header for these files. Default=None. If None, code uses
            self.primary_hdrext

        Returns
        -------
        dither_pattern, dither_id, offset_arcsec

        dither_pattern (str `numpy.ndarray`_):
            Array of dither pattern names
        dither_id (str `numpy.ndarray`_):
            Array of dither pattern IDs
        offset_arc (float `numpy.ndarray`_):
            Array of dither pattern offsets
        """
        nfiles = len(file_list)
        offset_arcsec = np.zeros(nfiles)
        dither_pattern = None
        dither_id = None
        for ifile, file in enumerate(file_list):
            hdr = fits.getheader(file, self.primary_hdrext if ext is None else ext)
            try:
                ra, dec = meta.convert_radec(self.get_meta_value(hdr, 'ra', no_fussing=True),
                                    self.get_meta_value(hdr, 'dec', no_fussing=True))
            except:
                msgs.warn('Encounter invalid value of your coordinates. Give zeros for both RA and DEC. Check that this does not cause problems with the offsets')
                ra, dec = 0.0, 0.0
            if ifile == 0:
                coord_ref = SkyCoord(ra*units.deg, dec*units.deg)
                offset_arcsec[ifile] = 0.0
                # ESOs position angle appears to be the negative of the canonical astronomical convention
                posang_ref = -(hdr['HIERARCH ESO INS SLIT POSANG']*units.deg)
                posang_ref_rad = posang_ref.to('radian').value
                # Unit vector pointing in direction of slit PA
                u_hat_slit = np.array([np.sin(posang_ref), np.cos(posang_ref)]) # [u_hat_ra, u_hat_dec]
            else:
                coord_this = SkyCoord(ra*units.deg, dec*units.deg)
                posang_this = coord_ref.position_angle(coord_this).to('deg')
                separation  = coord_ref.separation(coord_this).to('arcsec').value
                ra_off, dec_off = coord_ref.spherical_offsets_to(coord_this)
                u_hat_this  = np.array([ra_off.to('arcsec').value/separation, dec_off.to('arcsec').value/separation])
                dot_product = np.dot(u_hat_slit, u_hat_this)
                if not np.isclose(np.abs(dot_product),1.0, atol=1e-2):
                    msgs.error('The slit appears misaligned with the angle between the coordinates: dot_product={:7.5f}'.format(dot_product) + msgs.newline() +
                               'The position angle in the headers {:5.3f} differs from that computed from the coordinates {:5.3f}'.format(posang_this, posang_ref))
                offset_arcsec[ifile] = separation*np.sign(dot_product)

#            dither_id.append(hdr['FRAMEID'])
#            offset_arcsec[ifile] = hdr['YOFFSET']
        return dither_pattern, dither_id, offset_arcsec