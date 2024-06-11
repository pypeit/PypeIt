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
    supported = False
    comment = ' just getting started'


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
        par['calibrations']['arcframe']['exprng'] = [None, 1]
        par['calibrations']['standardframe']['exprng'] = [1, 61]
        #
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
        # TODO - FIX THIS!!
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
            # TODO -- CHECK THE FOLLOWING
            binning = parse.binning2string(binspec, binspatial)
            return binning
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
            # Flats and trace frames are typed together
            # TODO -- Are there internal flats?
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

        # These numbers are from the ESO FORS2 user manual at: 0
        # http://www.eso.org/sci/facilities/paranal/instruments/fors/doc/VLT-MAN-ESO-13100-1543_P01.1.pdf
        # They are for the MIT CCD (which is the default detector) for the high-gain, 100 khZ readout mode used for
        # spectroscopy. The other readout modes are not yet implemented. The E2V detector is not yet supported!!

        # CHIP1
        # TODO - UPDATE ALL OF THIS

        # TODO -- Deal with dataswec, oscansec
        detector_dict1 = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 2.1,  # e-/pixel/hour
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
            specaxis        = 0,
            specflip        = True,
            spatflip        = False,
            platescale      = 0.126,  # average between order 11 & 30, see manual
            darkcurr        = 1.4,  # e-/pixel/hour
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

        if self.get_meta_value(scifile, 'dispname') == 'SCFCGRMB01':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_subaru_focas_SCFCGRMB01.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'
            par['calibrations']['wavelengths']['stretch_func'] = 'quadratic'
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
        #return ['dispname', 'dispangle', 'decker', 'detector']
        # TODO -- Consider dispangle
        return ['dispname', 'decker', 'detector']


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
