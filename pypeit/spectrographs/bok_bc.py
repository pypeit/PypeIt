"""
Module for Bok/B&C specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit import io
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class BokBCSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle BOK specific code
    """
    ndet = 1
    telescope = telescopes.BokTelescopePar()
    name = 'bok_bc'
    camera = 'BC'
    url = 'http://james.as.arizona.edu/~psmith/90inch/90inch.html'
    comment = 'Bok B&C spectrometer'
    header_name = 'Bok B&C spectrometer'
    supported = True

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
        self.meta['decker'] = dict(ext=0, card='APERTURE')
        self.meta['binning'] = dict(card=None, compound=True)
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='DISPERSE')
        self.meta['dispangle'] = dict(ext=0, card='TILTPOS', rtol=1e-3)
        self.meta['idname'] = dict(ext=0, card='OBJECT')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        # used for arc and continuum lamps
        self.meta['lampstat01'] = dict(ext=0, card=None, compound=True)

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
            binspatial = headarr[0]['CCDBIN1']
            binspec = headarr[0]['CCDBIN2']
            return parse.binning2string(binspatial, binspec)
            #return parse.binning2string(binspec, binspatial)
        elif meta_key == 'mjd':
            """
            Need to combine 'DATE-OBS' and 'UT' headers and then use astropy to make an mjd.
            """
            date = headarr[0]['DATE-OBS']
            ut = headarr[0]['UT']
            ttime = Time(f"{date}T{ut}", format='isot')
            return ttime.mjd
        elif meta_key == 'lampstat01':
            """
            If the comparison mirror is in, there will be a 'COMPLAMP' header entry containing the lamps
            that are turned on. However, if the comparison mirror is out, then this header entry doesn't exist.
            So need to test for it and set to 'Off' if it's not there.
            """
            if 'COMPLAMP' in headarr[0]:
                return headarr[0]['COMPLAMP']
            else:
                return 'off'
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
        return ['dispname', 'decker', 'binning', 'dispangle']

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
        return ['DISPERSE', 'APERTURE', 'CCDBIN1', 'CCDBIN2', 'TILTPOS']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['slitwid']


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
        # TODO: Could this be detector dependent??
        binning = '1,1' if hdu is None \
                    else self.get_meta_value(self.get_headarr(hdu), 'binning')

        # Detector 1
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = False,
            spatflip        = False,
            #platescale      = 15.0/18.0,
            platescale      = 0.2,
            darkcurr        = 5.4,
            saturation      = 65535.,
            nonlinear       = 1.0,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(1.5),
            ronoise         = np.atleast_1d(3.0),
            datasec         = np.atleast_1d('[:,1:1200]')
            #datasec         = np.atleast_1d('[250:650,1:1200]'),
            )
        return detector_container.DetectorContainer(**detector_dict)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        # Turn off illumflat
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False,
                        use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)
        # TODO: Note this comment doesn't match up with what's actually done...
        # Require dark images to be subtracted from the flat images used for
        # tracing, pixelflats, and illumflats
        par['calibrations']['traceframe']['process']['use_darkimage'] = False
        par['calibrations']['pixelflatframe']['process']['use_darkimage'] = False
        par['calibrations']['illumflatframe']['process']['use_darkimage'] = False

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 50.0

        # JFH Is this correct?
        # Processing steps
        #turn_off = dict(use_overscan=False)
        #par.reset_all_processimages_par(**turn_off)

        # Turn off the overscan
        #for ftype in par['calibrations'].keys():
        #    try:
        #        par['calibrations'][ftype]['process']['overscan'] = 'none'
        #    except (TypeError, KeyError):
        #        pass
        par['scienceframe']['process']['use_overscan'] = False
        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = False
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'holy-grail'
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['lamps'] = ['NeI', 'ArI', 'ArII', 'HeI']
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.5
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['fwhm']= 5.0

        #par['calibrations']['wavelengths']['n_first'] = 3
        #par['calibrations']['wavelengths']['n_final'] = 5
        #par['calibrations']['wavelengths']['sigdetect'] = 10.0
        #par['calibrations']['wavelengths']['wv_cen'] = 4859.0
        #par['calibrations']['wavelengths']['disp'] = 0.2
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = False
        par['reduce']['skysub']['global_sky_std']  = False

        par['reduce']['extraction']['sn_gauss'] = 4.0
        par['reduce']['findobj']['snr_thresh'] = 5.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [5,5]

        # cosmic ray rejection parameters for science frames
        par['scienceframe']['process']['sigclip'] = 5.0
        par['scienceframe']['process']['objlim'] = 2.0

        # Sensitivity function parameters
        par['sensfunc']['polyorder'] = 7

        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

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
                Processed bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """

        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)

        if det == 1:
            msgs.info("Using hard-coded BPM for Bok B&C")

            bpm_img[:, -1] = 1

        else:
            msgs.error(f"Invalid detector number, {det}, for Bok B&C (only one detector).")

        return bpm_img

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

        # Wavelength calibrations
        if self.get_meta_value(scifile, 'dispname') == '300':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'bok_bc_300.fits'
            par['calibrations']['wavelengths']['method'] = 'full_template'

        return par

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
        if ftype in ['science']:
            science = []
            for obj in fitstbl['idname'].tolist():
                science.append(not (("Dome Flat" in obj) or ("STANDARD" in obj) or ("HIP05" in obj) or ("HZ" in obj) or ("G191" in obj) or ("PG0220" in obj)))
            return good_exp & (fitstbl['lampstat01'] == 'off') & science
        if ftype in ['standard']:
            standard = []
            for obj in fitstbl['idname'].tolist():
                standard.append(("STANDARD" in obj) or ("HIP05" in obj) or ("HZ" in obj) or ("G191" in obj) or ("PG0220" in obj) )
            return good_exp & (fitstbl['lampstat01'] == 'off') & (fitstbl['idname'] != 'Dome Flat') & standard
        if ftype == 'bias':
            bias = []
            for obj in fitstbl['idname'].tolist():
                bias.append(("BIAS" in obj) or ("Bias" in obj))
            return good_exp & (fitstbl['lampstat01'] == 'off') & bias
        if ftype in ['pixelflat', 'trace']:
            flat = []
            for obj in fitstbl['idname'].tolist():
                flat.append(("Dome Flat" in obj))
            return good_exp & (fitstbl['lampstat01'] == 'off') & flat
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['lampstat01'] != 'off')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)



