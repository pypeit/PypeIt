"""
Module for NTT EFOSC2

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from IPython import embed

class NTTEFOSC2Spectrograph(spectrograph.Spectrograph):
    """
    Child of Spectrograph to handle NTT/EFOSC2 specific code
    """
    ndet = 1  # Because each detector is written to a separate FITS file
    telescope = telescopes.NTTTelescopePar()
    name = 'ntt_efosc2'
    header_name = 'EFOSC'
    camera = 'EFOSC2'
    url = 'https://www.eso.org/sci/facilities/lasilla/instruments/efosc.html'
    supported = True
    comment = 'The ESO Faint Object Spectrograph and Camera version 2'

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA', required_ftypes=['science', 'standard'])
        self.meta['dec'] = dict(ext=0, card='DEC', required_ftypes=['science', 'standard'])
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['binning'] = dict(card=None, compound=True) #CDELT1 and CDELT2
        self.meta['mjd'] = dict(ext=0, card='MJD-OBS')
        
        self.meta['datasec'] = dict(card=None, compound=True)
        self.meta['oscansec'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        self.meta['decker'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='HIERARCH ESO INS GRIS1 NAME', required_ftypes=['science', 'standard'])
        self.meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
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
            binspatial = headarr[0]['CDELT1']
            binspec = headarr[0]['CDELT2']
            binning = parse.binning2string(int(binspec), int(binspatial))
            return binning
        elif meta_key == 'decker':
            try:  # Science
                decker = headarr[0]['HIERARCH ESO INS SLIT1 NAME']
            except KeyError:  # Standard!
                try:
                    decker = headarr[0]['HIERARCH ESO SEQ SPEC TARG']
                except KeyError:
                    return None
            return decker
        elif meta_key == 'datasec' or meta_key == 'oscansec':
            xbin =  int(headarr[0]['CDELT2'])
            data_x = int(headarr[0]['HIERARCH ESO DET OUT1 NX'] * headarr[0]['CDELT1']) #valid pixels along X
            data_y = int(headarr[0]['HIERARCH ESO DET OUT1 NY'] * headarr[0]['CDELT2']) #valid pixels along Y
            oscan_y = int(headarr[0]['HIERARCH ESO DET OUT1 OVSCY'] * headarr[0]['CDELT1']) #Overscan region in Y, no overscan in X
            pscan_x = int(headarr[0]['HIERARCH ESO DET OUT1 PRSCX'] * headarr[0]['CDELT2']) #Prescan region in X, no prescan in Y  
            pscan_y = int(headarr[0]['HIERARCH ESO DET OUT1 PRSCY'] * headarr[0]['CDELT2']) #Prescan region in Y
            oscan_x = int(headarr[0]['HIERARCH ESO DET OUT1 X'])  # X location of output;  Not binned
            max_x = int(headarr[0]['NAXIS1'] * headarr[0]['CDELT2']) # Maximum columns
            if meta_key == 'datasec':
                datasec = '[%s:%s,:%s]' % (pscan_y+1, pscan_y+data_y, 
                                           data_x)
                return datasec
            else:
                oscansec = '[%s:%s,%s:%s]' % (pscan_y+1, pscan_y+data_y,
                                              oscan_x+1*xbin, max_x-1*xbin) # Actually two overscan regions, here I only dealing with the region on x-axis
                return oscansec
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
        return ['dispname', 'decker', 'binning', 'datasec']

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
        return ['HIERARCH ESO INS GRIS1 NAME', 'HIERARCH ESO INS SLIT1 NAME',
                'HIERARCH ESO SEQ SPEC TARG', 'CDELT1', 'CDELT2',
                'HIERARCH ESO DET OUT1 NX', 'HIERARCH ESO DET OUT1 NY',
                'HIERARCH ESO DET OUT1 OVSCY', 'HIERARCH ESO DET OUT1 PRSCX',
                'HIERARCH ESO DET OUT1 PRSCY', 'HIERARCH ESO DET OUT1 X']

    def get_detector_par(self, det, hdu=None):
        """
        Return metadata for the selected detector.

        .. warning::

            Many of the necessary detector parameters are read from the file
            header, meaning the ``hdu`` argument is effectively **required** for
            NOT/EFOSC2.  The optional use of ``hdu`` is only viable for
            automatically generated documentation.

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
        if hdu is None:
            binning = '1,1'
            datasec = None
            oscansec = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            datasec = np.atleast_1d(self.get_meta_value(self.get_headarr(hdu), 'datasec'))
            oscansec = np.atleast_1d(self.get_meta_value(self.get_headarr(hdu), 'oscansec'))
        
        # Manual: https://www.eso.org/sci/facilities/lasilla/instruments/efosc/doc/manual/EFOSC2manual_v4.2.pdf
        # Instrument paper: http://articles.adsabs.harvard.edu/pdf/1984Msngr..38....9B
        detector_dict = dict(
            binning         = binning,
            det             = 1, # only one detector
            dataext         = 0,
            specaxis        = 0,
            specflip        = False,
            spatflip        = False,
            platescale      = 0.12, # Manual 2.2
            darkcurr        = 0.0,
            saturation      = 65535, # Maual Table 8
            nonlinear       = 0.80,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(0.91), # See fits header ['HIERARCH ESO DET OUT1 GAIN']
            ronoise         = np.atleast_1d(10.0), # manual page 108
            datasec         = datasec,
            oscansec         = oscansec,
            #suffix          = '_Thor',
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

        # Always correct for flexure, starting with default parameters
        par['flexure']['spec_method'] = 'boxcar'
        
        # Adjustments to slit and tilts for NIR
        par['calibrations']['traceframe']['process']['use_darkimage'] = False
        par['calibrations']['pixelflatframe']['process']['use_darkimage'] = False
        par['calibrations']['illumflatframe']['process']['use_darkimage'] = False
        
        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        #par['calibrations']['slitedges']['rm_slits'] = '1:500:120' # remove the fake slit due to bad pixels

        #edge parameters
        par['calibrations']['slitedges']['edge_thresh'] = 75.

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] = 25.0
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 4
        
        # Image processing
        # The overscan region might cause oversubtraction of the background, set it to False
        #par['scienceframe']['process']['use_overscan'] = False 

        # 1D wavelength solution
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['lamps'] = ['HeI', 'ArI']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 4.0
        par['calibrations']['wavelengths']['n_final'] = 4

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10
        
        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['skysub']['no_poly'] = True
        par['reduce']['skysub']['bspline_spacing'] = 0.6
        par['reduce']['skysub']['joint_fit'] = False
        par['reduce']['skysub']['global_sky_std']  = False

        par['reduce']['extraction']['sn_gauss'] = 4.0
        par['reduce']['skysub']['sky_sigrej'] = 5.0

        return par

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
        if self.get_meta_value(scifile, 'dispname') == 'Gr#6':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ntt_efosc2_Gr6.fits'
        elif self.get_meta_value(scifile, 'dispname') == 'Gr#5':
            par['calibrations']['wavelengths']['reid_arxiv'] = 'ntt_efosc2_Gr5.fits'
            # Fringes are affecting this Grism significantly, skip flat fielding
            par['scienceframe']['process']['use_pixelflat'] = False
            par['scienceframe']['process']['use_illumflat'] = False
            par['scienceframe']['process']['use_specillum'] = False

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
            return good_exp & ((fitstbl['target'] == 'BIAS')
                               |(fitstbl['target'] == 'DARK'))
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'FLAT')
                               | (fitstbl['target'] == 'SKY,FLAT')
                               | (fitstbl['target'] == 'DOME'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['target'] == 'WAVE'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
    
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

        msgs.info("Using hard-coded BPM for NTT EFOSC2")
        binning = self.get_meta_value(filename, 'binning')
        binspatial =  int(binning[0])
        binspec =  int(binning[2])
        bpm_img[int(232/binspec):, int(362/binspatial):int(366/binspatial)] = 1
        bpm_img[int(340/binspec):, int(1292/binspatial)] = 1
        #bpm_img[int(2050/binspec):, :] = 1

        return bpm_img
