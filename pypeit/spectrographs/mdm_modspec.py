"""
Module for MDM/Modspec specific methods.

.. include:: ../include/links.rst
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container

class MDMModspecEchelleSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle MDM Modspec Echelle instrument+detector
    """
    ndet = 1
    name = 'mdm_modspec_echelle'
    telescope = telescopes.KPNOTelescopePar()
    camera = 'Echelle'
    header_name = 'ModSpec'
    supported = True
    comment = 'MDM Modspec spectrometer'

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
        # Detector 1
        # See Echelle at 2.4m f/7.5 scale : http://mdm.kpno.noirlab.edu/mdm-ccds.html 
        gain = np.atleast_1d([1.3])      # Hardcoded in the header 
        ronoise = np.atleast_1d([7.90])    # Hardcoded in the header
        len1 = hdu[0].header['NAXIS1']
        len2 = hdu[0].header['NAXIS2']
    
        datasec = np.atleast_1d([
            '[{0:d}:{1:d},{2:d}:{3:d}]'.format(1+5, len1-5, 1, len2)])
        oscansec = np.atleast_1d([
            '[{0:d}:{1:d},{2:d}:{3:d}]'.format(1, 1+5, 1, len2),
            '[{0:d}:{1:d},{2:d}:{3:d}]'.format(len1-5, len1, 1, len2)
        ])
        if hdu is None:
            binning = '1,1'                 # Most common use mode
        else:
            binning = "{0},{1}".format(
                hdu[0].header['CCDBIN1'], hdu[0].header['CCDBIN2'])

        # Detector 1 continued
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,        # Native spectrum is along the x-axis
            specflip        = True,     # DeVeny CCD has blue at the right
            spatflip        = False,
            platescale      = 0.28,     # Arcsec / pixel
            darkcurr        = 0.0,      # Electrons per hour
            saturation      = 65535.,   # 16-bit ADC
            nonlinear       = 0.97,     # Linear to ~97% of saturation
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = gain,     # See above
            ronoise         = ronoise,  # See above
            # Data & Overscan Sections -- Edge tracing can handle slit edges
            datasec         = datasec,  # See above
            oscansec        = oscansec  # See above
            )
        # Return
        return detector_container.DetectorContainer(**detector_dict)
        
    @classmethod

    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of ``PypeIt`` methods.
        """
        par = super().default_pypeit_par()

        # Ignore PCA
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'mean'
        par['calibrations']['pixelflatframe']['process']['clip'] = True
        par['calibrations']['pixelflatframe']['process']['comb_sigrej'] = 3.0 
        par['calibrations']['pixelflatframe']['process']['n_lohi'] = [1, 1] #[nlow, nhigh]
        par['calibrations']['pixelflatframe']['process']['use_overscan'] = False
        
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'full_template' #more reliable than 'holy-grail', but requires an archived wavelength solution for the specific instrument/grating combination. See https://pypeit.readthedocs.io/en/latest/pypeit_par.html#wavelengthsolutionpar-keywords, also https://pypeit.readthedocs.io/en/latest/wave_calib.html#identify and https://pypeit.readthedocs.io/en/latest/master_edges.html and https://pypeit.readthedocs.io/en/latest/master_arc.html
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'XeI', 'NeI']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'mdm_modspec.fits' #abovementioned archived wavelength solution; need one for Echelle / Modspec
        ###|||||| do this one below ||||||###
        par['calibrations']['wavelengths']['sigdetect'] = 10.0 #Sigma threshold above fluctuations for arc-line detection
        
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['standardframe']['exprng'] = [10, 60]
        par['scienceframe']['exprng'] = [120, 600]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(card=None)
        self.meta['binning'] = dict(card=None, compound=True)
        #self.meta['binning'] = [dict(ext=0, card='CCDBIN1'), dict(ext=0, card='CCDBIN2')]

        self.meta['mjd'] = float(dict(ext=0, card='JD')) - 2400000.5
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        # in an example of what this code generates, see https://pypeit.readthedocs.io/en/latest/pypeit_file.html#pypeit-file
        ## on that note, 'dispname' is showing 600/4310 and that makes no sense to me 
        ## okay apparently according to https://watermark.silverchair.com/360-4-1281.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAuYwggLiBgkqhkiG9w0BBwagggLTMIICzwIBADCCAsgGCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMYDY3odQx18n61h_OAgEQgIICmRDa6yhMHVq4rhm9vxfUOxJJk38Ppwr37LtGSH08YK5TV_Y3sUaRQ6Pd_CRDP4HpCy0UBoRIEvJMaW77Tnq1k45akgXw26T5D-SsWWWX8v3JYjgD5wJRzcM_lbQOgljRuBgffDixM-kc4kFvlzAUGdPfSAKT-0kYp9Bj3X6SX2QUKJabemNm3kuApMDpgOQYj1_JXDxH1FAzqUco3abNAJSJwXp33t8vhNaK3hTI_CO0X66ul6PYx2aGEs1MG2AFLT0rnPt_Q1NSoXnGhRV4AJVMXJF58I7oZ_0bLuN_P-e87J9VNA1JBp9CNJzqYN9Rz0TOHEDhmC05NTYmfNpbbZUW6R4x9u6s8cQiaK0J5C6m_DtZgxokt5Rg7fSEtOZFzD42i-B1Aln4BqX57o423Uhi1QkYGfk5eeReq-fKRFEy-vrsuEiAPYhLE045xQR_OiaKeQkeQEMxumZcDA4FPLHcF9W_cIBc9-Qr24h-xdKXYpkCg994hOcdOzMvJon9nAFJyyW6CAENtDDmFewW9-Ht3EwPYCR0Remb3SuddgZPCxBXoijUcjf7YPw9PRqRLEJ26K0ag10B_eAgwjkcEAwThhHGDp9EkkB8cmlwjca5uPOlfrZ_lqen0y-UC-8wMh4bdcvqUCmsMg3GoPC4q8CVgTPdPvXOhROgyCSbJi_J53RZD8CkE1K7K9dfDN7UsFPwhb31qWOw1FoG5dAv4xWUJZG_zOe503hFhWvql8J9Go2J1IcdZXIA1eRBueW7GF9SffIp0YQddpg4e8oAwiMKD3Tqed0--lZ2Oqx6fbbubfuo5ppFTDZi83VrlRbxdPU-blJYdYJDSDH5DoQ6xQZ6QBOwU0M95TkmiNb7LgJWz5k6pdDCbDqu
        ## it is Grism/grating (lines mm^-1 / blaze), whatever that means
        #self.meta['dispname'] = 'Diffraction Grating'
        self.meta['dispname'] = '1200/5000' ## check on how this needs to be formatted
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPS')
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
            ## double-check the bin1 vs bin2 assignment
            binspatial = headarr[0]['CCDBIN1']
            binspec = headarr[0]['CCDBIN2']
            return parse.binning2string(binspec, binspatial)
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
        return ['dispname', 'decker', 'binning']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard ``PypeIt`` file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys() + ['slitwid']
    
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
            return good_exp & (fitstbl['imagetyp'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp & (fitstbl['imagetyp'] == 'zero')
        if ftype in ['arc','tilt']:
            return good_exp & np.array([ilamp in ['Ar','Xe'] for ilamp in fitstbl['lampstat01']]) & (fitstbl['idname'] == 'COMP')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
