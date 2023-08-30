"""
Module for MDM/Modspec specific methods.

.. include:: ../include/links.rst
"""

import array as arr
import numpy as np

from astropy.time import Time

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
    name = 'mdm_modspec'

    telescope = telescopes.HiltnerTelescopePar()

    camera = 'Echelle'
    header_name = 'Modspec'
    pypeline = 'MultiSlit'
    supported = True
    comment = 'MDM Modspec spectrometer; Only 1200l/mm disperser (so far)'

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

        # Allowing hdu=None is only needed for the automated documentation.
        # TODO: See if there's a better way to automatically create the detector
        # table for the docs.
        if hdu is None:
            lenSpat = None
            lenSpec = None
            datasec = None
            oscansec = None
            binning = '1,1'                 # Most common use mode
        else:
            # length of spatial axis, including overscan. Horizontal axis of
            # original .fits files
            lenSpat = hdu[0].header['NAXIS1']
            # length of spectral axis. Vertical axis of original .fits files
            lenSpec = hdu[0].header['NAXIS2']
            datasec = np.atleast_1d([f'[1:{lenSpec},1:3002]'])
            oscansec = np.atleast_1d([f'[1:{lenSpec},308:{lenSpat}]'])
            binning = self.compound_meta(self.get_headarr(hdu), 'binning')

        if binning != '1,1':
            msgs.error("Not ready for any binning except 1x1;  contact the developers")

        # Detector 1 continued
        detector_dict = dict(
            binning         = binning,
            det             = 1,
            dataext         = 0,
            specaxis        = 0,        # Native spectrum is along the x-axis 
            specflip        = True,     # Wavelength decreases as pixel number increases
            spatflip        = False,    # Spatial position increases as pixel number increases
            platescale      = 0.28,     # Arcsec / pixel
            darkcurr        = 0.0,      # e-/pixel/hour
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

        # Slit edge method
        par['calibrations']['slitedges']['sync_predict'] = 'nearest' # Ignore PCA
        par['calibrations']['slitedges']['bound_detector'] = True # Edges of slit fall off the detector, so assign edges of detector as the edges of the slit

        # Pixel flat method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'mean'
        par['calibrations']['pixelflatframe']['process']['clip'] = True
        par['calibrations']['pixelflatframe']['process']['comb_sigrej'] = 3.0 
        par['calibrations']['pixelflatframe']['process']['n_lohi'] = [1, 1] 
        par['calibrations']['pixelflatframe']['process']['use_overscan'] = True  
        
        # Wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'full_template' 
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'XeI', 'NeI']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'mdm_modspec_1200_5100.fits'
        par['calibrations']['wavelengths']['sigdetect'] = 5.0 # Sigma threshold above fluctuations for arc-line detection
        #par['calibrations']['wavelengths']['ech_fix_format'] = False
        par['calibrations']['wavelengths']['n_final'] = 9

        # Flat fielding
        par['calibrations']['flatfield']['slit_illum_finecorr'] = False
        
        # Bias method
        par['calibrations']['biasframe']['process']['overscan_method'] = 'median'
        
        # Arc method
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['arcframe']['process']['combine'] = 'mean'
        
        # Tilt method
        par['calibrations']['tiltframe']['process']['subtract_continuum'] = True
        par['calibrations']['tiltframe']['process']['clip'] = False
        par['calibrations']['tiltframe']['process']['combine'] = 'mean'

        
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['scienceframe']['exprng'] = [10, 600]

        return par

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the ``PypeIt``-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        # see https://pypeit.readthedocs.io/en/latest/setup.html?highlight=decker#overview
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        self.meta['decker'] = dict(card=None, compound=True)
        self.meta['binning'] = dict(card=None, compound=True)
        
        #self.meta['datasec'] = dict(ext=0, card='DATASEC') # possibly use local variable for this
        self.meta['filter1'] = dict(ext=0, card='FILTER') 

        
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        
        # Extras for config and frametyping
        self.meta['dispname'] = dict(card=None, compound=True) 
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
        self.meta['cenwave'] = dict(card=None, compound=True, rtol=2.0)

       
        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPS')
        self.meta['instrument'] = dict(ext=0, card='INSTRUME')
        
        # Mirror
        self.meta['mirror'] = dict(ext=0, card='MIRROR')
        
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
            return parse.binning2string(binspec, binspatial)
        if meta_key == 'mjd':
            return Time(headarr[0]['JD'], format='jd').mjd
        if meta_key == 'decker':
            return 'none'
        if meta_key == 'dispname':
            return '1200 l/mm'
        if meta_key == 'cenwave':
            return 5100.0
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
        return ['dispname', 'cenwave', 'filter1', 'binning']

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
                This table uses the Pypeit-specific metadata keywords, as defined
                under def init_meta(self).
            exprng (:obj:`list`, optional):
                Range in the allowed exposure time for a frame of type
                ``ftype``. See
                :func:`pypeit.core.framematch.check_frame_exptime`.

        Returns:
            `numpy.ndarray`_: Boolean array with the flags selecting the
            exposures in ``fitstbl`` that are ``ftype`` type frames.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science']:                # Standards and Sciences lumped together under 'science'
            return good_exp & (fitstbl['idname'] == 'Object') & (fitstbl['mirror'] == 'OUT')
        
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        
        if ftype in ['arc','tilt']:             # Lamps and Arcs
            return good_exp & np.array([target in ['Comp','Object'] for target in fitstbl['idname']]) & (fitstbl['mirror'] == 'IN')
        
        if ftype in ['pixelflat']:              # Internal Flats
            return good_exp & (fitstbl['idname'] == 'Flat') & (fitstbl['mirror'] == 'IN')
                
        if ftype in ['illumflat', 'trace']:     # Twilight Flats
            return good_exp & (fitstbl['idname'] == 'Flat') & (fitstbl['mirror'] == 'OUT')
        
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))

        return np.zeros(len(fitstbl), dtype=bool)
    