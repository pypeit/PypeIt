"""
Module for MDM/Modspec specific methods.

.. include:: ../include/links.rst
"""
import array as arr
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
    header_name = 'Modspec'
    pypeline = 'MultiSlit'
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
        len2 = hdu[0].header['NAXIS1']     ## switched with len1
        len1 = hdu[0].header['NAXIS2']      ## switched with len2
    
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
            specflip        = False,     ## ADD COMMENT
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
        # Edges of slit fall off the detector, so assign the edges of the detector as the edges of the slit
        par['calibrations']['slitedges']['bound_detector'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'mean'
        par['calibrations']['pixelflatframe']['process']['clip'] = True
        par['calibrations']['pixelflatframe']['process']['comb_sigrej'] = 3.0 
        par['calibrations']['pixelflatframe']['process']['n_lohi'] = [1, 1] #[nlow, nhigh]
        par['calibrations']['pixelflatframe']['process']['use_overscan'] = False
        
        # Wavelength calibration methods
        par['calibrations']['wavelengths']['method'] = 'full_template' #more reliable than 'holy-grail', but requires an archived wavelength solution for the specific instrument/grating combination. See https://pypeit.readthedocs.io/en/latest/pypeit_par.html#wavelengthsolutionpar-keywords, also https://pypeit.readthedocs.io/en/latest/wave_calib.html#identify and https://pypeit.readthedocs.io/en/latest/master_edges.html and https://pypeit.readthedocs.io/en/latest/master_arc.html
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'XeI', 'NeI']
        par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_mdm_modspec_echelle_20220714T1118.fits' # this is an example; this is based only on Xenon and the minimum files needed to run
        ###|||||| do this one below ||||||###
        par['calibrations']['wavelengths']['sigdetect'] = 5.0 #Sigma threshold above fluctuations for arc-line detection
        
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, None]  # Long arc exposures on this telescope
        par['calibrations']['arcframe']['process']['clip'] = False
        par['calibrations']['arcframe']['process']['subtract_continuum'] = True
        par['calibrations']['standardframe']['exprng'] = [10, 60]
        par['scienceframe']['exprng'] = [120, 600]

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
        self.meta['decker'] = dict(card=None, compound=True) # str -- ie long_1.0 -- Name of the decker or slit mask 
        self.meta['binning'] = dict(card=None, compound=True)
        
        ##self.meta['amp'] = '' # str -- ie SINGLE:B -- Name of the amplifier used to read the detector
        ##self.meta['arm'] = '' # str -- ie VIS -- Name of the spectrograph arm used to collect the data
        self.meta['datasec'] = dict(ext=0, card='DATASEC') # str -- ie [1:256,1:512] -- The science region of the detector
        ##self.meta['detector'] = '' # str -- ie CHIP1 -- Name of the detector
        ##self.meta['dichroic'] = '' # str -- ie 560 -- Name of the dichroic
        self.meta['filter1'] = dict(ext=0, card='FILTER') # str -- ie J -- Name of the order-sorting filter
        
        self.meta['mjd'] = dict(card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        
        # Extras for config and frametyping
        # in an example of what this code generates, see https://pypeit.readthedocs.io/en/latest/pypeit_file.html#pypeit-file
        self.meta['dispname'] = dict(card=None, compound=True) # str -- ie 830G -- Name of the dispersing element
        self.meta['dispangle'] = dict(card=None, compound=True) # float -- ie 7500.0 -- Central wavelength for the dispersing element at the observed angle
        self.meta['idname'] = dict(ext=0, card='IMAGETYP')
       
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
            ## double-check the bin1 vs bin2 assignment
            binspatial = headarr[0]['CCDBIN1']
            binspec = headarr[0]['CCDBIN2']
            return parse.binning2string(binspec, binspatial)
        if meta_key == 'mjd':
            return parse.Time(headarr[0]['JD'], format='jd').mjd
        if meta_key == 'decker':
            return 'none'
        if meta_key == 'dispname':
            return '1200/5100'
        if meta_key == 'dispangle':
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
        return ['dispname', 'decker', 'binning']
        ## since there is no 'decker' card, maybe replace this with a different parameter
        ## return ['dispname', 'dispangle', 'binning']
        ## for why 'dispangle' was chosen, see https://pypeit.readthedocs.io/en/latest/setup.html?highlight=decker#overview

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
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Comp')
        if ftype in ['pixelflat']: #Internal Flats
            return good_exp & (fitstbl['idname'] == 'Flat') & (fitstbl['mirror'] == 'IN')
        if ftype in ['illumflat', 'trace']: #Twilight Flats
            return good_exp & (fitstbl['idname'] == 'Flat') & (fitstbl['mirror'] == 'OUT')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        #msgs.warn('Cannot determine if frames are of type {0}. Frame idname and target are: {1}, {2}'.format(ftype, fitstbl['idname'], fitstbl['target']))
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
                Master bias frame used to identify bad pixels

        Returns:
            `numpy.ndarray`_: An integer array with a masked value set
            to 1 and an unmasked value set to 0.  All values are set to
            0.
        """
        # Validate the entered (list of) detector(s)
        nimg, _det = self.validate_det(det)
        _det = list(_det)

        # Call the base-class method to generate the empty bpm
        bpm_img = super().bpm(filename, det, shape=shape, msbias=msbias)
        # NOTE: expand_dims does *not* copy the array.  We can edit it directly
        # because we've created it inside this function.
        _bpm_img = np.expand_dims(bpm_img, 0) if nimg == 1 else bpm_img

        if 1 in _det:
            i = _det.index(1)
            _bpm_img[i,:,1052:1054] = 1
        if 2 in _det:
            i = _det.index(2)
            _bpm_img[i,:,0:4] = 1
            _bpm_img[i,:,376:381] = 1
            _bpm_img[i,:,489] = 1
            _bpm_img[i,:,1333:1335] = 1
            _bpm_img[i,:,2047] = 1
        if 3 in _det:
            i = _det.index(3)
            _bpm_img[i,:,0:4] = 1
            _bpm_img[i,:,221] = 1
            _bpm_img[i,:,260] = 1
            _bpm_img[i,:,366] = 1
            _bpm_img[i,:,816:819] = 1
            _bpm_img[i,:,851] = 1
            _bpm_img[i,:,940] = 1
            _bpm_img[i,:,1167] = 1
            _bpm_img[i,:,1280] = 1
            _bpm_img[i,:,1301:1303] = 1
            _bpm_img[i,:,1744:1747] = 1
            _bpm_img[i,:,-4:] = 1
        if 4 in _det:
            i = _det.index(4)
            _bpm_img[i,:,0:4] = 1
            _bpm_img[i,:,47] = 1
            _bpm_img[i,:,744] = 1
            _bpm_img[i,:,790:792] = 1
            _bpm_img[i,:,997:999] = 1
        if 5 in _det:
            i = _det.index(5)
            _bpm_img[i,:,25:27] = 1
            _bpm_img[i,:,128:130] = 1
            _bpm_img[i,:,1535:1539] = 1
        if 7 in _det:
            i = _det.index(7)
            _bpm_img[i,:,426:428] = 1
            _bpm_img[i,:,676] = 1
            _bpm_img[i,:,1176:1178] = 1
        if 8 in _det:
            i = _det.index(8)
            _bpm_img[i,:,440] = 1
            _bpm_img[i,:,509:513] = 1
            _bpm_img[i,:,806] = 1
            _bpm_img[i,:,931:934] = 1

        return _bpm_img[0] if nimg == 1 else _bpm_img
