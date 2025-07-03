"""
Module for P200/NGPS specific methods.

.. include:: ../include/links.rst
"""
from typing import List, Optional

import numpy as np

from astropy.io import fits
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.spectrographs import spectrograph
from pypeit.core import parse
from pypeit.images import detector_container


class P200NGPSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/NGPS specific code
    """
    ndet = 1 
    telescope = telescopes.P200TelescopePar()

    def init_meta(self):
        """
        Define how metadata are derived from the spectrograph files.

        That is, this associates the PypeIt-specific metadata keywords
        with the instrument-specific header cards using :attr:`meta`.
        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='TELRA', required_ftypes=['science', 'standard']) 
        self.meta['dec'] = dict(ext=0, card='TELDEC', required_ftypes=['science', 'standard']) 
        self.meta['target'] = dict(ext=0, card='NAME', compound=True, required_ftypes=['science', 'standard']) 

        self.meta['dispname'] = dict(card=None, compound=True, default='VPH')
        self.meta['decker'] = dict(ext=0, card='SLITW', rtol=1e-2) 
        self.meta['binning'] = dict(card=None, compound=True) # Got by compound meta function

        self.meta['mjd'] = dict(ext=0, card='MJD')
        self.meta['exptime'] = dict(ext=0, card='SHUTTIME') # SHUTTIME more accurate than EXPTIME
        self.meta['airmass'] = dict(ext=0, card='AIRMASS', required_ftypes=['science', 'standard'])

        # Extras for config and frametyping
        self.meta['dichroic'] = dict(card=None, compound=True) 
        self.meta['dispangle'] = dict(card=None, rtol=1e-2, compound=True) # Got by compound meta function
        self.meta['slitwid'] = dict(ext=0, card='SLITW', rtol=1e-2) 
        self.meta['idname'] = dict(ext=0, card='IMGTYPE') 
        self.meta['instrument'] = dict(ext=0, card='INSTRUME') 

        # Lamps
        self.meta['lampstat01'] = dict(ext=0, card='LAMPBLUC') # Blue Xe
        self.meta['lampstat02'] = dict(ext=0, card='LAMPFEAR') # FeAr
        self.meta['lampstat03'] = dict(ext=0, card='LAMPREDC') # Red Continuum
        self.meta['lampstat04'] = dict(ext=0, card='LAMPTHAR') # ThAr

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
            object
        """
        return ['binning']
    

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
        return ['GRATING', 'ANGLE', 'APERTURE']

    def pypeit_file_keys(self):
        """
        Define the list of keys to be output into a standard PypeIt file.

        Returns:
            :obj:`list`: The list of keywords in the relevant
            :class:`~pypeit.metadata.PypeItMetaData` instance to print to the
            :ref:`pypeit_file`.
        """
        return super().pypeit_file_keys()
    
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
            return good_exp & (fitstbl['idname'] == 'SCI')
        
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'BIAS')
        
        if ftype in ['pixelflat', 'trace', 'illumflat']:
            return ((good_exp & (fitstbl['idname'] == 'DOMEFLAT')) | (good_exp & (fitstbl['idname'] == 'CONT')))
        
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)

        if ftype in ['arc', 'tilt']:
            # return good_exp & ((fitstbl['idname'] == 'FEAR') | (fitstbl['idname'] == 'THAR'))  
            return good_exp & (fitstbl['idname'] == 'THAR') # Temporary fix, do not use FEAR arcs

        
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)


class P200NGPSSpectrograph_r(P200NGPSSpectrograph):
    """
    Child to handle P200/NGPS r-Channel specific code
    """
    name = 'p200_ngps_r'
    camera = 'NGPS_r'
    header_name = 'NGPS_r'
    supported = True
    comment = 'r-Channel'

    def get_rawimage(self, raw_file, det):
        """
        Read raw spectrograph image files and return data and relevant metadata
        needed for image processing.

        For P200/NGPS, the ``DATASEC`` and ``OSCANSEC`` regions are read
        directly from the file header and are automatically adjusted to account
        for the on-chip binning.  This is a simple wrapper for
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage` that
        sets ``sec_includes_binning`` to True.  See the base-class function for
        the detailed descriptions of the input parameters and returned objects.
        """

        return super().get_rawimage(raw_file, det=1, sec_includes_binning=True)
    


    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
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
        # Handle dispangle and mjd from superclass method
        retval = super().compound_meta(headarr, meta_key)
        
        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            return 0
            
        elif meta_key == 'binning':
            # Always the same binning for DET01 and DET02    
            binspat = headarr[1]['BINSPAT'] 
            binspec = headarr[1]['BINSPEC']
            return parse.binning2string(binspec, binspat)
        
        # If there is no target keyword, return image type
        elif meta_key == 'target':
            if 'TARGET' in headarr[0]:
                return headarr[0]['TARGET']
            else:
                return headarr[0]['IMGTYPE']
        elif meta_key == 'dichroic': 
            return None
        else:
            msgs.error(f"Not ready for this compound meta: {meta_key}")


    def get_detector_par(self, det: int, hdu: Optional[fits.HDUList] = None):
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

        if hdu is None:
            binning = '1,1'
            datasec = None
            oscansec = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            datasec = np.atleast_1d(parse.flip_fits_slice(hdu[1].header['DATASEC']))
            oscansec = np.atleast_1d(parse.flip_fits_slice(hdu[1].header['BIASSEC']))

        # Detector 1 (r Channel)
        detector_dict1 = dict(
            binning         = binning,
            det             = 1, # All r channel images assigned to extension 1
            dataext         = 1, # All r channel images assigned to extension 1
            specaxis        = 1,
            specflip        = False, 
            spatflip        = False, 
            platescale      = 0.5, 
            darkcurr        = 0.0,  # e-/pixel/hour (No dark current)
            saturation      = 65000., # ???
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1, 
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = datasec,
            oscansec        = oscansec,
        )

        return detector_container.DetectorContainer(**detector_dict1)


    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 50. # Lower edge tracing thresdhold to catch leftmost slit
        par['calibrations']['slitedges']['minimum_slit_length'] = 100 # Set minimum slit length 
        par['calibrations']['slitedges']['min_edge_side_sep'] = 1.0
#        #par['calibrations']['slitedges']['add_slits'] = ['1:1090:24:153']
        
        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] = 5.0

        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = True

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'

        par['calibrations']['wavelengths']['lamps'] = ['ThAr'] # FeAr and ThAR lamps for NGPS (ThAr)
        par['calibrations']['wavelengths']['method'] = 'full_template' # Use wavelength template
        par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1227.fits' #  Channel

        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 1.0 

        # Do not flux calibrate
        #par['fluxcalib'] = None
        par['sensfunc']['algorithm'] = 'UVIS' 

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        return par



class P200NGPSSpectrograph_i(P200NGPSSpectrograph):
    """
    Child to handle P200/NGPS i-Channel specific code
    """
    name = 'p200_ngps_i'
    camera = 'NGPS_i'
    header_name = 'NGPS_i'
    supported = True
    comment = 'i-Channel'


    def get_rawimage(self, raw_file, det):
        """
        Read raw spectrograph image files and return data and relevant metadata
        needed for image processing.

        For P200/NGPS, the ``DATASEC`` and ``OSCANSEC`` regions are read
        directly from the file header and are automatically adjusted to account
        for the on-chip binning.  This is a simple wrapper for
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.get_rawimage` that
        sets ``sec_includes_binning`` to True.  See the base-class function for
        the detailed descriptions of the input parameters and returned objects.
        """

        # Pull image from detector 2
        return super().get_rawimage(raw_file, det=2, sec_includes_binning=True)
    
    def compound_meta(self, headarr: List[fits.Header], meta_key: str):
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
        # Handle dispangle and mjd from superclass method
        retval = super().compound_meta(headarr, meta_key)
        
        # If superclass could not handle the meta key
        if retval is not None:
            return retval
        
        if meta_key == 'mjd':
            return Time(headarr[0]['UTSHUT']).mjd
        elif meta_key == 'dispangle':
            return 0
            
        elif meta_key == 'binning':
            # Always the same binning for DET01 and DET02    
            binspat = headarr[2]['BINSPAT'] 
            binspec = headarr[2]['BINSPEC']
            return parse.binning2string(binspec, binspat)
        
        # If there is no target keyword, return image type
        elif meta_key == 'target':
            if 'TARGET' in headarr[0]:
                return headarr[0]['TARGET']
            else:
                return headarr[0]['IMGTYPE']
        elif meta_key == 'dichroic': 
            return None
        else:
            msgs.error("Not ready for this compound meta: ", meta_key)


    def get_detector_par(self, det: int, hdu: Optional[fits.HDUList] = None):
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

        if hdu is None:
            binning = '1,1'
            datasec = None
            oscansec = None
        else:
            binning = self.get_meta_value(self.get_headarr(hdu), 'binning')
            datasec = np.atleast_1d(parse.flip_fits_slice(hdu[2].header['DATASEC']))
            oscansec = np.atleast_1d(parse.flip_fits_slice(hdu[2].header['BIASSEC']))

        # Detector 2 (i Channel)
        detector_dict2 = dict(
            binning         = binning,
            det             = 1, # All i channel images assigned to extension 2 ###################
            dataext         = 2, # All i channel images assigned to extension 2
            specaxis        = 1,
            specflip        = False, 
            spatflip        = False, 
            platescale      = 0.5, 
            darkcurr        = 0.0,  # e-/pixel/hour (No dark current)
            saturation      = 65000., # ???
            nonlinear       = 40./45.,
            mincounts       = -1e10, # check
            numamplifiers   = 1, # Updated
            gain            = np.atleast_1d(2.8),
            ronoise         = np.atleast_1d(8.5),
            datasec         = datasec,
            oscansec        = oscansec,
        )

        return detector_container.DetectorContainer(**detector_dict2)

    @classmethod
    def default_pypeit_par(cls):
        """
        Return the default parameters to use for this instrument.
        
        Returns:
            :class:`~pypeit.par.pypeitpar.PypeItPar`: Parameters required by
            all of PypeIt methods.
        """
        par = super().default_pypeit_par()

        par['calibrations']['slitedges']['sync_predict'] = 'nearest'
        par['calibrations']['slitedges']['edge_thresh'] = 50. # Lower edge tracing thresdhold to catch leftmost slit
        par['calibrations']['slitedges']['minimum_slit_length'] = 100 # Set minimum slit length 
        par['calibrations']['slitedges']['min_edge_side_sep'] = 1.0
#        #par['calibrations']['slitedges']['add_slits'] = ['1:1015:123:251']
        
        par['scienceframe']['process']['combine'] = 'median'
        par['calibrations']['standardframe']['process']['combine'] = 'median'

        par['scienceframe']['process']['use_overscan'] = True
        par['scienceframe']['process']['sigclip'] = 4.0
        par['scienceframe']['process']['objlim'] =  5.0

        # Make a bad pixel mask
        par['calibrations']['bpm_usebias'] = False

        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'

        par['calibrations']['wavelengths']['lamps'] = ['ThAr'] # FeAr and ThAR lamps for NGPS (ThAr)
        par['calibrations']['wavelengths']['method'] = 'full_template' # Use wavelength template
        par['calibrations']['wavelengths']['reid_arxiv'] = 'wvarxiv_p200_ngps_20250131T1354.fits' # I Channel Template

        par['calibrations']['wavelengths']['rms_thresh_frac_fwhm'] = 1.0 

        # Do not flux calibrate  
        #par['fluxcalib'] = None
        par['sensfunc']['algorithm'] = 'UVIS'

        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.001]
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]

        return par

    