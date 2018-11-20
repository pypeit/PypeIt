'''
Implements HIRES-specific functions, including reading in slitmask design files.
'''
from __future__ import absolute_import, division, print_function

import glob
import re
import numpy as np

from scipy import interpolate

from astropy.io import fits

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit.spectrographs.slitmask import SlitMask
from pypeit.spectrographs.opticalmodel import ReflectionGrating, OpticalModel, DetectorMap

from pypeit import debugger


class KECKHIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle KECK/HIRES specific code
    """
    def __init__(self):
        # Get it started
        super(KECKHIRESSpectrograph, self).__init__()
        self.spectrograph = 'keck_hires_base'
        self.telescope = telescopes.KeckTelescopePar()

    @property
    def pypeline(self):
        return 'Echelle'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for KECK HIRES reductions.
        """
        par = pypeitpar.PypeItPar()
        # Correct for flexure using the default approach
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def header_keys(self):
        hdr_keys = {}
        hdr_keys[0] = {}

        # The keyword that identifies the frame type (i.e. bias, flat, etc.)
        hdr_keys[0]['idname']  = 'OBSTYPE'
        # Header keyword for the name given by the observer to a given frame
        hdr_keys[0]['target']  = 'OBJECT'
        hdr_keys[0]['utc'] = 'UTC'
        # The UT date of the observation which is used for heliocentric
        # (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)
        hdr_keys[0]['date']    = 'DATE-OBS'
        # Right Ascension of the target
        hdr_keys[0]['ra']      = 'RA'
        # Declination of the target
        hdr_keys[0]['dec']     = 'DEC'
        # Airmass at start of observation
        hdr_keys[0]['airmass'] = 'AIRMASS'
        # Exposure time keyword
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['binning'] = 'BINNING'
        return hdr_keys

    def metadata_keys(self):
        return ['filename', 'date', 'frametype', 'idname', 'target', 'exptime', 'decker',
                'binning']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Object')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'Std') | (fitstbl['idname'] == 'Object'))
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'Bias')
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'Dark')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['idname'] == 'Flat') | (fitstbl['idname'] == 'IntFlat'))
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'Line')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_match_criteria(self):
        # TODO: Matching needs to be looked at...
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        #
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['binning'] = ''
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['binning'] = ''
        # Pixelflat
        match_criteria['pixelflat']['match'] = match_criteria['standard']['match'].copy()
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['standard']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = match_criteria['bias']['match'].copy()

        # Return
        return match_criteria


class KECKHIRESRSpectrograph(KECKHIRESSpectrograph):
    """
    Child to handle KECK/HIRES-R specific code
    """
    def __init__(self):
        # Get it started
        super(KECKHIRESRSpectrograph, self).__init__()
        self.spectrograph = 'keck_hires_red'
        self.camera = 'HIRES_R'
        self.detector = [
            # Detector 1 B
            pypeitpar.DetectorPar(dataext         = 1,
                        dispaxis        = 0,  # Device is fussed with by the image reader
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.191,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 1,
                        gain            = 0.78, # high gain, low gain 1.9
                        ronoise         = 2.8,
                        suffix          = '_01'
                        ),
            # Detector 2
            pypeitpar.DetectorPar(dataext         = 2,
                        dispaxis        = 0,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.191,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 1,
                        gain            = 0.86, # high gain, low gain 2.2
                        ronoise         = 3.1,
                        suffix          = '_02'
                        ),
            # Detector 3
            pypeitpar.DetectorPar(dataext         = 3,
                        dispaxis        = 0,
                        xgap            = 0.,
                        ygap            = 0.,
                        ysize           = 1.,
                        platescale      = 0.191,
                        darkcurr        = 0.0,
                        saturation      = 65535.,
                        nonlinear       = 0.86,
                        numamplifiers   = 1,
                        gain            = 0.84, # high gain, low gain 2.2
                        ronoise         = 3.1,
                        suffix          = '_03'
                        ),
        ]
        self.numhead = 4

    def default_pypeit_par(self):
        """
        Set default parameters for HIRES RED reductions.
        """
        par = KECKHIRESSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'keck_hires_red'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slits']['sigdetect'] = 600.
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['pcatype'] = 'pixel'
        par['calibrations']['tilts']['tracethresh'] = 20

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        # Reidentification parameters
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.json'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        par['calibrations']['standardframe']['exprng'] = [None, 600]
        par['scienceframe']['exprng'] = [600, None]

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an KECK/HIRES-R exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'HIRES: High Resolution Echelle Spectrometer',
                            '0.XDISPERS': 'RED'}
        super(KECKHIRESRSpectrograph, self).check_headers(headers,
                                                              expected_values=expected_values)

    def header_keys(self):
        hdr_keys = super(KECKHIRESRSpectrograph, self).header_keys()
        hdr_keys[0]['decker'] = 'DECKNAME'
        return hdr_keys


    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """
        Override parent bpm function with BPM specific to X-ShooterNIR.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """

        self.empty_bpm(shape=shape, filename=filename, det=det)
        return self.bpm_img


