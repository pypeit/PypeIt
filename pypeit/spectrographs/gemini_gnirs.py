""" Module for Gemini/GNIRS specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from pypeit import msgs
from pypeit.par.pypeitpar import DetectorPar
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import fsort

from pypeit import debugger

class GeminiGNIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GNIRS specific code
    """
    def __init__(self):
        # Get it started
        super(GeminiGNIRSSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gnirs'
        self.telescope = telescopes.GeminiNTelescopePar()
        self.camera = 'GNIRS'
        self.numhead = 1
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 1,
                            dispaxis        = 2,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.15,
                            saturation      = 90000.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 13.5,
                            ronoise         = 7.0,
                            datasec         = '[:,:]',#'[1:1024,1:1022]',
                            oscansec        = '[:,:]',#'[1:1024,1:1022]'
                            )
            ]
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Gemini GNIRS reductions.
        """
        par = pypeitpar.PypeItPar()
        # TODO: Make self.spectrograph a class attribute?
        # Use the ARMS pipeline
        par['rdx']['pipeline'] = 'ARMS'
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1
        # Bias
        par['calibrations']['biasframe']['useframe'] = 'overscan'
        # Set slits and tilts parameters
        par['calibrations']['tilts']['order'] = 2
        par['calibrations']['tilts']['tracethresh'] = [1000, 50, 50, 50, 50]
        par['calibrations']['slits']['polyorder'] = 5
        #par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['pcatype'] = 'order'
        par['calibrations']['slits']['sigdetect'] = 100
        #par['calibrations']['slits']['pcapar'] = [3, 2, 1,0]
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = None
        return par

    def gnirs_header_keys(self):
        def_keys = self.default_header_keys()

        def_keys[0]['target'] = 'OBJECT'
        def_keys[0]['time'] = 'MJD_OBS'
        def_keys[0]['decker'] = 'SLIT'
        def_keys[0]['dispname'] = 'GRATING'
        return def_keys

    def header_keys(self):
        head_keys = self.gnirs_header_keys()
        print(head_keys)
        return head_keys

    def get_match_criteria(self):
        """Set the general matching criteria for Shane Kast."""
        match_criteria = {}
        for key in fsort.ftype_list:
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}

        match_criteria['pixelflat']['match'] = {}

        match_criteria['trace']['match'] = {}

        match_criteria['arc']['match'] = {}

        return match_criteria


    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
        """
        Override parent bpm function with BPM specific to X-Shooter VIS.

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
        msgs.info("Custom bad pixel mask for GNIRS")
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[:, :20] = 1.
            self.bpm_img[:, 1000:] = 1.

        return self.bpm_img



    def setup_arcparam(self, arcparam, fitstbl=None, arc_idx=None,
                       msarc_shape=None, **null_kwargs):
        """

        Args:
            arcparam:
            disperser:
            fitstbl:
            arc_idx:
            msarc_shape:
            **null_kwargs:

        Returns:

        """
        # ToDo need to parse the sigdetect parameter to be here for detect_lines function in arc.py
        #      I force change sigdetect=5 for GNIRS.

        arcparam['lamps'] = ['OH_triplespec'] # Line lamps on
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_ampl'] = 10.       # Minimum amplitude
        arcparam['wvmnx'] = [8500., 24000.] # Guess at wavelength range
        arcparam['n_first'] = 1            # Order of polynomial for first fit
        arcparam['n_final'] = 1            # Order of polynomial for final fit
        arcparam['disp'] = 2.7  # Ang/unbinned pixel
        arcparam['func'] = 'legendre'  # Function for fitting
#        arcparam['llist'] = ''
#        arcparam['disp'] = 2.              # Ang/unbinned pixel
#        arcparam['b1'] = 0.                # Pixel fit term (binning independent)
#        arcparam['b2'] = 0.                # Pixel fit term
#        arcparam['wv_cen'] = 0.            # Estimate of central wavelength
#        arcparam['disp_toler'] = 0.1       # 10% tolerance
#        arcparam['match_toler'] = 3.       # Matching tolerance (pixels)
#        arcparam['func'] = 'legendre'      # Function for fitting
#        arcparam['n_first'] = 1            # Order of polynomial for first fit
#        arcparam['n_final'] = 3            # Order of polynomial for final fit
#        arcparam['nsig_rej'] = 2.          # Number of sigma for rejection
#        arcparam['nsig_rej_final'] = 2.0   # Number of sigma for rejection (final fit)
#        arcparam['Nstrong'] = 13           # Number of lines for auto-analysis

