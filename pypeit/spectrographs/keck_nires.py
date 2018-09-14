""" Module for Keck/NIRES specific codes
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

# TODO: KeckNIRESSPectrograph instead (i.e. two SS) ?
class KeckNIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRESpectrograph, self).__init__()
        self.spectrograph = 'keck_nires'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRES'
        self.numhead = 3
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = -1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.01,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 3.8,
                            ronoise         = 5.0,
                            datasec         = '[1:2048,1:1024]',
                            oscansec        = '[1:2048,980:1024]'
                            )
            ]
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast Blue reductions.
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
        par['calibrations']['tilts']['tracethresh'] = [50, 50, 60, 60, 2000]
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 3.
        par['calibrations']['slits']['pcatype'] = 'order'
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = None
        return par

    def nires_header_keys(self):
        def_keys = self.default_header_keys()

        def_keys[0]['target'] = 'OBJECT'
        def_keys[0]['exptime'] = 'ITIME'
        return def_keys

    def header_keys(self):
        head_keys = self.nires_header_keys()
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
        msgs.info("Custom bad pixel mask for NIRES")
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

        arcparam['lamps'] = ['OH_triplespec'] # Line lamps on
        arcparam['llist'] = ''
        arcparam['disp'] = 2.              # Ang/unbinned pixel
        arcparam['b1'] = 0.                # Pixel fit term (binning independent)
        arcparam['b2'] = 0.                # Pixel fit term
        arcparam['wv_cen'] = 0.            # Estimate of central wavelength
        arcparam['wvmnx'] = [9000., 25000.] # Guess at wavelength range
        arcparam['disp_toler'] = 0.1       # 10% tolerance
        arcparam['match_toler'] = 3.       # Matching tolerance (pixels)
        arcparam['min_ampl'] = 1000.       # Minimum amplitude
        arcparam['func'] = 'legendre'      # Function for fitting
        arcparam['n_first'] = 1            # Order of polynomial for first fit
        arcparam['n_final'] = 3            # Order of polynomial for final fit
        arcparam['nsig_rej'] = 2.          # Number of sigma for rejection
        arcparam['nsig_rej_final'] = 2.0   # Number of sigma for rejection (final fit)
        arcparam['Nstrong'] = 13           # Number of lines for auto-analysis


        """ OLD VERSION
        arcparam = dict(llist='',
                        disp=2.,  # Ang/unbinned pixel
                        b1=0.,  # Pixel fit term (binning independent)
                        b2=0.,  # Pixel fit term
                        lamps=['OH_triplespec'],  # Line lamps on
                        wv_cen=0.,  # Estimate of central wavelength
                        wvmnx=[9000., 25000.],  # Guess at wavelength range
                        disp_toler=0.1,  # 10% tolerance
                        match_toler=3.,  # Matching tolerance (pixels)
                        min_ampl=1000.,  # Minimum amplitude
                        func='legendre',  # Function for fitting
                        n_first=1,  # Order of polynomial for first fit
                        n_final=3,  # Order of polynomial for final fit
                        nsig_rej=2.,  # Number of sigma for rejection
                        nsig_rej_final=2.0,  # Number of sigma for rejection (final fit)
                        Nstrong=13)  # Number of lines for auto-analysis 
    """

        # NEW VERSION
