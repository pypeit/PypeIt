""" Module for Keck/NIRES specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class KeckNIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRESSpectrograph, self).__init__()
        self.spectrograph = 'keck_nires'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRES'
        self.numhead = 3
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dispaxis        = 1,
                            dispflip        = True,
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
                            )]
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @property
    def pypeline(self):
        return 'MultiSlit'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_nires'
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
        par['calibrations']['tilts']['tracethresh'] = [10, 10, 10, 10, 10]
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 3.
        par['calibrations']['slits']['pcatype'] = 'order'
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 20]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a Keck NIRES exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'NIRES',
                               '1.NAXIS': 2,
                              '1.NAXIS1': 2048,
                              '1.NAXIS2': 1024 }
        super(KeckNIRESSpectrograph, self).check_headers(headers, expected_values=expected_values)

    def header_keys(self):
        """
        Return a dictionary with the header keywords to read from the
        fits file.

        Returns:
            dict: A nested dictionary with the header keywords to read.
            The first level gives the extension to read and the second
            level gives the common name for header values that is passed
            on to the PypeItMetaData object.
        """
        hdr_keys = {}
        hdr_keys[0] = {}

        # Copied over defaults
        hdr_keys[0]['idname'] = 'OBSTYPE'
        hdr_keys[0]['time'] = 'MJD-OBS'
        hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['utc'] = 'UTC'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['exptime'] = 'ITIME'
        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'

        return hdr_keys

    def metadata_keys(self):
        return ['filename', 'date', 'frametype', 'target', 'exptime']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace']:
            return fitstbl['idname'] == 'domeflat'
        
        return (fitstbl['idname'] == 'object') \
                        & framematch.check_frame_exptime(fitstbl['exptime'], exprng)
  
    def get_match_criteria(self):
        """Set the general matching criteria for Shane Kast."""
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}

        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['binning'] = ''

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
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_nsig'] = 50.0         # Min significance for arc lines to be used
        arcparam['lowest_nsig'] = 10.0         # Min significance for arc lines to be used
        arcparam['wvmnx'] = [8000.,26000.]  # Guess at wavelength range
        # These parameters influence how the fts are done by pypeit.core.wavecal.fitting.iterative_fitting
        arcparam['match_toler'] = 3         # Matcing tolerance (pixels)
        arcparam['func'] = 'legendre'       # Function for fitting
        arcparam['n_first'] = 2             # Order of polynomial for first fit
        arcparam['n_final'] = 4             # Order of polynomial for final fit
        arcparam['nsig_rej'] = 2            # Number of sigma for rejection
        arcparam['nsig_rej_final'] = 3.0    # Number of sigma for rejection (final fit)


#        arcparam['llist'] = ''
#        arcparam['disp'] = 2.              # Ang/unbinned pixel
#        arcparam['b1'] = 0.                # Pixel fit term (binning independent)
#        arcparam['b2'] = 0.                # Pixel fit term
#        arcparam['wv_cen'] = 0.            # Estimate of central wavelength
#        arcparam['wvmnx'] = [9000., 25000.] # Guess at wavelength range
#        arcparam['disp_toler'] = 0.1       # 10% tolerance
#        arcparam['match_toler'] = 3.       # Matching tolerance (pixels)
#        arcparam['func'] = 'legendre'      # Function for fitting
#        arcparam['n_first'] = 1            # Order of polynomial for first fit
#        arcparam['n_final'] = 3            # Order of polynomial for final fit
#        arcparam['nsig_rej'] = 2.          # Number of sigma for rejection
#        arcparam['nsig_rej_final'] = 2.0   # Number of sigma for rejection (final fit)
#        arcparam['Nstrong'] = 13           # Number of lines for auto-analysis

