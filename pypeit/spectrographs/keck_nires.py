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
                            platescale      = 0.15,
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
        return 'Echelle'

    def default_pypeit_par(self):
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
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grating dependent..
#        par['calibrations']['wavelengths']['min_nsig'] = 10.0
        par['calibrations']['wavelengths']['sigdetect'] =5.0
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.json'
        par['calibrations']['wavelengths']['echelle'] = True


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
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = 'skip'
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
        #hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['utc'] = 'UTC'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['exptime'] = 'ITIME'
        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'
        hdr_keys[0]['binning'] = 1
        hdr_keys[0]['dispname'] = 'INSTR'  # Should be 'spec' if in the spectroscopy mode

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
        """Set the general matching criteria for NIRES"""
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'

        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'

        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'

        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'

        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'

        # OLD
        # Bias
        #match_criteria['bias']['match'] = {}
        #match_criteria['standard']['match'] = {}
        #match_criteria['pixelflat']['match'] = {}
        #match_criteria['trace']['match'] = {}
        #match_criteria['arc']['match'] = {}
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


    def slit2order(self,islit):

        """
        Parameters
        ----------
        islit: int, float, or string, slit number

        Returns
        -------
        order: int
        """

        orders = [3,4,5,6,7]
        return orders[int(islit)]

