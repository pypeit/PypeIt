""" Module for Shane/Kast specific codes
"""
from __future__ import absolute_import, division, print_function


import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class WHTISISSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):
        super(WHTISISSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_base'
        self.telescope = telescopes.WHTTelescopePar()

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
        hdr_keys[1] = {}

        # Copied over defaults
        hdr_keys[0]['idname'] = 'IMAGETYP'
        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['time'] = 'MJD-OBS'
        # hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['utc'] = 'UT'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['binning_x'] = 'CCDXBIN'
        hdr_keys[0]['binning_y'] = 'CCDYBIN'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['decker'] = 'SLITNAME'
        hdr_keys[0]['dichroic'] = 'DICHNAME'
        hdr_keys[0]['filter1'] = 'ISIFILTA'
        hdr_keys[0]['filter2'] = 'ISIFILTB'
        hdr_keys[0]['decker'] = 'ISISLITU'
        hdr_keys[0]['slitwid'] = 'ISISLITW'
        hdr_keys[0]['dichroic'] = 'ISIDICHR'
        hdr_keys[0]['dispname'] = 'ISIGRAT'
        hdr_keys[0]['dispangle'] = 'CENWAVE'
        hdr_keys[0]['lamps'] = 'CAGLAMPS'

        hdr_keys[1]['naxis1'] = 'NAXIS1'
        hdr_keys[1]['naxis0'] = 'NAXIS2'

        return hdr_keys

    def metadata_keys(self):
        return super(WHTISISSpectrograph, self).metadata_keys() \
               + ['binning', 'dichroic', 'dispangle']

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'wht_isis_blue'
        # Set pixel flat combination method
        par['calibrations']['pixelflatframe']['process']['combine'] = 'median'
        par['calibrations']['pixelflatframe']['process']['sig_lohi'] = [10., 10.]
        # Change the wavelength calibration method
        par['calibrations']['wavelengths']['method'] = 'arclines'
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Do not flux calibrate
        par['fluxcalib'] = None
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['arcframe']['exprng'] = [None, 120]
        par['calibrations']['standardframe']['exprng'] = [None, 120]
        par['scienceframe']['exprng'] = [90, None]
        return par


class WHTISISBlueSpectrograph(WHTISISSpectrograph):
    """
    Child to handle WHT/ISIS blue specific code
    """
    def __init__(self):
        # Get it started
        super(WHTISISBlueSpectrograph, self).__init__()
        self.spectrograph = 'wht_isis_blue'
        self.telescope = telescopes.WHTTelescopePar()
        self.camera = 'ISISb'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            dispaxis        = 0,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.225,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 1.2,
                            ronoise         = 5.0,
                            datasec         = '[:,2:4030]',
                            suffix          = '_blue'
                            )]
        self.numhead = 2
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    def check_headers(self, headers):
        """
        Check headers match expectations for an WHT ISIS Blue exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.DETECTOR': 'EEV12',
                              '0.ISIARM': 'Blue arm',
                               '1.NAXIS': 2 }
        super(WHTISISBlueSpectrograph, self).check_headers(headers,
                                                           expected_values=expected_values)

    def validate_metadata(self, fitstbl):
        fitstbl['binning'] = np.array(['{0},{1}'.format(bx,by) 
                                for bx,by in zip(fitstbl['binning_x'], fitstbl['binning_y'])])

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['lamps'] == 'Off') & (fitstbl['idname'] == 'object')
        if ftype == 'bias':
            return good_exp & (fitstbl['idname'] == 'zero')
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['lamps'] == 'W') & (fitstbl['idname'] == 'flat')
        if ftype in ['pinhole', 'dark']:
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & (fitstbl['lamps'] == 'CuNe+CuAr') & (fitstbl['idname'] == 'arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        # Standard
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'
#        match_criteria['standard']['match']['decker'] = ''
#        match_criteria['standard']['match']['dispangle'] = '|<=1'
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'
        # Pixelflat
        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'
        match_criteria['pixelflat']['match']['decker'] = ''
        match_criteria['pixelflat']['match']['dispangle'] = '|<=1'
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['pixelflat']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'
        match_criteria['arc']['match']['dispangle'] = '|<=1'

        return match_criteria

    def default_pypeit_par(self):
        """
        Set default parameters for Keck LRISr reductions.
        """
        par = WHTISISSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'wht_isis_blue'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['min_nsig'] = 10.0
        par['calibrations']['wavelengths']['lowest_nsig'] =10.0
        par['calibrations']['wavelengths']['lamps'] = ['CuI', 'NeI', 'ArI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 2

        return par
