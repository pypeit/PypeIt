""" Module for TNG/Dolores
"""
from __future__ import absolute_import, division, print_function


import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

class TNGDoloresSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """

    def __init__(self):
        super(TNGDoloresSpectrograph, self).__init__()
        self.spectrograph = 'tng_dolores'
        self.telescope = telescopes.TNGTelescopePar()
        self.camera = 'DOLORES'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            dispaxis        = 1,
                            dispflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.252,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 0.97,
                            ronoise         = 9.0,
                            datasec         = '[51:,1:2045]',
                            oscansec        = '[51:,2054:]',
                            suffix          = '_lrr'
                            )]
        self.numhead = 1
        # Uses default primary_hdrext
        self.timeunit = 'isot'
        # self.sky_file = ?

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for TNG Dolores reductions.
        """
        par = pypeitpar.PypeItPar()
        par['calibrations']['tilts']['params'] = [1,1,1]
        # Always sky subtract, starting with default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 0.1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['scienceframe']['exprng'] = [1, None]
        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an TNG Dolores exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.DET_ID': 'E2V4240',
                             '0.NAXIS': 2}
        super(TNGDoloresSpectrograph, self).check_headers(headers, expected_values=expected_values)

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
        hdr_keys[2] = {}
        hdr_keys[3] = {}
        hdr_keys[4] = {}

        # Copied over defaults
        hdr_keys[0]['idname'] = 'OBS-TYPE'
        hdr_keys[0]['target'] = 'OBJCAT'
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['time'] = 'DATE-OBS'
        hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'
        hdr_keys[0]['filter1'] = 'FLT_ID'
        hdr_keys[0]['dispname'] = 'GRM_ID'
        hdr_keys[0]['lamps'] = 'LMP_ID'

        return hdr_keys

    def metadata_keys(self):
        return ['filename', 'date', 'frametype', 'target', 'exptime', 'dispname', 'idname',
                'configuration', 'calib']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'OBJECT') & (fitstbl['lamps'] == 'Parking') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype == 'bias':
            return good_exp & (fitstbl['dispname'] == 'OPEN')
        if ftype == 'pixelflat' or ftype == 'trace':
            return good_exp & (fitstbl['idname'] == 'CALIB') & (fitstbl['lamps'] == 'Halogen') \
                        & (fitstbl['dispname'] != 'OPEN')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'arc') & (fitstbl['lamps'] == 'Ne+Hg') \
                        & (fitstbl['dispname'] != 'OPEN')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_match_criteria(self):
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}
        #
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'
        # Bias
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'
        # Pixelflat
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'
        match_criteria['pixelflat']['match']['dispname'] = ''
        # Traceflat
        match_criteria['trace']['match'] = match_criteria['pixelflat']['match'].copy()
        # Arc
        match_criteria['arc']['match'] = match_criteria['pixelflat']['match'].copy()

        # Return
        return match_criteria

    def setup_arcparam(self, arcparam, disperser=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place
            modify_dict: dict

        """
        arcparam['lamps'] = ['NeI', 'HgI']
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']

        if disperser == 'LR-R':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.61  # Ang per pixel (unbinned)
            arcparam['disp_toler'] = 0.1  # Ang per pixel (unbinned)
            arcparam['wvmnx'][0] = 4470.0
            arcparam['wvmnx'][1] = 10073.0
            arcparam['wv_cen'] = 7400.
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0] / binspectral
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

