""" Module for LBT/MODS specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np


from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

from pypeit import debugger

# ToDo: test MODS1B and MODS2B

class LBTMODSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODSSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods'
        self.telescope = telescopes.LBTTelescopePar()
        self.timeunit = 'isot'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast reductions.
        """
        par = pypeitpar.PypeItPar()
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 5
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1


        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, 60]
        par['calibrations']['standardframe']['exprng'] = [1, 200]
        par['scienceframe']['exprng'] = [200, None]
        return par

    def header_keys(self):
        """
        Provide the relevant header keywords
        """

        hdr_keys = {}
        hdr_keys[0] = {}

        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['idname'] = 'IMAGETYP'
        hdr_keys[0]['time'] = 'MJD-OBS'
        hdr_keys[0]['utc'] = 'UTC-OBS'
        hdr_keys[0]['date'] = 'DATE-OBS'
        hdr_keys[0]['ra'] = 'OBJRA'
        hdr_keys[0]['dec'] = 'OBJDEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['binning'] = 'CCDXBIN'
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['decker'] = 'MASKNAME'
        hdr_keys[0]['dichroic'] = 'FILTNAME'
        hdr_keys[0]['dispname'] = 'GRATNAME'
        hdr_keys[0]['spectrograph'] = 'INSTRUME'

        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'


        return hdr_keys

    # Uses parent metadata keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'bias':
            return good_exp  & (fitstbl['idname'] == 'BIAS')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp  & (fitstbl['idname'] == 'FLAT')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & ((fitstbl['idname'] == 'COMP') | (fitstbl['idname'] == 'OBJECT'))

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    def get_match_criteria(self):
        """Set the general matching criteria for LBT MODS."""
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
        match_criteria['pixelflat']['match']['decker'] = ''

        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'
        match_criteria['trace']['match']['decker'] = ''

        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'

        return match_criteria


class LBTMODS1RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS1RSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods1r'
        self.camera = 'MODS1R'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.4,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.5,
                            ronoise         = 4.2,
                            datasec='[1:8288,1:3088]',
                            oscansec='[8250:,:]', # ToDo: fix this
                            suffix          = '_mods1r'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS1R reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods1r'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1
        par['calibrations']['wavelengths']['n_final'] = 4

        # slit
        par['calibrations']['slits']['sigdetect'] = 300

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 5
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a LBT MODS1R exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.INSTRUME': 'MODS1R' }
        super(LBTMODS1RSpectrograph, self).check_headers(headers,
                                                             expected_values=expected_values)

class LBTMODS1BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS1BSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods1b'
        self.camera = 'MODS1B'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.120,
                            darkcurr        = 0.5,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.1,
                            ronoise         = 3.0,
                            datasec='[1:8288,1:3088]',
                            oscansec='[8250:,:]',
                            suffix          = '_mods1b'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS1B reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods1b'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1

        # slit
        par['calibrations']['slits']['sigdetect'] = 300

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a LBT MODS1B exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.INSTRUME': 'MODS1B' }
        super(LBTMODS1BSpectrograph, self).check_headers(headers,
                                                             expected_values=expected_values)


class LBTMODS2RSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS2RSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods2r'
        self.camera = 'MODS2R'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.4,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 1.7,
                            ronoise         = 2.8,
                            datasec='[1:8288,1:3088]',
                            oscansec='[8250:,:]',
                            suffix          = '_mods2r'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS2R reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods2r'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['fwhm'] = 10.
        #par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']]
        par['calibrations']['wavelengths']['lamps'] = ['OH_MODS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1
        par['calibrations']['wavelengths']['n_final'] = 4


        # slit
        par['calibrations']['slits']['sigdetect'] = 300

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a LBT MODS2R exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.INSTRUME': 'MODS2R' }
        super(LBTMODS2RSpectrograph, self).check_headers(headers,
                                                             expected_values=expected_values)

class LBTMODS2BSpectrograph(LBTMODSSpectrograph):
    """
    Child to handle LBT/MODS1R specific code
    """
    def __init__(self):
        # Get it started
        super(LBTMODS2RSpectrograph, self).__init__()
        self.spectrograph = 'lbt_mods2b'
        self.camera = 'MODS2B'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.120,
                            darkcurr        = 0.5,
                            saturation      = 65535.,
                            nonlinear       = 0.99,
                            numamplifiers   = 1,
                            gain            = 2.0,
                            ronoise         = 3.7,
                            datasec='[1:8288,1:3088]',
                            oscansec='[8250:,:]',
                            suffix          = '_mods2b'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for LBT MODS2B reductions.
        """
        par = LBTMODSSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'lbt_mods2b'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['sigdetect'] = 5.
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20
        par['calibrations']['wavelengths']['lamps'] = ['XeI','ArII','ArI','NeI','KrI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1

        # slit
        par['calibrations']['slits']['sigdetect'] = 300

        # Set wave tilts order
        par['calibrations']['tilts']['spat_order'] = 3
        par['calibrations']['tilts']['spec_order'] = 5
        par['calibrations']['tilts']['maxdev_tracefit'] = 0.02
        par['calibrations']['tilts']['maxdev2d'] = 0.02

        # reidentification stuff
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'shane_kast_blue_600_4310_d55.json'

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a LBT MODS2B exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.INSTRUME': 'MODS2B' }
        super(LBTMODS2BSpectrograph, self).check_headers(headers,
                                                             expected_values=expected_values)
