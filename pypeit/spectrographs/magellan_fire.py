""" Module for Magellan/FIRE specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels
from pypeit import debugger

class MagellanFIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Magellan/FIRE specific code

    .. note::
        For FIRE Echelle, we usually use high gain and SUTR read mode.
        The exposure time is usually around 900s. The detector
        parameters below are based on such mode. Standard star and
        calibrations are usually use Fowler 1 read mode in which case
        the read noise is ~20 electron.

    """
    def __init__(self):
        # Get it started
        super(MagellanFIRESpectrograph, self).__init__()
        self.spectrograph = 'magellan_fire'
        self.telescope = telescopes.MagellanTelescopePar()
        self.camera = 'FIRE'
        self.numhead = 1
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.01,
                            saturation      = 20000., # high gain mode, low gain is 32000
                            nonlinear       = 1.0, # high gain mode, low gain is 0.875
                            numamplifiers   = 1,
                            gain            = 1.2, # high gain mode, low gain is 3.8 e-/DN
                            ronoise         = 5.0, # for high gain mode and SUTR read modes with exptime ~ 900s
                            datasec         = '[1:2048,1:2048]',
                            oscansec        = '[:,:4]'
                            )]
        self.norders = 22
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @property
    def pypeline(self):
        return 'Echelle'

    # TODO: Remove dependency on self.  Non-linear counts does not need
    # to be a parameter.
    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_fire'
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
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['lamps'] = ['OH_XSHOOTER']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']

        # Reidentification parameters
        #par['calibrations']['wavelengths']['method'] = 'reidentify'
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire.json'
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

        # Set slits and tilts parameters
#        par['calibrations']['tilts']['order'] = 2
        par['calibrations']['tilts']['tracethresh'] = [10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 10]
        par['calibrations']['slits']['trace_npoly'] = 5
        par['calibrations']['slits']['sigdetect'] = 50
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['pcatype'] = 'pixel'
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Do not correct for flexure
        par['flexure'] = pypeitpar.FlexurePar()
        par['flexure']['method'] = 'skip'
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

#    def check_headers(self, headers):
#        """
#        Check headers match expectations for a Keck NIRES exposure.
#
#        See also
#        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.
#
#        Args:
#            headers (list):
#                A list of headers read from a fits file
#        """
#        expected_values = { '0.INSTRUME': 'FIRE',
#                               '0.NAXIS': 2,
#                              '0.NAXIS1': 2048,
#                              '0.NAXIS2': 2048 }
#        super(MagellanFIRESpectrograph, self).check_headers(headers, expected_values=expected_values)

#    def header_keys(self):
#        """
#        Return a dictionary with the header keywords to read from the
#        fits file.
#
#        Returns:
#            dict: A nested dictionary with the header keywords to read.
#            The first level gives the extension to read and the second
#            level gives the common name for header values that is passed
#            on to the PypeItMetaData object.
#        """
#        hdr_keys = {}
#        hdr_keys[0] = {}
#
#        # Copied over defaults
#        hdr_keys[0]['idname'] = 'OBSTYPE'
#        #hdr_keys[0]['time'] = 'MJD-OBS'
#        hdr_keys[0]['date'] = 'DATE-OBS'
#        hdr_keys[0]['utc'] = 'UT-TIME'
#        hdr_keys[0]['ra'] = 'RA'
#        hdr_keys[0]['dec'] = 'DEC'
#        hdr_keys[0]['airmass'] = 'AIRMASS'
#        hdr_keys[0]['exptime'] = 'EXPTIME'
#        hdr_keys[0]['target'] = 'OBJECT'
#        hdr_keys[0]['naxis0'] = 'NAXIS2'
#        hdr_keys[0]['naxis1'] = 'NAXIS1'
#        hdr_keys[0]['binning'] = 1
#        hdr_keys[0]['dispname'] = 'INSTR'  # Should be 'spec' if in the spectroscopy mode
#
#        return hdr_keys
#
#    def metadata_keys(self):
#        return ['filename', 'date', 'frametype', 'idname','target', 'exptime', 'setup', 'calib',
#                'obj_id', 'bkg_id']

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        self.meta = {}
        # Required (core)
        self.meta['ra'] = dict(ext=0, card='RA')
        self.meta['dec'] = dict(ext=0, card='DEC')
        self.meta['target'] = dict(ext=0, card='OBJECT')
        #TODO: Check decker is correct
        self.meta['decker'] = dict(ext=0, card='SLIT')
        self.meta['binning'] = dict(ext=0, card=None, default='1,1')
        self.meta['mjd'] = dict(ext=0, card=None, compound=True)
        self.meta['exptime'] = dict(ext=0, card='EXPTIME')
        self.meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        self.meta['dispname'] = dict(ext=0, card='INSTR')

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'mjd':
            time = headarr[0]['DATE']
            ttime = Time(time, format='isot')
            return ttime.mjd
        msgs.error("Not ready for this compound meta")

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['idname'] == 'PixFlat')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'Telluric')
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'Science')
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'Science')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

#    def parse_binning(self, inp, det=1):
#        return '1,1'

#    def get_match_criteria(self):
#        """Set the general matching criteria for FIRE"""
#        match_criteria = {}
#        for key in framematch.FrameTypeBitMask().keys():
#            match_criteria[key] = {}
#
#        match_criteria['standard']['match'] = {}
#        match_criteria['standard']['match']['naxis0'] = '=0'
#        match_criteria['standard']['match']['naxis1'] = '=0'
#
#        match_criteria['bias']['match'] = {}
#        match_criteria['bias']['match']['naxis0'] = '=0'
#        match_criteria['bias']['match']['naxis1'] = '=0'
#
#        match_criteria['pixelflat']['match'] = {}
#        match_criteria['pixelflat']['match']['naxis0'] = '=0'
#        match_criteria['pixelflat']['match']['naxis1'] = '=0'
#
#        match_criteria['trace']['match'] = {}
#        match_criteria['trace']['match']['naxis0'] = '=0'
#        match_criteria['trace']['match']['naxis1'] = '=0'
#
#        match_criteria['arc']['match'] = {}
#        match_criteria['arc']['match']['naxis0'] = '=0'
#        match_criteria['arc']['match']['naxis1'] = '=0'
#
#        return match_criteria

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
        msgs.info("Custom bad pixel mask for FIRE")
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[:, :4] = 1.

        return self.bpm_img

    @staticmethod
    def slitmask(tslits_dict, pad=None, binning=None):
        """
         Generic routine ton construct a slitmask image from a tslits_dict. Children of this class can
         overload this function to implement instrument specific slitmask behavior, for example setting
         where the orders on an echelle spectrograph end

         Parameters
         -----------
         tslits_dict: dict
            Trace slits dictionary with slit boundary information

         Optional Parameters
         pad: int or float
            Padding of the slit boundaries
         binning: tuple
            Spectrograph binning in spectral and spatial directions

         Returns
         -------
         slitmask: ndarray int
            Image with -1 where there are no slits/orders, and an integer where there are slits/order with the integer
            indicating the slit number going from 0 to nslit-1 from left to right.

         """

        # These lines are always the same
        pad = tslits_dict['pad'] if pad is None else pad
        slitmask = pixels.slit_pixels(tslits_dict['lcen'], tslits_dict['rcen'], tslits_dict['nspat'], pad=pad)

        spec_img = np.outer(np.arange(tslits_dict['nspec'], dtype=int), np.ones(tslits_dict['nspat'], dtype=int))  # spectral position everywhere along image

        order7bad = (slitmask == 0) & (spec_img < tslits_dict['nspec']/2)
        slitmask[order7bad] = -1
        return slitmask

    @staticmethod
    def slit2order(islit):

        """
        Parameters
        ----------
        islit: int, float, or string, slit number

        Returns
        -------
        order: int
        """

        if isinstance(islit, str):
            islit = int(islit)
        elif isinstance(islit, np.ndarray):
            islit = islit.astype(int)
        elif isinstance(islit, float):
            islit = int(islit)
        elif isinstance(islit, int):
            pass
        else:
            msgs.error('Unrecognized type for islit')

        orders = np.arange(32, 11, -1, dtype=int)
        return orders[islit]

    @staticmethod
    def order_platescale(binning = None):


        """
        Returns the plate scale in arcseconds for each order

        Parameters
        ----------
        None

        Optional Parameters
        --------------------
        binning: str

        Returns
        -------
        order_platescale: ndarray, float

        """

        # FIRE has no binning, but for an instrument with binning we would do this
        #binspatial, binspectral = parse.parse_binning(binning)
        return np.full(5, 0.15)

