""" Module for Gemini/GNIRS specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels


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
        self.numhead = 2
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 1,
                            specaxis        = 0,
                            specflip=True,
                            spatflip=True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.15,
                            saturation      = 90000.,
                            nonlinear       = 0.71,
                            numamplifiers   = 1,
                            gain            = 13.5,
                            ronoise         = 7.0,
                            datasec         = '[:,:]',#'[1:1024,1:1022]',
                            oscansec        = '[:,:]',#'[1:1024,1:1022]'
                            )]
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?
    @property
    def pypeline(self):
        return 'Echelle'

    @property
    def norders(self):
        return 6

    def default_pypeit_par(self):
        """
        Set default parameters for Gemini GNIRS reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'gemini_gnirs'
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1

        # Slits
        par['calibrations']['slits']['sigdetect'] = 50.
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 0.5

        # Wavelengths
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['lamps'] = ['OH_GNIRS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 2
        par['calibrations']['wavelengths']['n_final'] = [1,3,3,3,3,3]

        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.6
        par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs.json'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        # JFH This is provisional these IDs should be checked.
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Tilts
        par['calibrations']['tilts']['tracethresh'] = [5.0,10,10,10,10,10]
        par['calibrations']['tilts']['sig_neigh'] = 5.0
        par['calibrations']['tilts']['nfwhm_neigh'] = 2.0

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = False
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Extraction
        par['scienceimage']['bspline_spacing'] = 0.8
        par['scienceimage']['sn_gauss'] = 4.0

        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['calibrations']['standardframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [30, None]

        # Do not bias subtract
        par['scienceframe']['useframe'] ='overscan'
        # This is a hack for now until we can specify for each image type what to do. Bias currently
        # controls everything
        par['calibrations']['biasframe']['useframe'] = 'overscan'



        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an LRISb exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'GNIRS',
                               '1.NAXIS': 2 }
        super(GeminiGNIRSSpectrograph, self).check_headers(headers,
                                                           expected_values=expected_values)

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
        hdr_keys[0]['idname'] = 'OBSTYPE'
#        hdr_keys[0]['date'] = 'DATE'
        hdr_keys[0]['ut'] = 'UT'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['time'] = 'MJD_OBS'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['slit'] = 'SLIT'
        hdr_keys[0]['decker'] = 'DECKER'
        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['hatch'] = 'COVER'
        hdr_keys[0]['dispname'] = 'GRATING'
        hdr_keys[0]['dispangle'] = 'GRATTILT'
        hdr_keys[0]['wavecen'] = 'GRATWAVE'
        hdr_keys[0]['spectrograph'] = 'INSTRUME'
        hdr_keys[0]['binning'] = ' '

        return hdr_keys

    def configuration_keys(self):
        return ['decker', 'dispname']

    def metadata_keys(self):
        return super(GeminiGNIRSSpectrograph, self).metadata_keys() + ['dispangle', 'idname']

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype == 'science':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'standard':
            return good_exp & (fitstbl['idname'] == 'OBJECT')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype == 'pinhole' or ftype == 'dark' or ftype == 'bias':
            # Don't type pinhole, dark, or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & (fitstbl['idname'] == 'ARC')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)




    def order_platescale(self, binning=None):


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

        # Right now I just assume a simple linear trend
        return np.full(self.norders, 0.15)



    def slit2order(self, islit):

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
        elif isinstance(islit, (int,np.int64,np.int32,np.int)):
            pass
        else:
            msgs.error('Unrecognized type for islit')

        orders = np.arange(8,2,-1, dtype=int)
        return orders[islit]



    def slitmask(self, tslits_dict, pad=None, binning=None):
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

        nslits = tslits_dict['lcen'].shape[1]
        if nslits != self.norders:
            msgs.error('There is a problem with your slit bounadries. You have nslits={:d} orders, whereas GNIRS has norders={:d}'.format(nslits,self.norders))
        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        order_max = [1022,1022,1022,1022,1022,1022]
        order_min = [512,280, 0, 0, 0, 0]
        # TODO add binning adjustments to these
        for islit in range(nslits):
            orderbad = (slitmask == islit) & ((spec_img < order_min[islit]) | (spec_img > order_max[islit]))
            slitmask[orderbad] = -1
        return slitmask


    def get_match_criteria(self):

        """Set the general matching criteria for GNIRS. Copied from NIRES"""
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
#        match_criteria['standard']['match']['naxis0'] = '=0'
#        match_criteria['standard']['match']['naxis1'] = '=0'

        match_criteria['bias']['match'] = {}
#        match_criteria['bias']['match']['naxis0'] = '=0'
#        match_criteria['bias']['match']['naxis1'] = '=0'

        match_criteria['pixelflat']['match'] = {}
#        match_criteria['pixelflat']['match']['naxis0'] = '=0'
#        match_criteria['pixelflat']['match']['naxis1'] = '=0'

        match_criteria['trace']['match'] = {}
#        match_criteria['trace']['match']['naxis0'] = '=0'
#        match_criteria['trace']['match']['naxis1'] = '=0'

        match_criteria['arc']['match'] = {}
#        match_criteria['arc']['match']['naxis0'] = '=0'
#        match_criteria['arc']['match']['naxis1'] = '=0'

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



