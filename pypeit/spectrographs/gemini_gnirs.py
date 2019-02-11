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
                            saturation      = 150000.,
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
        par['calibrations']['slits']['trace_npoly'] = 5
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
        par['scienceimage']['sig_thresh'] = 5.0
        par['scienceimage']['bspline_spacing'] = 0.8
        par['scienceimage']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['scienceimage']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars
        par['scienceimage']['no_poly']  = True         # Do not use polynomial degree of freedom for global skysub


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

    def init_meta(self):
        """
        Generate the meta data dict
        Note that the children can add to this

        Returns:
            self.meta: dict (generated in place)

        """
        meta = {}
        # Required (core)
        meta['ra'] = dict(ext=0, card='RA')
        meta['dec'] = dict(ext=0, card='DEC')
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['decker'] = dict(ext=0, card='DECKER')

        meta['binning'] = dict(ext=0, card=None, default='1,1')
        meta['mjd'] = dict(ext=0, card='MJD_OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='GRATING')
        meta['hatch'] = dict(ext=0, card='COVER')
        meta['dispangle'] = dict(ext=0, card='GRATTILT', rtol=1e-4)
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta

    def configuration_keys(self):
        return ['decker', 'dispname', 'dispangle']

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

#    def parse_binning(self, inp, det=1):
#        return '1,1'

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


    def order_vec(self):
        return self.slit2order(np.arange(self.norders))


    def slit_minmax(self, nslits, binspectral=1):

        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        spec_max = np.asarray([1022,1022,1022,1022,1022,1022])
        spec_min = np.asarray([512,280, 0, 0, 0, 0])

        return spec_min, spec_max

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


    def wavegrid(self, binning=None):

        # Define the new wavelength grid for GNIRS
        ngrid = 5000
        dloglam = 0.000127888 # this is the average of the median dispersions
        logmin = 3.777
        osamp = 1.0
        loglam_grid = logmin + (dloglam / osamp) * np.arange(int(np.ceil(osamp * ngrid)))

        return np.power(10.0,loglam_grid)


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



