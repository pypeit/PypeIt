""" Module for Keck/NIRES specific codes
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
                            specaxis        = 1,
                            specflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.01,
                            saturation      = 1e6, # I'm not sure we actually saturate with the DITs???
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 3.8,
                            ronoise         = 5.0,
                            datasec         = '[1:2048,1:1024]',
                            oscansec        = '[1:2048,980:1024]'
                            )]
        self.norders = 5
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
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.json'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        #par['calibrations']['tilts']['spat_order'] =  3
        #par['calibrations']['tilts']['spec_order'] =  3

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = False

        # Extraction
        par['scienceimage']['bspline_spacing'] = 0.8
        par['scienceimage']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Do not bias subtract
        par['scienceframe']['useframe'] ='overscan'
        # This is a hack for now until we can specify for each image type what to do. Bias currently
        # controls everything
        par['calibrations']['biasframe']['useframe'] = 'overscan'




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
        meta['decker'] = dict(ext=0, card=None, default='default')
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='ITIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='INSTR')
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta

    def configuration_keys(self):
        return ['dispname']

    def pypeit_file_keys(self):
        pypeit_keys = super(KeckNIRESSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

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

#    def parse_binning(self, inp, det=1):
#        return '1,1'

#    def get_match_criteria(self):
#        """Set the general matching criteria for NIRES"""
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
#        # OLD
#        # Bias
#        #match_criteria['bias']['match'] = {}
#        #match_criteria['standard']['match'] = {}
#        #match_criteria['pixelflat']['match'] = {}
#        #match_criteria['trace']['match'] = {}
#        #match_criteria['arc']['match'] = {}
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
        msgs.info("Custom bad pixel mask for NIRES")
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[:, :20] = 1.
            self.bpm_img[:, 1000:] = 1.

        return self.bpm_img

    def slit_minmax(self, nslits, binspectral=1):

        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        spec_max = np.asarray([np.inf, np.inf, np.inf, np.inf, np.inf, np.inf])
        spec_min = np.asarray([1024, -np.inf, -np.inf, -np.inf, -np.inf, -np.inf])

        return spec_min, spec_max

    def slitmask(self, tslits_dict, pad=None):
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

        orders = np.arange(7, 2, -1, dtype=int)
        return orders[islit]

    def order_vec(self):
        return self.slit2order(np.arange(self.norders))


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

        # NIRES has no binning, but for an instrument with binning we would do this
        #binspatial, binspectral = parse.parse_binning(binning)
        return np.full(5, 0.15)

    # JFH This function should probably take a disperser name or something for other instruments??
    def wavegrid(self, binning=None):
        """
        For fixed format echelle spectrographs return the wavelength grid used for 2D coadds
        Args:
            binning:

        Returns:

        """

        # Define the grid for NIRES
        R = 2700.0 * 2.7
        dloglam = 1.0 / R / np.log(10.0)
        logmin = np.log10(9500.0)
        logmax = np.log10(26000)
        ngrid = int(np.ceil((logmax - logmin) / dloglam))
        osamp = 1.0
        loglam_grid = logmin + (dloglam / osamp) * np.arange(int(np.ceil(osamp * ngrid)))

        return np.power(10.0,loglam_grid)






