""" Module for Magellan/FIRE specific codes
"""
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
#        self.norders = 22
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @property
    def norders(self):
        return 22

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
        par['calibrations']['tilts']['tracethresh'] = [10, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 10]
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['edge_thresh'] = 50
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
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
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'Science')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

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

    @property
    def norders(self):
        return 22

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

