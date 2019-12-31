""" Module for Gemini/GNIRS specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.core.wavecal import wvutils
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
        # No overscan
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'

        # Slits
        par['calibrations']['slitedges']['edge_thresh'] = 20.
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['pca_order'] = 3

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
        par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs.fits'
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

        # Finding objects
        par['scienceimage']['skysub']['bspline_spacing'] = 0.8
        par['scienceimage']['findobj']['sig_thresh'] = 5.0
        par['scienceimage']['findobj']['find_trim_edge'] = [2,2]    # Slit is too short to trim 5,5 especially
        par['scienceimage']['findobj']['find_cont_fit'] = False     # Don't continuum fit objfind for narrow slits
        par['scienceimage']['findobj']['find_npoly_cont'] = 0       # Continnum order for determining thresholds
        # Extraction
        par['scienceimage']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        # Sky Subtraction
        par['scienceimage']['skysub']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars
        par['scienceimage']['skysub']['no_poly'] = True         # Do not use polynomial degree of freedom for global skysub


        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['pixelflatframe']['exprng'] = [None, 30]
        par['calibrations']['traceframe']['exprng'] = [None, 30]
        par['calibrations']['standardframe']['exprng'] = [None, 30]
        par['scienceframe']['exprng'] = [30, None]

        # Do not bias subtract
        #par['scienceframe']['useframe'] = 'overscan'
        # This is a hack for now until we can specify for each image type what to do. Bias currently
        # controls everything
        par['calibrations']['biasframe']['useframe'] = 'none'



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
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & (fitstbl['idname'] == 'FLAT')
        if ftype in ['pinhole', 'dark', 'bias']:
            # Don't type pinhole, dark, or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['idname'] == 'ARC')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

#    def parse_binning(self, inp, det=1):
#        return '1,1'

    def order_platescale(self, order_vec, binning=None):


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
        return np.full(order_vec.size, 0.15)


    @property
    def norders(self):
        return 6

    @property
    def order_spat_pos(self):
        ord_spat_pos = np.array([0.2955097 , 0.37635756, 0.44952223, 0.51935601, 0.59489503, 0.70210309])
        return ord_spat_pos

    @property
    def orders(self):
        return np.arange(8,2,-1, dtype=int)


    @property
    def spec_min_max(self):
        spec_max = np.asarray([1022,1022,1022,1022,1022,1022])
        spec_min = np.asarray([512,280, 0, 0, 0, 0])
        return np.vstack((spec_min, spec_max))

    @property
    def dloglam(self):
        dloglam = 0.000127888 # this is the average of the median dispersions
        return dloglam

    @property
    def loglam_minmax(self):
        return np.log10(7000), np.log10(26000)

    def wavegrid(self, binning=None, samp_fact=1.0, midpoint=False):

        # Define the grid for GNIRS
        logmin, logmax = self.loglam_minmax
        loglam_grid = wvutils.wavegrid(logmin, logmax, self.dloglam, samp_fact=samp_fact)
        if midpoint:
            loglam_grid = loglam_grid + self.dloglam/2.0
        return np.power(10.0,loglam_grid)


    def bpm(self, filename, det, shape=None):
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
        bpm_img = self.empty_bpm(filename, det, shape=shape)
        if det == 1:
            bpm_img[:, :20] = 1.
            bpm_img[:, 1000:] = 1.

        return bpm_img



