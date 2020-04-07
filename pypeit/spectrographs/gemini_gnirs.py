""" Module for Gemini/GNIRS specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.images import detector_container
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pkg_resources import resource_filename
from IPython import embed


from pypeit import debugger

class GeminiGNIRSSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Gemini/GNIRS specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(GeminiGNIRSSpectrograph, self).__init__()
        self.spectrograph = 'gemini_gnirs'
        self.telescope = telescopes.GeminiNTelescopePar()
        self.camera = 'GNIRS'
        self.dispname = None # TODO We need a model for setting setup specific parameters in spectrograph
        self.numhead = 2

    def get_detector_par(self, hdu, det):
        """
        Return a DectectorContainer for the current image

        Args:
            hdu (`astropy.io.fits.HDUList`):
                HDUList of the image of interest.
                Ought to be the raw file, or else..
            det (int):

        Returns:
            :class:`pypeit.images.detector_container.DetectorContainer`:

        """
        # Detector 1
        detector_dict = dict(
            binning         = '1,1',
            det             = 1,
            dataext         = 1,
            specaxis        = 0,
            specflip=True,
            spatflip=True,
            platescale      = 0.15,
            darkcurr        = 0.15,
            saturation      = 150000.,
            nonlinear       = 0.71,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(13.5),
            ronoise         = np.atleast_1d(7.0),
            datasec         = np.atleast_1d('[:,:]'),#'[1:1024,1:1022]',
            oscansec        = np.atleast_1d('[:,:]'),#'[1:1024,1:1022]'
        )

        return detector_container.DetectorContainer(**detector_dict)

    @property
    def pypeline(self):
        return 'Echelle'


    def default_pypeit_par(self):
        """
        Set default parameters for Gemini GNIRS reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'gemini_gnirs'
        # No overscan
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'

        # Flats
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10
        par['calibrations']['standardframe']['process']['illumflatten'] = False
        par['scienceframe']['process']['illumflatten'] = False

        # Finding objects
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['findobj']['sig_thresh'] = 5.0
        par['reduce']['findobj']['find_trim_edge'] = [2,2]    # Slit is too short to trim 5,5 especially
        par['reduce']['findobj']['find_cont_fit'] = False     # Don't continuum fit objfind for narrow slits
        par['reduce']['findobj']['find_npoly_cont'] = 0       # Continnum order for determining thresholds
        # Extraction
        par['reduce']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        # Sky Subtraction
        par['reduce']['skysub']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars
        par['reduce']['skysub']['no_poly'] = True         # Do not use polynomial degree of freedom for global skysub


        # Do not correct for flexure
        par['flexure']['spec_method'] = 'skip'

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

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

        return par

    def config_specific_par(self, scifile, inp_par=None):
        """
        Modify the PypeIt parameters to hard-wired values used for specific instrument configurations.


        Args:
           scifile (str):
                File to use when determining the configuration and how
                to adjust the input parameters.
           inp_par (:class:`pypeit.par.parset.ParSet`, optional):
                Parameter set used for the full run of PypeIt.  If None,
                use :func:`default_pypeit_par`.

        Returns:
           :class:`pypeit.par.parset.ParSet`: The PypeIt paramter set  adjusted for configuration specific parameter values.
        """

        par = self.default_pypeit_par() if inp_par is None else inp_par

        # TODO This is a hack for now until we figure out how to set dispname and other meta information in the
        # spectrograph class itself
        self.dispname = self.get_meta_value(scifile, 'dispname')
        # 32/mmSB_G5533 setup, covering XYJHK with short blue camera
        if '32/mm' in self.dispname:
            # Edges
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
            #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = [1, 3, 3, 3, 3, 3]

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
            par['calibrations']['tilts']['tracethresh'] = [5.0, 10, 10, 10, 10, 10]
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0
        # 10/mmLBSX_G5532 setup, covering YJHK with the long blue camera and SXD prism
        elif '10/mmLBSX' in self.dispname:
            # Edges
            par['calibrations']['slitedges']['edge_thresh'] = 20.
            par['calibrations']['slitedges']['trace_thresh'] = 10.
            par['calibrations']['slitedges']['fit_order'] = 2
            par['calibrations']['slitedges']['max_shift_adj'] = 0.5
            par['calibrations']['slitedges']['det_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['fit_min_spec_length'] = 0.20
            par['calibrations']['slitedges']['left_right_pca'] = True # Actually we need a parameter to disable PCA entirely
            par['calibrations']['slitedges']['sync_predict'] = 'nearest'

            # Wavelengths
            par['calibrations']['wavelengths']['rms_threshold'] = 1.0  # Might be grating dependent..
            par['calibrations']['wavelengths']['sigdetect'] = 5.0
            par['calibrations']['wavelengths']['lamps'] = ['Ar_IR_GNIRS']
            #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
            par['calibrations']['wavelengths']['n_first'] = 2
            par['calibrations']['wavelengths']['n_final'] = [3, 3, 3, 3]
            # Reidentification parameters
            par['calibrations']['wavelengths']['method'] = 'reidentify'
            par['calibrations']['wavelengths']['cc_thresh'] = 0.6
            par['calibrations']['wavelengths']['reid_arxiv'] = 'gemini_gnirs_10mm_LBSX.fits'
            par['calibrations']['wavelengths']['ech_fix_format'] = True
            # Echelle parameters
            par['calibrations']['wavelengths']['echelle'] = True
            par['calibrations']['wavelengths']['ech_nspec_coeff'] = 3
            par['calibrations']['wavelengths']['ech_norder_coeff'] = 3
            par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

            # Tilts
            par['calibrations']['tilts']['tracethresh'] = [10, 10, 10, 10]
            par['calibrations']['tilts']['sig_neigh'] = 5.0
            par['calibrations']['tilts']['nfwhm_neigh'] = 2.0
        else:
            msgs.erro('Unrecognized GNIRS dispname')

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
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.full(order_vec.size, 0.05)
        elif '32/mm' in self.dispname:
            return np.full(order_vec.size, 0.15)
        else:
            msgs.error('Unrecognized disperser')

    @property
    def norders(self):
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return 4
        elif '32/mm' in self.dispname:
            return 6
        else:
            msgs.error('Unrecognized disperser')

    @property
    def order_spat_pos(self):
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.array([0.050, 0.215, 0.442, 0.759])
        elif '32/mm' in self.dispname:
            return np.array([0.2955097 , 0.37635756, 0.44952223, 0.51935601, 0.59489503, 0.70210309])
        else:
            msgs.error('Unrecognized disperser')

    @property
    def orders(self):
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            return np.arange(6,2,-1, dtype=int)
        elif '32/mm' in self.dispname:
            return np.arange(8,2,-1,dtype=int)
        else:
            msgs.error('Unrecognized disperser')

    @property
    def spec_min_max(self):
        self.check_disperser()
        if '10/mmLBSX' in self.dispname:
            spec_max = np.asarray([1022, 1022, 1022, 1022])
            spec_min = np.asarray([450, 0, 0, 0])
            return np.vstack((spec_min, spec_max))
        elif '32/mm' in self.dispname:
            spec_max = np.asarray([1022, 1022, 1022, 1022, 1022, 1022])
            spec_min = np.asarray([512, 280, 0, 0, 0, 0])
            return np.vstack((spec_min, spec_max))
        else:
            msgs.error('Unrecognized disperser')



    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to Gemini/GNIRS

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

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        # JFH Changed. Dealing with detector scratch
        if det == 1:
            bpm_img[687:765,12:16] = 1.
            bpm_img[671:687,8:13] = 1.
        #    bpm_img[:, 1000:] = 1.

        return bpm_img



