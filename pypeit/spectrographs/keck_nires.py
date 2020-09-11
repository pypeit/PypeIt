""" Module for Keck/NIRES specific codes
"""
import numpy as np

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pkg_resources import resource_filename


class KeckNIRESSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(KeckNIRESSpectrograph, self).__init__()
        self.spectrograph = 'keck_nires'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRES'
        self.numhead = 3

    def get_detector_par(self, hdu, det):
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip=False,
            platescale      = 0.15,
            darkcurr        = 0.01,
            saturation      = 1e6, # I'm not sure we actually saturate with the DITs???
            nonlinear       = 0.76,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(3.8),
            ronoise         = np.atleast_1d(5.0),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = np.atleast_1d('[980:1024,:]')  # Is this a hack??
            )
        detector = detector_container.DetectorContainer(**detector_dict)
        return detector

    @property
    def pypeline(self):
        return 'Echelle'

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'keck_nires'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.20 #0.20  # Might be grating dependent..
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        #par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'keck_nires.fits'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.4
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fwhm_gaussian'] = 4.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0
        #par['calibrations']['tilts']['spat_order'] =  3
        #par['calibrations']['tilts']['spec_order'] =  3

        # Processing steps
        turn_off = dict(use_illumflat=False, use_biasimage=False, use_overscan=False, use_darkimage=False)
        par.reset_all_processimages_par(**turn_off)


        # Extraction
        par['reduce']['skysub']['bspline_spacing'] = 0.8
        par['reduce']['extraction']['sn_gauss'] = 4.0

        # Flexure
        par['flexure']['spec_method'] = 'skip'

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'
        par['reduce']['extraction']['boxcar_radius'] = 0.75  # arcsec


        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [100, None]
        par['calibrations']['tiltframe']['exprng'] = [100, None]
        par['calibrations']['darkframe']['exprng'] = [60, None]
        par['scienceframe']['exprng'] = [60, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

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
        """
        Add additional keys to determine the instrument configuration

        Returns:
            list:

        """
        return ['dispname']

    def pypeit_file_keys(self):
        """
        Add additional columns to the file block of the PypeIt file

        Returns:
            list:

        """
        pypeit_keys = super(KeckNIRESSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['pinhole', 'bias']:
            # No pinhole or bias frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        if ftype == 'dark':
            return good_exp & (fitstbl['idname'] == 'dark')
        if ftype in ['pixelflat', 'trace']:
            return fitstbl['idname'] == 'domeflat'
        if ftype in 'science':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        if ftype in ['arc', 'tilt']:
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        return np.zeros(len(fitstbl), dtype=bool)


    def bpm(self, filename, det, shape=None, msbias=None):
        """
        Override parent bpm function with BPM specific to X-Shooter VIS.

        .. todo::
            Allow for binning changes.

        Parameters
        ----------
        det : int, REQUIRED
        msbias : numpy.ndarray, required if the user wishes to generate a BPM based on a master bias
        **null_kwargs:
            Captured and never used

        Returns
        -------
        bpix : ndarray
          0 = ok; 1 = Mask

        """
        msgs.info("Custom bad pixel mask for NIRES")
        bpm_img = self.empty_bpm(filename, det, shape=shape)

        # Fill in bad pixels if a master bias frame is provided
        if msbias is not None:
            return self.bpm_frombias(msbias, det, bpm_img)

        if det == 1:
            bpm_img[:, :20] = 1.
            bpm_img[:, 1000:] = 1.

        return bpm_img



    @property
    def norders(self):
        return 5

    @property
    def order_spat_pos(self):
        ord_spat_pos = np.array([0.22773035, 0.40613574, 0.56009658,
                                   0.70260714, 0.86335914])
        return ord_spat_pos

    @property
    def orders(self):
        return np.arange(7, 2, -1, dtype=int)


    @property
    def spec_min_max(self):
        spec_max = np.asarray([np.inf]*self.norders)
        spec_min = np.asarray([1024, -np.inf, -np.inf, -np.inf, -np.inf])
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        NIRES has no binning

        Args:
            order_vec (np.ndarray):
            binning (optional):

        Returns:
            np.ndarray:

        """
        norders = order_vec.size
        return np.full(norders, 0.15)




