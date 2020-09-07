""" Module for P200/Triplespec specific codes
"""
import numpy as np
from astropy.time import Time

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.images import detector_container

from pkg_resources import resource_filename


class P200TSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle P200/TripleSpec specific code
    """
    ndet = 1

    def __init__(self):
        # Get it started
        super(P200TSPECSpectrograph, self).__init__()
        self.spectrograph = 'p200_tspec'
        self.telescope = telescopes.P200TelescopePar()
        self.camera = 'TSPEC'
        self.numhead = 1

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

        meta['mjd'] = dict(ext=0, card=None, compound=True)
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='FPA')
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'mjd':
            time = headarr[0]['UTSHUT']
            ttime = Time(time, format='isot')
            return ttime.mjd
        else:
            msgs.error("Not ready for this compound meta")

    def get_detector_par(self, hdu, det):
        # Detector 1
        detector_dict = dict(
            binning='1,1',
            det=1,
            dataext         = 0,
            specaxis        = 1,
            specflip        = True,
            spatflip=False,
            platescale      = 0.37,
            darkcurr        = 0.085,
            saturation      = 28000,
            nonlinear       = 0.9,
            mincounts       = -1e10,
            numamplifiers   = 1,
            gain            = np.atleast_1d(3.8),
            ronoise         = np.atleast_1d(3.5),
            datasec         = np.atleast_1d('[:,:]'),
            oscansec        = np.atleast_1d('[:,:]')
            )
        detector = detector_container.DetectorContainer(**detector_dict)
        return detector

    @property
    def pypeline(self):
        return 'Echelle'

    def default_pypeit_par(self):
        """
        Set default parameters for P200 TripleSpec reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'p200_tspec'
        # Wavelengths
        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.3
        par['calibrations']['wavelengths']['sigdetect']=5.0
        par['calibrations']['wavelengths']['fwhm']= 5.0
        par['calibrations']['wavelengths']['n_final']= [3,4,4,4,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_NIRES']
        par['calibrations']['wavelengths']['method'] = 'reidentify'

        # Reidentification parameters
        par['calibrations']['wavelengths']['reid_arxiv'] = 'p200_triplespec.fits'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        par['calibrations']['slitedges']['trace_thresh'] = 5.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.3
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['fwhm_gaussian'] = 4.0

        # Tilt parameters
        par['calibrations']['tilts']['tracethresh'] =  10.0

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
        par['calibrations']['darkframe']['exprng'] = [0, None]
        par['scienceframe']['exprng'] = [60, None]

        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')

        return par

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
        pypeit_keys = super(P200TSPECSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)

        if ftype in ['pinhole','bias']:
            # No pinhole frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'lamp_off')
        if ftype == 'standard':
            return good_exp & ((fitstbl['idname'] == 'object') | (fitstbl['idname'] == 'Object'))
        if ftype in ['pixelflat', 'trace']:
            return good_exp & (fitstbl['target'] == 'lamp_on')
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
        msgs.info("Custom bad pixel mask for TSPEC")
        bpm_img = self.empty_bpm(filename, det, shape=shape)
        # ToDo: Build a custom bad pixel mask.

        return bpm_img

    @property
    def norders(self):
        return 5

    @property
    def order_spat_pos(self):
        ord_spat_pos = np.array([0.3096, 0.4863, 0.6406, 0.7813, 0.9424])
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
        Args:
            order_vec (np.ndarray):
            binning (optional):

        Returns:
            np.ndarray:

        """
        norders = order_vec.size
        return np.full(norders, 0.37)
