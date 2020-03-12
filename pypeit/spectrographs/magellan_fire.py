"""
Module for Magellan/FIRE specific codes

Important Notes:

    - If you are reducing old FIRE data (before the broken happened
      in 2016), please change the ord_spat_pos array (see lines from
      ~220 to ~230)

"""
from pkg_resources import resource_filename
import numpy as np
from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph



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
        self.spectrograph = 'magellan_fire_base'
        self.telescope = telescopes.MagellanTelescopePar()

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for VLT XSHOOTER reductions.
        """
        par = pypeitpar.PypeItPar()
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
        meta['dichroic'] = dict(ext=0, card=None, default='default')
        meta['binning'] = dict(ext=0, card=None, default='1,1')

        meta['mjd'] = dict(ext=0, card='ACQTIME')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='AIRMASS')
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card='GRISM')
        meta['idname'] = dict(ext=0, card='OBSTYPE')

        # Ingest
        self.meta = meta



class MagellanFIREEchelleSpectrograph(MagellanFIRESpectrograph):
    """
    Child to handle Magellan/FIRE Echelle data

    .. note::
        For FIRE Echelle, we usually use high gain and SUTR read mode.
        The exposure time is usually around 900s. The detector
        parameters below are based on such mode. Standard star and
        calibrations are usually use Fowler 1 read mode in which case
        the read noise is ~20 electron.

    """
    def __init__(self):
        # Get it started
        super(MagellanFIREEchelleSpectrograph, self).__init__()
        #TODO Rename this magallen_fire_echelle??
        self.spectrograph = 'magellan_fire'
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
                            platescale      = 0.18,
                            darkcurr        = 0.01,
                            #saturation      = 20000., # high gain is 20000 ADU, low gain is 32000 ADU
                            saturation      = 100000., # This is an arbitrary value.
                            nonlinear       = 1.0, # high gain mode, low gain is 0.875
                            numamplifiers   = 1,
                            gain            = 1.2, # high gain mode, low gain is 3.8 e-/DN
                            ronoise         = 5.0, # for high gain mode and SUTR read modes with exptime ~ 900s
                            datasec         = '[5:2044,5:2044]',
                            oscansec        = '[5:2044,:5]'
                            #datasec         = '[:,:]',
                            #oscansec        = '[:,:]'
                            )]

    @property
    def pypeline(self):
        return 'Echelle'

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_fire'
        # No overscan
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'
        # Wavelengths
        # 1D wavelength solution with OH lines
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['sigdetect']=[5,10,10,10,10,20,30,30,30,30,30,10,30,30,60,30,30,10,20,30,10]
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=[3,3,3,2,4,4,4,3,4,4,4,3,4,4,4,4,4,4,6,6,4]
        par['calibrations']['wavelengths']['lamps'] = ['OH_FIRE_Echelle']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.35
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire_echelle.fits'
        par['calibrations']['wavelengths']['match_toler']=30.0

        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 6
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['slitedges']['edge_thresh'] = 10.
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_order'] = 5
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['pca_order'] = 3

        # Scienceimage default parameters
        par['reduce'] = pypeitpar.ReducePar()
        # Always flux calibrate, starting with default parameters
        #par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [20, None]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]

        # Sensitivity function parameters
        # Sensitivity function parameters
        par['sensfunc']['algorithm'] = 'IR'
        par['sensfunc']['polyorder'] = 8
        # place holder for telgrid file
        par['sensfunc']['IR']['telgridfile'] = resource_filename('pypeit', '/data/telluric/TelFit_MaunaKea_3100_26100_R20000.fits')


        return par


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

    @property
    def norders(self):
        return 21

    @property
    def order_spat_pos(self):
        # ToDo: We somehow need to automate this.
        ## For OLD data, i.e. before 2017
        #ord_spat_pos = np.array([0.06054688, 0.14160156, 0.17089844, 0.22753906, 0.27539062,
        #                         0.32128906, 0.36474609, 0.40673828, 0.45019531, 0.48974609,
        #                         0.52978516, 0.56054688, 0.59814453, 0.63378906, 0.66503906,
        #                         0.70019531, 0.7421875 , 0.77978516, 0.82763672, 0.87109375,
        #                         0.9296875])
        ## For NEW data
        ord_spat_pos = np.array([0.078125, 0.13769531, 0.19189453, 0.24414062, 0.29296875,
                                 0.34179688, 0.38330078, 0.42724609, 0.46582031, 0.50439453,
                                 0.54199219, 0.57763672, 0.61279297, 0.6484375 , 0.68457031,
                                 0.71875   , 0.75439453, 0.79443359, 0.83789062, 0.88671875,
                                 0.94091797])
        return ord_spat_pos

    @property
    def orders(self):
        return np.arange(31, 10, -1, dtype=int)

    @property
    def spec_min_max(self):
        spec_max = np.asarray([2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,2048,
                               2048,2048,2048,2048,2048])
        spec_min = np.asarray([ 500,   0,   0,   0,   0,   0,   0,    0,   0,   0,  0,   0,   0,   0,   0,   0,
                                  0,   0,   0,   0,   0])
        return np.vstack((spec_min, spec_max))

    def order_platescale(self, order_vec, binning=None):
        """
        FIRE has no binning

        Args:
            order_vec (np.ndarray):
            binning (optional):

        Returns:
            np.ndarray:

        """
        norders = order_vec.size
        return np.full(norders, 0.15)

    @property
    def dloglam(self):
        # This number was determined using the resolution and sampling quoted on the FIRE website
        R = 6000.0 * 2.7
        dloglam = 1.0 / R / np.log(10.0)
        return dloglam

    @property
    def loglam_minmax(self):
        return np.log10(8000.0), np.log10(25700)


class MagellanFIRELONGSpectrograph(MagellanFIRESpectrograph):
    """
    Child to handle Magellan/FIRE high-throughput data

    .. note::
        For FIRE longslit, science data are usually taken with SUTR readout mode with ~600s exposure
        (at least for quasar hunting people) and the readout noise is ~6 e-

    """
    def __init__(self):
        # Get it started
        super(MagellanFIRELONGSpectrograph, self).__init__()
        self.spectrograph = 'magellan_fire_long'
        self.camera = 'FIRE'
        self.numhead = 1
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.15,
                            darkcurr        = 0.01,
                            saturation      = 320000., #32000 for low gain, I set to a higher value to keep data in K-band
                            nonlinear       = 0.875,
                            numamplifiers   = 1,
                            gain            = 3.8,
                            ronoise         = 6.0, # SUTR readout mode with exposure~600s
                            datasec         = '[5:2044, 900:1250]',
                            oscansec        = '[:5, 900:1250]'
                            )]

    def default_pypeit_par(self):
        """
        Set default parameters.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'magellan_fire_long'

        # No overscan
        for key in par['calibrations'].keys():
            if 'frame' in key:
                par['calibrations'][key]['process']['overscan'] = 'none'
        # Wavelengths
        # 1D wavelength solution with arc lines
        par['calibrations']['wavelengths']['rms_threshold'] = 1.0
        par['calibrations']['wavelengths']['sigdetect']=3
        par['calibrations']['wavelengths']['fwhm'] = 20
        par['calibrations']['wavelengths']['n_first']=2
        par['calibrations']['wavelengths']['n_final']=4
        par['calibrations']['wavelengths']['lamps'] = ['ArI', 'ArII', 'ThAr', 'NeI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['method'] = 'full_template'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'magellan_fire_long.fits'
        par['calibrations']['wavelengths']['match_toler']=5.0

        # Set slits and tilts parameters
        par['calibrations']['tilts']['tracethresh'] = 5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['sync_predict'] = 'nearest'

        # Scienceimage parameters
        par['reduce']['findobj']['sig_thresh'] = 5
        #par['reduce']['maxnumber'] = 2
        par['reduce']['findobj']['find_trim_edge'] = [50,50]
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibratePar()
        # Do not correct for flexure
        par['flexure'] = None
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['standardframe']['exprng'] = [None, 60]
        par['calibrations']['arcframe']['exprng'] = [1, 50]
        par['calibrations']['darkframe']['exprng'] = [20, None]
        par['scienceframe']['exprng'] = [20, None]
        return par

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
            return good_exp & (fitstbl['idname'] == 'Arc')
        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

