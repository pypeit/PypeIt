""" Module for VLT X-Shooter
"""
import glob

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit import utils
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels
from IPython import embed

from pkg_resources import resource_filename

from pypeit import debugger

class VLTXShooterSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    def __init__(self):
        # Get it started
        super(VLTXShooterSpectrograph, self).__init__()
        self.spectrograph = 'vlt_xshooter_base'
        self.telescope = telescopes.VLTTelescopePar()

    @property
    def pypeline(self):
        return 'Echelle'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for VLT XSHOOTER reductions.
        """
        par = pypeitpar.PypeItPar()
        # Correct for flexure using the default approach
#        par['flexure'] = pypeitpar.FlexurePar()
        # Right now turn off flexure compensation
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
        meta['ra'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])  # Need to convert to : separated
        meta['dec'] = dict(card=None, compound=True, required_ftypes=['science', 'standard'])
        meta['target'] = dict(ext=0, card='OBJECT')
        meta['binning'] = dict(card=None, compound=True)

        meta['mjd'] = dict(ext=0, card='MJD-OBS')
        meta['exptime'] = dict(ext=0, card='EXPTIME')
        meta['airmass'] = dict(ext=0, card='HIERARCH ESO TEL AIRM START', required_ftypes=['science', 'standard'])
        # Extras for config and frametyping
        meta['dispname'] = dict(ext=0, card=None, default='default')
        meta['idname'] = dict(ext=0, card='HIERARCH ESO DPR CATG')
        meta['arm'] = dict(ext=0, card='HIERARCH ESO SEQ ARM')

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            if 'HIERARCH ESO DET WIN1 BINX' in headarr[0]:
                binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            else:
                binspatial = 1
            if 'HIERARCH ESO DET WIN1 BINY' in headarr[0]:
                binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            else:
                binspec = 1
            binning = parse.binning2string(binspec, binspatial)
            return binning
        elif meta_key in ['ra', 'dec']:
            try:  # Calibs do not have RA values
                coord = SkyCoord(ra=headarr[0]['RA'], dec=headarr[0]['DEC'], unit='deg')
            except:
                return None
            if meta_key == 'ra':
                return coord.ra.to_string(unit=units.hour,sep=':',pad=True,precision=2)
            else:
                return coord.dec.to_string(sep=':',pad=True,alwayssign=True,precision=1)
        else:
            msgs.error("Not ready for this compound meta")

    def configuration_keys(self):
        return ['arm']


    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        # TODO: Allow for 'sky' frame type, for now include sky in
        # 'science' category
        if ftype == 'science':
            return good_exp & ((fitstbl['idname'] == 'SCIENCE')
                                | (fitstbl['target'] == 'STD,TELLURIC')
                                | (fitstbl['target'] == 'STD,SKY'))
        if ftype == 'standard':
            return good_exp & (fitstbl['target'] == 'STD,FLUX')
        if ftype == 'bias':
            return good_exp & (fitstbl['target'] == 'BIAS')
        if ftype == 'dark':
            return good_exp & (fitstbl['target'] == 'DARK')
        if ftype in ['pixelflat', 'trace']:
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype in ['arc', 'tilt']:
            return good_exp & (fitstbl['target'] == 'LAMP,WAVE')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

    @property
    def norders(self):
        return None


class VLTXShooterNIRSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    def __init__(self):
        # Get it started
        super(VLTXShooterNIRSpectrograph, self).__init__()
        self.spectrograph = 'vlt_xshooter_nir'
        self.camera = 'XShooter_NIR'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 1,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.197, # average between order 11 & 30, see manual
                            darkcurr        = 0.0,
                            saturation      = 2.0e5, # I think saturation may never be a problem here since there are many DITs
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 2.12, #
                            ronoise         = 8.0, # ?? more precise value? #TODO the read noise is exposure time  dependent and should be grabbed from header
                            datasec         = '[4:2044,4:]', # These are all unbinned pixels
                            # EMA: No real overscan for XSHOOTER-NIR: 
                            # See Table 6 in http://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/VLT-MAN-ESO-14650-4942_P103v1.pdf
                            # The overscan region below contains only zeros
                            # ToDo should we just set it as empty?
                            oscansec        = '[4:2044,1:3]', # These are all unbinned pixels.
                            suffix          = '_NIR'
                            )]
        self.numhead = 1


    def default_pypeit_par(self):
        """
        Set default parameters for XSHOOTER NIR reductions.
        """
        par = VLTXShooterSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'vlt_xshooter_nir'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slitedges']['edge_thresh'] = 50.
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['fit_min_spec_length'] = 0.5
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # Tilt parameters
        par['calibrations']['tilts']['rm_continuum'] = True
        par['calibrations']['tilts']['tracethresh'] =  25.0
        par['calibrations']['tilts']['maxdev_tracefit'] =  0.04
        par['calibrations']['tilts']['maxdev2d'] =  0.04
        par['calibrations']['tilts']['spat_order'] =  3
        par['calibrations']['tilts']['spec_order'] =  4

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['OH_XSHOOTER']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.25
        par['calibrations']['wavelengths']['sigdetect'] = 10.0
        par['calibrations']['wavelengths']['fwhm'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = 4
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 5
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = False
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10


        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        # Is this needed below?
        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] = 'nothing'
        # TODO tune up LA COSMICS parameters here for X-shooter as tellurics are being excessively masked

        # Extraction
        par['scienceimage']['skysub']['bspline_spacing'] = 0.8
        par['scienceimage']['skysub']['global_sky_std']  = False # Do not perform global sky subtraction for standard stars
        par['scienceimage']['extraction']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        par['scienceimage']['findobj']['trace_npoly'] = 8
        par['scienceimage']['findobj']['find_npoly_cont'] = 0  # Continnum order for determining thresholds
        par['scienceimage']['findobj']['find_cont_fit'] = False  # Don't attempt to fit a continuum to the trace rectified image

        # The settings below enable X-shooter dark subtraction from the traceframe and pixelflatframe, but enforce
        # that this bias won't be subtracted from other images. It is a hack for now, because eventually we want to
        # perform this operation with the dark frame class, and we want to attach individual sets of darks to specific
        # images.
        par['calibrations']['biasframe']['useframe'] = 'bias'
        par['calibrations']['traceframe']['process']['bias'] = 'force'
        par['calibrations']['pixelflatframe']['process']['bias'] = 'force'
        par['calibrations']['arcframe']['process']['bias'] = 'skip'
        par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        par['calibrations']['standardframe']['process']['bias'] = 'skip'
        par['scienceframe']['process']['bias'] = 'skip'


        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for an VLT/XSHOOTER exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = { '0.INSTRUME': 'XSHOOTER',
                            '0.HIERARCH ESO SEQ ARM': 'NIR',
                            '0.NAXIS': 2 }
        super(VLTXShooterNIRSpectrograph, self).check_headers(headers,
                                                              expected_values=expected_values)

    def init_meta(self):
        """
        Meta data specific to VLT NIR

        Returns:

        """
        super(VLTXShooterNIRSpectrograph, self).init_meta()
        # No binning in the NIR
        self.meta['binning'] = dict(card=None, default='1,1')

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI5 NAME')

    def pypeit_file_keys(self):
        pypeit_keys = super(VLTXShooterSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

    def bpm(self, filename, det, shape=None):
        """
        Override parent bpm function with BPM specific to X-ShooterNIR.

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

        bpm_img = self.empty_bpm(filename, det, shape=shape)
        if det == 1:
            bpm_dir = resource_filename('pypeit', 'data/static_calibs/vlt_xshoooter/')
            try :
                bpm_loc = np.loadtxt(bpm_dir+'BP_MAP_RP_NIR.dat',usecols=(0,1))
            except IOError :
                msgs.warn('BP_MAP_RP_NIR.dat not present in the static database')
                bpm_fits = fits.open(bpm_dir+'BP_MAP_RP_NIR.fits.gz')
                # ToDo: this depends on datasec, biassec, specflip, and specaxis
                #       and should become able to adapt to these parameters.
                # Flipping and shifting BPM to match the PypeIt format
                y_shift = -2
                x_shift = 18
                bpm_data = np.flipud(bpm_fits[0].data)
                y_len = len(bpm_data[:,0])
                x_len = len(bpm_data[0,:])
                bpm_data_pypeit = np.full( ((y_len+abs(y_shift)),(x_len+abs(x_shift))) , 0)
                bpm_data_pypeit[:-abs(y_shift),:-abs(x_shift)] = bpm_data_pypeit[:-abs(y_shift),:-abs(x_shift)] + bpm_data
                bpm_data_pypeit = np.roll(bpm_data_pypeit,-y_shift,axis=0)
                bpm_data_pypeit = np.roll(bpm_data_pypeit,x_shift,axis=1)
                filt_bpm = bpm_data_pypeit[1:y_len,1:x_len]>100.
                y_bpm, x_bpm = np.where(filt_bpm)
                bpm_loc = np.array([y_bpm,x_bpm]).T
                np.savetxt(bpm_dir+'BP_MAP_RP_NIR.dat', bpm_loc, fmt=['%d','%d'])
            finally :
                bpm_img[bpm_loc[:,0].astype(int),bpm_loc[:,1].astype(int)] = 1.

        return bpm_img


    @property
    def norders(self):
        return 16

    @property
    def order_spat_pos(self):
        ord_spat_pos = np.array([0.08284662, 0.1483813 , 0.21158701, 0.27261607,
                                   0.33141317, 0.38813936, 0.44310197, 0.49637422,
                                   0.54839496, 0.59948157, 0.65005956, 0.70074477,
                                   0.75240745, 0.80622583, 0.86391259, 0.9280528 ])
        return ord_spat_pos

    @property
    def orders(self):
        return np.arange(26, 10, -1, dtype=int)


    @property
    def spec_min_max(self):
        spec_max = np.asarray([1467,1502,1540, 1580,1620,1665,1720, 1770,1825,1895, 1966, 2000,2000,2000,2000,2000])
        spec_min = np.asarray([420 ,390 , 370,  345, 315, 285, 248,  210, 165, 115,   63,   10,   0,   0,   0,   0])
        return np.vstack((spec_min, spec_max))


    def order_platescale(self, order_vec, binning=None):
        """
        Returns the spatial plate scale in arcseconds for each order

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

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # NIR: 52.4 pixels (0.210"/pix) at order 11 to 59.9 pixels (0.184"/pix) at order 26.

        # Right now I just assume a simple linear trend
        plate_scale = 0.184 + (order_vec - 26)*(0.184-0.210)/(26 - 11)
        return plate_scale


    @property
    def dloglam(self):
        # This number was computed by taking the mean of the dloglam for all the X-shooter orders. The specific
        # loglam across the orders deviates from this value by +-6% from this first to final order
        return 1.93724e-5

    @property
    def loglam_minmax(self):
        return np.log10(9500.0), np.log10(26000)



class VLTXShooterVISSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    def __init__(self):
        # Get it started
        super(VLTXShooterVISSpectrograph, self).__init__()
        self.spectrograph = 'vlt_xshooter_vis'
        self.camera = 'XShooter_VIS'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = False,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.16, # average from order 17 and order 30, see manual
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 0.595, # FITS format is flipped: PrimaryHDU  (2106, 4000) w/respect to Python
                            ronoise         = 3.1, # raw unbinned images are (4000,2106) (spec, spat)
                            datasec='[:,11:2058]',  # pre and oscan are in the spatial direction
                            oscansec='[:,2059:2106]',
                            suffix          = '_VIS'
                            )]
        self.numhead = 1



    def default_pypeit_par(self):
        """
        Set default parameters for VLT XSHOOTER VIS reductions.
        """
        par = VLTXShooterSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'vlt_xshooter_vis'

        # Adjustments to parameters for VIS
        par['calibrations']['arcframe']['process']['overscan'] = 'median'
        # X-SHOOTER arcs/tilts are also have different binning with bias frames
        par['calibrations']['arcframe']['process']['bias'] = 'skip'
        par['calibrations']['tiltframe']['process']['bias'] = 'skip'
        # Don't use the biases for the arcs or flats since it appears to be a different amplifier readout
        par['calibrations']['traceframe']['process']['overscan'] = 'median'

        par['calibrations']['slitedges']['edge_thresh'] = 8.0
        par['calibrations']['slitedges']['fit_order'] = 8
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        # These are the defaults
        par['calibrations']['tilts']['tracethresh'] = 15
        par['calibrations']['tilts']['spat_order'] =  3
        par['calibrations']['tilts']['spec_order'] =  5 # [5, 5, 5] + 12*[7] # + [5]

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_XSHOOTER_VIS']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.50 # This is for 1x1 binning. TODO GET BINNING SORTED OUT!!
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        par['calibrations']['wavelengths']['n_final'] = [3] + 13*[4] + [3]
        par['calibrations']['wavelengths']['fwhm'] = 11.0 # This is for 1x1 binning. Needs to be divided by binning for binned data!!
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # ToDo the arxived solution is for 1x1 binning. It needs to be generalized for different binning!
        #par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.json'
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.fits'
        par['calibrations']['wavelengths']['cc_thresh'] = 0.50
        par['calibrations']['wavelengths']['cc_local_thresh'] = 0.50
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 4
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # Flats
        par['calibrations']['flatfield']['illumflatten'] = True
        par['calibrations']['flatfield']['tweak_slits_thresh'] = 0.90
        par['calibrations']['flatfield']['tweak_slits_maxfrac'] = 0.10

        # Extraction
        par['scienceimage']['skysub']['bspline_spacing'] = 0.5
        par['scienceimage']['skysub']['global_sky_std'] = False
        par['scienceimage']['extraction']['model_full_slit'] = True # local sky subtraction operates on entire slit
        par['scienceimage']['findobj']['find_trim_edge'] = [3,3] # Mask 3 edges pixels since the slit is short, insted of default (5,5)
        par['scienceimage']['findobj']['find_npoly_cont'] = 0       # Continnum order for determining thresholds
        par['scienceimage']['findobj']['find_cont_fit'] = False # Don't attempt to fit a continuum to the trace rectified image

        # Right now we are using the overscan and not biases becuase the standards are read with a different read mode and we don't
        # yet have the option to use different sets of biases for different standards, or use the overscan for standards but not for science frames
        par['scienceframe']['useframe'] ='overscan'

        return par

    def init_meta(self):
        """
        Meta data specific to VLT NIR

        Returns:

        """
        super(VLTXShooterVISSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        # dispangle and filter1 are not defined for Shane Kast Blue

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI4 NAME')

    @property
    def norders(self):
        return 15

    @property
    def order_spat_pos(self):
        ord_spat_pos = np.array([0.13540436, 0.21055672, 0.2817009, 0.34907542,
                                 0.41289127, 0.4733839 , 0.53072208, 0.58509916,
                                 0.63671413, 0.685754, 0.73236772, 0.77676367,
                                 0.8191196 , 0.85968302, 0.89877932])
        return ord_spat_pos

    @property
    def orders(self):
        return np.arange(30, 15, -1, dtype=int)

    @property
    def spec_min_max(self):
        spec_max = np.asarray([4000]*14 + [3000])
        spec_min = np.asarray([2000,1000] + [0]*13)
        return np.vstack((spec_min, spec_max))


    def order_platescale(self, order_vec, binning=None):
        """
        Returns the plate scale in arcseconds for each order

        Args:
            order_vec (np.ndarray): Order numbers
            binning (optional):

        Returns:
            np.ndarray: Platescale

        """
        # VIS has no binning, but for an instrument with binning we would do this
        binspectral, binspatial = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # VIS: 65.9 pixels (0.167"/pix) at order 17 to 72.0 pixels (0.153"/pix) at order 30.

        # Right now I just assume a simple linear trend
        plate_scale = 0.153 + (order_vec - 30)*(0.153-0.167)/(30 - 17)
        return plate_scale*binspatial



    @property
    def dloglam(self):
        # This number was computed by taking the mean of the dloglam for all the X-shooter orders. The specific
        # loglam across the orders deviates from this value by +-7% from this first to final order. This is the
        # unbinned value. It was actually measured  to be 1.69207e-5  from a 2x1 data and then divided by two.
        return 8.46035e-06

    @property
    def loglam_minmax(self):
        return np.log10(5000.0), np.log10(11000)

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
        bpm_img = self.empty_bpm(filename, det, shape=shape)
        shape = bpm_img.shape
        #
        # ToDo Ema: This is just a workaround to deal with
        # different binning. I guess binspatial and binspectral
        # should be passed in.
        if shape[0]<3000.:
            binspectral_bpm=2
        else:
            binspectral_bpm=1
        if shape[1]<1500.:
            binspatial_bpm=2
        else:
            binspatial_bpm=1

        if det == 1:
            bpm_img[2912//binspectral_bpm:,842//binspatial_bpm:844//binspatial_bpm] = 1.
        return bpm_img


class VLTXShooterUVBSpectrograph(VLTXShooterSpectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    def __init__(self):
        # Get it started
        super(VLTXShooterUVBSpectrograph, self).__init__()
        self.spectrograph = 'vlt_xshooter_uvb'
        self.camera = 'XShooter_UVB'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            specaxis        = 0,
                            specflip        = True,
                            spatflip        = True,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            # average from order 14 and order 24, see manual
                            platescale      = 0.161,
                            darkcurr        = 0.0,
                            saturation      = 65000.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.61,
                            ronoise         = 2.60,
                            datasec         = '[:,49:2096]', # '[49:2000,1:2999]',
                            oscansec        = '[:,1:48]', # '[1:48, 1:2999]',
                            suffix          = '_UVB'
                            )]
        self.numhead = 1


#    @staticmethod
    def default_pypeit_par(self):
        """
        Set default parameters for VLT XSHOOTER UVB reductions.
        """
        par = VLTXShooterSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'vlt_xshooter_uvb'

        # Adjustments to slit and tilts for UVB
        par['calibrations']['slitedges']['edge_thresh'] = 8.
        par['calibrations']['slitedges']['max_shift_adj'] = 0.5
        par['calibrations']['slitedges']['trace_thresh'] = 10.
        par['calibrations']['slitedges']['left_right_pca'] = True
        par['calibrations']['slitedges']['length_range'] = 0.3

        par['calibrations']['arcframe']['process']['overscan'] = 'median'
        par['calibrations']['traceframe']['process']['overscan'] = 'median'
        # X-SHOOTER UVB arcs/tilts have different binning with bias frames
        par['calibrations']['arcframe']['process']['bias'] = 'skip'
        par['calibrations']['tiltframe']['process']['bias'] = 'skip'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['lamps'] = ['ThAr_XSHOOTER_UVB']
        # TODO: This is a KLUDGE; default_pypeit_par should be a
        # staticmethod meaning it should not depend on self
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['rms_threshold'] = 0.50 # This is for 1x1 binning. TODO GET BINNING SORTED OUT!!
        par['calibrations']['wavelengths']['sigdetect'] = 5.0
        # Reidentification parameters
        par['calibrations']['wavelengths']['method'] = 'reidentify'
        # ToDo the arxived solution is for 1x1 binning. It needs to be generalized for different binning!
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_uvb1x1_iraf.json'
        par['calibrations']['wavelengths']['ech_fix_format'] = True
        # Echelle parameters
        par['calibrations']['wavelengths']['echelle'] = True
        par['calibrations']['wavelengths']['ech_nspec_coeff'] = 4
        par['calibrations']['wavelengths']['ech_norder_coeff'] = 5
        par['calibrations']['wavelengths']['ech_sigrej'] = 3.0

        # TODO FIX THIS TO USE BIASES!!
        par['scienceframe']['useframe'] ='overscan'

        return par

    def init_meta(self):
        """
        Meta data specific to VLT NIR

        Returns:

        """
        super(VLTXShooterUVBSpectrograph, self).init_meta()
        # Add the name of the dispersing element
        # dispangle and filter1 are not defined for Shane Kast Blue

        # Required
        self.meta['decker'] = dict(ext=0, card='HIERARCH ESO INS OPTI3 NAME')

    def slit2order(self, islit, nslit):

        """
        Parameters
        ----------
        islit: int, float, or string, slit number

        Returns
        -------
        order: int
        """
        msgs.error("Refactor to use slit_spat_pos!!")

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

        orders = np.arange(24,12,-1, dtype=int)
        return orders[islit]


    def order_platescale(self, binning = None):
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
        msgs.error("REFACTOR")

        binspectral, binspatial = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # UVB: 65.9 pixels (0.167“/pix) at order 14 to 70.8 pixels (0.155”/pix) at order 24

        # Right now I just took the average
        return np.full(self.norders, 0.161)*binspatial

    def bpm(self, filename, det, shape=None):
        """
        Override parent bpm function with BPM specific to X-Shooter UVB.

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
        bpm_img = self.empty_bpm(filename, det, shape=shape)
        if det == 1:
            # TODO: This is for the 1x1 binning it should
            # change for other binning
            bpm_img[:2369,1326:1328] = 1.

        return bpm_img



