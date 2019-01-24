""" Module for VLT X-Shooter
"""
from __future__ import absolute_import, division, print_function

import glob

import numpy as np
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import units

from pypeit import msgs
from pypeit import telescopes
from pypeit.core import parse
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit.core import pixels

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

        # Ingest
        self.meta = meta

    def compound_meta(self, headarr, meta_key):
        if meta_key == 'binning':
            binspatial = headarr[0]['HIERARCH ESO DET WIN1 BINX']
            binspec = headarr[0]['HIERARCH ESO DET WIN1 BINY']
            binning = parse.binning2string(binspatial, binspec)
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
        return []

    def pypeit_file_keys(self):
        pypeit_keys = super(VLTXShooterSpectrograph, self).pypeit_file_keys()
        pypeit_keys += ['calib', 'comb_id', 'bkg_id']
        return pypeit_keys

#    def parse_binning(self, inp, **kwargs):
#        """
#        Get the pixel binning for an image.
#
#        Args:
#            inp (:obj:`str`, `astropy.io.fits.Header`):
#                String providing the file name to read, or the relevant
#                header object.
#
#        Returns:
#            str: String representation of the binning.  The ordering is
#            as provided in the header, regardless of which axis is
#            designated as the dispersion axis.  It is expected that this
#            be used with :func:`pypeit.core.parse.sec2slice` to setup
#            the data and overscane sections of the image data.
#
#        Raises:
#            PypeItError:
#                Raised if `inp` is not one of the accepted types.
#        """
#        # Get the header
#        if isinstance(inp, str):
#            hdu = fits.open(inp)
#            hdr = hdu[0].header
#        elif isinstance(inp, fits.Header):
#            hdr = inp
#        else:
#            msgs.error('Input must be a filename or fits.Header object')
#
#        # TODO: This is a hack.  These two keywords don't exist for the
#        # test NIR file in tests/test_load_images.py.  Can the NIR data
#        # be binned?  What are the appropriate keywords?
#        try:
#            return '{0},{1}'.format(hdr['HIERARCH ESO DET WIN1 BINX'],
#                                    hdr['HIERARCH ESO DET WIN1 BINY'])
#        except:
#            return '1,1'

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
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & ((fitstbl['target'] == 'LAMP,DFLAT')
                               | (fitstbl['target'] == 'LAMP,QFLAT')
                               | (fitstbl['target'] == 'LAMP,FLAT'))
        if ftype == 'pinhole':
            # Don't type pinhole
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & (fitstbl['target'] == 'LAMP,WAVE')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)

#    def get_match_criteria(self):
#        # TODO: Matching needs to be looked at...
#        match_criteria = {}
#        for key in framematch.FrameTypeBitMask().keys():
#            match_criteria[key] = {}
#
#        match_criteria['standard']['match'] = {}
#        # Bias
#        match_criteria['bias']['match'] = {}
#        match_criteria['bias']['match']['binning'] = ''
#        # Pixelflat
#        match_criteria['pixelflat']['match'] = {}
#        match_criteria['pixelflat']['match']['binning'] = ''
#        # Traceflat
#        match_criteria['trace']['match'] = {}
#        match_criteria['trace']['match']['binning'] = ''
#        # Arc
#        match_criteria['arc']['match'] = {}
#        # Return
#        return match_criteria


    @property
    def norders(self):
        return None

    def slit2order(self, islit):
        pass

    def order_vec(self):
        return self.slit2order(np.arange(self.norders))



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
                            gain            = 2.12,
                            ronoise         = 8.0, # ?? more precise value?
                            datasec         = '[4:,4:2044]',
                            # EMA: No real overscan for XSHOOTER-NIR: 
                            # See Table 6 in http://www.eso.org/sci/facilities/paranal/instruments/xshooter/doc/VLT-MAN-ESO-14650-4942_P103v1.pdf
                            # The overscan region below contains only zeros
                            # ToDo should we just set it as empty?
                            oscansec        = '[1:3,4:2044]',
                            suffix          = '_NIR'
                            )]
        self.numhead = 1

    @property
    def norders(self):
        return 16


    def default_pypeit_par(self):
        """
        Set default parameters for XSHOOTER NIR reductions.
        """
        par = VLTXShooterSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'vlt_xshooter_nir'

        # Adjustments to slit and tilts for NIR
        par['calibrations']['slits']['sigdetect'] = 120.
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['pcatype'] = 'order'

        # Tilt parameters
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
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_nir.json'
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

        # Extraction
        par['scienceimage']['bspline_spacing'] = 0.8
        par['scienceimage']['model_full_slit'] = True  # local sky subtraction operates on entire slit
        # Do not bias subtract
        par['scienceframe']['useframe'] ='none'
        # This is a hack for now until we can specify for each image type what to do. Bias currently
        # controls everything
        par['calibrations']['biasframe']['useframe'] = 'none'

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

    '''
    def header_keys(self):
        hdr_keys = super(VLTXShooterNIRSpectrograph, self).header_keys()
        hdr_keys[0]['decker'] = 'HIERARCH ESO INS OPTI5 NAME'
        hdr_keys[0]['utc'] = 'HIERARCH ESO DET EXP UTC'
        return hdr_keys
    '''

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


    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
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

        self.empty_bpm(shape=shape, filename=filename, det=det)
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
                self.bpm_img[bpm_loc[:,0].astype(int),bpm_loc[:,1].astype(int)] = 1.

        return self.bpm_img

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

        orders = np.arange(26,10,-1, dtype=int)
        return orders[islit]

    def order_platescale(self, binning = None):
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
        slit_vec = np.arange(self.norders)
        order_vec = self.slit2order(slit_vec)
        plate_scale = 0.184 + (order_vec - 26)*(0.184-0.210)/(26 - 11)
        return plate_scale

    def slit_minmax(self, nslits, binspectral=1):

        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        spec_max = np.asarray([1467,1502,1540, 1580,1620,1665,1720, 1770,1825,1895, 1966, 2000,2000,2000,2000,2000])
        spec_min = np.asarray([420 ,390 , 370,  345, 315, 285, 248,  210, 165, 115,   63,   10,   0,   0,   0,   0])

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

        nslits = tslits_dict['lcen'].shape[1]
        if nslits != self.norders:
            msgs.error('There is a problem with your slit bounadries. You have nslits={:d} orders, whereas NIR has norders={:d}'.format(nslits,self.norders))
        # These are the order boundaries determined by eye by JFH. 2025 is used as the maximum as the upper bit is not illuminated
        order_max = [1467,1502,1540, 1580,1620,1665,1720, 1770,1825,1895, 1966, 2000,2000,2000,2000,2000]
        order_min = [420 ,390 , 370,  345, 315, 285, 248,  210, 165, 115,   63,   10,   0,   0,   0,   0]
        for islit in range(nslits):
            orderbad = (slitmask == islit) & ((spec_img < order_min[islit]) | (spec_img > order_max[islit]))
            slitmask[orderbad] = -1
        return slitmask


    def wavegrid(self, binning=None):

        # Define the grid for VLT-XSHOOTER NIR
        dloglam = 1.93724e-5
        # This number was computed by taking the mean of the dloglam for all the X-shooter orders. The specific
        # loglam across the orders deviates from this value by +-6% from this first to final order
        logmin = np.log10(9500.0)
        logmax = np.log10(26000)
        ngrid = int(np.ceil((logmax - logmin) / dloglam))
        osamp = 1.0
        loglam_grid = logmin + (dloglam / osamp) * np.arange(int(np.ceil(osamp * ngrid)))

        return np.power(10.0,loglam_grid)




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
                            gain            = 0.595,
                            ronoise         = 3.1,
                            datasec         = '[11:2058,1:]', #'[29:1970,1:]',
                            oscansec        = '[2059:2106,1:]',
                            suffix          = '_VIS'
                            )]
        self.numhead = 1


    @property
    def norders(self):
        return 15

    def default_pypeit_par(self):
        """
        Set default parameters for VLT XSHOOTER VIS reductions.
        """
        par = VLTXShooterSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'vlt_xshooter_vis'

        # Adjustments to parameters for VIS
        par['calibrations']['arcframe']['process']['overscan'] = 'median'
        # Don't use the biases for the arcs or flats since it appears to be a different amplifier readout
        par['calibrations']['arcframe']['useframe']= 'overscan'
        par['calibrations']['traceframe']['process']['overscan'] = 'median'
        par['calibrations']['traceframe']['useframe']= 'overscan'
        par['calibrations']['biasframe']['useframe']= 'overscan'
        # TODO THIS IS STUPID. biasframe currently determines behvior for everyone. See Issue # 554

        par['calibrations']['slits']['sigdetect'] = 8.0
        par['calibrations']['slits']['pcatype'] = 'pixel'
        par['calibrations']['slits']['polyorder'] = 6
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['number'] = -1
        #par['calibrations']['slits']['fracignore'] = 0.01

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
        par['calibrations']['wavelengths']['reid_arxiv'] = 'vlt_xshooter_vis1x1.json'
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
        par['scienceimage']['bspline_spacing'] = 0.8
        par['scienceimage']['model_full_slit'] = True # local sky subtraction operates on entire slit
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

        orders = np.arange(30,15,-1, dtype=int)

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

        # VIS has no binning, but for an instrument with binning we would do this
        binspatial, binspectral = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # VIS: 65.9 pixels (0.167"/pix) at order 17 to 72.0 pixels (0.153"/pix) at order 30.

        # Right now I just assume a simple linear trend
        slit_vec = np.arange(self.norders)
        order_vec = self.slit2order(slit_vec)
        plate_scale = 0.153 + (order_vec - 30)*(0.153-0.167)/(30 - 17)
        return plate_scale*binspatial


    def slit_minmax(self, nslits, binspectral=1):

        spec_max = np.asarray([4000]*14 + [3000])//binspectral
        spec_min = np.asarray([2000,1000] + [0]*13)//binspectral

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
        nslits = tslits_dict['lcen'].shape[1]
        if nslits != self.norders:
            msgs.error('Not all the orders were identified!')

        slitmask = pixels.slit_pixels(tslits_dict['lcen'], tslits_dict['rcen'], tslits_dict['nspat'], pad=pad)
        spec_img = np.outer(np.arange(tslits_dict['nspec'], dtype=int), np.ones(tslits_dict['nspat'], dtype=int))  # spectral position everywhere along image

        binning = tslits_dict['binning']
        binspatial, binspectral = parse.parse_binning(binning)
        # These are the order boundaries determined by eye by JFH.
        order_max = np.asarray([4000]*14 + [3000])//binspectral
        order_min = np.asarray([2000,1000] + [0]*13)//binspectral
        # TODO add binning adjustments to these
        for iorder in range(self.norders):
            orderbad = (slitmask == iorder) & ((spec_img < order_min[iorder]) | (spec_img > order_max[iorder]))
            slitmask[orderbad] = -1

        return slitmask


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

        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[2912//binspectral_bpm:,842//binspatial_bpm:844//binspatial_bpm] = 1.
        return self.bpm_img


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
                            datasec         = '[49:2096,1:]', # '[49:2000,1:2999]',
                            oscansec        = '[1:48,1:]', # '[1:48, 1:2999]',
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
        par['calibrations']['slits']['sigdetect'] = 8.
        par['calibrations']['slits']['pcatype'] = 'pixel'
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['number'] = -1

        par['calibrations']['arcframe']['process']['overscan'] = 'median'
        par['calibrations']['traceframe']['process']['overscan'] = 'median'

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

        binspatial, binspectral = parse.parse_binning(binning)

        # ToDO Either assume a linear trend or measure this
        # X-shooter manual says, but gives no exact numbers per order.
        # UVB: 65.9 pixels (0.167“/pix) at order 14 to 70.8 pixels (0.155”/pix) at order 24

        # Right now I just took the average
        return np.full(self.norders, 0.161)*binspatial

    def bpm(self, shape=None, filename=None, det=None, **null_kwargs):
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
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            # TODO: This is for the 1x1 binning it should
            # change for other binning
            self.bpm_img[:2369,1326:1328] = 1.

        return self.bpm_img



