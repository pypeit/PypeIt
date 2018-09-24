""" Module for VLT X-Shooter
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np
from astropy.io import fits

from pypeit import msgs
from pypeit.core import parse
from pypeit.par.pypeitpar import DetectorPar
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import fsort

from pypeit import debugger

class VLTXShooterSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle VLT/XSHOOTER specific code
    """
    def __init__(self):
        # Get it started
        super(VLTXShooterSpectrograph, self).__init__()
        self.spectrograph = 'vlt_xshooter_vis'
        self.telescope = telescopes.VLTTelescopePar()

    def xshooter_header_keys(self):
        def_keys = self.default_header_keys()

        def_keys[0]['target']  = 'OBJECT'                     # Header keyword for the name given by the observer to a given frame
        def_keys[0]['idname']  = 'HIERARCH ESO DPR CATG'      # The keyword that identifies the frame type (i.e. bias, flat, etc.)
        def_keys[0]['time']    = 'MJD-OBS'                    # The time stamp of the observation (i.e. decimal MJD)
        def_keys[0]['date']    = 'DATE-OBS'                   # The UT date of the observation which is used for heliocentric (in the format YYYY-MM-DD  or  YYYY-MM-DDTHH:MM:SS.SS)
        def_keys[0]['ra']      = 'RA'                         # Right Ascension of the target
        def_keys[0]['dec']     = 'DEC'                        # Declination of the target
        def_keys[0]['airmass'] = 'HIERARCH ESO TEL AIRM START'# Airmass at start of observation
        def_keys[0]['binning'] = 'HIERARCH ESO DET WIN1 BINX' # Binning along X-axis (HIERARCH ESO DET WIN1 BINY along Y)
        def_keys[0]['exptime'] = 'EXPTIME'                    # Exposure time keyword
        def_keys[0]['naxis0']  = 'NAXIS2'
        def_keys[0]['naxis1']  = 'NAXIS1'
        def_keys[0]['decker']  = 'HIERARCH ESO INS OPTI3 NAME'# FOR UVB
        def_keys[0]['decker']  = 'HIERARCH ESO INS OPTI4 NAME'# FOR VIS
        def_keys[0]['decker']  = 'HIERARCH ESO INS OPTI5 NAME'# FOR NIR
        #
        def_keys[0]['binspatial']  = 'HIERARCH ESO DET WIN1 BINX' # Binning along X-axis (HIERARCH ESO DET WIN1 BINY along Y)
        def_keys[0]['binspectral'] = 'HIERARCH ESO DET WIN1 BINY' # Binning along X-axis (HIERARCH ESO DET WIN1 BINY along Y)

        ## def_keys[0]['utc'] = 'HIERARCH ESO DET FRAM UTC'


        # TODO: Should do something with the lamps

        return def_keys

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
                DetectorPar(dataext         = 0,
                            dispaxis        = 0,
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
                            datasec         = '[1:2000,10:2058]',
                            oscansec        = '[1:2000, 2060:2106]',
                            suffix          = '_VIS'
                            ),
            ]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        #self.sky_file = 'sky_LRISb_600.fits'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for VLT XSHOOTER reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'vlt_xshooter_vis'
        # Use the ARMED pipeline
        par['rdx']['pipeline'] = 'ARMED'
        # Set wave tilts order
        # par['calibrations']['slits']['polyorder'] = 5 # Might want 6 or 7
        par['calibrations']['slits']['maxshift'] = 0.5  # Trace crude
        # par['calibrations']['slits']['polyorder'] = 5
        # par['calibrations']['slits']['maxshift'] = 3.

        par['calibrations']['slits']['sigdetect'] = 310.
        par['calibrations']['slits']['pcatype'] = 'pixel'
        par['calibrations']['slits']['number'] = -1


        # par['calibrations']['tilts']['tracethresh'] = [100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100, 100]

        # par['calibrations']['tilts']['tracethresh'] = [50, 50, 60, 60, 2000]
        # par['calibrations']['traceframe']['process']['overscan'] = 'median'
        # par['calibrations']['arcframe']['process']['overscan'] = 'median'

        #par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # Always sky subtract, starting with default parameters
        # par['skysubtract'] = pypeitpar.SkySubtractionPar()
        # Always flux calibrate, starting with default parameters
        #par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def check_header(self):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[1]['HIERARCH.ESO.DET.CHIP1.NAME'] = 'MIT/LL CCID-20'
        return chk_dict

    def header_keys(self):
        head_keys = self.xshooter_header_keys()
        if self.spectrograph == 'vlt_xshooter_uvb':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI3 NAME'# FOR UVB
        if self.spectrograph == 'vlt_xshooter_vis':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI4 NAME'# FOR VIS
        if self.spectrograph == 'vlt_xshooter_nir':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI5 NAME'# FOR NIR
        return head_keys

    def get_match_criteria(self):
        """Set the general matching criteria for Xshooter vis."""
        match_criteria = {}
        for key in fsort.ftype_list:
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['pixelflat']['match'] = {}
        match_criteria['trace']['match'] = {}
        match_criteria['arc']['match'] = {}
        match_criteria['bias']['match'] = {}

        return match_criteria

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        debugger.set_trace() # THIS NEEDS TO BE DEVELOPED
        arcparam['lamps'] = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI', 'CdI', 'HgI']
        if disperser == '600/4000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        elif disperser == '400/3400':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.02
            arcparam['b1']= 2.72694493e-04
            arcparam['b2']= -5.30717321e-09
            arcparam['wvmnx'][1] = 6000.
        elif disperser == '300/5000':
            arcparam['n_first'] = 2
            arcparam['wv_cen'] = 4500.
            arcparam['disp'] = 1.43
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


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
        self.empty_bpm(shape=shape, filename=filename, det=det)
        if det == 1:
            self.bpm_img[1456:, 841:845] = 1.

        return self.bpm_img



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
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.197, # average between order 11 and order 30 see manual
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 2.12,
                            ronoise         = 8.0, # ?? more precise value?
                            datasec         = '[20:,4:2044]',
                            oscansec        = '[4:20,4:2044]',
                            suffix          = '_NIR'
                            ),
            ]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        #self.sky_file = 'sky_LRISb_600.fits'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Xshooter NIR reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'vlt_xshooter_nir'
        # Use the ARMED pipeline
        par['rdx']['pipeline'] = 'ARMED'
        # Set wave tilts order
        par['calibrations']['slits']['sigdetect'] = 500.
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 0.5
        par['calibrations']['slits']['pcatype'] = 'pixel'
        par['calibrations']['tilts']['tracethresh'] = [50,50,50,50,50,50,50,50,50, 50, 50, 60, 60, 2000,2000,6000]
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()

        par['scienceframe']['process']['sigclip'] = 20.0
        par['scienceframe']['process']['satpix'] ='nothing'


        #par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # Always sky subtract, starting with default parameters
        # par['skysubtract'] = pypeitpar.SkySubtractionPar()
        # Always flux calibrate, starting with default parameters
        #par['fluxcalib'] = pypeitpar.FluxCalibrationPar()

        return par

    def get_match_criteria(self):
        """Set the general matching criteria for Xshooter vis."""
        match_criteria = {}
        for key in fsort.ftype_list:
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['pixelflat']['match'] = {}
        match_criteria['trace']['match'] = {}
        match_criteria['arc']['match'] = {}
        match_criteria['bias']['match'] = {}

        return match_criteria

    def check_header(self):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[1]['HIERARCH.ESO.DET.CHIP1.NAME'] = 'MIT/LL CCID-20'
        return chk_dict

    def header_keys(self):
        head_keys = self.xshooter_header_keys()
        if self.spectrograph == 'vlt_xshooter_uvb':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI3 NAME'# FOR UVB
        if self.spectrograph == 'vlt_xshooter_vis':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI4 NAME'# FOR VIS
        if self.spectrograph == 'vlt_xshooter_nir':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI5 NAME'# FOR NIR
        return head_keys

    def setup_arcparam(self, arcparam, msarc_shape=None, 
                       disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        #debugger.set_trace() # THIS NEEDS TO BE DEVELOPED
        arcparam['lamps'] = ['OH_triplespec'] # Line lamps on
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_ampl'] = 1000.       # Minimum amplitude
        arcparam['wvmnx'] = [8000.,25000.]                     # Guess at wavelength range

#        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
#        arcparam['disp'] = 0.6                                 # Ang/unbinned pixel
#        arcparam['b1'] = 1./ arcparam['disp'] / msarc_shape[0] # Pixel fit term (binning independent)
#        arcparam['b2'] = 0.                                    # Pixel fit term
#        arcparam['lamps'] = ['OH_triplespec']                  # Line lamps on
#        arcparam['wv_cen']=17370.                              # Estimate of central wavelength
#        arcparam['disp_toler'] = 0.1                           # 10% tolerance
#        arcparam['match_toler'] = 3.                           # Matching tolerance (pixels)
#        arcparam['min_ampl'] = 1000.                           # Minimum amplitude
#        arcparam['func'] = 'legendre'                          # Function for fitting
#        arcparam['n_first'] = 1                                # Order of polynomial for first fit
#        arcparam['n_final'] = 3                                # Order of polynomial for final fit
#        arcparam['nsig_rej'] = 5.                              # Number of sigma for rejection
#        arcparam['nsig_rej_final'] = 5.0                       # Number of sigma for rejection (final fit)
#        arcparam['Nstrong'] = 20                               # Number of lines for auto-analysis



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
                DetectorPar(dataext         = 0,
                            dispaxis        = 0,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.161, # average from order 14 and order 24, see manual
                            darkcurr        = 0.0,
                            saturation      = 65000.,
                            nonlinear       = 0.86,
                            numamplifiers   = 1,
                            gain            = 1.61,
                            ronoise         = 2.60,
                            datasec         = '[1:2000,10:2058]',
                            oscansec        = '[1:2000, 2060:2106]',
                            suffix          = '_UVB'
                            ),
            ]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        #self.sky_file = 'sky_LRISb_600.fits'

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for VLT XSHOOTER reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'vlt_xshooter_uvb'
        # Use the ARMED pipeline
        par['rdx']['pipeline'] = 'ARMED'
        # Set wave tilts order
        par['calibrations']['slits']['polyorder'] = 5 # Might want 6 or 7
        par['calibrations']['slits']['maxshift'] = 0.5  # Trace crude
        par['calibrations']['slits']['polyorder'] = 5
        par['calibrations']['slits']['maxshift'] = 3.
        par['calibrations']['slits']['pcatype'] = 'order'
        #par['calibrations']['slits']['pcapar'] = [3,2,1,0]
        # Always sky subtract, starting with default parameters
        # par['skysubtract'] = pypeitpar.SkySubtractionPar()
        # Always flux calibrate, starting with default parameters
        #par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def check_header(self):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[1]['HIERARCH.ESO.DET.CHIP1.NAME'] = 'MIT/LL CCID-20'
        return chk_dict

    def header_keys(self):
        head_keys = self.xshooter_header_keys()
        if self.spectrograph == 'vlt_xshooter_uvb':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI3 NAME'# FOR UVB
        if self.spectrograph == 'vlt_xshooter_vis':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI4 NAME'# FOR VIS
        if self.spectrograph == 'vlt_xshooter_nir':
            head_keys[0]['decker'] = 'HIERARCH ESO INS OPTI5 NAME'# FOR NIR
        return head_keys

    def get_match_criteria(self):
        """Set the general matching criteria for Xshooter uvb."""
        match_criteria = {}
        for key in fsort.ftype_list:
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['pixelflat']['match'] = {}
        match_criteria['trace']['match'] = {}
        match_criteria['arc']['match'] = {}
        match_criteria['bias']['match'] = {}

        # binning
        # match_criteria['standard']['match']['binning'] = ''
        # match_criteria['bias']['match']['binning'] = ''


        return match_criteria

    def setup_arcparam(self, arcparam, disperser=None, **null_kwargs):
        """
        Setup the arc parameters

        Args:
            arcparam: dict
            disperser: str, REQUIRED
            **null_kwargs:
              Captured and never used

        Returns:
            arcparam is modified in place

        """
        debugger.set_trace() # THIS NEEDS TO BE DEVELOPED
        arcparam['lamps'] = ['Th', 'ArI']
        arcparam['n_first']=2 
        arcparam['disp']=0.2 # Ang per pixel (unbinned)
        arcparam['b1']= 0.
        arcparam['b2']= 0.
        arcparam['wvmnx'] = [2950.,5650.]
        arcparam['wv_cen'] = 4300.


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
            self.bpm_img[1456:, 841:845] = 1.

        return self.bpm_img



