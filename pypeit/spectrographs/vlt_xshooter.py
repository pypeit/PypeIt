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

        def_keys[0]['target'] = 'OBJECT'
        def_keys[0]['exptime'] = 'EXPTIME'
        #def_keys[0]['hatch'] = 'TRAPDOOR'
        def_keys[0]['idname'] = 'HIERARCH.ESO.DPR.TYPE'
        def_keys[0]['airmass'] = 'HIERARCH.ESO.TEL.AIRM.START'
        def_keys[0]['decker'] = 'HIERARCH.ESO.INS.OPTI4.NAME'
        def_keys[0]['naxis0'] = 'NAXIS2'
        def_keys[0]['naxis1'] = 'NAXIS1'

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
        par['calibrations']['slits']['polyorder'] = 5 # Might want 6 or 7
        par['calibrations']['slits']['maxshift'] = 0.5  # Trace crude
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
        return head_keys

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
        Set default parameters for Keck LRISb reductions.
        """
        par = pypeitpar.PypeItPar()
        par['rdx']['spectrograph'] = 'vlt_xshooter_nir'
        # Use the ARMS pipeline
        par['rdx']['pipeline'] = 'ARMED'
        # Set wave tilts order
        par['calibrations']['slits']['polyorder'] = 5 # Might want 6 or 7
        par['calibrations']['slits']['maxshift'] = 0.5  # Trace crude
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
        return head_keys

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
        #debugger.set_trace() # THIS NEEDS TO BE DEVELOPED
        arcparam['lamps'] = ['OH_triplespec']
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']



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


