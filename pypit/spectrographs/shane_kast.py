""" Module for Shane/Kast specific codes
"""
from __future__ import absolute_import, division, print_function

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import glob

import numpy as np

from astropy.io import fits

from pypit import msgs
from pypit.par.pypitpar import DetectorPar
from pypit.par.pypitpar import CalibrationsPar
from pypit.spectrographs import spectrograph
from pypit import telescopes
from pypit.core import arsort
from pypit.spectrographs import util

from pypit import ardebug as debugger


class ShaneKastSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Shane/Kast specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast'
        self.telescope = telescopes.ShaneTelescopePar()
        self.timeunit = 's'

    def kast_header_keys(self):
        def_keys = self.default_header_keys()
        # Update
        def_keys[0]['time'] = 'TSEC'   # A time stamp of the observation; used to find calibrations proximate to science frames. The units of this value are specified by fits+timeunit below
        def_keys[0]['naxis0'] = 'NAXIS2' # Number of pixels along the zeroth axis
        def_keys[0]['naxis1'] = 'NAXIS1' # Number of pixels along the first axis
        def_keys[0]['lampname01'] = 'LAMPNAM1' # Number of pixels along the first axis
        def_keys[0]['lampstat01'] = 'LAMPSTA1' # Number of pixels along the first axis
        def_keys[0]['lampname02'] = 'LAMPNAM2' # Number of pixels along the first axis
        def_keys[0]['lampstat02'] = 'LAMPSTA2' # Number of pixels along the first axis
        def_keys[0]['lampname03'] = 'LAMPNAM3' # Number of pixels along the first axis
        def_keys[0]['lampstat03'] = 'LAMPSTA3' # Number of pixels along the first axis
        def_keys[0]['lampname04'] = 'LAMPNAM4' # Number of pixels along the first axis
        def_keys[0]['lampstat04'] = 'LAMPSTA4' # Number of pixels along the first axis
        def_keys[0]['lampname05'] = 'LAMPNAM5' # Number of pixels along the first axis
        def_keys[0]['lampstat05'] = 'LAMPSTA5' # Number of pixels along the first axis
        def_keys[0]['lampname06'] = 'LAMPNAMA' # Number of pixels along the first axis
        def_keys[0]['lampstat06'] = 'LAMPSTAA' # Number of pixels along the first axis
        def_keys[0]['lampname07'] = 'LAMPNAMB' # Number of pixels along the first axis
        def_keys[0]['lampstat07'] = 'LAMPSTAB' # Number of pixels along the first axis
        def_keys[0]['lampname08'] = 'LAMPNAMC' # Number of pixels along the first axis
        def_keys[0]['lampstat08'] = 'LAMPSTAC' # Number of pixels along the first axis
        def_keys[0]['lampname09'] = 'LAMPNAMD' # Number of pixels along the first axis
        def_keys[0]['lampstat09'] = 'LAMPSTAD' # Number of pixels along the first axis
        def_keys[0]['lampname10'] = 'LAMPNAME' # Number of pixels along the first axis
        def_keys[0]['lampstat10'] = 'LAMPSTAE' # Number of pixels along the first axis
        def_keys[0]['lampname11'] = 'LAMPNAMF' # Number of pixels along the first axis
        def_keys[0]['lampstat11'] = 'LAMPSTAF' # Number of pixels along the first axis
        def_keys[0]['lampname12'] = 'LAMPNAMG' # Number of pixels along the first axis
        def_keys[0]['lampstat12'] = 'LAMPSTAG' # Number of pixels along the first axis
        def_keys[0]['lampname13'] = 'LAMPNAMH' # Number of pixels along the first axis
        def_keys[0]['lampstat13'] = 'LAMPSTAH' # Number of pixels along the first axis
        def_keys[0]['lampname14'] = 'LAMPNAMI' # Number of pixels along the first axis
        def_keys[0]['lampstat14'] = 'LAMPSTAI' # Number of pixels along the first axis
        def_keys[0]['lampname15'] = 'LAMPNAMJ' # Number of pixels along the first axis
        def_keys[0]['lampstat15'] = 'LAMPSTAJ' # Number of pixels along the first axis
        def_keys[0]['lampname16'] = 'LAMPNAMK' # Number of pixels along the first axis
        def_keys[0]['lampstat16'] = 'LAMPSTAK' # Number of pixels along the first axis
        #
        def_keys[0]['dichroic'] = 'BSPLIT_N' # Number of pixels along the first axis
        # Return
        return def_keys


class ShaneKastBlueSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast blue specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastBlueSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_blue'
        self.camera = 'KASTb'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.43,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 2,
                            gain            = [1.2, 1.2],
                            ronoise         = [3.7, 3.7],
                            datasec         = [ '[0:1024,:]', '[1024:2048,:]'],
                            oscansec        = [ '[2049:2080,:]', '[2080:2111,:]'],
                            suffix          = '_blue'
                            )
            ]
        self.numhead = 1
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        self.sky_file = 'sky_kastb_600.fits'

    def check_header(self, headers):
        chk_dict = {}
        chk_dict[1] = {}  # 1,2,3 indexing
        chk_dict[1]['NAXIS'] = 2                            # THIS IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['DSENSOR'] = 'Fairchild CCD 3041 2Kx2K' # Check the CCD name (replace any spaces with underscores)
        #

    def header_keys(self):
        """
        Header keys specific to shane_kast_blue

        Returns:

        """
        head_keys = self.kast_header_keys()
        head_keys[0]['decker'] = 'SLIT_N'  # Which decker is being used
        head_keys[0]['dispname'] = 'GRISM_N' # Number of pixels along the first axis
        #
        return head_keys

    def check_ftype(self, ftype, fitstbl):
        gd_chk = np.ones(len(fitstbl), dtype=bool)
        cond_dict = {}

        if ftype == 'science':
            cond_dict['condition1'] = 'lampstat01=off&lampstat02=off&lampstat03=off&lampstat04=off&lampstat05=off&lampstat06=off&lampstat07=off&lampstat08=off&lampstat09=off&lampstat10=off&lampstat11=off&lampstat12=off&lampstat13=off&lampstat14=off&lampstat15=off&lampstat16=off'
            cond_dict['condition2'] = 'exptime>1'
        elif ftype == 'bias':
            cond_dict['condition1'] = 'exptime<1'
        elif ftype == 'pixelflat':
            cond_dict['condition1'] = 'lampstat01=on|lampstat02=on|lampstat03=on|lampstat04=on|lampstat05=on'
            cond_dict['condition2'] = 'exptime>0'
        elif ftype == 'pinhole':
            gd_chk[:] = False
            return gd_chk
        elif ftype == 'trace':
            cond_dict['condition1'] = 'lampstat01=on|lampstat02=on|lampstat03=on|lampstat04=on|lampstat05=on'
            cond_dict['condition2'] = 'exptime>0'
        elif ftype == 'arc':
            cond_dict['condition1'] = 'lampstat06=on|lampstat07=on|lampstat08=on|lampstat09=on|lampstat10=on|lampstat11=on|lampstat12=on|lampstat13=on|lampstat14=on|lampstat15=on|lampstat16=on'
            cond_dict['condition2'] = 'exptime<=60'
        else:
            pass

        # Do it
        gd_chk &= arsort.chk_all_conditions(fitstbl, cond_dict)

        return gd_chk

    def _set_calib_par(self, user_supplied=None):
        self.calib_par = CalibrationsPar()

    def get_match_criteria(self):
        match_criteria = {}
        for key in arsort.ftype_list:
            match_criteria[key] = {}
        # Science
        match_criteria['science']['number'] = 1
        # Standard
        match_criteria['standard']['number'] = 1  # Can be over-ruled by flux calibrate = False
        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'
        # Bias
        match_criteria['bias']['number'] = 5  # Can be over-ruled by flux calibrate = False
        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'
        # Pixelflat
        match_criteria['pixelflat']['number'] = 5  # Can be over-ruled by flux calibrate = False
        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'
        match_criteria['pixelflat']['match']['decker'] = ''
        # Traceflat
        match_criteria['trace']['number'] = 5  # Can be over-ruled by flux calibrate = False
        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'
        match_criteria['trace']['match']['decker'] = ''
        # Arc
        match_criteria['arc']['number'] = 1
        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'

        # Return
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
        arcparam['lamps'] = ['CdI','HgI','HeI']
        if disperser == '600/4310':
            arcparam['disp']=1.02
            arcparam['b1']=6.88935788e-04
            arcparam['b2']=-2.38634231e-08
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4250.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


class ShaneKastRedSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast red specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastRedSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_red'
        self.camera = 'KASTr'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 0,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.43,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 2,
                            gain            = [1.9, 1.9],
                            ronoise         = [3.8, 3.8],
                            datasec         = ['[1:511,:]', '[512:525,:]'],
                            oscansec        = ['[526:625,:]', '[626:725,:]'],
                            suffix          = '_red'
                            )
            ]
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        # self.sky_file = ?

    def setup_arcparam(self, arcparam, disperser=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
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
        arcparam['lamps'] = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=1.30
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=0.63
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
            arcparam['wv_cen'] = 6600.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

class ShaneKastRedRetSpectrograph(ShaneKastSpectrograph):
    """
    Child to handle Shane/Kast red specific code
    """
    def __init__(self):
        # Get it started
        super(ShaneKastRedRetSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_red_ret'
        # WARNING: This is not unique wrt ShaneKastRed...
        self.camera = 'KASTr'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.774,
                            darkcurr        = 0.0,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 3.0,
                            ronoise         = 12.5,
                            oscansec        = '[1202:1232,:]',
                            suffix          = '_red'
                            )
            ]
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        # self.sky_file = ?

    def setup_arcparam(self, arcparam, disperser=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
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
        arcparam['lamps'] = ['NeI','HgI','HeI','ArI']
        if disperser == '600/7500':
            arcparam['disp']=2.35
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=1.17
            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

