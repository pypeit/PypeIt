""" Module for Shane/Kast specific codes
"""
from __future__ import absolute_import, division, print_function



import numpy as np


from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import fsort

from pypeit import debugger


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

    @staticmethod
    def kast_default_pypeit_par():
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = pypeitpar.PypeItPar()
        # TODO: Make self.spectrograph a class attribute?
        # Use the ARMS pipeline
        par['rdx']['pipeline'] = 'ARMS'
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 5
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 1
        # Set wave tilts order
        par['calibrations']['tilts']['order'] = 2
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always flux calibrate, starting with default parameters
        par['fluxcalib'] = pypeitpar.FluxCalibrationPar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()
        return par

    def kast_header_keys(self):
        """
        Provide the relevant header keywords
        """
        def_keys = self.default_header_keys()

        # A time stamp of the observation; used to find calibrations
        # proximate to science frames. The units of this value are
        # specified by fits+timeunit below
        def_keys[0]['time'] = 'TSEC'

        # Image size
        # TODO: Check ordering
        def_keys[0]['naxis0'] = 'NAXIS2'
        def_keys[0]['naxis1'] = 'NAXIS1'

        # Lamp names and statuses
        def_keys[0]['lampname01'] = 'LAMPNAM1'
        def_keys[0]['lampstat01'] = 'LAMPSTA1'
        def_keys[0]['lampname02'] = 'LAMPNAM2'
        def_keys[0]['lampstat02'] = 'LAMPSTA2'
        def_keys[0]['lampname03'] = 'LAMPNAM3'
        def_keys[0]['lampstat03'] = 'LAMPSTA3'
        def_keys[0]['lampname04'] = 'LAMPNAM4'
        def_keys[0]['lampstat04'] = 'LAMPSTA4'
        def_keys[0]['lampname05'] = 'LAMPNAM5'
        def_keys[0]['lampstat05'] = 'LAMPSTA5'
        def_keys[0]['lampname06'] = 'LAMPNAMA'
        def_keys[0]['lampstat06'] = 'LAMPSTAA'
        def_keys[0]['lampname07'] = 'LAMPNAMB'
        def_keys[0]['lampstat07'] = 'LAMPSTAB'
        def_keys[0]['lampname08'] = 'LAMPNAMC'
        def_keys[0]['lampstat08'] = 'LAMPSTAC'
        def_keys[0]['lampname09'] = 'LAMPNAMD'
        def_keys[0]['lampstat09'] = 'LAMPSTAD'
        def_keys[0]['lampname10'] = 'LAMPNAME'
        def_keys[0]['lampstat10'] = 'LAMPSTAE'
        def_keys[0]['lampname11'] = 'LAMPNAMF'
        def_keys[0]['lampstat11'] = 'LAMPSTAF'
        def_keys[0]['lampname12'] = 'LAMPNAMG'
        def_keys[0]['lampstat12'] = 'LAMPSTAG'
        def_keys[0]['lampname13'] = 'LAMPNAMH'
        def_keys[0]['lampstat13'] = 'LAMPSTAH'
        def_keys[0]['lampname14'] = 'LAMPNAMI'
        def_keys[0]['lampstat14'] = 'LAMPSTAI'
        def_keys[0]['lampname15'] = 'LAMPNAMJ'
        def_keys[0]['lampstat15'] = 'LAMPSTAJ'
        def_keys[0]['lampname16'] = 'LAMPNAMK'
        def_keys[0]['lampstat16'] = 'LAMPSTAK'

        # Dichroic and decker
        def_keys[0]['dichroic'] = 'BSPLIT_N'
        def_keys[0]['decker'] = 'SLIT_N'

        return def_keys

    def kast_cond_dict(self, ftype):
        """
        Create dict for image typing (from header)

        Args:
            ftype: str

        Returns:

        """
        cond_dict = {}

        if ftype == 'science':
            cond_dict['condition1'] = 'lampstat01=off&lampstat02=off&lampstat03=off' \
                                      '&lampstat04=off&lampstat05=off&lampstat06=off' \
                                      '&lampstat07=off&lampstat08=off&lampstat09=off' \
                                      '&lampstat10=off&lampstat11=off&lampstat12=off' \
                                      '&lampstat13=off&lampstat14=off&lampstat15=off' \
                                      '&lampstat16=off'
            cond_dict['condition2'] = 'exptime>1'
        elif ftype == 'bias':
            cond_dict['condition1'] = 'exptime<1'
        elif ftype == 'pixelflat':
            cond_dict['condition1'] = 'lampstat01=on|lampstat02=on|lampstat03=on' \
                                      '|lampstat04=on|lampstat05=on'
            cond_dict['condition2'] = 'exptime>0'
        elif ftype == 'pinhole':
            cond_dict['condition1'] = 'exptime>99999999'
        elif ftype == 'trace':
            cond_dict['condition1'] = 'lampstat01=on|lampstat02=on|lampstat03=on|lampstat04=on' \
                                      '|lampstat05=on'
            cond_dict['condition2'] = 'exptime>0'
        elif ftype == 'arc':
            cond_dict['condition1'] = 'lampstat06=on|lampstat07=on|lampstat08=on|lampstat09=on' \
                                      '|lampstat10=on|lampstat11=on|lampstat12=on|lampstat13=on' \
                                      '|lampstat14=on|lampstat15=on|lampstat16=on'
            cond_dict['condition2'] = 'exptime<=60'
        else:
            pass

        return cond_dict

    def check_ftype(self, ftype, fitstbl):
        """
        Check the frame type

        Args:
            ftype:
            fitstbl:

        Returns:

        """
        # Load up
        cond_dict = self.kast_cond_dict(ftype)

        # Do it
        gd_chk = fsort.chk_all_conditions(fitstbl, cond_dict)

        return gd_chk

    def get_match_criteria(self):
        """Set the general matching criteria for Shane Kast."""
        match_criteria = {}
        for key in fsort.ftype_list:
            match_criteria[key] = {}

        match_criteria['standard']['match'] = {}
        match_criteria['standard']['match']['naxis0'] = '=0'
        match_criteria['standard']['match']['naxis1'] = '=0'

        match_criteria['bias']['match'] = {}
        match_criteria['bias']['match']['naxis0'] = '=0'
        match_criteria['bias']['match']['naxis1'] = '=0'

        match_criteria['pixelflat']['match'] = {}
        match_criteria['pixelflat']['match']['naxis0'] = '=0'
        match_criteria['pixelflat']['match']['naxis1'] = '=0'
        match_criteria['pixelflat']['match']['decker'] = ''

        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'
        match_criteria['trace']['match']['decker'] = ''

        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'

        return match_criteria

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
                pypeitpar.DetectorPar(dataext         = 0,
                                     dispaxis        = 1,
                                     dispflip       = False,
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
                                     datasec         = [ '[1:1024,:]', '[1025:2048,:]'],
                                     oscansec        = [ '[2050:2080,:]', '[2081:2111,:]'],
                                     suffix          = '_blue'
                                     )]
        self.numhead = 1
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        self.sky_file = 'sky_kastb_600.fits'

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = self.kast_default_pypeit_par()
        par['rdx']['spectrograph'] = 'shane_kast_blue'
        return par

    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name (replace any spaces with underscores)
        chk_dict[1]['DSENSOR'] = 'Fairchild CCD 3041 2Kx2K'
        return chk_dict

    def header_keys(self):
        """
        Header keys specific to shane_kast_blue

        Returns:

        """
        head_keys = self.kast_header_keys()
        # Add the name of the dispersing element
        head_keys[0]['dispname'] = 'GRISM_N'
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
        arcparam['lamps'] = ['CdI','HgI','HeI']
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
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

        # TODO: NEED TO CHECK ORIENTATION OF DATASEC AND OSCANSEC ARE
        # CORRECT!!!!

        # Get it started
        super(ShaneKastRedSpectrograph, self).__init__()
        self.spectrograph = 'shane_kast_red'
        self.camera = 'KASTr'
        self.detector = [
                # Detector 1
                pypeitpar.DetectorPar(dataext         = 0,
                                     dispaxis        = 0,
                                     dispflip=False,
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
                                     datasec         = ['[:,2:511]', '[:,513:525]'],
                                     oscansec        = ['[:,527:625]', '[:,627:725]'],
                                     suffix          = '_red'
                                     )]
        self.numhead = 1
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        # self.sky_file = ?

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = self.kast_default_pypeit_par()
        #
        par['rdx']['spectrograph'] = 'shane_kast_red'
        #
        return par

    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[1]['DSENSOR'] = '2k x 4k Hamamatsu'
        return chk_dict

    def header_keys(self):
        """
        Header keys specific to shane_kast_red

        Returns:

        """
        head_keys = self.kast_header_keys()
        head_keys[0]['filter1'] = 'RDFILT_N'
        head_keys[0]['dispname'] = 'GRATING_N'
        head_keys[0]['dispangle'] = 'GRTILT_P'
        return head_keys

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
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_nsigl'] = 30.         # Minimum signififance
        arcparam['wvmnx'] = [3000.,11000.]  # Guess at wavelength range
        # These parameters influence how the fts are done by pypeit.core.wavecal.fitting.iterative_fitting
        arcparam['match_toler'] = 3         # Matcing tolerance (pixels)
        arcparam['func'] = 'legendre'       # Function for fitting
        arcparam['n_first'] = 2             # Order of polynomial for first fit
        arcparam['n_final'] = 4             # Order of polynomial for final fit
        arcparam['nsig_rej'] = 2            # Number of sigma for rejection
        arcparam['nsig_rej_final'] = 3.0    # Number of sigma for rejection (final fit)




#        if disperser == '600/7500':
#            arcparam['disp']=1.30
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
#            arcparam['wvmnx'][0] = 5000.
#            arcparam['n_first']=2 # Should be able to lock on
#        elif disperser == '1200/5000':
#            arcparam['disp']=0.63
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
#            arcparam['wvmnx'][0] = 5000.
#            arcparam['n_first']=2 # Should be able to lock on
#            arcparam['wv_cen'] = 6600.
#        else:
#            msgs.error('Not ready for this disperser {:s}!'.format(disperser))


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
                pypeitpar.DetectorPar(dataext         = 0,
                                     dispaxis        = 1,
                                     dispflip=False,
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
                                     oscansec        = '[1203:1232,:]',
                                     suffix          = '_red'
                                     )]
        self.numhead = 1
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        # self.sky_file = ?

    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = self.kast_default_pypeit_par()
        par['rdx']['spectrograph'] = 'shane_kast_red_ret'
        #
        par['calibrations']['pixelflatframe']['number'] = 3
        par['calibrations']['traceframe']['number'] = 3
        return par

    def header_keys(self):
        """
        Header keys specific to shane_kast_blue

        Returns:

        """
        head_keys = self.kast_header_keys()
        head_keys[0]['filter1'] = 'RDFILT_N'
        head_keys[0]['dispname'] = 'GRATNG_N'
        head_keys[0]['dispangle'] = 'GRTILT_P'
        #
        return head_keys

    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        # Check the CCD name
        chk_dict[1]['DSENSOR'] = 'Ret 400x1200'
        return chk_dict

    def get_match_criteria(self):
        match_criteria = super(ShaneKastRedRetSpectrograph, self).get_match_criteria()
        # Add more
        match_criteria['standard']['match']['dispangle'] = '|<=20'
        match_criteria['pixelflat']['match']['dispangle'] = '|<=20'
        match_criteria['arc']['match']['dispangle'] = '|<=10'
        return match_criteria

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
        arcparam['nonlinear_counts'] = self.detector[0]['nonlinear']*self.detector[0]['saturation']
        arcparam['min_nsigl'] = 30.         # Minimum signififance
        arcparam['wvmnx'] = [3000.,11000.]  # Guess at wavelength range
        # These parameters influence how the fts are done by pypeit.core.wavecal.fitting.iterative_fitting
        arcparam['match_toler'] = 3         # Matcing tolerance (pixels)
        arcparam['func'] = 'legendre'       # Function for fitting
        arcparam['n_first'] = 2             # Order of polynomial for first fit
        arcparam['n_final'] = 4             # Order of polynomial for final fit
        arcparam['nsig_rej'] = 2            # Number of sigma for rejection
        arcparam['nsig_rej_final'] = 3.0    # Number of sigma for rejection (final fit)



#        if disperser == '600/7500':
#            arcparam['disp']=2.35
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
#            arcparam['wvmnx'][0] = 5000.
#            arcparam['n_first']=2 # Should be able to lock on
#        elif disperser == '1200/5000':
#            arcparam['disp']=1.17
#            arcparam['b1']= 1./arcparam['disp']/msarc_shape[0] / binspectral
#            arcparam['wvmnx'][0] = 5000.
#            arcparam['n_first']=2 # Should be able to lock on
#        else:
#            msgs.error('Not ready for this disperser {:s}!'.format(disperser))

