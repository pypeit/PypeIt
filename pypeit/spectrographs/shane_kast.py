""" Module for Shane/Kast specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np


from pypeit import msgs
from pypeit import telescopes
from pypeit.core import framematch
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph

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
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast reductions.
        """
        par = pypeitpar.PypeItPar()
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
        # Set the default exposure time ranges for the frame typing
        par['calibrations']['biasframe']['exprng'] = [None, 1]
        par['calibrations']['darkframe']['exprng'] = [999999, None]     # No dark frames
        par['calibrations']['pinholeframe']['exprng'] = [999999, None]  # No pinhole frames
        par['calibrations']['pixelflatframe']['exprng'] = [0, None]
        par['calibrations']['traceframe']['exprng'] = [0, None]
        par['calibrations']['arcframe']['exprng'] = [None, 61]
        par['calibrations']['standardframe']['exprng'] = [1, 61]
        par['scienceframe']['exprng'] = [61, None]
        return par

    def header_keys(self):
        """
        Provide the relevant header keywords
        """

        hdr_keys = {}
        hdr_keys[0] = {}

        hdr_keys[0]['target'] = 'OBJECT'
        hdr_keys[0]['idname'] = 'OBSTYPE'
        hdr_keys[0]['time'] = 'TSEC'
        hdr_keys[0]['date'] = 'DATE'
        hdr_keys[0]['ra'] = 'RA'
        hdr_keys[0]['dec'] = 'DEC'
        hdr_keys[0]['airmass'] = 'AIRMASS'
        hdr_keys[0]['binning'] = 'BINNING'
        hdr_keys[0]['exptime'] = 'EXPTIME'
        hdr_keys[0]['decker'] = 'SLIT_N'
        hdr_keys[0]['dichroic'] = 'BSPLIT_N'
        
        #dispname is different for all three spectrographs

        hdr_keys[0]['naxis0'] = 'NAXIS2'
        hdr_keys[0]['naxis1'] = 'NAXIS1'

        lamp_names = [ '1', '2', '3', '4', '5',
                       'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J', 'K' ]

        for kk,lamp_name in enumerate(lamp_names):
            hdr_keys[0]['lampstat{:02d}'.format(kk+1)] = 'LAMPSTA{0}'.format(lamp_name)

        return hdr_keys

    # Uses parent metadata keys

    def check_frame_type(self, ftype, fitstbl, exprng=None):
        """
        Check for frames of the provided type.
        """
        good_exp = framematch.check_frame_exptime(fitstbl['exptime'], exprng)
        if ftype in ['science', 'standard']:
            return good_exp & self.lamps(fitstbl, 'off')
#            \
#                        & np.array([ t not in ['Arcs', 'Bias', 'Dome Flat']
#                                        for t in fitstbl['target']])
        if ftype == 'bias':
            return good_exp # & (fitstbl['target'] == 'Bias')
        if ftype == 'pixelflat' or ftype == 'trace':
            # Flats and trace frames are typed together
            return good_exp & self.lamps(fitstbl, 'dome') # & (fitstbl['target'] == 'Dome Flat')
        if ftype == 'pinhole' or ftype == 'dark':
            # Don't type pinhole or dark frames
            return np.zeros(len(fitstbl), dtype=bool)
        if ftype == 'arc':
            return good_exp & self.lamps(fitstbl, 'arcs')#  & (fitstbl['target'] == 'Arcs')

        msgs.warn('Cannot determine if frames are of type {0}.'.format(ftype))
        return np.zeros(len(fitstbl), dtype=bool)
  
    def lamps(self, fitstbl, status):
        """
        Check the lamp status.

        Args:
            fitstbl (:obj:`astropy.table.Table`):
                The table with the fits header meta data.
            status (:obj:`str`):
                The status to check.  Can be `off`, `arcs`, or `dome`.
        
        Returns:
            numpy.ndarray: A boolean array selecting fits files that
            meet the selected lamp status.

        Raises:
            ValueError:
                Raised if the status is not one of the valid options.
        """
        if status == 'off':
            # Check if all are off
            return np.all(np.array([ (fitstbl[k] == 'off') | (fitstbl[k] == 'None')
                                        for k in fitstbl.keys() if 'lampstat' in k]), axis=0)
        if status == 'arcs':
            # Check if any arc lamps are on
            arc_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(6,17) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in arc_lamp_stat]), axis=0)
        if status == 'dome':
            # Check if any dome lamps are on
            dome_lamp_stat = [ 'lampstat{0:02d}'.format(i) for i in range(1,6) ]
            return np.any(np.array([ fitstbl[k] == 'on' for k in fitstbl.keys()
                                            if k in dome_lamp_stat]), axis=0)
        raise ValueError('No implementation for status = {0}'.format(status))
        
    def get_match_criteria(self):
        """Set the general matching criteria for Shane Kast."""
        match_criteria = {}
        for key in framematch.FrameTypeBitMask().keys():
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
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            dispaxis        = 1,
                            dispflip        = False,
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

    #@staticmethod
    def default_pypeit_par(self):
        """
        Set default parameters for Shane Kast Blue reductions.
        """
        par = ShaneKastSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'shane_kast_blue'

        # 1D wavelength solution
        par['calibrations']['wavelengths']['rms_threshold'] = 0.15  # Might be grating dependent..
        par['calibrations']['wavelengths']['min_nsig'] = 30.
        par['calibrations']['wavelengths']['lowest_nsig'] = 10.
        par['calibrations']['wavelengths']['lamps'] = ['CdI','HgI','HeI']
        par['calibrations']['wavelengths']['nonlinear_counts'] = self.detector[0]['nonlinear'] * self.detector[0]['saturation']
        par['calibrations']['wavelengths']['n_first'] = 1

        return par

    def check_headers(self, headers):
        """
        Check headers match expectations for a Shane Kast blue exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.DSENSOR': 'Fairchild CCD 3041 2Kx2K' }
        super(ShaneKastBlueSpectrograph, self).check_headers(headers,
                                                             expected_values=expected_values)

    def header_keys(self):
        """
        Header keys specific to shane_kast_blue

        Returns:

        """
        hdr_keys = super(ShaneKastBlueSpectrograph, self).header_keys()
        # Add the name of the dispersing element
        # dispangle and filter1 are not defined for Shane Kast Blue
        hdr_keys[0]['dispname'] = 'GRISM_N'
        return hdr_keys

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
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            dispaxis        = 0,
                            dispflip        = False,
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

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast Red reductions.
        """
        par = ShaneKastSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'shane_kast_red'
        return par

    def check_header(self, headers):
        """
        Check headers match expectations for a Shane Kast red exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.DSENSOR': '2k x 4k Hamamatsu' }
        super(ShaneKastRedSpectrograph, self).check_headers(headers,
                                                            expected_values=expected_values)

    def header_keys(self):
        """
        Header keys specific to shane_kast_red

        Returns:

        """
        hdr_keys = super(ShaneKastRedSpectrograph, self).header_keys()
        hdr_keys[0]['dispname'] = 'GRATING_N'
        hdr_keys[0]['filter1'] = 'RDFILT_N'
        hdr_keys[0]['dispangle'] = 'GRTILT_P'
        return hdr_keys

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
        arcparam['min_nsig'] = 30.0         # Minimum signififance
        arcparam['lowest_nsig'] = 10.0      # Min significance for arc lines to be used
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
                pypeitpar.DetectorPar(
                            dataext         = 0,
                            dispaxis        = 1,
                            dispflip        = False,
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
        # TODO: Can we change suffix to be unique wrt ShaneKastRed?
        self.numhead = 1
        # Uses timeunit from parent class
        # Uses default primary_hdrext
        # self.sky_file = ?

    @staticmethod
    def default_pypeit_par():
        """
        Set default parameters for Shane Kast Red Ret reductions.
        """
        par = ShaneKastSpectrograph.default_pypeit_par()
        par['rdx']['spectrograph'] = 'shane_kast_red_ret'
        par['calibrations']['pixelflatframe']['number'] = 3
        par['calibrations']['traceframe']['number'] = 3
        return par

    def check_header(self, headers):
        """
        Check headers match expectations for a Shane Kast Red Ret
        exposure.

        See also
        :func:`pypeit.spectrographs.spectrograph.Spectrograph.check_headers`.

        Args:
            headers (list):
                A list of headers read from a fits file
        """
        expected_values = {   '0.NAXIS': 2,
                            '0.DSENSOR': 'Ret 400x1200' }
        super(ShaneKastRedRetSpectrograph, self).check_headers(headers,
                                                               expected_values=expected_values)
    
    def header_keys(self):
        """
        Header keys specific to shane_kast_red_ret

        Returns:

        """
        hdr_keys = super(ShaneKastRedRetSpectrograph, self).header_keys()
        hdr_keys[0]['dispname'] = 'GRATNG_N'
        hdr_keys[0]['filter1'] = 'RDFILT_N'
        hdr_keys[0]['dispangle'] = 'GRTILT_P'
        return hdr_keys

    def get_match_criteria(self):
        # Get the parent matching criteria ...
        match_criteria = super(ShaneKastRedRetSpectrograph, self).get_match_criteria()
        # ... add more
        match_criteria['standard']['match']['dispangle'] = '|<=20'
        match_criteria['pixelflat']['match']['dispangle'] = '|<=20'
        match_criteria['arc']['match']['dispangle'] = '|<=10'
        match_criteria['arc']['match']['decker'] = 'any'
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
        arcparam['min_nsig'] = 30.         # Minimum signififance
        arcparam['lowest_nsig'] = 10.0      # Min significance for arc lines to be used
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

