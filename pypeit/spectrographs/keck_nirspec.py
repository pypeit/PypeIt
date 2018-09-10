""" Module for Keck/NIRSPEC specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes
from pypeit.core import fsort

from pypeit import debugger

class KeckNIRSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRSPEC specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRSPECSpectrograph, self).__init__()
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRSPEC'
        self.detector = [
                # Detector 1
            pypeitpar.DetectorPar(dataext         = 0,
                            dispaxis        = 0,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.193,
                            darkcurr        = 0.8,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 5.8,
                            ronoise         = 23,
                            datasec         = '[:,:]',
                            oscansec        = '[:,:]'
                            )
            ]
        self.numhead = 1
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    @staticmethod
    def nirspec_default_pypeit_par():
        """
        Set default parameters for NIRSPEC redutions
        """
        par = pypeitpar.PypeItPar()
        # TODO: Make self.spectrograph a class attribute?
        # Frame numbers
        par['calibrations']['standardframe']['number'] = 1
        par['calibrations']['biasframe']['number'] = 0
        par['calibrations']['pixelflatframe']['number'] = 5
        par['calibrations']['traceframe']['number'] = 5
        par['calibrations']['arcframe']['number'] = 0
        # Set wave tilts order
        par['calibrations']['tilts']['order'] = 2
        # Scienceimage default parameters
        par['scienceimage'] = pypeitpar.ScienceImagePar()
        # Always correct for flexure, starting with default parameters
        par['flexure'] = pypeitpar.FlexurePar()

        return par

    def nirspec_header_keys(self):
        """
        Provide the relevant header keywords
        """
        def_keys = self.default_header_keys()

        # A time stamp of the observation; used to find calibrations
        # proximate to science frames. The units of this value are
        # specified by fits+timeunit below
        def_keys[0]['time'] = 'UTC'

        # Image size
        # TODO: Check ordering
        def_keys[0]['naxis0'] = 'NAXIS2'
        def_keys[0]['naxis1'] = 'NAXIS1'

        # Lamp names and statuses
        def_keys[0]['lampstat01'] = 'NEON'
        def_keys[0]['lampstat02'] = 'ARGON'
        def_keys[0]['lampstat03'] = 'KRYPTON'
        def_keys[0]['lampstat04'] = 'XENON'
        def_keys[0]['lampstat05'] = 'ETALON'
        def_keys[0]['lampstat06'] = 'FLAT'

        def_keys[0]['dispname'] = 'DISPERS'
        def_keys[0]['hatch'] = 'CALMPOS'
        def_keys[0]['slitwid'] = 'SLITWIDT'
        def_keys[0]['slitlen'] = 'SLITLEN'
        def_keys[0]['imagetype'] = 'IMAGETYP'

        return def_keys

    def nirspec_cond_dict(self, ftype):
        """
        Create dict for image typing (from header)

        Args:
            ftype: str

        Returns:

        """
        cond_dict = {}

        if ftype == 'science':
            cond_dict['condition1'] = 'lampstat01=0&lampstat02=0&lampstat03=0' \
                                      '&lampstat04=0&lampstat05=0&lampstat06=0' # 0 is 'off' for NIRSPEC
            cond_dict['condition2'] = 'exptime>=1'
            cond_dict['condition3'] = 'hatch=0'
            cond_dict['condition4'] = 'imagetype=object'
        elif ftype == 'bias':
            cond_dict['condition1'] = 'exptime<1'
            cond_dict['condition2'] = 'hatch=0'
            cond_dict['condition3'] = 'imagetype=dark'
        elif ftype == 'pixelflat':
            cond_dict['condition1'] = 'lampstat06=1' # This is the dome flat lamp for NIRSPEC; 1 is 'on'
            cond_dict['condition2'] = 'exptime>0'
            cond_dict['condition3'] = 'hatch=1'
            cond_dict['condition4'] = 'imagetype=flatlamp'
            condict
        elif ftype == 'pinhole':
            cond_dict['condition1'] = 'exptime>99999999'
        elif ftype == 'trace':
            cond_dict['condition1'] = 'lampstat06=1' # This is the dome flat lamp for NIRSPEC; 1 is 'on'
            cond_dict['condition2'] = 'exptime>0'
            cond_dict['conditoin3'] = 'hatch=1'
        elif ftype == 'arc':
            cond_dict['condition1'] = 'lampstat01=1|lampstat02=1|lampstat03=1|lampstat04=1|lampstat05=1'
            cond_dict['condition2'] = 'hatch=1'
            cond_dict['condition3'] = 'imagetype=arclamp'
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
        cond_dict = self.nirspec_cond_dict(ftype)

        # Do it
        gd_chk = fsort.chk_all_conditions(fitstbl, cond_dict)

        return gd_chk

    def get_match_criteria(self):
        """Set the general matching criteria for Keck NIRSPEC."""
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
        match_criteria['pixelflat']['match']['dispname'] = ''

        match_criteria['trace']['match'] = {}
        match_criteria['trace']['match']['naxis0'] = '=0'
        match_criteria['trace']['match']['naxis1'] = '=0'
        match_criteria['trace']['match']['dispname'] = ''

        match_criteria['arc']['match'] = {}
        match_criteria['arc']['match']['naxis0'] = '=0'
        match_criteria['arc']['match']['naxis1'] = '=0'

        return match_criteria


    # TODO: This function is unstable to shape...
    def bpm(self, shape=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED

        Returns
        -------
        badpix : ndarray

        """
        if shape is None:
            raise ValueError('Must provide shape for Keck NIRSPEC bpm.')
        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRSPEC")
        self.bpm_img = np.zeros(shape, dtype=np.int8)
        self.bpm_img[:, :20] = 1.
        self.bpm_img[:, 1000:] = 1.

        return self.bpm_img

class KeckNIRSPECLowSpectrograph(KeckNIRSPECSpectrograph):
    """
    Child to handle NIRSPEC low-dispersion specific code
    """

    def __init__(self):
        # Get it started
        super(KeckNIRSPECLowSpectrograph, self).__init__()
        self.spectrograph = 'keck_nirspec_low'


    def default_pypeit_par(self):
        """
        Set default parameters for NIRSPEC low-dispersion reductions
        """
        par = self.nirspec_default_pypeit_par()
        return par

    def check_header(self, headers):
        """Validate elements of the header."""
        chk_dict = {}
        # chk_dict is 1-indexed!
        chk_dict[1] = {}
        # THIS CHECK IS A MUST! It performs a standard check to make sure the data are 2D.
        chk_dict[1]['NAXIS'] = 2
        return chk_dict

    def header_keys(self):
        """
        Header keys specific to keck_nirspec

        Returns:

        """
        head_keys = self.nirspec_header_keys()
        # Add the name of the filter used
        head_keys[0]['filter'] = 'FILNAME'
        return head_keys

    def setup_arcparam(self, arcparam, disperser=None, msarc_shape=None,
                       binspectral=None, **null_kwargs):
        """

        Args:
            arcparam:
            disperser:
            fitstbl:
            arc_idx:
            msarc_shape:
            **null_kwargs:

        Returns:

        """
        arcparam['lamps'] = ['OH_R24000']
        if fitstbl['filter'][arc_idx] == 'NIRSPEC-1':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.1093  # Ang per pixel for Low-Res, NIRSPEC-1 filter
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0]
            arcparam['wvmnx'][0] = 9400.  # Min wavelength
            arcparam['wvmnx'][1] = 11300.  # Max wavelength
            arcparam['wv_cen'] = 10000.  # Central wavelength

