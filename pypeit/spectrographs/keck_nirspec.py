""" Module for Keck/NIRSPEC specific codes
"""
from __future__ import absolute_import, division, print_function


import numpy as np

from pypeit import msgs
from pypeit.par.pypeitpar import DetectorPar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes

from pypeit import debugger

class KeckNIRSPECSpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRSPEC specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRSPECSpectrograph, self).__init__()
        self.spectrograph = 'keck_nirspec'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRSPEC'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
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
        par['rdx']['spectrograph'] = 'keck_nirspec'
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
        def_keys[0]['lampname01'] = 'NEON'
        def_keys[0]['lampstat01'] = '01.NEON'
        def_keys[0]['lampname02'] = 'ARGON'
        def_keys[0]['lampstat02'] = '01.ARGON'
        def_keys[0]['lampname03'] = 'KRYPTON'
        def_keys[0]['lampstat03'] = '01.KRYPTON'
        def_keys[0]['lampname04'] = 'XENON'
        def_keys[0]['lampstat04'] = '01.XENON'
        def_keys[0]['lampname05'] = 'ETALON'
        def_keys[0]['lampstat05'] = '01.ETALON'
        def_keys[0]['lampname06'] = 'FLAT'
        def_keys[0]['lampstat06'] = '01.FLAT'


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
        self.spectrograph = 'keck_nirspec'


    def default_pypeit_par(self):
        """
        Set default parameters for NIRSPEC low-dispersion reductions
        """
        par = self.nirspec_default_pypeit_par()
        par['rdx']['pipeline'] = ARMS
        return par


    def header_keys(self):
        """
        Header keys specific to keck_nirspec low-dispersion

        Returns:

        """
        head_keys = self.nirspec_header_keys()
        # Add the name of the dispersing element
        head_keys[0]['dispname'] = '01.DISPERS'
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
        if fitstbl['filter1'][arc_idx] == 'NIRSPEC-1':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.1093  # Ang per pixel for Low-Res, NIRSPEC-1 filter
            arcparam['b1'] = 1. / arcparam['disp'] / msarc_shape[0]
            arcparam['wvmnx'][0] = 9400.  # Min wavelength
            arcparam['wvmnx'][1] = 11300.  # Max wavelength
            arcparam['wv_cen'] = 10000.  # Central wavelength

