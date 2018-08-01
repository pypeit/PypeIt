""" Module for Keck/NIRES specific codes
"""
from __future__ import absolute_import, division, print_function

import numpy as np

from pypeit import msgs
from pypeit.par.pypeitpar import DetectorPar
from pypeit.spectrographs import spectrograph
from pypeit import telescopes

from pypeit import debugger

# TODO: KeckNIRESSPectrograph instead (i.e. two SS) ?
class KeckNIRESpectrograph(spectrograph.Spectrograph):
    """
    Child to handle Keck/NIRES specific code
    """
    def __init__(self):
        # Get it started
        super(KeckNIRESpectrograph, self).__init__()
        self.spectrograph = 'keck_nires'
        self.telescope = telescopes.KeckTelescopePar()
        self.camera = 'NIRES'
        self.detector = [
                # Detector 1
                DetectorPar(dataext         = 0,
                            dispaxis        = 1,
                            xgap            = 0.,
                            ygap            = 0.,
                            ysize           = 1.,
                            platescale      = 0.123,
                            darkcurr        = 0.01,
                            saturation      = 65535.,
                            nonlinear       = 0.76,
                            numamplifiers   = 1,
                            gain            = 3.8,
                            ronoise         = 2,
                            datasec         = '[:,:]',
                            oscansec        = '[:,:]'
                            )
            ]
        # Uses default timeunit
        # Uses default primary_hdrext
        # self.sky_file = ?

    def bpm(self, shape=None, **null_kwargs):
        """ Generate a BPM

        Parameters
        ----------
        shape : tuple, REQUIRED

        Returns
        -------
        badpix : ndarray

        """
        # Edges of the detector are junk
        msgs.info("Custom bad pixel mask for NIRES")
        if shape is None and self.naxis is None:
            raise ValueError('Must define shape!')
        _shape = self.naxis if shape is None else shape
        self.bpm = np.zeros(_shape)
        #self.bpm[:, :20] = 1.
        #self.bpm[:, 1000:] = 1.

    def setup_arcparam(self, arcparam, fitstbl=None, arc_idx=None,
                       msarc_shape=None, **null_kwargs):
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
        arcparam = dict(llist='',
                        disp=2.,  # Ang/unbinned pixel
                        b1=0.,  # Pixel fit term (binning independent)
                        b2=0.,  # Pixel fit term
                        lamps=['OH_triplespec'],  # Line lamps on
                        wv_cen=0.,  # Estimate of central wavelength
                        wvmnx=[9000., 25000.],  # Guess at wavelength range
                        disp_toler=0.1,  # 10% tolerance
                        match_toler=3.,  # Matching tolerance (pixels)
                        min_ampl=1000.,  # Minimum amplitude
                        func='legendre',  # Function for fitting
                        n_first=1,  # Order of polynomial for first fit
                        n_final=3,  # Order of polynomial for final fit
                        nsig_rej=2.,  # Number of sigma for rejection
                        nsig_rej_final=2.0,  # Number of sigma for rejection (final fit)
                        Nstrong=13)  # Number of lines for auto-analysis

