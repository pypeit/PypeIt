# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload


from pypit import msgs
from pypit import ardebug as debugger
from pypit import masterframe
from pypit import ararc
from pypit import ararclines
from pypit.core import arsort

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

frametype = 'wv_calib'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.

#  See save_master() for the data model for output


class WaveCalib(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    Parameters
    ----------
    mstrace : ndarray
      Trace image
    pixlocn : ndarray
      Pixel location array
    binbpx : ndarray, optional
      Bad pixel mask
      If not provided, a dummy array with no masking is generated
    settings : dict, optional
      Settings for trace slits
    det : int, optional
      Detector number
    ednum : int, optional
      Edge number used for indexing

    Attributes
    ----------
    frametype : str
      Hard-coded to 'wv_calib'

    steps : list
      List of the processing steps performed
    """
    def __init__(self, msarc, spectrograph=None, settings=None, det=None, setup=None, fitstbl=None, sci_ID=None):

        # Required parameters (but can be None)
        self.msarc = msarc

        # Optional parameters
        self.det = det
        self.fitstbl = fitstbl
        self.setup = setup
        self.sci_ID = sci_ID
        self.settings = settings
        self.spectrograph = spectrograph

        # Attributes
        self.frametype = frametype
        self.steps = []

        # Main outputs
        self.wv_calib = {}

        # Key Internals

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

    def _load_arcparam(self, calibrate_lamps=None):
        """

        Parameters
        ----------
        calibrate_lamps : str, optional
           List of lamps used

        Returns
        -------

        """
        reload(ararc)
        reload(ararclines)
        # Setup arc parameters (e.g. linelist)
        arc_idx = arsort.ftype_indices(self.fitstbl, 'arc', self.sci_ID)
        self.arcparam = ararc.setup_param(self.spectrograph, self.msarc.shape, self.fitstbl, arc_idx[0],
                                          calibrate_lamps=calibrate_lamps)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arcparam


    def _extract_arcs(self, lordloc, rordloc, pixlocn):
        self.arccen, self.maskslit, _ = ararc.get_censpec(lordloc, rordloc, pixlocn,
                                                    self.msarc, self.det, self.settings,
                                                    gen_satmask=False)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arccen, self.maskslit

    def master(self, method, lordloc, rordloc, pixlocn, nonlinear=None):
        """ Main driver for wavelength calibration

        Parameters
        ----------

        Returns
        -------
        """
        # Attempt to load the Master Frame
        wv_calib, _, _ = self.load_master_frame(self, "wv_calib")
        if wv_calib is None:

            ###############
            # Extract an arc down each slit
            #   The settings here are settings.spect (saturation and nonlinear)
            _, _ = self._extract_arcs(lordloc, rordloc, pixlocn)

            # Load the arcparam
            _ = self._load_arcparam()


            #
            ok_mask = np.where(maskslit == 0)[0]

            # Fill up the calibrations
            wv_calib = {}
            for slit in ok_mask:
                ###############
                # Extract arc and identify lines
                #if settings.argflag['arc']['calibrate']['method'] == 'simple':
                #elif settings.argflag['arc']['calibrate']['method'] == 'arclines':
                if method == 'simple':
                    iwv_calib = ararc.simple_calib(self.det, self.msarc, arcparam, censpec=arccen[:, slit], slit=slit)
                elif method == 'arclines':
                    iwv_calib = ararc.calib_with_arclines(self.det, self.msarc, slit, arcparam, censpec=arccen[:, slit])
                wv_calib[str(slit)] = iwv_calib.copy()

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



