# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload

from matplotlib import pyplot as plt

from pypit import msgs
from pypit import ardebug as debugger
from pypit import masterframe
from pypit import ararc
from pypit import ararclines
from pypit import armasters
from pypit.core import arsort

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

frametype = 'wv_calib'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
default_settings = dict(calibrate={'nfitpix': 5,
                                   'IDpixels': None, # User input pixel values
                                   'IDwaves': None,  # User input wavelength values
                                   'lamps': None,
                                   'method': 'arclines',
                                   'detection':  6.,
                                   'numsearch': 20,
                                   }
                        )
#settings_spect[dnum]['saturation']*settings_spect[dnum]['nonlinear'])  -- For satmask (echelle)

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
        if settings is None:
            self.settings = default_settings.copy()
        else:
            self.settings = settings.copy()
            if 'calibrate' not in self.settings.keys():
                self.settings.update(default_settings)
        self.spectrograph = spectrograph

        # Attributes
        self.frametype = frametype
        self.steps = []

        # Main outputs
        self.wv_calib = {}

        # Key Internals

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

    def _build_wv_calib(self, method, skip_QA=False):
        # Loop
        self.wv_calib = {}
        ok_mask = np.where(self.maskslits == 0)[0]
        for slit in ok_mask:
            ###############
            # Extract arc and identify lines
            if method == 'simple':
                iwv_calib = ararc.simple_calib(self.msarc, self.arcparam,
                                               self.arccen[:, slit],
                                               nfitpix=self.settings['calibrate']['nfitpix'],
                                               IDpixels=self.settings['calibrate']['IDpixels'],
                                               IDwaves=self.settings['calibrate']['IDwaves'])
            elif method == 'arclines':
                iwv_calib = ararc.calib_with_arclines(self.arcparam, self.arccen[:, slit])
            self.wv_calib[str(slit)] = iwv_calib.copy()
            # QA
            if not skip_QA:
                ararc.arc_fit_qa(self.setup, iwv_calib, slit)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.wv_calib

    def calibrate_spec(self, slit, method='arclines'):
        reload(ararc)
        spec = self.wv_calib[str(slit)]['spec']
        if method == 'simple':
            iwv_calib = ararc.simple_calib(self.det, self.msarc, self.arcparam,
                                           censpec=self.arccen[:, slit], slit=slit)
        elif method == 'arclines':
            iwv_calib = ararc.calib_with_arclines(self.arcparam, spec)
        return iwv_calib

    def _extract_arcs(self, lordloc, rordloc, pixlocn):
        self.arccen, self.maskslits, _ = ararc.get_censpec(lordloc, rordloc, pixlocn,
                                                          self.msarc, self.det, self.settings,
                                                          gen_satmask=False)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arccen, self.maskslits

    def load_wv_calib(self, filename):
        self.wv_calib, _, _ =  armasters._load(filename, frametype=self.frametype)
        # Recast a few items as arrays -- MIGHT PUSH THIS INTO armasters._load
        for key in self.wv_calib.keys():
            if key in ['steps', 'arcparam']:  # This isn't really necessary
                continue
            for tkey in self.wv_calib[key].keys():
                if tkey in ['tcent', 'spec', 'xfit', 'yfit', 'xrej']:
                    self.wv_calib[key][tkey] = np.array(self.wv_calib[key][tkey])
        # arcparam
        if 'arcparam' in self.wv_calib.keys():
            self.arcparam = self.wv_calib['arcparam'].copy()

    def _load_arcparam(self, calibrate_lamps=None):
        """

        Parameters
        ----------
        calibrate_lamps : str, optional
           List of lamps used

        Returns
        -------

        """
        # Setup arc parameters (e.g. linelist)
        arc_idx = arsort.ftype_indices(self.fitstbl, 'arc', self.sci_ID)
        self.arcparam = ararc.setup_param(self.spectrograph, self.msarc.shape,
                                          self.fitstbl, arc_idx[0],
                                          calibrate_lamps=calibrate_lamps)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arcparam

    def _make_maskslits(self, shape):
        # Set mask based on wv_calib
        mask = np.array([True]*shape)
        for key in self.wv_calib.keys():
            if key in ['steps', 'arcparam']:
                continue
            #
            mask[int(key)] = False
        self.maskslits = mask
        return self.maskslits

    def master(self, shape):
        """ Main driver for wavelength calibration

        Parameters
        ----------

        Returns
        -------
        self.wv_calib : dict or None
        self.maskslits : dict or None
        """
        # Attempt to load the Master Frame
        self.wv_calib, _, _ = self.load_master_frame()
        if self.wv_calib is None:
            return None, None
        # Mask
        self.maskslits = self._make_maskslits(shape)
        # Finish
        return self.wv_calib, self.maskslits


    def run(self, lordloc, rordloc, pixlocn, nonlinear=None, skip_QA=False):
        """
        Main driver for wavelength calibration

        Parameters
        ----------
        lordloc
        rordloc
        pixlocn
        nonlinear
        skip_QA

        Returns
        -------

        """
        ###############
        # Extract an arc down each slit
        #   The settings here are settings.spect (saturation and nonlinear)
        _, _ = self._extract_arcs(lordloc, rordloc, pixlocn)

        # Load the arcparam
        _ = self._load_arcparam()

        # Fill up the calibrations and generate QA
        self.wv_calib = self._build_wv_calib(self.settings['calibrate']['method'], skip_QA=skip_QA)
        self.wv_calib['steps'] = self.steps
        sv_aparam = self.arcparam.copy()
        sv_aparam.pop('llist')
        self.wv_calib['arcparam'] = sv_aparam

        # Build mask
        self._make_maskslits(lordloc.shape[1])
        # Return
        return self.wv_calib, self.maskslits

    def show(self, item, slit=None):
        if item == 'spec':
            # spec
            spec = self.wv_calib[str(slit)]['spec']
            # tcent
            tcent = self.wv_calib[str(slit)]['tcent']
            yt = np.zeros_like(tcent)
            for jj,t in enumerate(tcent):
                it = int(np.round(t))
                yt[jj] = np.max(spec[it-1:it+1])
            # Plot
            plt.clf()
            ax=plt.gca()
            ax.plot(spec, drawstyle='steps-mid')
            ax.scatter(tcent, yt, color='red', marker='*')
            ax.set_xlabel('Pixel')
            ax.set_ylabel('Counts')
            plt.show()
        elif item == 'fit':
            ararc.arc_fit_qa(None, self.wv_calib[str(slit)], slit, outfile='show')

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: >'.format(self.__class__.__name__)
        return txt



