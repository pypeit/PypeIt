# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from matplotlib import pyplot as plt

from pypit import msgs
from pypit import ardebug as debugger
from pypit import masterframe
from pypit.core import ararc
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
    msarc : ndarray
      Arc image
    settings : dict, optional
      Settings for trace slits
    det : int, optional
      Detector number
    setup : str, optional
    fitstbl : Table, optional
      Used for arcparam
    sci_ID : int, optional
      Index of the science frame (also for arcparam)

    Attributes
    ----------
    frametype : str
      Hard-coded to 'wv_calib'
    steps : list
      List of the processing steps performed
    wv_calib : dict
      Primary output
      Keys
        0, 1, 2, 3 -- Solution for individual slit
        arcparam -- Parameters used
        steps
    arccen : ndarray (nwave, nslit)
      Extracted arc(s) down the center of the slit(s)
    maskslit : ndarray (nslit)
      Slits to ignore because they were not extacted
    arcparam : dict
      Arc parameter (instrument/disperser specific)
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
        self.arccen = None

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup, self.settings)

    def _build_wv_calib(self, method, skip_QA=False):
        """
        Main routine to generate the wavelength solutions in a loop over slits

        Wrapper to ararc.simple_calib or ararc.calib_with_arclines

        Parameters
        ----------
        method : str
          'simple' -- ararc.simple_calib
          'arclines' -- ararc.calib_with_arclines
        skip_QA : bool, optional

        Returns
        -------
        self.wv_calib : dict

        """
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
        """
        User method to calibrate a given spectrum from a chosen slit

        Wrapper to ararc.simple_calib or ararc.calib_with_arclines

        Parameters
        ----------
        slit : int
        method : str, optional
          'simple' -- ararc.simple_calib
          'arclines' -- ararc.calib_with_arclines

        Returns
        -------
        iwv_calib : dict
          Solution for that single slit

        """
        spec = self.wv_calib[str(slit)]['spec']
        if method == 'simple':
            iwv_calib = ararc.simple_calib(self.msarc, self.arcparam,
                                           self.arccen[:, slit])
        elif method == 'arclines':
            iwv_calib = ararc.calib_with_arclines(self.arcparam, spec)
        else:
            msgs.error("Not an allowed method")
        return iwv_calib

    def _extract_arcs(self, lordloc, rordloc, pixlocn):
        """
        Extract an arc down the center of each slit/order

        Wrapper to ararc.get_censpec

        Parameters
        ----------
        lordloc : ndarray
          Left edges (from TraceSlit)
        rordloc : ndarray
          Right edges (from TraceSlit)
        pixlocn : ndarray

        Returns
        -------
        self.arccen
          1D arc spectra from each slit
        self.maskslits

        """
        self.arccen, self.maskslits, _ = ararc.get_censpec(lordloc, rordloc, pixlocn,
                                                          self.msarc, self.det, self.settings,
                                                          gen_satmask=False)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arccen, self.maskslits

    def load_wv_calib(self, filename):
        """
        Load a full (all slit) wv_calib dict

        Includes converting the JSON lists of particular items into ndarray

        Parameters
        ----------
        filename : str

        Returns
        -------
        Fills:
          self.wv_calib
          self.arcparam

        """
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
        Load the arc parameters

        Wrapper to ararc.setup_param

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

    def _make_maskslits(self, nslit):
        """

        Parameters
        ----------
        nslit : int
          Number of slits/orders

        Returns
        -------
        self.maskslits : ndarray (bool)

        """
        # Set mask based on wv_calib
        mask = np.array([True]*nslit)
        for key in self.wv_calib.keys():
            if key in ['steps', 'arcparam']:
                continue
            #
            mask[int(key)] = False
        self.maskslits = mask
        return self.maskslits

    def run(self, lordloc, rordloc, pixlocn, nonlinear=None, skip_QA=False):
        """
        Main driver for wavelength calibration

        Code flow:
          1. Extract 1D arc spectra down the center of each slit/order
          2. Load the parameters guiding wavelength calibration
          3. Generate the 1D wavelength fits
          4. Generate a mask

        Parameters
        ----------
        lordloc : ndarray
          From a TraceSlit object
        rordloc : ndarray
          From a TraceSlit object
        pixlocn : ndarray
          From a TraceSlit object
        nonlinear : float, optional
          Would be passed to ararc.detect_lines but that routine is
          currently being run in arclines.holy
        skip_QA : bool, optional

        Returns
        -------
        self.wv_calib
        self.maskslits

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
        """
        Show one of the class internals

        Parameters
        ----------
        item : str
          'spec' -- Show the fitted points and solution;  requires slit
          'fit' -- Show fit QA; requires slit
        slit : int, optional

        Returns
        -------

        """
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
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt



