# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import os
import inspect
import numpy as np

#from importlib import reload

from matplotlib import pyplot as plt

from pypeit import msgs
from pypeit import masterframe
from pypeit.core import arc
from pypeit.core import masters
from pypeit.core import fsort
from pypeit.par import pypeitpar

from pypeit.spectrographs.spectrograph import Spectrograph
from pypeit.spectrographs.util import load_spectrograph

from pypeit import debugger



class WaveCalib(masterframe.MasterFrame):
    """Class to guide slit/order tracing

    .. todo::
        - Need to clean up these docs.

    Parameters
    ----------
    msarc : ndarray
      Arc image, created by the ArcImage class

    Optional Parameters
    --------------------
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

    # Frametype is a class attribute
    frametype = 'wv_calib'

    # ToDo This code will crash is spectrograph and det are not set. I see no reason why these should be optional
    # parameters since instantiating without them does nothing. Make them required
    def __init__(self, msarc, spectrograph=None, par=None, det=None, setup=None, master_dir=None,
                 mode=None, fitstbl=None, sci_ID=None, arcparam=None, redux_path=None):

        # Instantiate the spectograph
        if isinstance(spectrograph, str):
            self.spectrograph = load_spectrograph(spectrograph=spectrograph)
        elif isinstance(spectrograph, Spectrograph):
            self.spectrograph = spectrograph
        else:
            raise TypeError('Must provide a name or instance for the Spectrograph.')

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         master_dir=master_dir, mode=mode)

        # Required parameters (but can be None)
        self.msarc = msarc

        self.par = pypeitpar.WavelengthSolutionPar() if par is None else par

        # Optional parameters
        self.redux_path = redux_path
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID
        self.det = det
        self.setup = setup
        self.arcparam = arcparam

        # Attributes
        # Done by MasterFrame
        self.steps = []

        # Main outputs
        self.wv_calib = {}

        # Key Internals
        self.arccen = None


    def _build_wv_calib(self, method, skip_QA=False):
        """
        Main routine to generate the wavelength solutions in a loop over slits

        Wrapper to arc.simple_calib or arc.calib_with_arclines

        Parameters
        ----------
        method : str
          'simple' -- arc.simple_calib
          'arclines' -- arc.calib_with_arclines
        skip_QA : bool, optional

        Returns
        -------
        self.wv_calib : dict

        """
        # Obtain a list of good slits
        ok_mask = np.where(self.maskslits == 0)[0]

        # Obtain calibration for all slits
        if method == 'simple':
            self.wv_calib = arc.simple_calib_driver(self.msarc, self.arcparam, self.arccen, ok_mask,
                                                    nfitpix=self.par['nfitpix'],
                                                    IDpixels=self.par['IDpixels'],
                                                    IDwaves=self.par['IDwaves'])
        elif method == 'arclines':
            self.wv_calib = arc.calib_with_arclines(self.arcparam, self.arccen, ok_mask=ok_mask)

        # QA
        if not skip_QA:
            for slit in ok_mask:
                arc.arc_fit_qa(self.setup, self.wv_calib[str(slit)], slit, out_dir=self.redux_path)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.wv_calib

    def calibrate_spec(self, slit, method='arclines'):
        """
        TODO: Deprecate this function? It's only being used by the tests
        User method to calibrate a given spectrum from a chosen slit

        Wrapper to arc.simple_calib or arc.calib_with_arclines

        Parameters
        ----------
        slit : int
        method : str, optional
          'simple' -- arc.simple_calib
          'arclines' -- arc.calib_with_arclines

        Returns
        -------
        iwv_calib : dict
          Solution for that single slit

        """
        spec = self.wv_calib[str(slit)]['spec']
        if method == 'simple':
            iwv_calib = arc.simple_calib(self.msarc, self.arcparam, self.arccen[:, slit])
        elif method == 'arclines':
            iwv_calib = arc.calib_with_arclines(self.arcparam, spec.reshape((spec.size, 1)))
        else:
            msgs.error("Not an allowed method")
        return iwv_calib

    def _extract_arcs(self, lordloc, rordloc, pixlocn):
        """
        Extract an arc down the center of each slit/order

        Wrapper to arc.get_censpec

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

        # ToDO JFH As far as I can tell, passing in nonlinear_counts with gen_satmask = False does absolutely nothing.
        # Not sure what was intended here.
        nonlinear_counts = self.spectrograph.detector[self.det - 1]['saturation'] * \
                          self.spectrograph.detector[self.det - 1]['nonlinear']
        self.arccen, self.maskslits, _ \
                    = arc.get_censpec(lordloc, rordloc, pixlocn, self.msarc, self.det,
                                        nonlinear_counts=nonlinear_counts, gen_satmask=False)
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
        self.wv_calib, _, _ =  masters._load(filename, frametype=self.frametype)
        # Recast a few items as arrays
        # TODO -- Consider pushing into master
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

        Wrapper to arc.setup_param

        Parameters
        ----------
        calibrate_lamps : str, optional
           List of lamps used

        Returns
        -------

        """
        # Setup arc parameters (e.g. linelist)
        arc_idx = fsort.ftype_indices(self.fitstbl, 'arc', self.sci_ID)
        self.arcparam = arc.setup_param(self.spectrograph, self.msarc.shape, self.fitstbl,
                                          arc_idx[0], calibrate_lamps=calibrate_lamps)
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
          Would be passed to arc.detect_lines but that routine is
          currently being run in arclines.holy
        skip_QA : bool, optional

        Returns
        -------
        self.wv_calib
        self.maskslits

        """
        ###############
        # Extract an arc down each slit
        _, _ = self._extract_arcs(lordloc, rordloc, pixlocn)

        # Load the arcparam, if one was not passed in.
        if self.arcparam is None:
            _ = self._load_arcparam()
        # This call is unfortunate since it requires the fitstable in
        # the true PYPIT style of bloated overloaded and uncenessary
        # argument lists.  It is mainly here for backwards compatibility
        # with old methods for wavelength calibration that have been
        # superceded by holy grail. Holy grail only requires a linelist
        # in the arcparam dict, so if the user passes in an arcparam, no
        # need to run this.
        # KBW - What's stopping us from fixing this?

        # Fill up the calibrations and generate QA
        self.wv_calib = self._build_wv_calib(self.par['method'], skip_QA=skip_QA)
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
            arc.arc_fit_qa(None, self.wv_calib[str(slit)], slit, outfile='show')

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


