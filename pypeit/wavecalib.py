# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import os
import inspect
import numpy as np

#from importlib import reload

from matplotlib import pyplot as plt

from astropy.table import vstack

from pypeit import msgs
from pypeit import masterframe
from pypeit.core import arc
from pypeit.core import qa
from pypeit.core import masters
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
from pypeit.core.wavecal import autoid
from pypeit.core.wavecal import waveio

from pypeit import debugger

# Frametype is a class attribute
frametype = 'wv_calib'

# By moving this utility routine out of the class, we do not need to instantiate the class to read in the wv_calib dict
# from a file
def load_wv_calib(filename):
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
    wv_calib, _, _ = masters._load(filename, frametype=frametype)
    # Recast a few items as arrays
    # TODO -- Consider pushing into master
    for key in wv_calib.keys():
        if key in ['steps', 'par']:  # This isn't really necessary
            continue
        for tkey in wv_calib[key].keys():
            if tkey in ['tcent', 'spec', 'xfit', 'yfit', 'xrej']:
                wv_calib[key][tkey] = np.array(wv_calib[key][tkey])
    # parset
    if 'par' in wv_calib.keys():
        par = wv_calib['par'].copy()
    else:
        par = None

    return (wv_calib, par)

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

    # ToDo This code will crash is spectrograph and det are not set. I see no reason why these should be optional
    # parameters since instantiating without them does nothing. Make them required
    def __init__(self, msarc, spectrograph=None, par=None, det=None, setup=None, master_dir=None,
                 mode=None, fitstbl=None, sci_ID=None, redux_path=None, bpm = None):

        # Instantiate the spectograph
        self.spectrograph = load_spectrograph(spectrograph)

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         master_dir=master_dir, mode=mode)

        # Required parameters (but can be None)
        self.msarc = msarc
        self.bpm = bpm

        self.par = pypeitpar.WavelengthSolutionPar() if par is None else par

        # Optional parameters
        self.redux_path = redux_path
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID
        self.det = det
        self.setup = setup
        #self.arcparam = arcparam

        # Attributes
        # Done by MasterFrame
        self.steps = []

        # Main outputs
        self.wv_calib = {}

        # Key Internals
        self.arccen = None


    def _build_wv_calib(self, method, skip_QA=False, use_method='general'):
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
            # Should only run this on 1 slit
            #self.par['n_first'] = 2
            #self.par['n_final'] = 3
            #self.par['func'] = 'legendre'
            #self.par['sigrej_first'] = 2.
            #self.par['sigrej_final'] = 3.
            #self.par['match_toler'] = 3.
            #self.arcparam['Nstrong'] = 13

            CuI = waveio.load_line_list('CuI', use_ion=True, NIST=True)
            ArI = waveio.load_line_list('ArI', use_ion=True, NIST=True)
            ArII = waveio.load_line_list('ArII', use_ion=True, NIST=True)
            llist = vstack([CuI, ArI, ArII])
            self.arcparam['llist'] = llist

            self.wv_calib = arc.simple_calib_driver(self.msarc, self.par, self.arccen, ok_mask,
                                                    nfitpix=self.par['nfitpix'],
                                                    IDpixels=self.par['IDpixels'],
                                                    IDwaves=self.par['IDwaves'])
        elif method == 'arclines':
            if ok_mask is None:
                ok_mask = np.arange(self.arccen.shape[1])

            if use_method == "semi-brute":
                debugger.set_trace()  # THIS IS BROKEN
                final_fit = {}
                for slit in ok_mask:
                    # HACKS BY JXP
                    self.par['wv_cen'] = 8670.
                    self.par['disp'] = 1.524
                    # ToDO remove these hacks and use the parset in semi_brute
                    best_dict, ifinal_fit = autoid.semi_brute(self.arccen[:, slit],
                                                              self.par['lamps'], self.par['wv_cen'],
                                                              self.par['disp'],
                                                              match_toler=self.par['match_toler'], func=self.par['func'],
                                                              n_first=self.par['n_first'],sigrej_first=self.par['n_first'],
                                                              n_final=self.par['n_final'], sigrej_final=self.par['sigrej_final'],
                                                              min_nsig=self.par['min_nsig'],
                                                              nonlinear_counts= self.par['nonlinear_counts'])
                    final_fit[str(slit)] = ifinal_fit.copy()
            elif use_method == "basic":
                final_fit = {}
                for slit in ok_mask:
                    status, ngd_match, match_idx, scores, ifinal_fit = \
                        autoid.basic(self.arccen[:, slit], self.par['lamps'], self.par['wv_cen'], self.par['disp'],
                                     nonlinear_counts = self.par['nonlinear_counts'])
                    final_fit[str(slit)] = ifinal_fit.copy()
            else:
                # Now preferred
                arcfitter = autoid.General(self.arccen, par = self.par, ok_mask=ok_mask)
                patt_dict, final_fit = arcfitter.get_results()
            self.wv_calib = final_fit

        # QA
        if not skip_QA:
            for slit in ok_mask:
                outfile = qa.set_qa_filename(self.setup, 'arc_fit_qa', slit=(slit + 1), out_dir=self.redux_path)
                arc.arc_fit_qa(self.wv_calib[str(slit)], outfile)
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
            iwv_calib = arc.simple_calib(self.msarc, self.par, self.arccen[:, slit])
        elif method == 'arclines':
            iwv_calib = arc.calib_with_arclines(self.par, spec.reshape((spec.size, 1)))
        else:
            msgs.error("Not an allowed method")
        return iwv_calib

    def _extract_arcs(self, lordloc, rordloc, slitpix):
        """
        Extract an arc down the center of each slit/order

        Wrapper to arc.get_censpec

        Parameters
        ----------
        lordloc : ndarray
          Left edges (from TraceSlit)
        rordloc : ndarray
          Right edges (from TraceSlit)
        slitpix : ndarray

        Returns
        -------
        self.arccen
          1D arc spectra from each slit
        self.maskslits

        """
        inmask = (self.bpm == 0) if self.bpm is not None else None
        self.arccen, self.maskslits = arc.get_censpec(lordloc, rordloc, slitpix, self.msarc,
                                                      inmask=inmask)

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.arccen, self.maskslits

    def load_wv_calib_old(self, filename):
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
            if key in ['steps', 'par']:  # This isn't really necessary
                continue
            for tkey in self.wv_calib[key].keys():
                if tkey in ['tcent', 'spec', 'xfit', 'yfit', 'xrej']:
                    self.wv_calib[key][tkey] = np.array(self.wv_calib[key][tkey])
        # parset
        if 'par' in self.wv_calib.keys():
            self.par = self.wv_calib['par'].copy()



    def master(self):
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

        # This is code from masters._load. I hate the fact that the master loading is in a separate master._load method in
        # a different module and am moving it here.
        name = self.ms_name()
        # Check to see if file exists
        if not os.path.isfile(name):
            msgs.warn("Master frame does not exist: {:s}".format(name))
            if force:
                msgs.error("Crashing out because reduce-masters-force=True:" + msgs.newline() + name)
            return None, None, None

        msgs.info("Loading Master {0:s} frame:".format(self.frametype)+msgs.newline()+name)
        self.wv_calib = linetools.utils.loadjson(name)

        for key in self.wv_calib.keys():
            if key in ['steps', 'par']:  # This isn't really necessary
                continue
            for tkey in self.wv_calib[key].keys():
                if tkey in ['tcent', 'spec', 'xfit', 'yfit', 'xrej']:
                    self.wv_calib[key][tkey] = np.array(self.wv_calib[key][tkey])
        # parset
        if 'par' in self.wv_calib.keys():
            self.par = self.wv_calib['par'].copy()



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
        arc_idx = self.fitstbl.find_frames('arc', sci_ID=self.sci_ID, index=True)
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
            if key in ['steps', 'par']:
                continue
            #
            mask[int(key)] = False
        self.maskslits = mask
        return self.maskslits

    def run(self, lordloc, rordloc, slitpix, nonlinear=None, skip_QA=False):
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
        slitpix : ndarray
          slitmask from tslits_dict
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
        _, _ = self._extract_arcs(lordloc, rordloc, slitpix)


        #if self.arcparam is None:
        #    _ = self._load_arcparam()

        # Fill up the calibrations and generate QA
        self.wv_calib = self._build_wv_calib(self.par['method'], skip_QA=skip_QA)
        self.wv_calib['steps'] = self.steps
        sv_par = self.par.data.copy()
        #sv_par.pop('llist')
        self.wv_calib['par'] = sv_par

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


