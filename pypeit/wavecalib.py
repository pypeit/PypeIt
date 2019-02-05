# Module for guiding 1D Wavelength Calibration
from __future__ import absolute_import, division, print_function

import os
import inspect
import numpy as np

#from importlib import reload

from matplotlib import pyplot as plt

from astropy.table import vstack

import copy
from pypeit import msgs
from pypeit import masterframe
from pypeit.core import arc, wavecal, qa, pixels
from pypeit.par import pypeitpar
from pypeit.spectrographs.util import load_spectrograph
import linetools.utils


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
    binspectral: int, optional
      Binning of the Arc in the spectral dimension

    det : int, optional
      Detector number
    master_key : str, optional
    spectrograph : Spectrograph
      Optional so that one can load a solution without using the Spectrograph, e.g. waveio.load_reid_arxiv()

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
    maskslits : ndarray (nslit); bool
      Slits to ignore because they were not extracted
      WARNING: Outside of this Class, it is best to regenerate
      the mask using  make_maskslits()
    arcparam : dict
      Arc parameter (instrument/disperser specific)
    """

    # Frametype is a class attribute
    frametype = 'wv_calib'

    def __init__(self, msarc, tslits_dict, binspectral=None, spectrograph=None, par=None, det=1, master_key=None, master_dir=None,
                 reuse_masters=False, redux_path=None, bpm=None):

        # Instantiate the spectograph
        self.spectrograph = load_spectrograph(spectrograph)

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.frametype, master_key,
                                         master_dir=master_dir, reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.msarc = msarc
        self.tslits_dict = tslits_dict

        # Optional parameters
        self.bpm = bpm
        self.par = pypeitpar.WavelengthSolutionPar() if par is None else par
        self.binspectral = binspectral
        self.redux_path = redux_path
        self.det = det
        self.master_key = master_key

        # Attributes
        self.steps = []    # steps executed
        self.wv_calib = {} # main output
        self.arccen = None # central arc spectrum

        # TODO this code is duplicated verbatim in wavetilts. Should it be a function
        if spectrograph is not None:
            if self.spectrograph.detector is not None:
                self.nonlinear_counts = self.spectrograph.detector[self.det-1]['saturation']*self.spectrograph.detector[self.det-1]['nonlinear']
            else:
                self.nonlinear_counts=1e10
        else:
            self.nonlinear_counts=1e10
        # Set the slitmask and slit boundary related attributes that the code needs for execution. This also deals with
        # arcimages that have a different binning then the trace images used to defined the slits
        if self.tslits_dict is not None and self.msarc is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            inmask = (self.bpm == 0) if self.bpm is not None else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slit_righ'])
            self.slitcen   = arc.resize_slits2arc(self.shape_arc, self.shape_science, self.tslits_dict['slitcen'])
            self.slitmask  = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.inmask = arc.resize_mask2arc(self.shape_arc,inmask)
            # TODO -- Remove the following two lines if deemed ok
            if self.par['method'] != 'full_template':
                self.inmask &= self.msarc < self.nonlinear_counts
        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.nslits = 0
            self.slit_left = None
            self.slit_righ = None
            self.slitcen = None
            self.slitmask = None
            self.inmask = None


    def build_wv_calib(self, arccen, method, skip_QA=False):
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
        ok_mask = np.where(~self.maskslits)[0]

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

            #CuI = wavecal.waveio.load_line_list('CuI', use_ion=True, NIST=True)
            #ArI = wavecal.waveio.load_line_list('ArI', use_ion=True, NIST=True)
            #ArII = wavecal.waveio.load_line_list('ArII', use_ion=True, NIST=True)
            #llist = vstack([CuI, ArI, ArII])
            #self.arcparam['llist'] = llist
            lines = self.par['lamps']
            line_lists = wavecal.waveio.load_line_lists(lines)

            self.wv_calib = arc.simple_calib_driver(self.msarc, line_lists, arccen, ok_mask,
                                                    nfitpix=self.par['nfitpix'],

                                                    IDpixels=self.par['IDpixels'],
                                                    IDwaves=self.par['IDwaves'])
        elif method == 'semi-brute':
            debugger.set_trace()  # THIS IS BROKEN
            final_fit = {}
            for slit in ok_mask:
                # HACKS BY JXP
                self.par['wv_cen'] = 8670.
                self.par['disp'] = 1.524
                # ToDO remove these hacks and use the parset in semi_brute
                best_dict, ifinal_fit = wavecal.autoid.semi_brute(arccen[:, slit],
                                                                  self.par['lamps'], self.par['wv_cen'],
                                                                  (self)['disp'],match_toler=self.par['match_toler'],
                                                                  func=self.par['func'],n_first=self.par['n_first'],
                                                                  sigrej_first=self.par['n_first'],
                                                                  n_final=self.par['n_final'],
                                                                  sigrej_final=self.par['sigrej_final'],
                                                                  sigdetect=self.par['sigdetect'],
                                                                  nonlinear_counts= self.nonlinear_counts)
                final_fit[str(slit)] = ifinal_fit.copy()
        elif method == 'basic':
            final_fit = {}
            for slit in ok_mask:
                status, ngd_match, match_idx, scores, ifinal_fit = \
                    wavecal.autoid.basic(arccen[:, slit], self.par['lamps'], self.par['wv_cen'], self.par['disp'],
                                 nonlinear_counts = self.nonlinear_counts)
                final_fit[str(slit)] = ifinal_fit.copy()
                if status != 1:
                    self.maskslits[slit] = True
        elif method == 'holy-grail':
            # Sometimes works, sometimes fails
            arcfitter = wavecal.autoid.HolyGrail(arccen, par=self.par, ok_mask=ok_mask)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'reidentify':
            # Now preferred
            arcfitter = wavecal.autoid.ArchiveReid(arccen, par=self.par, ok_mask=ok_mask)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'full_template':
            # Now preferred
            if self.binspectral is None:
                msgs.error("You must specify binspectral for the full_template method!")
            final_fit = wavecal.autoid.full_template(arccen, self.par, ok_mask, self.det,
                                                     self.binspectral,
                                                     nsnippet=self.par['nsnippet'])

        else:
            msgs.error('Unrecognized wavelength calibration method: {:}'.format(method))

        self.wv_calib = final_fit

        # Remake mask (*mainly for the QA that follows*)
        self.maskslits = self.make_maskslits(len(self.maskslits))
        ok_mask = np.where(~self.maskslits)[0]

        # QA
        if not skip_QA:
            for slit in ok_mask:
                outfile = qa.set_qa_filename(self.master_key, 'arc_fit_qa', slit=slit, out_dir=self.redux_path)
                wavecal.qa.arc_fit_qa(self.wv_calib[str(slit)], outfile = outfile)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.wv_calib


    def echelle_2dfit(self, wv_calib,debug=False, skip_QA = False):
        """
        Evaluate 2-d wavelength solution for echelle data. Unpacks wv_calib for slits to be input into  arc.fit2darc

        Parameters
        ----------
        wv_calib: dict
           Wavelength calibration

        Optional Parameters
        -------------------
        debug: bool, default = False
           Show debugging info
        skip_QA: bool, default = False
          Not yet implemented

        Returns
        -------
        fit2d_dict: dict
           dictionary containing information from 2-d fit
        """

        msgs.info('Fitting 2-d wavelength solution for echelle....')
        all_wave = np.array([], dtype=float)
        all_pixel = np.array([], dtype=float)
        all_order = np.array([],dtype=float)

        # Obtain a list of good slits
        ok_mask = np.where(~self.maskslits)[0]
        nspec = self.msarc.shape[0]
        for islit in wv_calib.keys():
            if int(islit) not in ok_mask:
                continue
            iorder = self.spectrograph.slit2order(islit)
            mask_now = wv_calib[islit]['mask']
            all_wave = np.append(all_wave, wv_calib[islit]['wave_fit'][mask_now])
            all_pixel = np.append(all_pixel, wv_calib[islit]['pixel_fit'][mask_now])
            all_order = np.append(all_order, np.full_like(wv_calib[islit]['pixel_fit'][mask_now], float(iorder)))

        fit2d_dict = arc.fit2darc(all_wave, all_pixel, all_order, nspec, nspec_coeff=self.par['ech_nspec_coeff'],
                                  norder_coeff=self.par['ech_norder_coeff'],sigrej=self.par['ech_sigrej'],
                                  debug=debug)

        self.steps.append(inspect.stack()[0][3])

        # QA
        if not skip_QA:
            outfile_global = qa.set_qa_filename(self.master_key, 'arc_fit2d_global_qa', out_dir=self.redux_path)
            arc.fit2darc_global_qa(fit2d_dict, outfile=outfile_global)
            outfile_orders = qa.set_qa_filename(self.master_key, 'arc_fit2d_orders_qa', out_dir=self.redux_path)
            arc.fit2darc_orders_qa(fit2d_dict, outfile=outfile_orders)

        return fit2d_dict


    # JFH TODO this method is identical to the code in wavetilts. SHould we make it a separate function?
    def extract_arcs(self, slitcen, slitmask, msarc, inmask):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Returns
        -------
        (self.arccen, self.arc_maskslit_
           self.arccen: ndarray, (nspec, nslit)
              arc spectrum for all slits
            self.arc_maskslit: ndarray, bool (nsit)
              boolean array containing a mask indicating which slits are good

        """
        # Full template kludge
        if self.par['method'] == 'full_template':
            nonlinear = 1e10
        else:
            nonlinear = self.nonlinear_counts
        # Do it
        # TODO -- Consider *not* passing in nonlinear_counts;  Probably should not mask saturated lines at this stage
        arccen, arc_maskslit = arc.get_censpec(slitcen, slitmask, msarc, inmask=inmask, nonlinear_counts=nonlinear)
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arc_maskslit

    def load_master(self, filename):
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
          self.par

        """


        # Does the master file exist?
        if not os.path.isfile(filename):
            msgs.warn("No Master frame found of type {:s}: {:s}".format(self.frametype, filename))
            return None
        else:
            msgs.info("Loading Master {0:s} frame:".format(self.frametype) + msgs.newline() + filename)
            self.wv_calib = linetools.utils.loadjson(filename)
            # Recast a few items as arrays
            for key in self.wv_calib.keys():
                if key in ['steps', 'par']:  # This isn't really necessary
                    continue
                for tkey in self.wv_calib[key].keys():
                    if isinstance(self.wv_calib[key][tkey], list):
                        self.wv_calib[key][tkey] = np.array(self.wv_calib[key][tkey])
            # parset
            if 'par' in self.wv_calib.keys():
                self.par = self.wv_calib['par'].copy()
            return self.wv_calib

    def save_master(self, data, outfile=None, raw_files=None, steps=None, overwrite=True, extensions=None, names=None):

        _outfile = self.ms_name if outfile is None else outfile
        if os.path.exists(_outfile) and (not overwrite):
            msgs.warn("This file already exists.  Use overwrite=True to overwrite it")
            return
        #
        msgs.info("Saving master {0:s} frame as:".format(self.frametype) + msgs.newline() + _outfile)
        # Wavelength fit(s)

        # jsonify has the annoying property that it modifies the objects when it jsonifies them so make a copy,
        # which converts lists to arrays, so we make a copy
        data_for_json = copy.deepcopy(data)
        gddict = linetools.utils.jsonify(data_for_json)
        linetools.utils.savejson(_outfile, gddict, easy_to_read=True, overwrite=True)
        # Finish
        msgs.info("Master {0:s} frame saved successfully:".format(self.frametype) + msgs.newline() + _outfile)

        return

    def make_maskslits(self, nslit):
        """
        (re)Generate the mask for wv_calib based on its contents
        This is the safest way to go...

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
            if key in ['steps', 'par', 'fit2d']:
                continue
            if (self.wv_calib[key] is not None) and (len(self.wv_calib[key]) > 0):
                mask[int(key)] = False
        self.maskslits = mask
        return self.maskslits

    def run(self, skip_QA=False, debug=False):
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
        skip_QA : bool, optional

        Returns
        -------
        self.wv_calib
        self.maskslits

        """
        ###############
        # Extract an arc down each slit
        self.arccen, self.maskslits = self.extract_arcs(self.slitcen, self.slitmask, self.msarc, self.inmask)

        # Fill up the calibrations and generate QA
        self.wv_calib = self.build_wv_calib(self.arccen, self.par['method'], skip_QA=skip_QA)

        # Return
        if self.par['echelle'] is True:
            fit2d_dict = self.echelle_2dfit(self.wv_calib, skip_QA = skip_QA, debug=debug)
            self.wv_calib['fit2d'] = fit2d_dict

        # Build mask
        self.make_maskslits(self.nslits)

        # Pack up
        self.wv_calib['steps'] = self.steps
        sv_par = self.par.data.copy()
        self.wv_calib['par'] = sv_par

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
            wavecal.qa.arc_fit_qa(self.wv_calib[str(slit)])

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


def load_wv_calib(filename):
    """
    Utility function which enables one to load the wv_calib and parset from a master file one line
    of code without instantiating the class.

    Note:  This method instantiates without a Spectrograph

    Parameters
    ----------
    filename: str
       Master file name
           slit : int, optional

    Returns
    -------
    Fills:
      wv_calib
        wv_calib dict
      par
        parset
    """


    waveCalib = WaveCalib(None,None)
    wv_calib = waveCalib.load_master(filename)
    return (waveCalib.wv_calib, waveCalib.par)


