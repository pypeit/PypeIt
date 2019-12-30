"""
Module for guiding 1D Wavelength Calibration
"""
import os
import copy
import inspect

from IPython import embed

import numpy as np

from matplotlib import pyplot as plt

import linetools.utils

from pypeit import msgs
from pypeit import masterframe
from pypeit import edgetrace
from pypeit.core import arc, qa, pixels
from pypeit.core.wavecal import autoid, waveio, templates
from pypeit.core.gui import identify as gui_identify



class WaveCalib(masterframe.MasterFrame):
    """
    Class to guide wavelength calibration

    Args:
        msarc (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Arc image, created by the ArcImage class
        tslits_dict (dict or None):  TraceSlits dict
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.WaveSolutionPar` or None):
            The parameters used for the wavelength solution
        binspectral (int, optional): Binning of the Arc in the spectral dimension
        det (int, optional): Detector number
        master_key (str, optional)
        master_dir (str, optional): Path to master frames
        reuse_masters (bool, optional):  Load from disk if possible
        qa_path (str, optional):  For QA
        msbpm (ndarray, optional): Bad pixel mask image

    Attributes:
        frametype : str
            Hard-coded to 'wv_calib'
        steps : list
            List of the processing steps performed
        wv_calib : dict
            Primary output.  Keys 0, 1, 2, 3 are solution for individual
            slit steps
        arccen (ndarray):
            (nwave, nslit) Extracted arc(s) down the center of the slit(s)
        maskslits : ndarray (nslit); bool
            Slits to ignore because they were not extracted. WARNING:
            Outside of this Class, it is best to regenerate the mask
            using  make_maskslits()
        gpm (np.ndarray):
            Good pixel mask
            Eventually, we might attach this to self.msarc although that would then
            require that we write it to disk with self.msarc.image
    """
    # Frametype is a class attribute
    frametype = 'wv_calib'
    master_type = 'WaveCalib'

    def __init__(self, msarc, tslits_dict, spectrograph, par, binspectral=None, det=1,
                 master_key=None, master_dir=None, reuse_masters=False, qa_path=None,
                 msbpm=None):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, file_format='json',
                                         reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.msarc = msarc
        self.tslits_dict = tslits_dict
        self.spectrograph = spectrograph
        self.par = par

        # Optional parameters
        self.bpm = msbpm
        self.binspectral = binspectral
        self.qa_path = qa_path
        self.det = det
        self.master_key = master_key

        # Attributes
        self.steps = []     # steps executed
        self.wv_calib = {}  # main output
        self.arccen = None  # central arc spectrum

        # Get the non-linear count level
        self.nonlinear_counts = 1e10 if self.spectrograph is None \
            else self.spectrograph.nonlinear_counts(self.det)

        # --------------------------------------------------------------
        # TODO: Build another base class that does these things for both
        # WaveTilts and WaveCalib?

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits
        if self.tslits_dict is not None and self.msarc is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            gpm = self.bpm == 0 if self.bpm is not None \
                                    else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.image.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                  self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                  self.tslits_dict['slit_righ'])
            self.slitcen = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                  self.tslits_dict['slitcen'])
            self.slitmask  = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.gpm = arc.resize_mask2arc(self.shape_arc, gpm)
            # We want even the saturated lines in full_template for the cross-correlation
            #   They will be excised in the detect_lines() method on the extracted arc
            if self.par['method'] != 'full_template':
                self.gpm &= self.msarc.image < self.nonlinear_counts
            self.slit_spat_pos = edgetrace.slit_spat_pos(self.tslits_dict)

        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.nslits = 0
            self.slit_left = None
            self.slit_righ = None
            self.slitcen = None
            self.slitmask = None
            self.gpm = None
            self.nonlinear_counts = None

    def build_wv_calib(self, arccen, method, skip_QA=False):
        """
        Main routine to generate the wavelength solutions in a loop over slits
        Wrapper to arc.simple_calib or arc.calib_with_arclines

        self.maskslits is updated for slits that fail

        Args:
            method : str
              'simple' -- arc.simple_calib
              'arclines' -- arc.calib_with_arclines
              'holy-grail' -- wavecal.autoid.HolyGrail
              'reidentify' -- wavecal.auotid.ArchiveReid
              'identify' -- wavecal.identify.Identify
              'full_template' -- wavecal.auotid.full_template
            skip_QA (bool, optional)

        Returns:
            dict:  self.wv_calib
        """
        # Obtain a list of good slits
        ok_mask = np.where(~self.maskslits)[0]

        # Obtain calibration for all slits
        if method == 'simple':
            lines = self.par['lamps']
            line_lists = waveio.load_line_lists(lines)

            final_fit = arc.simple_calib_driver(line_lists, arccen, ok_mask,
                                                    n_final=self.par['n_final'],
                                                    sigdetect=self.par['sigdetect'],
                                                    IDpixels=self.par['IDpixels'],
                                                    IDwaves=self.par['IDwaves'])
        elif method == 'semi-brute':
            # TODO: THIS IS CURRENTLY BROKEN
            debugger.set_trace()
            final_fit = {}
            for slit in ok_mask:
                # HACKS BY JXP
                self.par['wv_cen'] = 8670.
                self.par['disp'] = 1.524
                # ToDO remove these hacks and use the parset in semi_brute
                best_dict, ifinal_fit \
                        = autoid.semi_brute(arccen[:, slit], self.par['lamps'], self.par['wv_cen'],
                                            self.par['disp'], match_toler=self.par['match_toler'],
                                            func=self.par['func'], n_first=self.par['n_first'],
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
                        autoid.basic(arccen[:, slit], self.par['lamps'], self.par['wv_cen'],
                                     self.par['disp'], nonlinear_counts=self.nonlinear_counts)
                final_fit[str(slit)] = ifinal_fit.copy()
                if status != 1:
                    self.maskslits[slit] = True
        elif method == 'holy-grail':
            # Sometimes works, sometimes fails
            arcfitter = autoid.HolyGrail(arccen, par=self.par, ok_mask=ok_mask)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'identify':
            # Manually identify lines
            msgs.info("Initializing the wavelength calibration tool")
            # Todo : Generalise to multislit case
            arcfitter = gui_identify.initialise(arccen, par=self.par)
            final_fit = arcfitter.get_results()
            slit = 0
            if final_fit[str(slit)] is not None:
                ans = 'y'
                # ans = ''
                # while ans != 'y' and ans != 'n':
                #     ans = input("Would you like to store this wavelength solution in the archive? (y/n): ")
                if ans == 'y' and final_fit[str(slit)]['rms'] < 0.1:
                    # Store the results in the user reid arxiv
                    specname = self.spectrograph.spectrograph
                    gratname = "UNKNOWN"  # input("Please input the grating name: ")
                    dispangl = "UNKNOWN"  # input("Please input the dispersion angle: ")
                    templates.pypeit_identify_record(final_fit[str(slit)], self.binspectral, specname, gratname, dispangl)
                    msgs.info("Your wavelength solution has been stored")
                    msgs.info("Please consider sending your solution to the PYPEIT team!")
        elif method == 'reidentify':
            # Now preferred
            # Slit positions
            arcfitter = autoid.ArchiveReid(arccen, self.spectrograph, self.par, ok_mask=ok_mask,
                                           slit_spat_pos=self.slit_spat_pos)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'full_template':
            # Now preferred
            if self.binspectral is None:
                msgs.error("You must specify binspectral for the full_template method!")
            final_fit = autoid.full_template(arccen, self.par, ok_mask, self.det,
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
                outfile = qa.set_qa_filename(self.master_key, 'arc_fit_qa', slit=slit,
                                             out_dir=self.qa_path)
                autoid.arc_fit_qa(self.wv_calib[str(slit)], outfile=outfile)

        # Return
        self.steps.append(inspect.stack()[0][3])
        return self.wv_calib

    def echelle_2dfit(self, wv_calib, debug=False, skip_QA=False):
        """
        Evaluate 2-d wavelength solution for echelle data. Unpacks
        wv_calib for slits to be input into  arc.fit2darc

        Args:
            wv_calib (dict): Wavelength calibration
            debug (bool, optional):  Show debugging info
            skip_QA (bool, optional): Skip QA

        Returns:
            dict: dictionary containing information from 2-d fit

        """
        msgs.info('Fitting 2-d wavelength solution for echelle....')
        all_wave = np.array([], dtype=float)
        all_pixel = np.array([], dtype=float)
        all_order = np.array([],dtype=float)

        # Obtain a list of good slits
        ok_mask = np.where(~self.maskslits)[0]
        nspec = self.msarc.image.shape[0]
        for islit in wv_calib.keys():
            if int(islit) not in ok_mask:
                continue
            iorder, iindx = self.spectrograph.slit2order(self.slit_spat_pos[int(islit)])
            mask_now = wv_calib[islit]['mask']
            all_wave = np.append(all_wave, wv_calib[islit]['wave_fit'][mask_now])
            all_pixel = np.append(all_pixel, wv_calib[islit]['pixel_fit'][mask_now])
            all_order = np.append(all_order, np.full_like(wv_calib[islit]['pixel_fit'][mask_now],
                                                          float(iorder)))

        # Fit
        fit2d_dict = arc.fit2darc(all_wave, all_pixel, all_order, nspec,
                                  nspec_coeff=self.par['ech_nspec_coeff'],
                                  norder_coeff=self.par['ech_norder_coeff'],
                                  sigrej=self.par['ech_sigrej'], debug=debug)

        self.steps.append(inspect.stack()[0][3])

        # QA
        if not skip_QA:
            outfile_global = qa.set_qa_filename(self.master_key, 'arc_fit2d_global_qa',
                                                out_dir=self.qa_path)
            arc.fit2darc_global_qa(fit2d_dict, outfile=outfile_global)
            outfile_orders = qa.set_qa_filename(self.master_key, 'arc_fit2d_orders_qa',
                                                out_dir=self.qa_path)
            arc.fit2darc_orders_qa(fit2d_dict, outfile=outfile_orders)

        return fit2d_dict

    # TODO: JFH this method is identical to the code in wavetilts.
    # SHould we make it a separate function?
    def extract_arcs(self):
        """
        Extract the arcs down each slit/order

        Wrapper to arc.get_censpec()

        Args:

        Returns:
            tuple: Returns the following:
                - self.arccen: ndarray, (nspec, nslit): arc spectrum for
                  all slits
                - self.arc_maskslit: ndarray, bool (nsit): boolean array
                  containing a mask indicating which slits are good

        """
        # Do it
        arccen, arccen_bpm, arc_maskslit = arc.get_censpec(
            self.slitcen, self.slitmask, self.msarc.image,
            gpm=self.gpm)  #, nonlinear_counts=nonlinear) -- Non-linear counts are already part of the gpm
        # Step
        self.steps.append(inspect.stack()[0][3])
        return arccen, arc_maskslit

    def save(self, outfile=None, overwrite=True):
        """
        Save the wavelength calibration data to a master frame.

        This is largely a wrapper for
        :func:`pypeit.core.wavecal.waveio.save_wavelength_calibration`.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return

        # Report and save

        # jsonify has the annoying property that it modifies the objects
        # when it jsonifies them so make a copy, which converts lists to
        # arrays, so we make a copy
        data_for_json = copy.deepcopy(self.wv_calib)
        gddict = linetools.utils.jsonify(data_for_json)
        linetools.utils.savejson(_outfile, gddict, easy_to_read=True, overwrite=True)
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load a full (all slit) wavelength calibration.

        This is largely a wrapper for
        :func:`pypeit.core.wavecal.waveio.load_wavelength_calibration`.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`master_file_path`.

        Returns:
            dict or None: self.wv_calib
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        # Read, save it to self, return
        msgs.info('Loading Master frame: {0}'.format(master_file))
        self.wv_calib = waveio.load_wavelength_calibration(master_file)
        return self.wv_calib

    def make_maskslits(self, nslit):
        """
        (re)Generate the mask for wv_calib based on its contents
        This is the safest way to go...

        Args:
            nslit (int): Number of slits/orders

        Returns:
            ndarray: self.maskslits, boolean array -- True = masked, i.e. do not use

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

        Args:
            skip_QA : bool, optional

        Returns:
            dict, ndarray:  wv_calib dict and maskslits bool array

        """
        ###############
        # Extract an arc down each slit
        self.arccen, self.maskslits = self.extract_arcs()

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

        Args:
            item (str):
              'spec' -- Show the fitted points and solution;  requires slit
              'fit' -- Show fit QA; requires slit
            slit (int, optional):

        Returns:

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
            autoid.arc_fit_qa(self.wv_calib[str(slit)])

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


