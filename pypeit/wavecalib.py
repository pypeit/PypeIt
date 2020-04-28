"""
Module for guiding 1D Wavelength Calibration

.. include:: ../links.rst
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
from pypeit.core import arc, qa
from pypeit.core.wavecal import autoid, waveio, templates
from pypeit.core.gui import identify as gui_identify
from pypeit import utils
from pypeit import datamodel

#class WaveCalib(datamodel.DataContainer):
#    # Peg the version of this class to that of PypeItImage
#    version = '1.0.0'
#
#    # I/O
#    output_to_disk = None
#    hdu_prefix = None
#
#    # Master fun
#    frametype = 'wv_calib'
#    master_type = 'WaveCalib'
#
#    # Data model
#    datamodel_v100 = {
#        'image': dict(otype=np.ndarray, atype=np.floating, desc='Main data image'),
#        'ivar': dict(otype=np.ndarray, atype=np.floating, desc='Main data inverse variance image'),
#    }
#
#    datamodel = datamodel_v100.copy()

class WaveCalib(object):
    """
    Class to guide wavelength calibration

    Args:
        msarc (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Arc image, created by the ArcImage class
        slits (:class:`pypeit.slittrace.SlitTraceSet`, None):
            Slit edges
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.WaveSolutionPar` or None):
            The parameters used for the wavelength solution
            Uses ['calibrations']['wavelengths']
        binspectral (int, optional): Binning of the Arc in the spectral dimension
        det (int, optional): Detector number
        msbpm (ndarray, optional): Bad pixel mask image
        qa_path (str, optional):  For QA
        master_key (:obj:`str`, optional):  For naming QA only

    Attributes:
        frametype : str
            Hard-coded to 'wv_calib'
        steps : list
            List of the processing steps performed
        wv_calib : dict
            Primary output.  Keys 0, 1, 2, 3 are solution for individual
            previously slit steps
        arccen (ndarray):
            (nwave, nslit) Extracted arc(s) down the center of the slit(s)
        maskslits : ndarray (nslit); bool
            Slits to ignore because they were not extracted. WARNING:
            Outside of this Class, it is best to regenerate the mask
            using  make_maskslits()
        gpm (`numpy.ndarray`_):
            Good pixel mask
            Eventually, we might attach this to self.msarc although that would then
            require that we write it to disk with self.msarc.image
        nonlinear_counts (float):
            Specifies saturation level for the arc lines
        wvc_bpm (`numpy.ndarray`_):  Mask for slits attempted to have a wv_calib solution
    """
    # Frametype is a class attribute
    frametype = 'wv_calib'
    master_type = 'WaveCalib'
    master_file_format = 'json'

    def __init__(self, msarc, slits, spectrograph, par, binspectral=None, det=1,
                 qa_path=None, msbpm=None, master_key=None):

        # Required parameters (but can be None)
        self.msarc = msarc
        self.slits = slits
        self.spectrograph = spectrograph
        self.par = par

        # Optional parameters
        self.bpm = msbpm
        if self.bpm is None:
            if msarc is not None:  # Can be None for load;  will remove this for DataContainer
                self.bpm = msarc.bpm
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
            else self.spectrograph.nonlinear_counts(self.msarc.detector)
            #else self.spectrograph.nonlinear_counts(self.det)

        # --------------------------------------------------------------
        # TODO: Build another base class that does these things for both
        # WaveTilts and WaveCalib?

        # Set the slitmask and slit boundary related attributes that the
        # code needs for execution. This also deals with arcimages that
        # have a different binning then the trace images used to defined
        # the slits
        if self.slits is not None and self.msarc is not None:
            # Load up slits
            # TODO -- Allow for flexure
            all_left, all_right, mask = self.slits.select_edges(initial=True, flexure=None)  # Grabs all, init slits + flexure
            self.spat_coo = self.slits.spatial_coordinates()  # All slits, even masked
            # Internal mask for failed wv_calib analysis
            # TODO -- Allow for an option to re-attempt those previously flagged as BADWVCALIB?
            self.wvc_bpm = np.invert(mask == 0)
            self.wvc_bpm_init = self.wvc_bpm.copy()
            # Slitmask -- Grabs only unmasked, initial slits
            self.slitmask_science = self.slits.slit_img(initial=True, flexure=None)
            # Resize
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msarc.image.shape
            # slitcen is padded to include slits that may be masked, for convenience in coding downstream
            self.slitcen = arc.resize_slits2arc(self.shape_arc, self.shape_science, (all_left+all_right)/2)
            self.slitmask = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            # Mask
            gpm = self.bpm == 0 if self.bpm is not None \
                else np.ones_like(self.slitmask_science, dtype=bool)
            self.gpm = arc.resize_mask2arc(self.shape_arc, gpm)
            # We want even the saturated lines in full_template for the cross-correlation
            #   They will be excised in the detect_lines() method on the extracted arc
            if self.par['method'] != 'full_template':
                self.gpm &= self.msarc.image < self.nonlinear_counts

        else:
            self.slitmask_science = None
            self.shape_science = None
            self.shape_arc = None
            self.nslits = 0
            self.slitcen = None
            self.slitmask = None
            self.gpm = None

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
        ok_mask_idx = np.where(np.invert(self.wvc_bpm))[0]

        # Obtain calibration for all slits
        if method == 'simple':
            lines = self.par['lamps']
            line_lists = waveio.load_line_lists(lines)

            final_fit = arc.simple_calib_driver(line_lists, arccen, ok_mask_idx,
                                                    n_final=self.par['n_final'],
                                                    sigdetect=self.par['sigdetect'],
                                                    IDpixels=self.par['IDpixels'],
                                                    IDwaves=self.par['IDwaves'])
#        elif method == 'basic':
#            final_fit = {}
#            for slit in ok_mask:
#                status, ngd_match, match_idx, scores, ifinal_fit = \
#                        autoid.basic(arccen[:, slit], self.par['lamps'], self.par['wv_cen'],
#                                     self.par['disp'], nonlinear_counts=self.nonlinear_counts)
#                final_fit[str(slit)] = ifinal_fit.copy()
#                if status != 1:
#                    self.maskslits[slit] = True
        elif method == 'holy-grail':
            # Sometimes works, sometimes fails
            arcfitter = autoid.HolyGrail(arccen, par=self.par, ok_mask=ok_mask_idx, nonlinear_counts=self.nonlinear_counts)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'identify':
            final_fit = {}
            # Manually identify lines
            msgs.info("Initializing the wavelength calibration tool")
            # TODO: Move this loop to the GUI initalise method
            embed()
            for slit_idx in ok_mask_idx:
                arcfitter = gui_identify.initialise(arccen, slit=slit_idx, par=self.par)
                final_fit[str(slit_idx)] = arcfitter.get_results()
                if final_fit[str(slit_idx)] is not None:
                    ans = 'y'
                    # ans = ''
                    # while ans != 'y' and ans != 'n':
                    #     ans = input("Would you like to store this wavelength solution in the archive? (y/n): ")
                    if ans == 'y' and final_fit[str(slit_idx)]['rms'] < self.par['rms_threshold']:
                        # Store the results in the user reid arxiv
                        specname = self.spectrograph.spectrograph
                        gratname = "UNKNOWN"  # input("Please input the grating name: ")
                        dispangl = "UNKNOWN"  # input("Please input the dispersion angle: ")
                        templates.pypeit_identify_record(final_fit[str(slit_idx)], self.binspectral, specname, gratname, dispangl)
                        msgs.info("Your wavelength solution has been stored")
                        msgs.info("Please consider sending your solution to the PYPEIT team!")

        elif method == 'reidentify':
            # Now preferred
            # Slit positions
            arcfitter = autoid.ArchiveReid(arccen, self.spectrograph, self.par, ok_mask=ok_mask_idx,
                                           slit_spat_pos=self.spat_coo,
                                           nonlinear_counts=self.nonlinear_counts)
            patt_dict, final_fit = arcfitter.get_results()
        elif method == 'full_template':
            # Now preferred
            if self.binspectral is None:
                msgs.error("You must specify binspectral for the full_template method!")
            final_fit = autoid.full_template(arccen, self.par, ok_mask_idx, self.det,
                                             self.binspectral, nonlinear_counts=self.nonlinear_counts,
                                             nsnippet=self.par['nsnippet'])
        else:
            msgs.error('Unrecognized wavelength calibration method: {:}'.format(method))

        # Convert keys to spatial system
        self.wv_calib = {}
        tmp = copy.deepcopy(final_fit)
        for idx in range(self.slits.nslits):
            if str(idx) in final_fit.keys():
                self.wv_calib[str(self.slits.slitord_id[idx])] = final_fit.pop(str(idx))

        # Update mask
        self.update_wvmask()

        #TODO For generalized echelle (not hard wired) assign order number here before, i.e. slits.ech_order

        # QA
        if not skip_QA:
            ok_mask_idx = np.where(np.invert(self.wvc_bpm))[0]
            for slit_idx in ok_mask_idx:
                outfile = qa.set_qa_filename(self.master_key, 'arc_fit_qa', slit=self.slits.slitord_id[slit_idx],
                                             out_dir=self.qa_path)
                autoid.arc_fit_qa(self.wv_calib[str(self.slits.slitord_id[slit_idx])], outfile=outfile)


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
        ok_mask_idx = np.where(np.invert(self.wvc_bpm))[0]
        ok_mask_order = self.slits.slitord_id[ok_mask_idx]
        nspec = self.msarc.image.shape[0]
        for iorder in wv_calib.keys():  # Spatial based
            if int(iorder) not in ok_mask_order:
                continue
            #try:
            #    iorder, iindx = self.spectrograph.slit2order(self.spat_coo[self.slits.spatid_to_zero(int(islit))])
            #except:
            #    embed()
            mask_now = wv_calib[iorder]['mask']
            all_wave = np.append(all_wave, wv_calib[iorder]['wave_fit'][mask_now])
            all_pixel = np.append(all_pixel, wv_calib[iorder]['pixel_fit'][mask_now])
            all_order = np.append(all_order, np.full_like(wv_calib[iorder]['pixel_fit'][mask_now],
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
                  True = masked (bad)

        """
        # Do it on the slits not masked in self.slitmask
        arccen, arccen_bpm, arc_maskslit = arc.get_censpec(
            self.slitcen, self.slitmask, self.msarc.image, gpm=self.gpm, slit_bpm=self.wvc_bpm)
        # Step
        self.steps.append(inspect.stack()[0][3])

        # Update the mask
        self.wvc_bpm |= arc_maskslit

        return arccen, self.wvc_bpm

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
        _outfile = outfile # self.master_file_path if outfile is None else outfile
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

    def load(self, ifile):
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
        self.wv_calib = waveio.load_wavelength_calibration(ifile)
        return self.wv_calib

    def update_wvmask(self):
        """
        (re)Generate the mask for wv_calib based on its contents
        This is the safest way to go...

        Args:
            nslit (int): Number of slits/orders

        Returns:
            `numpy.ndarray`_: self.wvc_bpm, boolean array -- True = masked, i.e. do not use

        """
        # Update mask based on wv_calib
        for key in self.wv_calib.keys():
            if key in ['steps', 'par', 'fit2d', 'bpm']:
                continue
            if (self.wv_calib[key] is None) or (len(self.wv_calib[key]) == 0):
                try:
                    idx = self.slits.spatid_to_zero(int(key))
                except:
                    embed(header='428 of wavecalib')
                self.wvc_bpm[idx] = True

    def run(self, skip_QA=False, debug=False):
        """
        Main driver for wavelength calibration

        Code flow:
          1. Extract 1D arc spectra down the center of each unmasked slit/order
          2. Load the parameters guiding wavelength calibration
          3. Generate the 1D wavelength fits
          4. Generate a mask

        Args:
            skip_QA : bool, optional

        Returns:
            dict:  wv_calib dict

        """
        ###############
        # Extract an arc down each slit
        self.arccen, self.wvc_bpm = self.extract_arcs()

        # Fill up the calibrations and generate QA
        self.wv_calib = self.build_wv_calib(self.arccen, self.par['method'], skip_QA=skip_QA)

        # Fit 2D?
        if self.par['echelle'] is True:
            fit2d_dict = self.echelle_2dfit(self.wv_calib, skip_QA = skip_QA, debug=debug)
            self.wv_calib['fit2d'] = fit2d_dict

        # Deal with mask
        self.update_wvmask()

        # Any masked during this analysis?
        wv_masked = np.where(np.invert(self.wvc_bpm_init) & self.wvc_bpm)[0]
        if len(wv_masked) > 0:
            self.slits.mask[wv_masked] = self.slits.bitmask.turn_on(
                    self.slits.mask[wv_masked], 'BADWVCALIB')

        # Pack up
        self.wv_calib['steps'] = self.steps
        sv_par = self.par.data.copy()
        self.wv_calib['par'] = sv_par

        return self.wv_calib

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



# TODO -- Move this as a method on a WaveCalib DataContainer
def build_waveimg(spectrograph, tilts, slits, wv_calib, spat_flexure=None):
    """
    Main algorithm to build the wavelength image

    Only applied to good slits, which means any non-flagged or flagged
     in the exclude_for_reducing list

    Args:
        spectrograph (:obj:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph object
        tilts (`numpy.ndarray`_):
            Image holding tilts
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
        wv_calib (dict):
        spat_flexure (float, optional):

    Returns:
        `numpy.ndarray`_: The wavelength image.
    """
    # Setup
    #ok_slits = slits.mask == 0
    bpm = slits.mask.astype(bool)
    bpm &= np.invert(slits.bitmask.flagged(slits.mask, flag=slits.bitmask.exclude_for_reducing))
    ok_slits = np.invert(bpm)
    #
    image = np.zeros_like(tilts)
    slitmask = slits.slit_img(flexure=spat_flexure, exclude_flag=slits.bitmask.exclude_for_reducing)

    par = wv_calib['par']
    slit_spat_pos = slits.spatial_coordinates(flexure=spat_flexure)

    # If this is echelle print out a status message and do some error checking
    if par['echelle']:
        msgs.info('Evaluating 2-d wavelength solution for echelle....')
        if len(wv_calib['fit2d']['orders']) != np.sum(ok_slits):
            msgs.error('wv_calib and ok_slits do not line up. Something is very wrong!')

    # Unpack some 2-d fit parameters if this is echelle
    for slit_spat in slits.spat_id[ok_slits]:
        thismask = (slitmask == slit_spat)
        if not np.any(thismask):
            msgs.error("Something failed in wavelengths or masking..")
        if par['echelle']:
            # TODO: Put this in `SlitTraceSet`?
            order, indx = spectrograph.slit2order(slit_spat_pos[slits.spatid_to_zero(slit_spat)])
            # evaluate solution
            image[thismask] = utils.func_val(wv_calib['fit2d']['coeffs'],
                                             tilts[thismask],
                                             wv_calib['fit2d']['func2d'],
                                             x2=np.ones_like(tilts[thismask])*order,
                                             minx=wv_calib['fit2d']['min_spec'],
                                             maxx=wv_calib['fit2d']['max_spec'],
                                             minx2=wv_calib['fit2d']['min_order'],
                                             maxx2=wv_calib['fit2d']['max_order'])
            image[thismask] /= order
        else:
            #iwv_calib = wv_calib[str(slit)]
            iwv_calib = wv_calib[str(slit_spat)]
            image[thismask] = utils.func_val(iwv_calib['fitc'], tilts[thismask],
                                             iwv_calib['function'],
                                             minx=iwv_calib['fmin'],
                                             maxx=iwv_calib['fmax'])
    # Return
    return image


