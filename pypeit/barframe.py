"""
Module for generating a Bar image to map constant spatial locations

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import copy
import inspect
import numpy as np
from IPython import embed

from pypeit import msgs
from pypeit import masterframe
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg, extract


class BarFrame(calibrationimage.CalibrationImage, masterframe.MasterFrame):
    """
    Class to generate/load the bar image

    This class is primarily designed to generate a Bar frame to map constant
    spatial locations.  It also contains I/O methods for the Master frames
    of PypeIt.  The build_master() method will return a simple command
    (str) if that is the specified parameter (`par['useframe']`).

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.

        files (:obj:`list`, optional):
            List of filenames to process.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):
            The parameters used to process the frames.  If None, set
            to::

                pypeitpar.FrameGroupPar('bar')

        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        reuse_masters (:obj:`bool`, optional):
            Load from disk if possible
    """

    # Frame type is a class attribute
    frametype = 'bar'
    master_type = 'Bar'

    @classmethod
    def from_master_file(cls, master_file, par=None):
        """
        Instantiate from a master file

        Args:
            master_file (str):
            par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):

        Returns:
            barframe.BarFrame:
                The PypeItImage is loaded into self.pypeitImage

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(spectrograph, par=par, master_dir=master_dir, master_key=master_key,
                  reuse_masters=True)
        slf.pypeitImage = slf.load(ifile=master_file)
        # Return
        return slf

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False, msbias=None):

        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.CalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)
        # Process steps
        self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['orient']
        self.process_steps += ['apply_gain']

    def save(self, outfile=None, overwrite=True):
        """
        Save the master bar data.

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
        #
        hdr = self.build_master_header(steps=self.process_steps, raw_files=self.file_list)
        self.pypeitImage.write(_outfile, hdr=hdr, iext='BAR')
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the bar frame data from a saved master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            Returns a `numpy.ndarray`_ with the bar master frame image.
            Also returns the primary header, if requested.
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        # Load it up
        self.pypeitImage = pypeitimage.PypeItImage.from_file(master_file)
        return self.pypeitImage


class BarProfile(masterframe.MasterFrame):
    """
    Class to guide the determination of the bar traces

    Args:
        msbar (:class:`pypeit.images.pypeitimage.PypeItImage` or None):
            Bar image, created by the BarFrame class
        tslits_dict (dict or None):  TraceSlits dict
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph` or None):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        par (:class:`pypeit.par.pypeitpar.BarPar` or None):
            The parameters used for the bar traces
        det (int, optional): Detector number
        master_key (str, optional)
        master_dir (str, optional): Path to master frames
        reuse_masters (bool, optional):  Load from disk if possible
        qa_path (str, optional):  For QA
        msbpm (ndarray, optional): Bad pixel mask image

    Attributes:
        frametype : str
            Hard-coded to 'bar_prof'
        steps : list
            List of the processing steps performed
    """
    # Frametype is a class attribute
    frametype = 'bar_prof'
    master_type = 'BarProfile'

    def __init__(self, msbar, tslits_dict, spectrograph, par, det=1,
                 master_key=None, master_dir=None, reuse_masters=False,
                 qa_path=None, msbpm=None):

        # MasterFrame
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, file_format='json',
                                         reuse_masters=reuse_masters)

        # Required parameters (but can be None)
        self.msbar = msbar
        self.tslits_dict = tslits_dict
        self.spectrograph = spectrograph
        self.par = par

        # Optional parameters
        self.bpm = msbar.mask if msbpm is None else msbpm
        self.qa_path = qa_path
        self.det = det
        self.master_key = master_key

        # Attributes
        self.steps = []  # steps executed

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
        if self.tslits_dict is not None and self.msbar is not None:
            self.slitmask_science = pixels.tslits2mask(self.tslits_dict)
            gpm = self.bpm == 0 if self.bpm is not None \
                else np.ones_like(self.slitmask_science, dtype=bool)
            self.shape_science = self.slitmask_science.shape
            self.shape_arc = self.msbar.image.shape
            self.nslits = self.tslits_dict['slit_left'].shape[1]
            self.slit_left = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                  self.tslits_dict['slit_left'])
            self.slit_righ = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                  self.tslits_dict['slit_righ'])
            self.slitcen = arc.resize_slits2arc(self.shape_arc, self.shape_science,
                                                self.tslits_dict['slitcen'])
            self.slitmask = arc.resize_mask2arc(self.shape_arc, self.slitmask_science)
            self.gpm = arc.resize_mask2arc(self.shape_arc, gpm)
            # We want even the saturated lines in full_template for the cross-correlation
            #   They will be excised in the detect_lines() method on the extracted arc
            if self.par['method'] != 'full_template':
                self.gpm &= self.msbar.image < self.nonlinear_counts
            self.slit_spat_pos = edgetrace.slit_spat_pos(self.tslits_dict['slit_left'],
                                                         self.tslits_dict['slit_righ'],
                                                         self.tslits_dict['nspec'],
                                                         self.tslits_dict['nspat'])

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

    def build_traces(self, image, show_peaks=False, show_trace=False, debug=False):
        """
        Main routine to generate the bar profile traces in all slits

        Args:
             image (np.ndarray):
             show_peaks (bool, optional):
               Generate QA showing peaks identified by object finding
             show_trace (bool, optional):
               Generate QA  showing traces identified. Requires an open ginga RC modules window
             debug (bool, optional):

        Returns:
            dict:  self.bar_prof
        """
        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        inmask = self.sciImg.mask == 0
        # Find objects
        specobj_dict = {'setup': self.setup, 'slitid': 999,  # 'orderindx': 999,
                        'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}
        # TODO This is a bad idea -- we want to find everything for standards
        # sig_thresh = 30.0 if std else self.redux_par['sig_thresh']
        sobjs_ech, _ = extract.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, self.tslits_dict['slit_left'],
            self.tslits_dict['slit_righ'], self.order_vec, self.maskslits,
            spec_min_max=np.vstack((self.tslits_dict['spec_min'],
                                    self.tslits_dict['spec_max'])),
            inmask=inmask, ir_redux=self.ir_redux, ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            hand_extract_dict=manual_extract_dict, plate_scale=plate_scale,
            std_trace=std_trace,
            specobj_dict=specobj_dict, sig_thresh=self.par['reduce']['findobj']['sig_thresh'],
            show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
            cont_fit=self.par['reduce']['findobj']['find_cont_fit'],
            npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
            fwhm=self.par['reduce']['findobj']['find_fwhm'],
            maxdev=self.par['reduce']['findobj']['find_maxdev'],
            max_snr=self.par['reduce']['findobj']['ech_find_max_snr'],
            min_snr=self.par['reduce']['findobj']['ech_find_min_snr'],
            nabove_min_snr=self.par['reduce']['findobj']['ech_find_nabove_min_snr'],
            show_trace=show_trace, debug=debug)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image * (self.sciImg.mask == 0), chname='ech_objfind', sobjs=sobjs_ech,
                      slits=False)

        return sobjs_ech, len(sobjs_ech), skymask
        # QA
        if not skip_QA:
            for slit in ok_mask:
                outfile = qa.set_qa_filename(self.master_key, 'arc_fit_qa', slit=slit,
                                             out_dir=self.qa_path)
                autoid.arc_fit_qa(self.wv_calib[str(slit)], outfile=outfile)

        # Return
        self.steps.append(inspect.stack()[0][3])
        return self.bar_prof

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
        all_order = np.array([], dtype=float)

        # Obtain a list of good slits
        ok_mask = np.where(np.invert(self.maskslits))[0]
        nspec = self.msbar.image.shape[0]
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
            self.slitcen, self.slitmask, self.msbar.image, gpm=self.gpm)
        # , nonlinear_counts=nonlinear) -- Non-linear counts are already part of the gpm
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

    def run(self, skip_QA=False, debug=False):
        """
        Main driver for bar profile tracing

        Args:
            skip_QA : bool, optional

        Returns:
            dict, ndarray:  wv_calib dict and maskslits bool array

        """
        ###############
        # Fill up the calibrations and generate QA
        self.bar_prof = self.build_traces(skip_QA=skip_QA)

        # Pack up
        self.bar_prof['steps'] = self.steps
        sv_par = self.par.data.copy()
        self.bar_prof['par'] = sv_par

        return self.bar_prof

    def show(self, item, slit=None):
        """
        Show one of the class internals

        Args:
            item (str):
            slit (int, optional):

        Returns:

        """
        return

    def __repr__(self):
        # Generate sets string
        txt = '<{:s}: '.format(self.__class__.__name__)
        if len(self.steps) > 0:
            txt += ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2] + ']'  # Trim the trailing comma
        txt += '>'
        return txt
