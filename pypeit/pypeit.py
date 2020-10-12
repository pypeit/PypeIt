"""
Main driver class for PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import time
import os
import numpy as np
import copy
from astropy.io import fits
from pypeit import msgs
from pypeit import calibrations
from pypeit.images import buildimage
from pypeit.display import display
from pypeit import reduce
from pypeit import spec2dobj
from pypeit.core import qa
from pypeit.core import extract
from pypeit import specobjs
from pypeit.spectrographs.util import load_spectrograph
from pypeit import slittrace

from configobj import ConfigObj
from pypeit.par.util import parse_pypeit_file
from pypeit.par import PypeItPar
from pypeit.metadata import PypeItMetaData

from IPython import embed

class PypeIt(object):
    """
    This class runs the primary calibration and extraction in PypeIt

    .. todo::
        Fill in list of attributes!

    Args:
        pypeit_file (:obj:`str`):
            PypeIt filename.
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be:

                - 0: No output
                - 1: Minimal output (default)
                - 2: All output

        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        reuse_masters (:obj:`bool`, optional):
            Reuse any pre-existing calibration files
        logname (:obj:`str`, optional):
            The name of an ascii log file with the details of the
            reduction.
        show: (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via ``ginga --modules=RC &``
        redux_path (:obj:`str`, optional):
            Over-ride reduction path in PypeIt file (e.g. Notebook usage)
        calib_only: (:obj:`bool`, optional):
            Only generate the calibration files that you can

    Attributes:
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  PypeIt files have a
            specific set of valid formats. A description can be found
            :ref:`pypeit_file`.
        fitstbl (:obj:`pypeit.metadata.PypeItMetaData`): holds the meta info

    """
#    __metaclass__ = ABCMeta

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_masters=False, logname=None,
                 show=False, redux_path=None, calib_only=False):

        # Load
        cfg_lines, data_files, frametype, usrdata, setups \
                = parse_pypeit_file(pypeit_file, runtime=True)
        self.pypeit_file = pypeit_file
        self.calib_only = calib_only

        # Spectrograph
        cfg = ConfigObj(cfg_lines)
        spectrograph_name = cfg['rdx']['spectrograph']
        self.spectrograph = load_spectrograph(spectrograph_name)
        msgs.info('Loaded spectrograph {0}'.format(self.spectrograph.spectrograph))

        # --------------------------------------------------------------
        # Get the full set of PypeIt parameters
        #   - Grab a science or standard file for configuration specific parameters

        config_specific_file = None
        for idx, row in enumerate(usrdata):
            if ('science' in row['frametype']) or ('standard' in row['frametype']):
                config_specific_file = data_files[idx]
        # search for arcs, trace if no scistd was there
        if config_specific_file is None:
            for idx, row in enumerate(usrdata):
                if ('arc' in row['frametype']) or ('trace' in row['frametype']):
                    config_specific_file = data_files[idx]
        if config_specific_file is not None:
            msgs.info(
                'Setting configuration-specific parameters using {0}'.format(os.path.split(config_specific_file)[1]))
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(config_specific_file).to_config()

        #   - Build the full set, merging with any user-provided
        #     parameters
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)
        msgs.info('Built full PypeIt parameter set.')

        # Check the output paths are ready
        if redux_path is not None:
            self.par['rdx']['redux_path'] = redux_path

        # TODO: Write the full parameter set here?
        # --------------------------------------------------------------

        # --------------------------------------------------------------
        # Build the meta data
        #   - Re-initilize based on the file data
        msgs.info('Compiling metadata')
        self.fitstbl = PypeItMetaData(self.spectrograph, self.par, files=data_files,
                                      usrdata=usrdata, strict=True)
        #   - Interpret automated or user-provided data from the PypeIt
        #   file
        self.fitstbl.finalize_usr_build(frametype, setups[0])
        # --------------------------------------------------------------
        #   - Write .calib file (For QA naming amongst other things)
        calib_file = pypeit_file.replace('.pypeit', '.calib')
        self.fitstbl.write_calib(calib_file)

        # Other Internals
        self.logname = logname
        self.overwrite = overwrite

        # Currently the runtime argument determines the behavior for
        # reuse_masters.
        self.reuse_masters = reuse_masters
        self.show = show

        # Set paths
        self.calibrations_path = os.path.join(self.par['rdx']['redux_path'], self.par['calibrations']['master_dir'])

        # Check for calibrations
        if not self.calib_only:
            calibrations.check_for_calibs(self.par, self.fitstbl,
                                          raise_error=self.par['calibrations']['raise_chk_error'])

        # Report paths
        msgs.info('Setting reduction path to {0}'.format(self.par['rdx']['redux_path']))
        msgs.info('Master calibration data output to: {0}'.format(self.calibrations_path))
        msgs.info('Science data output to: {0}'.format(self.science_path))
        msgs.info('Quality assessment plots output to: {0}'.format(self.qa_path))
        # TODO: Is anything written to the qa dir or only to qa/PNGs?
        # Should we have separate calibration and science QA
        # directories?
        # An html file wrapping them all too

#        # Instantiate Calibrations class
#        if self.spectrograph.pypeline in ['MultiSlit', 'Echelle']:
#            self.caliBrate \
#                = calibrations.MultiSlitCalibrations(self.fitstbl, self.par['calibrations'],
#                                                     self.spectrograph, self.calibrations_path,
#                                                     qadir=self.qa_path,
#                                                     reuse_masters=self.reuse_masters,
#                                                     show=self.show,
#                                                     slitspat_num=self.par['rdx']['slitspatnum'])
#        elif self.spectrograph.pypeline in ['IFU']:
#            self.caliBrate \
#                = calibrations.IFUCalibrations(self.fitstbl, self.par['calibrations'],
#                                               self.spectrograph,
#                                               self.calibrations_path,
#                                               qadir=self.qa_path,
#                                               reuse_masters=self.reuse_masters,
#                                               show=self.show)
#        else:
#            msgs.error("No calibration available to support pypeline: {0:s}".format(self.spectrograph.pypeline))

        # Init
        self.verbosity = verbosity
        # TODO: I don't think this ever used

        self.det = None

        self.tstart = None
        self.basename = None
        self.sciI = None
        self.obstime = None

    @property
    def science_path(self):
        """Return the path to the science directory."""
        return os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])

    @property
    def qa_path(self):
        """Return the path to the top-level QA directory."""
        return os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])

    def build_qa(self):
        """
        Generate QA wrappers
        """
        qa.gen_mf_html(self.pypeit_file, self.qa_path)
        qa.gen_exp_html()

    # TODO: This should go in a more relevant place
    def spec_output_file(self, frame, twod=False):
        """
        Return the path to the spectral output data file.
        
        Args:
            frame (:obj:`int`):
                Frame index from :attr:`fitstbl`.
            twod (:obj:`bool`):
                Name for the 2D output file; 1D file otherwise.
        
        Returns:
            :obj:`str`: The path for the output file
        """
        return os.path.join(self.science_path, 'spec{0}d_{1}.fits'.format('2' if twod else '1',
                                                    self.fitstbl.construct_basename(frame)))

    def outfile_exists(self, frame):
        """
        Check whether the 2D outfile of a given frame already exists

        Args:
            frame (int): Frame index from fitstbl

        Returns:
            bool: True if the 2d file exists, False if it does not exist
        """
        return os.path.isfile(self.spec_output_file(frame, twod=True))

    def get_std_outfile(self, standard_frames):
        """
        Grab the output filename from an input list of standard_frame indices

        If more than one index is provided, the first is taken

        Args:
            standard_frames (list): List of indices corresponding to standard stars

        Returns:
            str: Full path to the standard spec1d output file
        """
        # TODO: Need to decide how to associate standards with
        # science frames in the case where there is more than one
        # standard associated with a given science frame.  Below, I
        # just use the first standard

        std_outfile = None
        std_frame = None if len(standard_frames) == 0 else standard_frames[0]
        # Prepare to load up standard?
        if std_frame is not None:
            std_outfile = self.spec_output_file(std_frame) \
                            if isinstance(std_frame, (int,np.integer)) else None
        if std_outfile is not None and not os.path.isfile(std_outfile):
            msgs.error('Could not find standard file: {0}'.format(std_outfile))
        return std_outfile

    def calib_all(self):
        """
        Create calibrations for all setups

        This will not crash if not all of the standard set of files are not provided


        """

        self.tstart = time.time()

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))
        for i in range(self.fitstbl.n_calib_groups):
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(i)
            grp_frames = frame_indx[in_grp]

            # Find the detectors to reduce
            detectors = PypeIt.select_detectors(detnum=self.par['rdx']['detnum'],
                                            ndet=self.spectrograph.ndet)
            # Loop on Detectors
            for self.det in detectors:
                # Instantiate Calibrations class
                self.caliBrate = calibrations.Calibrations.get_instance(
                    self.fitstbl, self.par['calibrations'], self.spectrograph,
                    self.calibrations_path, qadir=self.qa_path, reuse_masters=self.reuse_masters,
                    show=self.show, slitspat_num=self.par['rdx']['slitspatnum'])
                # Do it
                self.caliBrate.set_config(grp_frames[0], self.det, self.par['calibrations'])
                self.caliBrate.run_the_steps()

        # Finish
        self.print_end_time()

    def reduce_all(self):
        """
        Main driver of the entire reduction

        Calibration and extraction via a series of calls to reduce_exposure()

        """
        # Validate the parameter set
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.time()

        # Find the standard frames
        is_standard = self.fitstbl.find_frames('standard')

        # Find the science frames
        is_science = self.fitstbl.find_frames('science')

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))

        # Iterate over each calibration group and reduce the standards
        for i in range(self.fitstbl.n_calib_groups):

            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(i)

            # Find the indices of the standard frames in this calibration group:
            grp_standards = frame_indx[is_standard & in_grp]

            # Reduce all the standard frames, loop on unique comb_id
            u_combid_std= np.unique(self.fitstbl['comb_id'][grp_standards])
            for j, comb_id in enumerate(u_combid_std):
                frames = np.where(self.fitstbl['comb_id'] == comb_id)[0]
                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]
                if not self.outfile_exists(frames[0]) or self.overwrite:
                    std_spec2d, std_sobjs = self.reduce_exposure(frames, bg_frames=bg_frames)
                    # TODO come up with sensible naming convention for save_exposure for combined files
                    self.save_exposure(frames[0], std_spec2d, std_sobjs, self.basename)
                else:
                    msgs.info('Output file: {:s} already exists'.format(self.fitstbl.construct_basename(frames[0])) +
                              '. Set overwrite=True to recreate and overwrite.')

        # Iterate over each calibration group again and reduce the science frames
        for i in range(self.fitstbl.n_calib_groups):
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(i)

            # Find the indices of the science frames in this calibration group:
            grp_science = frame_indx[is_science & in_grp]
            # Associate standards (previously reduced above) for this setup
            std_outfile = self.get_std_outfile(frame_indx[is_standard])
            # Reduce all the science frames; keep the basenames of the science frames for use in flux calibration
            science_basename = [None]*len(grp_science)
            # Loop on unique comb_id
            u_combid = np.unique(self.fitstbl['comb_id'][grp_science])
            for j, comb_id in enumerate(u_combid):
                frames = np.where(self.fitstbl['comb_id'] == comb_id)[0]
                # Find all frames whose comb_id matches the current frames bkg_id.
                bg_frames = np.where((self.fitstbl['comb_id'] == self.fitstbl['bkg_id'][frames][0]) &
                                     (self.fitstbl['comb_id'] >= 0))[0]
                # JFH changed the syntax below to that above, which allows frames to be used more than once
                # as a background image. The syntax below would require that we could somehow list multiple
                # numbers for the bkg_id which is impossible without a comma separated list
#                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]
                if not self.outfile_exists(frames[0]) or self.overwrite:
                    # TODO -- Should we reset/regenerate self.slits.mask for a new exposure
                    sci_spec2d, sci_sobjs = self.reduce_exposure(frames, bg_frames=bg_frames,
                                                    std_outfile=std_outfile)
                    science_basename[j] = self.basename
                    # TODO come up with sensible naming convention for save_exposure for combined files
                    self.save_exposure(frames[0], sci_spec2d, sci_sobjs, self.basename)
                else:
                    msgs.warn('Output file: {:s} already exists'.format(self.fitstbl.construct_basename(frames[0])) +
                              '. Set overwrite=True to recreate and overwrite.')

            msgs.info('Finished calibration group {0}'.format(i))

        # Check if this is an IFU reduction. If so, make a datacube
        if self.spectrograph.pypeline == "IFU" and self.par['reduce']['cube']['make_cube']:
            msgs.work("Generate datacube")

        # Finish
        self.print_end_time()

    # This is a static method to allow for use in coadding script 
    @staticmethod
    def select_detectors(detnum=None, ndet=1, slitspatnum=None):
        """
        Return the 1-indexed list of detectors to reduce.

        Args:
            detnum (:obj:`int`, :obj:`list`, optional):
                One or more detectors to reduce.  If None, return the
                full list for the provided number of detectors (`ndet`).
            ndet (:obj:`int`, optional):
                The number of detectors for this instrument.  Only used
                if `detnum is None`.

        Returns:
            list:  List of detectors to be reduced

        """
        if detnum is not None and slitspatnum is not None:
            msgs.error("You cannot specify both detnum and slitspatnum.  Too painful for over-writing SpecObjs")
        if detnum is None and slitspatnum is None:
            return np.arange(1, ndet+1).tolist()
        elif detnum is not None:
            return np.atleast_1d(detnum).tolist()
        else:
            return slittrace.parse_slitspatnum(slitspatnum)[0].tolist()

    def reduce_exposure(self, frames, bg_frames=None, std_outfile=None):
        """
        Reduce a single exposure

        Args:
            frame (:obj:`int`):
                0-indexed row in :attr:`fitstbl` with the frame to
                reduce.
            bg_frames (:obj:`list`, optional):
                List of frame indices for the background.
            std_outfile (:obj:`str`, optional):
                File with a previously reduced standard spectrum from
                PypeIt.

        Returns:
            dict: The dictionary containing the primary outputs of
            extraction.

        """

        # TODO:
        # - change doc string to reflect that more than one frame can be
        #   provided

        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            # TODO: Put this in a try/except block?
            display.clear_all()

        has_bg = True if bg_frames is not None and len(bg_frames) > 0 else False

        # Is this an IR reduction?
        # TODO: Why specific to IR?
        self.ir_redux = True if has_bg else False

        # Container for all the Spec2DObj
        all_spec2d = spec2dobj.AllSpec2DObj()
        all_spec2d['meta']['ir_redux'] = self.ir_redux

        # TODO -- Should we reset/regenerate self.slits.mask for a new exposure

        all_specobjs = specobjs.SpecObjs()

        # Print status message
        msgs_string = 'Reducing target {:s}'.format(self.fitstbl['target'][frames[0]]) + msgs.newline()
        # TODO: Print these when the frames are actually combined,
        # backgrounds are used, etc?
        msgs_string += 'Combining frames:' + msgs.newline()
        for iframe in frames:
            msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
        msgs.info(msgs_string)
        if has_bg:
            bg_msgs_string = ''
            for iframe in bg_frames:
                bg_msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
            bg_msgs_string = msgs.newline() + 'Using background from frames:' + msgs.newline() + bg_msgs_string
            msgs.info(bg_msgs_string)

        # Find the detectors to reduce
        detectors = PypeIt.select_detectors(detnum=self.par['rdx']['detnum'],
                                            slitspatnum=self.par['rdx']['slitspatnum'],
                                            ndet=self.spectrograph.ndet)
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([ str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Loop on Detectors
        # TODO: Attempt to put in a multiprocessing call here?
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            # Instantiate Calibrations class
            self.caliBrate = calibrations.Calibrations.get_instance(
                self.fitstbl, self.par['calibrations'], self.spectrograph,
                self.calibrations_path, qadir=self.qa_path, reuse_masters=self.reuse_masters,
                show=self.show, slitspat_num=self.par['rdx']['slitspatnum'])
            # These need to be separate to accomodate COADD2D
            self.caliBrate.set_config(frames[0], self.det, self.par['calibrations'])
            self.caliBrate.run_the_steps()
            # Extract
            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            all_spec2d[self.det], tmp_sobjs \
                    = self.reduce_one(frames, self.det, bg_frames, std_outfile=std_outfile)
            # Hold em
            if tmp_sobjs.nobj > 0:
                all_specobjs.add_sobj(tmp_sobjs)
            # JFH TODO write out the background frame?

            # TODO -- Save here?  Seems like we should.  Would probably need to use update_det=True

        # Return
        return all_spec2d, all_specobjs

    def get_sci_metadata(self, frame, det):
        """
        Grab the meta data for a given science frame and specific detector

        Args:
            frame (int): Frame index
            det (int): Detector index

        Returns:
            5 objects are returned::
                - str: Object type;  science or standard
                - str: Setup string from master_key()
                - astropy.time.Time: Time of observation
                - str: Basename of the frame
                - str: Binning of the detector

        """

        # Set binning, obstime, basename, and objtype
        binning = self.fitstbl['binning'][frame]
        obstime  = self.fitstbl.construct_obstime(frame)
        basename = self.fitstbl.construct_basename(frame, obstime=obstime)
        objtype  = self.fitstbl['frametype'][frame]
        if 'science' in objtype:
            objtype_out = 'science'
        elif 'standard' in objtype:
            objtype_out = 'standard'
        else:
            msgs.error('Unrecognized objtype')
        setup = self.fitstbl.master_key(frame, det=det)
        return objtype_out, setup, obstime, basename, binning

    def get_std_trace(self, std_redux, det, std_outfile):
        """
        Returns the trace of the standard if it is applicable to the current reduction

        Args:
            std_redux (bool): If False, proceed
            det (int): Detector index
            std_outfile (str): Filename for the standard star spec1d file

        Returns:
            ndarray or None: Trace of the standard star on input detector

        """
        if std_redux is False and std_outfile is not None:
            sobjs = specobjs.SpecObjs.from_fitsfile(std_outfile)
            # Does the detector match?
            # TODO Instrument specific logic here could be implemented with the parset. For example LRIS-B or LRIS-R we
            # we would use the standard from another detector
            this_det = sobjs.DET == det
            if np.any(this_det):
                sobjs_det = sobjs[this_det]
                sobjs_std = sobjs_det.get_std()
                # No standard extracted on this detector??
                if sobjs_std is None:
                    return None
                std_trace = sobjs_std.TRACE_SPAT
                # flatten the array if this multislit
                if 'MultiSlit' in self.spectrograph.pypeline:
                    std_trace = std_trace.flatten()
                elif 'Echelle' in self.spectrograph.pypeline:
                    std_trace = std_trace.T
                else:
                    msgs.error('Unrecognized pypeline')
            else:
                std_trace = None
        else:
            std_trace = None

        return std_trace

    def reduce_one(self, frames, det, bg_frames, std_outfile=None):
        """
        Reduce + Extract a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Args:
            frames (:obj:`list`):
                List of frames to extract; stacked if more than one
                is provided
            det (:obj:`int`):
                Detector number (1-indexed)
            bg_frames (:obj:`list`):
                List of frames to use as the background. Can be
                empty.
            std_outfile (:obj:`str`, optional):
                Filename for the standard star spec1d file. Passed
                directly to :func:`get_std_trace`.

        Returns:
            tuple: Returns six `numpy.ndarray`_ objects and a
            :class:`pypeit.specobjs.SpecObjs` object with the
            extracted spectra from this exposure/detector pair. The
            six `numpy.ndarray`_ objects are (1) the science image,
            (2) its inverse variance, (3) the sky model, (4) the
            object model, (5) the model inverse variance, and (6) the
            mask.

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning \
                = self.get_sci_metadata(frames[0], det)
        msgs.info("Extraction begins for {} on det={}".format(self.basename, det))
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype
        if self.std_redux:
            frame_par = self.par['calibrations']['standardframe']
        else:
            frame_par = self.par['scienceframe']
        # Get the standard trace if need be
        std_trace = self.get_std_trace(self.std_redux, det, std_outfile)

        # Build Science image
        sci_files = self.fitstbl.frame_paths(frames)
        sciImg = buildimage.buildimage_fromlist(
            self.spectrograph, det, frame_par,
            sci_files, bias=self.caliBrate.msbias, bpm=self.caliBrate.msbpm,
            dark=self.caliBrate.msdark,
            flatimages=self.caliBrate.flatimages,
            slits=self.caliBrate.slits,  # For flexure correction
            ignore_saturation=False)

        # Background Image?
        if len(bg_frames) > 0:
            bg_file_list = self.fitstbl.frame_paths(bg_frames)
            sciImg = sciImg.sub(
                buildimage.buildimage_fromlist(
                self.spectrograph, det, frame_par,bg_file_list,
                bpm=self.caliBrate.msbpm, bias=self.caliBrate.msbias,
                dark=self.caliBrate.msdark,
                flatimages=self.caliBrate.flatimages,
                slits=self.caliBrate.slits,  # For flexure correction
                ignore_saturation=False), frame_par['process'])

        # Instantiate Reduce object
        # Required for pypeline specific object
        # At instantiaton, the fullmask in self.sciImg is modified
        self.redux = reduce.Reduce.get_instance(sciImg, self.spectrograph,
                                                self.par, self.caliBrate,
                                                self.objtype,
                                                ir_redux=self.ir_redux,
                                                std_redux=self.std_redux,
                                                setup=self.setup,
                                                show=self.show,
                                                det=det, binning=self.binning,
                                                std_outfile=std_outfile,
                                                basename=self.basename)
        # Show?
        if self.show:
            self.redux.show('image', image=sciImg.image, chname='processed',
                            slits=True, clear=True)

        # Do it
        skymodel, objmodel, ivarmodel, outmask, sobjs, scaleImg, waveImg, tilts = self.redux.run(
            std_trace=std_trace, show_peaks=self.show,
            ra=self.fitstbl["ra"][frames[0]], dec=self.fitstbl["dec"][frames[0]],
            obstime=self.obstime)

        # TODO -- Save the slits yet again?


        # TODO -- Do this upstream
        # Tack on detector
        for sobj in sobjs:
            sobj.DETECTOR = sciImg.detector

        # Construct the Spec2DObj
        spec2DObj = spec2dobj.Spec2DObj(det=self.det,
                                        sciimg=sciImg.image,
                                        ivarraw=sciImg.ivar,
                                        skymodel=skymodel,
                                        objmodel=objmodel,
                                        ivarmodel=ivarmodel,
                                        scaleimg=scaleImg,
                                        waveimg=waveImg,
                                        bpmmask=outmask,
                                        detector=sciImg.detector,
                                        sci_spat_flexure=sciImg.spat_flexure,
                                        tilts=tilts,
                                        slits=copy.deepcopy(self.caliBrate.slits))
        spec2DObj.process_steps = sciImg.process_steps

        # Return
        return spec2DObj, sobjs

    def save_exposure(self, frame, all_spec2d, all_specobjs, basename):
        """
        Save the outputs from extraction for a given exposure

        Args:
            frame (:obj:`int`):
                0-indexed row in the metadata table with the frame
                that has been reduced.
            all_spec2d(:class:`pypeit.spec2dobj.AllSpec2DObj`):
            sci_dict (:obj:`dict`):
                Dictionary containing the primary outputs of
                extraction
            basename (:obj:`str`):
                The root name for the output file.

        Returns:
            None or SpecObjs:  All of the objects saved to disk

        """
        # TODO: Need some checks here that the exposure has been reduced?

        # Determine the headers
        row_fitstbl = self.fitstbl[frame]
        # Need raw file header information
        rawfile = self.fitstbl.frame_paths(frame)
        head2d = fits.getheader(rawfile, ext=self.spectrograph.primary_hdrext)

        # Check for the directory
        if not os.path.isdir(self.science_path):
            os.makedirs(self.science_path)

        subheader = self.spectrograph.subheader_for_spec(row_fitstbl, head2d)
        # 1D spectra
        if all_specobjs.nobj > 0:
            # Spectra
            outfile1d = os.path.join(self.science_path, 'spec1d_{:s}.fits'.format(basename))
            all_specobjs.write_to_fits(subheader, outfile1d,
                                       update_det=self.par['rdx']['detnum'],
                                       slitspatnum=self.par['rdx']['slitspatnum'])
            # Info
            outfiletxt = os.path.join(self.science_path, 'spec1d_{:s}.txt'.format(basename))
            all_specobjs.write_info(outfiletxt, self.spectrograph.pypeline)

        # 2D spectra
        outfile2d = os.path.join(self.science_path, 'spec2d_{:s}.fits'.format(basename))
        # Build header
        pri_hdr = all_spec2d.build_primary_hdr(head2d, self.spectrograph,
                                               redux_path=self.par['rdx']['redux_path'],
                                               master_key_dict=self.caliBrate.master_key_dict,
                                               master_dir=self.caliBrate.master_dir,
                                               subheader=subheader)
        # Write
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr, update_det=self.par['rdx']['detnum'])


    def msgs_reset(self):
        """
        Reset the msgs object
        """

        # Reset the global logger
        msgs.reset(log=self.logname, verbosity=self.verbosity)
        msgs.pypeit_file = self.pypeit_file

    def print_end_time(self):
        """
        Print the elapsed time
        """
        # Capture the end time and print it to user
        tend = time.time()
        codetime = tend-self.tstart
        if codetime < 60.0:
            msgs.info('Execution time: {0:.2f}s'.format(codetime))
        elif codetime/60.0 < 60.0:
            mns = int(codetime/60.0)
            scs = codetime - 60.0*mns
            msgs.info('Execution time: {0:d}m {1:.2f}s'.format(mns, scs))
        else:
            hrs = int(codetime/3600.0)
            mns = int(60.0*(codetime/3600.0 - hrs))
            scs = codetime - 60.0*mns - 3600.0*hrs
            msgs.info('Execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs))

    # TODO: Move this to fitstbl?
    def show_science(self):
        """
        Simple print of science frames
        """
        indx = self.fitstbl.find_frames('science')
        print(self.fitstbl[['target','ra','dec','exptime','dispname']][indx])

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


