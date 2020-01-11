"""
Main driver class for PypeIt run
"""
import time
import os
import numpy as np
from collections import OrderedDict
from astropy.io import fits
from pypeit import msgs
from pypeit import calibrations
from pypeit.images import scienceimage
from pypeit import ginga
from pypeit import reduce
from pypeit.core import qa
from pypeit.core import wave
from pypeit.core import save
from pypeit import specobjs
from pypeit.core import pixels
from pypeit.spectrographs.util import load_spectrograph

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

    Attributes:
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  PypeIt files have a
            specific set of valid formats. A description can be found
            :ref:`pypeit_file`.
        fitstbl (:obj:`pypit.metadata.PypeItMetaData`): holds the meta info

    """
#    __metaclass__ = ABCMeta

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_masters=False, logname=None,
                 show=False, redux_path=None):

        # Load
        cfg_lines, data_files, frametype, usrdata, setups \
                = parse_pypeit_file(pypeit_file, runtime=True)
        self.pypeit_file = pypeit_file

        # Spectrograph
        cfg = ConfigObj(cfg_lines)
        spectrograph_name = cfg['rdx']['spectrograph']
        self.spectrograph = load_spectrograph(spectrograph_name, ifile=data_files[0])
        msgs.info('Loaded spectrograph {0}'.format(self.spectrograph.spectrograph))

        # --------------------------------------------------------------
        # Get the full set of PypeIt parameters
        #   - Grab a science or standard file for configuration specific parameters
        scistd_file = None
        for idx, row in enumerate(usrdata):
            if ('science' in row['frametype']) or ('standard' in row['frametype']):
                scistd_file = data_files[idx]
                break
        #   - Configuration specific parameters for the spectrograph
        if scistd_file is not None:
            msgs.info('Setting configuration-specific parameters using {0}'.format(
                      os.path.split(scistd_file)[1]))
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(scistd_file).to_config()
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
        if self.par['calibrations']['caldir'] == 'default':
            self.calibrations_path = os.path.join(self.par['rdx']['redux_path'], 'Masters')
        else:
            self.calibrations_path = self.par['calibrations']['caldir']

        # Report paths
        msgs.info('Setting reduction path to {0}'.format(self.par['rdx']['redux_path']))
        msgs.info('Master calibration data output to: {0}'.format(self.calibrations_path))
        msgs.info('Science data output to: {0}'.format(self.science_path))
        msgs.info('Quality assessment plots output to: {0}'.format(self.qa_path))
        # TODO: Is anything written to the qa dir or only to qa/PNGs?
        # Should we have separate calibration and science QA
        # directories?

        # Instantiate Calibrations class
        self.caliBrate \
            = calibrations.MultiSlitCalibrations(self.fitstbl, self.par['calibrations'],
                                                 self.spectrograph,
                                                 caldir=self.calibrations_path,
                                                 qadir=self.qa_path,
                                                 reuse_masters=self.reuse_masters,
                                                 show=self.show)
        # Init
        self.verbosity = verbosity
        # TODO: I don't think this ever used

        self.frame = None
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

    def reduce_all(self):
        """
        Main driver of the entire reduction

        Calibration and extraction via a series of calls to reduce_exposure()

        """
        # Validate the parameter set
        required = ['rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib']
        can_be_None = ['flexure', 'fluxcalib']
        self.par.validate_keys(required=required, can_be_None=can_be_None)

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
                    std_dict = self.reduce_exposure(frames, bg_frames=bg_frames)
                    # TODO come up with sensible naming convention for save_exposure for combined files
                    self.save_exposure(frames[0], std_dict, self.basename)
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
                    sci_dict = self.reduce_exposure(frames, bg_frames=bg_frames,
                                                    std_outfile=std_outfile)
                    science_basename[j] = self.basename
                    # TODO come up with sensible naming convention for save_exposure for combined files
                    self.save_exposure(frames[0], sci_dict, self.basename)
                else:
                    msgs.warn('Output file: {:s} already exists'.format(self.fitstbl.construct_basename(frames[0])) +
                              '. Set overwrite=True to recreate and overwrite.')

            msgs.info('Finished calibration group {0}'.format(i))

        # Finish
        self.print_end_time()

    # This is a static method to allow for use in coadding script 
    @staticmethod
    def select_detectors(detnum=None, ndet=1):
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
        if detnum is None:
            return np.arange(1, ndet+1).tolist()
        return [detnum] if isinstance(detnum, int) else detnum

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
        # - bg_frames should be None by default
        # - change doc string to reflect that more than one frame can be
        #   provided

        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            # TODO: Put this in a try/except block?
            ginga.clear_all()

        has_bg = True if bg_frames is not None and len(bg_frames) > 0 else False

        # Is this an IR reduction?
        # TODO: Why specific to IR?
        self.ir_redux = True if has_bg else False

        # TODO: JFH Why does this need to be ordered?
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['ir_redux'] = self.ir_redux

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
                                            ndet=self.spectrograph.ndet)
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([ str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Loop on Detectors
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            sci_dict[self.det] = {}
            # Calibrate
            #TODO Is the right behavior to just use the first frame?
            self.caliBrate.set_config(frames[0], self.det, self.par['calibrations'])
            self.caliBrate.run_the_steps()
            # Extract
            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            sci_dict[self.det]['sciimg'], sci_dict[self.det]['sciivar'], \
                sci_dict[self.det]['skymodel'], sci_dict[self.det]['objmodel'], \
                sci_dict[self.det]['ivarmodel'], sci_dict[self.det]['outmask'], \
                sci_dict[self.det]['specobjs'], \
                        = self.extract_one(frames, self.det, bg_frames,
                                           std_outfile=std_outfile)
            # JFH TODO write out the background frame?

        # Return
        return sci_dict

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
            ndarray: Trace of the standard star on input detector

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

    def extract_one(self, frames, det, bg_frames, std_outfile=None):
        """
        Extract a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Args:
            frames (list):
                List of frames to extract;  stacked if more than one is provided
            det (int):
            bg_frames (list):
                List of frames to use as the background
                Can be empty
            std_outfile (str, optional):

        Returns:
            seven objects are returned::
                - ndarray: Science image
                - ndarray: Science inverse variance image
                - ndarray: Model of the sky
                - ndarray: Model of the object
                - ndarray: Model of inverse variance
                - ndarray: Mask
                - :obj:`pypeit.specobjs.SpecObjs`: spectra

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning = self.get_sci_metadata(frames[0], det)
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype
        # Get the standard trace if need be
        std_trace = self.get_std_trace(self.std_redux, det, std_outfile)

        # Build Science image
        sci_files = self.fitstbl.frame_paths(frames)
        self.sciImg = scienceimage.build_from_file_list(
            self.spectrograph, det, self.par['scienceframe']['process'],
            self.caliBrate.msbpm, sci_files, self.caliBrate.msbias,
            self.caliBrate.mspixelflat, illum_flat=self.caliBrate.msillumflat)

        # Background Image?
        if len(bg_frames) > 0:
            bg_file_list = self.fitstbl.frame_paths(bg_frames)
            self.sciImg = self.sciImg - scienceimage.build_from_file_list(
                self.spectrograph, det, self.par['scienceframe']['process'],
                self.caliBrate.msbpm, bg_file_list, self.caliBrate.msbias,
                self.caliBrate.mspixelflat, illum_flat=self.caliBrate.msillumflat)

        # Update mask for slitmask
        slitmask = pixels.tslits2mask(self.caliBrate.tslits_dict)
        self.sciImg.update_mask_slitmask(slitmask)

        # For QA on crash
        msgs.sciexp = self.sciImg

        # Instantiate Reduce object
        self.maskslits = self.caliBrate.tslits_dict['maskslits'].copy()
        # Required for pypeline specific object
        # TODO -- caliBrate should be replaced by the ~3 primary Objects needed
        #   once we have the data models in place.
        self.redux = reduce.instantiate_me(self.sciImg, self.spectrograph,
                                           self.par, self.caliBrate,
                                           maskslits=self.maskslits,
                                           ir_redux=self.ir_redux,
                                           std_redux=self.std_redux,
                                           objtype=self.objtype,
                                           setup=self.setup,
                                           show=self.show,
                                           det=det, binning=self.binning)
        # Show?
        if self.show:
            self.redux.show('image', image=self.sciImg.image, chname='processed',
                            slits=True, clear=True)

        # Prep for manual extraction (if requested)
        manual_extract_dict = self.fitstbl.get_manual_extract(frames, det)

        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = self.redux.run(
            std_trace=std_trace, manual_extract_dict=manual_extract_dict, show_peaks=self.show,
            basename=self.basename, ra=self.fitstbl["ra"][frames[0]], dec=self.fitstbl["dec"][frames[0]],
            obstime=self.obstime)

        # Return
        return self.sciImg.image, self.sciImg.ivar, self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    # TODO: Why not use self.frame?
    def save_exposure(self, frame, sci_dict, basename):
        """
        Save the outputs from extraction for a given exposure

        Args:
            frame (:obj:`int`):
              0-indexed row in the metadata table with the frame that
              has been reduced.
            sci_dict (:obj:`dict`):
              Dictionary containing the primary outputs of extraction
            basename (:obj:`str`):
                The root name for the output file.

        Returns:
            None or SpecObjs:  All of the objects saved to disk

        """
        # TODO: Need some checks here that the exposure has been reduced

        # Determine the headers
        head1d = self.fitstbl[frame]
        # Need raw file header information
        rawfile = self.fitstbl.frame_paths(frame)
        head2d = fits.getheader(rawfile, ext=self.spectrograph.primary_hdrext)
        refframe = 'pixel' if self.caliBrate.par['wavelengths']['reference'] == 'pixel' else \
            self.caliBrate.par['wavelengths']['frame']

        # Determine the paths/filenames
        save.save_all(sci_dict, self.caliBrate.master_key_dict, self.caliBrate.master_dir,
                      self.spectrograph, head1d, head2d, self.science_path, basename,
                      update_det=self.par['rdx']['detnum'], binning=self.fitstbl['binning'][frame])

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


