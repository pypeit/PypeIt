"""
Main driver class for PypeIt run
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
#from abc import ABCMeta
import os
import datetime
import numpy as np
from collections import OrderedDict

from astropy.io import fits
from pypeit import msgs
from pypeit import calibrations
from pypeit import scienceimage
from pypeit import specobjs
from pypeit import ginga
from pypeit import reduce
from pypeit.core import paths
from pypeit.core import qa
from pypeit.core import wave
from pypeit.core import save
from pypeit.core import load
from pypeit.spectrographs.util import load_spectrograph
from linetools import utils as ltu



from configobj import ConfigObj
from pypeit.par.util import parse_pypeit_file
from pypeit.par import PypeItPar
from pypeit.metadata import PypeItMetaData

from pypeit import debugger

class PypeIt(object):
    """
    This class runs the primary calibration and extraction in PypeIt

    Args:
        pypeit_file (:obj:`str`):  PypeIt filename
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be::
                - 0: No output
                - 1: Minimal output (default)
                - 2: All output
        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        reuse_masters (bool, optional): Reuse any pre-existing calibration files
        logname (:obj:`str`, optional):
            The name of an ascii log file with the details of the
            reduction.
        redux_path (:obj:`str`, optional):
            Over-ride reduction path in PypeIt file (e.g. Notebook usage)
        show: (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via "ginga --modules=RC &"

    Attributes:
        pypeit_file (:obj:`str`):
            Name of the pypeit file to read.  PypeIt files have a specific
            set of valid formats. A description can be found `here`_
            (include doc link).
        fitstbl (:obj:`pypit.metadata.PypeItMetaData`): holds the meta info
    """
#    __metaclass__ = ABCMeta

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_masters=False, logname=None,
                 show=False, redux_path=None):

        # Load
        cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(pypeit_file, runtime=True)
        self.pypeit_file = pypeit_file

        # Spectrograph
        cfg = ConfigObj(cfg_lines)
        spectrograph_name = cfg['rdx']['spectrograph']
        self.spectrograph = load_spectrograph(spectrograph_name)

        # Par
        # Defaults
        spectrograph_def_par = self.spectrograph.default_pypeit_par()
        # Grab a science file for configuration specific parameters
        for idx, row in enumerate(usrdata):
            if 'science' in row['frametype']:
                sci_file = data_files[idx]
                break
        # Set
        spectrograph_cfg_lines = self.spectrograph.config_specific_par(spectrograph_def_par, sci_file).to_config()
        self.par = PypeItPar.from_cfg_lines(cfg_lines=spectrograph_cfg_lines, merge_with=cfg_lines)

        # Fitstbl
        self.fitstbl = PypeItMetaData(self.spectrograph, self.par, file_list=data_files,
                                      usrdata=usrdata, strict=True)

        # The following could be put in a prepare_to_run() method in PypeItMetaData
        if 'setup' not in self.fitstbl.keys():
            self.fitstbl['setup'] = setups[0]
        self.fitstbl.get_frame_types(user=frametype)  # This sets them using the user inputs
        self.fitstbl.set_defaults()  # Only does something if values not set in PypeIt file
        self.fitstbl._set_calib_group_bits()
        self.fitstbl._check_calib_groups()
        # Write .calib file (For QA naming amongst other things)
        calib_file = pypeit_file.replace('.pypeit', '.calib')
        self.fitstbl.write_calib(calib_file)

        # Other Internals
        self.logname = logname
        self.overwrite = overwrite
        # Currently the runtime argument determines the behavior for reuse_masters. There is also a reuse_masters
        # parameter in the parset but it is currently ignored.
        self.reuse_masters=reuse_masters
        self.show = show

        # Make the output directories
        self.par['rdx']['redux_path'] = os.getcwd() if redux_path is None else redux_path
        msgs.info("Setting reduction path to {:s}".format(self.par['rdx']['redux_path']))
        paths.make_dirs(self.spectrograph.spectrograph, self.par['calibrations']['caldir'],
                        self.par['rdx']['scidir'], self.par['rdx']['qadir'],
                        overwrite=self.overwrite, redux_path=self.par['rdx']['redux_path'])

        # Instantiate Calibrations class
        self.caliBrate \
            = calibrations.MultiSlitCalibrations(self.fitstbl, self.par['calibrations'], self.spectrograph,
                                                 redux_path=self.par['rdx']['redux_path'],
                                                 reuse_masters=self.reuse_masters,
                                                 save_masters=True, write_qa=True,
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

    def build_qa(self):
        """
        Generate QA wrappers
        """
        qa.gen_mf_html(self.pypeit_file)
        qa.gen_exp_html()

    def outfile_exists(self, frame):
        """
        Check whether the 2D outfile of a given frame already exists

        Args:
            frame (int): Frame index from fitstbl

        Returns:
            bool: True if the 2d file exists
                 False if it does not exist
        """
        # Check if the 2d output file exists
        scidir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])
        basename = self.fitstbl.construct_basename(frame)
        outfile = scidir + '/spec2d_{:s}.fits'.format(basename)
        return os.path.isfile(outfile)

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
            std_outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],
            'spec1d_{:s}.fits'.format(self.fitstbl.construct_basename(std_frame))) \
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
                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]
                if not self.outfile_exists(frames[0]) or self.overwrite:
                    sci_dict = self.reduce_exposure(frames, bg_frames=bg_frames, std_outfile=std_outfile)
                    science_basename[j] = self.basename
                    # TODO come up with sensible naming convention for save_exposure for combined files
                    self.save_exposure(frames[0], sci_dict, self.basename)
                else:
                    msgs.info('Output file: {:s} already exists'.format(self.fitstbl.construct_basename(frames[0])) +
                              '. Set overwrite=True to recreate and overwrite.')

            msgs.info('Finished calibration group {0}'.format(i))

        # Finish
        self.print_end_time()


    def select_detectors(self):
        """
        Return the 1-indexed list of detectors to reduce.

        Returns:
            list:  List of detectors to be reduced

        """
        if self.par['rdx']['detnum'] is None:
            return np.arange(self.spectrograph.ndet)+1
        return [self.par['rdx']['detnum']] if isinstance(self.par['rdx']['detnum'], int) \
                    else self.par['rdx']['detnum']

    def reduce_exposure(self, frames, bg_frames=[], std_outfile=None):
        """
        Reduce a single exposure

        Args:
            frame (:obj:`int`):
                0-indexed row in :attr:`fitstbl` with the frame to
                reduce
            bgframes (:obj:`list`, optional):
                List of frame indices for the background
            std_outfile (:obj:`str`, optional):
                the name of a file with a previously PypeIt-reduced standard spectrum.

        Returns:
            dict: The dictionary containing the primary outputs of extraction
        """
        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            ginga.clear_all()

        # Save the frame
        self.frames = frames
        self.bg_frames = bg_frames
        # Is this an IR reduction?
        self.ir_redux = True if len(bg_frames) > 0 else False

        # JFH Why does this need to be ordered?
        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        sci_dict['meta']['ir_redux'] = self.ir_redux

        # Print status message
        msgs_string = 'Reducing target {:s}'.format(self.fitstbl['target'][self.frames[0]]) + msgs.newline()
        msgs_string += 'Combining frames:' + msgs.newline()
        for iframe in self.frames:
            msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
        msgs.info(msgs_string)
        if len(bg_frames) > 0:
            bg_msgs_string = ''
            for iframe in self.bg_frames:
                bg_msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
            bg_msgs_string = msgs.newline() + 'Using background from frames:' + msgs.newline() + bg_msgs_string
            msgs.info(bg_msgs_string)

        # Find the detectors to reduce
        detectors = self.select_detectors()
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([ str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Loop on Detectors
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            sci_dict[self.det] = {}
            # Calibrate
            #TODO Is the right behavior to just use the first frame?
            self.caliBrate.set_config(self.frames[0], self.det, self.par['calibrations'])
            self.caliBrate.run_the_steps()
            # Extract
            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            sci_dict[self.det]['sciimg'], sci_dict[self.det]['sciivar'], sci_dict[self.det]['skymodel'], \
                sci_dict[self.det]['objmodel'], sci_dict[self.det]['ivarmodel'], sci_dict[self.det]['outmask'], \
                sci_dict[self.det]['specobjs'], vel_corr \
                    = self.extract_one(self.frames, self.det, bg_frames = self.bg_frames, std_outfile = std_outfile)
            if vel_corr is not None:
                sci_dict['meta']['vel_corr'] = vel_corr

            # JFH TODO write out the background frame?

        # Return
        return sci_dict

    def flexure_correct(self, sobjs, maskslits):
        """
        Correct for flexure

        Spectra are modified in place (wavelengths are shifted)

        Args:
            sobjs (SpecObjs):
            maskslits (ndarray): Mask of SpecObjs

        """

        if self.par['flexure']['method'] != 'skip':
            flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                         self.par['flexure']['spectrum'],
                                         mxshft=self.par['flexure']['maxshift'])
            # QA
            wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,
                            out_dir=self.par['rdx']['redux_path'])
        else:
            msgs.info('Skipping flexure correction.')

    def helio_correct(self, sobjs, maskslits, frame, obstime):
        """
        Perform a heliocentric correction on a set of spectra

        Args:
            sobjs (pypeit.specobjs.SpecObjs): Spectra
            maskslits (ndarray): Slits that are masked
            frame (int): Frame to use for meta info
            obstime (astropy.time.Time):

        Returns:
            astropy.units.Quantity: Velocity correction in km/s

        """
        # Helio, correct Earth's motion
        if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                and (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
            # TODO change this keyword to refframe instead of frame
            msgs.info("Performing a {0} correction".format(self.caliBrate.par['wavelengths']['frame']))
            vel, vel_corr = wave.geomotion_correct(sobjs, maskslits, self.fitstbl, frame, obstime,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   self.caliBrate.par['wavelengths']['frame'])
        else:
            msgs.info('A wavelength reference-frame correction will not be performed.')
            vel_corr = None

        return vel_corr

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
            sobjs, hdr_std = load.load_specobjs(std_outfile)
            # Does the detector match?
            # TODO Instrument specific logic here could be implemented with the parset. For example LRIS-B or LRIS-R we
            # we would use the standard from another detector
            this_det = sobjs.det == det
            if np.any(this_det):
                sobjs_det = sobjs[this_det]
                sobjs_std = sobjs_det.get_std()
                std_trace = sobjs_std.trace_spat
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

    def extract_one(self, frames, det, bg_frames=[], std_outfile=None):
        """
        Extract a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Args:
            frames (list):  List of frames to extract;  stacked if more than one is provided
            det (int):
            bg_frames (list, optional): List of frames to use as the background
            std_outfile (str, optional):

        Returns:
            eight objects are returned::
                - ndarray: Science image
                - ndarray: Science inverse variance image
                - ndarray: Model of the sky
                - ndarray: Model of the object
                - ndarray: Model of inverse variance
                - ndarray: Mask
                - :obj:`pypeit.specobjs.SpecObjs`: spectra
                - astropy.units.Quantity: velocity correction

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning = self.get_sci_metadata(frames[0], det)
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype
        # Get the standard trace if need be
        std_trace = self.get_std_trace(self.std_redux, det, std_outfile)
        # Instantiate ScienceImage for the files we will reduce
        sci_files = self.fitstbl.frame_paths(frames)
        self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_files,
                                              bg_file_list=self.fitstbl.frame_paths(bg_frames),
                                              ir_redux = self.ir_redux,
                                              par=self.par['scienceframe'],
                                              det=det,
                                              binning=self.binning)
        # For QA on crash.
        msgs.sciexp = self.sciI

        # Process images (includes inverse variance image, rn2 image, and CR mask)
        self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = \
            self.sciI.proc(self.caliBrate.msbias, self.caliBrate.mspixflatnrm.copy(),
                           self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                           show=self.show)
        # Object finding, first pass on frame without sky subtraction
        self.maskslits = self.caliBrate.maskslits.copy()

        self.redux = reduce.instantiate_me(self.spectrograph, self.caliBrate.tslits_dict,
                                           self.mask, self.par,
                                           ir_redux = self.ir_redux,
                                           objtype=self.objtype, setup=self.setup,
                                           det=det, binning=self.binning)

        # Prep for manual extraction (if requested)
        manual_extract_dict = self.fitstbl.get_manual_extract(frames, det)

        # Do one iteration of object finding, and sky subtract to get initial sky model
        self.sobjs_obj, self.nobj, skymask_init = \
            self.redux.find_objects(self.sciimg, self.sciivar, std=self.std_redux, ir_redux=self.ir_redux,
                                    std_trace=std_trace,maskslits=self.maskslits,
                                    show=self.show & (not self.std_redux),
                                    manual_extract_dict=manual_extract_dict)

        # Global sky subtraction, first pass. Uses skymask from object finding step above
        self.initial_sky = \
            self.redux.global_skysub(self.sciimg, self.sciivar, self.caliBrate.tilts_dict['tilts'], skymask=skymask_init,
                                    std=self.std_redux, maskslits=self.maskslits, show=self.show)

        if not self.std_redux:
            # Object finding, second pass on frame *with* sky subtraction. Show here if requested
            self.sobjs_obj, self.nobj, self.skymask = \
                self.redux.find_objects(self.sciimg - self.initial_sky, self.sciivar, std=self.std_redux, ir_redux=self.ir_redux,
                                  std_trace=std_trace,maskslits=self.maskslits,show=self.show,
                                        manual_extract_dict=manual_extract_dict)

        # If there are objects, do 2nd round of global_skysub, local_skysub_extract, flexure, geo_motion
        if self.nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            self.global_sky = self.initial_sky if self.std_redux else \
                self.redux.global_skysub(self.sciimg, self.sciivar, self.caliBrate.tilts_dict['tilts'],
                skymask=self.skymask, maskslits=self.maskslits, show=self.show)

            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
            self.redux.local_skysub_extract(self.sciimg, self.sciivar, self.caliBrate.tilts_dict['tilts'], self.caliBrate.mswave,
                                            self.global_sky, self.rn2img, self.sobjs_obj,
                                            model_noise=(not self.ir_redux),std = self.std_redux,
                                            maskslits=self.maskslits, show_profile=self.show,show=self.show)

            # Purge out the negative objects if this was a near-IR reduction.
            # TODO should we move this purge call to local_skysub_extract??
            if self.ir_redux:
                self.sobjs.purge_neg()

            # Flexure correction if this is not a standard star
            if not self.std_redux:
                self.redux.flexure_correct(self.sobjs, self.basename)

            # Grab coord
            radec = ltu.radec_to_coord((self.fitstbl["ra"][frames[0]], self.fitstbl["dec"][frames[0]]))
            self.vel_corr = self.redux.helio_correct(self.sobjs, radec, self.obstime)

        else:
            # Print status message
            msgs_string = 'No objects to extract for target {:s}'.format(self.fitstbl['target'][frames[0]]) + msgs.newline()
            msgs_string += 'On frames:' + msgs.newline()
            for iframe in frames:
                msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
            msgs.warn(msgs_string)
            # set to first pass global sky
            self.skymodel = self.initial_sky
            self.objmodel = np.zeros_like(self.sciimg)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciivar)
            # Set to inmask in case on objects were found
            self.outmask = self.mask
            # empty specobjs object from object finding
            if self.ir_redux:
                self.sobjs_obj.purge_neg()
            self.sobjs = self.sobjs_obj
            self.vel_corr = None

        return self.sciimg, self.sciivar, self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs, self.vel_corr

    # TODO: Why not use self.frame?
    def save_exposure(self, frame, sci_dict, basename, only_1d=False):
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
            only_1d (:obj:`bool`, optional):
              Save only the 1D spectra?

        Returns:
            None or SpecObjs:  All of the objects saved to disk

        """
        # TODO: Need some checks here that the exposure has been reduced

        # Determine the headers
        head1d = self.fitstbl[frame]
        # Need raw file header information
        rawfile = self.fitstbl.frame_paths(frame)
        head2d = fits.getheader(rawfile, ext=self.spectrograph.primary_hdrext,)
        refframe = 'pixel' if self.caliBrate.par['wavelengths']['reference'] == 'pixel' else \
            self.caliBrate.par['wavelengths']['frame']

        # Determine the paths/filenames
        scipath = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])

        save.save_all(sci_dict, self.caliBrate.master_key_dict, self.caliBrate.master_dir, self.spectrograph,
                      head1d, head2d, scipath, basename, only_1d=only_1d, refframe=refframe,
                      update_det=self.par['rdx']['detnum'], binning=self.fitstbl['binning'][frame])

        return



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


