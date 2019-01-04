from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import time
from abc import ABCMeta
import os
import datetime
import numpy as np
from collections import OrderedDict

from astropy.time import Time
from pypeit import msgs
from pypeit import pypeitsetup
from pypeit import calibrations
from pypeit import scienceimage
from pypeit import specobjs
from pypeit import fluxspec
from pypeit import ginga
from pypeit.core import paths
from pypeit.core import qa
#from pypeit.core import pypsetup
from pypeit.core import wave
from pypeit.core import save
from pypeit.core import load
from pypeit.spectrographs.util import load_spectrograph


class PypeIt(object):
    """
    This class is designed to run PypeIt

    .. todo::
        Improve docstring...

    Args:
        pypeit_file (:obj:`str`,
            Filename
        verbosity (:obj:`int`, optional):
            Verbosity level of system output.  Can be::
                - 0: No output
                - 1: Minimal output (default)
                - 2: All output
        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
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
    """
    __metaclass__ = ABCMeta

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_masters=False, logname=None, show=False,
                 redux_path=None):

        # Setup
        self.pypeit_file = pypeit_file
        ps = pypeitsetup.PypeItSetup.from_pypeit_file(self.pypeit_file)
        ps.run(setup_only=False)
        # Only need the parameters, spectrograph, and metadata for the remainder
        self.par = ps.par
        # self.spectrograph = ps.spectrograph
        self.fitstbl = ps.fitstbl

        self.pypeitSetup = ps

        # Other Internals
        self.logname = logname
        self.overwrite = overwrite
        # Currently the runtime argument determines the behavior for reuse_masters. There is also a reuse_masters
        # parameter in the parset but it is currently ignored.
        self.reuse_masters=reuse_masters
        self.show = show


        # Spectrometer class
        self.spectrograph = load_spectrograph(ps.spectrograph)

        # Make the output directories
        self.par['rdx']['redux_path'] = os.getcwd() if redux_path is None else redux_path
        msgs.info("Setting reduction path to {:s}".format(self.par['rdx']['redux_path']))
        paths.make_dirs(self.spectrograph.spectrograph, self.par['calibrations']['caldir'],
                        self.par['rdx']['scidir'], self.par['rdx']['qadir'],
                        overwrite=self.overwrite, redux_path=self.par['rdx']['redux_path'])

        # Instantiate Calibrations class
        self.caliBrate \
            = calibrations.MultiSlitCalibrations(self.fitstbl, spectrograph=self.spectrograph,
                                                 par=self.par['calibrations'],
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

        Returns:

        """
        qa.gen_mf_html(self.pypeit_file)
        qa.gen_exp_html()

    def outfile_exists(self, frame):
        """
         Returns: True if the 2d file exists
                 False if it does not exist

        """
        # Check if the 2d output file exists
        scidir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])
        basename = self.fitstbl.construct_basename(frame)
        outfile = scidir + '/spec2d_{:s}.fits'.format(basename)
        return os.path.isfile(outfile)

    def get_std_outfile(self, standard_frames):
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
        Reduce all of the science exposures
        Generate all needed calibration files

        Returns:

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

            # Apply the flux calibration for this calibration group
            # TODO: I don't think this function is written yet...
            #self.flux_calibrate()

            msgs.info('Finished calibration group {0}'.format(i))

        # Finish
        self.print_end_time()


    def select_detectors(self):
        """
        Return the 1-indexed list of detectors to reduce.
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
            bgframe (:obj:`int`, optional):
                0-indexed row in :attr:`fitstbl` with the matching background frame
            std_frame (:obj:`int`, :obj:`str`, optional):
                0-indexed row in :attr:`fitstbl` with a standard frame
                associated with the frame to reduce, or the name of a
                file with a previously PypeIt-reduced standard spectrum.

        Returns:
            dict: The dictionary containing the primary outputs of
            extraction
        """
        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            ginga.clear_all()

        # Save the frame
        self.frames = frames
        self.bg_frames = bg_frames

        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.

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

    def flexure_correct(self,sobjs,maskslits):
        """ Correct for flexure """

        if self.par['flexure']['method'] != 'skip':
            sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
            flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                         sky_spectrum, sky_file=sky_file,
                                         mxshft=self.par['flexure']['maxshift'])
            # QA
            wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,out_dir=self.par['rdx']['redux_path'])
        else:
            msgs.info('Skipping flexure correction.')


    def helio_correct(self, sobjs, maskslits, frame, obstime):
        """ Perform a heliocentric correction """
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

        # Set binning, obstime, basename, and objtype
        try:
            binning = self.fitstbl['binning'][frame]
        except:
            binning = None
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

        Returns:
            sciimg
            sciivar
            skymodel
            objmodel
            ivarmodel
            outmask
            sobjs
            vel_corr

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning = self.get_sci_metadata(frames[0], det)
        # Is this an IR reduction
        self.ir_redux = True if len(bg_frames) > 0 else False
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype
        # Get the standard trace if need be
        std_trace = self.get_std_trace(self.std_redux, det, std_outfile)
        # Instantiate ScienceImage for the files we will reduce
        self.sciI = scienceimage.ScienceImage(self.caliBrate.tslits_dict, self.spectrograph,
                                              self.fitstbl.frame_paths(frames),
                                              bg_file_list=self.fitstbl.frame_paths(bg_frames),
                                              ir_redux = self.ir_redux,
                                              par=self.par['scienceimage'],
                                              frame_par=self.par['scienceframe'],
                                              objtype=self.objtype,
                                              det=det,
                                              binning=self.binning,
                                              setup=self.setup)
        # For QA on crash
        msgs.sciexp = self.sciI

        # Process images (includes inverse variance image, rn2 image, and CR mask)
        self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = \
            self.sciI.proc(self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
                           self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                           show=self.show)
        # Object finding, first pass on frame without sky subtraction
        self.maskslits = self.caliBrate.maskslits.copy()
        # Do one iteration of object finding, and sky subtract to get initial sky model
        self.sobjs_obj, self.nobj, skymask_init = \
            self.find_objects(self.sciimg, std=self.std_redux, ir_redux=self.ir_redux,
                              std_trace=std_trace, snr_trim=False,maskslits=self.maskslits,
                              show = (not self.std_redux))

        # Global sky subtraction, first pass. Uses skymask from object finding step above
        self.initial_sky = \
            self.sciI.global_skysub(self.caliBrate.tilts_dict['tilts'], skymask=skymask_init,
                                    std=self.std_redux, maskslits=self.maskslits, show=self.show)

        if not self.std_redux:
            # Object finding, second pass on frame *with* sky subtraction. Show here if requested
            self.sobjs_obj, self.nobj, self.skymask = \
                self.find_objects(self.sciimg - self.initial_sky, std=self.std_redux, ir_redux=self.ir_redux,
                                  std_trace=std_trace, snr_trim=True,
                                  maskslits=self.maskslits,show=self.show)

        # If there are objects, do 2nd round of global_skysub, local_skysub_extract, flexure, geo_motion
        if self.nobj > 0:
            if not self.std_redux:
                # Global sky subtraction second pass. Uses skymask from object finding
                self.global_sky = self.sciI.global_skysub(self.caliBrate.tilts_dict['tilts'],
                                                     skymask=self.skymask, maskslits=self.maskslits, show=self.show)
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
                self.sciI.local_skysub_extract(self.sobjs_obj, self.caliBrate.mswave, model_noise=(not self.ir_redux),
                                               std = self.std_redux,maskslits=self.maskslits, show_profile=self.show,
                                               show=self.show)

            # Purge out the negative objects if this was a near-IR reduction
            if self.ir_redux:
                self.sobjs.purge_neg()

            # Flexure correction if this is not a standard star
            if not self.std_redux:
                self.flexure_correct(self.sobjs, self.maskslits)
            vel_corr = self.helio_correct(self.sobjs, self.maskslits, frames[0], self.obstime)

        else:
            # Print status message
            msgs_string = 'No objects to extract for target {:s}'.format(self.fitstbl['target'][frames[0]]) + msgs.newline()
            msgs_string += 'On frames:' + msgs.newline()
            for iframe in frames:
                msgs_string += '{0:s}'.format(self.fitstbl['filename'][iframe]) + msgs.newline()
            msgs.warn(msgs_string)
            # set to first pass global sky
            skymodel = self.initial_sky
            objmodel = np.zeros_like(sciimg)
            # Set to sciivar. Could create a model but what is the point?
            ivarmodel = np.copy(self.sciivar)
            # Set to inmask in case on objects were found
            outmask = self.mask
            # empty specobjs object from object finding
            sobjs = self.sobjs_obj
            vel_corr = None

        return self.sciimg, self.sciivar, self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs, vel_corr

    def find_objects(self, image, std=False, ir_redux=False, std_trace=None, snr_trim=False, maskslits=None,
                          show_peaks=False, show_fits=False, show_trace=False, show=False):
        """
        Dummy method for object finding. Overloaded by class specific object finding.

        Returns:

        """

        return None, None, None

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
            DOC!

            This can return None or all_specobjs

        """
        # TODO: Need some checks here that the exposure has been reduced

        # Build the final list of specobjs and vel_corr
        all_specobjs = specobjs.SpecObjs()

        vel_corr = 0.  # This will not be set for Standard stars, which is fine
        for key in sci_dict:
            if key in ['meta']:
                vel_corr = sci_dict['meta']['vel_corr']
                continue
            #
            try:
                all_specobjs.add_sobj(sci_dict[key]['specobjs'])
            except KeyError:  # No object extracted
                continue

        if len(all_specobjs) == 0:
            msgs.warn('No objects to save!')
            return

        # Write 1D spectra
        outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],'spec1d_{:s}.fits'.format(basename))
        helio_dict = dict(refframe='pixel' if self.caliBrate.par['wavelengths']['reference'] == 'pixel' else \
            self.caliBrate.par['wavelengths']['frame'],vel_correction=vel_corr)
        # Did the user re-run a single detector?
        save.save_1d_spectra_fits(all_specobjs, self.fitstbl[frame], self.spectrograph.pypeline, outfile,
                                  helio_dict=helio_dict, telescope=self.spectrograph.telescope,
                                  update_det=self.par['rdx']['detnum'])
        # 1D only?
        if only_1d:
            return
        # Obj info
        scidir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])
        save.save_obj_info(all_specobjs, self.spectrograph, basename, scidir, binning=self.fitstbl['binning'][frame])
        # Write 2D images for the Science Frame
        # Need raw file header information
        # TODO: Why is the raw file header needed?  Can the data be
        # gotten from fitstbl?  If not, is it worth adding the relevant
        # information to fitstbl?
        rawfile = self.fitstbl.frame_paths(frame)
        # TODO: Make sure self.det is correct!
        #master_key = self.fitstbl.master_key(frame, det=self.det)
        outfile = scidir + '/spec2d_{:s}.fits'.format(basename)
        save.save_2d_images(sci_dict, rawfile, self.spectrograph.primary_hdrext,
                            self.caliBrate.master_key_dict, self.caliBrate.master_dir, outfile,
                            update_det=self.par['rdx']['detnum'])
        return all_specobjs



    def msgs_reset(self):
        """
        Reset the msgs object

        Returns:

        """

        # Reset the global logger
        msgs.reset(log=self.logname, verbosity=self.verbosity)
        msgs.pypeit_file = self.pypeit_file

    def print_end_time(self):
        """
        Print the elapsed time

        Returns:

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

        Returns:

        """
        indx = self.fitstbl.find_frames('science')
        print(self.fitstbl[['target','ra','dec','exptime','dispname']][indx])

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


class MultiSlit(PypeIt):
    """
    Child of PypeIt for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, **kwargs):
        super(MultiSlit, self).__init__(spectrograph, **kwargs)

        # WARNING: Defining these here means it can't be used by the
        # parent class...

        # TODO: Will you need these in Echelle as well?  Do you think
        # effectively *any* PypeIt object would need them?  If so, we
        # should move them into the parent object.
        self.std_idx = None
        self.std_dict = None
        self.std_basename = None
        self.stdI = None


    def find_objects(self, image, std=False, ir_redux=False, std_trace=None, snr_trim=False, maskslits=None,
                          show_peaks=False, show_fits=False, show_trace=False, show=False):

        sobjs_obj_init, nobj_init, skymask_pos = \
            self.sciI.find_objects(image, std=std, std_trace=std_trace, maskslits=maskslits,
                                   show_peaks = show_peaks, show_fits = show_fits, show_trace = show_trace)

        if ir_redux:
            sobjs_obj_init_neg, nobj_init_neg, skymask_neg = \
                self.sciI.find_objects(-image, std=std, std_trace=std_trace, maskslits=maskslits,
                show_peaks = show_peaks, show_fits = show_fits, show_trace = show_trace)
            skymask = skymask_pos & skymask_neg
            sobjs_obj_init.append_neg(sobjs_obj_init_neg)
        else:
            skymask = skymask_pos

        if show:
            self.sciI.show('image', image=image*(self.mask == 0), chname='objfind',sobjs=sobjs_obj_init, slits=True)

        return sobjs_obj_init, len(sobjs_obj_init), skymask

    # TODO: I don't think this function is written yet...
    def flux_calibrate(self):
        """
        Doc it
        """
        # Flux?
        if self.par['fluxcalib'] is None or len(self.std_dict) == 0:
            msgs.info('Flux calibration is not performed.')
            return
        elif self.par['fluxcalib'] is None and len(self.std_dict) > 0:
            msgs.info('Flux calibration parameters not provided.  Standards not used.')
            return

        # Standard star (is this a calibration, e.g. goes above?)
        msgs.info("Taking one star per detector mosaic")
        msgs.info("Waited until very end to work on it")
        msgs.warn("You should probably consider using the pypeit_flux_spec script anyhow...")

        # Get the sensitivity function
        if self.par['fluxcalib']['sensfunc'] is None:
            # Take the first standard
            std_idx = list(self.std_dict.keys())[0]
            # Build the list of stdobjs
            #all_std_objs = []
            #for det in self.std_dict[std_idx].keys():
            #    all_std_objs += self.std_dict[std_idx][det]['specobjs']
            # Need the Header for RA/DEC
            std_header = {}
            for key in ['ra', 'dec', 'airmass', 'exptime']:
                std_header[key.upper()] = self.fitstbl[std_idx][key]
            # Go
            # TODO: This going to be wrong.  How des FluxSpec iterate
            # through detectors, and then how does it use "setup"
            setup = self.fitstbl.master_key(self.frame, det=self.det)
            FxSpec = fluxspec.FluxSpec(std_specobjs=std_spec_objs.specobjs,
                                       spectrograph=self.spectrograph, setup=self.setup,
                                       master_dir=self.caliBrate.master_dir,
                                       std_header=std_header,
                                       mode=self.par['calibrations']['masters'])
            sens_dict = FxSpec.get_sens_dict(self.fitstbl[std_idx])
        else:
            # User provided it
            FxSpec = fluxspec.FluxSpec(sens_file=self.par['fluxcalib']['sensfunc'],
                                       spectrograph=self.spectrograph,
                                       master_dir=self.caliBrate.master_dir,
                                       mode=self.par['calibrations']['masters'])
            sens_dict = FxSpec.sens_dict

        # Apply the flux calibration
        msgs.info("Fluxing with {:s}".format(sens_dict['std_name']))
        save_format = 'fits'
        for kk, sci_ID in enumerate(all_sci_ID):
            # Load from disk (we zero'd out the object to free memory)
            if save_format == 'fits':
                sci_spec1d_file = os.path.join(self.par['rdx']['scidir'],
                                               'spec1d_{:s}.fits'.format(basenames[kk]))

            # Load
            sci_specobjs, sci_header = load.load_specobj(sci_spec1d_file)
            # TODO: (KBW) This is dangerous.  We want FluxSpec to check
            # that its internals make sense and this bypasses any of
            # that checking.  You should be able to instantiate FluxSpec
            # from a specobj...
            FxSpec.sci_specobjs = sci_specobjs
            FxSpec.sci_header = sci_header

            # Flux
            FxSpec.flux_science()
            # Over-write
            FxSpec.write_science(sci_spec1d_file)

#    def _extract_std(self):
#        self._extract_one(std=True)

    # THESE ARENT USED YET BUT WE SHOULD CONSIDER IT
    @staticmethod
    def default_sci_find_obj_steps():
        return ['process', 'find_objects', 'global_skysub', 'find_objects']

    @staticmethod
    def default_std_find_obj_steps():
        return ['process', 'global_skysub', 'find_objects']



class Echelle(PypeIt):
    """
    Child of PypeIt for Multislit and Longslit reductions

    """
    def __init__(self, spectrograph, **kwargs):
        super(Echelle, self).__init__(spectrograph, **kwargs)


    def find_objects(self, image, std=False, ir_redux=False, std_trace=None, snr_trim=False, maskslits=None,
                          show_peaks=False, show_fits=False, show_trace=False, show=False):

        sobjs_obj_init, nobj_init, skymask_pos = \
            self.sciI.find_objects_ech(image, std=std, std_trace=std_trace, snr_trim=snr_trim,
                                   show_peaks = show_peaks, show_fits = show_fits, show_trace = show_trace)

        if ir_redux:
            sobjs_obj_init_neg, nobj_init_neg, skymask_neg = \
                self.sciI.find_objects_ech(-image, std=std, std_trace=std_trace, snr_trim=snr_trim,
                show_peaks = show_peaks, show_fits = show_fits, show_trace = show_trace)
            skymask = skymask_pos & skymask_neg
            sobjs_obj_init.append_neg(sobjs_obj_init_neg)
        else:
            skymask = skymask_pos

        if show:
            self.sciI.show('image', image=image*(self.mask == 0), chname='ech_objfind',
                      sobjs=sobjs_obj_init, slits=True)

        return sobjs_obj_init, len(sobjs_obj_init), skymask


    # THESE ARENT USED YET BUT WE SHOULD CONSIDER IT
    @staticmethod
    def default_sci_find_obj_steps():
        return ['process', 'find_objects', 'global_skysub', 'find_objects']

    @staticmethod
    def default_std_find_obj_steps():
        return ['process', 'global_skysub', 'find_objects']


def instantiate_me(spectrograph, pypeit_file, **kwargs):
    """
    Instantiate the PypeIt subclass appropriate for the provided
    spectrograph.

    The class must be subclassed from PypeIt.  See :class:`PypeIt` for
    the description of the valid keyword arguments.

    Args:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The instrument used to collect the data to be reduced.

    Returns:
        :class:`PypeIt`: One of the classes with :class:`PypeIt` as its
        base.
    """
    indx = [ c.__name__ == spectrograph.pypeline for c in PypeIt.__subclasses__() ]
    if not np.any(indx):
        msgs.error('Pipeline {0} is not defined!'.format(spectrograph.pypeline))
    return PypeIt.__subclasses__()[np.where(indx)[0][0]](pypeit_file, **kwargs)

