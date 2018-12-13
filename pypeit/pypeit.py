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
from pypeit.core import pypsetup
from pypeit.core import wave
from pypeit.core import save
from pypeit.core import load
from pypeit.spectrographs.util import load_spectrograph
from pypeit.scripts import run_pypeit

from pypeit import debugger

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

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, logname=None, show=False,
                 redux_path=None):

        # Setup
        self.pypeit_file = pypeit_file
        ps = pypeitsetup.PypeItSetup.from_pypeit_file(self.pypeit_file)
        ps.run(setup_only=False)
        '''
        ps.build_fitstbl(strict=True)
        ps.get_frame_types()

        # Determine the configurations and assign each frame to the
        # specified configuration
        cfgs = self.fitstbl.unique_configurations(ignore_frames=['bias', 'dark'])
        self.fitstbl.set_configurations(cfgs)

        # Assign frames to calibration groups
        self.fitstbl.set_calibration_groups(global_frames=['bias', 'dark'])
        '''
        # Only need the parameters, spectrograph, and metadata for the remainder
        self.par = ps.par
        # self.spectrograph = ps.spectrograph
        self.fitstbl = ps.fitstbl

        self.pypeitSetup = ps

        # Other Internals
        self.logname = logname
        self.overwrite = overwrite
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

#    def calibrate_one(self, frame, det=1):
#        """
#        Calibrate a science exposure / detector pair
#
#        Args:
#            frame (:obj:`int`):
#                0-indexed index of the row in :attr:`fitstbl` to calibrate
#            det (:obj:`int`, optional):
#                1-indexed detector on this frame to calibrate
#        """
#        self.caliBrate.reset(frame, det, par=self.par['calibrations'])
#        self.caliBrate.run_the_steps()

#    def _chk_for_std(self):
#        # Can only reduce these frames if the mask is the same
#        std_idx = self.fitstbl.find_frames('standard', sci_ID=self.sci_ID, index=True)
#        if len(std_idx) > 0:
#            std_idx = std_idx[0]
#            return std_idx
#        else:
#            msgs.info("No standard star associated with this science frame")
#            return -1

    def reduce_all(self, reuse_masters=False):
        """
        Reduce all of the science exposures
        Generate all needed calibration files

        Args:
            reuse_masters (:obj:`bool`, optional):
                Use the master frames if available (same as setting
                par['calibrations']['masters'] = 'reuse'.

        Returns:

        """
        # Validate the parameter set
        required = ['rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib']
        can_be_None = ['flexure', 'fluxcalib']
        self.par.validate_keys(required=required, can_be_None=can_be_None)

        # TODO: add arguments and return values here to make the control
        # flow understandable
        self.tstart = time.time()

        # Find the standard frames
        is_standard = self.fitstbl.find_frames('standard')

        # Find the science frames
        is_science = self.fitstbl.find_frames('science')

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))

        # Iterate through each calibration group
        for i in range(self.fitstbl.n_calib_groups):

            # TODO: Could put everything in this loop into a new
            # function: reduce_calibgroup(i)
            
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(i)

            # Find the indices of the standard frames in this
            # calibration group:
            grp_standards = frame_indx[is_standard & in_grp]

            # TODO -- Turn standards back on!
            '''
            # Reduce all the standard frames
            for frame in grp_standards:
                # This sets: frame, sciI, obstime, basename
                # reduce_exposure(filename, group, std=False,
                std_dict = self.reduce_exposure(frame, reuse_masters=reuse_masters)
                self.save_exposure(frame, std_dict, self.basename)
            '''

            # Find the indices of the science frames in this calibration
            # group:
            grp_science = frame_indx[is_science & in_grp]

            # TODO: Need to decide how to associate standards with
            # science frames in the case where there is more than one
            # standard associated with a given science frame.  Below, I
            # just use the first standard
            std_frame = None if len(grp_standards) == 0 else grp_standards[0]

            # Reduce all the science frames; keep the basenames of the
            # science frames for use in flux calibration
            science_basename = [None]*len(grp_science)
            for j,frame in enumerate(grp_science):
                # This sets: frame, sciI, obstime, basename
                sci_dict = self.reduce_exposure(frame, std_frame=std_frame,
                                                reuse_masters=reuse_masters)
                science_basename[j] = self.basename
                self.save_exposure(frame, sci_dict, science_basename)

            # Apply the flux calibration for this calibration group
            # TODO: I don't think this function is written yet...
            #self.flux_calibrate(reuse_masters=False)

            msgs.info('Finished calibration group {0}'.format(i))

        # Finish
        self.print_end_time()


#    def reduce_all_old(self, reuse_masters=False):
#        """
#        Reduce all of the science exposures
#        Generate all needed calibration files
#
#        Args:
#            reuse_masters (:obj:`bool`, optional):
#                Use the master frames if available (same as setting
#                par['calibrations']['masters'] = 'reuse'.
#
#        Returns:
#
#        """
#
#        self.tstart = time.time()
#        self.std_dict = {}
#        # Science IDs are in a binary system: 1,2,4,8, etc.
#        all_sci_ID = self.fitstbl['sci_ID'][self.fitstbl.find_frames('science')]
#        numsci = len(all_sci_ID)
#        basenames = [None]*numsci  # For fluxing at the very end
#
#        # Check par
#        required = ['rdx', 'calibrations', 'scienceframe', 'scienceimage', 'flexure', 'fluxcalib']
#        can_be_None = ['flexure', 'fluxcalib']
#        self.par.validate_keys(required=required, can_be_None=can_be_None)
#
#        # Save
#        for kk,sci_ID in enumerate(all_sci_ID):
#            sci_dict = self.reduce_exposure(sci_ID, reuse_masters=reuse_masters)
#            scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
#            self.save_exposure(scidx, sci_dict, self.basename)
#            basenames[kk] = self.basename
#
#        # Standard stars
#        for std_idx in self.std_dict.keys():
#            # Basename
#            ikey = list(self.std_dict[std_idx].keys())[0]  # Any will do, so take the first
#            std_spec_objs = self.save_exposure(std_idx, self.std_dict[std_idx],
#                                               self.std_dict[std_idx][ikey]['basename'])
#
#        # Flux?
#        if self.par['fluxcalib'] is None or len(self.std_dict) == 0:
#            msgs.info('Flux calibration is not performed.')
#        elif self.par['fluxcalib'] is None and len(self.std_dict) > 0:
#            msgs.info('Flux calibration parameters not provided.  Standards not used.')
#        else:
#            # Standard star (is this a calibration, e.g. goes above?)
#            msgs.info("Taking one star per detector mosaic")
#            msgs.info("Waited until very end to work on it")
#            msgs.warn("You should probably consider using the pypeit_flux_spec script anyhow...")
#
#            # Get the sensitivity function
#            if self.par['fluxcalib']['sensfunc'] is None:
#                # Take the first standard
#                std_idx = list(self.std_dict.keys())[0]
#                # Build the list of stdobjs
#                #all_std_objs = []
#                #for det in self.std_dict[std_idx].keys():
#                #    all_std_objs += self.std_dict[std_idx][det]['specobjs']
#                # Need the Header for RA/DEC
#                std_header = {}
#                for key in ['ra', 'dec', 'airmass', 'exptime']:
#                    std_header[key.upper()] = self.fitstbl[std_idx][key]
#                # Go
#                FxSpec = fluxspec.FluxSpec(std_specobjs=std_spec_objs.specobjs, spectrograph=self.spectrograph,
#                                           setup=self.setup, master_dir=self.caliBrate.master_dir, std_header=std_header, mode=self.par['calibrations']['masters'])
#                sens_dict = FxSpec.get_sens_dict(self.fitstbl[std_idx])
#            else:
#                # User provided it
#                FxSpec = fluxspec.FluxSpec(sens_file=self.par['fluxcalib']['sensfunc'],
#                                           spectrograph=self.spectrograph, master_dir=self.caliBrate.master_dir,
#                                           mode=self.par['calibrations']['masters'])
#                sens_dict = FxSpec.sens_dict
#
#            # Apply the flux calibration
#            msgs.info("Fluxing with {:s}".format(sens_dict['std_name']))
#            save_format = 'fits'
#            for kk, sci_ID in enumerate(all_sci_ID):
#                # Load from disk (we zero'd out the object to free memory)
#                if save_format == 'fits':
#                    sci_spec1d_file = os.path.join(self.par['rdx']['scidir'],
#                                                   'spec1d_{:s}.fits'.format(basenames[kk]))
#
#                # Load
#                sci_specobjs, sci_header = load.load_specobj(sci_spec1d_file)
#                # TODO: (KBW) I'm wary of this kind of approach.  We want
#                # FluxSpec to check that its internals make sense and this
#                # bypasses any of that checking.
#                FxSpec.sci_specobjs = sci_specobjs
#                FxSpec.sci_header = sci_header
#
#                # Flux
#                FxSpec.flux_science()
#                # Over-write
#                FxSpec.write_science(sci_spec1d_file)
#
#        # Finish
#        self.print_end_time()

    def select_detectors(self):
        """
        Return the 1-indexed list of detectors to reduce.
        """
        if self.par['rdx']['detnum'] is None:
            return np.arange(self.spectrograph.ndet)+1
        return [self.par['rdx']['detnum']] if isinstance(self.par['rdx']['detnum'], int) \
                    else self.par['rdx']['detnum']

    # JFH ToDO this operates on a file, and takes bgframe as an input, stdframe, is_std = False
    def reduce_exposure(self, frame, std_frame=None, reuse_masters=False):
        """
        Reduce a single exposure

        Args:
            frame (:obj:`int`):
                0-indexed row in :attr:`fitstbl` with the frame to
                reduce
            std_frame (:obj:`int`, :obj:`str`, optional):
                0-indexed row in :attr:`fitstbl` with a standard frame
                associated with the frame to reduce, or the name of a
                file with a previously PypeIt-reduced standard spectrum.
            reuse_masters (:obj:`bool`, optional):
                Re-use MasterFrame files when available

        Returns:
            dict: The dictionary containing the primary outputs of
            extraction
        """
        # Prepare to load up standard?
        if std_frame is not None:
            std_outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],
                                       'spec1d_{:s}.fits'.format(
                                            self.fitstbl.construct_basename(std_frame))) \
                                if isinstance(std_frame, int) else std_frame
            if std_outfile is not None and not os.path.isfile(std_outfile):
                msgs.error('Could not open standard file: {0}'.format(std_outfile))

        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            ginga.clear_all()

        # Save the frame
        self.frame = frame

        # Insist on re-using MasterFrames where applicable
        if reuse_masters:
            self.par['calibrations']['masters'] = 'reuse'

        sci_dict = OrderedDict()  # This needs to be ordered
        sci_dict['meta'] = {}
        sci_dict['meta']['vel_corr'] = 0.
        #
        msgs.info("Reducing file {0:s}, target {1:s}".format(self.fitstbl['filename'][self.frame],
                                                             self.fitstbl['target'][self.frame]))

        # Check if the frame is a standard
        is_standard = self.frame in self.fitstbl.find_frames('standard', index=True)

        # Find the detectors to reduce
        detectors = self.select_detectors()
        if len(detectors) != self.spectrograph.ndet:
            msgs.warn('Not reducing detectors: {0}'.format(' '.join([ str(d) for d in 
                                set(np.arange(self.spectrograph.ndet))-set(detectors)])))

        # Loop on Detectors
        for self.det in detectors:
            msgs.info("Working on detector {0}".format(self.det))
            sci_dict[self.det] = {}
            setup = self.fitstbl.master_key(self.frame, det=self.det)

            # Calibrate
            self.caliBrate.set_config(self.frame, self.det, self.par['calibrations'])
            self.caliBrate.run_the_steps()

            # Initialize the time and output file root
            #   - This sets frame, det, sciI, obstime, basename
            self.init_one_science(self.frame, det=self.det)

            # Extract
            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            sci_dict[det]['sciimg'], sci_dict[det]['sciivar'], sci_dict[det]['skymodel'], \
                sci_dict[det]['objmodel'], sci_dict[det]['ivarmodel'], sci_dict[det]['outmask'], \
                sci_dict[det]['specobjs'], vel_corr \
                    = self._extract_one(scifiles, bgframes, std=is_standard,
                                        std_outfile=std_outfile)
            if vel_corr is not None:
                sci_dict['meta']['vel_corr'] = vel_corr

            # JFH TODO write out the background frame

        # Return
        return sci_dict


#    def reduce_exposure_old(self, frame, reuse_masters=False):
#        """
#        Reduce a single science exposure
#
#        Args:
#            sci_ID: int
#              binary flag indicating the science frame
#            reuse_masters: bool, optional
#              Reuse MasterFrame files (where available)
#
#
#        Returns:
#            sci_dict: dict
#              dict containing the primary outputs of extraction
#
#        """
#        self.sci_ID = sci_ID
#
#        # Insist on re-using MasterFrames where applicable
#        if reuse_masters:
#            self.par['calibrations']['masters'] = 'reuse'
#
#        sci_dict = OrderedDict()  # This needs to be ordered
#        sci_dict['meta'] = {}
#        sci_dict['meta']['vel_corr'] = 0.
#        #
#        scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
#        msgs.info("Reducing file {0:s}, target {1:s}".format(self.fitstbl['filename'][scidx],
#                                                             self.fitstbl['target'][scidx]))
#
#        # Loop on Detectors
#        for kk in range(self.spectrograph.ndet):
#            det = kk + 1  # Detectors indexed from 1
#            self.det = det
#            if self.par['rdx']['detnum'] is not None:
#                detnum = [self.par['rdx']['detnum']] if isinstance(self.par['rdx']['detnum'],int) else self.par['rdx']['detnum']
#                if det not in map(int, detnum):
#                    msgs.warn("Skipping detector {:d}".format(det))
#                    continue
#                else:
#                    msgs.warn("Restricting the reduction to detector {:d}".format(det))
#            # Setup
#            msgs.info("Working on detector {0}".format(det))
#            sci_dict[det] = {}
#
#            # Calibrate
#            self.calibrate_one(frame, det)
#
#            # Init ScienceImage class
#            self.init_one_science(sci_ID, det)
#            # Extract
#            sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr = self._extract_one()
#
#            # Save for outputing (after all detectors are done)
#            sci_dict[det]['sciimg'] = sciimg
#            sci_dict[det]['sciivar'] = sciivar
#            sci_dict[det]['skymodel'] = skymodel
#            sci_dict[det]['objmodel'] = objmodel
#            sci_dict[det]['ivarmodel'] = ivarmodel
#            sci_dict[det]['outmask'] = outmask
#            sci_dict[det]['specobjs'] = sobjs   #utils.unravel_specobjs([specobjs])
#            if vel_corr is not None:
#                sci_dict['meta']['vel_corr'] = vel_corr
#
#            # Standard star
#            # TODO -- Make this more modular
#            self.std_idx = self._chk_for_std()
#            if self.std_idx is not -1:
#                self._extract_std()
#
#        # Return
#        return sci_dict

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
        save_format = 'fits'
        if save_format == 'fits':
            outfile = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'],
                                   'spec1d_{:s}.fits'.format(basename))
            helio_dict = dict(refframe='pixel'
            if self.caliBrate.par['wavelengths']['reference'] == 'pixel'
            else self.caliBrate.par['wavelengths']['frame'],
                              vel_correction=vel_corr)
            # Did the user re-run a single detector?
            save.save_1d_spectra_fits(all_specobjs, self.fitstbl[frame], outfile,
                                      helio_dict=helio_dict, telescope=self.spectrograph.telescope,
                                      update_det=self.par['rdx']['detnum'])
        #        elif save_format == 'hdf5':
        #            debugger.set_trace()  # NEEDS REFACTORING
        #            arsave.save_1d_spectra_hdf5(None)
        else:
            msgs.error(save_format + ' is not a recognized output format!')
        # 1D only?
        if only_1d:
            return
        # Obj info
        scidir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['scidir'])
        save.save_obj_info(all_specobjs, self.fitstbl, self.spectrograph, basename, scidir)
        # Write 2D images for the Science Frame
        # Need raw file header information
        # TODO: Why is the raw file header needed?  Can the data be
        # gotten from fitstbl?  If not, is it worth adding the relevant
        # information to fitstbl?
        rawfile = self.fitstbl.frame_paths(frame)
        # TODO: Make sure self.det is correct!
        setup = self.fitstbl.master_key(frame, det=self.det)
        save.save_2d_images(sci_dict, rawfile, self.spectrograph.primary_hdrext,
                            setup, self.caliBrate.master_dir, scidir, basename,
                            update_det=self.par['rdx']['detnum'])
        return all_specobjs

    def _extract_one(self, std=False):
        """
        Dummy method for object extraction. Overloaded by class specific extraction.

        Returns:

        """
        assert False

#    def _extract_std(self):
#        """
#        Dummy method for std extraction
#
#        Returns:
#
#        """
#        assert False

# This is no longer required

#    def _init_calibrations(self):
#        """
#        Instantiate the Calibrations class
#        Returns:
#
#        """
#        # TODO -- Need to make save_masters and write_qa optional
#        # Init calib dict
#        self.caliBrate \
#                = calibrations.MultiSlitCalibrations(self.fitstbl, spectrograph=self.spectrograph,
#                                                     par=self.par['calibrations'],
#                                                     redux_path=self.par['rdx']['redux_path'],
#                                                     save_masters=True, write_qa=True,
#                                                     show=self.show)

    def init_one_science(self, frame, det=1):
        """
        Instantiate ScienceImage class and run the first step with it

        Args:
            frame (:obj:`int`):
                0-indexed index of the row in :attr:`fitstbl` to calibrate
            det (:obj:`int`, optional):
                1-indexed detector on this frame to calibrate
        """
        self.frame = frame
        self.det = det

        sci_image_files = self.fitstbl.frame_paths(self.frame)
        self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_image_files, det=self.det,
                                               binning=self.fitstbl['binning'][self.frame],
                                               objtype=self.fitstbl['frametype'][self.frame],
                                               scidx=self.frame,
                                               setup=self.fitstbl.master_key(self.frame, det=det),
                                               par=self.par['scienceimage'],
                                               frame_par=self.par['scienceframe'])
        # For QA on crash
        msgs.sciexp = self.sciI
        self.obstime = self.fitstbl.construct_obstime(self.frame)
        self.basename = self.fitstbl.construct_basename(self.frame, obstime=self.obstime)

#    def init_one_science_old(self, sci_ID, det):
#        """
#        Instantiate ScienceImage class and run the first step with it
#
#        Args:
#            sci_ID: int
#              binary flag indicating the science frame
#            det: int
#              detector index
#
#        Returns:
#            self.obstime : Time
#            self.basename : str
#        """
#        self.sci_ID = sci_ID
#        self.det = det
#
#        sci_image_files = self.fitstbl.find_frame_files('science', sci_ID=sci_ID)
#        scidx = self.fitstbl.find_frames('science', sci_ID=sci_ID, index=True)[0]
#        try:
#            binning = self.fitstbl['binning'][scidx]
#        except:
#            binning = (1,1)
#
#        self.sciI = scienceimage.ScienceImage(self.spectrograph, sci_image_files, det=det,
#                                              binning = binning,
#                                              objtype='science', scidx=scidx, setup=self.setup,
#                                              par=self.par['scienceimage'],
#                                              frame_par=self.par['scienceframe'])
#        msgs.sciexp = self.sciI  # For QA on crash
#
#        # Names and time
#        self.obstime, self.basename = self.init_time_names(self.fitstbl, scidx)
#        # Return
#        return self.obstime, self.basename  # For fluxing

    # Move out of the scienceimage class to here. It is more appropriate here where fitstable exists
    # Moved to fitstbl
#    def _init_time_names(self):
#        """
#        Setup the basename (for file output mainly)
#        and time objects (for heliocentric)
#
#        Parameters
#        ----------
#        camera : str
#          Taken from settings['mosaic']['camera']
#        timeunit : str
#          mjd
#
#        Returns
#        -------
#        self.time : Time
#        self.basename : str
#
#        """
#
#        timeunit = self.spectrograph.timeunit
#        camera = self.spectrograph.camera
#
#        self.fitstbl = fitstbl
#
#        # TODO: Given that we just read the file header to get the
#        # datasec_img in the init function above, I don't see why I need
#        # to access the fits table for exptime and binning. This
#        # information is also in the headers. By simply pulling the
#        # stuff from the header, we would remove the fitstbl entirely.
#        # Another option would be to put the datasec_img stuff in the
#        # fitstbl for each detector
#
#        self.exptime = self.fitstbl['exptime'][scidx]
#        try:
#            self.binning = self.fitstbl['binning'][scidx]
#        except:
#            self.binning = (1,1)
#
#        # This should have been set when we construct the fitstbl.
#        #
#        # JFH These time routines need to exit cleanly with warnings rather than crashing the code
#        # until we get the fitstbl working in as stable way.
#        try:
#            tval = Time(fitstbl['time'][scidx], format='mjd')#'%Y-%m-%dT%H:%M:%S.%f')
#        except:
#            msgs.warn('There is no time in the fitstbl.' + msgs.newline() +
#            'The time and heliocentric corrections will be off!!' + msgs.newline() +
#            'This is a bad idea. Continuing with a dummy time value')
#            tval = '2010-01-01'
#        # Time
#        tiso = Time(tval, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
#        dtime = datetime.datetime.strptime(tiso.value, '%Y-%m-%dT%H:%M:%S.%f')
#        self.time = tval
#        # Basename
#        self.inst_name = camera
#        self.target_name = self.fitstbl['target'][scidx].replace(" ", "")
#        self.basename = self.target_name+'_'+self.inst_name+'_'+ \
#                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
#                         tiso.value.split("T")[1].replace(':','')
#        # Return
#        return self.time, self.basename


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

    def _extract_one(self, std=False):
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
        # Standard star specific
        if std:
            # Dict
            msgs.info("Processing standard star")
            # TODO: Where is self.std_idx defined
            if self.std_idx in self.std_dict.keys():
                if self.det in self.std_dict[self.std_idx].keys():
                    return
            else:
                self.std_dict[self.std_idx] = {}

            # Files
            # TODO: (KBW) I don't understand why you're selecting all
            # standards here, but this is the new way to do it.
            is_standard = self.fitstbl.find_frames('standard')
            std_image_files = self.fitstbl.frame_paths(is_standard)
            if self.par['calibrations']['standardframe'] is None:
                msgs.warn('No standard frame parameters provided.  Using default parameters.')

            # Instantiate for the Standard
            # TODO: Uses the same trace and extraction parameter sets used for the science
            # frames.  Should these be different for the standards?
            setup = self.fitstbl.master_key(self.frame, det=self.det)
            self.stdI = scienceimage.ScienceImage(self.spectrograph, file_list=std_image_files,
                                          frame_par=self.par['calibrations']['standardframe'],
                                          det=self.det,
                                          binning=self.fitstbl['binning'][self.std_idx],
                                          setup=setup, scidx=self.std_idx, objtype='standard',
                                          par=self.par['scienceimage'])
            # Names and time
            self.std_basename = self.fitstbl.construct_basename(self.std_idx)
            sciI = self.stdI
        else:
            sciI = self.sciI

        # Process images (includes inverse variance image, rn2 image,
        # and CR mask)
        sciimg, sciivar, rn2img, crmask \
                = sciI.proc(self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
                               self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                               apply_gain=True, trim=self.caliBrate.par['trim'], show=self.show)

        # Object finding, first pass on frame without sky subtraction
        maskslits = self.caliBrate.maskslits.copy()
        if not std:
            sobjs_obj0, nobj0 = sciI.find_objects(self.caliBrate.tslits_dict, skysub=False,
                                                   maskslits=maskslits)

        # Global sky subtraction, first pass. Uses skymask from object
        # finding
        global_sky0 = sciI.global_skysub(self.caliBrate.tslits_dict,
                                         self.caliBrate.tilts_dict['tilts'],
                                         use_skymask=True,maskslits=maskslits, show=self.show)

        # Object finding, second pass on frame *with* sky subtraction.
        # Show here if requested
        sobjs_obj, nobj = sciI.find_objects(self.caliBrate.tslits_dict, skysub=True,
                                            maskslits=maskslits, show_peaks=self.show)

        if std:
            if nobj == 0:
                msgs.warn('No objects to extract for standard frame' + msgs.newline()
                          + self.fitstbl['filename'][self.sciI.scidx])
                return
            # Extract
            skymodel, objmodel, ivarmodel, outmask, sobjs \
                    = self.stdI.local_skysub_extract(sobjs_obj, self.caliBrate.mswave,
                                                     maskslits=maskslits, show_profile=self.show,
                                                     show=self.show)

            # Save for fluxing and output later
            self.std_dict[self.std_idx][self.det] = {}
            self.std_dict[self.std_idx][self.det]['basename'] = self.std_basename
            self.std_dict[self.std_idx][self.det]['specobjs'] = sobjs
            # Done
            return

        # If there are objects, do 2nd round of global_skysub,
        # local_skysub_extract, flexure, geo_motion
        vel_corr = None
        if nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            global_sky = sciI.global_skysub(self.caliBrate.tslits_dict,
                                            self.caliBrate.tilts_dict['tilts'], use_skymask=True,
                                            maskslits=maskslits, show=self.show)

            skymodel, objmodel, ivarmodel, outmask, sobjs \
                    = sciI.local_skysub_extract(sobjs_obj, self.caliBrate.mswave,
                                                maskslits=maskslits, show_profile=self.show,
                                                show=self.show)

            # Flexure correction?
            if self.par['flexure']['method'] != 'skip':
                sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                             sky_spectrum, sky_file=sky_file,
                                             mxshft=self.par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,
                                out_dir=self.par['rdx']['redux_path'])

            # Helio
            # Correct Earth's motion
            # vel_corr = -999999.9
            if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                    and (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
                if sobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                                                    self.caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr \
                            = wave.geomotion_correct(sobjs, maskslits, self.fitstbl,
                                                     self.sciI.scidx, self.obstime,
                                                     self.spectrograph.telescope['longitude'],
                                                     self.spectrograph.telescope['latitude'],
                                                     self.spectrograph.telescope['elevation'],
                                                     self.caliBrate.par['wavelengths']['frame'])
                else:
                    msgs.info('There are no objects on detector {0} to perform a '.format(self.det)
                              + '{1} correction'.format(self.caliBrate.par['wavelengths']['frame']))
            else:
                msgs.info('A wavelength reference-frame correction will not be performed.')

        else:
            msgs.warn('No objects to extract for science frame' + msgs.newline()
                      + self.fitstbl['filename'][self.sciI.scidx])
            # set to first pass global sky
            skymodel = global_sky0
            objmodel = np.zeros_like(sciimg)
            # Set to sciivar. Could create a model but what is the point?
            ivarmodel = np.copy(sciivar)
            # Set to inmask in case on objects were found
            outmask = sciI.mask
            # empty specobjs object from object finding
            sobjs = sobjs_obj

        return sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr


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

    # JFH Take background files as an argument
    def _extract_one(self,std=False):
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

        # TODO Move the init_one_science functinon lines here. No three line functions!

        sciI = self.sciI

        # Process images (includes inverse variance image, rn2 image, and CR mask)
        # JFH ToDO this should take a file list!
        sciimg, sciivar, rn2img, crmask \
                = sciI.proc(self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
                            self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
                            apply_gain=True, trim=self.caliBrate.par['trim'], show=self.show)

        # if bgfiles is not None:
        # sciimg, sciivar, rn2img, crmask = sciI.proc(bgfiles, self.caliBrate.msbias, self.caliBrate.mspixflatnrm,
        #                             self.caliBrate.msbpm, illum_flat=self.caliBrate.msillumflat,
        #                             apply_gain=True, trim=self.caliBrate.par['trim'], show=self.show)


        # Object finding, first pass on frame without sky subtraction
        maskslits = self.caliBrate.maskslits.copy()

        # Do one iteration of object finding, and sky subtract to get initial sky model
        initial_sky = sciI.get_init_sky(self.caliBrate.tslits_dict, self.caliBrate.tilts_dict['tilts'], show = self.show)

        sobjs_ech, nobj = sciI.get_ech_objects(self.caliBrate.tslits_dict, show=self.show)


        # If there are objects, do 2nd round of global_skysub,
        # local_skysub_extract, flexure, geo_motion
        vel_corr = None
        if nobj > 0:
            skymodel, objmodel, ivarmodel, outmask, sobjs \
                    = sciI.local_skysub_extract(sobjs_ech, self.caliBrate.mswave, maskslits=maskslits, std = std,
                                                show_profile=self.show, show=self.show)

            # Flexure correction?
            if self.par['flexure']['method'] != 'skip':
                sky_file, sky_spectrum = self.spectrograph.archive_sky_spectrum()
                flex_list = wave.flexure_obj(sobjs, maskslits, self.par['flexure']['method'],
                                             sky_spectrum, sky_file=sky_file,
                                             mxshft=self.par['flexure']['maxshift'])
                # QA
                wave.flexure_qa(sobjs, maskslits, self.basename, self.det, flex_list,
                                out_dir=self.par['rdx']['redux_path'])

            # Helio
            # Correct Earth's motion
            # vel_corr = -999999.9
            if (self.caliBrate.par['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                    and (self.caliBrate.par['wavelengths']['reference'] != 'pixel'):
                if sobjs is not None:
                    msgs.info("Performing a {0} correction".format(
                                                    self.caliBrate.par['wavelengths']['frame']))

                    vel, vel_corr \
                            = wave.geomotion_correct(sobjs, maskslits, self.fitstbl, self.sciI.scidx,
                                                     self.obstime,
                                                     self.spectrograph.telescope['longitude'],
                                                     self.spectrograph.telescope['latitude'],
                                                     self.spectrograph.telescope['elevation'],
                                                     self.caliBrate.par['wavelengths']['frame'])
                else:
                    msgs.info('There are no objects on detector {0} to perform a '.format(self.det)
                              + '{1} correction'.format(self.caliBrate.par['wavelengths']['frame']))
            else:
                msgs.info('A wavelength reference-frame correction will not be performed.')

        else:
            msgs.warn('No objects to extract for science frame' + msgs.newline()
                      + self.fitstbl['filename'][self.sciI.scidx])
            # set to first pass global sky
            skymodel = global_sky0
            objmodel = np.zeros_like(sciimg)
            # Set to sciivar. Could create a model but what is the point?
            ivarmodel = np.copy(sciivar)
            # Set to inmask in case on objects were found
            outmask = sciI.mask
            # empty specobjs object from object finding
            sobjs = sobjs_obj

        return sciimg, sciivar, skymodel, objmodel, ivarmodel, outmask, sobjs, vel_corr

#    def _extract_std(self):
#        self._extract_one(std=True)

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

