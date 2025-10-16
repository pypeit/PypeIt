"""
Main driver class for PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import copy
import datetime
import os
from pathlib import Path
import time

from astropy.io import fits
from astropy.table import Table
import numpy as np

from pypeit import io
from pypeit import inputfiles
from pypeit.calibframe import CalibFrame
from pypeit.core import parse, wave, qa
from pypeit import msgs
from pypeit import calibrations
from pypeit.images import buildimage
from pypeit.images import pypeitimage
from pypeit.display import display
from pypeit import find_objects
from pypeit import extraction
from pypeit import spec2dobj
from pypeit import specobjs
#from pypeit.spectrographs.util import load_spectrograph
from pypeit import slittrace
from pypeit import utils
from pypeit.history import History
#from pypeit.par import PypeItPar
#from pypeit.par.pypeitpar import ql_is_on
from pypeit.metadata import PypeItMetaData
from pypeit.manual_extract import ManualExtractionObj
from pypeit.core import skysub
from pypeit.core import procimg
from pypeit.core import flat

from linetools import utils as ltu

#TODO FIX THIS
try:
    from jwst import datamodels
except ModuleNotFoundError:
    raise ModuleNotFoundError('Unable to import jwst.  Install pypeit with the jwst '
                              'option to reduce jwst data.')
from IPython import embed


class PypeIt:
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
        reuse_calibs (:obj:`bool`, optional):
            Reuse any pre-existing calibration files
        logname (:obj:`str`, optional):
            The name of an ascii log file with the details of the
            reduction.
        show: (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via ``ginga --modules=RC,SlitWavelength &``
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

    @classmethod
    def get_instance(cls, pypeit_file, verbosity=2, overwrite=True, reuse_calibs=False, logname=None,
                 show=False, redux_path=None, calib_only=False):
        """
        Instantiate and return the :class:`PypeIt` subclass appropriate for
        the provided spectrograph.

        For argument descriptions, see :class:`PypeIt`.

        Returns:
            :class:`PypeIt`: One of the subclasses with :class:`PypeIt` as its base."""

        spectrograph = inputfiles.PypeItFile.from_file(pypeit_file).get_spectrograph()

        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (spectrograph.pypeline + 'PypeIt'))(pypeit_file,
                                                                        verbosity=verbosity,
                                                                        overwrite=overwrite,
                                                                        reuse_calibs=reuse_calibs,
                                                                        logname=logname,
                                                                        show=show,
                                                                        redux_path=redux_path,
                                                                        calib_only=calib_only)

    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_calibs=False, logname=None,
                 show=False, redux_path=None, calib_only=False):

        # Set up logging
        self.logname = logname
        self.verbosity = verbosity
        self.pypeit_file = pypeit_file
        
        self.msgs_reset()
        
        # Load up PypeIt file
        self.pypeItFile = inputfiles.PypeItFile.from_file(pypeit_file)
        self.calib_only = calib_only

        # Build the spectrograph and the parameters
        self.spectrograph, self.par, config_specific_file = self.pypeItFile.get_pypeitpar()
        msgs.info(f'Loaded spectrograph {self.spectrograph.name}')
        msgs.info('Setting configuration-specific parameters using '
                  f'{os.path.split(config_specific_file)[1]}.')

        # Check the output paths are ready
        if redux_path is not None:
            self.par['rdx']['redux_path'] = redux_path

        # Write the full parameter set here
        # --------------------------------------------------------------
        par_file = pypeit_file.replace(
            '.pypeit', f"_UTC_{datetime.datetime.now(datetime.UTC).date()}.par")
        self.par.to_config(par_file, include_descr=False)

        # --------------------------------------------------------------
        # Build the meta data
        #   - Re-initilize based on the file data
        msgs.info('Compiling metadata')
        self.fitstbl = PypeItMetaData(self.spectrograph, self.par, 
                                      files=self.pypeItFile.filenames,
                                      usrdata=self.pypeItFile.data, 
                                      strict=True)
        #   - Interpret automated or user-provided data from the PypeIt
        #   file
        self.fitstbl.finalize_usr_build(
            self.pypeItFile.frametypes, 
            self.pypeItFile.setup_name)

        # Other Internals
        self.overwrite = overwrite

        # Currently the runtime argument determines the behavior for
        # reusing calibrations
        self.reuse_calibs = reuse_calibs
        self.show = show

        # Set paths
        self.calibrations_path = Path(self.par['rdx']['redux_path']) / self.par['calibrations']['calib_dir']

        # Check for calibrations
        if not self.calib_only:
            calibrations.check_for_calibs(self.par, self.fitstbl,
                                          raise_error=self.par['calibrations']['raise_chk_error'])

        # --------------------------------------------------------------
        #   - Write .calib file (For QA naming amongst other things)
        calib_file = pypeit_file.replace('.pypeit', '.calib')
        calibrations.Calibrations.association_summary(calib_file, self.fitstbl, self.spectrograph,
                                                      self.calibrations_path, overwrite=True)

        # Report paths
        msgs.info('Setting reduction path to {0}'.format(self.par['rdx']['redux_path']))
        msgs.info('Calibration frames saved to: {0}'.format(self.calibrations_path))
        msgs.info('Science data output to: {0}'.format(self.science_path))
        msgs.info('Quality assessment plots output to: {0}'.format(self.qa_path))

        # Init
        self.det = None
        self.tstart = None
        self.basename = None
        self.obstime = None

    @property
    def science_path(self) -> Path:
        """Return the path to the science directory."""
        return Path(self.par['rdx']['redux_path']) / self.par['rdx']['scidir']

    @property
    def qa_path(self) -> str:
        """Return the path to the top-level QA directory."""
        return os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])

    def build_qa(self):
        """
        Generate QA wrappers
        """
        msgs.qa_path = self.qa_path
        qa.gen_qa_dir(self.qa_path)
        qa.gen_mf_html(self.pypeit_file, self.qa_path)
        qa.gen_exp_html()

    # TODO: This should go in a more relevant place
    def spec_output_file(self, frame:int, twod:bool=False) -> Path:
        """
        Return the path to the spectral output data file.
        
        Args:
            frame (:obj:`int`):
                Frame index from :attr:`fitstbl`.
            twod (:obj:`bool`), optional:
                Name for the 2D output file; 1D file otherwise.
        
        Returns:
            `Path`_: The path for the output file
        """
        return self.get_spec_file_name(self.science_path, self.fitstbl.construct_basename(frame),
                                       twod=twod)

    @staticmethod
    def get_spec_file_name(science_path:Path, basename:str, twod:bool=False) -> Path:
        """
        Get the spectrum filename

        Args:
            science_path (`Path`_):
                Path to the science files
            basename (:obj:`str`):
                Base name for this frame
            twod (:obj:`bool`), optional:
                Is this a 2D science frame?

        Returns:
            `Path`_: The spectrum filename
        """
        return science_path / f'spec{"2" if twod else "1"}d_{basename}.fits'

    def outfile_exists(self, frame:int) -> bool:
        """
        Check whether the 2D outfile of a given frame already exists

        Args:
            frame (:obj:`int`): Frame index from fitstbl

        Returns:
            :obj:`bool`: True if the 2d file exists, False if it does not exist
        """
        return self.spec_output_file(frame, twod=True).is_file()

    def get_std_outfile(self, standard_frames):
        """
        Return the spec1d file name for a reduced standard to use as a tracing
        crutch.

        The file is either constructed using the provided standard frame indices
        or it is directly pulled from the
        :class:`~pypeit.par.pypeitpar.FindObjPar` parameters in :attr:`par`.
        The latter takes precedence.  If more than one row is provided by
        ``standard_frames``, the first index is used.

        Args:
            standard_frames (array-like):
                Set of rows in :attr:`fitstbl` with standards.

        Returns:
            :obj:`str`: Full path to the standard spec1d output file to use.
        """
        # NOTE: I'm not sure if this is the best place to put this, but it does
        # isolate where the name of the standard-star spec1d file is defined.
        std_outfile = self.par['reduce']['findobj']['std_spec1d']
        if std_outfile is not None:
            if not self.par['reduce']['findobj']['use_std_trace']:
                msgs.error('If you provide a standard star spectrum for tracing, you must set use_std_trace=True')
            elif not Path(std_outfile).absolute().exists():
                msgs.error(f'Provided standard spec1d file does not exist: {std_outfile}')
            return std_outfile

        # TODO: Need to decide how to associate standards with
        # science frames in the case where there is more than one
        # standard associated with a given science frame.  Below, I
        # just use the first standard

        std_frame = None if (len(standard_frames) == 0 or not self.par['reduce']['findobj']['use_std_trace']) \
            else standard_frames[0]
        # Prepare to load up standard?
        if std_frame is not None:
            std_outfile = self.spec_output_file(std_frame) \
                            if isinstance(std_frame, (int,np.integer)) else None
        if std_outfile is not None and not std_outfile.is_file():
            msgs.error(f'Could not find standard file: {std_outfile}')
        return std_outfile

    def calib_all(self):
        """
        Process all calibration frames.

        Provides an avenue to reduce a dataset without (or omitting) any
        science/standard frames.
        """
        self.tstart = time.perf_counter()

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))
        for calib_ID in self.fitstbl.calib_groups:
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)
            if not any(in_grp):
                continue
            grp_frames = frame_indx[in_grp]

            # Find the detectors to reduce
            detectors = self.select_detectors(self.spectrograph, self.par['rdx']['detnum'],
                                              slitspatnum=self.par['rdx']['slitspatnum'])
            msgs.info(f'Detectors to work on: {detectors}')

            # Loop on Detectors
            for self.det in detectors:
                msgs.info(f'Working on detector {self.det}')

                self.caliBrate = self.calib_one(grp_frames, self.det)
                if not self.caliBrate.success:
                    msgs.warn(f'Calibrations for detector {self.det} were unsuccessful!  The step '
                              f'that failed was {self.caliBrate.failed_step}.  Continuing to next '
                              f'detector.')

        # Finish
        self.print_end_time()

    def reduce_all(self):
        """
        Main driver of the entire reduction

        Calibration and extraction via a series of calls to
        :func:`reduce_exposure`.

        """
        # Validate the parameter set
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

        # Find the standard frames
        is_standard = self.fitstbl.find_frames('standard')
        if np.any(is_standard):
            msgs.info(f'Found {np.sum(is_standard)} standard frames to reduce.')

        # Find the science frames
        is_science = self.fitstbl.find_frames('science')
        if np.any(is_science):
            msgs.info(f'Found {np.sum(is_science)} science frames to reduce.')

        # This will give an error to alert the user that no reduction will be
        # run if there are no science/standard frames and `run_pypeit` is run
        # without -c flag
        if not np.any(is_science) and not np.any(is_standard):
            msgs.error('No science/standard frames provided. Add them to your PypeIt file '
                       'if this is a standard run! Otherwise run calib_only reduction using -c flag')

        # Frame indices
        frame_indx = np.arange(len(self.fitstbl))

        # Standard Star(s) Loop
        # Iterate over each calibration group and reduce the standards
        for calib_ID in self.fitstbl.calib_groups:

            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)

            if not np.any(is_standard & in_grp):
                continue

            # Find the indices of the standard frames in this calibration group:
            grp_standards = frame_indx[is_standard & in_grp]

            msgs.info(f'Found {len(grp_standards)} standard frames in calibration group '
                      f'{calib_ID}.')

            # Reduce all the standard frames, loop on unique comb_id
            u_combid_std = np.unique(self.fitstbl['comb_id'][grp_standards])
            for j, comb_id in enumerate(u_combid_std):
                frames = np.where(self.fitstbl['comb_id'] == comb_id)[0]
                # Find all frames whose comb_id matches the current frames
                # bkg_id (same as for science frames).
                bg_frames = np.where((self.fitstbl['comb_id'] == self.fitstbl['bkg_id'][frames][0])
                                     & (self.fitstbl['comb_id'] >= 0))[0]
                if not self.outfile_exists(frames[0]) or self.overwrite:
                    # # Build history to document what contributed to the reduced
                    # # exposure
                    # history = History(self.fitstbl.frame_paths(frames[0]))
                    # history.add_reduce(calib_ID, self.fitstbl, frames, bg_frames)
                    self.reduce_exposure(frames, bg_frames=bg_frames)
                    # # TODO come up with sensible naming convention for save_exposure for combined files
                    # self.save_exposure(frames[0], std_spec2d, std_sobjs, self.basename, history)
                else:
                    msgs.info('Output file: {:s} already exists'.format(self.fitstbl.construct_basename(frames[0])) +
                              '. Set overwrite=True to recreate and overwrite.')

        # Science Frame(s) Loop
        # Iterate over each calibration group again and reduce the science frames
        for calib_ID in self.fitstbl.calib_groups:
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)

            if not np.any(is_science & in_grp):
                continue

            # Find the indices of the science frames in this calibration group:
            grp_science = frame_indx[is_science & in_grp]
            msgs.info(f'Found {len(grp_science)} science frames in calibration group {calib_ID}.')

            # Associate standards (previously reduced above) for this setup
            std_outfile = self.get_std_outfile(frame_indx[is_standard])
            # Reduce all the science frames; keep the basenames of the science
            # frames for use in flux calibration
            science_basename = [None]*len(grp_science)
            # Loop on unique comb_id
            u_combid = np.unique(self.fitstbl['comb_id'][grp_science])
        
            for j, comb_id in enumerate(u_combid):
                # TODO: This was causing problems when multiple science frames
                # were provided to quicklook and the user chose *not* to stack
                # the frames.  But this means it now won't skip processing the
                # B-A pair when the background image(s) are defined.  Punting
                # for now...
#                # Quicklook mode?
#                if self.par['rdx']['quicklook'] and j > 0:
#                    msgs.warn('PypeIt executed in quicklook mode.  Only reducing science frames '
#                              'in the first combination group!')
#                    break
                #
                frames = np.where(self.fitstbl['comb_id'] == comb_id)[0]
                # Find all frames whose comb_id matches the current frames bkg_id.
                #TODO TESTING!!!!
                bg_frames = np.where((np.isin(self.fitstbl['comb_id'], [int(b) for b in self.fitstbl['bkg_id'][frames][0].split(',')]))
                                     & (self.fitstbl['comb_id'] >= 0))[0]
                # JFH changed the syntax below to that above, which allows
                # frames to be used more than once as a background image. The
                # syntax below would require that we could somehow list multiple
                # numbers for the bkg_id which is impossible without a comma
                # separated list
#                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]
                if not self.outfile_exists(frames[0]) or self.overwrite:

                    # # Build history to document what contributd to the reduced
                    # # exposure
                    # history = History(self.fitstbl.frame_paths(frames[0]))
                    # history.add_reduce(calib_ID, self.fitstbl, frames, bg_frames)

                    # TODO -- Should we reset/regenerate self.slits.mask for a new exposure
                    self.reduce_exposure(frames, bg_frames=bg_frames, std_outfile=std_outfile)
                    science_basename[j] = self.basename

                    # # TODO: come up with sensible naming convention for
                    # # save_exposure for combined files
                    # if len(sci_spec2d.detectors) > 0:
                    #     self.save_exposure(frames[0], sci_spec2d, sci_sobjs, self.basename, history)
                    # else:
                    #     msgs.warn('No spec2d and spec1d saved to file because the '
                    #               'calibration/reduction was not successful for all the detectors')

                else:
                    msgs.warn(f'Output file: {self.fitstbl.construct_basename(frames[0])} already '
                              'exists. Set overwrite=True to recreate and overwrite.')

            msgs.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()

    @staticmethod
    def select_detectors(spectrograph, detnum, slitspatnum=None):
        """
        Get the set of detectors to be reduced.

        This is mostly a wrapper for
        :func:`~pypeit.spectrographs.spectrograph.Spectrograph.select_detectors`,
        except that it applies any limitations set by the
        :class:`~pypeit.par.pypeitpar.ReduxPar` parameters.

        Args:
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
                Spectrograph instance that defines the allowed
                detectors/mosaics.
            detnum (:obj:`int`, :obj:`tuple`):
                The detectors/mosaics to parse
            slitspatnum (:obj:`str`, optional):
                Used to restrict the reduction to a specified slit.  See
                :class:`~pypeit.par.pypeitpar.ReduxPar`.

        Returns:
            :obj:`list`: List of unique detectors or detector mosaics to be
            reduced.
        """
        return spectrograph.select_detectors(subset=detnum if slitspatnum is None else slitspatnum)

    def reduce_exposure(self, frames, bg_frames=None, std_outfile=None):
        """
        Reduce a single exposure

        Args:
            frames (:obj:`list`):
                List of 0-indexed rows in :attr:`fitstbl` with the frames to
                reduce.
            bg_frames (:obj:`list`, optional):
                List of frame indices for the background.
            std_outfile (:obj:`str`, optional):
                File with a previously reduced standard spectrum from
                PypeIt.

        """

        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            # TODO: Put this in a try/except block?
            display.clear_all(allow_new=True)

        has_bg = True if bg_frames is not None and len(bg_frames) > 0 else False
        # Is this an b/g subtraction reduction?
        if has_bg:
            self.bkg_redux = True
            # The default is to find_negative objects if the bg_frames are
            # classified as "science", and to not find_negative objects if the
            # bg_frames are classified as "sky". This can be explicitly
            # overridden if par['reduce']['findobj']['find_negative'] is set to
            # something other than the default of None.
            self.find_negative = (('science' in self.fitstbl['frametype'][bg_frames[0]]) |
                                  ('standard' in self.fitstbl['frametype'][bg_frames[0]])) \
                            if self.par['reduce']['findobj']['find_negative'] is None else \
                                self.par['reduce']['findobj']['find_negative']
        else:
            self.bkg_redux = False
            self.find_negative= False

        # Container for all the Spec2DObj
        all_spec2d = spec2dobj.AllSpec2DObj()
        all_spec2d['meta']['bkg_redux'] = self.bkg_redux
        all_spec2d['meta']['find_negative'] = self.find_negative
        # TODO -- Should we reset/regenerate self.slits.mask for a new exposure

        # container for specobjs during first loop (objfind)
        all_specobjs_objfind = specobjs.SpecObjs()
        # container for specobjs during second loop (extraction)
        all_specobjs_extract = specobjs.SpecObjs()
        # list of global_sky obtained during objfind and used in extraction
        initial_sky_list = []
        # list of sciImg
        sciImg_list = []
        # list of bkg_redux_sciimg
        bkg_redux_sciimg_list = []
        # List of detectors with successful calibration
        calibrated_det = []
        # list of successful slits calibrations to be used in the extraction loop
        calib_slits = []
        # List of objFind objects
        objFind_list = []

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
        detectors = self.select_detectors(self.spectrograph, self.par['rdx']['detnum'],
                                          slitspatnum=self.par['rdx']['slitspatnum'])
        msgs.info(f'Detectors to work on: {detectors}')

        # Loop on Detectors -- Calibrate, process image, find objects
        # TODO: Attempt to put in a multiprocessing call here?
        for self.det in detectors:
            msgs.info(f'Reducing detector {self.det}')
            # run calibration
            self.caliBrate = self.calib_one(frames, self.det)
            if not self.caliBrate.success:
                msgs.warn(f'Calibrations for detector {self.det} were unsuccessful!  The step '
                          f'that failed was {self.caliBrate.failed_step}.  Continuing by '
                          f'skipping this detector.')
                continue

            # we save only the detectors that had a successful calibration,
            # and we use only those in the extract loop below
            calibrated_det.append(self.det)
            # we also save the successful slit calibrations because they are used and modified
            # in the slitmask stuff in between the two loops
            calib_slits.append(self.caliBrate.slits)
            # global_sky, skymask and sciImg are needed in the extract loop
            initial_sky, sobjs_obj, sciImg, bkg_redux_sciimg, objFind = self.objfind_one(
                frames, self.det, bg_frames=bg_frames, std_outfile=std_outfile)
            if len(sobjs_obj)>0:
                all_specobjs_objfind.add_sobj(sobjs_obj)
            initial_sky_list.append(initial_sky)
            sciImg_list.append(sciImg)
            bkg_redux_sciimg_list.append(bkg_redux_sciimg)
            objFind_list.append(objFind)

        # slitmask stuff
        if len(calibrated_det) > 0 and self.par['reduce']['slitmask']['assign_obj']:
            # get object positions from slitmask design and slitmask offsets for all the detectors
            spat_flexure = np.array([ss.spat_flexure for ss in sciImg_list])
            # Grab platescale with binning
            bin_spec, bin_spat = parse.parse_binning(self.binning)
            platescale = np.array([ss.detector.platescale*bin_spat for ss in sciImg_list])
            # get the dither offset if available and if desired
            dither_off = None
            if self.par['reduce']['slitmask']['use_dither_offset']:
                if 'dithoff' in self.fitstbl.keys():
                    dither_off = self.fitstbl['dithoff'][frames[0]]

            calib_slits = slittrace.get_maskdef_objpos_offset_alldets(
                all_specobjs_objfind, calib_slits, spat_flexure, platescale,
                self.par['calibrations']['slitedges']['det_buffer'],
                self.par['reduce']['slitmask'], dither_off=dither_off)
            # determine if slitmask offsets exist and compute an average offsets over all the detectors
            calib_slits = slittrace.average_maskdef_offset(
                calib_slits, platescale[0], self.spectrograph.list_detectors(mosaic='MSC' in calib_slits[0].detname))
            # slitmask design matching and add undetected objects
            all_specobjs_objfind = slittrace.assign_addobjs_alldets(
                all_specobjs_objfind, calib_slits, spat_flexure, platescale,
                self.par['reduce']['slitmask'], self.par['reduce']['findobj']['find_fwhm'])

        # Extract
        for i, self.det in enumerate(calibrated_det):
            # re-run (i.e., load) calibrations
            self.caliBrate = self.calib_one(frames, self.det)
            self.caliBrate.slits = calib_slits[i]

            detname = sciImg_list[i].detector.name

            # TODO: pass back the background frame, pass in background
            # files as an argument. extract one takes a file list as an
            # argument and instantiates science within
            if all_specobjs_objfind.nobj > 0:
                all_specobjs_on_det = all_specobjs_objfind[all_specobjs_objfind.DET == detname]
            else:
                all_specobjs_on_det = all_specobjs_objfind

            # Extract
            all_spec2d[detname], tmp_sobjs \
                    = self.extract_one(frames, self.det, sciImg_list[i], bkg_redux_sciimg_list[i], objFind_list[i],
                                       initial_sky_list[i], all_specobjs_on_det)
            # Hold em
            if tmp_sobjs.nobj > 0:
                all_specobjs_extract.add_sobj(tmp_sobjs)

            # Add calibration associations to the SpecObjs object
            all_specobjs_extract.calibs = calibrations.Calibrations.get_association(
                                    self.fitstbl, self.spectrograph, self.calibrations_path,
                                    self.fitstbl[frames[0]]['setup'],
                                    self.fitstbl.find_frame_calib_groups(frames[0])[0], self.det,
                                    must_exist=True, proc_only=True)

            # JFH TODO write out the background frame?

            # TODO -- Save here?  Seems like we should.  Would probably need to use update_det=True

        if len(all_spec2d.detectors) > 0:
            # Build history to document what contributd to the reduced
            # exposure
            history = History(self.fitstbl.frame_paths(frames[0]))
            history.add_reduce(self.fitstbl.find_frame_calib_groups(frames[0])[0], self.fitstbl, frames, bg_frames)
            self.save_exposure(frames[0], all_spec2d, all_specobjs_extract, self.basename, history,
                               skip_write_2d=self.par['scienceframe']['process']['skip_write_2d'])
        else:
            msgs.warn('No spec2d and spec1d saved to file because the '
                      'calibration/reduction was not successful for all the detectors')


    def get_sci_metadata(self, frame, det):
        """
        Grab the meta data for a given science frame and specific detector

        Args:
            frame (int): Frame index
            det (int): Detector index

        Returns:
            5 objects are returned::
                - str: Object type;  science or standard
                - str: Setup/configuration string
                - `astropy.time.Time`_: Time of observation
                - str: Basename of the frame
                - str: Binning of the detector

        """

        # Set binning, obstime, basename, and objtype
        binning = self.fitstbl['binning'][frame]
        obstime  = self.fitstbl.construct_obstime(frame)
        basename = self.fitstbl.construct_basename(frame, obstime=obstime)
        types  = self.fitstbl['frametype'][frame].split(',')
        if 'science' in types:
            objtype_out = 'science'
        elif 'standard' in types:
            objtype_out = 'standard'
        else:
            msgs.error('get_sci_metadata() should only be run on standard or science frames.  '
                       f'Types of this frame are: {types}')
        calib_key = CalibFrame.construct_calib_key(self.fitstbl['setup'][frame],
                                                   self.fitstbl['calib'][frame],
                                                   self.spectrograph.get_det_name(det))
        return objtype_out, calib_key, obstime, basename, binning

    def calib_one(self, frames, det, stop_at_step:str=None):
        """
        Run Calibration for a single exposure/detector pair

        Args:
            frames (:obj:`list`):
                List of frames (rows) to calibrate
                Only used to idetify the setup and calibration group
            det (:obj:`int`):
                Detector number (1-indexed)
            stop_at_step (:obj:`str`, optional):
                Run only up to this calibration step.
                

        Returns:
            caliBrate (:class:`pypeit.calibrations.Calibrations`)

        """

        msgs.info(f'Building/loading calibrations for detector {det}')
        # Instantiate Calibrations class
        user_slits = slittrace.merge_user_slit(self.par['rdx']['slitspatnum'],
                                               self.par['rdx']['maskIDs'])
        caliBrate = calibrations.Calibrations.get_instance(
            self.fitstbl, self.par['calibrations'], self.spectrograph,
            self.calibrations_path, qadir=self.qa_path,
            reuse_calibs=self.reuse_calibs, show=self.show, user_slits=user_slits,
            chk_version=self.par['rdx']['chk_version'])

        # Check
        if stop_at_step is not None and stop_at_step not in caliBrate.steps:
            msgs.error(f"Requested stop_at_step={stop_at_step} is not a valid calibration step.\n Allowed steps are: {caliBrate.steps}")
            
        # These need to be separate to accomodate COADD2D
        caliBrate.set_config(frames[0], det, self.par['calibrations'])

        # Run
        caliBrate.run_the_steps(stop_at_step=stop_at_step)

        return caliBrate

    def objfind_one(self, frames, det, bg_frames=None, std_outfile=None):
        """
        Reduce + Find Objects in a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Parameters
        ----------
        frames : :obj:`list`
            List of frames to extract; stacked if more than one is provided
        det : :obj:`int`
            Detector number (1-indexed)
        bg_frames : :obj:`list`, optional
            List of frames to use as the background. Can be empty.
        std_outfile : :obj:`str`, optional
            Filename for the standard star spec1d file. Passed directly to
            :func:`~pypeit.specobjs.get_std_trace`.

        Returns
        -------
        global_sky : `numpy.ndarray`_
            Initial global sky model
        sobjs_obj : :class:`~pypeit.specobjs.SpecObjs`
            List of objects found
        sciImg : :class:`~pypeit.images.pypeitimage.PypeItImage`
            Science image
        bkg_redux_sciimg : :class:`~pypeit.images.pypeitimage.PypeItImage`
            Science image before background subtraction
            if self.bkg_redux is True, otherwise None.
            It's used to generate a global sky model without bkg subtraction.
        objFind : :class:`~pypeit.find_objects.FindObjects`
            Object finding speobject

        """
        # Grab some meta-data needed for the reduction from the fitstbl
        self.objtype, self.setup, self.obstime, self.basename, self.binning \
                = self.get_sci_metadata(frames[0], det)

        msgs.info("Object finding begins for {} on det={}".format(self.basename, det))

        # Is this a standard star?
        self.std_redux = self.objtype == 'standard'
        frame_par = self.par['calibrations']['standardframe'] \
                        if self.std_redux else self.par['scienceframe']
        # Get the standard trace if need be

        if self.std_redux is False and std_outfile is not None:
            std_trace = specobjs.get_std_trace(self.spectrograph.get_det_name(det), std_outfile)
        else:
            std_trace = None

        # Build Science image
        sci_files = self.fitstbl.frame_paths(frames)
        sciImg = buildimage.buildimage_fromlist(
            self.spectrograph, det, frame_par,
            sci_files, bias=self.caliBrate.msbias, bpm=self.caliBrate.msbpm,
            dark=self.caliBrate.msdark,
            scattlight=self.caliBrate.msscattlight,
            flatimages=self.caliBrate.flatimages,
            slits=self.caliBrate.slits,  # For flexure correction
            ignore_saturation=False)

        # get no bkg subtracted sciImg to generate a global sky model without bkg subtraction.
        # it's a dictionary with only `image` and `ivar` keys if bkg_redux=False, otherwise it's None
        bkg_redux_sciimg = None

        # Background Image?
        if bg_frames is not None and len(bg_frames) > 0:
            # get no bkg subtracted sciImg
            bkg_redux_sciimg = sciImg
            # Build the background image
            bg_file_list = self.fitstbl.frame_paths(bg_frames)
            # TODO I think we should create a separate self.par['bkgframe'] parameter set to hold the image
            # processing parameters for the background frames.  This would allow the user to specify different
            # parameters for the background frames than for the science frames.  
            bgimg = buildimage.buildimage_fromlist(self.spectrograph, det, frame_par, bg_file_list,
                                                   bpm=self.caliBrate.msbpm,
                                                   bias=self.caliBrate.msbias,
                                                   dark=self.caliBrate.msdark,
                                                   scattlight=self.caliBrate.msscattlight,
                                                   flatimages=self.caliBrate.flatimages,
                                                   slits=self.caliBrate.slits,
                                                   ignore_saturation=False)

            # NOTE: If the spatial flexure exists for sciImg, the subtraction
            # function propagates that to the subtracted image, ignoring any
            # spatial flexure determined for the background image.
            sciImg = bkg_redux_sciimg.sub(bgimg)

        # Flexure
        spat_flexure = None
        # use the flexure correction in the "shift" column
        manual_flexure = self.fitstbl[frames[0]]['shift']
        if (self.objtype == 'science' and self.par['scienceframe']['process']['spat_flexure_correct']) or \
                (self.objtype == 'standard' and self.par['calibrations']['standardframe']['process']['spat_flexure_correct']) or \
                    manual_flexure:
            if (manual_flexure or manual_flexure == 0) and not (np.issubdtype(self.fitstbl[frames[0]]["shift"], np.integer)):
                msgs.info(f'Implementing manual flexure of {manual_flexure}')
                spat_flexure = np.float64(manual_flexure)
                sciImg.spat_flexure = spat_flexure
            else:
                msgs.info(f'Using auto-computed flexure')
                spat_flexure = sciImg.spat_flexure
        msgs.info(f'Flexure being used is: {spat_flexure}')
        # Build the initial sky mask
        initial_skymask = self.load_skyregions(initial_slits=self.spectrograph.pypeline != 'SlicerIFU',
                                               scifile=sciImg.files[0], frame=frames[0], spat_flexure=spat_flexure)

        # Deal with manual extraction
        row = self.fitstbl[frames[0]]
        manual_obj = ManualExtractionObj.by_fitstbl_input(
            row['filename'], row['manual'], self.spectrograph) if len(row['manual'].strip()) > 0 else None

        # Instantiate Reduce object
        # Required for pypeline specific object
        # At instantiaton, the fullmask in self.sciImg is modified
        objFind = find_objects.FindObjects.get_instance(sciImg, self.caliBrate.slits,
                                                        self.spectrograph, self.par, self.objtype,
                                                        wv_calib=self.caliBrate.wv_calib,
                                                        waveTilts=self.caliBrate.wavetilts,
                                                        initial_skymask=initial_skymask,
                                                        bkg_redux=self.bkg_redux,
                                                        manual=manual_obj,
                                                        find_negative=self.find_negative,
                                                        std_redux=self.std_redux,
                                                        show=self.show,
                                                        basename=self.basename)

        # Do it
        initial_sky, sobjs_obj = objFind.run(std_trace=std_trace, show_peaks=self.show)

        # Return
        return initial_sky, sobjs_obj, sciImg, bkg_redux_sciimg, objFind

    def load_skyregions(self, initial_slits=False, scifile=None, frame=None, spat_flexure=None):
        """
        Generate or load sky regions, if defined by the user.

        Sky regions are defined by the internal provided parameters; see
        ``user_regions`` in :class:`~pypeit.par.pypeitpar.SkySubPar`.  If
        included in the pypeit file like so,

        .. code-block:: ini

            [reduce]
                [[skysub]]
                    user_regions = :25,75:

        The first and last 25% of all slits are used as sky.  If the user has
        used the ``pypeit_skysub_regions`` GUI to generate a sky mask for a
        given frame, this can be searched for and loaded by setting the
        parameter to ``user``:

        .. code-block:: ini

            [reduce]
                [[skysub]]
                    user_regions = user

        Parameters
        ----------
        initial_slits : :obj:`bool`, optional
            Flag to use the initial slits before any tweaking based on the
            slit-illumination profile; see
            :func:`~pypeit.slittrace.SlitTraceSet.select_edges`.
        scifile : :obj:`str`, optional
            The file name used to define the user-based sky regions.  Only used
            if ``user_regions = user``.
        frame : :obj:`int`, optional
            The index of the frame used to construct the calibration key.  Only
            used if ``user_regions = user``.
        spat_flexure : :obj:`float`, None, optional
            The spatial flexure (measured in pixels) of the science frame relative to the trace frame.

        Returns
        -------
        skymask : `numpy.ndarray`_
            A boolean array used to select sky pixels; i.e., True is a pixel
            that corresponds to a sky region.  If the ``user_regions`` parameter
            is not set (or an empty string), the returned value is None.
        """
        if self.par['reduce']['skysub']['user_regions'] in [None, '']:
            return None

        # First priority given to user_regions first
        if self.par['reduce']['skysub']['user_regions'] == 'user':
            # Build the file name
            calib_key = CalibFrame.construct_calib_key(
                                self.fitstbl['setup'][frame],
                                CalibFrame.ingest_calib_id(self.fitstbl['calib'][frame]),
                                self.spectrograph.get_det_name(self.det))
            regfile = buildimage.SkyRegions.construct_file_name(calib_key,
                                                                calib_dir=self.calibrations_path,
                                                                basename=io.remove_suffix(scifile))
            regfile = Path(regfile).absolute()
            if not regfile.exists():
                msgs.error(f'Unable to find SkyRegions file: {regfile} . Create a SkyRegions '
                           'frame using pypeit_skysub_regions, or change the user_regions to '
                           'the percentage format.  See documentation.')
            msgs.info(f'Loading SkyRegions file: {regfile}')
            return buildimage.SkyRegions.from_file(regfile).image.astype(bool)

        skyregtxt = self.par['reduce']['skysub']['user_regions']
        if isinstance(skyregtxt, list):
            skyregtxt = ",".join(skyregtxt)
        msgs.info(f'Generating skysub mask based on the user defined regions: {skyregtxt}')
        # NOTE : Do not include spatial flexure here!
        #        It is included when generating the mask in the return statement below
        slits_left, slits_right, _ \
            = self.caliBrate.slits.select_edges(initial=initial_slits, flexure=None)

        maxslitlength = np.max(slits_right-slits_left)
        # Get the regions
        status, regions = skysub.read_userregions(skyregtxt, self.caliBrate.slits.nslits, maxslitlength)
        if status == 1:
            msgs.error("Unknown error in sky regions definition. Please check the value:" + msgs.newline() + skyregtxt)
        elif status == 2:
            msgs.error("Sky regions definition must contain a percentage range, and therefore must contain a ':'")
        # Generate and return image
        return skysub.generate_mask(self.spectrograph.pypeline, regions, self.caliBrate.slits,
                                    slits_left, slits_right, spat_flexure=spat_flexure)

    def extract_one(self, frames, det, sciImg, bkg_redux_sciimg, objFind, initial_sky, sobjs_obj):
        """
        Extract Objects in a single exposure/detector pair

        sci_ID and det need to have been set internally prior to calling this method

        Args:
            frames (:obj:`list`):
                List of frames to extract; stacked if more than one
                is provided
            det (:obj:`int`):
                Detector number (1-indexed)
            sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Data container that holds a single image from a
                single detector and its related images (e.g. ivar, mask)
            bkg_redux_sciimg (:class:`~pypeit.images.pypeitimage.PypeItImage`, optional):
                Data container that holds a single image from a
                single detector and its related images (e.g. ivar, mask)
                before background subtraction if self.bkg_redux is True,
                otherwise None. It's used to generate a global sky
                model without bkg subtraction.
            objFind : :class:`~pypeit.find_objects.FindObjects`
                Object finding object
            initial_sky (`numpy.ndarray`_):
                Initial global sky model
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                List of objects found during `run_objfind`

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
        # Is this a standard star?
        self.std_redux = 'standard' in self.objtype

        ## TODO JFH I think all of this about determining the final global sky should be moved out of this method
        ## and preferably into the FindObjects class. I see why we are doing it like this since for multislit we need
        # to find all of the objects first using slitmask meta data,  but this comes at the expense of a much more complicated
        # control structure.

        # Update the global sky
        skymask = None
        if 'standard' in self.fitstbl['frametype'][frames[0]] or \
                self.par['reduce']['findobj']['skip_skysub'] or \
                self.par['reduce']['findobj']['skip_final_global'] or \
                self.par['reduce']['skysub']['user_regions'] is not None:
            final_global_sky = initial_sky
        else:
            # Update the skymask
            skymask = objFind.create_skymask(sobjs_obj)
            final_global_sky = objFind.global_skysub(previous_sky=initial_sky,
                                                     skymask=skymask, show=self.show,
                                                     reinit_bpm=False)
        # get the bkg_redux_global_sky
        bkg_redux_global_sky = None
        if self.bkg_redux:
            skymask = objFind.create_skymask(sobjs_obj) if skymask is None else skymask
            # DO NOT reinit_bpm, nor update_crmask
            bkg_redux_global_sky = objFind.global_skysub(skymask=skymask, bkg_redux_sciimg=bkg_redux_sciimg,
                                                     reinit_bpm=False, update_crmask=False, show=self.show)

        scaleImg = objFind.scaleimg

        # Each spec2d file includes the slits object with unique flagging
        #  for extraction failures.  So we make a copy here before those flags
        #  are modified.
        maskdef_designtab = self.caliBrate.slits.maskdef_designtab
        slits = copy.deepcopy(self.caliBrate.slits)
        slits.maskdef_designtab = None

        # update here slits.mask since global_skysub modify reduce_bpm and we need to propagate it into extraction
        flagged_slits = np.where(objFind.reduce_bpm)[0]
        if len(flagged_slits) > 0:
            slits.mask[flagged_slits] = \
                slits.bitmask.turn_on(slits.mask[flagged_slits], 'BADSKYSUB')

        if not self.par['reduce']['extraction']['skip_extraction']:
            msgs.info(f"Extraction begins for {self.basename} on det={det}")
            # Instantiate Reduce object
            # Required for pipeline specific object
            # At instantiation, the fullmask in self.sciImg is modified
            # TODO Are we repeating steps in the init for FindObjects and Extract??
            self.exTract = extraction.Extract.get_instance(
                sciImg, slits, sobjs_obj, self.spectrograph,
                self.par, self.objtype, global_sky=final_global_sky, bkg_redux_global_sky=bkg_redux_global_sky,
                waveTilts=self.caliBrate.wavetilts, wv_calib=self.caliBrate.wv_calib, flatimages=self.caliBrate.flatimages,
                bkg_redux=self.bkg_redux, return_negative=self.par['reduce']['extraction']['return_negative'],
                std_redux=self.std_redux, basename=self.basename, show=self.show)
            # Perform the extraction
            skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts, slits = self.exTract.run()
            slitgpm = np.logical_not(self.exTract.extract_bpm)
            slitshift = self.exTract.slitshift
        else:
            msgs.info(f"Extraction skipped for {self.basename} on det={det}")
            # Since the extraction was not performed, fill the arrays with the best available information
            skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts = \
                final_global_sky, \
                bkg_redux_global_sky, \
                np.zeros_like(objFind.sciImg.image), \
                np.copy(objFind.sciImg.ivar), \
                objFind.sciImg.fullmask, \
                sobjs_obj, \
                objFind.waveimg, \
                objFind.tilts
            slitgpm = (slits.mask == 0)
            slitshift = objFind.slitshift
            # If waveImg has not yet been created, make it now
            if waveImg is None:
                waveImg = self.caliBrate.wv_calib.build_waveimg(tilts, slits, spat_flexure=objFind.spat_flexure_shift)

        # Apply a reference frame correction to each object and the waveimg
        vel_corr, waveImg = self.refframe_correct(slits, self.fitstbl["ra"][frames[0]], self.fitstbl["dec"][frames[0]],
                                                  self.obstime, slitgpm=slitgpm, waveimg=waveImg, sobjs=sobjs)

        # TODO -- Do this upstream
        # Tack on wavelength RMS
        for sobj in sobjs:
            iwv = np.where(self.caliBrate.wv_calib.spat_ids == sobj.SLITID)[0][0]
            sobj.WAVE_RMS =self.caliBrate.wv_calib.wv_fits[iwv].rms

        # Construct table of spectral flexure
        spec_flex_table = Table()
        spec_flex_table['spat_id'] = slits.spat_id
        spec_flex_table['sci_spec_flexure'] = slitshift

        # Construct the Spec2DObj
        spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                        ivarraw=sciImg.ivar,
                                        skymodel=skymodel,
                                        bkg_redux_skymodel=bkg_redux_skymodel,
                                        objmodel=objmodel,
                                        ivarmodel=ivarmodel,
                                        scaleimg=scaleImg,
                                        waveimg=waveImg,
                                        bpmmask=outmask,
                                        detector=sciImg.detector,
                                        sci_spat_flexure=sciImg.spat_flexure,
                                        sci_spec_flexure=spec_flex_table,
                                        vel_corr=vel_corr,
                                        vel_type=self.par['calibrations']['wavelengths']['refframe'],
                                        tilts=tilts,
                                        slits=slits,
                                        wavesol=self.caliBrate.wv_calib.wave_diagnostics(print_diag=False),
                                        maskdef_designtab=maskdef_designtab)
        spec2DObj.process_steps = sciImg.process_steps

        spec2DObj.calibs = calibrations.Calibrations.get_association(
                                    self.fitstbl, self.spectrograph, self.calibrations_path,
                                    self.fitstbl[frames[0]]['setup'],
                                    self.fitstbl.find_frame_calib_groups(frames[0])[0], det,
                                    must_exist=True, proc_only=True)

        # QA
        spec2DObj.gen_qa()

        # Return
        return spec2DObj, sobjs

    # TODO :: Should this be moved outside of this class?
    def refframe_correct(self, slits, ra, dec, obstime, slitgpm=None, waveimg=None, sobjs=None):
        """ Correct the calibrated wavelength to the user-supplied reference frame

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit trace set object
            ra (float, str):
                Right Ascension
            dec (float, str):
                Declination
            obstime (`astropy.time.Time`_):
                Observation time
            slitgpm (`numpy.ndarray`_, None, optional):
                1D boolean array indicating the good slits (True). If None, the gpm will be taken from slits
            waveimg (`numpy.ndarray`_, optional)
                Two-dimensional image specifying the wavelength of each pixel
            sobjs (:class:`~pypeit.specobjs.SpecObjs`, None, optional):
                Spectrally extracted objects

        """
        if slitgpm is None:
            slitgpm = (slits.mask == 0)
        # Correct Telescope's motion
        refframe = self.par['calibrations']['wavelengths']['refframe']
        vel_corr = 0.0
        if refframe in ['heliocentric', 'barycentric'] \
                and self.par['calibrations']['wavelengths']['reference'] != 'pixel':
            msgs.info("Performing a {0} correction".format(self.par['calibrations']['wavelengths']['refframe']))
            # Calculate correction
            radec = ltu.radec_to_coord((ra, dec))
            vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   refframe)
            # Apply correction to objects
            msgs.info('Applying {0} correction = {1:0.5f} km/s'.format(refframe, vel))
            if (sobjs is not None) and (sobjs.nobj != 0):
                # Loop on slits to apply
                gd_slitord = slits.slitord_id[slitgpm]
                for slitord in gd_slitord:
                    indx = sobjs.slitorder_indices(slitord)
                    this_specobjs = sobjs[indx]
                    # Loop on objects
                    for specobj in this_specobjs:
                        if specobj is None:
                            continue
                        specobj.apply_helio(vel_corr, refframe)

            # Apply correction to wavelength image
            if waveimg is not None:
                waveimg *= vel_corr
        else:
            msgs.info('A wavelength reference frame correction will not be performed.')

        # Return the value of the correction and the corrected wavelength image
        return vel_corr, waveimg

    def save_exposure(self, frame:int, all_spec2d:spec2dobj.AllSpec2DObj,
                      all_specobjs:specobjs.SpecObjs, basename:str, history:History=None,
                      skip_write_2d:bool=False):
        """
        Save the outputs from extraction for a given exposure

        Args:
            frame (:obj:`int`):
                0-indexed row in the metadata table with the frame
                that has been reduced.
            all_spec2d(:class:`~pypeit.spec2dobj.AllSpec2DObj`):
                The 2D reduced spectrum objects.
            all_specobjs (:class:`~pypeit.specobjs.SpecObjs`):
                The 1D spectral extraction objects.
            basename (:obj:`str`):
                The root name for the output file.
            history (:class:`~pypeit.history.History`), optional:
                History entries to be added to fits header
            skip_write_2d (:obj:`bool`), optional:
                Skip writing the 2D spectrum to disk
        """
        # TODO: Need some checks here that the exposure has been reduced?

        # Determine the headers
        row_fitstbl = self.fitstbl[frame]
        # Need raw file header information
        rawfile = self.fitstbl.frame_paths(frame)
        head2d = fits.getheader(rawfile, ext=self.spectrograph.primary_hdrext)

        # Check for the directory
        if not self.science_path.is_dir():
            self.science_path.mkdir()

        # NOTE: There are some gymnastics here to keep from altering
        # self.par['rdx']['detnum'].  I.e., I can't just set update_det =
        # self.par['rdx']['detnum'] because that can alter the latter if I don't
        # deepcopy it...
        if self.par['rdx']['detnum'] is None:
            update_det = None
        elif isinstance(self.par['rdx']['detnum'], list):
            update_det = [self.spectrograph.allowed_mosaics.index(d)+1 
                            if isinstance(d, tuple) else d for d in self.par['rdx']['detnum']]
        else:
            update_det = self.par['rdx']['detnum']

        subheader = self.spectrograph.subheader_for_spec(row_fitstbl, head2d)
        # 1D spectra
        if all_specobjs.nobj > 0 and not self.par['reduce']['extraction']['skip_extraction']:
            # Spectra
            outfile1d = self.science_path / f'spec1d_{basename}.fits'
            # TODO
            #embed(header='deal with the following for maskIDs;  713 of pypeit')
            all_specobjs.write_to_fits(subheader, outfile1d,
                                       update_det=update_det,
                                       slitspatnum=self.par['rdx']['slitspatnum'],
                                       history=history)
            # Info
            outfiletxt = self.science_path / f'spec1d_{basename}.txt'
            # TODO: Note we re-read in the specobjs from disk to deal with situations where
            # only a single detector is run in a second pass but in the same reduction directory.
            # Thiw was to address Issue #1116 in PR #1154. Slightly inefficient, but only other
            # option is to re-work write_info to also "append"
            sobjs = specobjs.SpecObjs.from_fitsfile(outfile1d, chk_version=False)
            sobjs.write_info(outfiletxt, self.spectrograph.pypeline)
            #all_specobjs.write_info(outfiletxt, self.spectrograph.pypeline)

        if skip_write_2d:
            return

        # 2D spectra
        outfile2d = self.science_path / f'spec2d_{basename}.fits'
        # Build header
        pri_hdr = all_spec2d.build_primary_hdr(head2d, self.spectrograph,
                                               redux_path=self.par['rdx']['redux_path'],
                                               calib_dir=self.caliBrate.calib_dir,
                                               subheader=subheader,
                                               history=history)

        # Write
        all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr,
                                 update_det=update_det,
                                 slitspatnum=self.par['rdx']['slitspatnum'])

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
        msgs.info(utils.get_time_string(time.perf_counter()-self.tstart))

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


class NIRSpecSlitPypeIt(PypeIt):
    """
    PypeIt class for performing JWST NIRSpec calibration.
    See :class:`PypeIt` for arguments.
    """

    def get_sci_metadata(self, frame, det, slit_name=None):
        """
        Get the metadata for a science frame

        Args:
            frame (:obj:`int`):
                0-indexed row in the metadata table with the frame
                that has been reduced.
            det (:obj:`int`):
                Detector number (1-indexed)
            slit_name (:obj:`str`, optional):
                Name of the slit, if applicable.

        Returns:
            tuple: The object type, setup, obstime, basename, and binning
            for the science frame.
        """

        objtype, setup, obstime, basename, binning = super().get_sci_metadata(frame,det)
        # update basename for NIRSpec
        basename = basename.replace('_assign_wcs', '')
        basename = basename.replace('_nrs1', '')
        basename = basename.replace('_nrs2', '')
        if slit_name is not None:
            basename = basename.replace(self.spectrograph.camera, f'{slit_name}_{self.spectrograph.camera}')
            setup = setup + f'_{slit_name}'
        return objtype, setup, obstime, basename, binning

    @staticmethod
    def get_slit_slice(cal_data, slit_name):
        """
        Get the slice for a slit from the provided calibration data.

        Args:
            cal_data (:class:`~jwst.datamodels.MultiSlitModel`):
                The input NIRSpec MultiSlit data model from the _cal.fits output file
            slit_name (:obj:`str'):
                The name of the slit to get the slice for, e.g. 'S200A1' or 'S200B1'.

        Returns:
            :obj:`slice`: The slice for the specified slit.
        """
        # Check if the slit name is present in the cal_data
        slit_names = np.array([slit.name for slit in cal_data.slits])
        if slit_name not in slit_names:
            msgs.error(f'Slit name {slit_name} not found in the calibration data {cal_data.meta.filename}')

        # Get the index of the slit
        indx = np.where(slit_names == slit_name)[0][0]

        # get the slice for the slit
        this_slit = cal_data.slits[indx]
        spec_lo = this_slit.xstart -1
        spec_hi = spec_lo + this_slit.xsize
        spat_lo = this_slit.ystart - 1
        spat_hi = spat_lo + this_slit.ysize
        return slice(spec_lo, spec_hi), slice(spat_lo, spat_hi)

    def build_image(self, frames, _det, _detectors, _sci_data, _slit_slices, kludge_err=1.2):
        # ronoise = np.array([copy.deepcopy(d['ronoise']) for d in _detectors])
        # # Get the non-linear count level
        # if self.rawflatimg.is_mosaic:
        #     # if this is a mosaic we take the maximum value among all the detectors
        #     nonlinear_counts = np.max(
        #         [rawdets.nonlinear_counts() for rawdets in self.rawflatimg.detector.detectors])
        # else:
        #     nonlinear_counts = self.rawflatimg.detector.nonlinear_counts()

        headarr = copy.deepcopy(self.spectrograph.get_headarr(fits.open(self.fitstbl.frame_paths(frames)[0])))
        ronoise = np.array([copy.deepcopy(d['ronoise']) for d in _detectors])
        exptime = _sci_data[0].meta.exposure.effective_exposure_time
        # dark = np.array([np.prod(parse.parse_binning(d.binning)) * d.darkcurr* exptime / 3600. for d in _detectors])
        bpm = self.caliBrate.msbpm
        nimg = self.caliBrate.nimg
        img_scale = utils.inverse(self.caliBrate.flatimages.pixelflat_norm)

        if nimg == 1:
            sci = _sci_data[0].data.T[_slit_slices[0]] * exptime
            var_rnoise = kludge_err**2 * _sci_data[0].var_rnoise.T[_slit_slices[0]] * exptime ** 2
            var_poisson = kludge_err**2 * _sci_data[0].var_poisson.T[_slit_slices[0]] * exptime ** 2
            datasec = np.full_like(sci, 1, dtype=int)
            rn2img = procimg.rn2_frame(datasec, ronoise[0], units='e-', gain=1.)
            # darkim = np.repeat(dark, np.prod(sci.shape)).reshape(sci.shape)

        else:
            sci_list = [_sci_data[d].data.T[_slit_slices[d]] * exptime for d in range(_sci_data.size)]
            sci = self.spectrograph.make_mosaic(sci_list, _det, _slit_slices)
            var_rnoise_list = [_sci_data[d].var_rnoise.T[_slit_slices[d]] * exptime ** 2 for d in range(_sci_data.size)]
            var_rnoise = self.spectrograph.make_mosaic(var_rnoise_list, _det, _slit_slices)
            var_poisson_list = [_sci_data[d].var_poisson.T[_slit_slices[d]] * exptime ** 2 for d in
                                range(_sci_data.size)]
            var_poisson = self.spectrograph.make_mosaic(var_poisson_list, _det, _slit_slices)
            datasec_list = [np.full_like(sci_list[d], 1, dtype=int) for d in range(_sci_data.size)]
            datasec = self.spectrograph.make_mosaic(datasec_list, _det, _slit_slices)
            rn2img_list = [procimg.rn2_frame(datasec_list[d], ronoise[d], units='e-', gain=1.) for d in
                           range(_sci_data.size)]
            rn2img = self.spectrograph.make_mosaic(rn2img_list, _det, _slit_slices)
            # darkim_list = [np.repeat(dark[d], np.prod(sci_list[d])).reshape(sci_list[d].shape) for d in range(_sci_data.size)]
            # darkim = self.spectrograph.make_mosaic(darkim_list, _det, _slit_slices)

        sci2, flat_bpm = flat.flatfield(sci, self.caliBrate.flatimages.pixelflat_norm)
        var_rnoise2, _ = flat.flatfield(var_rnoise, self.caliBrate.flatimages.pixelflat_norm ** 2)
        var_poisson2, _ = flat.flatfield(var_poisson, self.caliBrate.flatimages.pixelflat_norm ** 2)
        # bpm |= flat_bpm
        base_var = procimg.base_variance(var_rnoise2, count_scale=None)
        var = procimg.variance_model(base_var, counts=var_poisson2, count_scale=None,
                                     noise_floor=self.par['scienceframe']['process']['noise_floor'])

        # make zero the values outside the slit and the nans
        sci2[(self.caliBrate.flatimages.pixelflat_norm==1) | np.logical_not(np.isfinite(sci2))] = 0.
        var[(self.caliBrate.flatimages.pixelflat_norm==1) | np.logical_not(np.isfinite(var))] = 0.
        base_var[(self.caliBrate.flatimages.pixelflat_norm==1) | np.logical_not(np.isfinite(base_var))] = 0.
        img_scale[(self.caliBrate.flatimages.pixelflat_norm==1) | np.logical_not(np.isfinite(img_scale))] = 0.
        bpm[(self.caliBrate.flatimages.pixelflat_norm==1) | np.logical_not(np.isfinite(sci2))] = 1

        pypeitImage = pypeitimage.PypeItImage(sci2, ivar=utils.inverse(var), amp_img=datasec,
                                              det_img=datasec, rn2img=rn2img, base_var=base_var,
                                              img_scale=img_scale,
                                              detector=self.spectrograph.get_mosaic_par(_det) if nimg == 2 else
                                              _detectors[0],
                                              PYP_SPEC=self.spectrograph.name,
                                              units='e-', exptime=exptime,
                                              noise_floor=self.par['scienceframe']['process']['noise_floor'],
                                              shot_noise=self.par['scienceframe']['process']['shot_noise'],
                                              bpm=bpm.astype(bool),
                                              filename=None)
        pypeitImage.rawheadlist = headarr

        return pypeitImage

    def reduce_exposure(self, frames, bg_frames=None, std_outfile=None):
        """
        Reduce a single exposure

        Args:
            frames (:obj:`list`):
                List of 0-indexed rows in :attr:`fitstbl` with the frames to
                reduce.
            bg_frames (:obj:`list`, optional):
                List of frame indices for the background.
            std_outfile (:obj:`str`, optional):
                File with a previously reduced standard spectrum from
                PypeIt.

        """

        # if show is set, clear the ginga channels at the start of each new sci_ID
        if self.show:
            # TODO: Put this in a try/except block?
            display.clear_all(allow_new=True)

        has_bg = True if bg_frames is not None and len(bg_frames) > 0 else False
        if has_bg and len(bg_frames) % len(frames) != 0:
            msgs.error('Background frames must be provided for all science frames. '
                       'Please provide a number of background frames that is a multiple of the number of science frames.')
        # Is this an b/g subtraction reduction?
        if has_bg:
            self.bkg_redux = True
            # The default is to find_negative objects if the bg_frames are
            # classified as "science", and to not find_negative objects if the
            # bg_frames are classified as "sky". This can be explicitly
            # overridden if par['reduce']['findobj']['find_negative'] is set to
            # something other than the default of None.
            self.find_negative = (('science' in self.fitstbl['frametype'][bg_frames[0]]) |
                                  ('standard' in self.fitstbl['frametype'][bg_frames[0]])) \
                            if self.par['reduce']['findobj']['find_negative'] is None else \
                                self.par['reduce']['findobj']['find_negative']
        else:
            self.bkg_redux = False
            self.find_negative= False

        # Find the detectors to reduce
        detnums = self.par['rdx']['detnum']
        msgs.info(f'Detectors to work on: {detnums}')

        # book keeping for successful calibrations
        if len(detnums) == 2:
            # check that both detectors are present, i.e., both files exist
            both_exist = np.all(np.any(s in fname for fname in self.fitstbl[frames]['filename'].data) for s in ['_nrs1', '_nrs2'])
            if not both_exist:
                msgs.error('Reduction requested for both detectors, but they are not both present in the input files.')
        elif len(detnums) == 1 and detnums[0] == 1:
            # check that the first detector is present, i.e., the file exists
            if not np.any('_nrs1' in fname for fname in self.fitstbl[frames]['filename'].data):
                msgs.error('Reduction requested for detector 1, but it is not present in the input files.')
        elif len(detnums) == 1 and detnums[0] == 2:
            # check that the second detector is present, i.e., the file exists
            if not np.any('_nrs2' in fname for fname in self.fitstbl[frames]['filename'].data):
                msgs.error('Reduction requested for detector 2, but it is not present in the input files.')

        # TODO: use jwst.exp_to_source.exp_to_source to parse the jwst models
        # Create arrays to hold JWST spec2
        # NOTE: all this implies that the number of frames is the same as the number of detectors,
        # i.e., we are not combining different exposures for the main science frame but we can combine the bg frames
        sci_data = np.array(datamodels.open(self.fitstbl.frame_paths(frames)))
        if self.bkg_redux:
            sci_data_bkg = np.array(datamodels.open(self.fitstbl.frame_paths(bg_frames)))
        # find the calibration files for this frame
        # NOTE: we assume that calib_id is the same for all `frames` and `bg_frames`
        calib_grps = self.fitstbl.find_frame_calib_groups(frames[0])
        cal_files = self.fitstbl.find_frame_files('trace', calib_ID=calib_grps)
        if len(cal_files) != frames.size:
            msgs.error(f'Number of _cal files ({len(cal_files)}) does not match number of science frames ({frames.size}) for calib_ID={calib_grps}.')
        cal_data = np.array(datamodels.open(cal_files))
        # find the flat field files for this frame
        flat_files = self.fitstbl.find_frame_files('pixelflat', calib_ID=calib_grps)
        if len(flat_files) != frames.size:
            msgs.error(f'Number of _flat files ({len(flat_files)}) does not match number of science frames ({frames.size}) for calib_ID={calib_grps}.')
        flat_data = np.array(datamodels.open(flat_files))
        # If the flat field files have interpolatedflat_fs slits, append them to the flat_data
        for i, _flat_data in enumerate(flat_data):
            fs_path = flat_files[i].replace('_interpolatedflat', '_interpolatedflat_fs')
            if Path(fs_path).exists():
                msgs.info(f'Appending interpolatedflat FS slits into MOS output for {Path(flat_files[i]).name}')
                _flat_data_fs = datamodels.open(fs_path)
                for slit in _flat_data_fs.slits:
                    _flat_data.slits.append(slit)


        # If the background frames are provided, we need to skip the sky-subtraction
        if self.bkg_redux:
            self.par['reduce']['findobj']['skip_skysub'] = True
            self.par['reduce']['extraction']['skip_optimal'] = True  # Skip local_skysubtraction and profile fitting

        # get slits and sources names
        # we slit slit_names between NRS1 and NRS2 because we need them later to determine which detector to use,
        # but we don't need to do that for the sources
        slit_names_nrs1 = np.array([[slit.name for slit in cal_data[i].slits] for i in range(cal_data.size) if
                                    cal_data[i] is not None and cal_data[i].meta.instrument.detector == 'NRS1']).flatten()
        slit_names_nrs2 = np.array([[slit.name for slit in cal_data[i].slits] for i in range(cal_data.size) if
                                    cal_data[i] is not None and cal_data[i].meta.instrument.detector == 'NRS2']).flatten()
        slit_names = np.hstack([slit_names_nrs1, slit_names_nrs2])
        source_names = np.hstack([[slit.source_name for slit in cal_data[i].slits] for i in range(cal_data.size) if cal_data[i] is not None])
        source_ids = np.hstack([[slit.source_id for slit in cal_data[i].slits] for i in range(cal_data.size) if cal_data[i] is not None])
        source_aliases = np.hstack([[slit.source_alias for slit in cal_data[i].slits] for i in range(cal_data.size) if cal_data[i] is not None])
        # Find the unique slit names and the unique sources aligned with those slits
        slit_names_uni, uni_indx = np.unique(slit_names, return_index=True)
        source_names_uni = source_names[uni_indx]
        source_ids_uni = source_ids[uni_indx]
        source_aliases_uni = source_aliases[uni_indx]
        slit_sources_uni = [(slit, source_name, source_id, source_alias)
                            for slit, source_name, source_id, source_alias in
                            zip(slit_names_uni, source_names_uni, source_ids_uni, source_aliases_uni)]


        if self.par['rdx']['maskIDs'] is not None:
            maskIDs = [str(mid).strip() for mid in self.par['rdx']['maskIDs']]
            gd_slits_sources = [(slt, src_n, src_id, src_alias)
                                for slt, src_n, src_id, src_alias in slit_sources_uni
                                if str(slt).strip() in maskIDs or str(src_id) in maskIDs or str(src_alias) in maskIDs]
            if not gd_slits_sources:
                msgs.warn(f'No slits or sources found for maskIDs={self.par["rdx"]["maskIDs"]}. '
                          'Skipping reduction.')
                return
            # find if there are maskIDs that are not present in the slit_sources_uni
            elif len(gd_slits_sources) < len(maskIDs):
                missing_maskIDs = [mid for mid in maskIDs if mid not in [slt for slt, _, _, _ in gd_slits_sources] and
                                  mid not in [src_id for _, _, src_id, _ in gd_slits_sources] and
                                  mid not in [src_alias for _, _, _, src_alias in gd_slits_sources]]
                msgs.warn(f'The following maskIDs were not found: {", ".join(missing_maskIDs)}. '
                          'Reduction will be performed on the available slits and sources.')

        else:
            gd_slits_sources = slit_sources_uni

        # MSGS info on which slits and sources are being reduced
        msgs.info(f'Reducing the following (slit_name, src_name): '
                  f'{", ".join([f"({slt}, {src_n})" for slt, src_n, _, _ in gd_slits_sources])}')

        for ii, (islit, isource, isource_id, isource_alias) in enumerate(gd_slits_sources):
            # Print status message
            add_to_msgs = f'Slit name: {islit}' if isource is None else f'SRC name: {isource}'
            msgs_string = f'Reducing target {self.fitstbl['target'][frames[0]]} - {add_to_msgs}' + msgs.newline()

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

            # Get the detector, cal, flat, science data, and slit slices
            if islit in slit_names_nrs1 and islit in slit_names_nrs2:
                msgs.error('Needto fix this part of the code to reduce slits that are on both NRS1 and NRS2 and have multiple background frames.')
                # # this is a mosaic
                # _det = (1, 2)
                # _detectors = self.spectrograph.get_mosaic_par(_det).detectors
                # _cal_data = np.array([cal_data[[cal_d.meta.instrument.detector == 'NRS1' for cal_d in cal_data]][0],
                #              cal_data[[cal_d.meta.instrument.detector == 'NRS2' for cal_d in cal_data]][0]])
                # _slit_slices = [self.get_slit_slice(_cal_data[0], islit), self.get_slit_slice(_cal_data[1], islit)]
                # _flat_data = np.array([flat_data[[flat_d.meta.instrument.detector == 'NRS1' for flat_d in flat_data]][0],
                #                 flat_data[[flat_d.meta.instrument.detector == 'NRS2' for flat_d in flat_data]][0]])
                # _sci_data = np.array([sci_data[[sci_d.meta.instrument.detector == 'NRS1' for sci_d in sci_data]][0],
                #                 sci_data[[sci_d.meta.instrument.detector == 'NRS2' for sci_d in sci_data]][0]])
                # if self.bkg_redux:
                #     _sci_data_bkg = np.array([sci_data_bkg[[sci_d.meta.instrument.detector == 'NRS1' for sci_d in sci_data_bkg]][0],
                #                     sci_data_bkg[[sci_d.meta.instrument.detector == 'NRS2' for sci_d in sci_data_bkg]][0]])
            elif islit in slit_names_nrs1:
                # this is NRS1
                _det = 1
                _detectors = np.array([self.spectrograph.get_detector_par(_det)])
                _cal_data = cal_data[[cal_d.meta.instrument.detector == 'NRS1' for cal_d in cal_data]]
                _slit_slices = [self.get_slit_slice(_cal_data[0], islit)]
                _flat_data = flat_data[[flat_d.meta.instrument.detector == 'NRS1' for flat_d in flat_data]]
                _sci_data = sci_data[[sci_d.meta.instrument.detector == 'NRS1' for sci_d in sci_data]]
                if self.bkg_redux:
                    _sci_data_bkg = sci_data_bkg[[sci_d.meta.instrument.detector == 'NRS1' for sci_d in sci_data_bkg]]
            elif islit in slit_names_nrs2:
                # this is NRS2
                _det = 2
                _detectors = np.array([self.spectrograph.get_detector_par(_det)])
                _cal_data = cal_data[[cal_d.meta.instrument.detector == 'NRS2' for cal_d in cal_data]]
                _slit_slices = [self.get_slit_slice(_cal_data[0], islit)]
                _flat_data = flat_data[[flat_d.meta.instrument.detector == 'NRS2' for flat_d in flat_data]]
                _sci_data = sci_data[[sci_d.meta.instrument.detector == 'NRS2' for sci_d in sci_data]]
                if self.bkg_redux:
                    _sci_data_bkg = sci_data_bkg[[sci_d.meta.instrument.detector == 'NRS2' for sci_d in sci_data_bkg]]
            else:
                msgs.error(f'Slit {islit} not found in the slit names of the calibration files. '
                           'Please check the slit names in the calibration files or the maskIDs provided.')
            msgs.info(f'Reducing detector {_det}')

            self.objtype, self.setup, self.obstime, self.basename, self.binning = self.get_sci_metadata(frames[0], _det, slit_name=islit)

            # Loop on detectors to get the calibrations
            # from pypeitdev.jwst.jwst_utils import NIRSpecSlitCalibrations, jwst_mosaic, jwst_reduce
            from pypeit.images import combineimage
            _calibrate = []
            # for d, _det in enumerate(detectors):
            #     _calibrate.append(NIRSpecSlitCalibrations(self.spectrograph.get_detector_par(_det), cal_data[d], flat_data[d],
            #                                     islit, f070_f100_rescale=None))
            # THIS IS JUST TO NOT MAKE THE CODE COMPLAIN
            self.caliBrate = calibrations.Calibrations.get_instance(self.fitstbl, self.par['calibrations'],
                                                                    self.spectrograph,self.calibrations_path,
                                                                    qadir=self.qa_path, reuse_calibs=self.reuse_calibs,
                                                                    show=self.show, user_slits=islit,
                                                                    chk_version=self.par['rdx']['chk_version'],
                                                                    cal_data=_cal_data, flat_data=_flat_data,
                                                                    slit_slices=_slit_slices)
            self.caliBrate.set_config(frames[0], _det, self.par['calibrations'])
            self.caliBrate.run_the_steps()

            # build the science PypeItImage
            sciImg = self.build_image(frames, _det, _detectors, _sci_data, _slit_slices)
            if self.bkg_redux:
                sciImg_bkg_list = []
                for _bg in _sci_data_bkg:
                    sciImg_bkg = self.build_image(bg_frames, _det, _detectors, [_bg], _slit_slices)
                    sciImg_bkg_list.append(sciImg_bkg)
                # combine the background images
                combineImage = combineimage.CombineImage(sciImg_bkg_list, self.par['scienceframe']['process'])
                comb_sciImg_bkg = combineImage.run(ignore_saturation=True)
                sciImg = sciImg.sub(comb_sciImg_bkg)

            # REDUCE
            # Is this a standard star?
            self.std_redux = self.objtype == 'standard'
            frame_par = self.par['calibrations']['standardframe'] \
                if self.std_redux else self.par['scienceframe']

            # NO FLEXURE FOR JWST (?)
            # Build the initial sky mask
            initial_skymask = self.load_skyregions(initial_slits=self.spectrograph.pypeline != 'SlicerIFU',
                                                   scifile=None, frame=frames[0], spat_flexure=None)

            # Deal with manual extraction
            row = self.fitstbl[frames[0]]
            manual_obj = ManualExtractionObj.by_fitstbl_input(
                row['filename'], row['manual'], self.spectrograph) if len(row['manual'].strip()) > 0 else None

            # Instantiate Reduce object
            # Required for pypeline specific object
            # At instantiaton, the fullmask in self.sciImg is modified
            objFind = find_objects.FindObjects.get_instance(sciImg, self.caliBrate.slits,
                                                            self.spectrograph, self.par, self.objtype,
                                                            tilts=self.caliBrate.tilts,
                                                            initial_skymask=initial_skymask,
                                                            bkg_redux=self.bkg_redux,
                                                            manual=manual_obj,
                                                            find_negative=self.find_negative,
                                                            std_redux=self.std_redux,
                                                            show=self.show,
                                                            basename=self.basename)

            # Do it
            initial_sky, sobjs_obj = objFind.run(std_trace=None, show_peaks=self.show)

            if 'standard' in self.fitstbl['frametype'][frames[0]] or \
                    self.par['reduce']['findobj']['skip_skysub'] or \
                    self.par['reduce']['findobj']['skip_final_global'] or \
                    self.par['reduce']['skysub']['user_regions'] is not None:
                final_global_sky = initial_sky
            else:
                # Update the skymask
                skymask = objFind.create_skymask(sobjs_obj)
                final_global_sky = objFind.global_skysub(previous_sky=initial_sky,
                                                         skymask=skymask, show=self.show,
                                                         reinit_bpm=False)
            scaleImg = objFind.scaleimg

            # Each spec2d file includes the slits object with unique flagging
            #  for extraction failures.  So we make a copy here before those flags
            #  are modified.
            maskdef_designtab = self.caliBrate.slits.maskdef_designtab
            slits = copy.deepcopy(self.caliBrate.slits)
            slits.maskdef_designtab = None

            # update here slits.mask since global_skysub modify reduce_bpm and we need to propagate it into extraction
            flagged_slits = np.where(objFind.reduce_bpm)[0]
            if len(flagged_slits) > 0:
                slits.mask[flagged_slits] = \
                    slits.bitmask.turn_on(slits.mask[flagged_slits], 'BADSKYSUB')

            if not self.par['reduce']['extraction']['skip_extraction']:
                msgs.info(f"Extraction begins for {self.basename} on det={_det}")
                # Instantiate Reduce object
                # Required for pipeline specific object
                # At instantiation, the fullmask in self.sciImg is modified
                # TODO Are we repeating steps in the init for FindObjects and Extract??
                self.exTract = extraction.Extract.get_instance(
                    sciImg, slits, sobjs_obj, self.spectrograph,
                    self.par, self.objtype, global_sky=final_global_sky, bkg_redux_global_sky=None,
                    tilts=self.caliBrate.tilts, waveimg=self.caliBrate.waveimg,
                    flatimages=self.caliBrate.flatimages,
                    bkg_redux=self.bkg_redux, return_negative=self.par['reduce']['extraction']['return_negative'],
                    std_redux=self.std_redux, basename=self.basename, show=self.show)
                # Perform the extraction
                skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts, slits = self.exTract.run()
                slitshift = self.exTract.slitshift
            else:
                msgs.info(f"Extraction skipped for {self.basename} on det={_det}")
                # Since the extraction was not performed, fill the arrays with the best available information
                skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs, waveImg, tilts = \
                    final_global_sky, \
                        None, \
                        np.zeros_like(objFind.sciImg.image), \
                        np.copy(objFind.sciImg.ivar), \
                        objFind.sciImg.fullmask, \
                        sobjs_obj, \
                        objFind.waveimg, \
                        objFind.tilts
                slitshift = objFind.slitshift

            # Construct table of spectral flexure
            spec_flex_table = Table()
            spec_flex_table['spat_id'] = slits.spat_id
            spec_flex_table['sci_spec_flexure'] = slitshift

            # Construct the Spec2DObj
            spec2DObj = spec2dobj.Spec2DObj(sciimg=sciImg.image,
                                            ivarraw=sciImg.ivar,
                                            skymodel=skymodel,
                                            bkg_redux_skymodel=bkg_redux_skymodel,
                                            objmodel=objmodel,
                                            ivarmodel=ivarmodel,
                                            scaleimg=scaleImg,
                                            waveimg=waveImg,
                                            bpmmask=outmask,
                                            detector=sciImg.detector,
                                            sci_spat_flexure=sciImg.spat_flexure,
                                            sci_spec_flexure=spec_flex_table,
                                            vel_corr=None,
                                            vel_type=self.par['calibrations']['wavelengths']['refframe'],
                                            tilts=tilts,
                                            slits=slits,
                                            wavesol=None,
                                            maskdef_designtab=maskdef_designtab)
            spec2DObj.process_steps = sciImg.process_steps

            spec2DObj.calibs = calibrations.Calibrations.get_association(
                self.fitstbl, self.spectrograph, self.calibrations_path,
                self.fitstbl[frames[0]]['setup'],
                self.fitstbl.find_frame_calib_groups(frames[0])[0], _det,
                must_exist=True, proc_only=True)

            # QA
            spec2DObj.gen_qa()

            # Container for all the Spec2DObj
            all_spec2d = spec2dobj.AllSpec2DObj()
            all_spec2d['meta']['bkg_redux'] = self.bkg_redux
            all_spec2d['meta']['find_negative'] = self.find_negative
            # Container for the specobjs
            all_specobjs = specobjs.SpecObjs()

            # fill the spec2DObj
            all_spec2d[sciImg.detector.name] = spec2DObj
            # fill the specobjs
            if sobjs.nobj > 0:
                all_specobjs.add_sobj(sobjs)

            if len(all_spec2d.detectors) > 0:
                # Build history to document what contributd to the reduced
                # exposure
                history = History(self.fitstbl.frame_paths(frames[0]))
                history.add_reduce(self.fitstbl.find_frame_calib_groups(frames[0])[0], self.fitstbl, frames, bg_frames)
                self.save_exposure(frames[0], all_spec2d, all_specobjs, self.basename, history,
                                   skip_write_2d=self.par['scienceframe']['process']['skip_write_2d'])
            else:
                msgs.warn('No spec2d and spec1d saved to file because the '
                          'calibration/reduction was not successful for all the detectors')

