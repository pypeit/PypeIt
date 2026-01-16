"""
Main driver class for PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path
import time
import os
import datetime

from IPython import embed

import numpy as np


from pypeit import inputfiles
from pypeit.core import qa
from pypeit import msgs
from pypeit import calibrations
from pypeit import utils
from pypeit.history import History
from pypeit.metadata import PypeItMetaData
from pypeit import outputfiles
from pypeit import exposure

from pypeit import pypeit_steps

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
    def __init__(self, pypeit_file, verbosity=2, overwrite=True, reuse_calibs=False, logname=None,
                 show=False, redux_path=None, calib_only=False):

        # Set up logging
        self.logname = logname
        self.verbosity = verbosity
        self.pypeit_file = pypeit_file

        # State
        #self.run_state = state.RunPypeItState(pypeit_file=pypeit_file, 
        #                                      current_step='init',
        #                                      current_det=-1,
        #                                      current_calibID=-1)
        #self.run_state = self.run_state.load()
        self.run_state = None
        
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
        #self.basename = None
        self.obstime = None

    @property
    def science_path(self) -> Path:
        """Return the path to the science directory."""
        return outputfiles.science_path(self.par)

    @property
    def qa_path(self) -> str:
        """Return the path to the top-level QA directory."""
        return os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])

    def build_qa(self):
        """
        Generate QA wrappers

        Called by run_pypeit.py
        """
        msgs.qa_path = self.qa_path
        qa.gen_qa_dir(self.qa_path)
        qa.gen_mf_html(self.pypeit_file, self.qa_path)
        qa.gen_exp_html()

    def calib_all(self):
        """
        Process all calibration frames.

        Provides an avenue to process the calibrations for a dataset 
        without (or omitting) any science/standard frames.
        """
        self.tstart = time.perf_counter()

        # Frame indices
        for calib_ID in self.fitstbl.calib_groups:
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)
            if not any(in_grp):
                continue
            # Find the detectors to reduce
            detectors = self.spectrograph.select_detectors(subset=self.par['rdx']['detnum'] if self.par['rdx']['slitspatnum'] is None 
                                              else self.par['rdx']['slitspatnum'])
            msgs.info(f'Detectors to work on: {detectors}')

            # Loop on Detectors
            for self.det in detectors:
                msgs.info(f'Working on detector {self.det}')

                caliBrate = pypeit_steps.calib_one(self.spectrograph, self.fitstbl, self.par,
                                       self.det, calib_ID, self.calibrations_path)
                                       

        # Finish
        self.print_end_time()

    def reduce_all(self):
        """
        Main driver of the end-to-end reduction

        Calibration and extraction via a series of calls to
        :func:`reduce_exposure`.

        """
        # Validate the parameter set
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

#        # Find the standard frames
#        is_standard = self.fitstbl.find_frames('standard')
#        if np.any(is_standard):
#            msgs.info(f'Found {np.sum(is_standard)} standard frames to reduce.')
#
#        # Find the science frames
#        is_science = self.fitstbl.find_frames('science')
#        if np.any(is_science):
#            msgs.info(f'Found {np.sum(is_science)} science frames to reduce.')
#
#        # This will give an error to alert the user that no reduction will be
#        # run if there are no science/standard frames and `run_pypeit` is run
#        # without -c flag
#        if not np.any(is_science) and not np.any(is_standard):
#            msgs.error('No science/standard frames provided. Add them to your PypeIt file '
#                       'if this is a standard run! Otherwise run calib_only reduction using -c flag')
#
#        # Frame indices
#        frame_indx = np.arange(len(self.fitstbl))

        # ############################################################################
        # Standard Star(s) Loop
        # ############################################################################
        # Iterate over each calibration group and reduce the standards
        for calib_ID in self.fitstbl.calib_groups:

            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                           calib_ID, self.calibrations_path,
                           reduce_standard=True, overwrite=self.overwrite,
                           show=self.show, 
                           run_state=self.run_state,
                           reuse_calibs=self.reuse_calibs)

        # ############################################################################
        # Science Frame(s) Loop
        # ############################################################################
        # Iterate over each calibration group again and reduce the science frames
        for calib_ID in self.fitstbl.calib_groups:
            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                                        calib_ID, self.calibrations_path,
                                        reduce_standard=False, overwrite=self.overwrite,
                                        show=self.show, run_state=self.run_state,
                                        reuse_calibs=self.reuse_calibs)
            msgs.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()


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

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


def reduce_calibID(spectrograph, par, fitstbl, calib_ID:str, 
                   calibrations_path:str,
                   reduce_standard:bool=False, overwrite:bool=False,
                   show:bool=False,
                   run_state=None,
                   reuse_calibs:bool=True):

    """
    Reduce all the frames in a given calibration group

    Outputs are written to disk.

    Calls :func:`~pypeit.exposure.reduce_exposure` to do the
    actual reduction.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph object for the instrument being reduced.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction, including slitmask and
            object finding parameters.
        fitstbl (:class:`pypeit.metadata.PypeItMetaData`):
            The metadata table for the current reduction.
        calib_ID (:obj:`str`):
            The calibration group ID to reduce.
        calibrations_path (:obj:`str`):
            The path to the calibration files.
        reduce_standard (:obj:`bool`, optional):
            Reduce the standard frames if True; science frames if
            False.
        overwrite (:obj:`bool`, optional):
            Overwrite any existing files.
        show (:obj:`bool`, optional):
            Show reduction steps via plots (which will block further
            execution until clicked on) and outputs to ginga. Requires
            remote control ginga session via
            ``ginga --modules=RC,SlitWavelength &``
        run_state (:class:`~pypeit.state.RunPypeItState`, optional):
            The current state of the reduction.
        reuse_calibs (:obj:`bool`, optional):
            Reuse any pre-existing calibration files
    """

    if reduce_standard:
        is_this = fitstbl.find_frames('standard')
        rtype = 'standard'
    else:
        is_this = fitstbl.find_frames('science')
        rtype = 'science'

    # Frame indices
    frame_indx = np.arange(len(fitstbl))

    # Find all the frames in this calibration group
    in_grp = fitstbl.find_calib_group(calib_ID)

    if not np.any(is_this & in_grp):
        return

    # Find the indices of the science frames in this calibration group:
    grp_this = frame_indx[is_this & in_grp]
    msgs.info(f'Found {len(grp_this)} {rtype} frames in calibration group {calib_ID}.')

    # Associate standards (previously reduced above) for this setup
    if not reduce_standard:
        is_standard = fitstbl.find_frames('standard')
        std_outfile = outputfiles.get_std_outfile(fitstbl, par, frame_indx[is_standard])
    else:
        std_outfile = None

    # Loop on unique comb_id
    u_combid = np.unique(fitstbl['comb_id'][grp_this])

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
        frames = np.where(fitstbl['comb_id'] == comb_id)[0]
        # Find all frames whose comb_id matches the current frames bkg_id.
        bg_frames = np.where((fitstbl['comb_id'] == fitstbl['bkg_id'][frames][0])
                                & (fitstbl['comb_id'] >= 0))[0]
        # JFH changed the syntax below to that above, which allows
        # frames to be used more than once as a background image. The
        # syntax below would require that we could somehow list multiple
        # numbers for the bkg_id which is impossible without a comma
        # separated list
#                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]

        outfile2d = outputfiles.spec_output_file(fitstbl, par,
                                            frames[0], twod=True)
        if not outfile2d.is_file() or overwrite:

            # Build history to document what contributd to the reduced
            # exposure
            history = History(fitstbl.frame_paths(frames[0]))
            history.add_reduce(calib_ID, fitstbl, frames, bg_frames)

            # TODO -- Should we reset/regenerate self.slits.mask for a new exposure
            #sci_spec2d, sci_sobjs = self.reduce_exposure(
            #    frames, calib_ID, bg_frames=bg_frames, 
            #    std_outfile=std_outfile)

            this_spec2d, this_sobjs = exposure.reduce_exposure(
                spectrograph, fitstbl, par, frames, calib_ID, 
                calibrations_path, bg_frames=bg_frames,
                reuse_calibs=reuse_calibs, run_state=run_state,
                show=show,
                std_outfile=std_outfile)

            # TODO: come up with sensible naming convention for
            # save_exposure for combined files
            if len(this_spec2d.detectors) > 0:
                exposure.save_exposure(spectrograph,
                                    fitstbl, par, frames[0], 
                                    this_spec2d, this_sobjs, calibrations_path,
                                    history=history,
                                    skip_write_2d=par['scienceframe']['process']['skip_write_2d'])
            else:
                msgs.warn('No spec2d and spec1d saved to file because the '
                            'calibration/reduction was not successful for all the detectors')
        else:
            msgs.warn(f'Output file: {fitstbl.construct_basename(frames[0])} already '
                        'exists. Set overwrite=True to recreate and overwrite.')
