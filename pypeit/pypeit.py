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
from pypeit import qa
from pypeit import log
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
        overwrite (:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
        reuse_calibs (:obj:`bool`, optional):
            Reuse any pre-existing calibration files
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

    @staticmethod
    def get_instance(pypeit_file, **kwargs):
        """
        Instantiate and return the :class:`PypeIt` subclass appropriate for
        the provided spectrograph.

        Args:
            pypeit_file (:obj:`str`):
                PypeIt filename.
            **kwargs:
                Passed to the constructor.

        Returns:
            :class:`PypeIt`: The appropriate subclass instance.
        """
        spectrograph = inputfiles.PypeItFile.from_file(pypeit_file).get_spectrograph()
        if spectrograph.pypeline == 'NIRSpecSlit':
            return NIRSpecSlitPypeIt(pypeit_file, **kwargs)
        return PypeIt(pypeit_file, **kwargs)

    def __init__(
        self, pypeit_file, overwrite=True, reuse_calibs=False, show=False, redux_path=None,
        calib_only=False
    ):

        # Set up logging
        self.pypeit_file = pypeit_file

        # State
        #self.run_state = state.RunPypeItState(pypeit_file=pypeit_file,
        #                                      current_step='init',
        #                                      current_det=-1,
        #                                      current_calibID=-1)
        #self.run_state = self.run_state.load()
        # TODO: Implement RunPypeItState.
        self.run_state = None
        
        # Load up PypeIt file
        self.pypeItFile = inputfiles.PypeItFile.from_file(pypeit_file)
        self.calib_only = calib_only

        # Build the spectrograph and the parameters
        self.spectrograph, self.par, config_specific_file = self.pypeItFile.get_pypeitpar()
        log.info(f'Loaded spectrograph {self.spectrograph.name}')
        log.info('Setting configuration-specific parameters using '
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
        log.info('Compiling metadata')
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
            self._check_calibrations()

        # --------------------------------------------------------------
        #   - Write .calib file (For QA naming amongst other things)
        calib_file = pypeit_file.replace('.pypeit', '.calib')
        calibrations.Calibrations.association_summary(calib_file, self.fitstbl, self.spectrograph,
                                                      self.calibrations_path, overwrite=True)

        # Report paths
        log.info('Setting reduction path to {0}'.format(self.par['rdx']['redux_path']))
        log.info('Calibration frames saved to: {0}'.format(self.calibrations_path))
        log.info('Science data output to: {0}'.format(self.science_path))
        log.info('Quality assessment plots output to: {0}'.format(self.qa_path))

        # Init
        self.det = None
        self.tstart = None
        #self.basename = None
        self.obstime = None

    def _check_calibrations(self):
        """Check that required calibration frames are available."""
        calibrations.check_for_calibs(self.par, self.fitstbl,
                                      raise_error=self.par['calibrations']['raise_chk_error'])

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
#        log.qa_path = self.qa_path
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

        log.info('Starting calibration only reduction')

        # Frame indices
        for calib_ID in self.fitstbl.calib_groups:
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)
            if not any(in_grp):
                continue
            log.info(f'Working on calibration group {calib_ID}')
            # Find the detectors to reduce
            detectors = self.spectrograph.select_detectors(subset=self.par['rdx']['detnum'] if self.par['rdx']['slitspatnum'] is None
                                              else self.par['rdx']['slitspatnum'])
            log.info(f'Detectors to work on: {detectors}')

            # Loop on Detectors
            for self.det in detectors:
                log.info(f'Working on detector {self.det}')

                caliBrate = pypeit_steps.calib_one(self.spectrograph, self.fitstbl, self.par,
                                       self.det, calib_ID, self.calibrations_path)

            log.info(f'Finished calibration group {calib_ID}')


        # Finish
        self.print_end_time()

    def reduce_all(self):
        """
        Main driver of the end-to-end reduction

        Calibration and extraction via a series of calls to
        :func:`reduce_calibID`.
        """
        # Validate the parameter set
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

        # ############################################################################
        # Standard Star(s) Loop
        # ############################################################################
        if np.any(self.fitstbl.find_frames('standard')):
            log.info('Starting standard star reduction')
        # Iterate over each calibration group and reduce the standards
        for calib_ID in self.fitstbl.calib_groups:
            # check first if there is any standard stars frame in this group
            if not np.any(self.fitstbl.find_frames('standard', calib_ID)):
                continue
            log.info(f'Working on calibration group {calib_ID}')
            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                           calib_ID, self.calibrations_path,
                           reduce_standard=True, overwrite=self.overwrite,
                           show=self.show,
                           run_state=self.run_state,
                           reuse_calibs=self.reuse_calibs)
            log.info(f'Finished calibration group {calib_ID}')

        # ############################################################################
        # Science Frame(s) Loop
        # ############################################################################
        if np.any(self.fitstbl.find_frames('science')):
            log.info('Starting science frame reduction')
        # Iterate over each calibration group again and reduce the science frames
        for calib_ID in self.fitstbl.calib_groups:
            # check first if there is any science frames frame in this group
            if not np.any(self.fitstbl.find_frames('science', calib_ID)):
                continue
            log.info(f'Working on calibration group {calib_ID}')
            reduce_calibID(self.spectrograph, self.par, self.fitstbl,
                                        calib_ID, self.calibrations_path,
                                        reduce_standard=False, overwrite=self.overwrite,
                                        show=self.show, run_state=self.run_state,
                                        reuse_calibs=self.reuse_calibs)
            log.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()

    def print_end_time(self):
        """
        Print the elapsed time
        """
        # Capture the end time and print it to user
        log.info(utils.get_time_string(time.perf_counter()-self.tstart))

    def __repr__(self):
        # Generate sets string
        return '<{:s}: pypeit_file={}>'.format(self.__class__.__name__, self.pypeit_file)


class NIRSpecSlitPypeIt(PypeIt):
    """
    Child of :class:`PypeIt` for JWST NIRSpec slit-by-slit reductions.

    This class skips the standard calibration checks that don't apply to
    NIRSpec and overrides :func:`reduce_all` to call
    :func:`reduce_calibID` for each slit/calibration-group pair, passing
    ``slitname`` and ``detnum`` to select the NIRSpec reduction path.

    The slit loop lives here (outermost), wrapping around the calibration
    group loops.  :func:`reduce_calibID` calls the unified
    :func:`~pypeit.exposure.reduce_exposure` which handles both
    standard and NIRSpec pipelines via the ``slitname`` parameter.

    See :class:`PypeIt` for arguments.
    """

    def _check_calibrations(self):
        """NIRSpec uses a different calibration workflow -- skip standard checks."""
        pass

    def calib_all(self):
        """
        Process all calibration frames for all NIRSpec slits.

        Mirrors the slit-loop structure of :func:`reduce_all` but calls
        :func:`~pypeit.pypeit_steps.calib_one` for each slit instead of
        running the full science reduction.  The slit list and detector
        mapping are obtained once per calibration group via
        :meth:`~pypeit.spectrographs.jwst_nirspec.JWSTNIRSpecSpectrograph.get_nirspec_slits`.

        Calibration files are written to disk for each slit independently.
        """
        self.tstart = time.perf_counter()

        log.info('Starting calibration only reduction')

        for calib_ID in self.fitstbl.calib_groups:
            log.info(f'Working on calibration group {calib_ID}')
            gd_slits_sources, slit_det_map = self.spectrograph.get_nirspec_slits(
                self.fitstbl, calib_ID, self.par)
            if gd_slits_sources is None:
                log.warning(f'No slits found for calibration group {calib_ID}. Skipping.')
                continue

            for ii, (islit, isource, isource_id, isource_alias) in enumerate(gd_slits_sources):
                det = slit_det_map[islit]
                log.info(f'Working on detector {det}, slit {islit} ({ii + 1}/{len(gd_slits_sources)})')
                pypeit_steps.calib_one(
                    self.spectrograph, self.fitstbl, self.par,
                    det, calib_ID, self.calibrations_path,
                    slitname=islit, reuse_calibs=self.reuse_calibs,
                    qa_path=self.qa_path, show=self.show,
                    run_state=self.run_state)
            log.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()

    def reduce_all(self):
        """
        NIRSpec-specific reduction driver.

        Determines slits to reduce via
        :func:`~pypeit.spectrographs.jwst_nirspec.JWSTNIRSpecSpectrograph.get_nirspec_slits`,
        then iterates over slits.  For each slit the standard-star and
        science calibration-group loops call
        :func:`reduce_calibID` (which internally calls the unified
        :func:`~pypeit.exposure.reduce_exposure` with ``slitname``
        and ``detnum`` set for the NIRSpec slit-by-slit path).
        """
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

        # ################################################################
        # Standard Star(s) Loop
        # ################################################################
        if np.any(self.fitstbl.find_frames('standard')):
            log.info('Starting standard star reduction')
        for calib_ID in self.fitstbl.calib_groups:
            # check first if there is any standard stars frame in this group
            if not np.any(self.fitstbl.find_frames('standard', calib_ID)):
                continue
            log.info(f'Working on calibration group {calib_ID}')
            # get ids of the slits that we want to reduce, could be all (but still the reduction will be separated per slit),
            # or user input slits passed through the parameter self.par['rdx']['maskIDs']
            gd_slits_sources, slit_det_map = self.spectrograph.get_nirspec_slits(self.fitstbl, calib_ID, self.par)
            if gd_slits_sources is None:
                continue
            for ii, (islit, isource, isource_id, isource_alias) in enumerate(gd_slits_sources):
                obj_str = f'Slit name: {islit}' if isource is None \
                    else f'SRC name: {isource}'
                log.info(f'Processing {obj_str} ({ii + 1}/{len(gd_slits_sources)})')

                reduce_calibID(
                    self.spectrograph, self.par, self.fitstbl,
                    calib_ID, self.calibrations_path,
                    slitname=islit, detnum=slit_det_map[islit],
                    reduce_standard=True, overwrite=self.overwrite,
                    show=self.show, run_state=self.run_state, reuse_calibs=self.reuse_calibs)
            log.info(f'Finished calibration group {calib_ID}')

        # ################################################################
        # Science Frame(s) Loop
        # ################################################################
        if np.any(self.fitstbl.find_frames('science')):
            log.info('Starting science frame reduction')
        for calib_ID in self.fitstbl.calib_groups:
            # check first if there is any science frames frame in this group
            if not np.any(self.fitstbl.find_frames('science', calib_ID)):
                continue
            log.info(f'Working on calibration group {calib_ID}')
            # get ids of the slits that we want to reduce, could be all (but still the reduction will be separated per slit),
            # or user input slits passed through the parameter self.par['rdx']['maskIDs']
            gd_slits_sources, slit_det_map = self.spectrograph.get_nirspec_slits(self.fitstbl, calib_ID, self.par)
            if gd_slits_sources is None:
                continue
            for ii, (islit, isource, isource_id, isource_alias) in enumerate(gd_slits_sources):
                obj_str = f'Slit name: {islit}' if isource is None \
                    else f'SRC name: {isource}'
                log.info(f'Processing {obj_str} ({ii + 1}/{len(gd_slits_sources)})')

                reduce_calibID(
                    self.spectrograph, self.par, self.fitstbl,
                    calib_ID, self.calibrations_path,
                    slitname=islit, detnum=slit_det_map[islit],
                    reduce_standard=False, overwrite=self.overwrite,
                    show=self.show, run_state=self.run_state, reuse_calibs=self.reuse_calibs)
            log.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()


def reduce_calibID(spectrograph, par, fitstbl, calib_ID:str,
                   calibrations_path:str,
                   reduce_standard:bool=False, overwrite:bool=False,
                   show:bool=False,
                   run_state=None,
                   reuse_calibs:bool=True,
                   slitname:str=None, detnum=None):
    """
    Reduce all the frames in a given calibration group.

    This function drives the reduction for both the standard pipeline
    (MultiSlit, Echelle, SlicerIFU) and the JWST NIRSpec slit-by-slit
    pipeline.  When ``slitname`` is ``None`` (the default), the standard
    pipeline path is used.  When ``slitname`` is provided, the
    NIRSpec slit-by-slit path is used.  In both cases,
    :func:`~pypeit.exposure.reduce_exposure` is called, which internally
    dispatches based on ``slitname``.

    Outputs are written to disk.

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
        slitname (:obj:`str`, optional):
            NIRSpec slit name.  When provided, the NIRSpec slit-by-slit
            reduction path is used.  Default is None (standard pipeline).
        detnum (:obj:`int` or :obj:`tuple`, optional):
            Detector number(s) to reduce.  Used only for NIRSpec to
            restrict the reduction to the detector hosting the
            requested slit.  Default is None.
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
    log.info(f'Found {len(grp_this)} {rtype} frames in calibration group {calib_ID}.')

    # Associate standards (previously reduced above) for this setup
    if not reduce_standard:
        is_standard = fitstbl.find_frames('standard')
        std_outfile = outputfiles.get_std_outfile(fitstbl, par, frame_indx[is_standard],
                                                  slitname=slitname)
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
#                    log.warning('PypeIt executed in quicklook mode.  Only reducing science frames '
#                              'in the first combination group!')
#                    break
        #
        frames = np.where(fitstbl['comb_id'] == comb_id)[0]
        # Find all frames whose comb_id matches the current frames bkg_id.
        # bg_frames = np.where((fitstbl['comb_id'] == fitstbl['bkg_id'][frames][0])
        #                         & (fitstbl['comb_id'] >= 0))[0]
        # NOTE: The comma-separated parsing handles both single bkg_id values
        # (standard pipeline) and comma-separated lists (NIRSpec).
        bg_frames = np.where(
            (np.isin(fitstbl['comb_id'],
                     [int(b) for b in str(fitstbl['bkg_id'][frames][0]).split(',')]))
            & (fitstbl['comb_id'] >= 0))[0]
        # JFH changed the syntax below to that above, which allows
        # frames to be used more than once as a background image. The
        # syntax below would require that we could somehow list multiple
        # numbers for the bkg_id which is impossible without a comma
        # separated list
#                bg_frames = np.where(self.fitstbl['bkg_id'] == comb_id)[0]

        outfile2d = outputfiles.spec_output_file(fitstbl, par,
                                            frames[0], slitname=slitname, twod=True)
        if not outfile2d.is_file() or overwrite:

            # Build history to document what contributed to the reduced
            # exposure
            history = History(fitstbl.frame_paths(frames[0]))
            history.add_reduce(calib_ID, fitstbl, frames, bg_frames)

            # TODO -- Should we reset/regenerate self.slits.mask for a new exposure

            # Reduce the exposure (standard and NIRSpec paths are both
            # handled by reduce_exposure via the slitname parameter)
            this_spec2d, this_sobjs = exposure.reduce_exposure(
                spectrograph, fitstbl, par, frames, calib_ID,
                calibrations_path, bg_frames=bg_frames,
                reuse_calibs=reuse_calibs, run_state=run_state,
                show=show, std_outfile=std_outfile,
                slitname=slitname, detnum=detnum)

            # TODO: come up with sensible naming convention for
            # save_exposure for combined files
            if this_spec2d is not None and len(this_spec2d.detectors) > 0:
                exposure.save_exposure(spectrograph,
                                    fitstbl, par, frames[0],
                                    this_spec2d, this_sobjs, calibrations_path,
                                    slitname=slitname,
                                    history=history,
                                    skip_write_2d=par['scienceframe']['process']['skip_write_2d'])
            else:
                log.warning(
                    'No spec2d and spec1d saved to file because the calibration/reduction was '
                    'not successful for all the detectors'
                )
        else:
            log.warning(
                f'Output file: {fitstbl.construct_basename(frames[0], slitname=slitname)} '
                'already exists. Set overwrite=True to recreate and overwrite.'
            )






