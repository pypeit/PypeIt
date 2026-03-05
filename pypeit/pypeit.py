"""
Main driver class for PypeIt run

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
from pathlib import Path
import copy
import time
import os
import datetime

from astropy.io import fits
from astropy.table import Table
from IPython import embed

import numpy as np

from pypeit import inputfiles
from pypeit.core import qa
from pypeit import log
from pypeit import calibrations
from pypeit import utils
from pypeit.history import History
from pypeit.metadata import PypeItMetaData
from pypeit import outputfiles
from pypeit import exposure
from pypeit import pypeit_steps
from pypeit import PypeItError


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
        calib_only=False):

        # Set up logging
        self.pypeit_file = pypeit_file

        # State
        #self.run_state = state.RunPypeItState(pypeit_file=pypeit_file,
        #                                      current_step='init',
        #                                      current_det=-1,
        #                                      current_calibID=-1)
        #self.run_state = self.run_state.load()
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

        # Frame indices
        for calib_ID in self.fitstbl.calib_groups:
            # Find all the frames in this calibration group
            in_grp = self.fitstbl.find_calib_group(calib_ID)
            if not any(in_grp):
                continue
            # Find the detectors to reduce
            detectors = self.spectrograph.select_detectors(subset=self.par['rdx']['detnum'] if self.par['rdx']['slitspatnum'] is None
                                              else self.par['rdx']['slitspatnum'])
            log.info(f'Detectors to work on: {detectors}')

            # Loop on Detectors
            for self.det in detectors:
                log.info(f'Working on detector {self.det}')

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

    This class overrides :func:`reduce_all` to use the NIRSpec-specific
    ``reduce_calibID_nirspec`` function, and skips the standard calibration
    checks that don't apply to NIRSpec.

    See :class:`PypeIt` for arguments.
    """

    def _check_calibrations(self):
        """NIRSpec uses a different calibration workflow -- skip standard checks."""
        pass

    def reduce_all(self):
        """
        NIRSpec-specific reduction driver.

        Follows the same pattern as :func:`PypeIt.reduce_all` but calls
        :func:`reduce_calibID_nirspec` instead of :func:`reduce_calibID`.
        """
        self.par.validate_keys(required=['rdx', 'calibrations', 'scienceframe', 'reduce',
                                         'flexure'])
        self.tstart = time.perf_counter()

        # ############################################################################
        # Standard Star(s) Loop
        # ############################################################################
        for calib_ID in self.fitstbl.calib_groups:
            reduce_calibID_nirspec(self.spectrograph, self.par, self.fitstbl,
                                   calib_ID, self.calibrations_path,
                                   reduce_standard=True, overwrite=self.overwrite,
                                   show=self.show, reuse_calibs=self.reuse_calibs,
                                   qa_path=self.qa_path)

        # ############################################################################
        # Science Frame(s) Loop
        # ############################################################################
        for calib_ID in self.fitstbl.calib_groups:
            reduce_calibID_nirspec(self.spectrograph, self.par, self.fitstbl,
                                   calib_ID, self.calibrations_path,
                                   reduce_standard=False, overwrite=self.overwrite,
                                   show=self.show, reuse_calibs=self.reuse_calibs,
                                   qa_path=self.qa_path)
            log.info(f'Finished calibration group {calib_ID}')

        # Finish
        self.print_end_time()


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
    log.info(f'Found {len(grp_this)} {rtype} frames in calibration group {calib_ID}.')

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
#                    log.warning('PypeIt executed in quicklook mode.  Only reducing science frames '
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
                log.warning(
                    'No spec2d and spec1d saved to file because the calibration/reduction was '
                    'not successful for all the detectors'
                )
        else:
            log.warning(
                f'Output file: {fitstbl.construct_basename(frames[0])} already exists. Set '
                'overwrite=True to recreate and overwrite.'
            )


def reduce_calibID_nirspec(spectrograph, par, fitstbl, calib_ID, calibrations_path,
                           reduce_standard=False, overwrite=False, show=False,
                           reuse_calibs=True, qa_path=None):
    """
    Reduce all frames in a given calibration group for JWST NIRSpec,
    performing slit-by-slit calibration and extraction using JWST data models.

    This function parallels :func:`reduce_calibID` but handles the NIRSpec
    slit-by-slit workflow: loading JWST data models, iterating over individual
    slits, running :class:`~pypeit.calibrations.NIRSpecSlitCalibrations`,
    building science images via :class:`~pypeit.images.rawimage.NIRSpecRawImage`,
    and calling the standard find-objects and extraction machinery.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph object.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            The parameter set for the reduction.
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            The metadata table.
        calib_ID (:obj:`str`):
            The calibration group ID.
        calibrations_path (:obj:`str`):
            Path to the calibration files.
        reduce_standard (:obj:`bool`, optional):
            Reduce standard frames if True; science frames if False.
        overwrite (:obj:`bool`, optional):
            Overwrite existing files.
        show (:obj:`bool`, optional):
            Show plots interactively.
        reuse_calibs (:obj:`bool`, optional):
            Reuse pre-existing calibration files.
        qa_path (:obj:`str`, optional):
            Path for QA output.
    """
    try:
        from jwst import datamodels
    except ModuleNotFoundError:
        raise PypeItError('Unable to import jwst. Install pypeit with the jwst '
                          'option to reduce jwst data.')

    from pypeit import find_objects
    from pypeit import extraction
    from pypeit import spec2dobj
    from pypeit import specobjs
    from pypeit.images import combineimage
    from pypeit.images.rawimage import NIRSpecRawImage
    from pypeit.calibframe import CalibFrame

    if reduce_standard:
        is_this = fitstbl.find_frames('standard')
        rtype = 'standard'
    else:
        is_this = fitstbl.find_frames('science')
        rtype = 'science'

    frame_indx = np.arange(len(fitstbl))
    in_grp = fitstbl.find_calib_group(calib_ID)

    if not np.any(is_this & in_grp):
        return

    grp_this = frame_indx[is_this & in_grp]
    log.info(f'Found {len(grp_this)} {rtype} frames in calibration group {calib_ID}.')

    u_combid = np.unique(fitstbl['comb_id'][grp_this])

    for comb_id in u_combid:
        frames = np.where(fitstbl['comb_id'] == comb_id)[0]
        # bg_frames = np.where(
        #     (fitstbl['comb_id'] == fitstbl['bkg_id'][frames][0])
        #     & (fitstbl['comb_id'] >= 0))[0]
        # Find all frames whose comb_id matches the current frames bkg_id.
        # TODO TESTING!!!!
        bg_frames = \
        np.where((np.isin(fitstbl['comb_id'], [int(b) for b in fitstbl['bkg_id'][frames][0].split(',')]))
                 & (fitstbl['comb_id'] >= 0))[0]

        has_bg = len(bg_frames) > 0
        bkg_redux = has_bg
        find_negative = False
        if has_bg:
            find_negative = (
                ('science' in fitstbl['frametype'][bg_frames[0]])
                | ('standard' in fitstbl['frametype'][bg_frames[0]])
            ) if par['reduce']['findobj']['find_negative'] is None \
                else par['reduce']['findobj']['find_negative']

        # Load JWST data models
        sci_data = np.array(datamodels.open(fitstbl.frame_paths(frames)))
        sci_data_bkg = np.array(datamodels.open(
            fitstbl.frame_paths(bg_frames))) if bkg_redux else None

        # Load calibration and flat data
        calib_grps = fitstbl.find_frame_calib_groups(frames[0])
        cal_files = fitstbl.find_frame_files('trace', calib_ID=calib_grps)
        cal_data = np.array(datamodels.open(cal_files))
        flat_files = fitstbl.find_frame_files('pixelflat', calib_ID=calib_grps)
        flat_data = np.array(datamodels.open(flat_files))

        # Append FS slits if present
        for i, _flat_data in enumerate(flat_data):
            fs_path = flat_files[i].replace('_interpolatedflat',
                                            '_interpolatedflat_fs')
            if Path(fs_path).exists():
                log.info(f'Appending FS slits for {Path(flat_files[i]).name}')
                _flat_data_fs = datamodels.open(fs_path)
                for slit in _flat_data_fs.slits:
                    _flat_data.slits.append(slit)

        # Make a working copy of par so we can modify it for bkg_redux
        _par = copy.deepcopy(par)
        if bkg_redux:
            _par['reduce']['findobj']['skip_skysub'] = True
            _par['reduce']['extraction']['skip_optimal'] = True

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


        if _par['rdx']['maskIDs'] is not None:
            maskIDs = [str(mid).strip() for mid in _par['rdx']['maskIDs']]
            gd_slits_sources = [(slt, src_n, src_id, src_alias)
                                for slt, src_n, src_id, src_alias in slit_sources_uni
                                if str(slt).strip() in maskIDs or str(src_id) in maskIDs or str(src_alias) in maskIDs]
            if not gd_slits_sources:
                log.warning(f'No slits or sources found for maskIDs={_par["rdx"]["maskIDs"]}. '
                          'Skipping reduction.')
                return
            # find if there are maskIDs that are not present in the slit_sources_uni
            elif len(gd_slits_sources) < len(maskIDs):
                missing_maskIDs = [mid for mid in maskIDs if mid not in [slt for slt, _, _, _ in gd_slits_sources] and
                                  mid not in [src_id for _, _, src_id, _ in gd_slits_sources] and
                                  mid not in [src_alias for _, _, _, src_alias in gd_slits_sources]]
                log.warning(f'The following maskIDs were not found: {", ".join(missing_maskIDs)}. '
                          'Reduction will be performed on the available slits and sources.')

        else:
            gd_slits_sources = slit_sources_uni

        # MSGS info on which slits and sources are being reduced
        log.info(f'Reducing the following (slit_name, src_name): '
                  f'{", ".join([f"({slt}, {src_n})" for slt, src_n, _, _ in gd_slits_sources])}')

        # Loop over individual slits
        for ii, (islit, isource, isource_id, isource_alias) in enumerate(gd_slits_sources):
            # Print status message
            add_to_msgs = f'Slit name: {islit}' if isource is None else f'SRC name: {isource}'
            msgs_string = f'Reducing target {fitstbl["target"][frames[0]]} - {add_to_msgs}\n'

            msgs_string += 'Combining frames:\n'
            for iframe in frames:
                msgs_string += '{0:s}'.format(fitstbl['filename'][iframe]) + '\n'
            log.info(msgs_string)
            if has_bg:
                bg_msgs_string = ''
                for iframe in bg_frames:
                    bg_msgs_string += '{0:s}'.format(fitstbl['filename'][iframe]) + '\n'
                bg_msgs_string = '\nUsing background from frames:\n' + bg_msgs_string
                log.info(bg_msgs_string)

            # Determine detector for this slit
            if islit in slit_names_nrs1 and islit not in slit_names_nrs2:
                _det = 1
            elif islit in slit_names_nrs2 and islit not in slit_names_nrs1:
                _det = 2
            else:
                # TODO: implement mosaic reduction for slits spanning both detectors
                log.warning(f'Slit {islit} spans both NRS1 and NRS2. '
                            'Mosaic reduction not yet implemented. Skipping.')
                continue

            _cal_data = cal_data[[
                cd.meta.instrument.detector == f'NRS{_det}'
                for cd in cal_data]]
            _slit_slices = [
                spectrograph.get_slit_slice(_cal_data[0], islit)]
            _flat_data = flat_data[[
                fd.meta.instrument.detector == f'NRS{_det}'
                for fd in flat_data]]
            _sci_data = sci_data[[
                sd.meta.instrument.detector == f'NRS{_det}'
                for sd in sci_data]]
            if bkg_redux:
                _sci_data_bkg = sci_data_bkg[[
                    sd.meta.instrument.detector == f'NRS{_det}'
                    for sd in sci_data_bkg]]

            # Build basename and metadata
            objtype, setup, obstime, basename, binning \
                = _get_nirspec_metadata(fitstbl, spectrograph, frames[0], _det,
                                        slit_name=islit)

            log.info(f'Reducing detector {_det}')

            # --- Calibrations ---
            caliBrate = calibrations.Calibrations.get_instance(
                fitstbl, _par['calibrations'],
                spectrograph, calibrations_path,
                calib_ID, frames[0], _det,
                qadir=qa_path,
                reuse_calibs=reuse_calibs,
                show=show,
                user_slits=islit,
                jwst_cal_data=_cal_data,
                jwst_flat_data=_flat_data,
                slit_slices=_slit_slices)
            caliBrate.run_the_steps()

            if not caliBrate.success:
                log.warning(f'Calibrations failed for slit {islit}. Skipping.')
                continue

            # --- Build science image ---
            headarr = copy.deepcopy(
                spectrograph.get_headarr(
                    fits.open(fitstbl.frame_paths(frames)[0])))

            nirspec_raw = NIRSpecRawImage(
                spectrograph, _det, _sci_data, _slit_slices,
                caliBrate.flatimages, caliBrate.msbpm, headarr)
            sciImg = nirspec_raw.process(_par['scienceframe']['process'])

            # --- Background subtraction ---
            if bkg_redux:
                sciImg_bkg_list = []
                for _bg in _sci_data_bkg:
                    bkg_raw = NIRSpecRawImage(
                        spectrograph, _det, [_bg], _slit_slices,
                        caliBrate.flatimages, caliBrate.msbpm, headarr)
                    sciImg_bkg_list.append(
                        bkg_raw.process(_par['scienceframe']['process']))
                combImg = combineimage.CombineImage(
                    sciImg_bkg_list,
                    _par['scienceframe']['process'])
                comb_bkg = combImg.run(ignore_saturation=True)
                sciImg = sciImg.sub(comb_bkg)

            # --- Object finding ---
            std_redux = objtype == 'standard'
            objFind = find_objects.FindObjects.get_instance(
                sciImg, caliBrate.slits,
                spectrograph, _par, objtype,
                tilts=caliBrate.tilts,
                initial_skymask=None,
                bkg_redux=bkg_redux,
                find_negative=find_negative,
                std_redux=std_redux,
                show=show,
                basename=basename)

            initial_sky, sobjs_obj = objFind.run(
                std_trace=None, show_peaks=show)

            # --- Final global sky ---
            if (std_redux
                    or _par['reduce']['findobj']['skip_skysub']
                    or _par['reduce']['findobj']['skip_final_global']
                    or _par['reduce']['skysub']['user_regions'] is not None):
                final_global_sky = initial_sky
            else:
                skymask = objFind.create_skymask(sobjs_obj)
                final_global_sky = objFind.global_skysub(
                    previous_sky=initial_sky,
                    skymask=skymask, show=show,
                    reinit_bpm=False)

            scaleImg = objFind.scaleimg
            slits = copy.deepcopy(caliBrate.slits)
            slits.maskdef_designtab = None

            flagged_slits = np.where(objFind.reduce_bpm)[0]
            if len(flagged_slits) > 0:
                slits.mask[flagged_slits] = \
                    slits.bitmask.turn_on(
                        slits.mask[flagged_slits], 'BADSKYSUB')

            # --- Extraction ---
            if not _par['reduce']['extraction']['skip_extraction']:
                log.info(f'Extraction begins for {basename} on det={_det}')
                exTract = extraction.Extract.get_instance(
                    sciImg, slits, sobjs_obj, spectrograph,
                    _par, objtype,
                    global_sky=final_global_sky,
                    bkg_redux_global_sky=None,
                    tilts=caliBrate.tilts,
                    waveimg=caliBrate.waveimg,
                    flatimages=caliBrate.flatimages,
                    bkg_redux=bkg_redux,
                    return_negative=_par['reduce']['extraction']['return_negative'],
                    std_redux=std_redux,
                    basename=basename,
                    show=show)
                (skymodel, bkg_redux_skymodel, objmodel,
                 ivarmodel, outmask, sobjs, waveImg,
                 tilts_out, slits_out) = exTract.run()
                slitshift = exTract.slitshift
            else:
                log.info(f'Extraction skipped for {basename}')
                skymodel = final_global_sky
                bkg_redux_skymodel = None
                objmodel = np.zeros_like(objFind.sciImg.image)
                ivarmodel = np.copy(objFind.sciImg.ivar)
                outmask = objFind.sciImg.fullmask
                sobjs = sobjs_obj
                waveImg = objFind.waveimg
                tilts_out = objFind.tilts
                slits_out = slits
                slitshift = objFind.slitshift

            # --- Build outputs (same pattern as reduce_calibID) ---
            spec_flex_table = Table()
            spec_flex_table['spat_id'] = slits.spat_id
            spec_flex_table['sci_spec_flexure'] = slitshift

            spec2DObj = spec2dobj.Spec2DObj(
                sciimg=sciImg.image,
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
                vel_type=_par['calibrations']['wavelengths']['refframe'],
                tilts=tilts_out,
                slits=slits_out,
                wavesol=None,
                maskdef_designtab=None)
            spec2DObj.process_steps = sciImg.process_steps

            all_spec2d = spec2dobj.AllSpec2DObj()
            all_spec2d['meta']['bkg_redux'] = bkg_redux
            all_spec2d['meta']['find_negative'] = find_negative
            all_spec2d[sciImg.detector.name] = spec2DObj

            all_specobjs = specobjs.SpecObjs()
            if sobjs.nobj > 0:
                all_specobjs.add_sobj(sobjs)

            # --- Save (parallel to exposure.save_exposure) ---
            if len(all_spec2d.detectors) > 0:
                history = History(fitstbl.frame_paths(frames[0]))
                history.add_reduce(calib_ID, fitstbl, frames, bg_frames)
                _save_nirspec_exposure(
                    spectrograph, fitstbl, _par, calibrations_path,
                    frames[0], all_spec2d, all_specobjs,
                    basename, history)
            else:
                log.warning(f'No spec2d/spec1d saved for slit {islit}.')


def _get_nirspec_metadata(fitstbl, spectrograph, frame, det, slit_name=None):
    """
    Get metadata for a NIRSpec science frame.

    This is a module-level helper paralleling the metadata extraction
    done in :func:`reduce_calibID` / :func:`exposure.reduce_exposure`.

    Args:
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            Metadata table.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph instance.
        frame (:obj:`int`): Frame index.
        det (:obj:`int`): Detector number.
        slit_name (:obj:`str`, optional): Slit name.

    Returns:
        :obj:`tuple`: objtype, setup, obstime, basename, binning
    """
    from pypeit.calibframe import CalibFrame

    binning = fitstbl['binning'][frame]
    obstime = fitstbl.construct_obstime(frame)
    basename = fitstbl.construct_basename(frame, obstime=obstime)
    types = fitstbl['frametype'][frame].split(',')
    objtype = 'standard' if 'standard' in types else 'science'

    # Clean up NIRSpec-specific substrings from basename
    basename = basename.replace('_assign_wcs', '')
    basename = basename.replace('_nrs1', '')
    basename = basename.replace('_nrs2', '')

    setup = CalibFrame.construct_calib_key(
        fitstbl['setup'][frame],
        fitstbl['calib'][frame],
        spectrograph.get_det_name(det))

    if slit_name is not None:
        basename = basename.replace(
            spectrograph.camera,
            f'{slit_name}_{spectrograph.camera}')
        setup = setup + f'_{slit_name}'

    return objtype, setup, obstime, basename, binning


def _save_nirspec_exposure(spectrograph, fitstbl, par, calibrations_path,
                           frame, all_spec2d, all_specobjs, basename, history):
    """
    Save NIRSpec reduction outputs.

    This is a module-level helper paralleling :func:`exposure.save_exposure`.

    Args:
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph instance.
        fitstbl (:class:`~pypeit.metadata.PypeItMetaData`):
            Metadata table.
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Parameter set.
        calibrations_path (:obj:`str`):
            Path to calibrations.
        frame (:obj:`int`): Frame index.
        all_spec2d (:class:`~pypeit.spec2dobj.AllSpec2DObj`): 2D spectra.
        all_specobjs (:class:`~pypeit.specobjs.SpecObjs`): 1D spectra.
        basename (:obj:`str`): Output basename.
        history (:class:`~pypeit.history.History`): History object.
    """
    from pypeit import specobjs as specobjs_mod

    row_fitstbl = fitstbl[frame]
    rawfile = fitstbl.frame_paths(frame)
    head2d = fits.getheader(rawfile, ext=spectrograph.primary_hdrext)

    science_path = outputfiles.science_path(par)
    if not science_path.is_dir():
        science_path.mkdir(parents=True)

    subheader = spectrograph.subheader_for_spec(row_fitstbl, head2d)

    # 1D spectra
    if all_specobjs.nobj > 0 and not par['reduce']['extraction']['skip_extraction']:
        outfile1d = science_path / f'spec1d_{basename}.fits'
        all_specobjs.write_to_fits(subheader, outfile1d, history=history)
        outfiletxt = science_path / f'spec1d_{basename}.txt'
        sobjs = specobjs_mod.SpecObjs.from_fitsfile(outfile1d, chk_version=False)
        sobjs.write_info(outfiletxt, spectrograph.pypeline)

    if par['scienceframe']['process']['skip_write_2d']:
        return

    # 2D spectra
    outfile2d = science_path / f'spec2d_{basename}.fits'
    pri_hdr = all_spec2d.build_primary_hdr(
        head2d, spectrograph,
        redux_path=par['rdx']['redux_path'],
        calib_dir=calibrations_path,
        subheader=subheader, history=history)
    all_spec2d.write_to_fits(outfile2d, pri_hdr=pri_hdr)

