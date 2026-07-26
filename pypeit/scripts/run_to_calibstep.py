"""
Script to run to a single calibration step for an input frame

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase

class RunToCalibStep(scriptbase.ScriptBase):

    valid_steps = ['align', 'arc', 'bias', 'bpm', 'dark', 'flats', 'scattlight', 
                   'slits', 'tiltimg', 'tilts', 'wv_calib']


    @classmethod
    def get_parser(cls, width=None):
        parser = super().get_parser(description='Run PypeIt to a single calibration step for an input frame',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')
        parser.add_argument('step', type=str, help=f"Calibration step to perform.  Valid steps are: {', '.join(cls.valid_steps)}")

        parser.add_argument('--science_frame', type=str, help='Raw science frame to reduce as listed in your PypeIt file, e.g. b28.fits.gz. Either this or the calib_group must be provided')
        parser.add_argument('--calib_group', type=str, help='Calibration group ID to reduce. Either this or the frame must be provided')
        parser.add_argument('--det', type=str, help='Detector to reduce')

        # TODO -- Grab these from run_pypeit.py ?
        parser.add_argument('-r', '--redux_path', default=None,
                            help='Path to directory for the reduction.  Only advised for testing')
        parser.add_argument('-s', '--show', default=False, action='store_true',
                            help='Show reduction steps via plots (which will block further '
                                 'execution until clicked on) and outputs to ginga. Requires '
                                 'remote control ginga session via '
                                 '"ginga --modules=RC,SlitWavelength &"')

        return parser

    @classmethod
    def main(cls, args):

        import numpy as np
        from IPython import embed
        from pathlib import Path

        from pypeit import pypeit
        from pypeit import pypeit_steps
        from pypeit import log
        from pypeit import PypeItError
        from pypeit.core import parse

        # Set a default log file based on the name of the pypeit file, not the
        # name of the script
        if args.log_file == 'default':
            _pypeit_file = Path(args.pypeit_file)
            if _pypeit_file.suffix != '.pypeit':
                raise PypeItError('Input file must have a .pypeit extension!')
            args.log_file = _pypeit_file.with_suffix('.log')

        # Initialize the log
        cls.init_log(args)

        # Check for the frame or calib_group
        if args.science_frame is None and args.calib_group is None:
            raise PypeItError('Must provide either a science frame or a calibration group ID')
        elif args.science_frame is not None and args.calib_group is not None:
            log.warning("Both science_frame and calib_group ID provided.  Will use the science_frame")

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(
            args.pypeit_file, redux_path=args.redux_path, show=args.show, calib_only=True
        )
        pypeIt.reuse_calibs = True

        # Find the detectors to reduce
        if args.det is None:
            dets = pypeIt.par['rdx']['detnum']
        else:
            dets = parse.eval_detectors(args.det)
        # NOTE: dets *can be* None

        detectors = pypeIt.spectrograph.select_detectors(dets if pypeIt.par['rdx']['slitspatnum'] is None else dets)

        # Find the row of the frame.  Both branches leave ``row`` as an array
        # of matching indices, so the int(row[0]) below selects the first.
        if args.science_frame is not None:
            row = np.where(pypeIt.fitstbl['filename'] == args.science_frame)[0]
            if len(row) != 1:
                raise PypeItError(f"Frame {args.science_frame} not found or not unique")
        elif args.calib_group is not None:
            row = np.where((pypeIt.fitstbl['calib'].data.astype(str) == args.calib_group))[0]
            if len(row) == 0:
                raise PypeItError(f"Calibration group {args.calib_group} not found")
        row = int(row[0])
        calib_id = pypeIt.fitstbl.find_frame_calib_groups(row)[0]

        # Preserve any existing reduction state — especially the science
        # entries (a science step-build via pypeit_reduce_by_step, or the
        # Dashboard's planned frames).  PypeIt starts with a fresh run_state
        # (no science); without this, writing the state after a calibration
        # build would blank out the science portion of *_state.json.
        pypeIt.run_state.merge_from_disk()

        # Calibrations?
        # Pass run_state so the reduction state file (*_state.json) is written
        # as each step runs -- without it the state is never updated, so e.g.
        # the PypeIt Dashboard would show stale statuses after a (re)build
        # (and could not monitor the (re)build live).
        for det in detectors:
            pypeit_steps.calib_one(pypeIt.spectrograph, pypeIt.fitstbl, pypeIt.par, det, calib_id, pypeIt.calibrations_path, stop_at_step=args.step, run_state=pypeIt.run_state)

        # QA HTML
        log.info('Generating QA HTML')
        pypeIt.build_qa()

        # Refresh the *full* reduction state from disk and persist it.  The
        # stop-at-step run above only visited steps up to args.step, so steps
        # after it (and any in other calibration groups) would otherwise be left
        # at their initial status even when their files exist on disk.  A
        # status-only reload pass re-derives every step's status from the
        # Calibrations/ outputs (present -> success, missing -> pending); since
        # status_only does not auto-write per step, persist it explicitly.
        pypeIt.calib_all(status_only=True, reload_only=True)
        pypeIt.run_state.write()

        return 0

