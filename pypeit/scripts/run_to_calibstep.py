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
        import argparse

        parser = super().get_parser(description='Run PypeIt to a single calibration step for an input frame',
                                    width=width, formatter=scriptbase.SmartFormatter)
        parser.add_argument('pypeit_file', type=str,
                            help='PypeIt reduction file (must have .pypeit extension)')
        parser.add_argument('step', type=str, help=f"Calibration step to perform.  Valid steps are: {', '.join(cls.valid_steps)}")

        parser.add_argument('--science_frame', type=str, help='Raw science frame to reduce as listed in your PypeIt file, e.g. b28.fits.gz. Either this or the calib_group must be provided')
        parser.add_argument('--calib_group', type=str, help='Calibration group ID to reduce. Either this or the frame must be provided')
        parser.add_argument('--det', type=str, help='Detector to reduce')

        # TODO -- Grab these from run_pypeit.py ?
        parser.add_argument('-v', '--verbosity', type=int, default=2,
                            help='Verbosity level between 0 [none] and 2 [all]')

        parser.add_argument('-r', '--redux_path', default=None,
                            help='Path to directory for the reduction.  Only advised for testing')
        parser.add_argument('-s', '--show', default=False, action='store_true',
                            help='Show reduction steps via plots (which will block further '
                                 'execution until clicked on) and outputs to ginga. Requires '
                                 'remote control ginga session via '
                                 '"ginga --modules=RC,SlitWavelength &"')

        return parser

    @staticmethod
    def main(args):

        import ast
        import numpy as np
        from IPython import embed
        from pathlib import Path

        from pypeit import pypeit
        from pypeit import msgs

        # Load options from command line
        _pypeit_file = Path(args.pypeit_file).absolute()
        if _pypeit_file.suffix != '.pypeit':
            msgs.error(f'Input file {_pypeit_file} must have a .pypeit extension!')
        logname = _pypeit_file.parent / f'{_pypeit_file.stem}.log'

        # Check for the frame or calib_group
        if args.science_frame is None and args.calib_group is None:
            msgs.error('Must provide either a science frame or a calibration group ID')
        elif args.science_frame is not None and args.calib_group is not None:
            msgs.warn("Both science_frame and calib_group ID provided.  Will use the science_frame")

        # Instantiate the main pipeline reduction object
        pypeIt = pypeit.PypeIt(args.pypeit_file, verbosity=args.verbosity,
                               redux_path=args.redux_path, 
                               logname=logname, show=args.show, calib_only=True)
        pypeIt.reuse_calibs = True

        # Find the detectors to reduce
        dets = pypeIt.par['rdx']['detnum'] if args.det is None else ast.literal_eval(args.det)
        detectors = pypeIt.select_detectors(
            pypeIt.spectrograph, dets,
            slitspatnum=pypeIt.par['rdx']['slitspatnum'])

        # Find the row of the frame
        if args.science_frame is not None:
            row = np.where(pypeIt.fitstbl['filename'] == args.science_frame)[0]
            if len(row) != 1:
                msgs.error(f"Frame {args.frame} not found or not unique")
        elif args.calib_group is not None:
            rows = np.where((pypeIt.fitstbl['calib'].data.astype(str) == args.calib_group))[0] 
            if len(rows) == 0:
                msgs.error(f"Calibration group {args.calib_group} not found")
            row = rows[0]
        row = int(row)

        # Calibrations?
        for det in detectors:
            pypeIt.calib_one([row], det, stop_at_step=args.step)
        
        # QA HTML
        msgs.info('Generating QA HTML')
        pypeIt.build_qa()
        msgs.close()

        return 0

