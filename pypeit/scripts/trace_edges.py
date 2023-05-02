"""
Trace slit edges for a set of images.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from pypeit.scripts import scriptbase


class TraceEdges(scriptbase.ScriptBase):

    @classmethod
    def get_parser(cls, width=None):
        from pypeit.spectrographs import available_spectrographs

        parser = super().get_parser(description='Trace slit edges', width=width)

        # Require either a pypeit file or a fits file
        inp = parser.add_mutually_exclusive_group(required=True)
        inp.add_argument('-f', '--pypeit_file', default=None, type=str,
                         help='PypeIt reduction file')
        inp.add_argument('-t', '--trace_file', default=None, type=str, help='Image to trace')

        parser.add_argument('-g', '--group', default=None,
                            help='If providing a pypeit file, use the trace images for this '
                                 'calibration group.  If None, use the first calibration group.')
        parser.add_argument('-d', '--detector', type=str, default=None, nargs='*',
                            help='Detector(s) to process.  If more than one, the list of detectors '
                                 'must be one of the allowed mosaics hard-coded for the selected '
                                 'spectrograph.  Using "mosaic" for gemini_gmos, keck_deimos, or '
                                 'keck_lris will use the default mosaic.')
        parser.add_argument('-s', '--spectrograph', default=None, type=str,
                            help='A valid spectrograph identifier, which is only used if providing'
                                 ' files directly: {0}'.format(', '.join(available_spectrographs)))
        parser.add_argument('-b', '--binning', default=None, type=str,
                            help='Image binning in spectral and spatial directions.  Only used if '
                                 'providing files directly; default is 1,1.')
        parser.add_argument('-p', '--redux_path', type=str, default='current working directory',
                            help='Path to top-level output directory.')
        parser.add_argument('-c', '--calib_dir', default='Calibrations',
                            help='Name for directory in output path for calibration file(s) '
                                 'relative to the top-level directory.')
        parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                            help='Overwrite any existing files/directories')

        parser.add_argument('--debug', default=False, action='store_true',
                            help='Run in debug mode.')
        parser.add_argument('--show', default=False, action='store_true',
                            help='Show the stages of trace refinements (only for the new code).')
        parser.add_argument('-v', '--verbosity', type=int, default=1,
                            help='Verbosity level between 0 [none] and 2 [all]. Default: 1. '
                                 'Level 2 writes a log with filename trace_edges_YYYYMMDD-HHMM.log')

        return parser

    @staticmethod
    def main(args):

        import time
        from pathlib import Path
        import numpy as np
        from pypeit import msgs
        from pypeit.spectrographs.util import load_spectrograph
        from pypeit import edgetrace
        from pypeit.pypeit import PypeIt
        from pypeit.images import buildimage

        from IPython import embed

        # Set the verbosity, and create a logfile if verbosity == 2
        msgs.set_logfile_and_verbosity('trace_edges', args.verbosity)

        if args.pypeit_file is not None:
            pypeit_file = Path(args.pypeit_file).resolve()
            if not pypeit_file.exists():
                msgs.error(f'File does not exist: {pypeit_file}')
            redux_path = pypeit_file.parent if args.redux_path is None \
                            else Path(args.redux_path).resolve()

            rdx = PypeIt(str(pypeit_file), redux_path=str(redux_path))
            detectors = rdx.par['rdx']['detnum'] if args.detector is None else args.detector
            # Save the spectrograph
            spec = rdx.spectrograph
            # Get the calibration group to use
            group = np.unique(rdx.fitstbl['calib'])[0] if args.group is None else args.group
            if group not in np.unique(rdx.fitstbl['calib']):
                msgs.error(f'Invalid calibration group: {group}')
            # Find the rows in the metadata table with trace frames in the
            # specified calibration group
            tbl_rows = rdx.fitstbl.find_frames('trace', calib_ID=int(group), index=True)
            setup = rdx.fitstbl['setup'][tbl_rows[0]]
            calib_id = rdx.fitstbl['calib'][tbl_rows[0]]
            # Save the binning
            binning = rdx.fitstbl['binning'][tbl_rows[0]]
            # Save the full file paths
            files = rdx.fitstbl.frame_paths(tbl_rows)
            # Trace image processing parameters
            proc_par = rdx.par['calibrations']['traceframe']
            # Slit tracing parameters
            trace_par = rdx.par['calibrations']['slitedges']

            # Get the bias files, if requested
            bias_rows = rdx.fitstbl.find_frames('bias', calib_ID=int(group), index=True)
            bias_files = rdx.fitstbl.frame_paths(bias_rows)
            bias_par = rdx.par['calibrations']['biasframe']
            if len(bias_files) == 0:
                bias_files = None

            # Get the dark files, if requested
            dark_rows = rdx.fitstbl.find_frames('dark', calib_ID=int(group), index=True)
            dark_files = rdx.fitstbl.frame_paths(dark_rows)
            dark_par = rdx.par['calibrations']['darkframe']
            if len(dark_files) == 0:
                dark_files = None

            # Set the QA path
            qa_path = rdx.qa_path
        else:
            detectors = args.detector
            spec = load_spectrograph(args.spectrograph)
            setup = 'A'     # Dummy value
            calib_id = '1'  # Dummy value
            binning = '1,1' if args.binning is None else args.binning
            trace_file = Path(args.trace_file).resolve()
            if not trace_file.exists():
                msgs.error(f'File does not exist: {trace_file}')
            files = [str(trace_file)]
            redux_path = trace_file.parent if args.redux_path is None \
                            else Path(args.redux_path).resolve()
            par = spec.default_pypeit_par()
            proc_par = par['calibrations']['traceframe']
            trace_par = par['calibrations']['slitedges']
            bias_files = None
            bias_par = None

            dark_files = None
            dark_par = None

            # Set the QA path
            qa_path = redux_path / 'QA'

        if detectors is None:
            detectors = np.arange(spec.ndet)+1
        elif isinstance(detectors, (int, tuple)):
            detectors = [detectors]
        elif any([isinstance(d,str) for d in detectors]):
            detectors = [eval(d) for d in detectors]

        calib_dir = redux_path / args.calib_dir
        for det in detectors:
            # Get the bias frame if requested
            if bias_files is None:
                proc_par['process']['use_biasimage'] = False
                msbias = None
            else:
                msbias = buildimage.buildimage_fromlist(spec, det, bias_par, bias_files)

            # Get the dark frame if requested
            if dark_files is None:
                proc_par['process']['use_darkimage'] = False
                msdark = None
            else:
                msdark = buildimage.buildimage_fromlist(spec, det, dark_par, dark_files)

            msbpm = spec.bpm(files[0], det)

            # Build the trace image
            traceImage = buildimage.buildimage_fromlist(spec, det, proc_par, files, bias=msbias,
                                                        bpm=msbpm, dark=msdark, setup=setup,
                                                        calib_id=calib_id, calib_dir=calib_dir)
            # Trace the slit edges
            t = time.perf_counter()
            edges = edgetrace.EdgeTraceSet(traceImage, spec, trace_par, auto=True,
                                           debug=args.debug, show_stages=args.show,
                                           qa_path=qa_path)

            print('Tracing for detector {0} finished in {1} s.'.format(det, time.perf_counter()-t))
            # Write the two calibration frames
            edges.to_file()
            edges.get_slits().to_file()

        return 0

