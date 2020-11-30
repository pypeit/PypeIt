#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
Trace slit edges for a set of images.
"""

def parse_args(options=None, return_parser=False):

    import argparse
    from pypeit.spectrographs import available_spectrographs

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Require either a pypeit file or a fits file 
    inp = parser.add_mutually_exclusive_group(required=True)
    inp.add_argument('-f', '--pypeit_file', default=None, type=str, help='PypeIt reduction file')
    inp.add_argument('-t', '--trace_file', default=None, type=str, help='Image to trace')

    parser.add_argument('-g', '--group', default=None,
                        help='If providing a pypeit file, use the trace images for this '
                             'calibration group.  If None, use the first calibration group.')
    parser.add_argument('-d', '--detector', default=None, type=int,
                        help='Only analyze the specified detector; otherwise analyze all or '
                             'detectors selected by the pypeit file, if provided.')
    parser.add_argument('-s', '--spectrograph', default=None, type=str,
                        help='A valid spectrograph identifier, which is only used if providing '
                             'files directly: {0}'.format(', '.join(available_spectrographs)))
    parser.add_argument('-b', '--binning', default=None, type=str,
                        help='Image binning in spectral and spatial directions.  Only used if '
                             'providing files directly; default is 1,1.')
    parser.add_argument('-p', '--redux_path', default=None,
                        help='Path to top-level output directory.  '
                             'Default is the current working directory.')
    parser.add_argument('-m', '--master_dir', default='Masters',
                        help='Name for directory in output path for Master file(s) relative '
                             'to the top-level directory.')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')

    parser.add_argument('--debug', default=False, action='store_true', help='Run in debug mode.')
    parser.add_argument('--show', default=False, action='store_true',
                        help='Show the stages of trace refinements (only for the new code).')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import time
    import os
    import numpy as np
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit import edgetrace
    from pypeit import slittrace
    from pypeit.pypeit import PypeIt
    from pypeit.images import buildimage
    from pypeit import masterframe

    from IPython import embed

    if args.pypeit_file is not None:
        pypeit_file = args.pypeit_file
        if not os.path.isfile(pypeit_file):
            raise FileNotFoundError('File does not exist: {0}'.format(pypeit_file))
        pypeit_file = os.path.abspath(pypeit_file)
        redux_path = os.path.abspath(os.path.split(pypeit_file)[0]
                                     if args.redux_path is None else args.redux_path)

        rdx = PypeIt(pypeit_file, redux_path=redux_path)
        detectors = rdx.par['rdx']['detnum'] if args.detector is None else args.detector
        # Save the spectrograph
        spec = rdx.spectrograph
        # Get the calibration group to use
        group = np.unique(rdx.fitstbl['calib'])[0] if args.group is None else args.group
        if group not in np.unique(rdx.fitstbl['calib']):
            raise ValueError('Not a valid calibration group: {0}'.format(group))
        # Find the rows in the metadata table with trace frames in the
        # specified calibration group
        tbl_rows = rdx.fitstbl.find_frames('trace', calib_ID=int(group), index=True)
        # Master keyword
        master_key_base = '_'.join(rdx.fitstbl.master_key(tbl_rows[0]).split('_')[:2])
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
        master_key_base = 'A_1'
        binning = '1,1' if args.binning is None else args.binning
        if not os.path.isfile(args.trace_file):
            raise FileNotFoundError('File does not exist: {0}'.format(args.trace_file))
        files = [os.path.abspath(args.trace_file)]
        redux_path = os.path.abspath(os.path.split(files[0])[0]
                                     if args.redux_path is None else args.redux_path)
        par = spec.default_pypeit_par()
        proc_par = par['calibrations']['traceframe']
        trace_par = par['calibrations']['slitedges']
        bias_files = None
        bias_par = None

        dark_files = None
        dark_par = None

        # Set the QA path
        qa_path = os.path.join(os.path.abspath(os.path.split(files[0])[0]), 'QA')

    if detectors is None: 
        detectors = np.arange(spec.ndet)+1
    if isinstance(detectors, int):
        detectors = [detectors]
    master_dir = os.path.join(redux_path, args.master_dir)
    for det in detectors:
        # Master keyword for output file name
        master_key = '{0}_{1}'.format(master_key_base, str(det).zfill(2))

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
            msdark = buildimage.buildimage_fromlist(spec, det, dark_par, dark_files) #, bias=msbias)

        msbpm = spec.bpm(files[0], det)

        # Build the trace image
        traceImage = buildimage.buildimage_fromlist(spec, det, proc_par, files, bias=msbias,
                                                    bpm=msbpm, dark=msdark)

        # Trace the slit edges
        t = time.perf_counter()
        edges = edgetrace.EdgeTraceSet(traceImage, spec, trace_par, bpm=msbpm, auto=True,
                                       debug=args.debug, show_stages=args.show,
                                       qa_path=qa_path)

        print('Tracing for detector {0} finished in {1} s.'.format(det, time.perf_counter()-t))
        # Write the MasterEdges file
        edge_masterframe_name = masterframe.construct_file_name(edgetrace.EdgeTraceSet,
                                                                master_key, master_dir=master_dir)
        edges.to_master_file(edge_masterframe_name)

        # Write the MasterSlits file
        slit_masterframe_name = masterframe.construct_file_name(slittrace.SlitTraceSet,
                                                                master_key, master_dir=master_dir)
        edges.get_slits().to_master_file(slit_masterframe_name)

    return 0



