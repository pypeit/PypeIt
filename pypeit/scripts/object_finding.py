#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.
"""

import os
import argparse
import numpy as np

from astropy.table import Table
from astropy.io import fits

from pypeit.core.gui import object_find as gui_object_find
from pypeit import msgs
from pypeit.core.parse import get_dnum
from pypeit.traceslits import TraceSlits
from pypeit import edgetrace
from pypeit.masterframe import MasterFrame
from pypeit.core import trace_slits


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in the'
                                                 'interactive object finding GUI.  Run above'
                                                 'the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PYPEIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?",
                        action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument("--old", default=False, action="store_true", help="Used old slit tracing")

    return parser.parse_args() if options is None else parser.parse_args(options)


def parse_traces(hdulist_1d, det_nm):
    """Extract the relevant trace information
    """
    traces = dict(traces=[], fwhm=[])
    pkflux = []
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['TRACE']
            fwhm = tbl['FWHM']
            obj_id = hdu.name.split('-')[0]
            traces['traces'].append(trace.copy())
            traces['fwhm'].append(np.median(fwhm))
            pkflux.append(np.median(tbl['BOX_COUNTS']))
    traces['pkflux'] = np.array(pkflux)
    return traces


def main(args):

    # List only?
    hdu = fits.open(args.file)
    head0 = hdu[0].header
    if args.list:
        hdu.info()
        return

    # Init
    sdet = get_dnum(args.det, prefix=False)

    # One detector, sky sub for now
    names = [hdu[i].name for i in range(len(hdu))]

    try:
        exten = names.index('DET{:s}-PROCESSED'.format(sdet))
    except:  # Backwards compatability
        msgs.error('Requested detector {:s} was not processed.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))
    sciimg = hdu[exten].data
    try:
        exten = names.index('DET{:s}-SKY'.format(sdet))
    except:  # Backwards compatability
        msgs.error('Requested detector {:s} has no sky model.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))
    skymodel = hdu[exten].data
    try:
        exten = names.index('DET{:s}-MASK'.format(sdet))
    except ValueError:  # Backwards compatability
        msgs.error('Requested detector {:s} has no bit mask.\n'
                   'Maybe you chose the wrong one to view?\n'
                   'Set with --det= or check file contents with --list'.format(sdet))

    mask = hdu[exten].data
    frame = (sciimg - skymodel) * (mask == 0)

    mdir = head0['PYPMFDIR']
    if not os.path.exists(mdir):
        mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
        msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
        mdir = mdir_base

    trace_key = '{0}_{1:02d}'.format(head0['TRACMKEY'], args.det)
    trc_file = os.path.join(mdir, MasterFrame.construct_file_name('Trace', trace_key))

    # TODO -- Remove this once the move to Edges is complete
    if args.old:
        tslits_dict = TraceSlits.load_from_file(trc_file)[0]
    else:
        trc_file = trc_file.replace('Trace', 'Edges')+'.gz'
        tslits_dict = edgetrace.EdgeTraceSet.from_file(trc_file).convert_to_tslits_dict()
    shape = (tslits_dict['nspec'], tslits_dict['nspat'])
    slit_ids = [trace_slits.get_slitid(shape, tslits_dict['slit_left'], tslits_dict['slit_righ'], ii)[0]
                for ii in range(tslits_dict['slit_left'].shape[1])]

    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')

    det_nm = 'DET{:s}'.format(sdet)
    if os.path.isfile(spec1d_file):
        hdulist_1d = fits.open(spec1d_file)
    else:
        hdulist_1d = []
        msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '                          No objects were extracted.')
    tslits_dict['objtrc'] = parse_traces(hdulist_1d, det_nm)

    # TODO :: Need to include standard star trace in the spec2d files
    std_trace = None

    # Extract some trace models
    fwhm = 2  # Start with some default value
    # Brightest object on slit
    trace_model_obj = None
    trace_model_dict = dict()
    if len(tslits_dict['objtrc']['pkflux']) > 0:
        smash_peakflux = tslits_dict['objtrc']['pkflux']
        ibri = smash_peakflux.argmax()
        trace_model_obj = tslits_dict['objtrc']['traces'][ibri]
        fwhm = tslits_dict['objtrc']['fwhm'][ibri]
    trace_model_dict['object'] = dict(trace_model=trace_model_obj, fwhm=fwhm)
    # Standard star trace
    trace_model_dict['std'] = dict(trace_model=std_trace, fwhm=fwhm)
    # Trace of the slit edge
    trace_model_dict['slit'] = dict(trace_model=tslits_dict['slit_left'].copy(), fwhm=fwhm)
    tslits_dict['trace_model'] = trace_model_dict

    # Finally, initialise the GUI
    gui_object_find.initialise(args.det, frame, tslits_dict, None, printout=True, slit_ids=slit_ids)
