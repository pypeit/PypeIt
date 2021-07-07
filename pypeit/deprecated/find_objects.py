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

from pypeit import msgs
from pypeit import io
from pypeit import slittrace
from pypeit.core import gui
from pypeit.core.parse import get_dnum


def parse_args(options=None, return_parser=False):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in the'
                                                 'interactive object finding GUI.  Run above'
                                                 'the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PYPEIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?",
                        action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument("--old", default=False, action="store_true", help="Used old slit tracing")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def parse_traces(hdulist_1d, det_nm):
    """Extract the relevant trace information
    """
    traces = dict(traces=[], fwhm=[])
    pkflux = []
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['TRACE_SPAT']
            fwhm = tbl['FWHMFIT']
            obj_id = hdu.name.split('-')[0]
            traces['traces'].append(trace.copy())
            traces['fwhm'].append(np.median(fwhm))
            pkflux.append(np.median(tbl['BOX_COUNTS']))
    traces['pkflux'] = np.array(pkflux)
    return traces


def main(args):

    raise NotImplementedError('This script is currently out of date.')

    # List only?
    hdu = io.fits_open(args.file)
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
    mkey = head0['FRAMMKEY']
    mast_key = '{0}_{1:02d}'.format(mkey, args.det)
    if not os.path.exists(mdir):
        mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
        msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
        mdir = mdir_base

    # Assumes a MasterSlit file has been written
    #slits = slittrace.SlitTraceSet.from_master('{0}_{1:02d}'.format(head0['TRACMKEY'], args.det),
    #                                           mdir)
    # Load the slits information
    slits = slittrace.SlitTraceSet.from_master(mast_key, mdir)

    # Object traces
    left, right, mask = slits.select_edges()
    msgs.error("You need to choose which slits you care about here")

    # Get object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    if os.path.isfile(spec1d_file):
        hdulist_1d = io.fits_open(spec1d_file)
    else:
        hdulist_1d = []
        msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '                          No objects were extracted.')

    msgs.error("This code needs to be refactored since tslits_dict was removed...")
    import pdb
    pdb.set_trace()
    tslits_dict['objtrc'] = parse_traces(hdulist_1d, det_nm)
    obj_trace = parse_traces(hdulist_1d, 'DET{:s}'.format(sdet))

    # TODO :: Need to include standard star trace in the spec2d files
    std_trace = None

    # Extract some trace models
    fwhm = 2  # Start with some default value
    # TODO: Dictionaries like this are a pet peeve of mine.  I'd prefer
    # either individual objects or a class with a well-formed data model.
    # TODO: Why do all of these dictionary elements need fwhm?  Can they
    # be different?
    trace_models = dict()
    # Brightest object on slit
    trace_models['object'] = dict(trace_model=None, fwhm=fwhm)
    if len(obj_trace['pkflux']) > 0:
        smash_peakflux = obj_trace['pkflux']
        ibri = smash_peakflux.argmax()
        trace_models['object']['trace_model'] = obj_trace['traces'][ibri]
        trace_models['object']['fwhm'] = obj_trace['fwhm'][ibri]
    # Standard star trace
    trace_models['std'] = dict(trace_model=std_trace, fwhm=trace_models['object']['fwhm'])
    # Trace of the slit edge
    # TODO: Any particular reason to use the lefts?
    trace_models['slit'] = dict(trace_model=left.copy(), fwhm=trace_models['object']['fwhm'])

    # Finally, initialise the GUI
    gui.object_find.initialise(args.det, frame, left, right, obj_trace, trace_models, None,
                               printout=True, slit_ids=slits.id)

    ofgui = gui_object_find.initialise(args.det, frame, tslits_dict, None, printout=True, slit_ids=slits.id)
