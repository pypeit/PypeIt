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
from pypeit.masterframe import MasterFrame
from pypeit import slittrace


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
    mkey = head0['FRAMMKEY']
    mast_key = '{0}_{1:02d}'.format(mkey, args.det)
    if not os.path.exists(mdir):
        mdir_base = os.path.join(os.getcwd(), os.path.basename(mdir))
        msgs.warn('Master file dir: {0} does not exist. Using {1}'.format(mdir, mdir_base))
        mdir = mdir_base

    # Load the slits information
    slits = slittrace.SlitTraceSet.from_master(mast_key, mdir)

    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')

    det_nm = 'DET{:s}'.format(sdet)
    if os.path.isfile(spec1d_file):
        hdulist_1d = fits.open(spec1d_file)
    else:
        hdulist_1d = []
        msgs.warn('Could not find spec1d file: {:s}'.format(spec1d_file) + msgs.newline() +
                  '                          No objects were extracted.')

    msgs.error("This code needs to be refactored since tslits_dict was removed...")
    import pdb
    pdb.set_trace()
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
    gui_object_find.initialise(args.det, frame, tslits_dict, None, printout=True,
                               slit_ids=slits.id)
