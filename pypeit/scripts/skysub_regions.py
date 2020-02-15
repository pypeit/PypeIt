#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script enables the user to view a 2D FITS file
and define the sky background regions interactively.
Run above the Science/ folder.
"""

import os
import argparse
import numpy as np

from astropy.table import Table
from astropy.io import fits

from pypeit.core.gui import skysub_regions as gui_skysub_regions
from pypeit import msgs
from pypeit.core.parse import get_dnum
from pypeit import edgetrace
from pypeit.masterframe import MasterFrame


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display a spec2d image to interactively define the'
                                                 'sky regions GUI.  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PYPEIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?",
                        action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")

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
    trc_file = trc_file.replace('Trace', 'Edges')+'.gz'
    tslits_dict = edgetrace.EdgeTraceSet.from_file(trc_file).convert_to_tslits_dict()
    shape = (tslits_dict['nspec'], tslits_dict['nspat'])

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
    gui_skysub_regions.initialise(args.det, frame, tslits_dict, None, printout=True)

"""

def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser()
    parser = argparse.ArgumentParser(description='Launch PypeIt identify tool, display extracted'
                                                 'MasterArc, and load linelist.'
                                                 'Run above the Masters/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type=str, default=None, help='PypeIt MasterArc file')
    parser.add_argument("--lamps", default=None, help="Comma separated list of calibration lamps (no spaces)",
                        action="store_true")
    parser.add_argument("--wmin", default=3000.0, help="Minimum wavelength range",
                        action="store_true")
    parser.add_argument("--wmax", default=10000.0, help="Maximum wavelength range",
                        action="store_true")
    parser.add_argument("--slit", default=0, help="Slit number to wavelength calibrate",
                        action="store_true")
    parser.add_argument("--det", default=1, help="Detector index",
                        action="store_true")

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import pdb
    import os
    import sys
    import astropy.io.fits as fits
    from pypeit.masterframe import MasterFrame
    from pypeit.spectrographs.util import load_spectrograph
    from pypeit.core import parse
    from pypeit.core.gui import identify as gui_identify
    from pypeit.core.wavecal import waveio, templates
    from pypeit.wavecalib import WaveCalib
    from pypeit import arcimage, edgetrace
    from pypeit.images import pypeitimage

    from IPython import embed

    # Load the MasterArc file
    if os.path.exists(args.file):
        arcfil = args.file
    else:
        try:
            arcfil = "Masters/{0:s}".format(args.file)
        except FileNotFoundError:
            print("Could not find MasterArc file.")
            sys.exit()
    msarc = pypeitimage.PypeItImage.from_file(arcfil)

    mdir = msarc.head0['MSTRDIR']
    mkey = msarc.head0['MSTRKEY']

    # Load the spectrograph
    specname = msarc.head0['PYP_SPEC']
    spec = load_spectrograph(specname)
    par = spec.default_pypeit_par()['calibrations']['wavelengths']

    # Get the lamp list
    if args.lamps is None:
        lamplist = par['lamps']
        if lamplist is None:
            print("ERROR :: Cannot determine the lamps")
            sys.exit()
    else:
        lamplist = args.lamps.split(",")
    par['lamps'] = lamplist

    # Load the tslits_dict
    trc_file = os.path.join(mdir, MasterFrame.construct_file_name('Edges', mkey, file_format='fits.gz'))
    tslits_dict = edgetrace.EdgeTraceSet.from_file(trc_file).convert_to_tslits_dict()

    # Check if a solution exists
    solnname = os.path.join(mdir, MasterFrame.construct_file_name('WaveCalib', mkey, file_format='json'))
    wv_calib = waveio.load_wavelength_calibration(solnname) if os.path.exists(solnname) else None

    # Load the MasterFrame (if it exists and is desired)?
    wavecal = WaveCalib(msarc, tslits_dict, spec, par)
    arccen, arc_maskslit = wavecal.extract_arcs()

    binspec, binspat = parse.parse_binning(spec.get_meta_value(msarc.head0['F1'], 'binning'))

    # Launch the identify window
    arcfitter = gui_identify.initialise(arccen, slit=int(args.slit), par=par, wv_calib_all=wv_calib)
    final_fit = arcfitter.get_results()

    # Ask the user if they wish to store the result in PypeIt calibrations
    ans = ''
    while ans != 'y' and ans != 'n':
        ans = input("Would you like to store this wavelength solution in the archive? (y/n):")
    if ans == 'y' and final_fit['rms'] < 0.1:
        gratname = fits.getheader(msarc.head0['F1'])[spec.meta['dispname']['card']].replace("/", "_")
        dispangl = "UNKNOWN"
        templates.pypeit_identify_record(final_fit, binspec, specname, gratname, dispangl)
        print("Your wavelength solution has been stored")
        print("Please consider sending your solution to the PypeIt team!")
"""