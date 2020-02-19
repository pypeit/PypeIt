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
import pdb
import argparse
import numpy as np
from configobj import ConfigObj

from astropy.table import Table
from astropy.io import fits

from pypeit.core.gui import skysub_regions as gui_skysub_regions
from pypeit import msgs
from pypeit.core.parse import get_dnum
from pypeit import edgetrace
from pypeit.par.util import parse_pypeit_file
from pypeit.masterframe import MasterFrame
from pypeit.spectrographs.util import load_spectrograph


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display a Raw science image and interactively define'
                                                 'the sky regions using a GUI. Run in the same folder'
                                                 'as your .pypeit file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PypeIt file')
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


def get_science_frame(usrdata):
    """Find all of the indices that correspond to science frames
    """
    sciidx = np.array([], dtype=np.int)
    cntr = 0
    for tt in range(len(usrdata)):
        ftype = usrdata['frametype'][tt].split(",")
        for ff in range(len(ftype)):
            if ftype[ff] == "science":
                sciidx = np.append(sciidx, tt)
                print("  ({0:d}) {1:s}".format(cntr+1, usrdata['filename'][tt]))
                cntr += 1
                break
    # Determine which science frame the user wants
    if cntr == 1:
        msgs.info("Only one science frame listed in .pypeit file - using that frame")
        idx = sciidx[0]
    else:
        ans = ''
        while True:
            ans = input(" Which frame would you like to select [1-{0:d}]: ".format(cntr))
            try:
                ans = int(ans)
                if 1 <= ans <= cntr:
                    idx = sciidx[ans-1]
                    break
            except ValueError:
                pass
            msgs.info("That is not a valid option!")
    return idx


def main_blah(args):

    # Load fits file
    cfg_lines, data_files, frametype, usrdata, setups = parse_pypeit_file(args.file, runtime=False)

    # Get the raw fits file name
    sciIdx = get_science_frame(usrdata)
    fname = data_files[sciIdx]
    setup = setups[0]

    # Load the spectrograph
    cfg = ConfigObj(cfg_lines)
    spectrograph_name = cfg['rdx']['spectrograph']
    spectrograph = load_spectrograph(spectrograph_name, ifile=data_files[0])
    msgs.info('Loaded spectrograph {0}'.format(spectrograph.spectrograph))

    pdb.set_trace()
    frame, _, _, datasec, _ = spectrograph.get_rawimage(fname, args.det)

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

    # Derive an appropriate output filename
    prefix = None
    outname = "{0:s}_skyregions.fits".format(prefix)

    # Finally, initialise the GUI
    skyreg = gui_skysub_regions.initialize(args.det, frame, tslits_dict, outname=outname, runtime=False, printout=True)

    # Get the results
    skyreg.get_results()


def main(args):

    # Load fits file
    hdu = fits.open(args.file)
    head0 = hdu[0].header

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
    frame = (sciimg) * (mask == 0)

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

    # Derive an appropriate output filename
    prefix = "blah"
    outname = "{0:s}_skyregions.fits".format(prefix)

    # Finally, initialise the GUI
    skyreg = gui_skysub_regions.initialize(args.det, frame, tslits_dict, outname=outname, runtime=False, printout=True)

    # Get the results
    skyreg.get_result()

