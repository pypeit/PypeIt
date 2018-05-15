#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
from pypit import traceslits

def parser(options=None):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in a Ginga viewer.  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'PYPIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):

    import os

    import pdb as debugger

    from astropy.io import fits
    from astropy.table import Table

    import pypit.ginga
    from pypit import pyputils
    from pypit import armasters
    from pypit.arparse import get_dnum
    from pypit.arspecobj import get_slitid

    # List only?
    hdu = fits.open(args.file)
    head0 = hdu[0].header
    if args.list:
        print(hdu.info())
        return

    # Setup for PYPIT imports
    msgs = pyputils.get_dummy_logger()

    # Init
    sdet = get_dnum(args.det, prefix=False)

    # One detector, sky sub for now
    names = [hdu[i].name for i in range(len(hdu))]
    try:
        exten = names.index('DET{:s}-SKYSUB'.format(sdet))
    except ValueError:  # Backwards compatability
        try:
            exten = names.index('DET{:d}-SKYSUB'.format(args.det))
        except ValueError:
            raise IOError("Requested detector {:s} was not processed.\n Maybe you chose the wrong one to view?\n  Set with --det=".format(sdet))
    skysub = hdu[exten].data

    # Show Image
    cwd = os.getcwd()
    wcs_img = cwd+'/'+head0['PYPMFDIR']+'/MasterWave_'+'{:s}_{:02d}_{:s}.fits'.format(head0['PYPCNFIG'], args.det, head0['PYPCALIB'])
    viewer, ch = pypit.ginga.show_image(skysub, chname='DET{:s}'.format(sdet), wcs_img=wcs_img)

    # Add slits
    testing = False
    if testing:
        mdir = 'MF_keck_lris_blue/'
        setup = 'A_{:s}_aa'.format(sdet)
    else:
        mdir = head0['PYPMFDIR']+'/'
        setup = '{:s}_{:s}_{:s}'.format(head0['PYPCNFIG'], sdet, head0['PYPCALIB'])

    # Load Tslits
    trc_file = armasters.master_name('trace', setup, mdir=mdir)
    Tslits = traceslits.TraceSlits.from_master_files(trc_file)
    lordloc = Tslits.lcen
    rordloc = Tslits.rcen

    # Get slit ids
    stup = (Tslits.mstrace.shape, lordloc, rordloc)
    slit_ids = [get_slitid(stup, None, ii)[0] for ii in range(lordloc.shape[1])]
    pypit.ginga.show_slits(viewer, ch, lordloc, rordloc, slit_ids)#, args.det)

    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    hdulist_1d = fits.open(spec1d_file)
    det_nm = 'D{:s}'.format(sdet)
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['obj_trace']
            obj_id = hdu.name.split('-')[0]
            pypit.ginga.show_trace(viewer, ch, trace, obj_id, color='green') #hdu.name)

