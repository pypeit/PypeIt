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
from astropy.table import Table
from pypeit import traceslits
from pypeit import ginga


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in a Ginga viewer.  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'PYPIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")
#    parser.add_argument('--resid', default=False, help="Set to show model residuals map: chi = (image - sky - obj)/noise",
#                        action = "store_true",)
#    parser.add_argument('--sky_resid', default=False, help="Set to show sky subtraction residuals map : chi = (image - sky)/noise",
#                        action = "store_true")
    parser.add_argument('--showmask', default=False, help="Overplot masked pixels", action = "store_true")
    parser.add_argument('--embed', default=False, help="Upong completetion embed in ipython shell", action = "store_true")


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def show_trace(hdulist_1d, det_nm, viewer, ch):

    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['TRACE']
            obj_id = hdu.name.split('-')[0]
            ginga.show_trace(viewer, ch, trace, obj_id, color='orange') #hdu.name)



def main(args):

    import os
    import numpy as np
    import pdb as debugger
    import IPython

    from astropy.io import fits
    from astropy.stats import sigma_clipped_stats

    from pypeit import msgs
    from pypeit.core import masters
    from pypeit.core.parse import get_dnum
    from pypeit.core import trace_slits
    from pypeit import traceslits


    # List only?
    hdu = fits.open(args.file)
    head0 = hdu[0].header
    if args.list:
        hdu.info()
        return

    # Setup for PYPIT imports
    msgs.reset(verbosity=2)

    # Init
    sdet = get_dnum(args.det, prefix=False)

    # One detector, sky sub for now
    names = [hdu[i].name for i in range(len(hdu))]

    try:
        exten = names.index('DET{:s}-PROCESSED'.format(sdet))
    except:  # Backwards compatability
        msgs.error("Requested detector {:s} was not processed.\n Maybe you chose the wrong one to view? "
                   "Set with --det= or check file contents with --list".format(sdet))
    sciimg = hdu[exten].data
    try:
        exten = names.index('DET{:s}-SKY'.format(sdet))
    except:  # Backwards compatability
        msgs.error("Requested detector {:s} has no sky model.\n Maybe you chose the wrong one to view? "
                   "Set with --det= or check file contents with --list".format(sdet))
    skymodel = hdu[exten].data
    try:
        exten = names.index('DET{:s}-MASK'.format(sdet))
    except ValueError:  # Backwards compatability
        msgs.error("Requested detector {:s} has no bit mask.\n Maybe you chose the wrong one to view?\n" +
                   "Set with --det= or check file contents with --list".format(sdet))
    bitmask = hdu[exten].data
    try:
        exten = names.index('DET{:s}-IVARMODEL'.format(sdet))
    except ValueError:  # Backwards compatability
        msgs.error("Requested detector {:s} has no IVARMODEL.\n Maybe you chose the wrong one to view?\n" +
                   "Set with --det= or check file contents with --list".format(sdet))
    ivarmodel = hdu[exten].data
    # Read in the object model for residual map
    try:
        exten = names.index('DET{:s}-OBJ'.format(sdet))
    except ValueError:  # Backwards compatability
        msgs.error("Requested detector {:s} has no object model.\n Maybe you chose the wrong one to view?\n" +
                   "Set with --det= or check file contents with --list".format(sdet))
    objmodel = hdu[exten].data
    # Get waveimg
    cwd = os.getcwd()
    wcs_img = cwd+'/'+ os.path.basename(os.path.normpath(head0['PYPMFDIR'])) +\
              '/MasterWave_'+'{:s}_{:02d}_{:s}.fits'.format(head0['PYPCNFIG'], args.det, head0['PYPCALIB'])
    # Load Tslits
    mdir = head0['PYPMFDIR']+'/'
    setup = '{:s}_{:s}_{:s}'.format(head0['PYPCNFIG'], sdet, head0['PYPCALIB'])
    trc_file = masters.master_name('trace', setup, mdir)
    Tslits = traceslits.TraceSlits.from_master_files(trc_file)
    slit_ids = [trace_slits.get_slitid(Tslits.mstrace.shape, Tslits.lcen, Tslits.rcen, ii)[0] for ii in range(Tslits.lcen.shape[1])]
    # Show the bitmask?
    if args.showmask:
        bitmask_in = bitmask
    else:
        bitmask_in = None
    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    hdulist_1d = fits.open(spec1d_file)
    det_nm = 'DET{:s}'.format(sdet)

    # Now show each image to a separate channel

    # SKYSUB
    image = (sciimg - skymodel) * (bitmask == 0)  # sky subtracted image
    (mean, med, sigma) = sigma_clipped_stats(image[bitmask == 0], sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    viewer, ch = ginga.show_image(image, chname='skysub-det{:s}'.format(sdet).format(sdet), wcs_img=wcs_img,
                                  bitmask=bitmask_in) #, cuts=(cut_min, cut_max))
                                  # JFH For some reason Ginga crashes when I try to put cuts in here.
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, Tslits.lcen, Tslits.rcen, slit_ids)#, args.det)

    # SKRESIDS
    image = (sciimg - skymodel) * np.sqrt(ivarmodel) * (bitmask == 0)  # sky residual map
    viewer, ch = ginga.show_image(image, chname='sky_resid-det{:s}'.format(sdet).format(sdet), wcs_img=wcs_img,
                                  cuts=(-5.0, 5.0), bitmask = bitmask_in)
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, Tslits.lcen, Tslits.rcen, slit_ids)#, args.det)

    # RESIDS
    image = (sciimg - skymodel - objmodel) * np.sqrt(ivarmodel) * (bitmask == 0)  # full model residual map
    viewer, ch = ginga.show_image(image, chname='resid-det{:s}'.format(sdet).format(sdet), wcs_img=wcs_img,
                                  cuts = (-5.0, 5.0), bitmask = bitmask_in)
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, Tslits.lcen, Tslits.rcen, slit_ids)#, args.det)

    # After displaying all the images since up the images with WCS_MATCH



    if args.embed:
        IPython.embed()

