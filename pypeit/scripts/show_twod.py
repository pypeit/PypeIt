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
from pypeit import traceslits

def parser(options=None):

    parser = argparse.ArgumentParser(description='Display sky subtracted, spec2d image in a Ginga viewer.  Run above the Science/ folder',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default = None, help = 'PYPIT spec2d file')
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('--resid', default=False, help="Set to show model residuals map: chi = (image - sky - obj)/noise",
                        action = "store_true",)
    parser.add_argument('--sky_resid', default=False, help="Set to show sky subtraction residuals map : chi = (image - sky)/noise",
                        action = "store_true")
    parser.add_argument('--showmask', default=False, help="Overplot masked pixels", action = "store_true")
    parser.add_argument('--embed', default=False, help="Upong completetion embed in ipython shell", action = "store_true")


    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):

    import os
    import numpy as np
    import pdb as debugger
    import IPython

    from astropy.io import fits
    from astropy.table import Table
    from astropy.stats import sigma_clipped_stats

    from pypeit import msgs
    from pypeit import ginga
    from pypeit.core import masters
    from pypeit.core.parse import get_dnum
    from pypeit.core import trace_slits
    from pypeit import traceslits

    # Check that arguments make sense
    if (args.resid is True & args.sky_resid is True):
        msgs.error('You cannot set both --resid and --sky_resid.\n '
                   'EITHER you show the model residual map (--resid) OR sky subtraction residuals (--sky_resid)')

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

    # Determine which images to read in and show
    if(args.resid is False and args.sky_resid is False):
        image = (sciimg - skymodel)*(bitmask == 0)  # sky subtracted image
        (mean, med, sigma) = sigma_clipped_stats(image[bitmask ==0],sigma_lower=5.0,sigma_upper=5.0)
        cut_min = mean -1.0*sigma
        cut_max = mean +4.0*sigma
        msgs.info('Showing sky subtracted image')
    else:
        # Cut levels for chi are always set to (-5,5)
        cut_min = -5.0
        cut_max = 5.0
        # Read in the noise model for residual maps
        try:
            exten = names.index('DET{:s}-IVARMODEL'.format(sdet))
        except ValueError:  # Backwards compatability
            msgs.error("Requested detector {:s} has no IVARMODEL.\n Maybe you chose the wrong one to view?\n" +
                       "Set with --det= or check file contents with --list".format(sdet))
        ivarmodel = hdu[exten].data
        if(args.sky_resid is True):
            image = (sciimg - skymodel)*np.sqrt(ivarmodel)*(bitmask == 0) # sky subtraction residual map
            msgs.info('Showing sky subtracted residual map')
        elif(args.resid is True):
            msgs.info('Showing full model residual map')
            # Read in the object model for residual map
            try:
                exten = names.index('DET{:s}-OBJ'.format(sdet))
            except ValueError:  # Backwards compatability
                msgs.error("Requested detector {:s} has no object model.\n Maybe you chose the wrong one to view?\n" +
                           "Set with --det= or check file contents with --list".format(sdet))
            objmodel = hdu[exten].data
            image = (sciimg - skymodel - objmodel) * np.sqrt(ivarmodel)*(bitmask == 0) # full model residual map


    # Show Image
    cwd = os.getcwd()
    wcs_img = cwd+'/'+ os.path.basename(os.path.normpath(head0['PYPMFDIR'])) +\
              '/MasterWave_'+'{:s}_{:02d}_{:s}.fits'.format(head0['PYPCNFIG'], args.det, head0['PYPCALIB'])


    if args.showmask:
        bitmask_in = bitmask
    else:
        bitmask_in = None
    viewer, ch = ginga.show_image(image, chname='DET{:s}'.format(sdet), wcs_img=wcs_img,
                                  cuts = (cut_min, cut_max), bitmask = bitmask_in)

    mdir = head0['PYPMFDIR']+'/'
    setup = '{:s}_{:s}_{:s}'.format(head0['PYPCNFIG'], sdet, head0['PYPCALIB'])

    # Load Tslits
    trc_file = masters.master_name('trace', setup, mdir)
    Tslits = traceslits.TraceSlits.from_master_files(trc_file)

    # Get slit ids
    #stup = (Tslits.mstrace.shape, Tslits.lcen, Tslits.rcen)
    slit_ids = [trace_slits.get_slitid(Tslits.mstrace.shape, Tslits.lcen, Tslits.rcen, ii)[0] for ii in range(Tslits.lcen.shape[1])]
    ginga.show_slits(viewer, ch, Tslits.lcen, Tslits.rcen, slit_ids)#, args.det)

    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    hdulist_1d = fits.open(spec1d_file)
    det_nm = 'DET{:s}'.format(sdet)
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['TRACE']
            obj_id = hdu.name.split('-')[0]
            ginga.show_trace(viewer, ch, trace, obj_id, color='orange') #hdu.name)


    if args.embed:
        IPython.embed()

