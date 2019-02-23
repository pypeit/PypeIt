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
from pypeit import ginga
from pypeit.spectrographs import util
from pypeit.processimages import ProcessImagesBitMask as bitmask
from pypeit.core import pixels

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
    from pypeit import masterframe
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
    mask = hdu[exten].data
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
    mdir = head0['PYPMFDIR']+'/'
    if not os.path.exists(mdir):
        mdir_base = os.path.basename(os.path.dirname(mdir)) + '/'
        msgs.warn('Master file dir: {:s}'.format(mdir) + ' does not exist. Using ./{:s}'.format(mdir_base))
        mdir=mdir_base

    wave_key = '{:s}'.format(head0['ARCMKEY']) +  '_{:02d}'.format(args.det)
    waveimg = masterframe.master_name('wave', wave_key, mdir)

    trace_key = '{:s}'.format(head0['TRACMKEY']) + '_{:02d}'.format(args.det)
    trc_file = masterframe.master_name('trace', trace_key, mdir)
    tslits_dict, _ = traceslits.load_tslits(trc_file)
    spectrograph = util.load_spectrograph(tslits_dict['spectrograph'])
    slitmask = pixels.tslits2mask(tslits_dict)
    shape = (tslits_dict['nspec'], tslits_dict['nspat'])
    slit_ids = [trace_slits.get_slitid(shape, tslits_dict['slit_left'], tslits_dict['slit_righ'], ii)[0]
                for ii in range(tslits_dict['slit_left'].shape[1])]
    # Show the bitmask?
    if args.showmask:
        mask_in = mask
    else:
        mask_in = None
    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    hdulist_1d = fits.open(spec1d_file)
    det_nm = 'DET{:s}'.format(sdet)

    # Unpack the bitmask
    bitMask = bitmask()
    bpm, crmask, satmask, minmask, offslitmask, nanmask, ivar0mask, ivarnanmask, extractmask = bitMask.unpack(mask)

    # Now show each image to a separate channel

    # SCIIMG

    image = sciimg # Raw science image
    (mean, med, sigma) = sigma_clipped_stats(image[mask == 0], sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    chname_skysub='sciimg-det{:s}'.format(sdet)
    # Clear all channels at the beginning
    viewer, ch = ginga.show_image(image, chname=chname_skysub, waveimg=waveimg,
                                  bitmask=mask_in, clear=True) #, cuts=(cut_min, cut_max), wcs_match=True)
    ginga.show_slits(viewer, ch, tslits_dict['slit_left'], tslits_dict['slit_righ'], slit_ids)#, args.det)

    # SKYSUB
    image = (sciimg - skymodel) * (mask == 0)  # sky subtracted image
    (mean, med, sigma) = sigma_clipped_stats(image[mask == 0], sigma_lower=5.0, sigma_upper=5.0)
    cut_min = mean - 1.0 * sigma
    cut_max = mean + 4.0 * sigma
    chname_skysub='skysub-det{:s}'.format(sdet)
    # Clear all channels at the beginning
    viewer, ch = ginga.show_image(image, chname=chname_skysub, waveimg=waveimg,
                                  bitmask=mask_in) #, cuts=(cut_min, cut_max),wcs_match=True)
                                  # JFH For some reason Ginga crashes when I try to put cuts in here.
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, tslits_dict['slit_left'], tslits_dict['slit_righ'], slit_ids)#, args.det)

    # SKRESIDS
    chname_skyresids = 'sky_resid-det{:s}'.format(sdet)
    image = (sciimg - skymodel) * np.sqrt(ivarmodel) * (mask == 0)  # sky residual map
    viewer, ch = ginga.show_image(image, chname_skyresids, waveimg=waveimg,
                                  cuts=(-5.0, 5.0), bitmask = mask_in) #,wcs_match=True)
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, tslits_dict['slit_left'], tslits_dict['slit_righ'], slit_ids)#, args.det)

    # RESIDS
    chname_resids = 'resid-det{:s}'.format(sdet)
    image = (sciimg - skymodel - objmodel) * np.sqrt(ivarmodel) * (mask == 0)  # full model residual map
    viewer, ch = ginga.show_image(image, chname=chname_resids, waveimg=waveimg,
                                  cuts = (-5.0, 5.0), bitmask = mask_in) #,wcs_match=True)
    show_trace(hdulist_1d, det_nm, viewer, ch)
    ginga.show_slits(viewer, ch, tslits_dict['slit_left'], tslits_dict['slit_righ'], slit_ids)#, args.det)


    # After displaying all the images sync up the images with WCS_MATCH
    shell = viewer.shell()
    out = shell.start_global_plugin('WCSMatch')
    out = shell.call_global_plugin_method('WCSMatch', 'set_reference_channel', [chname_resids], {})

    if args.embed:

        IPython.embed()
        # Playing with some mask stuff
        #out = shell.start_operation('TVMask')
        #maskfile = '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Shane_Kast_blue/600_4310_d55/shane_kast_blue_setup_A/crmask.fits'
        #out = shell.call_local_plugin_method(chname_resids, 'TVMask', 'load_file', [maskfile], {})






