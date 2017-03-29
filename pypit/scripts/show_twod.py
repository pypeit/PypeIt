#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-

"""
This script enables the viewing of a processed FITS file
with extras.  Run above the Science/ folder.
"""

def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser(description='Display spec2d image in a Ginga viewer.  Run above the Science/ folder',
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

    # List only?
    from astropy.io import fits
    from astropy.table import Table
    hdu = fits.open(args.file)
    head0 = hdu[0].header
    if args.list:
        print(hdu.info())
        return

    # Setup for PYPIT imports
    from pypit import pyputils
    msgs = pyputils.get_dummy_logger()
    from pypit import arparse as settings  # Has to come after the logger
    from pypit import ginga as pyp_ginga
    from pypit import armasters
    import pdb

    # One detector, sky sub for now
    names = [hdu[i].name for i in range(len(hdu))]
    exten = names.index('DET{:d}-SKYSUB'.format(args.det))
    skysub = hdu[exten].data

    # Show Image
    viewer, ch = pyp_ginga.show_image(skysub, chname='DET-{:02d}'.format(args.det))

    # Add slits
    '''
    testing = False
    if testing:
        mdir = 'MF_lris_blue/'
        setup = 'A_{:02d}_aa'.format(args.det)
    else:
        mdir = head0['PYPMFDIR']+'/'
        setup = '{:s}_{:02d}_{:s}'.format(head0['PYPCNFIG'], args.det, head0['PYPCALIB'])
    trc_file = armasters.master_name('trace', setup, mdir=mdir)
    trc_hdu = fits.open(trc_file)
    lordloc = trc_hdu[1].data  # Should check name
    rordloc = trc_hdu[2].data  # Should check name
    pyp_ginga.show_slits(viewer, ch, lordloc, rordloc)#, args.det)

    # Object traces
    spec1d_file = args.file.replace('spec2d', 'spec1d')
    hdulist_1d = fits.open(spec1d_file)
    det_nm = 'D{:02d}'.format(args.det)
    for hdu in hdulist_1d:
        if det_nm in hdu.name:
            tbl = Table(hdu.data)
            trace = tbl['obj_trace']
            pyp_ginga.show_trace(viewer, ch, trace, hdu.name)
    '''

