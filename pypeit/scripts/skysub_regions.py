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

from pypeit.core.gui.skysub_regions import SkySubGUI
from pypeit.core import flexure
from pypeit.scripts import utils


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display a Raw science image and interactively define'
                                                 'the sky regions using a GUI. Run in the same folder'
                                                 'as your .pypeit file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PypeIt file')
    parser.add_argument('--det', default=1, type=int, help="Detector")
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    parser.add_argument('-i', '--initial', default=False, action='store_true',
                        help='Use initial slit edges?')
    parser.add_argument('-f', '--flexure', default=False, action='store_true',
                        help='Use flexure corrected slit edges?')

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    # Generate a utilities class
    info = utils.Utilities(args.file, args.det)

    # Interactively select a science frame
    sciIdx = info.select_science_frame()

    # Load the spectrograph and parset
    info.load_spectrograph_parset(sciIdx)

    # Get the master key and directory
    mdir, mkey = info.get_master_dirkey()

    # Load the image data
    frame = info.load_frame(sciIdx)

    # Load the slits information
    slits = utils.get_slits(mkey, mdir)
    spat_flexure = None
    if args.flexure:
        spat_flexure = flexure.spat_flexure_shift(frame, slits)

    # Derive an appropriate output filename
    file_base = info.get_basename(sciIdx)
    prefix = os.path.splitext(file_base)
    if prefix[1] == ".gz":
        outname = os.path.splitext(prefix[0])[0]
    else:
        outname = prefix[0]
    outname = "{0:s}/MasterSkyRegions_{1:s}_{2:s}.fits.gz".format(mdir, mkey, outname)

    # Finally, initialise the GUI
    skyreg = SkySubGUI.initialize(args.det, frame, slits, info.spectrograph.pypeline, outname=outname, overwrite=args.overwrite,
                                  runtime=False, printout=True, initial=args.initial, flexure=spat_flexure)

    # Get the results
    skyreg.get_result()
