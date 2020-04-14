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

from pypeit import slittrace
from pypeit.core.gui.skysub_regions import SkySubGUI
from pypeit.scripts import utils


def parser(options=None):

    parser = argparse.ArgumentParser(description='Display a Raw science image and interactively define'
                                                 'the sky regions using a GUI. Run in the same folder'
                                                 'as your .pypeit file',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('file', type = str, default=None, help='PypeIt file')
    parser.add_argument('--det', default=1, type=int, help="Detector")

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    # Generate a utilities class
    info = utils.Utilities(args)

    # Interactively select a science frame
    sciIdx = info.select_science_frame()

    # Load the spectrograph and parset
    info.load_spectrograph_parset(sciIdx)

    # Get the master key and directory
    mdir, mkey = info.get_master_dirkey(args.file)

    # Load the image data
    frame = info.load_frame(sciIdx)

    # Load the slits
    slits = slittrace.SlitTraceSet.from_master(mkey, mdir)

    # Derive an appropriate output filename
    file_base = info.get_basename(sciIdx)
    prefix = os.path.splitext(file_base)
    if prefix[1] == ".gz":
        prefix = os.path.splitext(prefix[0])[0]
    else:
        prefix = prefix[0]
    outname = "{0:s}_skyregions.fits".format(prefix)

    # Finally, initialise the GUI
    skyreg = SkySubGUI.initialize(args.det, frame, slits, outname=outname, runtime=False, printout=True)

    # Get the results
    skyreg.get_result()
