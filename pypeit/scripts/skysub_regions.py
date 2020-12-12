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

def parse_args(options=None, return_parser=False):
    import argparse

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
    parser.add_argument('-s', '--standard', default=False, action='store_true',
                        help='List standard stars as well?')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os

    from pypeit.core.gui.skysub_regions import SkySubGUI
    from pypeit.core import flexure
    from pypeit.scripts import utils
    from pypeit import masterframe
    from pypeit.images import buildimage

    # Generate a utilities class
    info = utils.Utilities(args.file, args.det)

    # Interactively select a science frame
    sciIdx = info.select_science_frame(standard=args.standard)

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
    ext = buildimage.SkyRegions.master_file_format
    regfile = masterframe.construct_file_name(buildimage.SkyRegions, master_key=mkey, master_dir=mdir)
    regfile = regfile.replace(".{0:s}".format(ext), "_{0:s}.{1:s}".format(outname, ext))
    #outname = "{0:s}/MasterSkyRegions_{1:s}_{2:s}.fits.gz".format(mdir, mkey, outname)

    # Finally, initialise the GUI
    skyreg = SkySubGUI.initialize(args.det, frame, slits, info.spectrograph.pypeline,
                                  info.spectrograph.name, outname=regfile, overwrite=args.overwrite,
                                  runtime=False, printout=True, initial=args.initial,
                                  flexure=spat_flexure)

    # Get the results
    skyreg.get_result()
