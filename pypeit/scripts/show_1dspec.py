#!/usr/bin/env python

"""
Wrapper to the linetools XSpecGUI
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument("--exten", type=int, default=1, help="FITS extension")
    parser.add_argument("--obj", type=str, help="Object name in lieu of extension, e.g. O424-S1466-D02-I0013")
    parser.add_argument("--extract", type=str, default='OPT', help="Extraction method. Default is optimal. ['BOX', 'OPT']")
    parser.add_argument("--flux", default=False, action="store_true", help="Show fluxed spectrum?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False):
    """ Runs the XSpecGui on an input file
    """

    import sys
    import pdb

    from astropy.io import fits
    from PyQt5.QtWidgets import QApplication

    from linetools.guis.xspecgui import XSpecGui

    from pypeit.core import load

    # List only?
    if args.list:
        print("Showing object names for input file...")
        hdu = fits.open(args.file)
        for ii in range(1,len(hdu)):
            name = hdu[ii].name
            print("EXT{:07d} = {}".format(ii, name))
        return

    # Load spectrum
    spec = load.load_1dspec(args.file, exten=args.exten, extract=args.extract,
                              objname=args.obj, flux=args.flux)
    if unit_test is False:
        app = QApplication(sys.argv)
        # Screen dimensions
        width = app.desktop().screenGeometry().width()
        scale = 2. * (width/3200.)

    gui = XSpecGui(spec, unit_test=unit_test, screen_scale=scale)
    if unit_test is False:
        gui.show()
        app.exec_()

