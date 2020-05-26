#!/usr/bin/env python
"""
Wrapper to the linetools XSpecGUI
"""
import argparse
import sys
from linetools.guis.xspecgui import XSpecGui
from PyQt5.QtWidgets import QApplication

from pypeit import specobjs
from pypeit import msgs
import numpy as np
from IPython import embed

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument("--exten", type=int, default=1, help="FITS extension")
    parser.add_argument("--obj", type=str, help="Object name in lieu of extension, e.g. SPAT0424-SLIT0000-DET01")
    parser.add_argument("--extract", type=str, default='OPT', help="Extraction method. Default is OPT. ['BOX', 'OPT']")
    parser.add_argument("--flux", default=False, action="store_true", help="Show fluxed spectrum?")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args):
    """ Runs the XSpecGui on an input file
    """

    sobjs = specobjs.SpecObjs.from_fitsfile(args.file)
    # List only?
    if args.list:
        print("Showing object names for input file...")
        for ii in range(len(sobjs)):
            name = sobjs[ii].NAME
            print("EXT{:07d} = {}".format(ii+1, name))
        return

    if args.obj is not None:
        exten = np.where(sobjs.NAME == args.obj)[0][0]
        if exten < 0:
            msgs.error("Bad input object name: {:s}".format(args.obj))
    else:
        exten = args.exten-1 # 1-index in FITS file

    # Check Extraction
    if args.extract == 'OPT':
        if sobjs[exten]['OPT_WAVE'] is None: #not in sobjs[exten]._data.keys():
                msgs.error("Spectrum not extracted with OPT.  Try --extract BOX")

    spec = sobjs[exten].to_xspec1d(extraction=args.extract, fluxed=args.flux)

    # Setup
    app = QApplication(sys.argv)
    # Screen dimensions
    width = app.desktop().screenGeometry().width()
    scale = 2. * (width/3200.)

    # Launch
    gui = XSpecGui(spec, screen_scale=scale)
    gui.show()
    app.exec_()

