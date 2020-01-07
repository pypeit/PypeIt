#!/usr/bin/env python
"""
Wrapper to the linetools XSpecGUI
"""
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

    from pypeit import specobjs
    from pypeit import msgs

    try:
        sobjs = specobjs.SpecObjs.from_fitsfile(args.file)

        # List only?
        if args.list:
            print("Showing object names for input file...")
            for ii in range(len(sobjs)):
                name = sobjs[ii].NAME
                print("EXT{:07d} = {}".format(ii+1, name))
            return

        # Load spectrum
        #spec = load.load_1dspec(args.file, exten=args.exten, extract=args.extract,
        #                          objname=args.obj, flux=args.flux)
        if args.obj is not None:
            exten = sobjs.name.index(args.obj)
            if exten < 0:
                msgs.error("Bad input object name: {:s}".format(args.obj))
        else:
            exten = args.exten-1 # 1-index in FITS file

        spec = sobjs[exten].to_xspec1d(extraction=args.extract, fluxed=args.flux)
    except:
        # place holder for coadd data model
        import numpy as np
        from pypeit import utils
        from linetools.spectra.xspectrum1d import XSpectrum1D
        from pypeit.core.telluric import general_spec_reader
        wave, counts, counts_ivar, counts_mask, meta_spec, head = general_spec_reader(args.file, ret_flam=False)
        spec = XSpectrum1D.from_tuple((wave, counts, np.sqrt(utils.inverse(counts_ivar))), masking='none')

    if unit_test is False:
        app = QApplication(sys.argv)
        # Screen dimensions
        width = app.desktop().screenGeometry().width()
        scale = 2. * (width/3200.)

    gui = XSpecGui(spec, unit_test=unit_test, screen_scale=scale)
    if unit_test is False:
        gui.show()
        app.exec_()

