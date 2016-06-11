#!/usr/bin/env python

"""
Wrapper to the linetools XSpecGUI
"""


def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--list", default=False, help="List the extensions only?", action="store_true")
    parser.add_argument("--exten", type=int, help="FITS extension")
    parser.add_argument("--optimal", default=False,
                        help="Show Optimal? Default is boxcar", action="store_true")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False):
    """ Runs the XSpecGui on an input file
    """
    import sys

    # List only?
    if args.list:
        from astropy.io import fits
        hdu = fits.open(args.file)
        print(hdu.info())
        return

    from linetools.guis.xspecgui import XSpecGui

    # Extension
    exten = (args.exten if hasattr(args, 'exten') else 0)

    # Read spec keywords
    rsp_kwargs = {}
    if args.optimal:
        rsp_kwargs['wave_tag'] = 'opt_wave'
        rsp_kwargs['flux_tag'] = 'opt_counts'
        rsp_kwargs['var_tag'] = 'opt_var'
    else:
        rsp_kwargs['wave_tag'] = 'box_wave'
        rsp_kwargs['flux_tag'] = 'box_counts'
        rsp_kwargs['var_tag'] = 'box_var'

    if unit_test is False:
        from PyQt4 import QtGui
        app = QtGui.QApplication(sys.argv)

    gui = XSpecGui(args.file, exten=exten, rsp_kwargs=rsp_kwargs, unit_test=unit_test)
    if unit_test is False:
        gui.show()
        app.exec_()
