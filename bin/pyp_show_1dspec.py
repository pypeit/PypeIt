#!/usr/bin/env python

"""
Plot a spectrum with an interactive QT GUI
"""
import pdb
import sys


# Script to run XSpec from the command line or ipython
def main(*args, **kwargs):
    """ Runs the XSpecGui on an input file
    """
    import argparse

    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="Spectral file")
    parser.add_argument("--exten", type=int, help="FITS extension")
    parser.add_argument("--optim", default=False,
                        help="Show Optimal? Default is boxcar", action="store_true")

    pargs = parser.parse_args()

    from PyQt4 import QtGui
    from linetools.guis.xspecgui import XSpecGui

    # Extension
    exten = (pargs.exten if hasattr(pargs, 'exten') else 0)

    # Read spec keywords
    rsp_kwargs = {}
    if pargs.optim:
        rsp_kwargs['wave_tag'] = 'opt_wave'
        rsp_kwargs['flux_tag'] = 'opt_counts'
        rsp_kwargs['var_tag'] = 'opt_var'
    else:
        rsp_kwargs['wave_tag'] = 'box_wave'
        rsp_kwargs['flux_tag'] = 'box_counts'
        rsp_kwargs['var_tag'] = 'box_var'

    app = QtGui.QApplication(sys.argv)

    gui = XSpecGui(pargs.file, exten=exten, rsp_kwargs=rsp_kwargs)
    gui.show()
    app.exec_()

if __name__ == '__main__':
    main()
