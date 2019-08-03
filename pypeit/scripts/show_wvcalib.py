#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""
import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="WaveCalib JSON file")
    parser.add_argument("slit", type=str)

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs, unit_test=False):
    """ Shows the spectrum
    """

    import sys
    import pdb

    from matplotlib import pyplot as plt
    from linetools import utils as ltu

    wvcalib = ltu.loadjson(pargs.file)

    # Grab it
    spec = wvcalib[pargs.slit]['spec']

    plt.clf()
    ax = plt.gca()
    ax.plot(spec)
    plt.show()
