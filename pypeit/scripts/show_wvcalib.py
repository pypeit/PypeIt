#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""

def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Show the result of wavelength calibration',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", type=str, help="WaveCalib JSON file")
    parser.add_argument("slit", type=str)

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(pargs, unit_test=False):
    """ Shows the spectrum
    """

    from matplotlib import pyplot as plt
    from linetools import utils as ltu

    wvcalib = ltu.loadjson(pargs.file)

    # Grab it
    spec = wvcalib[pargs.slit]['spec']

    plt.clf()
    ax = plt.gca()
    ax.plot(spec)
    plt.show()
