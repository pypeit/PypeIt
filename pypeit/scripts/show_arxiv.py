#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""
import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="WaveCalib JSON file")
    parser.add_argument('--det', default=1, type=int, help='Detector number')

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(pargs):
    """ Shows the spectrum
    """

    import sys
    import pdb

    from matplotlib import pyplot as plt
    from pypeit.core.wavecal import waveio

    wave, flux, binspec = waveio.load_template(pargs.file, pargs.det)

    plt.clf()
    ax = plt.gca()
    ax.plot(wave, flux)
    plt.show()
