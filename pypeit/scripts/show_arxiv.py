#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""
import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("file", type=str, help="WaveCalib JSON file")
    parser.add_argument('--det', default=1, type=int, help='Detector number')

    return parser.parse_args() if options is None else parser.parse_args(options)


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
