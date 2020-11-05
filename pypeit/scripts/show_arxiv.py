#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""

def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Show an archived arc spectrum',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", type=str, help="WaveCalib JSON file")
    parser.add_argument('--det', default=1, type=int, help='Detector number')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Shows the spectrum
    """

    from matplotlib import pyplot as plt
    from pypeit.core.wavecal import waveio

    wave, flux, binspec = waveio.load_template(args.file, args.det)

    plt.clf()
    ax = plt.gca()
    ax.plot(wave, flux)
    plt.show()
