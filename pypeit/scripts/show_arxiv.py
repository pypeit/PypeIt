#!/usr/bin/env python
"""
Wrapper to matplotlib to show an arc spectrum
"""

def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Show an archived arc spectrum in pypeit/data/arc_liens/reid_arxiv',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("file", type=str, help="Arxiv filename, e.g. gemini_gmos_r831_ham.fits")
    parser.add_argument('--det', default=1, type=int, help='Detector number')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):
    """ Shows the spectrum
    """
    import os
    from pkg_resources import resource_filename

    from matplotlib import pyplot as plt
    from pypeit.core.wavecal import waveio

    # Path
    if os.path.basename(args.file) == args.file:
        args.file = os.path.join(resource_filename('pypeit', 'data'),
                                 'arc_lines', 'reid_arxiv', args.file)

    wave, flux, binspec = waveio.load_template(args.file, args.det)

    plt.clf()
    ax = plt.gca()
    ax.plot(wave, flux)
    plt.show()


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
