#!/usr/bin/env python

"""
Plots an extracted sky spectrum with an archived one
  Probably most useful for exploring sky spectra in the blue
"""

def parse_args(options=None, return_parser=False):
    import argparse
    parser = argparse.ArgumentParser(description='Compare the extracted sky spectrum against an '
                                                 'archived sky model maintained by PypeIt.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('file', type=str, help='Spectral file')
    parser.add_argument('skyfile', type=str,
                        help='Archived PypeIt sky file (e.g. paranal_sky.fits)')
    parser.add_argument('--exten', type=int, help='FITS extension')
    parser.add_argument('--optimal', default=False, action='store_true',
                        help='Show Optimal? Default is boxcar')
    parser.add_argument('--scale_user', default=1., type=float,
                        help='Scale user spectrum by a factor')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)

# Script to run XSpec from the command line or ipython
def main(args):

    import os
    from pkg_resources import resource_filename

    from matplotlib import pyplot as plt

    from linetools.spectra.io import readspec

    # Path to archived sky spectra
    sky_path = os.path.join(resource_filename('pypeit', 'data'), 'sky_spec')

    # Extension
    exten = args.exten if hasattr(args, 'exten') else 0

    # Read spec keywords
    ikwargs = {}
    if args.optimal:
        ikwargs['wave_tag'] = 'opt_wave'
        ikwargs['flux_tag'] = 'opt_sky'
    else:
        ikwargs['wave_tag'] = 'box_wave'
        ikwargs['flux_tag'] = 'box_sky'

    # Load user file
    user_sky = readspec(args.file, exten=exten, **ikwargs)
    # Load sky spec
    arx_sky = readspec(sky_path+args.skyfile)

    # Plot
    plt.clf()
    plt.plot(user_sky.wavelength, user_sky.flux*args.scale_user, 'k-', label='user')
    plt.plot(arx_sky.wavelength, arx_sky.flux, 'b-', label='archive')
    legend = plt.legend(loc='upper left', scatterpoints=1, borderpad=0.3,
                        handletextpad=0.3, fontsize='small', numpoints=1)
    plt.show()

