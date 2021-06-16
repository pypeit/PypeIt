#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script displays the Arc image from and the traces
in an RC Ginga window (must be previously launched)

.. todo::

    This script is out of date and will not function with the current
    pypeit master branch.

"""

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Display MasterArc image in a previously '
                                                 'launched RC Ginga viewer with tilts',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('option', type = str, default = None,
                        help='Item to show [fweight, model, tilts, final_tilts]')
    parser.add_argument('master_file', type = str, default = None,
                        help='path to Master file, e.g. Masters/MasterTilts_C_1_03.fits')
    parser.add_argument('--slit', type=int, default=None, help='Slit/Order [0,1,2..]')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    from pypeit import wavetilts
    from pypeit.display import display

    # Load up
    wTilts = wavetilts.WaveTilts.from_file(args.master_file)

    # Connect to the ginga viewer
    # TODO: I don't think raise_err needs to be specified here, but
    # this is what's done in show_2dspec
    display.connect_to_ginga(raise_err=True, allow_new=True)

    # Show
    cname = None if args.slit is None else 'Slit{:03d}'.format(args.slit)
    wTilts.show(args.option, slit=args.slit, cname=cname)


def entry_point():
    main(parse_args())


if __name__ == '__main__':
    entry_point()
