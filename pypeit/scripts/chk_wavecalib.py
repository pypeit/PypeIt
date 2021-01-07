"""
This script displays the flat images in an RC Ginga window.
"""

def parse_args(options=None, return_parser=False):

    import argparse

    parser = argparse.ArgumentParser(description='Print QA on Wavelength Calib to the screen',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('master_file', type=str,
                        help='PypeIt MasterWaveCalib file [e.g. MasterWaveCalib_A_1_01.fits]')
    #parser.add_argument('--try_old', default=False, action='store_true',
    #                    help='Attempt to load old datamodel versions.  A crash may ensue..')

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    from pypeit import wavecalib

    # Load
    waveCalib = wavecalib.WaveCalib.from_file(args.master_file)#, chk_version=(not args.try_old))

    # Do it
    waveCalib.print_diagnostics()

