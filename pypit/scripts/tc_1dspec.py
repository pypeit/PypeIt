#!/usr/bin/env python

from pypit import pyputils
from pdb as debugger

msgs = pyputils.get_dummy_logger()

def parser(options=None):i

    description = (
                  'Script to telluric correct a spec!D file. '
                  'Currently only works for LRIS 800/10000 grating.'
                  )
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("infile", type=str, help="Input file (YAML)")
    parser.add_argument("--debug", default=False, i
                        action='store_true', help="Turn debugging on")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args

def main(args, unit_test=False, path=''):
