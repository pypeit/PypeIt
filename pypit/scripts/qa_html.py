#!/usr/bin/env python

"""
Built HTML for PYPIT QA
"""
from pypit import pyputils
import pdb as debugger

msgs = pyputils.get_dummy_logger()
from numpy import isnan

def parser(options=None):

    import argparse

    parser = argparse.ArgumentParser(description='Script to build HTML files for PYPIT QA. [v1.0]')
    parser.add_argument("pypit_file", type=str, help="PYPIT file")
    parser.add_argument("type", type=str, help="QA Type (MF, Exposures, all)")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args, unit_test=False, path=''):
    """ Builds the HTML files
    path : str, optional
      Mainly for running the unit test
    """
    import yaml, glob
    from pypit import arqa
    from pypit.pypit import load_input

    # Parse PYPIT file (for setup)
    pyp_dict = load_input(args.pypit_file, msgs)
    debugger.set_trace()

    if args.type == 'MF':
        arqa.init_mf_html()

