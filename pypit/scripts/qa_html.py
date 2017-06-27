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
    parser.add_argument("type", type=str, help="QA Type (MF, exp, all)")

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
    import numpy as np
    import glob
    from pypit import arqa

    # Flags
    flg_MF, flg_exp = False, False
    if args.type == 'MF':
        flg_MF = True
    elif args.type == 'exp':
        flg_exp = True
    elif args.type == 'all':
        flg_exp, flg_MF = True, True

    # Master Frame
    if flg_MF:
        arqa.gen_mf_html(args.pypit_file)

    # Exposures
    if flg_exp:
        arqa.gen_exp_html()

