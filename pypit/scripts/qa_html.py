#!/usr/bin/env python
"""
Built HTML for PYPIT QA
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse
import pdb as debugger

from pypit import arqa
from pypit import pyputils
msgs = pyputils.get_dummy_logger()

def parser(options=None):
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

