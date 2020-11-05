#!/usr/bin/env python
"""
Built HTML for PYPIT QA
"""

def parse_args(options=None, return_parser=False):
    import argparse

    parser = argparse.ArgumentParser(description='Script to build HTML files for PYPIT QA. [v1.0]',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("pypeit_file", type=str, help="PYPIT file")
    parser.add_argument("type", type=str, help="QA Type (MF, exp, all)")
    parser.add_argument("--qapath", type=str, default='QA/', help="Path the QA folder including QA/)")

    if return_parser:
        return parser

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args, unit_test=False, path=''):
    """ Builds the HTML files
    path : str, optional
      Mainly for running the unit test
    """

    from pypeit.core import qa

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
        qa.gen_mf_html(args.pypeit_file, args.qapath)

    # Exposures
    if flg_exp:
        qa.gen_exp_html()

