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
        # Read calib file
        calib_file = args.pypit_file.replace('.pypit', '.calib')
        with open(calib_file, 'r') as infile:
            calib_dict = yaml.load(infile)
        # Parse
        setup = list(calib_dict.keys())[0]
        dets, cbsets = [], []
        for key in calib_dict[setup].keys():
            if key == '--':
                continue
            try:
                dets.append(int(key))
            except ValueError:
                cbsets.append(key)
        #
        MF_filename = 'QA/MF_{:s}.html'.format(setup)
        links = '<h2>Quick Links</h2>\n'
        links += '<ul>\n'
        body = ''
        with open(MF_filename,'w') as f:
            # Start
            arqa.html_mf_init(f, 'QA  Setup {:s}: MasterFrame files'.format(setup))
            # Loop on calib_sets
            for cbset in cbsets:
                for det in dets:
                    # Run
                    new_links, new_body = arqa.html_mf_pngs(setup, cbset, det)
                    # Save
                    links += new_links
                    body += new_body
            # Write links
            f.write(links)
            f.write('</ul>\n')
            f.write('<hr>\n')
            # Write rest of body
            f.write(body)
            # Finish
            f.write(arqa.html_end())
        #
        print("Wrote: {:s}".format(MF_filename))

