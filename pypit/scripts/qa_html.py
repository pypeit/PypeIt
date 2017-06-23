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
        # Generate MF file
        MF_filename = 'QA/MF_{:s}.html'.format(setup)
        body = ''
        with open(MF_filename,'w') as f:
            # Start
            links = arqa.html_init(f, 'QA  Setup {:s}: MasterFrame files'.format(setup))
            # Loop on calib_sets
            for cbset in cbsets:
                for det in dets:
                    # Run
                    new_links, new_body = arqa.html_mf_pngs(setup, cbset, det)
                    # Save
                    links += new_links
                    body += new_body
            # End
            arqa.html_end(f, body, links)
        #
        print("Wrote: {:s}".format(MF_filename))

    # Exposures
    if flg_exp:
        # Find all obj_trace files -- Not fool proof but ok
        obj_files = glob.glob('QA/PNGs/*obj_trace.png')
        # Parse for names
        names = []
        for obj_file in obj_files:
            i0 = obj_file.rfind('/')+1
            i1 = obj_file.rfind('_D')
            name = obj_file[i0:i1]
            names.append(name)
        uni_names = np.unique(names)
        # Loop
        for uni_name in uni_names:
            # Generate MF file
            exp_filename = 'QA/{:s}.html'.format(uni_name)
            body = ''
            with open(exp_filename,'w') as f:
                # Start
                links = arqa.html_init(f, 'QA for {:s}'.format(uni_name))
               # Loop on detector
                for det in range(1,99):
                    # Run
                    new_links, new_body = arqa.html_exp_pngs(uni_name, det)
                    # Save
                    links += new_links
                    body += new_body
                # End
                arqa.html_end(f, body, links)
            print("Wrote: {:s}".format(exp_filename))
