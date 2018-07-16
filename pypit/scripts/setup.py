#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-


"""
This script generates files to setup a PYPIT run
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

import argparse

def parser(options=None):
    parser = argparse.ArgumentParser(description="Script to setup a PYPIT run [v2]")
    parser.add_argument("files_root", type=str, help="File path+root, e.g. /data/Kast/b ")
    parser.add_argument("spectrograph", type=str, help="Name of spectrograph")
    parser.add_argument("-v", "--verbosity", type=int, default=2, help="(2) Level of verbosity (0-2)")
    parser.add_argument("-d", "--develop", default=False, action='store_true', help="Turn develop debugging on")
    parser.add_argument("--extension", default='.fits',
                        help="Extension for data files.  Note any extension for compression (e.g. .gz) is not required.")
    parser.add_argument("--pypit_file", default=False, action='store_true', help='Input is the .pypit file')
    parser.add_argument("--redux_path", default='./', help='Path to reduction folder (Mainly for tests)')
    parser.add_argument("-c", "--custom", default=False, action='store_true', help='Generate custom folders and pypit files?')
    #parser.add_argument("-q", "--quick", default=False, help="Quick reduction", action="store_true")
    #parser.add_argument("-c", "--cpus", default=False, help="Number of CPUs for parallel processing", action="store_true")
    #parser.print_help()

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    #
    return pargs


def main(args):

    import os
    import datetime
    import pdb as debugger

    from pypit import msgs
    from pypit import armeta
    from pypit.core import arsetup
    from pypit.par.util import make_pypit_file, parse_pypit_file
    from pypit.scripts import run_pypit

    # Check that input spectrograph is supported
    # TODO: This is outdated!
    if args.spectrograph not in armeta.instr_list():
        print("-------------------------------------------------------------")
        print("Input instrument {:s} is not supported by PYPIT".format(args.spectrograph))
        print("Here is the list of options: ")
        print(armeta.instr_list())
        raise IOError("Consult the documentation for further info.")

    # setup_files dir
    outdir = args.redux_path+'/setup_files'
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # Generate a dummy .pypit file
    date = str(datetime.date.today().strftime('%Y-%b-%d'))
    root = args.spectrograph+'_'+date
    pypit_file = outdir+'/'+root+'.pypit'
    # Generate
    dfname = "{:s}*{:s}*".format(args.files_root, args.extension)
    # configuration lines
    cfg_lines = ['[rdx]']
    cfg_lines += ['    spectrograph = {0}'.format(args.spectrograph)]
    cfg_lines += ['    sortroot = {0}'.format(root)]
    make_pypit_file(pypit_file, args.spectrograph, [dfname], cfg_lines=cfg_lines, setup_mode=True)
    print("Wrote {:s}".format(pypit_file))

    # Parser
    pinp = [pypit_file]
    if args.develop:
        pinp += ['-d', '-o']
    pargs = run_pypit.parser(pinp)
    sorted_file = pypit_file.replace('.pypit', '.sorted')

    # Run
    run_pypit.main(pargs)

    # #####################
    # Generate custom .pypit files
    if not args.custom:
        return

    msgs.reset(verbosity=0)

    # Read master file
    cfg_lines, data_files, frametype, setups = parse_pypit_file(filename)
    
    pyp_dict = load_input(pypit_file, msgs)
    parlines, datlines, spclines, dfnames = [pyp_dict[ii] for ii in ['par','dat','spc','dfn']]

    # Get paths
    paths = []
    for datafile in data_files:
        islsh = datline.rfind('/')
        path = datline[:islsh+1]
        if path not in paths:
            paths.append(path)

    # Generate .pypit files and sub-folders
    all_setups, all_setuplines, all_setupfiles = arsetup.load_sorted(sorted_file)
    for setup, setup_lines, sorted_files in zip(all_setups, all_setuplines, all_setupfiles):
        root = args.spectrograph+'_setup_'
        # Make the dir
        newdir = args.redux_path+root+setup
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        # Now the file
        pypit_file = newdir+'/'+root+setup+'.pypit'
        # Modify parlines
        for kk in range(len(cfg_lines)):
            if 'sortroot' in cfg_lines[kk]:
                cfg_lines[kk] = ['    sortroot = {0}'.format(root+setup)]

        make_pypit_file(pypit_file, args.spectrograph, [], cfg_lines=cfg_lines,
                        setup_lines=setup_lines, sorted_files=sorted_files, paths=None)
        print("Wrote {:s}".format(pypit_file))

