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
    parser.add_argument("-v", "--verbosity", type=int, default=2,
                        help="(2) Level of verbosity (0-2)")
    parser.add_argument("-d", "--develop", default=False, action='store_true',
                        help="Turn develop debugging on")
    parser.add_argument("--extension", default='.fits',
                        help='File extension; compression indicators (e.g. .gz) not required.')
    parser.add_argument("--pypit_file", default=False, action='store_true',
                        help='Input is the .pypit file')
    parser.add_argument("--redux_path", default=None,
                        help='Path to reduction folder.  Default is current working directory.')
    parser.add_argument("-c", "--custom", default=False, action='store_true',
                        help='Generate custom folders and pypit files?')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
#    parser.add_argument("-q", "--quick", default=False, help="Quick reduction",
#                        action="store_true")
#    parser.add_argument("-c", "--cpus", default=False,
#                        help="Number of CPUs for parallel processing", action="store_true")
#    parser.print_help()

    return parser.parse_args() if options is None else parser.parse_args(options)


def main(args):

    import os
    import datetime
    import pdb as debugger

    from pypit import msgs
    from pypit.spectrographs.util import valid_spectrographs
    from pypit.core import arsetup
    from pypit.par.util import make_pypit_file, parse_pypit_file
    from pypit.scripts import run_pypit

    # Check that input spectrograph is supported
    instruments_served = valid_spectrographs()
    if args.spectrograph not in instruments_served:
        raise ValueError('Instrument \'{0}\' unknown to PYPIT.\n'.format(args.spectrograph)
                         + '\tAvailable options are: {0}\n'.format(', '.join(instruments_served))
                         + '\tSelect an available instrument or consult the documentation '
                         + 'on how to add a new instrument.')

    # setup_files dir
    redux_path = os.getcwd() if args.redux_path is None else args.redux_path
    outdir = os.path.join(redux_path, 'setup_files')
    msgs.info('Setup files will be written to: {0}'.format(outdir))
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
    msgs.info('Wrote template pypit file: {0}'.format(pypit_file))

    # Parser
    pinp = [pypit_file, '-p', '-s {0}'.format(root) ]
    if args.overwrite:
        pinp += ['-o']
    if args.develop:
        pinp += ['-d']
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
    _, data_files, frametype, setups = parse_pypit_file(pypit_file)

    # Get paths
    paths = []
    for data_file in data_files:
        islsh = data_file.rfind('/')
        path = data_file[:islsh+1]
        if path not in paths:
            paths.append(path)

    # Generate .pypit files and sub-folders
    all_setups, all_setuplines, all_setupfiles = arsetup.load_sorted(sorted_file)
    for setup, setup_lines, sorted_files in zip(all_setups, all_setuplines, all_setupfiles):
        root = args.spectrograph+'_setup_'
        # Make the dir
        newdir = os.path.join(redux_path, root+setup)
        if not os.path.exists(newdir):
            os.mkdir(newdir)
        # Now the file
        pypit_file = os.path.join(newdir, root+setup+'.pypit')
        # Modify parlines
        for kk in range(len(cfg_lines)):
            if 'sortroot' in cfg_lines[kk]:
                cfg_lines[kk] = '    sortroot = {0}'.format(root+setup)

        make_pypit_file(pypit_file, args.spectrograph, [], cfg_lines=cfg_lines,
                        setup_lines=setup_lines, sorted_files=sorted_files, paths=paths)
        print("Wrote {:s}".format(pypit_file))

