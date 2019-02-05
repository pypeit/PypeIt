#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PypeIt
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import argparse

from pypeit import msgs


def parser(options=None):

    parser = argparse.ArgumentParser(description=msgs.usage('PypeIt'),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pypeit_file', type=str,
                        help='PypeIt reduction file (must have .pypeit extension)')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='Verbosity level between 0 [none] and 2 [all]')
    # JFH TODO Are the -t and -r keyword still valid given that run_pypeit no longer runs setup?
    parser.add_argument('-t', '--hdrframetype', default=False, action='store_true',
                        help='Use file headers and the instument-specific keywords to determine'
                             'the type of each frame')
    parser.add_argument('-r', '--sort_dir', default=None,
                        help='Directory used to store the sorted files.  Default is to omit '
                             'writing these files.')
    parser.add_argument('-m', '--use_masters', default=False, action='store_true',
                        help='Load previously generated MasterFrames')
    parser.add_argument('-s', '--show', default=False, action='store_true',
                        help='Show reduction steps via plots (which will block further execution until clicked on) '
                             'and outputs to ginga. Requires remote control ginga session via "ginga --modules=RC &"')
    # JFH Should the default now be true with the new definition.
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    group = parser.add_mutually_exclusive_group()
#    group.add_argument('-p', '--prep_setup', default=False, action='store_true',
#                       help='Run pypeit to prepare the setup only')
#    group.add_argument('-c', '--calcheck', default=False, action='store_true',
#                       help='Run pypeit only as a check on the calibrations')
    group.add_argument('-d', '--detector', default=None, help='Detector to limit reductions on.  If the output files exist and -o is used, the outputs for the input detector will be replaced.')

#    parser.add_argument('-q', '--quick', default=False, help='Quick reduction',
#                        action='store_true')
#    parser.add_argument('-c', '--cpus', default=False, action='store_true',
#                         help='Number of CPUs for parallel processing')
#    parser.print_help()

    if options is None:
        pargs = parser.parse_args()
    else:
        pargs = parser.parse_args(options)
    return pargs


def main(args):

    import os
    import sys
    import traceback

    from pypeit import check_requirements
    from pypeit import pypeit
    from pypeit import pypeitsetup
    from pypeit import debugger

    # Initiate logging for bugs and command line help
    # These messages will not be saved to a log file
    # Set the default variables
    qck = False
    cpu = 1
    #vrb = 2

    # Load options from command line
    splitnm = os.path.splitext(args.pypeit_file)
    if splitnm[1] != '.pypeit':
        msgs.error("Bad extension for PypeIt reduction file."+msgs.newline()+".pypeit is required")
    logname = splitnm[0] + ".log"

    # Load PypeIt file to get the spectrograph (might happen twice but that is ok)
    pypeitSetup = pypeitsetup.PypeItSetup.from_pypeit_file(args.pypeit_file)

    pypeIt = pypeit.PypeIt(args.pypeit_file,verbosity=args.verbosity,reuse_masters=args.use_masters,
                           overwrite=args.overwrite, logname=logname, show=args.show)

    # JFH I don't see why this is an optional argument here. We could allow the user to modify an infinite number of parameters
    # from the command line? Why do we have the PypeIt file then? This detector can be set in the pypeit file.
    # Detector?
    if args.detector is not None:
        msgs.info("Restricting reductions to detector={}".format(args.detector))
        pypeIt.par['rdx']['detnum'] = int(args.detector)

    pypeIt.reduce_all()
    msgs.info('Data reduction complete')
    # QA HTML
    msgs.info('Generating QA HTML')
    pypeIt.build_qa()

    return 0

