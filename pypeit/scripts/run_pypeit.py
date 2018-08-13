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
    parser.add_argument('-t', '--hdrframetype', default=False, action='store_true',
                        help='Use file headers and the instument-specific keywords to determine'
                             'the type of each frame')
    parser.add_argument('-s', '--sort_dir', default=None,
                        help='Directory used to store the sorted files.  Default is to omit '
                             'writing these files.')
    parser.add_argument('-m', '--use_masters', default=False, action='store_true',
                        help='Load previously generated MasterFrames')
#    parser.add_argument('--devtest', default=False, action='store_true',
#                        help='Running development tests')
    parser.add_argument('--debug_arc', default=False, action='store_true',
                        help='Turn wavelength/arc debugging on')
    parser.add_argument('-o', '--overwrite', default=False, action='store_true',
                        help='Overwrite any existing files/directories')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-p', '--prep_setup', default=False, action='store_true',
                       help='Run pypeit to prepare the setup only')
    group.add_argument('-c', '--calcheck', default=False, action='store_true',
                       help='Run pypeit only as a check on the calibrations')

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

    pypeit.PypeIt(args.pypeit_file, setup_only=args.prep_setup, calibration_check=args.calcheck,
                  use_header_frametype=args.hdrframetype, sort_dir=args.sort_dir, overwrite=args.overwrite,
                  verbosity=args.verbosity, use_masters=args.use_masters, logname=logname)

    return 0

#    # Execute the reduction, and catch any bugs for printout
#    if debug['develop']:
#        pypeit.PYPIT(args.pypeit_file, progname=pypeit.__file__, quick=qck, ncpus=cpu,
#                    verbosity=args.verbosity, use_masters=args.use_masters, devtest=args.devtest,
#                    logname=logname)
#    else:
#        try:
#            pypeit.PYPIT(args.pypeit_file, progname=pypeit.__file__, quick=qck, ncpus=cpu,
#                        verbosity=args.verbosity, use_masters=args.use_masters,
#                        devtest=args.devtest, logname=logname)
#        except:
#            # There is a bug in the code, print the file and line number of the error.
#            et, ev, tb = sys.exc_info()
#            filename, line_no = "<filename>", "<line_no>"
#            while tb:
#                co = tb.tb_frame.f_code
#                filename = str(co.co_filename)
#                try:
#                    line_no = str(traceback.tb_lineno(tb))
#                except AttributeError:  # Python 3
#                    line_no = 'UNDEFINED'
#                tb = tb.tb_next
#            filename = filename.split('/')[-1]
#            if str(ev) != "":
#                msgs.bug("There appears to be a bug on Line " + line_no + " of " + filename
#                         + " with error:" + msgs.newline() + str(ev) + msgs.newline()
#                         + "---> please contact the authors")
#            msgs.close()
#            return 1
#    return 0


