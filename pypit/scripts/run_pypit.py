#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-
"""
This script runs PYPIT
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
from __future__ import unicode_literals

import os
import sys
import argparse
import traceback

from pypit import pypit
from pypit import msgs
from pypit import ardebug
#from pypit.scripts import qa_html

# Globals
#debug = ardebug.init()
#debug['develop'] = True
#debug['arc'] = True
#debug['sky_sub'] = True
#debug['trace'] = True
#debug['obj_profile'] = True
#debug['trace_obj'] = True
#debug['tilts'] = True
#debug['flexure'] = True
#debug['no_qa'] = True

#from pypit.armsgs import Messages as Initmsg
#initmsgs = Initmsg(None, debug, 1)


def parser(options=None):

    parser = argparse.ArgumentParser(description=msgs.usage('PYPIT'),
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('pypit_file', type=str,
                        help='PYPIT reduction file (must have .pypit extension)')
    parser.add_argument('-v', '--verbosity', type=int, default=2,
                        help='(2) Level of verbosity (0-2)')
    parser.add_argument('-m', '--use_masters', default=False, action='store_true',
                        help='Load previously generated MasterFrames')
    parser.add_argument('-d', '--develop', default=False, action='store_true',
                        help='Turn develop debugging on')
    parser.add_argument('--devtest', default=False, action='store_true',
                        help='Running development tests')
    parser.add_argument('--debug_arc', default=False, action='store_true',
                        help='Turn wavelength/arc debugging on')
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

    # Initiate logging for bugs and command line help
    # These messages will not be saved to a log file
    # Set the default variables
    qck = False
    cpu = 1
    #vrb = 2

    # Load options from command line
    debug = ardebug.init()
    debug['develop'] = debug['develop'] or args.develop
    debug['arc'] = debug['arc'] or args.debug_arc
    splitnm = os.path.splitext(args.pypit_file)
    if splitnm[1] != '.pypit':
        msgs.error("Bad extension for PYPIT reduction file."+msgs.newline()+".pypit is required")
    logname = splitnm[0] + ".log"

    # Execute the reduction, and catch any bugs for printout
    if debug['develop']:
        pypit.PYPIT(args.pypit_file, progname=pypit.__file__, quick=qck, ncpus=cpu,
                    verbosity=args.verbosity, use_masters=args.use_masters, devtest=args.devtest,
                    logname=logname, debug=debug)
    else:
        try:
            pypit.PYPIT(args.pypit_file, progname=pypit.__file__, quick=qck, ncpus=cpu,
                        verbosity=args.verbosity, use_masters=args.use_masters,
                        devtest=args.devtest, logname=logname, debug=debug)
        except:
            # There is a bug in the code, print the file and line number of the error.
            et, ev, tb = sys.exc_info()
            filename, line_no = "<filename>", "<line_no>"
            while tb:
                co = tb.tb_frame.f_code
                filename = str(co.co_filename)
                try:
                    line_no = str(traceback.tb_lineno(tb))
                except AttributeError:  # Python 3
                    line_no = 'UNDEFINED'
                tb = tb.tb_next
            filename = filename.split('/')[-1]
            if str(ev) != "":
                msgs.bug("There appears to be a bug on Line " + line_no + " of " + filename
                         + " with error:" + msgs.newline() + str(ev) + msgs.newline()
                         + "---> please contact the authors")
            msgs.close()
#            # Get armsgs instance to terminate
#            from pypit.armsgs import get_logger
#            # TODO: Close logger and close QA
#            get_logger().close()
            return 1
    return 0
