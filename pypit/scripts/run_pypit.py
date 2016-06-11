#!/usr/bin/env python
#
# See top-level LICENSE file for Copyright information
#
# -*- coding: utf-8 -*-


"""
This script pushes a FITS file to ginga
"""
from __future__ import (print_function, absolute_import, division,
                        unicode_literals)

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

def parser(options=None):
    import argparse

    parser = argparse.ArgumentParser(description='Parse')
    parser.add_argument("pypit_file", type=str, help="PYPIT reduction file")
    parser.add_argument("-v", "--verbose", default=2, action="store_true",
                        help="(2) Level of verbosity (0-2)")
    parser.add_argument("-m", "--use_masters", default=False, help="Load previously generated MasterFrames", action="store_true")
    parser.add_argument("-d", "--develop", default=False, help="Turn develop debugging on", action="store_true")
    #parser.add_argument("-q", "--quick", default=False, help="Quick reduction", action="store_true")
    #parser.add_argument("-c", "--cpus", default=False, help="Number of CPUs for parallel processing", action="store_true")

    if options is None:
        args = parser.parse_args()
    else:
        args = parser.parse_args(options)
    return args


def main(args):

    import sys, os
    from pypit import pypit
    import traceback
    from pypit.armsgs import Messages as Initmsg

    # Import PYPIT routines
    from pypit import ardebug
    debug = ardebug.init()
    #debug['develop'] = True
    #debug['arc'] = True
    #debug['sky_sub'] = True
    #debug['trace'] = True
    #debug['obj_profile'] = True
    #debug['tilts'] = True
    #debug['flexure'] = True

    # Initiate logging for bugs and command line help
    # These messages will not be saved to a log file
    initmsgs = Initmsg(None, debug, 1)
    # Set the default variables
    qck = False
    cpu = 1
    #vrb = 2
    #use_masters = False

    #if len(sys.argv) < 2:
    #    initmsgs.usage(None)

    # Load options from command line
    """
        for o, a in opt:
            elif o in ('-q', '--quick'):
                qck = True
            elif o in ('-c', '--cpus'):
                cpu = int(a)
            elif o in ('-v', '--verbose'):
                vrb = int(a)
            elif o in ('-m', '--use_masters'):
                use_masters=True
            elif o in ('-d', '--develop'):
                debug['develop'] = True
        red = arg[0]
    except getopt.GetoptError, err:
        initmsgs.error(err.msg, usage=True)
    """
    debug['develop'] = debug['develop'] or args.develop
    splitnm = os.path.splitext(args.pypit_file)
    if splitnm[1] != '.pypit':
        initmsgs.error("Bad extension for PYPIT reduction file."+initmsgs.newline()+".pypit is required")
    logname = splitnm[0] + ".log"

    # Execute the reduction, and catch any bugs for printout
    if debug['develop']:
        pypit.PYPIT(args.pypit_file, progname=pypit.__file__, quick=qck, ncpus=cpu, verbose=args.verbose,
              use_masters=args.use_masters, logname=logname, debug=debug)
    else:
        try:
            pypit.PYPIT(args.pypit_file, progname=pypit.__file__, quick=qck, ncpus=cpu, verbose=args.verbose,
                  use_masters=args.use_masters, logname=logname, debug=debug)
        except:
            # There is a bug in the code, print the file and line number of the error.
            et, ev, tb = sys.exc_info()
            filename, line_no = "<filename>", "<line_no>"
            while tb:
                co = tb.tb_frame.f_code
                filename = str(co.co_filename)
                line_no = str(traceback.tb_lineno(tb))
                tb = tb.tb_next
            filename = filename.split('/')[-1]
            if str(ev) != "":
                initmsgs.bug("There appears to be a bug on Line " + line_no + " of " + filename + " with error:" +
                             initmsgs.newline() + str(ev) + initmsgs.newline() +
                             "---> please contact the authors")
            # Get armsgs instance to terminate
            from pypit.armsgs import get_logger
            get_logger().close()
