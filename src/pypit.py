#import matplotlib
#matplotlib.use('Agg')  # For Travis

import os
import sys
import getopt
from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
import traceback

# Import PYPIT routines
debug = True
last_updated = "26 November 2015"
version = '0.3'

# Init logger
from armsgs import Messages as Messages
import armsgs
msgs = armsgs.get_logger((None, debug, last_updated, version))

import arload

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except:
    pass

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger



def PYPIT(argflag, quick=False):
    """
    Main driver of the PYPIT code. Default settings and
    user-specified changes are made, and passed to the
    appropriate code for data reduction.

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    quick : bool
      If True, a quick reduction (but possibly less
      accurate) will be performed. This flag is most
      useful for observing at a telescope, but not
      for publication quality results.
    ---------------------------------------------------
    """
    # First send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
    sigsignal(SIGINT, msgs.signal_handler)

    # Ignore all warnings given by python
    resetwarnings()
    simplefilter("ignore")

    # Record the starting time
    tstart = time()

    # Load the Input file
    argflag, parlines, datlines, spclines = arload.load_input(argflag)

    # If a quick reduction has been requested, make sure the requested pipeline
    # is the quick implementation (if it exists), otherwise run the standard pipeline.
    if quick:
        # Change to a "quick" settings file
        msgs.work("TO BE DONE")

    # Load the Spectrograph settings
    spect = arload.load_spect(argflag)

    # Load any changes to the spectrograph settings
    spect = arload.load_spect(argflag, spect=spect, lines=spclines)

    # Load the important information from the fits headers
    fitsdict = arload.load_headers(argflag, spect, datlines)

    # Reduce the data!
    status = 0
    msgs.work("Make appropriate changes to quick reduction")
    if quick:
        msgs.work("define what is needed here")
    # Send the data away to be reduced
    if spect['mosaic']['reduction'] == 'ARMLSD':
        msgs.info("Data reduction will be performed using PYPIT-ARMLSD")
        import armlsd
        status = armlsd.ARMLSD(argflag, spect, fitsdict)
    elif spect['mosaic']['reduction'] == 'ARMED':
        msgs.info("Data reduction will be performed using PYPIT-ARMED")
        import armed
        status = armed.ARMED(argflag, spect, fitsdict, msgs)
    # Check for successful reduction
    if status == 0:
        msgs.info("Data reduction complete")
    else:
        msgs.error("Data reduction failed with status ID {0:d}".format(status))
    # Capture the end time and print it to user
    tend = time()
    codetime = tend-tstart
    if codetime < 60.0:
        msgs.info("Data reduction execution time: {0:.2f}s".format(codetime))
    elif codetime/60.0 < 60.0:
        mns = int(codetime/60.0)
        scs = codetime - 60.0*mns
        msgs.info("Data reduction execution time: {0:d}m {1:.2f}s".format(mns, scs))
    else:
        hrs = int(codetime/3600.0)
        mns = int(60.0*(codetime/3600.0 - hrs))
        scs = codetime - 60.0*mns - 3600.0*hrs
        msgs.info("Data reduction execution time: {0:d}h {1:d}m {2:.2f}s".format(hrs, mns, scs))
    return


if __name__ == "__main__":
    argflag = dict({})
    prognm = sys.argv[0]
    debug = True
    quick = False
    cpus = 1
    verbose = 2

    # Init logger
    #msgs = armsgs.get_logger((None, debug, last_updated, version))

    if len(sys.argv) < 2:
        msgs.usage(None)

    # Load options from command line
    try:
        opt, arg = getopt.getopt(sys.argv[1:], 'hqc:v:', ['help',
                                                          'quick',
                                                          'cpus',
                                                          'use_master',
                                                          'verbose'])
        for o, a in opt:
            if o in ('-h', '--help'):
                msgs.usage(None)
            elif o in ('-q', '--quick'):
                quick = True
            elif o in ('-c', '--cpus'):
                cpus = int(a)
            elif o in ('-v', '--verbose'):
                verbose = int(a)
            elif o in ('--use_master'):
                verbose = str(a)
        lname = os.path.splitext(arg[0])[0] + ".log"
        argflag = arload.optarg(sys.argv)
        argflag['run']['ncpus'] = cpus
        argflag['out']['verbose'] = verbose
    except getopt.GetoptError, err:
        msgs.error(err.msg, usage=True)

    if debug:
        PYPIT(argflag, quick=quick)
    else:
        try:
            PYPIT(argflag, quick=quick)
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
            msgs.bug("There appears to be a bug on Line " + line_no + " of " + filename + " with error:" +
                     msgs.newline() + str(ev) + msgs.newline() +
                     "---> please contact the authors")
