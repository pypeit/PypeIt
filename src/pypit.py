#import matplotlib
#matplotlib.use('Agg')  # For Travis

import os
import sys
import getopt
from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
import traceback
import arparse

# Import PYPIT routines
import ardebug
debug = ardebug.init()
#debug['develop'] = True
#debug['arc'] = True
#debug['sky_sub'] = True
#debug['trace'] = True
#debug['obj_profile'] = True
#debug['tilts'] = True
#debug['flexure'] = True
last_updated = "2 May 2016"
version = '0.6'

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except ImportError:
    pass

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger


def PYPIT(redname, progname=__file__, quick=False, ncpus=1, verbose=1,
          logname=None, use_masters=False):
    """
    Main driver of the PYPIT code. Default settings and
    user-specified changes are made, and passed to the
    appropriate code for data reduction.

    Parameters
    ----------
    redname : string
      Input reduction script
    progname : string
      Name of the program
    quick : bool
      If True, a quick reduction (but possibly less
      accurate) will be performed. This flag is most
      useful for observing at a telescope, but not
      for publication quality results.
    ncpus : int
      Number of CPUs to use for multiprocessing the
      data reduction (sometimes not used)
    verbose : int (0,1,2)
      Level of verbosity:
        0 = No output
        1 = Minimal output (default - suitable for the average user)
        2 = All output
    use_masters : bool, optional
      Load calibration files from MasterFrames directory, if they exist
    logname : string
      The name of an ascii log file which is used to
      save the output details of the reduction
    ---------------------------------------------------
    """
    # Init logger
    import armsgs
    msgs = armsgs.get_logger((logname, debug, last_updated, version, verbose))
    import arload

    # First send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
    sigsignal(SIGINT, msgs.signal_handler)

    # Ignore all warnings given by python
    resetwarnings()
    simplefilter("ignore")

    # Record the starting time
    tstart = time()

    # Load the input file
    parlines, datlines, spclines = arload.load_input(redname)

    # Initialize the arguments and flags
    argflag = arload.argflag_init()
    argflag['run']['ncpus'] = ncpus
    argflag['out']['verbose'] = verbose

    # Determine the name of the spectrograph
    specname = None
    for i in range(len(parlines)):
        parspl = parlines[i].split()
        if len(parspl) < 3:
            msgs.error("There appears to be a missing argument on the following input line" + msgs.newline() +
                       parlines[i])
        if (parspl[0] == 'run') and (parspl[1] == 'spectrograph'):
            specname = parspl[2]
            break
    if specname is None:
        msgs.error("Please specify the spectrograph settings to be used with the command" + msgs.newline() +
                   "run spectrograph <name>")
    msgs.info("Reducing data from the {0:s} spectrograph".format(specname))

    # Determine the type of reduction used for this spectrograph
    redtype = None
    # Get the software path
    prgn_spl = progname.split('/')
    fname = "/".join(prgn_spl[:-1])
    fname += '/settings.'+specname
    spl = open(fname, 'r').readlines()
    for i in range(len(spl)):
        parspl = spl[i].split()
        if len(parspl) < 3:
            continue
        if (parspl[0] == 'mosaic') and (parspl[1] == 'reduction'):
            redtype = parspl[2]
            break
    if redtype is None:
        msgs.bug("The {0:s} instrument settings file must contain the reduction type".format(specname))
        msgs.error("Please specify the reduction type with the command" + msgs.newline() +
                   "mosaic reduction <type>")

    # Load default reduction arguments/flags, and set any command line arguments
    arparse.get_argflag((redtype.upper(), ".".join(redname.split(".")[:-1])))
    assert(False)

    # Load the default settings
    prgn_spl = progname.split('/')
    tfname = ""
    for i in range(0,len(prgn_spl)-2): tfname += prgn_spl[i]+"/"
    fname = tfname + prgn_spl[-2] + '/settings.' + spect['mosaic']['reduction'].lower()
    argflag = arload.load_settings(fname, argflag)
    argflag['run']['prognm'] = progname
    argflag['run']['pypitdir'] = tfname

    # Now update the settings based on the user input file
    argflag = arload.set_params(parlines, argflag, setstr="Input ")

    # Load the Spectrograph settings
    spect = arload.load_spect(progname, specname)

    # Load any changes to the spectrograph settings based on the user input file
    spect = arload.load_spect(progname, specname, spect=spect, lines=spclines)

    # Command line arguments
    if use_masters:
        argflag['masters']['use'] = True

    # If a quick reduction has been requested, make sure the requested pipeline
    # is the quick implementation (if it exists), otherwise run the standard pipeline.
    if quick:
        # Change to a "quick" settings file
        msgs.work("QUICK REDUCTION TO STILL BE DONE")

    # Load the important information from the fits headers
    fitsdict = arload.load_headers(argflag, spect, datlines)

    # Reduce the data!
    status = 0
    msgs.work("Make appropriate changes to quick reduction")
    if quick:
        msgs.work("define what is needed here for quick reduction")
    # Send the data away to be reduced
    if spect['mosaic']['reduction'] == 'ARMLSD':
        msgs.info("Data reduction will be performed using PYPIT-ARMLSD")
        import armlsd
        status = armlsd.ARMLSD(argflag, spect, fitsdict)
    elif spect['mosaic']['reduction'] == 'ARMED':
        msgs.info("Data reduction will be performed using PYPIT-ARMED")
        import armed
        status = armed.ARMED(argflag, spect, fitsdict)
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
    # Initiate logging for bugs and comand line help
    # These messages will not be saved to a log file
    from armsgs import Messages as Initmsg
    initmsgs = Initmsg(None, debug, last_updated, version, 1)
    # Set the default variables
    red = "script.pypit"
    qck = False
    cpu = 1
    vrb = 2
    use_masters = False

    if len(sys.argv) < 2:
        initmsgs.usage(None)

    # Load options from command line
    try:
        opt, arg = getopt.getopt(sys.argv[1:], 'hmqc:v:', ['help',
                                                           'use_masters',
                                                          'quick',
                                                          'cpus',
                                                          'verbose'])
        for o, a in opt:
            if o in ('-h', '--help'):
                initmsgs.usage(None)
            elif o in ('-q', '--quick'):
                qck = True
            elif o in ('-c', '--cpus'):
                cpu = int(a)
            elif o in ('-v', '--verbose'):
                vrb = int(a)
            elif o in ('-m', '--use_masters'):
                use_masters=True
        splitnm = os.path.splitext(arg[0])
        if splitnm[1] != '.pypit':
            initmsgs.error("Bad extension for PYPIT reduction file."+initmsgs.newline()+".pypit is required")
        lnm = splitnm[0] + ".log"
        red = arg[0]
    except getopt.GetoptError, err:
        initmsgs.error(err.msg, usage=True)

    # Execute the reduction, and catch any bugs for printout
    if debug['develop']:
        PYPIT(red, progname=sys.argv[0], quick=qck, ncpus=cpu, verbose=vrb,
              logname=lnm, use_masters=use_masters)
    else:
        try:
            PYPIT(red, progname=sys.argv[0], quick=qck, ncpus=cpu, verbose=vrb,
                  logname=lnm, use_masters=use_masters)
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
            from armsgs import get_logger
            get_logger().close()
