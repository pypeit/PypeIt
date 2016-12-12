from __future__ import absolute_import, division, print_function

from signal import SIGINT, signal as sigsignal
from warnings import resetwarnings, simplefilter
from time import time
from pypit import armsgs
from pypit import archeck
import glob
import numpy as np

# Import PYPIT routines

try:
    from linetools.spectra.xspectrum1d import XSpectrum1D
except ImportError:
    pass

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger


def PYPIT(redname, debug=None, progname=__file__, quick=False, ncpus=1, verbosity=1,
          use_masters=False, logname=None):
    """
    Main driver of the PYPIT code. Default settings and
    user-specified changes are made, and passed to the
    appropriate code for data reduction.

    Parameters
    ----------
    redname : string
      Input reduction script
    debug : dict, optional
      Debug dict
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
    verbosity : int (0,1,2)
      Level of verbosity:
        0 = No output
        1 = Minimal output (default - suitable for the average user)
        2 = All output
    use_masters : bool, optional
      Load calibration files from MasterFrames directory, if they exist
    logname : str or None
          The name of an ascii log file which is used to
          save the output details of the reduction
        debug : dict
          A PYPIT debug dict (from ardebug.init)
        version : str
        last_updated : str
    ---------------------------------------------------
    """
    from pypit import ardebug
    # Init logger
    if debug is None:
        debug = ardebug.init()
    msgs = armsgs.get_logger((logname, debug, verbosity))

    # This needs to be loaded after msgs
    from pypit import arparse

    # version checking
    try:
        archeck.version_check()
    except archeck.VersionError as err:
        msgs.error(err.message)
        
    # First send all signals to messages to be dealt with (i.e. someone hits ctrl+c)
    sigsignal(SIGINT, msgs.signal_handler)

    # Ignore all warnings given by python
    resetwarnings()
    simplefilter("ignore")

    # Record the starting time
    tstart = time()

    # Load the input file
    pyp_dict = load_input(redname, msgs)
    parlines, datlines, spclines = [pyp_dict[ii] for ii in ['par','dat','spc']]

    # Initialize the arguments and flags
#    argflag = arload.argflag_init()
#    settings.argflag['run']['ncpus'] = ncpus
#    settings.argflag['output']['verbosity'] = verbosity

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
    tfname = "/".join(prgn_spl[:-1]) + "/"
    fname = tfname + 'settings/settings.' + specname
    try:
        spl = open(fname, 'r').readlines()
    except IOError:
        msgs.error("The following instrument settings file cannot be found:" + msgs.newline() + fname + msgs.newline() +
                   "Please check the settings file exists, and that the instrument name is spelt correctly.")
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
    argf = arparse.get_argflag_class((redtype.upper(), ".".join(redname.split(".")[:-1])))
    lines = argf.load_file()
    argf.set_param('run pypitdir {0:s}'.format(tfname))
    argf.set_param('run progname {0:s}'.format(progname))
    argf.set_param('run redname {0:s}'.format(redname))
    argf.set_paramlist(lines)
    # Load user changes to the arguments/flags
    plines = argf.load_lines(parlines)
    argf.set_paramlist(plines)
    # If the user wishes to load a settings file, do that now
    if argf.__dict__['_argflag']['run']['load']['settings'] is not None:
        lines = argf.load_file(argf.__dict__['_argflag']['run']['load']['settings'])
        argf.set_paramlist(lines)

    # Load default spectrograph settings
    spect = arparse.get_spect_class((redtype.upper(), specname, ".".join(redname.split(".")[:-1])))
    lines = spect.load_file()
    spect.set_paramlist(lines)
    # Load user changes to the arguments/flags
    plines = spect.load_lines(spclines)
    spect.set_paramlist(plines)
    if argf.__dict__['_argflag']['run']['load']['spect'] is not None:
        lines = spect.load_file(argf.__dict__['_argflag']['run']['load']['spect'])
        spect.set_paramlist(lines)
    # If the instrument settings file sets some argflag settings, implement those changes now
    if len(spect.__dict__['_settings']) != 0:
        argf.set_paramlist(spect.__dict__['_settings'])
    # Load command line changes
    argf.set_param('run ncpus {0:d}'.format(ncpus))
    argf.set_param('output verbosity {0:d}'.format(verbosity))
    if use_masters:
        argf.set_param('reduce masters reuse True')
    msgs.work("Make appropriate changes to quick reduction")
    if quick:
        # If a quick reduction has been requested, make sure the requested pipeline
        # is the quick implementation (if it exists), otherwise run the standard pipeline.
        msgs.work("QUICK REDUCTION TO STILL BE DONE")
    # Finally, save the arguments/flags and spectrograph settings used for this reduction
    argf.save()
    spect.save()

    # Now that all of the relevant settings are loaded, globalize the settings
    arparse.init(argf, spect)

    '''
    # Test that a maximum of one .setup files is present
    from pypit import arsort
    setup_file, nexist = arsort.get_setup_file()
    if nexist == 1:
        msgs.info("Found setup_file: {:s}".format(setup_file))
        msgs.info("Will use this to guide the data reduction.")
    '''

    # Load the important information from the fits headers
    from pypit import arload
    fitsdict = arload.load_headers(datlines)

    # If the dispersion direction is 1, flip the axes
    if arparse.argflag['trace']['dispersion']['direction'] == 1:
        # Update the keywords of all fits files
        for ff in range(len(fitsdict['naxis0'])):
            temp = fitsdict['naxis0'][ff]
            fitsdict['naxis0'][ff] = fitsdict['naxis1'][ff]
            fitsdict['naxis1'][ff] = temp
        # Update the spectrograph settings for all detectors in the mosaic
        for dd in range(arparse.spect['mosaic']['ndet']):
            ddnum = arparse.get_dnum(dd+1)
            # Change the user-specified (x,y) pixel sizes
            tmp = arparse.spect[ddnum]['xgap']
            arparse.spect[ddnum]['xgap'] = arparse.spect[ddnum]['ygap']
            arparse.spect[ddnum]['ygap'] = tmp
            arparse.spect[ddnum]['ysize'] = 1.0 / arparse.spect[ddnum]['ysize']
            # Update the amplifier/data/overscan sections
            for i in range(arparse.spect[ddnum]['numamplifiers']):
                # Flip the order of the sections
                arparse.spect[ddnum]['datasec{0:02d}'.format(i + 1)] = arparse.spect[ddnum][
                                                                            'datasec{0:02d}'.format(i + 1)][::-1]
                arparse.spect[ddnum]['oscansec{0:02d}'.format(i + 1)] = arparse.spect[ddnum][
                                                                             'oscansec{0:02d}'.format(i + 1)][::-1]
    # Reduce the data!
    status = 0
    # Send the data away to be reduced
    if spect.__dict__['_spect']['mosaic']['reduction'] == 'ARMLSD':
        msgs.info("Data reduction will be performed using PYPIT-ARMLSD")
        from pypit import armlsd
        status = armlsd.ARMLSD(fitsdict)
    elif spect.__dict__['_spect']['mosaic']['reduction'] == 'ARMED':
        msgs.info("Data reduction will be performed using PYPIT-ARMED")
        from pypit import armed
        status = armed.ARMED(fitsdict)
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


def load_input(redname, msgs):
    """
    Load user defined input .pypit reduction file. Updates are
    made to the argflag dictionary.

    Parameters
    ----------
    redname : string
      Name of reduction script
    msgs : Messages
      logger for PYPIT

    Returns
    -------
    pyp_dict : dict
      Contains the following keys --
      'par'
        parlines : list
          Input (uncommented) lines specified by the user.
          parlines is used in this routine to update the
          argflag dictionary
      'dat'
        datlines : list
          Input (uncommented) lines specified by the user.
          datlines contains the full data path to every
          raw exposure listed by the user
      'spc'
        spclines : list
          Input (uncommented) lines specified by the user.
          spclines contains a list of user-specified changes
          that should be made to the default spectrograph
          settings.
      'dfn'
        dfnames : list
          Input data lines
      'setup'
        dict of setup info
          'name' list of setups
          'lines' list of lines in the setup block
      'ftype'
        dict of filename: frametype
    """
    import os
    # Read in the model file
    msgs.info("Loading the input file")
    try:
        infile = open(redname, 'r')
    except IOError:
        msgs.error("The filename does not exist -" + msgs.newline() + redname)
    lines = infile.readlines()
    parlines = []
    datlines = []
    skip_files = []
    spclines = []
    dfnames = []
    setuplines = []
    paths = []
    rddata, rdspec, rdsetup, rdsfiles = 0, 0, 0, -1
    setups = []
    ftype_dict = {}
    ftype_col = -1
    for i in range(len(lines)):
        if lines[i].strip() == '': continue
        linspl = lines[i].split()
        if rddata == 1: # Read datafile(s)
            if linspl[0] == 'data' and linspl[1] == 'end':
                rddata += 1
                # Deal with skip files
                if len(skip_files) > 0:
                    keep = np.array([True]*len(datlines))
                    for skip_file in skip_files:
                        for kk,datfile in enumerate(datlines):
                            if skip_file in datfile:
                                keep[kk] = False
                                msgs.warn("Skipping file {:s}".format(skip_file))
                    # Save
                    datlines = np.array(datlines)[keep].tolist()
                continue
            #
            dfname = lines[i].rstrip('\n').strip()
            if rdsfiles == -1:
                if 'path' in dfname[0:5]:
                    rdsfiles = 1
                else:
                    rdsfiles = 0
            if rdsfiles == 0:
                # is there a comment?
                aux = dfname.split('#')
                if len(aux) > 1:  # yes, there is a comment
                    dfname = aux[0].strip()
                if len(dfname) == 0:  # line is fully commented out
                    continue
                elif dfname[0] == '~':
                    dfname = os.path.expanduser(dfname)
                    print(dfname)
                elif dfname[:4] == 'skip':
                    skip_files.append(dfname.split(' ')[1])
                elif dfname[0] != '/':
                    msgs.error("You must specify the full datapath for the file:" + msgs.newline() + dfname)
                elif len(dfname.split()) != 1:
                    msgs.error("There must be no spaces when specifying the datafile:" + msgs.newline() + dfname)
                dfnames.append(dfname)
                listing = glob.glob(dfname)
                for lst in listing: datlines.append(lst)
            else:
                if 'path' in dfname[0:5]:
                    paths.append(linspl[1])
                else:  # Grab filename and frametype
                    if ftype_col == -1:  # Identify columns for frametype
                        ftype_col = np.where(np.array(linspl) == 'frametype')[0]
                        dfile_col = np.where(np.array(linspl) == 'filename')[0]
                    else:
                        # Find datafile and update ftype dict
                        for path in paths:
                            if os.path.isfile(path+linspl[dfile_col]):
                                datlines.append(path+linspl[dfile_col])
                                ftype_dict[linspl[dfile_col]] = linspl[ftype_col]
            continue
        elif rddata == 0 and linspl[0] == 'data' and linspl[1] == 'read': # Begin data read block
            rddata += 1
            continue
        if rdsetup == 1:  # Read setup command
            if linspl[0] == 'setup' and linspl[1] == 'end':
                rdsetup += 1
                continue
            if 'Setup' in lines[i]:
                setups.append(lines[i][6:].strip())
            setuplines.append(lines[i])
            continue
        elif rdsetup == 0 and linspl[0] == 'setup' and linspl[1] == 'read':  # Begin setup read block
            rdsetup += 1
            continue
        if rdspec == 1:  # Read spect command
            if linspl[0] == 'spect' and linspl[1] == 'end':
                rdspec += 1
                continue
            spclines.append(lines[i])
            continue
        elif rdspec == 0 and linspl[0] == 'spect' and linspl[1] == 'read':  # Begin spect read block
            rdspec += 1
            continue
        if lines[i].lstrip()[0] == '#': continue
        parlines.append(lines[i])
    # Do some quick checks
    if rddata == 0:
        msgs.error("You haven't specified any data!")
    elif rddata == 1:
        msgs.error("Missing 'data end' in " + redname)
    if rddata == 0:
        msgs.info("Using Default spectrograph parameters")
    elif rddata != 2:
        msgs.error("Missing 'spect end' in " + redname)
    # Check there are no duplicate inputs
    if len(datlines) != len(set(datlines)):
        msgs.error("There are duplicate files in the list of data.")
    if len(datlines) == 0:
        msgs.error("There are no raw data frames" + msgs.newline() +
                   "Perhaps the path to the data is incorrect?")
    else:
        msgs.info("Found {0:d} raw data frames".format(len(datlines)))
    msgs.info("Input file loaded successfully")
    # Let's return a dict
    pypit_dict = dict(par=parlines, dat=datlines, spc=spclines,
                      dfn=dfnames, setup={'name': setups, 'lines': setuplines},
                    ftype=ftype_dict)
    return pypit_dict # parlines, datlines, spclines, dfnames, setup, setuplines, ftype_dict



