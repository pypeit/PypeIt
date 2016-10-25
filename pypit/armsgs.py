from __future__ import absolute_import, division, print_function

import sys
from os.path import dirname, basename
from textwrap import wrap as wraptext
from inspect import currentframe, getouterframes
from glob import glob

pypit_logger = None

class Messages:
    """
    Create coloured text for messages printed to screen.

    For further details on colours see the following example:
    http://ascii-table.com/ansi-escape-sequences.php
    """

    def __init__(self, log, debug, verbosity, colors=True):
        """
        Initialize the Message logging class

        Parameters
        ----------
        log : str or None
          Name of saved log file (no log will be saved if log=="")
        debug : dict
          dict used for debugging.
          'LOAD', 'BIAS', 'ARC', 'TRACE'
        verbosity : int (0,1,2)
          Level of verbosity:
            0 = No output
            1 = Minimal output (default - suitable for the average user)
            2 = All output
        colors : bool
          If true, the screen output will have colors, otherwise
          normal screen output will be displayed
        """
        # Version
        from pypit import pyputils
        version, last_updated = pyputils.get_version()
        # Import for version
        import scipy
        import numpy
        import astropy

        # Initialize the log
        if log is not None:
            self._log = open(log, 'w')
        else:
            self._log = log
        # Initialize other variables
        self._debug = debug
        self._last_updated = last_updated
        self._version = version
        self._verbosity = verbosity
        self.sciexp = None
        # Save the version of the code including last update information to the log file
        if self._log:
            self._log.write("------------------------------------------------------\n\n")
            self._log.write("PYPIT was last updated {0:s}\n".format(last_updated))
            self._log.write("This log was generated with version {0:s} of PYPIT\n\n".format(version))
            self._log.write("You are using scipy version={:s}\n".format(scipy.__version__))
            self._log.write("You are using numpy version={:s}\n".format(numpy.__version__))
            self._log.write("You are using astropy version={:s}\n\n".format(astropy.__version__))
            self._log.write("------------------------------------------------------\n\n")
        # Use colors?
        self._start, self._end = "", ""
        self._black_CL, self._yellow_CL, self._blue_CL, self._green_CL, self._red_CL = "", "", "", "", ""
        self._white_RD, self._white_GR, self._white_BK = "", "", ""
        self._white_BL, self._black_YL, self._yellow_BK = "", "", ""
        if colors:
            self.enablecolors()
        else:
            self.disablecolors()

    # Headers and usage
    def armedheader(self, prognm):
        """
        Get the info header for ARMED
        """
        header = "##  "
        header += self._start + self._white_GR + "ARMED : "
        header += "Automated Reduction and Modelling of Echelle Data v{0:s}".format(self._version) + self._end + "\n"
        header += "##  "
        header += "Usage : "
        header += "python %s [options] filelist".format(prognm)
        return header

    def pypitheader(self, prognm):
        """
        Get the info header for PYPIT
        """
        header = "##  "
        header += self._start + self._white_GR + "PYPIT : "
        header += "The Python Spectroscopic Data Reduction Pipeline v{0:s}".format(self._version) + self._end + "\n"
        header += "##  "
        #header += "Usage : "
        #if prognm is None:
        #    header += "pypit [options] filename.red"
        #else:
        #    header += "python %s [options] filename.red".format(prognm)
        return header

    def usage(self, prognm):
        stgs_arm = glob(dirname(__file__)+"/settings/settings.arm*")
        stgs_all = glob(dirname(__file__)+"/settings/settings.*")
        stgs_spc = list(set(stgs_arm) ^ set(stgs_all))
        armlist = basename(stgs_arm[0]).split(".")[-1]
        for i in range(1, len(stgs_arm)):
            armlist += ", " + basename(stgs_arm[i]).split(".")[-1]
        spclist = basename(stgs_spc[0]).split(".")[-1]
        for i in range(1, len(stgs_spc)):
            spclist += ", " + basename(stgs_spc[i]).split(".")[-1]
        spcl = wraptext(spclist, width=60)
        #print("\n#################################################################")
        #print(self.pypitheader(prognm))
        descs = self.pypitheader(prognm)
        #print("##  -------------------------------------------------------------")
        #print("##  Options: (default values in brackets)")
        #print("##   -c or --cpus      : (all) Number of cpu cores to use")
        #print("##   -h or --help      : Print this message")
        #print("##   -v or --verbosity : (2) Level of verbosity (0-2)")
        #print("##   -m or --use_masters : Use files in MasterFrames for reduction")
        #print("##   -d or --develop   : Turn develop debugging on")
        #print("##  -------------------------------------------------------------")
        descs += "\n##  Available pipelines include:"
        #print("##  Available pipelines include:")
        descs += "\n##   " + armlist
        #print("##  Available spectrographs include:")
        descs += "\n##  Available spectrographs include:"
        for i in spcl:
            descs += "\n##   " + i
            #print("##   " + i)
        #print("##  -------------------------------------------------------------")
        #print("##  Last updated: {0:s}".format(self._last_updated))
        descs += "\n##  Last updated: {0:s}".format(self._last_updated)
        #print("#################################################################\n")
        #sys.exit()
        return descs

    def debugmessage(self):
        if self._debug['develop']:
            info = getouterframes(currentframe())[2]
            dbgmsg = self._start+self._blue_CL+info[1].split("/")[-1]+" "+str(info[2])+" "+info[3]+"()"+self._end+" - "
        else:
            dbgmsg = ""
        return dbgmsg

    def close(self):
        """
        Close the log file and QA PDFs before the code exits
        """
        # Close PDFs
        try:
            self.sciexp._qa.close()
        except AttributeError:
            pass
        else:
            if self._debug['develop']:
                from pypit import armasters
                armasters.save_masters(self.sciexp, self.sciexp.det,
                                       self.sciexp._argflag['reduce']['masters']['setup'])
        # Close log
        if self._log:
            self._log.close()

    def signal_handler(self, signalnum, handler):
        """
        Handle signals sent by the keyboard during code execution
        """
        if signalnum == 2:
            self.info("Ctrl+C was pressed. Ending processes...")
            self.close()
            sys.exit()
        return

    def error(self, msg, usage=False):
        """
        Print an error message
        """
        dbgmsg = self.debugmessage()
        premsg = "\n"+self._start + self._white_RD + "[ERROR]   ::" + self._end + " "
        print(premsg+dbgmsg+msg, file=sys.stderr)
        if self._log:
            self._log.write(self.cleancolors(premsg+dbgmsg+msg)+"\n")
        # Close PDFs and log file
        self.close()
        # Print command line usage
        if usage:
            self.usage(None)
        sys.exit()

    def info(self, msg):
        """
        Print an information message
        """
        dbgmsg = self.debugmessage()
        premsg = self._start + self._green_CL + "[INFO]    ::" + self._end + " "
        print(premsg+dbgmsg+msg, file=sys.stderr)
        if self._log:
            self._log.write(self.cleancolors(premsg+dbgmsg+msg)+"\n")
        return

    def info_update(self, msg, last=False):
        """
        Print an information message that needs to be updated
        """
        dbgmsg = self.debugmessage()
        premsg = "\r" + self._start + self._green_CL + "[INFO]    ::" + self._end + " "
        if last:
            print(premsg+dbgmsg+msg, file=sys.stderr)
            if self._log:
                self._log.write(self.cleancolors(premsg+dbgmsg+msg)+"\n")
        else:
            print(premsg+dbgmsg+msg, file=sys.stderr)
            if self._log:
                self._log.write(self.cleancolors(premsg+dbgmsg+msg))
        return

    def test(self, msg):
        """
        Print a test message
        """
        if self._verbosity == 2:
            dbgmsg = self.debugmessage()
            premsg = self._start + self._white_BL + "[TEST]    ::" + self._end + " "
            print(premsg+dbgmsg+msg, file=sys.stderr)
            if self._log:
                self._log.write(self.cleancolors(premsg+dbgmsg+msg)+"\n")
        return

    def warn(self, msg):
        """
        Print a warning message
        """
        dbgmsg = self.debugmessage()
        premsg = self._start + self._red_CL + "[WARNING] ::" + self._end + " "
        print(premsg+dbgmsg+msg, file=sys.stderr)
        if self._log:
            self._log.write(self.cleancolors(premsg+dbgmsg+msg)+"\n")
        return

    def bug(self, msg):
        """
        Print a bug message
        """
        dbgmsg = self.debugmessage()
        premsg = self._start + self._white_BK + "[BUG]     ::" + self._end + " "
        print(premsg+dbgmsg+msg, file=sys.stderr)
        if self._log:
            self._log.write(self.cleancolors(premsg+dbgmsg+msg))
        return

    def work(self, msg):
        """
        Print a work in progress message
        """
        if self._verbosity == 2:
            dbgmsg = self.debugmessage()
            premsgp = self._start + self._black_CL + "[WORK IN ]::" + self._end + "\n"
            premsgs = self._start + self._yellow_CL + "[PROGRESS]::" + self._end + " "
            print(premsgp+premsgs+dbgmsg+msg, file=sys.stderr)
            if self._log:
                self._log.write(self.cleancolors(premsgp+premsgs+dbgmsg+msg)+"\n")
        return

    def prindent(self, msg):
        """
        Print an indent
        """
        premsg = "             "
        print(premsg+msg, file=sys.stderr)
        if self._log:
            self._log.write(self.cleancolors(premsg+msg)+"\n")
        return

    def input(self):
        """
        Return a text string to be used to display input required from the user
        """
        premsg = self._start + self._blue_CL + "[INPUT]   ::" + self._end + " "
        return premsg

    @staticmethod
    def newline():
        """
        Return a text string containing a newline to be used with messages
        """
        return "\n             "

    @staticmethod
    def indent():
        """
        Return a text string containing an indent to be used with messages
        """
        return "             "

    # Set the colors
    def enablecolors(self):
        """
        Enable colored output text
        """

        # Start and end coloured text
        self._start = "\x1B["
        self._end = "\x1B[" + "0m"

        # Clear Backgrounds
        self._black_CL = "1;30m"
        self._yellow_CL = "1;33m"
        self._blue_CL = "1;34m"
        self._green_CL = "1;32m"
        self._red_CL = "1;31m"

        # Coloured Backgrounds
        self._white_RD = "1;37;41m"
        self._white_GR = "1;37;42m"
        self._white_BK = "1;37;40m"
        self._white_BL = "1;37;44m"
        self._black_YL = "1;37;43m"
        self._yellow_BK = "1;33;40m"

    def cleancolors(self, msg):
        cols = [self._end, self._start,
                self._black_CL, self._yellow_CL, self._blue_CL, self._green_CL, self._red_CL,
                self._white_RD, self._white_GR, self._white_BK, self._white_BL, self._black_YL, self._yellow_BK]
        for i in cols:
            msg = msg.replace(i, "")
        return msg

    def disablecolors(self):
        """
        Disable colored output text
        """

        # Start and end coloured text
        self._start = ""
        self._end = ""

        # Clear Backgrounds
        self._black_CL = ""
        self._yellow_CL = ""
        self._blue_CL = ""
        self._green_CL = ""
        self._red_CL = ""

        # Coloured Backgrounds
        self._white_RD = ""
        self._white_GR = ""
        self._white_BK = ""
        self._white_BL = ""
        self._black_YL = ""
        self._yellow_BK = ""


def get_logger(init=None):
    """ Logger
    Parameters
    ----------
    init : tuple
      For instantiation
      (log, debug, verbosity)

    Returns
    -------
    msgs : Messages
    """
    global pypit_logger

    # Instantiate??
    if init is not None:
        pypit_logger = Messages(init[0], init[1], init[2])

    return pypit_logger


