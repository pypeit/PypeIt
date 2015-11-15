import numpy as np
# Import PYPIT routines
import armsgs as msgs

class ScienceExposure:

    def __init__(self, scidx, spect):
        """
        A Science Exposure class that carries all information for a given science exposure
        """
        #############################
        # Set some universal parameters
        self._argflag = argflag   # Arguments and Flags
        self._transpose = False   # Determine if the frames need to be transposed
        # Arrays to store the name for the frames that have already been combined
        self._done_bias, self._name_bias = [], []
        self._done_flat, self._name_flat = [], []
        self._done_arcs, self._name_arcs = [], []

        # Set indices used for frame combination
        self._idx_sci = scidx
        self._idx_arcs = None
        self._idx_bias = None
        self._idx_flat = None
        # Initialize Variables
        self._dispaxis = None
        self._pixlocn = None
        self._lordloc = None
        self._rordloc = None
        self._pixcen  = None
        self._pixwid  = None
        self._lordpix = None
        self._rordpix = None
        self._tilts   = None
        self._satmask = None
        self._arcparam = None
        self._wvcalib = None
        # Initialize the Master Calibration frames
        self._bpix   = None    # Bad Pixel Mask
        self._msarc = None     # Master Arc
        self._msbias = None    # Master Bias
        self._mstrace = None   # Master Trace
        self._mspixflat = None
        self._mspixflatnrm = None
        self._msblaze = None
        # Initialize the Master Calibration frame names
        self._msarc_name = None   # Master Arc Name
        self._msbias_name = None  # Master Bias Name
        self._mstrace_name = None  # Master Trace Name
        self._mspixflat_name = None  # Master Pixel Flat Name
        # Initialize the science, variance, and background frames
        self._sciframe = None
        self._varframe = None
        self._bgframe  = None
        # Initialize some extraction products
        self._ext_boxcar = None
        self._ext_optimal = None


    # Setters
    def SetMasterFrame(self, frame, ftype):
        if ftype == "arc": self._msarc = frame.copy()
        elif ftype == "badpix": self._bpix = frame.copy()
        elif ftype == "bias": self._msbias = frame.copy()
        elif ftype == "normpixflat": self._mspixflatnrm = frame.copy()
        elif ftype == "pixflat": self._mspixflat = frame.copy()
        elif ftype == "trace": self._mstrace = frame.copy()
        else:
            msgs.bug("I could not set master frame of type: {0:s}".format(ftype))
            msgs.error("Please contact the authors")
        return

    # Getters
    def GetMasterFrame(self, ftype):
        if ftype == "arc": return self._msarc.copy()
        elif ftype == "badpix": return self._bpix.copy()
        elif ftype == "bias": return self._msbias.copy()
        elif ftype == "normpixflat": return self._mspixflatnrm.copy()
        elif ftype == "pixflat": return self._mspixflat.copy()
        elif ftype == "trace": return self._mstrace.copy()
        else:
            msgs.bug("I could not set master frame of type: {0:s}".format(ftype))
        return None

