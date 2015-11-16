import numpy as np
# Import PYPIT routines
import armsgs as msgs
from matplotlib.backends.backend_pdf import PdfPages

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

        # Initialize the QA for this science exposure
        qafn = "testname.pdf"
        self._qa = pp = PdfPages(qafn)
        # Set indices used for frame combination
        self._idx_sci = scidx
        self._idx_arcs = None
        self._idx_bias = None
        self._idx_flat = None
        # Initialize Variables
        self._dispaxis = None  # Which direction is the predominant spectral (dispersion) axis
        self._pixlocn = None   # Physical locations of each pixel on the detector
        self._lordloc = None   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc = None   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen  = None   # Central slit traces in apparent pixel coordinates
        self._pixwid  = None   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix = None   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix = None   # Array of slit traces (right side) in apparent pixel coordinates
        self._tilts   = None   # Array of spectral tilts at each position on the detector
        self._satmask = None   # Array of Arc saturation streaks
        self._arcparam = None  #
        self._wvcalib = None   #
        self._resnarr = None   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = None          # Bad Pixel Mask
        self._msarc = None         # Master Arc
        self._msbias = None        # Master Bias
        self._mstrace = None       # Master Trace
        self._mspixflat = None     # Master pixel flat
        self._mspixflatnrm = None  # Normalized Master pixel flat
        self._msblaze = None       # Blaze function
        # Initialize the Master Calibration frame names
        self._msarc_name = None      # Master Arc Name
        self._msbias_name = None     # Master Bias Name
        self._mstrace_name = None    # Master Trace Name
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

