
class MasterFrames:

    def __init__(self, ndet):
        """
        A Master Calibrations class that carries the information generated for the master calibration
        """

        self._nspec    = [None for all in xrange(ndet)]   # Number of spectral pixels
        self._nspat    = [None for all in xrange(ndet)]   # Number of spatial pixels
        self._ampsec   = [None for all in xrange(ndet)]   # Locations of the amplifiers on each detector
        self._pixlocn  = [None for all in xrange(ndet)]   # Physical locations of each pixel on the detector
        self._lordloc  = [None for all in xrange(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc  = [None for all in xrange(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen   = [None for all in xrange(ndet)]   # Central slit traces in apparent pixel coordinates
        self._pixwid   = [None for all in xrange(ndet)]   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix  = [None for all in xrange(ndet)]   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix  = [None for all in xrange(ndet)]   # Array of slit traces (right side) in apparent pixel coordinates
        self._tilts    = [None for all in xrange(ndet)]   # Array of spectral tilts at each position on the detector
        self._satmask  = [None for all in xrange(ndet)]   # Array of Arc saturation streaks
        self._arcparam = [None for all in xrange(ndet)]   #
        self._wvcalib  = [None for all in xrange(ndet)]   #
        self._resnarr  = [None for all in xrange(ndet)]   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = [None for all in xrange(ndet)]          # Bad Pixel Mask
        self._msarc = [None for all in xrange(ndet)]         # Master Arc
        self._msbias = [None for all in xrange(ndet)]        # Master Bias
        self._mstrace = [None for all in xrange(ndet)]       # Master Trace
        self._mspixflat = [None for all in xrange(ndet)]     # Master pixel flat
        self._mspixflatnrm = [None for all in xrange(ndet)]  # Normalized Master pixel flat
        self._msblaze = [None for all in xrange(ndet)]       # Blaze function
        # Initialize the Master Calibration frame names
        self._msarc_name = [None for all in xrange(ndet)]      # Master Arc Name
        self._msbias_name = [None for all in xrange(ndet)]     # Master Bias Name
        self._mstrace_name = [None for all in xrange(ndet)]    # Master Trace Name
        self._mspixflat_name = [None for all in xrange(ndet)]  # Master Pixel Flat Name
