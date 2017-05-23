from __future__ import (print_function, absolute_import, division, unicode_literals)

from pypit import armsgs
from pypit import arparse as settings
from pypit import arsave

try:
    basestring
except NameError:  # For Python 3
    basestring = str

# Logging
msgs = armsgs.get_logger()

from pypit import ardebug as debugger

class MasterFrames:

    def __init__(self, ndet):
        """
        A Master Calibrations class that carries the information generated for the master calibration
        """

        self._nspec    = [None for all in range(ndet)]   # Number of spectral pixels
        self._nspat    = [None for all in range(ndet)]   # Number of spatial pixels
        self._ampsec   = [None for all in range(ndet)]   # Locations of the amplifiers on each detector
        self._pixlocn  = [None for all in range(ndet)]   # Physical locations of each pixel on the detector
        self._lordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._rordloc  = [None for all in range(ndet)]   # Array of slit traces (left side) in physical pixel coordinates
        self._pixcen   = [None for all in range(ndet)]   # Central slit traces in apparent pixel coordinates
        self._pixwid   = [None for all in range(ndet)]   # Width of slit (at each row) in apparent pixel coordinates
        self._lordpix  = [None for all in range(ndet)]   # Array of slit traces (left side) in apparent pixel coordinates
        self._rordpix  = [None for all in range(ndet)]   # Array of slit traces (right side) in apparent pixel coordinates
        self._tilts    = [None for all in range(ndet)]   # Array of spectral tilts at each position on the detector
        self._satmask  = [None for all in range(ndet)]   # Array of Arc saturation streaks
        self._arcparam = [None for all in range(ndet)]   #
        self._wvcalib  = [None for all in range(ndet)]   #
        self._resnarr  = [None for all in range(ndet)]   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = [None for all in range(ndet)]          # Bad Pixel Mask
        self._msarc = [None for all in range(ndet)]         # Master Arc
        self._msbias = [None for all in range(ndet)]        # Master Bias
        self._mstrace = [None for all in range(ndet)]       # Master Trace
        self._mspixelflat = [None for all in range(ndet)]     # Master pixel flat
        self._mspixelflatnrm = [None for all in range(ndet)]  # Normalized Master pixel flat
        self._msblaze = [None for all in range(ndet)]       # Blaze function
        # Initialize the Master Calibration frame names
        self._msarc_name = [None for all in range(ndet)]      # Master Arc Name
        self._msbias_name = [None for all in range(ndet)]     # Master Bias Name
        self._mstrace_name = [None for all in range(ndet)]    # Master Trace Name
        self._mspixelflat_name = [None for all in range(ndet)]  # Master Pixel Flat Name


def master_name(ftype, setup, mdir=None):
    """ Default filenames
    Parameters
    ----------
    ftype
    mdir : str, optional
      Master directory; usually taken from settings

    Returns
    -------
    """
    if mdir is None:
        mdir = settings.argflag['run']['directory']['master']+'_'+settings.argflag['run']['spectrograph']
    name_dict = dict(bias='{:s}/MasterBias_{:s}.fits'.format(mdir, setup),
                     badpix='{:s}/MasterBadPix_{:s}.fits'.format(mdir, setup),
                     trace='{:s}/MasterTrace_{:s}.fits'.format(mdir, setup),
                     normpixelflat='{:s}/MasterFlatField_{:s}.fits'.format(mdir, setup),
                     arc='{:s}/MasterArc_{:s}.fits'.format(mdir, setup),
                     wave='{:s}/MasterWave_{:s}.fits'.format(mdir, setup),
                     wave_calib='{:s}/MasterWaveCalib_{:s}.json'.format(mdir, setup),
                     tilts='{:s}/MasterTilts_{:s}.fits'.format(mdir, setup),
                     slitprof='{:s}/MasterSlitProfile_{:s}.fits'.format(mdir, setup),
                     )
    return name_dict[ftype]

'''
def load_masters(slf, det, setup):
    """ Load master frames
    Parameters
    ----------
    slf
    det
    setup
    Returns
    -------
    """
    def load_master(file, exten=0):
        hdu = pyfits.open(file)
        data = hdu[exten].data
        return data
    # Bias
    slf._msbias[det-1] = load_master(master_name('bias', setup))
    # Bad Pixel
    slf._bpix[det-1] = load_master(master_name('badpix', setup))
    # Trace
    slf._mstrace[det-1] = load_master(master_name('trace', setup))
    slf._pixcen[det-1] = load_master(master_name('trace', setup), exten=1)
    slf._pixwid[det-1] = load_master(master_name('trace', setup), exten=2)
    slf._lordpix[det-1] = load_master(master_name('trace', setup), exten=3)
    slf._rordpix[det-1] = load_master(master_name('trace', setup), exten=4)
    # Flat
    slf._mspixelflatnrm[det-1] = load_master(master_name('normpixelflat', setup))
    # Arc/wave
    slf._msarc[det-1] = load_master(master_name('arc', setup))
    slf._mswave[det-1] = load_master(master_name('wave', setup))
    # Tilts
    slf._tilts[det-1] = load_master(master_name('tilts', setup))
'''


def save_masters(slf, det, setup):
    """ Save Master Frames
    Parameters
    ----------
    slf
    setup
    Returns
    -------
    """
    from linetools import utils as ltu
    import io, json

    transpose = bool(settings.argflag['trace']['dispersion']['direction'])

    # Bias
    if 'bias'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        if not isinstance(slf._msbias[det-1], (basestring)):
            arsave.save_master(slf, slf._msbias[det-1],
                               filename=master_name('bias', setup),
                               frametype='bias')
    # Bad Pixel
    if 'badpix'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        arsave.save_master(slf, slf._bpix[det-1],
                               filename=master_name('badpix', setup),
                               frametype='badpix')
    # Trace
    if 'trace'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        extensions = [slf._lordloc[det-1], slf._rordloc[det-1],
                      slf._pixcen[det-1], slf._pixwid[det-1],
                      slf._lordpix[det-1], slf._rordpix[det-1],
                      slf._slitpix[det-1]]
        names = ['LeftEdges_det', 'RightEdges_det', 'SlitCentre', 'SlitLength', 'LeftEdges_pix', 'RightEdges_pix', 'SlitPixels']
        arsave.save_master(slf, slf._mstrace[det-1],
                           filename=master_name('trace', setup),
                           frametype='trace', extensions=extensions, names=names)
    # Pixel Flat
    if 'normpixelflat'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        arsave.save_master(slf, slf._mspixelflatnrm[det-1],
                           filename=master_name('normpixelflat', setup),
                           frametype='normpixelflat')
    # Arc/Wave
    if 'arc'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        arsave.save_master(slf, slf._msarc[det-1],
                           filename=master_name('arc', setup),
                           frametype='arc', keywds=dict(transp=transpose))
    if 'wave'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        # Wavelength image
        arsave.save_master(slf, slf._mswave[det-1],
                           filename=master_name('wave', setup),
                           frametype='wave')
        # Wavelength fit
        gddict = ltu.jsonify(slf._wvcalib[det-1])
        json_file = master_name('wave_calib', setup)
        if gddict is not None:
            ltu.savejson(json_file, gddict, easy_to_read=True, overwrite=True)
        else:
            msgs.warn("The master wavelength solution has not been saved")
    if 'tilts'+settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        arsave.save_master(slf, slf._tilts[det-1],
                           filename=master_name('tilts', setup),
                           frametype='tilts')

    # Spatial slit profile
    if 'slitprof' + settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        arsave.save_master(slf, slf._slitprof[det - 1],
                           filename=master_name('slitprof', setup),
                           frametype='slit profile')


def user_master_name(mdir, input_name):
    """ Convert user-input filename for master into full name
    Mainly used to append MasterFrame directory

    Parameters
    ----------
    mdir : str
    input_name : str

    Returns
    -------
    full_name : str

    """
    islash = input_name.find('/')
    if islash >= 0:
        full_name = input_name
    else:
        full_name = mdir+'/'+input_name
    # Return
    return full_name
