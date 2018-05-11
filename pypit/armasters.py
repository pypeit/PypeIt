from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import numpy as np
import yaml

import linetools.utils

from pypit import msgs
from pypit import arload
from pypit import arparse as settings
from pypit import arsave
from pypit import arutils
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
        self._wvcalib  = [None for all in range(ndet)]   # List of dicts
        self._resnarr  = [None for all in range(ndet)]   # Resolution array
        # Initialize the Master Calibration frames
        self._bpix = [None for all in range(ndet)]          # Bad Pixel Mask
        self._msarc = [None for all in range(ndet)]         # Master Arc
        self._msbias = [None for all in range(ndet)]        # Master Bias
        self._mstrace = [None for all in range(ndet)]       # Master Trace
        self._mspinhole = [None for all in range(ndet)]       # Master Pinhole
        self._mspixelflat = [None for all in range(ndet)]     # Master pixel flat
        self._mspixelflatnrm = [None for all in range(ndet)]  # Normalized Master pixel flat
        self._msblaze = [None for all in range(ndet)]       # Blaze function
        self._sensfunc = [None for all in range(ndet)]       # Sensitivity function
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
    setup : str
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
                     pinhole='{:s}/MasterPinhole_{:s}.fits'.format(mdir, setup),
                     normpixelflat='{:s}/MasterFlatField_{:s}.fits'.format(mdir, setup),
                     arc='{:s}/MasterArc_{:s}.fits'.format(mdir, setup),
                     wave='{:s}/MasterWave_{:s}.fits'.format(mdir, setup),
                     wv_calib='{:s}/MasterWaveCalib_{:s}.json'.format(mdir, setup),
                     tilts='{:s}/MasterTilts_{:s}.fits'.format(mdir, setup),
                     slitprof='{:s}/MasterSlitProfile_{:s}.fits'.format(mdir, setup),
                     sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.yaml'.format(mdir, setup[0], setup[-2:]),
                     )
    return name_dict[ftype]


def get_master_frame(slf, mftype, det=None):
    ret, head = core_get_master_frame(mftype, settings.argflag, slf.setup)
    if ret is None:
        return None
    elif mftype == 'arc':
        slf._transpose = head['transp']
        if slf._transpose:  # Need to setup for flipping
            settings.argflag['trace']['dispersion']['direction'] = 1
        else:
            settings.argflag['trace']['dispersion']['direction'] = 0
    elif mftype == 'trace':
        slf.SetFrame(slf._lordloc, ret['lordloc'], det)
        slf.SetFrame(slf._rordloc, ret['rordloc'], det)
        slf.SetFrame(slf._pixcen, ret['pixcen'].astype(np.int), det)
        slf.SetFrame(slf._pixwid, ret['pixwid'].astype(np.int), det)
        slf.SetFrame(slf._lordpix, ret['lordpix'].astype(np.int), det)
        slf.SetFrame(slf._rordpix, ret['rordpix'].astype(np.int), det)
        slf.SetFrame(slf._slitpix, ret['slitpix'].astype(np.int), det)
        # Mask -- It is assumed that all slits loaded are ok
        slf._maskslits[det-1] = np.array([False] * slf._lordloc[det-1].shape[1])
    return ret


def core_get_master_frame(mftype, settings, setup):
    """ If a MasterFrame exists, load it

    Parameters
    ----------
    mftype : str
    settings : dict
    setup : str

    Returns
    -------
    msfile : ndarray or dict or None
    head : Header or None

    """
    #setup = settings['reduce']['masters']['setup']
    # Were MasterFrames even desired?
    #  TODO -- Consider asking this if statement *before* calling get master_frame
    if (settings['reduce']['masters']['reuse']) or (settings['reduce']['masters']['force']):
        ms_name = master_name(mftype, setup)
        msfile, head = arload.load_master(ms_name, frametype=mftype)
        if msfile is None:
            msgs.warn("No Master frame found of type {:s}: {:s}".format(mftype,ms_name))
            return None, None
        else:  # Extras
            if mftype == 'trace':
                tdict = {}
                tdict['lordloc'], _ = arload.load_master(ms_name, frametype="trace", exten=1)
                tdict['rordloc'], _ = arload.load_master(ms_name, frametype="trace", exten=2)
                tdict['pixcen'], _ = arload.load_master(ms_name, frametype="trace", exten=3)
                tdict['pixwid'], _ = arload.load_master(ms_name, frametype="trace", exten=4)
                tdict['lordpix'], _ = arload.load_master(ms_name, frametype="trace", exten=5)
                tdict['rordpix'], _ = arload.load_master(ms_name, frametype="trace", exten=6)
                tdict['slitpix'], _ = arload.load_master(ms_name, frametype="trace", exten=7)
                msfile = tdict  # For returning
            # Append as loaded
            settings['reduce']['masters']['loaded'].append(mftype+setup)
            return msfile, head
    else:
        return None, None

def save_masters(slf, det, mftype='all'):
    """ Save Master Frames

    Parameters
    ----------
    slf
    mftype : str
      'all' -- Save them all

    """
    setup = slf.setup
    transpose = bool(settings.argflag['trace']['dispersion']['direction'])

    # Bias
    if (mftype == 'bias'):
        msgs.error("Should not get here anymore.  Save the bias in the BiasPrep class")
        #and ('bias'+setup not in settings.argflag['reduce']['masters']['loaded']):
        #arsave.core_save_master(object, filename=master_name('bias', setup),
        #                   frametype='bias', raw_files=raw_files)
    # Bad Pixel
    if (mftype in ['badpix', 'all']) and ('badpix'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._bpix[det-1],
                               filename=master_name('badpix', setup),
                               frametype='badpix')
    # Trace
    if (mftype in ['trace', 'all']) and ('trace'+setup not in settings.argflag['reduce']['masters']['loaded']):
        extensions = [slf._lordloc[det-1], slf._rordloc[det-1],
                      slf._pixcen[det-1], slf._pixwid[det-1],
                      slf._lordpix[det-1], slf._rordpix[det-1],
                      slf._slitpix[det-1]]
        names = ['LeftEdges_det', 'RightEdges_det', 'SlitCentre', 'SlitLength', 'LeftEdges_pix', 'RightEdges_pix', 'SlitPixels']
        arsave.save_master(slf, slf._mstrace[det-1],
                           filename=master_name('trace', setup),
                           frametype='trace', extensions=extensions, names=names)
    # Pixel Flat
    if (mftype in ['normpixelflat', 'all']) and ('normpixelflat'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._mspixelflatnrm[det-1],
                           filename=master_name('normpixelflat', setup),
                           frametype='normpixelflat')
    # Pinhole Flat
    if (mftype in ['pinhole', 'all']) and ('pinhole'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._mspinhole[det-1],
                           filename=master_name('pinhole', setup),
                           frametype='pinhole')
    # Arc/Wave
    if (mftype in ['arc', 'all']) and ('arc'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._msarc[det-1],
                           filename=master_name('arc', setup),
                           frametype='arc', keywds=dict(transp=transpose))
    if (mftype in ['wave', 'all']) and ('wave'+setup not in settings.argflag['reduce']['masters']['loaded']):
        # Wavelength image
        arsave.save_master(slf, slf._mswave[det-1],
                           filename=master_name('wave', setup),
                           frametype='wave')
    if (mftype in ['wv_calib', 'all']) and ('wv_calib'+setup not in settings.argflag['reduce']['masters']['loaded']):
        # Wavelength fit
        gddict = linetools.utils.jsonify(slf._wvcalib[det-1])
        json_file = master_name('wv_calib', setup)
        if gddict is not None:
            linetools.utils.savejson(json_file, gddict, easy_to_read=True, overwrite=True)
        else:
            msgs.warn("The master wavelength solution has not been saved")
    # Tilts
    if (mftype in ['tilts', 'all']) and ('tilts'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._tilts[det-1],
                           filename=master_name('tilts', setup),
                           frametype='tilts')
    # Spatial slit profile
    if (mftype in ['slitprof', 'all']) and ('slitprof'+setup not in settings.argflag['reduce']['masters']['loaded']):
        arsave.save_master(slf, slf._slitprof[det - 1],
                           filename=master_name('slitprof', setup),
                           frametype='slit profile')


def save_sensfunc(slf, setup):
    """ Make YAML friendly and write to disk
    Separate routine as this process is detector independent
    
    Parameters
    ----------
    slf
    setup : str
    """
    # Sensitivity Function
    if 'sensfunc' + settings.argflag['reduce']['masters']['setup'] not in settings.argflag['reduce']['masters']['loaded']:
        # yamlify
        ysens = arutils.yamlify(slf._sensfunc)
        with open(master_name('sensfunc', setup), 'w') as yamlf:
            yamlf.write(yaml.dump(ysens))


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
