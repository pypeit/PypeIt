""" Routines related to MasterFrames"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    basestring
except NameError:  # For Python 3
    basestring = str

import numpy as np
import os
import yaml

from astropy.io import fits
from astropy import units

import linetools.utils

from pypit import msgs
from pypit import arparse as settings
from pypit import arutils
from pypit import traceslits

from pypit import ardebug as debugger

'''
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
'''


def master_name(ftype, setup, mdir=None):
    if mdir is None:
        mdir = settings.argflag['run']['directory']['master']+'_'+settings.argflag['run']['spectrograph']
    return core_master_name(ftype, setup, mdir)


def core_master_name(ftype, setup, mdir):
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
    name_dict = dict(bias='{:s}/MasterBias_{:s}.fits'.format(mdir, setup),
                     badpix='{:s}/MasterBadPix_{:s}.fits'.format(mdir, setup),
                     trace='{:s}/MasterTrace_{:s}'.format(mdir, setup),   # Just a root as FITS+JSON are generated
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


def load_master_frame(slf, mftype, det=None):
    # Were MasterFrames even desired?
    if (settings.argflag['reduce']['masters']['reuse']) or (settings.argflag['reduce']['masters']['force']):
        ret, head, _ = core_load_master_frame(mftype, slf.setup,
                                           settings.argflag['run']['directory']['master']+'_'+settings['run']['spectrograph'],
                                           force=settings.argflag['reduce']['masters']['force'])
    else:
        return None, None
    if ret is None:
        return None
    elif mftype == 'arc':
        slf._transpose = head['transp']
        if slf._transpose:  # Need to setup for flipping
            settings.argflag['trace']['dispersion']['direction'] = 1
        else:
            settings.argflag['trace']['dispersion']['direction'] = 0
    elif mftype == 'trace':
        Tslits = ret[0]
        Tslits._make_pixel_arrays()
        #
        slf.SetFrame(slf._lordloc, Tslits.lcen, det)
        slf.SetFrame(slf._rordloc, Tslits.rcen, det)
        slf.SetFrame(slf._pixcen, Tslits.pixcen, det)
        slf.SetFrame(slf._pixwid, Tslits.pixwid, det)
        slf.SetFrame(slf._lordpix, Tslits.lordpix, det)
        slf.SetFrame(slf._rordpix, Tslits.rordpix, det)
        slf.SetFrame(slf._slitpix, Tslits.slitpix, det)
        # Mask -- It is assumed that all slits loaded are ok
        slf._maskslits[det-1] = np.array([False] * slf._lordloc[det-1].shape[1])
        # We only want to send back the mstrace image (for now)
        #     This should change when slf is Refactored
        ret = Tslits.mstrace
    # Append as loaded
    settings.argflag['reduce']['masters']['loaded'].append(mftype+slf.setup)
    return ret


def core_load_master_frame(mftype, setup, mdir, force=False):
    """ If a MasterFrame exists, load it

    Parameters
    ----------
    mftype : str
    setup : str
    mdir : str

    Returns
    -------
    msfile : ndarray or dict or None
    head : Header or None
    file_list : list (or None)
      Typically the files used the generate the master frame (may be incomplete or None)

    """
    # Name
    ms_name = core_master_name(mftype, setup, mdir)
    # Load
    msframe, head, file_list = _core_load(ms_name, exten=0, frametype=mftype, force=force)
    # Check
    if msframe is None:
        msgs.warn("No Master frame found of type {:s}: {:s}".format(mftype,ms_name))
        return None, None, None
    #else:  # Extras?
    #    if mftype == 'trace':
    #        tdict = {}
    #        tdict['mstrace'] = msframe0.copy()
    #        tdict['lordloc'], _, _ = _core_load(ms_name, frametype="trace", exten=1)
    #        tdict['rordloc'], _, _ = _core_load(ms_name, frametype="trace", exten=2)
    #        tdict['pixcen'], _, _ = _core_load(ms_name, frametype="trace", exten=3)
    #        tdict['pixwid'], _, _ = _core_load(ms_name, frametype="trace", exten=4)
    #        tdict['lordpix'], _, _ = _core_load(ms_name, frametype="trace", exten=5)
    #        tdict['rordpix'], _, _ = _core_load(ms_name, frametype="trace", exten=6)
    #        tdict['slitpix'], _, _ = _core_load(ms_name, frametype="trace", exten=7)
    #        msframe = tdict  # Just for returning
    #    else:
    #        msframe = msframe0
    # Return
    return msframe, head, file_list


def load_master(name, exten=0, frametype='<None>'):
    # TODO -Deprecate
    return _core_load(name, exten=exten, frametype=frametype,
                         force=settings.argflag['reduce']['masters']['force'])

def _core_load(name, exten=0, frametype='<None>', force=False):
    """
    Low level load method for master frames
      Should only be called by core_load_master_frame

    Parameters
    ----------
    name : str
      Name of the master calibration file to be loaded
    exten : int, optional
    frametype : str, optional
      The type of master calibration frame being loaded.
      This keyword is only used for terminal print out.

    Returns
    -------
    frame : ndarray or dict or TraceSlits or None
      The data from the master calibration frame
    head : str (or None)
    file_list : list (or None)
    """
    # Check to see if file exists
    if not os.path.isfile(name):
        msgs.warn("Master frame does not exist: {:s}".format(name))
        if force:
            msgs.error("Crashing out because reduce-masters-force=True:"+msgs.newline()+name)
        return None, None, None
    #
    if frametype == 'wv_calib':
        msgs.info("Loading Master {0:s} frame:".format(frametype)+msgs.newline()+name)
        ldict = linetools.utils.loadjson(name)
        return ldict, None, [name]
    elif frametype == 'sensfunc':
        with open(name, 'r') as f:
            sensfunc = yaml.load(f)
        sensfunc['wave_max'] = sensfunc['wave_max']*units.AA
        sensfunc['wave_min'] = sensfunc['wave_min']*units.AA
        return sensfunc, None, [name]
    elif frametype == 'trace':
        Tslits = traceslits.TraceSlits.from_master_files(name)
        return Tslits, None None
    else:
        msgs.info("Loading a pre-existing master calibration frame")
        hdu = fits.open(name)
        msgs.info("Master {0:s} frame loaded successfully:".format(hdu[0].header['FRAMETYP'])+msgs.newline()+name)
        head0 = hdu[0].header
        data = hdu[exten].data.astype(np.float)
        # List of files used to generate the Master frame (e.g. raw file frames)
        file_list = []
        for key in head0:
            if 'FRAME' in key:
                file_list.append(head0[key])
        return data, head0, file_list


def save_masters(slf, det, mftype='all'):
    """ Save Master Frames

    Parameters
    ----------
    slf
    mftype : str
      'all' -- Save them all

    """
    # TODO - Deprecate
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
        save_master(slf, slf._bpix[det-1],
                               filename=master_name('badpix', setup),
                               frametype='badpix')
    # Trace
    if (mftype in ['trace', 'all']) and ('trace'+setup not in settings.argflag['reduce']['masters']['loaded']):
        extensions = [slf._lordloc[det-1], slf._rordloc[det-1],
                      slf._pixcen[det-1], slf._pixwid[det-1],
                      slf._lordpix[det-1], slf._rordpix[det-1],
                      slf._slitpix[det-1]]
        names = ['LeftEdges_det', 'RightEdges_det', 'SlitCentre', 'SlitLength', 'LeftEdges_pix', 'RightEdges_pix', 'SlitPixels']
        save_master(slf, slf._mstrace[det-1],
                           filename=master_name('trace', setup),
                           frametype='trace', extensions=extensions, names=names)
    # Pixel Flat
    if (mftype in ['normpixelflat', 'all']) and ('normpixelflat'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._mspixelflatnrm[det-1],
                           filename=master_name('normpixelflat', setup),
                           frametype='normpixelflat')
    # Pinhole Flat
    if (mftype in ['pinhole', 'all']) and ('pinhole'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._mspinhole[det-1],
                           filename=master_name('pinhole', setup),
                           frametype='pinhole')
    # Arc/Wave
    if (mftype in ['arc', 'all']) and ('arc'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._msarc[det-1],
                           filename=master_name('arc', setup),
                           frametype='arc', keywds=dict(transp=transpose))
    if (mftype in ['wave', 'all']) and ('wave'+setup not in settings.argflag['reduce']['masters']['loaded']):
        # Wavelength image
        save_master(slf, slf._mswave[det-1],
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
        save_master(slf, slf._tilts[det-1],
                           filename=master_name('tilts', setup),
                           frametype='tilts')
    # Spatial slit profile
    if (mftype in ['slitprof', 'all']) and ('slitprof'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._slitprof[det - 1],
                           filename=master_name('slitprof', setup),
                           frametype='slit profile')


def save_master(slf, data, filename="temp.fits", frametype="<None>", ind=[],
                        extensions=None, keywds=None, names=None):
    """ Wrapper to core_save_master
    Will be Deprecated

    Parameters
    ----------
    slf
    data
    filename
    frametype
    ind
    extensions
    keywds
    names

    Returns
    -------

    """
    if len(ind) > 0:
        raw_files=slf._fitsdict['filename']
    else:
        raw_files=None
    core_save_master(data, filename=filename, frametype=frametype,
                     extensions=extensions, keywds=keywds, names=names,
                     raw_files=raw_files)

def core_save_master(data, filename="temp.fits", frametype="<None>",
                extensions=None, keywds=None, names=None, raw_files=None,
                     overwrite=True):
    """ Core function to write a MasterFrame

    Parameters
    ----------
    data : ndarray
    filename : str (optional)
    frametype : str (optional)
    extensions : list, optional
      Additional data images to write
    names : list, optional
      Names of the extensions
    keywds : Additional keywords for the Header
    raw_files : list or ndarray
      Names of the raw files used to generate the image

    Returns
    -------

    """
    # Check for existing
    if os.path.exists(filename) and (not overwrite):
        msgs.warn("This file already exists.  Use overwrite=True to overwrite it")
        return
    #
    msgs.info("Saving master {0:s} frame as:".format(frametype)+msgs.newline()+filename)
    hdu = fits.PrimaryHDU(data)
    hlist = [hdu]
    # Extensions
    if extensions is not None:
        for kk,exten in enumerate(extensions):
            hdu = fits.ImageHDU(exten)
            if names is not None:
                hdu.name = names[kk]
            hlist.append(hdu)
    # HDU list
    hdulist = fits.HDUList(hlist)
    # Header
    msgs.info("Writing header information")
    if raw_files is not None:
        for i in range(len(raw_files)):
            hdrname = "FRAME{0:03d}".format(i+1)
            hdulist[0].header[hdrname] = (raw_files[i], 'PYPIT: File used to generate Master {0:s}'.format(frametype))
    hdulist[0].header["FRAMETYP"] = (frametype, 'PYPIT: Master calibration frame type')
    if keywds is not None:
        for key in keywds.keys():
            hdulist[0].header[key] = keywds[key]
    # Write the file to disk
    if os.path.exists(filename):
        msgs.warn("Overwriting file:"+msgs.newline()+filename)

    hdulist.writeto(filename, overwrite=True)
    msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    return


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
