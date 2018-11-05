""" Routines related to MasterFrames"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
# import yaml
import json

from astropy.io import fits
from astropy import units

import linetools.utils

from pypeit.par import pypeitpar
from pypeit import msgs

from pypeit import debugger

def set_master_dir(redux_path, spectrograph, par):
    """
    Set the master directory auto-magically

    Args:
        redux_path: str or None
        spectrograph: Spectrograph or None
        par: ParSet or None

    Returns:
        master_dir : str
          Path of the MasterFrame directory

    """
    # Parameters
    if par is None:
        tmppar = pypeitpar.CalibrationsPar()
    else:
        if 'caldir' not in par.keys():
            tmppar = pypeitpar.CalibrationsPar()
        else:
            tmppar = par
    # Redux path
    if redux_path is None:
        redux_path = os.getcwd()
    master_dir = os.path.join(redux_path, tmppar['caldir'])
    # Spectrograph
    if spectrograph is not None:
        master_dir += '_'+spectrograph.spectrograph
    # Return
    return master_dir

def master_name(ftype, setup, mdir):
    """ Default filenames for MasterFrames

    Parameters
    ----------
    ftype : str
      Frame type
    setup : str
      Setup name
    mdir : str, optional
      Master directory

    Returns
    -------
    msname : str
    """
    name_dict = dict(bias='{:s}/MasterBias_{:s}.fits'.format(mdir, setup),
                     badpix='{:s}/MasterBadPix_{:s}.fits'.format(mdir, setup),
                     trace='{:s}/MasterTrace_{:s}'.format(mdir, setup),   # Just a root as FITS+JSON are generated
                     pinhole='{:s}/MasterPinhole_{:s}.fits'.format(mdir, setup),
                     pixelflat='{:s}/MasterPixelFlat_{:s}.fits'.format(mdir, setup),
                     illumflat='{:s}/MasterIllumFlat_{:s}.fits'.format(mdir, setup),
                     arc='{:s}/MasterArc_{:s}.fits'.format(mdir, setup),
                     wave='{:s}/MasterWave_{:s}.fits'.format(mdir, setup),
                     wv_calib='{:s}/MasterWaveCalib_{:s}.json'.format(mdir, setup),
                     tilts='{:s}/MasterTilts_{:s}.fits'.format(mdir, setup),
                     # sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.yaml'.format(mdir, setup[0], setup[-2:]),
                     sensfunc='{:s}/MasterSensFunc_{:s}_{:s}.fits'.format(mdir, setup[0], setup[-2:]),
                     )
    return name_dict[ftype]

'''
def load_master_frame(slf, mftype, det=None):
    """
    Mainly a wrapper on core_load_master_frame
       This method will be deprecated by load methods in the MasterFrame classes

    Parameters
    ----------
    slf
    mftype
    det

    Returns
    -------

    """
    # TODO -- This method will be deprecated by load methods in the MasterFrame classes
    #   Presently there are only 4 calls to this (tilts, mswave, wavecalib)
    # Were MasterFrames even desired?
    if (settings.argflag['reduce']['masters']['reuse']) or (settings.argflag['reduce']['masters']['force']):
        ret, head, _ = core_load_master_frame(mftype, slf.setup,
                                           settings.argflag['run']['directory']['master']+'_'+settings.argflag['run']['spectrograph'],
                                           force=settings.argflag['reduce']['masters']['force'])
    else:
        return None
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
'''


def load_master_frame(mftype, setup, mdir, force=False):
    """ If a MasterFrame exists, load it
    Will soon replace the load_master_frame above

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
    ms_name = master_name(mftype, setup, mdir)
    # Load
    msframe, head, file_list = _load(ms_name, exten=0, frametype=mftype, force=force)
    # Check
    if msframe is None:
        msgs.warn("No Master frame found of type {:s}: {:s}".format(mftype,ms_name))
        return None, None, None
    # Return
    return msframe, head, file_list


def _load(name, exten=0, frametype='<None>', force=False):
    """
    Low level load method for master frames
      Should mainly be called by core_load_master_frame

    Parameters
    ----------
    name : str
      Name of the master calibration file to be loaded
    exten : int, optional
    frametype : str, optional
      The type of master calibration frame being loaded.
      This keyword is only used for terminal print out.
    force : bool, optional
      Crash out if the file does not exist!

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
        msgs.error('Load from the class not this method')
        #msgs.info("Loading Master {0:s} frame:".format(frametype)+msgs.newline()+name)
        #ldict = linetools.utils.loadjson(name)
        #return ldict, None, [name]
    elif frametype == 'sensfunc':
        msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(frametype) + " from filename: {:}".format(name))

        hdu = fits.open(name)
        head = hdu[0].header
        tbl = hdu['SENSFUNC'].data
        sens_dict = {}
        sens_dict['wave'] = tbl['WAVE']
        sens_dict['sensfunc'] = tbl['SENSFUNC']
        for key in ['wave_min','wave_max','exptime','airmass','std_file','std_ra','std_dec','std_name','cal_file']:
            try:
                sens_dict[key] = head[key.upper()]
            except:
                pass
        return sens_dict, head, [name]
    elif frametype == 'trace':
        msgs.error('Load from the class not this method')
    elif frametype == 'tilts':
        msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(frametype) + " from filename: {:}".format(name))
        hdu = fits.open(name)
        head0 = hdu[0].header
        tilts = hdu[0].data
        head1 = hdu[1].header
        coeffs = hdu[1].data
        tilts_dict = {'tilts':tilts,'coeffs':coeffs,'func2D': head1['FUNC2D']} # This is the tilts_dict
        return tilts_dict, head0, [name]
    else:
        msgs.info("Loading a pre-existing master calibration frame of type: {:}".format(frametype) + " from filename: {:}".format(name))
        hdu = fits.open(name)
        #msgs.info("Master {0:s} frame loaded successfully:".format(hdu[0].header['FRAMETYP'])+msgs.newline()+name)
        head0 = hdu[0].header
        data = hdu[exten].data.astype(np.float)
        # List of files used to generate the Master frame (e.g. raw file frames)
        file_list = []
        for key in head0:
            if 'FRAME' in key:
                file_list.append(head0[key])
        return data, head0, file_list


'''
def save_masters(slf, det, mftype='all'):
    """ Save Master Frames
    THIS WILL BE DEPRECATED BIT BY BIT

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
        msgs.error("Should not get here anymore.  Save the bias in the BiasFrame class")
    # Bad Pixel
    if (mftype in ['badpix', 'all']) and ('badpix'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the trace in the TraceSlits class")
        #save_master(slf, slf._bpix[det-1],
        #                       filename=master_name('badpix', setup),
        #                       frametype='badpix')
    # Trace
    if (mftype in ['trace', 'all']) and ('trace'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the trace in the TraceSlits class")
    # Pixel Flat
    if (mftype in ['pixelflat', 'all']) and ('pixelflat'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the trace in the TraceSlits class")
        #save_master(slf, slf._mspixelflatnrm[det-1],
        #                   filename=master_name('normpixelflat', setup),
        #                   frametype='normpixelflat')
    # Pinhole Flat
    if (mftype in ['pinhole', 'all']) and ('pinhole'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._mspinhole[det-1],
                           filename=master_name('pinhole', setup),
                           frametype='pinhole')
    # Arc
    if (mftype in ['arc', 'all']) and ('arc'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the arc in the ArcImage class")
    # Wavelength image
    if (mftype in ['wave', 'all']) and ('wave'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._mswave[det-1],
                           filename=master_name('wave', setup),
                           frametype='wave')
    if (mftype in ['wv_calib', 'all']) and ('wv_calib'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the arc in the ArcImage class")
        # Wavelength fit
        #gddict = linetools.utils.jsonify(slf._wvcalib[det-1])
        #json_file = master_name('wv_calib', setup)
        #if gddict is not None:
        #    linetools.utils.savejson(json_file, gddict, easy_to_read=True, overwrite=True)
        #else:
        #    msgs.warn("The master wavelength solution has not been saved")
    # Tilts
    if (mftype in ['tilts', 'all']) and ('tilts'+setup not in settings.argflag['reduce']['masters']['loaded']):
        msgs.error("Should not get here anymore.  Save the arc in the ArcImage class")
        #save_master(slf, slf._tilts[det-1],
        #                   filename=master_name('tilts', setup),
        #                   frametype='tilts')
    # Spatial slit profile
    if (mftype in ['slitprof', 'all']) and ('slitprof'+setup not in settings.argflag['reduce']['masters']['loaded']):
        save_master(slf, slf._slitprof[det - 1],
                           filename=master_name('slitprof', setup),
                           frametype='slit profile')
'''

'''
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
'''

def save_master(data, filename="temp.fits", frametype="<None>",
                extensions=None, keywds=None, names=None, raw_files=None,
                     overwrite=True):
    """ Core function to write a MasterFrame image
    More sophisticated MasterFrame objects may be written by their own Class, e.g. TraceSlits

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
    if frametype == 'wv_calib':
        # Wavelength fit(s)
        gddict = linetools.utils.jsonify(data)
        linetools.utils.savejson(filename, gddict, easy_to_read=True, overwrite=True)
    else:  # 2D Image
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
    # Finish
    msgs.info("Master {0:s} frame saved successfully:".format(frametype)+msgs.newline()+filename)
    return


def save_sensfunc(sens_dict, outfile):
    """ Write Sensitivity function to disk as FITS table

    Args:
        sens_dict: dict
        outfile: str

    Returns:

    """
    prihdu = fits.PrimaryHDU()
    hdus = [prihdu]
    # Add critical keys from sens_dict to header
    for key in ['wave_min', 'wave_max', 'exptime', 'airmass', 'std_file', 'std_ra',
                'std_dec', 'std_name', 'cal_file']:
        try:
            prihdu.header[key.upper()] = sens_dict[key].value
        except AttributeError:
            prihdu.header[key.upper()] = sens_dict[key]
        except KeyError:
            pass # Will not require all of these

    cols = []
    cols += [fits.Column(array=sens_dict['wave'], name=str('WAVE'), format=sens_dict['wave'].dtype)]
    cols += [
        fits.Column(array=sens_dict['sensfunc'], name=str('SENSFUNC'), format=sens_dict['sensfunc'].dtype)]
    # Finish
    coldefs = fits.ColDefs(cols)
    tbhdu = fits.BinTableHDU.from_columns(coldefs)
    tbhdu.name = 'SENSFUNC'
    hdus += [tbhdu]
    # Finish
    hdulist = fits.HDUList(hdus)
    hdulist.writeto(outfile, overwrite=True)


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
