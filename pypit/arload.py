from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import astropy.io.fits as pyfits
from astropy.time import Time
import numpy as np
from pypit import arparse as settings
from pypit import armsgs
from pypit import arproc
from pypit import arlris
#from multiprocessing import Pool as mpPool
#from multiprocessing.pool import ApplyResult
#import arutils

try:
    basestring
except NameError:
    basestring = str

try:
    from xastropy.xutils import xdebug as debugger
except ImportError:
    import pdb as debugger

# Logging
msgs = armsgs.get_logger()


def load_headers(datlines):
    """
    Load the header information for each fits file

    Parameters
    ----------
    datlines : list
      Input (uncommented) lines specified by the user.
      datlines contains the full data path to every
      raw exposure listed by the user.

    Returns
    -------
    fitsdict : dict
      The relevant header information of all fits files
    """
    chks = settings.spect['check'].keys()
    keys = settings.spect['keyword'].keys()
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    whddict = dict({})
    for k in keys:
        fitsdict[k]=[]
    headarr = [None for k in range(settings.spect['fits']['numhead'])]
    numfiles = len(datlines)
    for i in range(numfiles):
        # Try to open the fits file
        try:
            for k in range(settings.spect['fits']['numhead']):
                headarr[k] = pyfits.getheader(datlines[i], ext=settings.spect['fits']['headext{0:02d}'.format(k+1)])
                whddict['{0:02d}'.format(settings.spect['fits']['headext{0:02d}'.format(k+1)])] = k
        except:
            msgs.error("Error reading header from extension {0:d} of file:".format(settings.spect['fits']['headext{0:02d}'.format(k+1)])+msgs.newline()+datlines[i])
        # Perform checks on each fits files, as specified in the settings.instrument file.
        skip = False
        for ch in chks:
            tfrhd = int(ch.split('.')[0])-1
            kchk  = '.'.join(ch.split('.')[1:])
            frhd  = whddict['{0:02d}'.format(tfrhd)]
            if settings.spect['check'][ch] != str(headarr[frhd][kchk]).strip():
                print(ch, frhd, kchk)
                print(settings.spect['check'][ch], str(headarr[frhd][kchk]).strip())
                msgs.warn("The following file:"+msgs.newline()+datlines[i]+msgs.newline()+"is not taken with the settings.{0:s} detector".format(settings.argflag['run']['spectrograph'])+msgs.newline()+"Remove this file, or specify a different settings file.")
                msgs.warn("Skipping the file..")
                skip = True
        if skip:
            numfiles -= 1
            continue
        # Now set the key values for each of the required keywords
        dspl = datlines[i].split('/')
        fitsdict['directory'].append('/'.join(dspl[:-1])+'/')
        fitsdict['filename'].append(dspl[-1])
        # Attempt to load a UTC
        utcfound = False
        for k in range(settings.spect['fits']['numhead']):
            if 'UTC' in headarr[k].keys():
                utc = headarr[k]['UTC']
                utcfound = True
                break
            elif 'UT' in headarr[k].keys():
                utc = headarr[k]['UT']
                utcfound = True
                break
        if utcfound:
            fitsdict['utc'].append(utc)
        else:
            fitsdict['utc'].append(None)
            msgs.warn("UTC is not listed as a header keyword in file:"+msgs.newline()+datlines[i])
        # Read binning-dependent detector properties here? (maybe read speed too)
        #if settings.argflag['run']['spectrograph'] in ['lris_blue']:
        #    arlris.set_det(fitsdict, headarr[k])
        # Now get the rest of the keywords
        for kw in keys:
            if settings.spect['keyword'][kw] is None:
                value = str('None')  # This instrument doesn't have/need this keyword
            else:
                ch = settings.spect['keyword'][kw]
                try:
                    tfrhd = int(ch.split('.')[0])-1
                except ValueError:
                    value = ch  # Keyword given a value. Only a string allowed for now
                else:
                    frhd = whddict['{0:02d}'.format(tfrhd)]
                    kchk = '.'.join(ch.split('.')[1:])
                    try:
                        value = headarr[frhd][kchk]
                    except KeyError: # Keyword not found in header
                        msgs.warn("{:s} keyword not in header. Setting to None".format(kchk))
                        value=str('None')
            # Convert the input time into hours
            if kw == 'time':
                if settings.spect['fits']['timeunit']   == 's'  : value = float(value)/3600.0    # Convert seconds to hours
                elif settings.spect['fits']['timeunit'] == 'm'  : value = float(value)/60.0      # Convert minutes to hours
                elif settings.spect['fits']['timeunit'] in Time.FORMATS.keys() : # Astropy time format
                    if settings.spect['fits']['timeunit'] in ['mjd']:
                        ival = float(value)
                    else:
                        ival = value
                    tval = Time(ival, scale='tt', format=settings.spect['fits']['timeunit'])
                    # dspT = value.split('T')
                    # dy,dm,dd = np.array(dspT[0].split('-')).astype(np.int)
                    # th,tm,ts = np.array(dspT[1].split(':')).astype(np.float64)
                    # r=(14-dm)/12
                    # s,t=dy+4800-r,dm+12*r-3
                    # jdn = dd + (153*t+2)/5 + 365*s + s/4 - 32083
                    # value = jdn + (12.-th)/24 + tm/1440 + ts/86400 - 2400000.5  # THIS IS THE MJD
                    value = tval.mjd * 24.0 # Put MJD in hours
                else:
                    msgs.error('Bad time unit')
            # Put the value in the keyword
            typv = type(value)
            if typv is int or typv is np.int_:
                fitsdict[kw].append(value)
            elif typv is float or typv is np.float_:
                fitsdict[kw].append(value)
            elif isinstance(value, basestring) or typv is np.string_:
                fitsdict[kw].append(value.strip())
            else:
                debugger.set_trace()
                msgs.bug("I didn't expect useful headers to contain type {0:s}".format(typv).replace('<type ','').replace('>',''))

        msgs.info("Successfully loaded headers for file:"+msgs.newline()+datlines[i])
    del headarr
    # Convert the fitsdict arrays into numpy arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])
    msgs.info("Headers loaded for {0:d} files successfully".format(numfiles))
    if numfiles != len(datlines):
        msgs.warn("Headers were not loaded for {0:d} files".format(len(datlines) - numfiles))
    if numfiles == 0:
        msgs.error("The headers could not be read from the input data files." + msgs.newline() +
                   "Please check that the settings file matches the data.")
    return fitsdict


def load_frames(fitsdict, ind, det, frametype='<None>', msbias=None, trim=True):
    """
    Load data frames, usually raw.
    Bias subtract (if not msbias!=None) and trim (if True)

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    ind : list or array
      integers of indices
    det : int
      Detector number, starts at 1

    Returns
    -------
    frames : ndarray (3 dimensional)
      One image per ind
    """
    def load_indfr(name,ext):
        msgs.work("Trim and overscan has not been applied")
        temp = pyfits.getdata(name, ext)
        return temp

    msgs.info("Loading individual {0:s} frames".format(frametype))
    if np.size(ind) == 0:
        msgs.warn("No {0:s} frames to load".format(frametype))
        return None
    msgs.work("Implement multiprocessing here (better -- at the moment it's slower than not) to speed up data reading")
    for i in range(np.size(ind)):
        # Instrument specific read
        if settings.argflag['run']['spectrograph'] in ['lris_blue', 'lris_red']:
            temp, head0, _ = arlris.read_lris(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], det=det)
        else:
            temp = pyfits.getdata(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], settings.spect['fits']['dataext'])
        temp = temp.astype(np.float)  # Let us avoid uint16
        if settings.argflag['trace']['dispersion']['direction'] == 1:
            temp = temp.T
        if msbias is not None:
            if type(msbias) is np.ndarray:
                temp -= msbias  # Subtract the master bias frame
            elif type(msbias) is str:
                if msbias == "overscan":
                    arproc.sub_overscan(temp, det)
                else:
                    msgs.error("Could not subtract bias level when loading {0:s} frames".format(frametype))
            if trim:
                temp = arproc.trim(temp, det)
        if i == 0:
            frames = np.zeros((temp.shape[0], temp.shape[1], np.size(ind)))
            frames[:,:,i] = temp.copy()
        else:
            frames[:,:,i] = temp.copy()
        del temp
#	pool = mpPool(processes=np.min([settings.argflag['run']['ncpus'],np.size(ind)]))
#	async_results = []
#	for i in range(np.size(ind)):
#		async_results.append(pool.apply_async(pyfits.getdata, (fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], settings.spect['fits']['dataext'])))
#	pool.close()
#	pool.join()
#	map(ApplyResult.wait, async_results)
#	for j in range(np.size(ind)):
#		if j == 0:
#			temp = async_results[j].get()
#			frames = np.zeros((temp.shape[0], temp.shape[1], np.size(ind)))
#			if msbias is None:
#				frames[:,:,i] = temp
#			else:
#				frames[:,:,i] = temp - msbias
#			del temp
#		else:
#			if msbias is None:
#				frames[:,:,i] = async_results[j].get()
#			else:
#				frames[:,:,i] = async_results[j].get() - msbias
    if np.size(ind) == 1:
        msgs.info("Loaded {0:d} {1:s} frame successfully".format(np.size(ind), frametype))
    else:
        msgs.info("Loaded {0:d} {1:s} frames successfully".format(np.size(ind), frametype))
    return frames


def load_extraction(name, frametype='<None>', wave=True):
    msgs.info("Loading a pre-existing {0:s} extraction frame:".format(frametype)+msgs.newline()+name)
    props_savas = dict({"ORDWN":"ordwnum"})
    props = dict({})
    props_allow = props_savas.keys()
    infile = pyfits.open(name)
    sciext = np.array(infile[0].data, dtype=np.float)
    sciwav = -999999.9*np.ones((sciext.shape[0],sciext.shape[1]))
    hdr = infile[0].header
    norders = hdr["NUMORDS"]
    pxsz    = hdr["PIXSIZE"]
    props = dict({})
    for o in range(norders):
        hdrname = "CDELT{0:03d}".format(o+1)
        cdelt = hdr[hdrname]
        hdrname = "CRVAL{0:03d}".format(o+1)
        crval = hdr[hdrname]
        hdrname = "CLINV{0:03d}".format(o+1)
        clinv = hdr[hdrname]
        hdrname = "CRPIX{0:03d}".format(o+1)
        crpix = hdr[hdrname]
        hdrname = "CNPIX{0:03d}".format(o+1)
        cnpix = hdr[hdrname]
        sciwav[:cnpix,o] = 10.0**(crval + cdelt*(np.arange(cnpix)-crpix))
        #sciwav[:cnpix,o] = clinv * 10.0**(cdelt*(np.arange(cnpix)-crpix))
        #sciwav[:cnpix,o] = clinv * (1.0 + pxsz/299792.458)**np.arange(cnpix)
    for k in props_allow:
        prsav = np.zeros(norders)
        try:
            for o in range(norders):
                hdrname = "{0:s}{1:03d}".format(k,o+1)
                prsav[o] = hdr[hdrname]
            props[props_savas[k]] = prsav.copy()
        except:
            pass
    del infile, hdr, prsav
    if wave is True:
        return sciext, sciwav, props
    else:
        return sciext, props


def load_master(name, exten=0, frametype='<None>'):
    """
    Load a pre-existing master calibration frame

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
    frame : ndarray
      The data from the master calibration frame
    """
    if frametype is None:
        msgs.info("Loading a pre-existing master calibration frame")
        try:
            hdu = pyfits.open(name)
        except:
            msgs.error("Master calibration file does not exist:"+msgs.newline()+name)
        msgs.info("Master {0:s} frame loaded successfully:".format(hdu[0].header['FRAMETYP'])+msgs.newline()+name)
        head = hdu[0].header
        data = hdu[exten].data.astype(np.float)
        return data, head
        #return np.array(infile[0].data, dtype=np.float)
    else:
        from linetools import utils as ltu
        msgs.info("Loading Master {0:s} frame:".format(frametype)+msgs.newline()+name)
        if frametype == 'wv_calib':
            ldict = ltu.loadjson(name)
            return ldict
        else:
            # Load
            hdu = pyfits.open(name)
            head = hdu[0].header
            data = hdu[exten].data.astype(np.float)
            return data, head
        #return np.array(pyfits.getdata(name, 0), dtype=np.float)


def load_ordloc(fname):
    # Load the files
    mstrace_bname, mstrace_bext = os.path.splitext(fname)
    lname = mstrace_bname+"_ltrace"+mstrace_bext
    rname = mstrace_bname+"_rtrace"+mstrace_bext
    # Load the order locations
    ltrace = np.array(pyfits.getdata(lname, 0),dtype=np.float)
    msgs.info("Loaded left order locations for frame:"+msgs.newline()+fname)
    rtrace = np.array(pyfits.getdata(rname, 0),dtype=np.float)
    msgs.info("Loaded right order locations for frame:"+msgs.newline()+fname)
    return ltrace, rtrace


def load_tilts(fname):
    # Load the files
    msarc_bname, msarc_bext = os.path.splitext(fname)
    tname = msarc_bname+"_tilts"+msarc_bext
    sname = msarc_bname+"_satmask"+msarc_bext
    # Load the order locations
    tilts = np.array(pyfits.getdata(tname, 0),dtype=np.float)
    msgs.info("Loaded order tilts for frame:"+msgs.newline()+fname)
    satmask = np.array(pyfits.getdata(sname, 0),dtype=np.float)
    msgs.info("Loaded saturation mask for frame:"+msgs.newline()+fname)
    return tilts, satmask


def load_1dspec(fname, exten=1, extract='opt'):
    """
    Parameters
    ----------
    fname : str
      Name of the file
    exten : int
      Extension of the spectrum
    extract : str, optional
      Extraction type ('opt', 'box')

    Returns
    -------
    spec : XSpectrum1D

    """
    from linetools.spectra.xspectrum1d import XSpectrum1D
    # Keywords for Table
    rsp_kwargs = {}
    rsp_kwargs['wave_tag'] = '{:s}_wave'.format(extract)
    rsp_kwargs['flux_tag'] = '{:s}_counts'.format(extract)
    rsp_kwargs['var_tag'] = '{:s}_var'.format(extract)
    # Load
    spec = XSpectrum1D.from_file(fname, exten=exten, **rsp_kwargs)
    # Return
    return spec

def waveids(fname):
    infile = pyfits.open(fname)
    pixels=[]
    msgs.info("Loading fitted arc lines")
    try:
        o = 1
        while True:
            pixels.append(infile[o].data.astype(np.float))
            o+=1
    except:
        pass
    return pixels
