from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import sys
import copy
import glob
import getopt
import astropy.io.fits as pyfits
from astropy.time import Time
import numpy as np
from pypit import armsgs
from pypit import arparse
from pypit import arproc
from pypit import arlris
from multiprocessing import cpu_count
#from multiprocessing import Pool as mpPool
#from multiprocessing.pool import ApplyResult
#import arutils

try:
    basestring
except NameError:
    basestring = str

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger

# Logging and settings
msgs = armsgs.get_logger()
argflag = arparse.get_argflag().__dict__['argflag']
spect = arparse.get_spect().__dict__['spect']


def load_input(redname):
    """
    Load user defined input reduction file. Updates are
    made to the argflag dictionary.

    Parameters
    ----------
    redname : string
      Name of reduction script

    Returns
    -------
    parlines : list
      Input (uncommented) lines specified by the user.
      parlines is used in this routine to update the
      argflag dictionary
    datlines : list
      Input (uncommented) lines specified by the user.
      datlines contains the full data path to every
      raw exposure listed by the user
    spclines : list
      Input (uncommented) lines specified by the user.
      spclines contains a list of user-specified changes
      that should be made to the default spectrograph
      settings.
    """
    # Read in the model file
    msgs.info("Loading the input file")
    try:
        infile = open(redname, 'r')
    except IOError:
        msgs.error("The filename does not exist -"+msgs.newline()+redname)
        sys.exit()
    lines = infile.readlines()
    parlines = []
    datlines = []
    spclines = []
    rddata, rdspec = 0, 0
    for i in range(len(lines)):
        if lines[i].strip() == '': continue
        linspl = lines[i].split()
        if rddata == 1:
            if linspl[0] == 'data' and linspl[1] == 'end':
                rddata += 1
                continue
            dfname = lines[i].rstrip('\n').strip()
            # is there a comment?
            aux = dfname.split('#')
            if len(aux) > 1:  # yes, there is a comment
                dfname = aux[0].strip()
            if len(dfname) == 0:  # line is fully commented out
                continue
            elif dfname[0] == '~':
                dfname = os.path.expanduser(dfname)
            elif dfname[0] != '/':
                msgs.error("You must specify the full datapath for the file:"+msgs.newline()+dfname)
            elif len(dfname.split()) != 1:
                msgs.error("There must be no spaces when specifying the datafile:"+msgs.newline()+dfname)
            listing = glob.glob(dfname)
            for lst in listing: datlines.append(lst)
            continue
        elif rddata == 0 and linspl[0] == 'data' and linspl[1] == 'read':
            rddata += 1
            continue
        if rdspec == 1:
            if linspl[0] == 'spect' and linspl[1] == 'end':
                rdspec += 1
                continue
            spclines.append(lines[i])
            continue
        elif rdspec == 0 and linspl[0] == 'spect' and linspl[1] == 'read':
            rdspec += 1
            continue
        if lines[i].lstrip()[0] == '#': continue
        parlines.append(lines[i])
    # Do some quick checks
    if rddata == 0:
        msgs.error("You haven't specified any data!")
    elif rddata == 1:
        msgs.error("Missing 'data end' in "+redname)
    if rddata == 0:
        msgs.info("Using Default spectrograph parameters")
    elif rddata != 2:
        msgs.error("Missing 'spect end' in "+redname)
    # Check there are no duplicate inputs
    if len(datlines) != len(set(datlines)):
        msgs.error("There are duplicate files in the list of data.")
    if len(datlines) == 0: msgs.error("There are no raw data frames" + msgs.newline() +
                                      "Perhaps the path to the data is incorrect?")
    else: msgs.info("Found {0:d} raw data frames".format(len(datlines)))
    msgs.info("Input file loaded successfully")
    return parlines, datlines, spclines


def load_spect(progname, specname, spect=None, lines=None):
    """
    Load spectrograph settings

    Parameters
    ----------
    progname : string
      Name of the program
    specname : string
      Name of spectrograph settings file
    spect : dict, optional
      Properties of the spectrograph.
      If None, spect will be created, otherwise spect
      will be updated.
    lines : list, optional
      Input (uncommented) lines specified by the user.
      lines contains a list of user-specified changes
      that should be made to the default spectrograph
      settings.

    Returns
    -------
    spect : dict
      Loaded or updated properties of the spectrograph
    """
    def initialise():
        msc = dict({'ndet': 0, 'latitude': 0.0, 'longitude': 0.0, 'elevation': 0.0, 'minexp': 0., 'reduction': 'ARMLSD', 'camera': 'UNKNWN'})
        # det starts as a dict but becomes a list of dicts in set_params
        ddet = dict({'xgap': 0.0, 'ygap': 0.0, 'ysize': 1.0, 'darkcurr': 0.0, 'ronoise': [1.0], 'gain': [1.0], 'saturation': 65536.0, 'nonlinear': 1.0, 'numamplifiers': 1, 'suffix': ""})
        #
        chk = dict({})
        stf = dict({'science': [], 'standard': [], 'bias': [], 'pixflat': [], 'blzflat': [], 'arc': [], 'trace': [], 'dark': [], 'readnoise': []})
        kyw = dict({'target': '01.OBJECT', 'idname': '01.OBSTYPE', 'time': '01.MJD', 'date': '', 'equinox': '', 'ra': '', 'dec': '', 'airmass': '', 'naxis0': '01.NAXIS2', 'naxis1': '01.NAXIS1', 'exptime': '01.EXPTIME', 'hatch': '01.TRAPDOOR', 'filter1': '01.FILTNAME', 'filter2': None, 'lamps': '01.LAMPNAME', 'decker': '01.DECKNAME', 'slitwid': '01.SLITWIDTH', 'slitlen': '01.SLITLENGTH', 'detrot': '01.DETECTORROTATION', 'cdangle': '01.XDISPANGLE', 'echangle': '01.ECHELLEANGLE', 'crossdisp': '01.XDISPERS', 'dichroic': '', 'disperser': '', 'binning': ''})
        fts = dict({'numhead': 1, 'numlamps':1, 'dataext':0, 'calwin':12.0, 'timeunit':'mjd'})
        sci = dict({'index': [], 'check': dict({}), 'idname': 'OBJECT', 'canbe': None})
        std = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({})})
        pfl = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        bfl = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        trc = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        arc = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        bia = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}) })
        rn = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}) })
        drk = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}) })
        spectt = dict({'mosaic': msc, 'det': ddet, 'check': chk, 'set': stf, 'keyword': kyw, 'fits': fts, 'science': sci, 'standard': std, 'pixflat': pfl, 'blzflat': bfl, 'trace': trc, 'arc': arc, 'bias': bia, 'dark': drk, 'readnoise': rn})
        return spectt
    if lines is None:
        # Read in the default settings
        # Get the software path
        prgn_spl = progname.split('/')
        fname = ""
        for i in range(0, len(prgn_spl)-1): fname += prgn_spl[i]+"/"
        fname += 'settings/settings.'+specname
        msgs.info("Loading the "+specname+" settings")
        spect = initialise()
        with open(fname, 'r') as infile:
            lines = infile.readlines()
        spect = set_params(lines, spect, setstr="Default "+specname+" ")
    else:
        if spect is not None:
            spect = set_params(lines, spect, setstr="Infile "+specname+" ")
    return spect


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
    spect : dict
      Loaded or updated properties of the spectrograph
    """
    chks = spect['check'].keys()
    keys = spect['keyword'].keys()
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    whddict = dict({})
    for k in keys:
        fitsdict[k]=[]
    headarr = [None for k in range(spect['fits']['numhead'])]
    for i in range(len(datlines)):
        # Try to open the fits file
        try:
            for k in range(spect['fits']['numhead']):
                headarr[k] = pyfits.getheader(datlines[i], ext=spect['fits']['headext{0:02d}'.format(k+1)])
                whddict['{0:02d}'.format(spect['fits']['headext{0:02d}'.format(k+1)])] = k
        except:
            msgs.error("Error reading header from extension {0:d} of file:".format(spect['fits']['headext{0:02d}'.format(k+1)])+msgs.newline()+datlines[i])
        # Perform checks on each fits files, as specified in the settings.instrument file.
        skip = False
        for ch in chks:
            tfrhd = int(ch.split('.')[0])-1
            kchk  = '.'.join(ch.split('.')[1:])
            frhd  = whddict['{0:02d}'.format(tfrhd)]
            if spect['check'][ch] != str(headarr[frhd][kchk]).strip():
                #print ch, frhd, kchk
                #print spect['check'][ch], str(headarr[frhd][kchk]).strip()
                msgs.warn("The following file:"+msgs.newline()+datlines[i]+msgs.newline()+"is not taken with the settings.{0:s} detector".format(argflag['run']['spectrograph'])+msgs.newline()+"Remove this file, or specify a different settings file.")
                msgs.warn("Skipping the file..")
                skip = True
        if skip:
            continue
        # Now set the key values for each of the required keywords
        dspl = datlines[i].split('/')
        fitsdict['directory'].append('/'.join(dspl[:-1])+'/')
        fitsdict['filename'].append(dspl[-1])
        # Attempt to load a UTC
        utcfound = False
        for k in range(spect['fits']['numhead']):
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
        #if argflag['run']['spectrograph'] in ['lris_blue']:
        #    arlris.set_det(fitsdict, headarr[k])
        # Now get the rest of the keywords
        for kw in keys:
            if spect['keyword'][kw] is None:
                value = str('None')  # This instrument doesn't have/need this keyword
            else:
                ch = spect['keyword'][kw]
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
                if spect['fits']['timeunit']   == 's'  : value = float(value)/3600.0    # Convert seconds to hours
                elif spect['fits']['timeunit'] == 'm'  : value = float(value)/60.0      # Convert minutes to hours
                elif spect['fits']['timeunit'] in Time.FORMATS.keys() : # Astropy time format
                    if spect['fits']['timeunit'] in ['mjd']:
                        ival = float(value)
                    else:
                        ival = value
                    tval = Time(ival, scale='tt', format=spect['fits']['timeunit'])
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
    msgs.info("Headers loaded for {0:d} files successfully".format(len(datlines)))
    return fitsdict


def load_frames(fitsdict, ind, det, frametype='<None>', msbias=None,
                trim=True, transpose=False):
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
        if argflag['run']['spectrograph'] in ['lris_blue', 'lris_red']:
            temp, head0, _ = arlris.read_lris(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], det=det)
        else:
            temp = pyfits.getdata(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], spect['fits']['dataext'])
        temp = temp.astype(float)  # Let us avoid uint16
        if transpose: temp = temp.T
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
#	pool = mpPool(processes=np.min([argflag['run']['ncpus'],np.size(ind)]))
#	async_results = []
#	for i in range(np.size(ind)):
#		async_results.append(pool.apply_async(pyfits.getdata, (fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], spect['fits']['dataext'])))
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
