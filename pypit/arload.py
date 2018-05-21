from __future__ import (print_function, absolute_import, division, unicode_literals)
from future.utils import iteritems

try:
    basestring
except NameError:
    basestring = str

import os

import numpy as np

from astropy import units
from astropy.time import Time
from astropy.io import fits
from astropy.table import Table
import yaml

import linetools.utils
from linetools.spectra.xspectrum1d import XSpectrum1D

from pypit import msgs
from pypit import arparse as settings
from pypit import arproc
from pypit import arlris
from pypit import arspecobj
from pypit import ardeimos
from pypit import ardebug as debugger


def load_headers(datlines, settings_spect, settings_argflag):
    """ Load the header information for each fits file
    The cards of interest are specified in the instrument settings file
    A check of specific cards is performed if specified in settings

    Parameters
    ----------
    datlines : list
      Input (uncommented) lines specified by the user.
      datlines contains the full data path to every
      raw exposure provided by the user.

    Returns
    -------
    fitstbl : Table
      The relevant header information of all fits files
    keylst : list
    """
    def generate_updates(dct, keylst, keys, whddict, headarr):
        """ Generate a list of settings to be updated
        """
        for (key, value) in iteritems(dct):
            keys += [str(key)]
            if isinstance(value, dict):
                generate_updates(value, keylst, keys, whddict, headarr)
            else:
                try:
                    tfrhd = int(value.split('.')[0]) - 1
                    kchk = '.'.join(value.split('.')[1:])
                    frhd = whddict['{0:02d}'.format(tfrhd)]
                    hdrval = headarr[frhd][kchk]
                    if keys[0] not in ["check", "keyword"]:
                        keylst += [str(' ').join(keys) + str(" ") +
                                   str("{0}\n".format(hdrval).replace(" ", ""))]
                        keylst[-1] = keylst[-1].split()
                except (AttributeError, ValueError, KeyError):
                    pass
            del keys[-1]

    chks = list(settings_spect['check'].keys())
    # FITS dict/table keys
    keys = list(settings_spect['keyword'].keys())
    # Init
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    headdict = {}
    for k in range(settings_spect['fits']['numhead']):
        headdict[k] = []
    whddict = dict({})
    for k in keys:
        fitsdict[k]=[]
    numfiles = len(datlines)
    # Loop on files
    for i in range(numfiles):
        # Try to open the fits file
        headarr = ['None' for k in range(settings_spect['fits']['numhead'])]
        try:
            for k in range(settings_spect['fits']['numhead']):
                headarr[k] = fits.getheader(datlines[i], ext=settings_spect['fits']['headext{0:02d}'.format(k+1)])
                whddict['{0:02d}'.format(settings_spect['fits']['headext{0:02d}'.format(k+1)])] = k
        except:
            if settings_argflag['run']['setup']:
                msgs.warn("Bad header in extension {0:d} of file:".format(settings_spect['fits']['headext{0:02d}'.format(k+1)])+msgs.newline()+datlines[i])
                msgs.warn("Proceeding on the hopes this was a calibration file, otherwise consider removing.")
            else:
                msgs.error("Error reading header from extension {0:d} of file:".format(settings_spect['fits']['headext{0:02d}'.format(k+1)])+msgs.newline()+datlines[i])
        # Save the headers into its dict
        for k in range(settings_spect['fits']['numhead']):
            headdict[k].append(headarr[k].copy())
            #tmp = [head.copy() for head in headarr]
            #allhead.append(tmp)
        # Perform checks on each FITS file, as specified in the settings instrument file.
        skip = False
        for ch in chks:
            tfrhd = int(ch.split('.')[0])-1
            kchk = '.'.join(ch.split('.')[1:])
            frhd = whddict['{0:02d}'.format(tfrhd)]
            # JFH changed to in instead of !=
            if ((settings_spect['check'][ch] in str(headarr[frhd][kchk]).strip()) == False):
                print(ch, frhd, kchk)
                print(settings_spect['check'][ch], str(headarr[frhd][kchk]).strip())
                msgs.warn("The following file:"+msgs.newline()+datlines[i]+msgs.newline()+"is not taken with the settings.{0:s} detector".format(settings_argflag['run']['spectrograph'])+msgs.newline()+"Remove this file, or specify a different settings file.")
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
        for k in range(settings_spect['fits']['numhead']):
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
            fitsdict['utc'].append('None') # Changed from None so it writes to disk
            msgs.warn("UTC is not listed as a header keyword in file:"+msgs.newline()+datlines[i])
        # Read binning-dependent detector properties here? (maybe read speed too)
        #if settings_argflag['run']['spectrograph'] in ['keck_lris_blue']:
        #    arlris.set_det(fitsdict, headarr[k])
        # Now get the rest of the keywords
        for kw in keys:
            if settings_spect['keyword'][kw] is None:
                value = str('None')  # This instrument doesn't have/need this keyword
            else:
                ch = settings_spect['keyword'][kw]
                # Parse the header extension holding the key
                try:
                    tfrhd = int(ch.split('.')[0])-1
                except ValueError:
                    # Keyword given a value. Only a string allowed for now
                    value = ch
                else:
                    # Load up the header
                    frhd = whddict['{0:02d}'.format(tfrhd)]
                    kchk = '.'.join(ch.split('.')[1:])
                    try:
                        value = headarr[frhd][kchk]
                    except KeyError: # Keyword not found in header
                        msgs.warn("{:s} keyword not in header. Setting to None".format(kchk))
                        value=str('None')
            # Convert the input time into hours -- Should we really do this here??
            if kw == 'time':
                if settings_spect['fits']['timeunit']   == 's'  : value = float(value)/3600.0    # Convert seconds to hours
                elif settings_spect['fits']['timeunit'] == 'm'  : value = float(value)/60.0      # Convert minutes to hours
                elif settings_spect['fits']['timeunit'] in Time.FORMATS.keys() : # Astropy time format
                    if settings_spect['fits']['timeunit'] in ['mjd']:
                        ival = float(value)
                    else:
                        ival = value
                    tval = Time(ival, scale='tt', format=settings_spect['fits']['timeunit'])
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
            elif typv is bool or typv is np.bool_:
                fitsdict[kw].append(value)
            else:
                msgs.bug("I didn't expect a useful header ({0:s}) to contain type {1:s}".format(kw, typv).replace('<type ','').replace('>',''))

        msgs.info("Successfully loaded headers for file:" + msgs.newline() + datlines[i])

    # Check if any other settings require header values to be loaded
    msgs.info("Checking spectrograph settings for required header information")
    # Just use the header info from the last file
    keylst = []
    generate_updates(settings_spect.copy(), keylst, [], whddict, headarr)

    # Convert the fitsdict arrays into numpy arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])
    #
    msgs.info("Headers loaded for {0:d} files successfully".format(numfiles))
    if numfiles != len(datlines):
        msgs.warn("Headers were not loaded for {0:d} files".format(len(datlines) - numfiles))
    if numfiles == 0:
        msgs.error("The headers could not be read from the input data files." + msgs.newline() +
                   "Please check that the settings file matches the data.")
    #  Might have to carry the headers around separately
    #    as packing them into a table could be problematic..
    #for key in headdict.keys():
    #    fitsdict['head{:d}'.format(key)] = headdict[key]
    # Return after creating a Table
    fitstbl = Table(fitsdict)
    fitstbl.sort('time')

    # TODO -- Remove the following (RC has an idea)
    # Instrument specific
    if settings_argflag['run']['spectrograph'] == 'keck_deimos':
        # Handle grating position
        for gval in [3,4]:
            gmt = fitstbl['gratepos'] == gval
            fitstbl['dispangle'][gmt] = fitstbl['g3tltwav'][gmt]
    return fitstbl, keylst


def load_frames(fitsdict, ind, det, frametype='<None>', msbias=None, trim=True):
    """  Now a wrapper on several core methods.  This might well get broken
    down further in a future Refactor.

    Load data frames, usually raw.
    Bias subtract (if not msbias!=None) and trim (if True)

    Parameters
    ----------

    fitsdict : dict
        Contains relevant information from fits header files
    ind : list or array
        integers of indices
    det : int
    msbias : ndarray, str (optional)
    trim : bool (optional)

    Returns
    -------

    frames : ndarray
      3D with the 3rd index corresponding to the frames returned
    """
    # Wrap me
    dnum = settings.get_dnum(det)
    spectrograph = settings.argflag['run']['spectrograph']
    if 'dataext01' in settings.spect[dnum].keys():
        dataexto01 = settings.spect[dnum]['dataext01']
    else:
        dataexto01 = None
    disp_dir = settings.argflag['trace']['dispersion']['direction']
    numamplifiers = settings.spect[dnum]['numamplifiers']
    # Build datasecs, oscansec
    datasecs, oscansecs = [], []
    for jj in range(numamplifiers):
        datasec = "datasec{0:02d}".format(jj+1)
        datasecs.append(settings.spect[dnum][datasec])
        oscansec = "oscansec{0:02d}".format(jj+1)
        oscansecs.append(settings.spect[dnum][oscansec])

    # Now run
    for i in range(len(ind)):
        raw_file = fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]]
        temp, head0 = load_raw_frame(spectrograph, raw_file, det,
                              dataext01=dataexto01, disp_dir=disp_dir)

        # TODO -- Take these next two steps out and put in a arproc.proc_image() method
        # Bias subtract?
        if msbias is not None:
            temp = arproc.bias_subtract(temp, msbias, numamplifiers=numamplifiers,
                                        datasec=datasecs, oscansec=oscansecs)

        if trim:
            # Trim
            temp = arproc.trim(temp, numamplifiers, datasecs)

        # Save
        if i == 0:
            frames = np.zeros((temp.shape[0], temp.shape[1], len(ind)))
            frames[:,:,i] = temp.copy()
        else:
            frames[:,:,i] = temp.copy()
        del temp

    # Finish
    if len(ind) == 1:
        msgs.info("Loaded {0:d} {1:s} frame successfully".format(len(ind), frametype))
    else:
        msgs.info("Loaded {0:d} {1:s} frames successfully".format(len(ind), frametype))
    return frames


def load_raw_frame(spectrograph, raw_file, det, dataext01=None, disp_dir=0):
    """
    Load data frames, usually raw.
    Bias subtract (if not msbias!=None) and trim (if True)

    Parameters
    ----------
    raw_file : str
       Full path to raw_file
    det : int
      Detector number requested, starts at 1
    disp_dir : int, optional
      if 1, Transpose the image to align spectral dimension with columns

    Returns
    -------
    frame : ndarray
      the raw_frame
    """
    msgs.info("Loading raw_file: {:s}".format(raw_file))
    # Get detector number
    msgs.work("Implement multiprocessing here (better -- at the moment it's slower than not) to speed up data reading")
    # Instrument specific read
    if spectrograph in ['keck_lris_blue', 'keck_lris_red']:
        #temp, head0, _ = arlris.read_lris(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], det=det)
        temp, head0, _ = arlris.read_lris(raw_file, det=det)
    elif spectrograph in ['keck_deimos']:
        temp, head0, _ = ardeimos.read_deimos(raw_file)
        #temp, head0, _ = ardeimos.read_deimos(fitsdict['directory'][ind[i]] + fitsdict['filename'][ind[i]])
    else:
        #hdulist = fits.open(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]])
        hdulist = fits.open(raw_file)
        #temp = hdulist[settings.spect[dnum]['dataext01']].data
        temp = hdulist[dataext01].data
        head0 = hdulist[0].header
    temp = temp.astype(np.float)  # Let us avoid uint16
    #if settings.argflag['trace']['dispersion']['direction'] == 1:
    if disp_dir == 1:
        temp = temp.T
    return temp, head0


def load_extraction(name, frametype='<None>', wave=True):
    msgs.info("Loading a pre-existing {0:s} extraction frame:".format(frametype)+msgs.newline()+name)
    props_savas = dict({"ORDWN":"ordwnum"})
    props = dict({})
    props_allow = props_savas.keys()
    infile = fits.open(name)
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

    Ret        # HAS NOT BEEN DEVELOPED SINCE THE SetupClass refactor;  no test case..urns
    -------
    frame : ndarray or dict
      The data from the master calibration frame
    head : str (or None)
    """
    if frametype == 'wv_calib':
        msgs.info("Loading Master {0:s} frame:".format(frametype)+msgs.newline()+name)
        ldict = linetools.utils.loadjson(name)
        return ldict, None
    elif frametype == 'sensfunc':
        with open(name, 'r') as f:
            sensfunc = yaml.load(f)
        sensfunc['wave_max'] = sensfunc['wave_max']*units.AA
        sensfunc['wave_min'] = sensfunc['wave_min']*units.AA
        return sensfunc, None
    else:
        msgs.info("Loading a pre-existing master calibration frame")
        try:
            hdu = fits.open(name)
        except IOError:
            if settings.argflag['reduce']['masters']['force']:
                msgs.error("Master calibration file does not exist:"+msgs.newline()+name)
            else:
                msgs.warn("Could not read Master calibration file:"+msgs.newline()+name)
                raise IOError
        msgs.info("Master {0:s} frame loaded successfully:".format(hdu[0].header['FRAMETYP'])+msgs.newline()+name)
        head = hdu[0].header
        data = hdu[exten].data.astype(np.float)
        return data, head


def load_ordloc(fname):
    # Load the files
    mstrace_bname, mstrace_bext = os.path.splitext(fname)
    lname = mstrace_bname+"_ltrace"+mstrace_bext
    rname = mstrace_bname+"_rtrace"+mstrace_bext
    # Load the order locations
    ltrace = np.array(fits.getdata(lname, 0),dtype=np.float)
    msgs.info("Loaded left order locations for frame:"+msgs.newline()+fname)
    rtrace = np.array(fits.getdata(rname, 0),dtype=np.float)
    msgs.info("Loaded right order locations for frame:"+msgs.newline()+fname)
    return ltrace, rtrace


def load_specobj(fname):
    """ Load a spec1d file into a list of SpecObjExp objects
    Parameters
    ----------
    fname : str

    Returns
    -------
    specobjs : list of SpecObjExp
    """
    speckeys = ['wave', 'sky', 'mask', 'flam', 'flam_var', 'var', 'counts']
    #
    specobjs = []
    hdulist = fits.open(fname)
    for hdu in hdulist:
        if hdu.name == 'PRIMARY':
            continue
        # Parse name
        objp = hdu.name.split('-')
        # Load data
        spec = Table(hdu.data)
        shape = (len(spec), 1024)  # 2nd number is dummy
        # Init
        specobj = arspecobj.SpecObjExp(shape, 'dum_config', int(objp[-1][1:]),
            int(objp[-2][1:]), [float(objp[1][1:])/10000.]*2, 0.5,
            float(objp[0][1:])/1000., 'unknown')
        # Add spectrum
        if 'box_counts' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.boxcar[skey] = spec['box_{:s}'.format(skey)].data
                except KeyError:
                    pass
            # Add units on wave
            specobj.boxcar['wave'] = specobj.boxcar['wave'] * units.AA

        if 'opt_counts' in spec.keys():
            for skey in speckeys:
                try:
                    specobj.optimal[skey] = spec['opt_{:s}'.format(skey)].data
                except KeyError:
                    pass
        # Append
        specobjs.append(specobj)
    # Return
    return specobjs


def load_tilts(fname):
    # Load the files
    msarc_bname, msarc_bext = os.path.splitext(fname)
    tname = msarc_bname+"_tilts"+msarc_bext
    sname = msarc_bname+"_satmask"+msarc_bext
    # Load the order locations
    tilts = np.array(fits.getdata(tname, 0),dtype=np.float)
    msgs.info("Loaded order tilts for frame:"+msgs.newline()+fname)
    satmask = np.array(fits.getdata(sname, 0),dtype=np.float)
    msgs.info("Loaded saturation mask for frame:"+msgs.newline()+fname)
    return tilts, satmask


def load_1dspec(fname, exten=None, extract='opt', objname=None, flux=False):
    """
    Parameters
    ----------
    fname : str
      Name of the file
    exten : int, optional
      Extension of the spectrum
      If not given, all spectra in the file are loaded
    extract : str, optional
      Extraction type ('opt', 'box')
    objname : str, optional
      Identify extension based on input object name
    flux : bool, optional
      Return fluxed spectra?

    Returns
    -------
    spec : XSpectrum1D

    """
    # Keywords for Table
    rsp_kwargs = {}
    rsp_kwargs['wave_tag'] = '{:s}_wave'.format(extract)
    if flux:
        rsp_kwargs['flux_tag'] = '{:s}_flam'.format(extract)
        rsp_kwargs['var_tag'] = '{:s}_flam_var'.format(extract)
    else:
        rsp_kwargs['flux_tag'] = '{:s}_counts'.format(extract)
        rsp_kwargs['var_tag'] = '{:s}_var'.format(extract)
    # Identify extension from objname?
    if objname is not None:
        hdulist = fits.open(fname)
        hdu_names = [hdu.name for hdu in hdulist]
        exten = hdu_names.index(objname)
        if exten < 0:
            msgs.error("Bad input object name: {:s}".format(objname))
    # Load
    spec = XSpectrum1D.from_file(fname, exten=exten, **rsp_kwargs)
    # Return
    return spec

def waveids(fname):
    infile = fits.open(fname)
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
