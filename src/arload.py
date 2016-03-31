import os
import sys
import copy
import glob
import getopt
import astropy.io.fits as pyfits
from astropy.time import Time
import numpy as np
import armsgs
import arproc
import arlris
from multiprocessing import cpu_count
#from multiprocessing import Pool as mpPool
#from multiprocessing.pool import ApplyResult
#import arutils

try:
    from xastropy.xutils.xdebug import set_trace
#    from xastropy.xutils import xdebug as xdb
except ImportError:
    from pdb import set_trace

# Logging
msgs = armsgs.get_logger()


def argflag_init():
    """
    Initialise the default settings called argflag
    """
    rna = dict({'prognm':'pypit.py', 'redname':'filelist.red', 'spectrograph':'hamspec', 'masterdir':'MasterFrames', 'plotsdir':'Plots', 'scidir':'Science', 'ncpus':-1, 'nsubpix':5, 'calcheck':False, 'qcontrol':True, 'preponly':False, 'stopcheck':False, 'use_idname':False})
    red = dict({'locations':None, 'nlcorr':False, 'trim':True, 'badpix':True, 'usebias':'bias', 'usetrace':'trace', 'usearc':'arc', 'usewave':'wave', 'useflat':'pixflat', 'subdark':False, 'flatfield':True, 'FlatMethod':'SpatialFit', 'FlatParams':[0], 'bgsubtraction':True, 'arcmatch':2.0, 'flatmatch':2.0, 'calibrate':True, 'fluxcalibrate':True, 'extraction':'2D', 'oscanMethod':'polynomial', 'oscanParams':[1], 'heliocorr':True, 'pixelsize':2.5})
    csq = dict({'atol':1.0E-3, 'xtol':1.0E-10, 'gtol':1.0E-10, 'ftol':1.0E-10, 'fstep':2.0})
    opa = dict({'verbose':2, 'sorted':None, 'plots':True, 'overwrite':False})
    sci = dict({'load':dict({'extracted':False}),
                'extraction':dict({'method':'2D', 'profile':'gaussian', 'centorder':1, 'widthorder':1, 'function':'legendre', 'pcacent':[1,0], 'pcawidth':[1,0], 'bintrace':10})
                })
    pfl = dict({'comb':dict({'method':None, 'rej_cosmicray':50.0, 'rej_lowhigh':[0,0], 'rej_level':[3.0,3.0], 'sat_pix':'reject', 'set_allrej':'median'}) })
    bfl = dict({'comb':dict({'method':None, 'rej_cosmicray':50.0, 'rej_lowhigh':[0,0], 'rej_level':[3.0,3.0], 'sat_pix':'reject', 'set_allrej':'median'}) })
    trc = dict({'comb':dict({'method':'weightmean', 'rej_cosmicray':50.0, 'rej_lowhigh':[0,0], 'rej_level':[3.0,3.0], 'sat_pix':'reject', 'set_allrej':'maxnonsat'}),
                'disp':dict({'window':None, 'direction':None}),
                'orders':dict({'tilts':'trace', 'pcatilt':[2,1,0], 'tiltorder':1, 'tiltdisporder':2, 'function':'polynomial', 'polyorder':2, 'diffpolyorder':2, 'fracignore':0.6, 'sigdetect':3.0, 'pca':[3,2,1,0,0,0], 'pcxpos':3, 'pcxneg':3}) })
    arc = dict({'comb':dict({'method':'weightmean', 'rej_cosmicray':50.0, 'rej_lowhigh':[0,0], 'rej_level':[3.0,3.0], 'sat_pix':'reject', 'set_allrej':'maxnonsat'}),
                'extract':dict({'binby':1.0}),
                'load':dict({'extracted':False, 'calibrated':False}),
                'calibrate':dict({'cwpolyorder':2, 'threshold':3.0, 'polyorderpri':4, 'polyordersec':8, 'pcapri':[4,2,2,0], 'pcasec':[5,4,3,2,1,1], 'detection':6.0, 'method':'simple', 'nfitpix':7, 'idfile':'wave_ThAr_3100-11000.npy', 'linelist':'arclist.ThAr', 'numsearch':20, 'sigmacut':2.0}),
                })
    bia = dict({'comb':dict({'method':'mean', 'rej_cosmicray':20.0, 'rej_lowhigh':[0,0], 'rej_level':[3.0,3.0], 'sat_pix':'reject', 'set_allrej':'median'}) })
    drk = dict({})
    argflag = dict({'run':rna, 'reduce':red, 'science':sci, 'pixflat':pfl, 'blzflat':bfl, 'trace':trc, 'arc':arc, 'bias':bia, 'dark':drk, 'chisq':csq, 'out':opa})
    return argflag


def cpucheck(ncpu, curcpu=0):
    cpucnt = cpu_count()
    if ncpu == 'all':
        ncpu = cpucnt  # Use all available cpus
        if cpucnt != curcpu: msgs.info("Setting {0:d} CPUs".format(ncpu))
    elif ncpu is None:
        ncpu = cpucnt-1  # Use all but 1 available cpus
        if ncpu != curcpu:
            msgs.info("Setting {0:d} CPUs".format(ncpu))
    else:
        try:
            ncpu = int(ncpu)
            if ncpu > cpucnt:
                msgs.warn("You don't have {0:d} CPUs!".format(ncpu))
                ncpu = cpucnt
            elif ncpu < 0:
                ncpu += cpucnt
            if ncpu != curcpu:
                msgs.info("Setting {0:d} CPUs".format(ncpu))
        except:
            msgs.error("Incorrect argument given for number of CPUs"+msgs.newline()+"Please choose from -"+msgs.newline()+"all, 1..."+str(cpucnt))
            if cpucnt == 1:
                if cpucnt != curcpu:
                    msgs.info("Setting 1 CPU")
                ncpu = 1
            else:
                ncpu = cpu_count()-1
                if ncpu != curcpu:
                    msgs.info("Setting {0:d} CPUs".format(ncpu))
    return ncpu


def optarg(argflag, argv, pypname):
    """
    Load the command line options and arguments

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    argv : list
      command line arguments
    pypname : string
      Name of the pipeline settings to be loaded

    Returns
    -------
    argflag : dict
      Arguments and flags used for reduction
    """

    # Load the default settings
    prgn_spl = argv[0].split('/')
    tfname = ""
    for i in range(0,len(prgn_spl)-2): tfname += prgn_spl[i]+"/"
    fname = tfname + prgn_spl[-2] + '/settings.' + pypname
    argflag = load_settings(fname, argflag)
    argflag['run']['prognm'] = argv[0]
    argflag['run']['pypitdir'] = tfname
    # Load options from command line
    opt, arg = None, None
    try:
        opt, arg = getopt.getopt(argv[1:], 'hc:v:', ['help',
                                                     'cpus',
                                                     'verbose',
                                                     ])
    except getopt.GetoptError, err:
        msgs.error(err.msg)
        msgs.usage(None)
    for o, a in opt:
        if o in ('-h', '--help'): msgs.usage(None)
        elif o in ('-c', '--cpus'): argflag['run']['ncpus'] = a
        elif o in ('-v', '--verbose'): argflag['out']['verbose'] = int(a)

    #######################
    # Now do some checks: #
    #######################

    # Check requested CPUs
    argflag['run']['ncpus'] = cpucheck(argflag['run']['ncpus'])

    # Assign filelist:
    argflag['run']['redname'] = arg[0]
    return argflag


def set_params_wtype(tvalue, svalue, lines="", setstr="", argnum=3):
    try:
        if type(tvalue) is int:
            tvalue = int(svalue)
        elif type(tvalue) is str:
            if svalue.lower() == 'none': tvalue = None
            elif svalue[0] == '[' and svalue[-1] == ']' and ',' in svalue and len(svalue.split(':')) == 3: tvalue = load_sections(svalue, strtxt=setstr)
            else: tvalue = svalue
        elif type(tvalue) is float:
            tvalue = float(svalue)
        elif type(tvalue) is list:
            if ',' in svalue:
                temp = svalue.lstrip('([').rstrip(')]').split(',')
                addarr = []
                # Find the type of the array elements
                for i in temp:
                    if i.lower() == 'none': # None type
                        addarr += [None]
                    elif i.lower() == 'true' or i.lower() == 'false': # bool type
                        addarr += [i.lower() in ['true']]
                    elif ',' in i: # a list
                        addarr += i.lstrip('([').rstrip('])').split(',')
                        msgs.bug("list in a list could cause trouble if elements are not strings!")
                    elif '.' in i: # Might be a float
                        try: addarr += [float(i)]
                        except: addarr += [i] # Must be a string
                    else:
                        try: addarr += [int(i)] # Could be an integer
                        except: addarr += [i] # Must be a string
                tvalue = addarr
            else:
                tvalue.append(svalue)
        elif type(tvalue) is bool:
            tvalue = svalue.lower() in ['true']
        elif tvalue is None: # If it was None, it may not be anymore. We now need to find the new type
            if svalue.lower() == 'none': # None type
                tvalue = None
            elif svalue.lower() == 'true' or svalue.lower() == 'false': # bool type
                tvalue = svalue.lower() in ['true']
            elif ',' in svalue: # a list
                tvalue = svalue.split(',')
            elif '.' in svalue: # Might be a float
                try: tvalue = float(svalue)
                except: tvalue = svalue # Must be a string
            else:
                try: tvalue = int(svalue) # Could be an integer
                except: tvalue = svalue # Must be a string
        elif type(tvalue) is dict:
            nvalue = lines.split()[argnum]
            tvalue[svalue] = set_params_wtype(tvalue[svalue], nvalue, argnum=argnum+1)
        else:
            msgs.bug("Type not found for:"+msgs.newline()+lines.split('#')[0].strip())
    except:
        msgs.error(setstr + "Settings contains bad line (arg {0:d}):".format(argnum)+msgs.newline()+lines.split('#')[0].strip())
    return tvalue


def set_params(lines, indict, setstr=""):
    """
    Adjust settings parameters.
    lines    : an array of settings with the same format as the default 'settings.armed'
    indict  : a dictionary generated by initialise that contains all settings
    setstr   : a string argument for error messages that tells the user which file the error occured in.
    """
    for i in range(len(lines)):
        if lines[i].strip() == '' or lines[i].strip() == '\n': continue
        if lines[i].strip()[0] == '#': continue
        tline = lines[i].strip().split("#")[0]
        linspl = tline.split()
        if len(linspl) <= 2:
            msgs.error("Not enough parameters given on line:"+msgs.newline()+lines[i])
        if linspl[0] == 'check':
            text = str(linspl[2]).strip().replace('_', ' ')
            if ',' in text:  # There are multiple possibilities
                indict[linspl[0]][linspl[1]] += text.split(',')
            else:
                indict[linspl[0]][linspl[1]] = text
        elif linspl[0] in indict.keys():
            if linspl[1] in ['check', 'match', 'combsame']:
                text = str(linspl[3]).strip().replace('_', ' ')
                if ',' in text and text[0:2] != '%,':  # There are multiple possibilities - split the infile
                    indict[linspl[0]][linspl[1]][linspl[2]] += text.split(',')
                else:
                    indict[linspl[0]][linspl[1]][linspl[2]] = text
            elif linspl[1][:6] == 'ndet':  # Mosaic of Detectors
                indict[linspl[0]][linspl[1]] = int(linspl[2])
                tmp = []
                for ii in range(indict['mosaic']['ndet']):  # List
                    tmpi = copy.deepcopy(indict['det'])
                    tmpi['suffix'] = str(ii)
                    tmp.append(tmpi)
                indict['det'] = tmp
            elif linspl[1][:7] == 'headext':  # Header Sections
                try:
                    null = np.int(linspl[1][7:])
                except ValueError:
                    msgs.error("keyword headext must contain an integer suffix")
                indict[linspl[0]][linspl[1]] = int(linspl[2])
            elif linspl[1][:8] == 'lampname':  # Lamp names
                try:
                    null = np.int(linspl[1][8:])
                except ValueError:
                    msgs.error("keyword lampname must contain an integer suffix")
                indict[linspl[0]][linspl[1]] = linspl[2]
            elif linspl[1][:8] == 'lampstat': # Lamp status
                try:
                    null = np.int(linspl[1][8:])
                except ValueError:
                    msgs.error("keyword lampstat must contain an integer suffix")
                indict[linspl[0]][linspl[1]] = linspl[2]
            elif linspl[1] in indict[linspl[0]].keys():
                indict[linspl[0]][linspl[1]] = set_params_wtype(indict[linspl[0]][linspl[1]], linspl[2], lines=tline, setstr=setstr)
            else:
                set_trace()
                msgs.error(setstr + "Settings contains bad line (arg 2):"+msgs.newline()+lines[i].split('#')[0].strip())
        elif linspl[0][:3] == 'det': # Detector parameters
            try:
                didx = np.int(linspl[0][4:]) - 1 
            except ValueError:
                msgs.error("keyword det must contain an integer suffix")
            else:
                linspl[0] = 'det'
            if linspl[1][:6] == 'ampsec': # Amplifier Sections
                try:
                    null = np.int(linspl[1][6:])
                except ValueError:
                    msgs.error("keyword ampsec must contain an integer suffix")
                indict[linspl[0]][didx][linspl[1]] = load_sections(linspl[2], strtxt=linspl[1])
            elif linspl[1][:7] == 'datasec': # Data Sections
                try:
                    null = np.int(linspl[1][7:])
                except ValueError:
                    msgs.error("keyword datasec must contain an integer suffix")
                indict[linspl[0]][didx][linspl[1]] = load_sections(linspl[2], strtxt=linspl[1])
            elif linspl[1][:8] == 'oscansec': # Overscan Sections
                try:
                    null = np.int(linspl[1][8:])
                except ValueError:
                    msgs.error("keyword oscansec must contain an integer suffix")
                indict[linspl[0]][didx][linspl[1]] = load_sections(linspl[2], strtxt=linspl[1])
            else:  # Read value
                indict[linspl[0]][didx][linspl[1]] = set_params_wtype(indict[linspl[0]][didx][linspl[1]], linspl[2], lines=tline,setstr=setstr)
        else: 
            msgs.error(setstr + "Settings contains bad line (arg 1):"+msgs.newline()+lines[i].split('#')[0].strip())
    return indict


def load_sections(string, strtxt="<not specified>"):
    """
    From the input string, return the coordinate sections

    Parameters
    ----------
    string : str
      character string of the form [x1:x2,y1:y2]
      x1 = left pixel
      x2 = right pixel
      y1 = bottom pixel
      y2 = top pixel

    Returns
    -------
    sections : list (or None)
      the detector sections
    """
    try:
        xyrng = string.strip('[]()').split(',')
        if xyrng[0] == ":": xyarrX = [0,0]
        else: xyarrX = xyrng[0].split(':')
        if xyrng[1] == ":": xyarrY = [0,0]
        else: xyarrY = xyrng[1].split(':')
        return [[np.int(xyarrX[0]), np.int(xyarrX[1])], [np.int(xyarrY[0]), np.int(xyarrY[1])]]
    except:
        msgs.error("Keyword value {0:s} must be of the form:".format(strtxt)+msgs.newline()+"[x1:x2,y1:y2]")
    return None


def load_settings(fname, argflag):
    # Read in the default settings
    msgs.info("Loading the default settings")
    with open(fname, 'r') as infile:
        lines = infile.readlines()
    argflag = set_params(lines, argflag, setstr="Default ")
    return argflag


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
            if dfname[0] == "#":
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


def check_argflag(argflag):
    curcpu = argflag['run']['ncpus']  # Store the current number of CPUs
    # Check requested CPUs
    argflag['run']['ncpus'] = cpucheck(argflag['run']['ncpus'], curcpu=curcpu)
    # Perform some checks on the input parameters:
    if argflag['chisq']['fstep'] < 1.0: msgs.error("Setting 'fstep' in family 'chisq' must be >= 1.0")
    if argflag['out']['verbose'] not in [0, 1, 2]: msgs.error("Setting 'verbose' in family 'out' must equal 0, 1, or 2.")
    if argflag['reduce']['usebias'].lower() not in ['bias', 'overscan', 'dark', 'none']:
        if not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usebias']) and not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usebias']+".fits") and not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usebias']+".fit"):
            msgs.warn("The following file does not exist:"+msgs.newline()+"{0:s}/{1:s}".format(argflag['run']['masterdir'],argflag['reduce']['usebias']))
            msgs.error("Setting 'usebias' in family 'reduce' must be one of:"+msgs.newline()+"'bias', 'overscan', 'dark', 'none'"+msgs.newline()+"or the name of a fits file.")
    if argflag['reduce']['usetrace'].lower() not in ['trace', 'blzflat', 'science']:
        if not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usetrace']) and not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usetrace']+".fits") and not os.path.exists(argflag['run']['masterdir']+'/'+argflag['reduce']['usetrace']+".fit"):
            msgs.warn("The following file does not exist:" + msgs.newline() +
                      "{0:s}/{1:s}".format(argflag['run']['masterdir'], argflag['reduce']['usetrace']))
            msgs.error("Setting 'usetrace' in family 'reduce' must be one of:" + msgs.newline() +
                       "'trace', 'blzflat', 'science'" + msgs.newline() +
                       "or the name of a fits file.")
    # Check that the sorted data will be output in VOTable format:
    osrtspl = argflag['out']['sorted'].split('.')
    if len(osrtspl) != 1:
        if osrtspl[-1] != 'xml': msgs.error("The output format for 'sorted' is .xml, not .{0:s}".format(osrtspl[-1]))
    return


def load_spect(progname, specname, spect=None, lines=None):
    """
    Load spectrograph settings

    Parameters
    ----------
    progname : string
      Name of the program
    specname : string
      Name of spectrograph settings file
    spect : dict
      Properties of the spectrograph.
      If None, spect will be created, otherwise spect
      will be updated.
    lines : list
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
        msc = dict({'ndet': 0, 'latitude': 0.0, 'longitude': 0.0, 'elevation': 0.0, 'minexp': 0., 'reduction': 'ARMLSD'})
        # det starts as a dict but becomes a list of dicts in set_params
        ddet = dict({'xgap': 0.0, 'ygap': 0.0, 'ysize': 1.0, 'darkcurr': 0.0, 'ronoise': 1.0, 'gain': 1.0, 'saturation': 65536.0, 'nonlinear': 1.0, 'numamplifiers': 1, 'suffix': ""})
        #
        chk = dict({})
        stf = dict({'science': [], 'standard': [], 'bias': [], 'pixflat': [], 'blzflat': [], 'arc': [], 'trace': [], 'dark': []})
        kyw = dict({'target': '01.OBJECT', 'idname': '01.OBSTYPE', 'time': '01.MJD', 'date': '', 'equinox': '', 'ra': '', 'dec': '', 'airmass': '', 'naxis0': '01.NAXIS2', 'naxis1': '01.NAXIS1', 'exptime': '01.EXPTIME', 'hatch': '01.TRAPDOOR', 'filter1': '01.FILTNAME', 'filter2': None, 'lamps': '01.LAMPNAME', 'decker': '01.DECKNAME', 'slitwid': '01.SLITWIDTH', 'slitlen': '01.SLITLENGTH', 'detrot': '01.DETECTORROTATION', 'cdangle': '01.XDISPANGLE', 'echangle': '01.ECHELLEANGLE', 'crossdisp': '01.XDISPERS', 'dichroic': '', 'disperser': '', 'binning': ''})
        fts = dict({'numhead': 1, 'numlamps':1, 'dataext':0, 'calwin':12.0, 'timeunit':'mjd'})
        sci = dict({'index': [], 'check': dict({}), 'idname': 'OBJECT', 'canbe': None})
        std = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({})})
        pfl = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        bfl = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        trc = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        arc = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 1, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}), 'lscomb': False})
        bia = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}) })
        drk = dict({'index': [], 'check': dict({}), 'match': dict({}), 'number': 5, 'idname': 'OBJECT', 'canbe': None, 'combsame': dict({}) })
        spectt = dict({'mosaic': msc, 'det': ddet, 'check': chk, 'set': stf, 'keyword': kyw, 'fits': fts, 'science': sci, 'standard': std, 'pixflat': pfl, 'blzflat': bfl, 'trace': trc, 'arc': arc, 'bias': bia, 'dark': drk})
        return spectt
    if lines is None:
        # Read in the default settings
        # Get the software path
        prgn_spl = progname.split('/')
        fname = ""
        for i in range(0, len(prgn_spl)-1): fname += prgn_spl[i]+"/"
        fname += 'settings.'+specname
        msgs.info("Loading the "+specname+" settings")
        spect = initialise()
        with open(fname, 'r') as infile:
            lines = infile.readlines()
        spect = set_params(lines, spect, setstr="Default "+specname+" ")
    else:
        if spect is not None:
            spect = set_params(lines, spect, setstr="Infile "+specname+" ")
    return spect


def load_headers(argflag, spect, datlines):
    """
    Load the header information for each fits file

    Parameters
    ----------
    argflag : dict
      Arguments and flags used for reduction
    spect : dict
      Properties of the spectrograph.
      If None, spect will be created, otherwise spect
      will be updated.
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
    for k in keys: fitsdict[k]=[]
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
        for ch in chks:
            tfrhd = int(ch.split('.')[0])-1
            kchk  = '.'.join(ch.split('.')[1:])
            frhd  = whddict['{0:02d}'.format(tfrhd)]
            if spect['check'][ch] != str(headarr[frhd][kchk]).strip():
                #set_trace()
                #print ch, frhd, kchk
                #print spect['check'][ch], str(headarr[frhd][kchk]).strip()
                msgs.error("The following file:"+msgs.newline()+datlines[i]+msgs.newline()+"is not taken with the settings.{0:s} detector".format(argflag['run']['spectrograph'])+msgs.newline()+"Remove this file, or specify a different settings file.")
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
            if spect['keyword'][kw] is None: value = 'None'  # This instrument doesn't have/need this keyword
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
                        value='None'
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
            elif typv is str or typv is np.string_:
                fitsdict[kw].append(value.strip())
            else:
                msgs.bug("I didn't expect useful headers to contain type {0:s}".format(typv).replace('<type ','').replace('>',''))

        if argflag['out']['verbose'] == 2: msgs.info("Successfully loaded headers for file:"+msgs.newline()+datlines[i])
    del headarr
    # Convert the fitsdict arrays into numpy arrays
    for k in fitsdict.keys(): fitsdict[k] = np.array(fitsdict[k])
    msgs.info("Headers loaded for {0:d} files successfully".format(len(datlines)))
    return fitsdict


def load_frames(slf, fitsdict, ind, det, frametype='<None>', msbias=None, trim=True, transpose=False):
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
        if slf._argflag['run']['spectrograph'] in ['lris_blue']:
#            set_trace()
            temp, head0, _ = arlris.read_lris(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], det=det)
        else:
            temp = pyfits.getdata(fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], slf._spect['fits']['dataext'])
        temp = temp.astype(float)  # Let us avoid uint16
        if transpose: temp = temp.T
        if msbias is not None:
            if type(msbias) is np.ndarray:
                temp -= msbias  # Subtract the master bias frame
            elif type(msbias) is str:
                if msbias == "overscan":
                    arproc.sub_overscan(slf, det, temp)
                else:
                    msgs.error("Could not subtract bias level when loading {0:s} frames".format(frametype))
            if trim: 
                temp = arproc.trim(slf, temp, det)
        if i == 0:
            frames = np.zeros((temp.shape[0], temp.shape[1], np.size(ind)))
            frames[:,:,i] = temp.copy()
        else:
            frames[:,:,i] = temp.copy()
        del temp
#	pool = mpPool(processes=np.min([slf._argflag['run']['ncpus'],np.size(ind)]))
#	async_results = []
#	for i in range(np.size(ind)):
#		async_results.append(pool.apply_async(pyfits.getdata, (fitsdict['directory'][ind[i]]+fitsdict['filename'][ind[i]], slf._spect['fits']['dataext'])))
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


def load_master(name, frametype='<None>'):
    """
    Load a pre-existing master calibration frame

    Parameters
    ----------
    name : str
      Name of the master calibration file to be loaded
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
            infile = pyfits.open(name)
        except:
            msgs.error("Master calibration file does not exist:"+msgs.newline()+name)
        msgs.info("Master {0:s} frame loaded successfully:".format(infile[0].header['FRAMETYP'])+msgs.newline()+name)
        return np.array(infile[0].data, dtype=np.float)
    else:
        msgs.info("Loading Master {0:s} frame:".format(frametype)+msgs.newline()+name)
        return np.array(pyfits.getdata(name, 0), dtype=np.float)


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