import os
import re
import sys
import shutil
import numpy as np
import armsgs as msgs
import arutils
import arcyutils
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
#from scipy.stats import chi2 as chisq

from xastropy.xutils import xdebug as xdb

def sort_data(slf):
    from arvcorr import radec_to_decdeg
    msgs.bug("There appears to be a bug with the assignment of arc frames when only one science frame is supplied")
    msgs.info("Sorting files")
    numfiles = slf._fitsdict['filename'].size
    # A dictionary of filetypes
    ftag = dict({'science':np.array([], dtype=np.int),
                 'standard':np.array([], dtype=np.int),
                 'bias':np.array([], dtype=np.int),
                 'pixflat':np.array([], dtype=np.int),
                 'blzflat':np.array([], dtype=np.int),
                 'trace':np.array([], dtype=np.int),
                 'arc':np.array([], dtype=np.int)})
    fkey = np.array(ftag.keys())
    # Create an array where 1 means it is a certain type of frame and 0 means it isn't.
    filarr = np.zeros((len(fkey),numfiles), dtype=np.int)
    # Identify the frames:
    for i in range(len(fkey)):
        w = np.where(slf._fitsdict['idname']==slf._spect[fkey[i]]['idname'])[0]
        n = np.arange(numfiles)
        n = np.intersect1d(n,w)
        # Perform additional checks in order to make sure this identification is true
        chkk = slf._spect[fkey[i]]['check'].keys()
        xdb.set_trace()
        for ch in chkk:
            if ch[0:9]=='condition':
                # Deal with a conditional argument
                conds = re.split("(\||\&)",slf._spect[fkey[i]]['check'][ch])
                tcond = conds[0].split("=")
                ntmp = (slf._fitsdict[tcond[0]]==tcond[1])
                for cn in range((len(conds)-1)/2):
                    tcond = conds[2*cn+2].split("=")
                    if conds[2*cn+1]=="|":
                        ntmp = ntmp | (slf._fitsdict[tcond[0]]==tcond[1])
                    elif conds[2*cn+1]=="&":
                        ntmp = ntmp & (slf._fitsdict[tcond[0]]==tcond[1])
                w = np.where(ntmp)[0]
            else:
                w = np.where(slf._fitsdict[ch]==slf._spect[fkey[i]]['check'][ch])[0]
            n = np.intersect1d(n,w)
        # Assign these filetypes
        filarr[i,:][n] = 1
        # Check if these files can also be another type
        if slf._spect[fkey[i]]['canbe'] is not None:
            for cb in slf._spect[fkey[i]]['canbe']:
                # Assign these filetypes
                fa = np.where(fkey==cb)[0]
                if np.size(fa) == 1: filarr[fa[0],:][n] = 1
                else: msgs.error("Unknown type for argument 'canbe': {0:s}".format(cb))
#		# Check for filetype clashes
#		bdf=np.where(np.sum(filarr,axis=0)[n] != 0)[0]
#		if np.size(bdf) != 0:
#			msgs.info("Clash with file identifications when assigning frames as type {0:s}:".format(fkey[i]))
#			# Check if this clash is allowed:
#			clashfound=False
#			for b in range(np.size(bdf)):
#				# For each file with a clash, get all frames that have been assigned to it
#				tarr = np.where(filarr[:,n[bdf[b]]]==1)[0]
#				for a in tarr:
#					if slf._spect[fkey[i]]['canbe'] is None:
#						print "{0:s} can only be of type: {1:s}, and is marked as {2:s}".format(slf._fitsdict['filename'][n[bdf[b]]],fkey[i],fkey[a])
#						clashfound=True
#					elif fkey[a] not in slf._spect[fkey[i]]['canbe']:
#						print "{0:s}  current filetype: {1:s}".format(slf._fitsdict['filename'][n[bdf[b]]],fkey[a])
#						clashfound=True
#			if clashfound: msgs.error("Check these files and your settings.{0:s} file before continuing.".format(slf._argflag['run']['spectrograph'])+msgs.newline()+"You can use the 'canbe' option to allow one frame to have multiple types.")
#			else: msgs.info("Clash permitted")
    # Identify the standard stars
    # Find the nearest standard star to each science frame
    wscistd = np.where(filarr[np.where(fkey=='standard')[0],:].flatten() == 1)[0]
    for i in range(wscistd.size):
        raval, decval = radec_to_decdeg(slf._fitsdict['ra'][wscistd[i]], slf._fitsdict['dec'][wscistd[i]])
        offset = arutils.calc_offset(15.0*raval, decval, slf._standardStars["RA"], slf._standardStars["DEC"], distance=True)
        # If an object exists within 1 arcmin of a listed standard, then it must be a standard star
        if (np.min(offset) < 60.0):
            filarr[np.where(fkey=='science')[0],wscistd[i]]=0
        else:
            filarr[np.where(fkey=='standard')[0],wscistd[i]]=0
    # Check that all files have an identification
    badfiles=np.where(np.sum(filarr,axis=0) == 0)[0]
    if np.size(badfiles) != 0:
        msgs.info("Couldn't identify the following files:")
        for i in range(np.size(badfiles)): print slf._fitsdict['filename'][badfiles[i]]
        msgs.error("Check these files and your settings.{0:s} file before continuing".format(slf._argflag['run']['spectrograph']))
    # Now identify the dark frames
    wdark = np.where((filarr[np.where(fkey=='bias')[0],:]==1).flatten() & (slf._fitsdict['exptime'].astype(np.float64) > 0.0))[0]
    ftag['dark'] = wdark
    # Make any forced changes
    msgs.info("Making forced file identification changes")
    skeys = slf._spect['set'].keys()
    for sk in skeys:
        for j in slf._spect['set'][sk]:
            w = np.where(slf._fitsdict['filename']==j)[0]
            filarr[:,w]=0
            filarr[np.where(fkey==sk)[0],w]=1
            del w
    # Store the frames in the ftag array
    for i in range(len(fkey)):
        ftag[fkey[i]] = np.where(filarr[i,:]==1)[0]
    # Finally check there are no duplicates (the arrays will automatically sort with np.unique)
    msgs.info("Finalising frame sorting, and removing duplicates")
    for key in ftag.keys():
        ftag[key] = np.unique(ftag[key])
        if np.size(ftag[key])==1:
            msgs.info("Found {0:d} {1:s} frame".format(np.size(ftag[key]),key))
        else:
            msgs.info("Found {0:d} {1:s} frames".format(np.size(ftag[key]),key))
    # Return ftag!
    msgs.info("Sorting completed successfully")
    return ftag

def sort_write(slf,space=3):
    """
    Write out an ascii file that contains the details of the sorting.
    By default, the filename is printed first, followed by the filetype.
    After these, all parameters listed in the 'keyword' item in the
    settings file will be printed
    -----------------------
    space : keyword to set how many blank spaces to place between keywords.
    """
    msgs.info("Preparing to write out the data sorting details")
    nfiles=slf._fitsdict['filename'].size
    # Specify which keywords to print after 'filename' and 'filetype'
    prord = ['filename', 'frametype', 'target', 'exptime', 'naxis0', 'naxis1', 'filter1', 'filter2']
    prdtp = ["char",     "char",      "char",   "double",  "int",    "int",    "char",     "char"]
    # Now insert the remaining keywords:
    fkey = slf._spect['keyword'].keys()
    for i in fkey:
        if i not in prord:
            prord.append(i)
            # Append the type of value this keyword holds
            typv = type(slf._fitsdict[i][0])
            if typv is int or typv is np.int_:
                prdtp.append("int")
            elif typv is str or typv is np.string_:
                prdtp.append("char")
            elif typv is float or typv is np.float_:
                prdtp.append("double")
            else:
                msgs.bug("I didn't expect useful headers to contain type {0:s}".format(typv).replace('<type ','').replace('>',''))
    # Open a VOTable for writing
    votable = VOTableFile()
    resource = Resource()
    votable.resources.append(resource)
    table = Table(votable)
    resource.tables.append(table)
    # Define VOTable fields
    tabarr=[]
    # Insert the filename and filetype first
    for i in range(len(prord)): tabarr.append(Field(votable, name=prord[i], datatype=prdtp[i], arraysize="*"))
    table.fields.extend(tabarr)
    table.create_arrays(nfiles)
    filtyp = slf._filesort.keys()
    for i in range(nfiles):
        values = ()
        for pr in prord:
            if pr=='frametype':
                addval = ""
                for ft in filtyp:
                    if i in slf._filesort[ft]:
                        if len(addval) != 0: addval += ","
                        addval += ft
                addval = (addval,)
            else: addval = (slf._fitsdict[pr][i],)
            values = values + addval
        table.array[i] = values
    osspl = slf._argflag['out']['sorted'].split('.')
    if len(osspl) > 1:
        fname = slf._argflag['out']['sorted']
    else:
        fname = slf._argflag['out']['sorted']+'.xml'
    votable.to_xml(fname)
    msgs.info("Successfully written sorted data information file:"+msgs.newline()+"{0:s}".format(fname))
    return

def match_science(slf):
    """
    For a given set of identified data, match frames to science frames
    """
    msgs.info("Matching calibrations to Science frames")
    ftag = ['standard','bias','dark','pixflat','blzflat','trace','arc']
    nfiles = slf._fitsdict['filename'].size
    iSCI = slf._filesort['science']
    iSTD = slf._filesort['standard']
    iBIA = slf._filesort['bias']
    iDRK = slf._filesort['dark']
    iPFL = slf._filesort['pixflat']
    iBFL = slf._filesort['blzflat']
    iTRC = slf._filesort['trace']
    iARC = slf._filesort['arc']
    iARR = [iSTD,iBIA,iDRK,iPFL,iBFL,iTRC,iARC]
    nSCI = iSCI.size
#	nBIA = iBIA.size
#	nDRK = iDRK.size
#	nFLT = iFLT.size
#	nTRC = iTRC.size
#	nARC = iARC.size
    i=0
    while i<nSCI:
        msgs.info("Matching calibrations to {0:s}".format(slf._fitsdict['target'][iSCI[i]]))
        slf._spect['science']['index'].append(np.array([iSCI[i]]))
        #find nearby calibration frames
        for ft in range(len(ftag)):
            # Some checks first to make sure we need to find matching frames
            if ftag[ft] == 'dark' and slf._argflag['reduce']['usebias'] != 'dark':
                msgs.info("  Dark frames not required")
                continue
            if ftag[ft] == 'bias' and slf._argflag['reduce']['usebias'] != 'bias' and not slf._argflag['reduce']['badpix']:
                msgs.info("  Bias frames not required")
                continue
            # Now go ahead and match the frames
            n = np.arange(nfiles)
            chkk = slf._spect[ftag[ft]]['match'].keys()
            for ch in chkk:
                tmtch = slf._spect[ftag[ft]]['match'][ch]
                if tmtch == "''":
                    w = np.where(slf._fitsdict[ch]==slf._fitsdict[ch][iSCI[i]])[0]
                elif tmtch[0] == '=':
                    mtch = np.float64(slf._fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                    w = np.where((slf._fitsdict[ch]).astype(np.float64)==mtch)[0]
                elif tmtch[0] == '<':
                    if tmtch[1] == '=':
                        mtch = np.float64(slf._fitsdict[ch][iSCI[i]]) + np.float64(tmtch[2:])
                        w = np.where((slf._fitsdict[ch]).astype(np.float64)<=mtch)[0]
                    else:
                        mtch = np.float64(slf._fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                        w = np.where((slf._fitsdict[ch]).astype(np.float64)<mtch)[0]
                elif tmtch[0] == '>':
                    if tmtch[1] == '=':
                        mtch = np.float64(slf._fitsdict[ch][iSCI[i]]) + np.float64(tmtch[2:])
                        w = np.where((slf._fitsdict[ch]).astype(np.float64)>=mtch)[0]
                    else:
                        mtch = np.float64(slf._fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                        w = np.where((slf._fitsdict[ch]).astype(np.float64)>mtch)[0]
                elif tmtch[0] == '|':
                    if tmtch[1] == '=':
                        mtch = np.float64(tmtch[2:])
                        w = np.where(np.abs((slf._fitsdict[ch]).astype(np.float64)-np.float64(slf._fitsdict[ch][iSCI[i]]))==mtch)[0]
                    elif tmtch[1] == '<':
                        if tmtch[2] == '=':
                            mtch = np.float64(tmtch[3:])
                            w = np.where(np.abs((slf._fitsdict[ch]).astype(np.float64)-np.float64(slf._fitsdict[ch][iSCI[i]]))<=mtch)[0]
                        else:
                            mtch = np.float64(tmtch[2:])
                            w = np.where(np.abs((slf._fitsdict[ch]).astype(np.float64)-np.float64(slf._fitsdict[ch][iSCI[i]]))<mtch)[0]
                    elif tmtch[1] == '>':
                        if tmtch[2] == '=':
                            mtch = np.float64(tmtch[3:])
                            w = np.where(np.abs((slf._fitsdict[ch]).astype(np.float64)-np.float64(slf._fitsdict[ch][iSCI[i]]))>=mtch)[0]
                        else:
                            mtch = np.float64(tmtch[2:])
                            w = np.where(np.abs((slf._fitsdict[ch]).astype(np.float64)-np.float64(slf._fitsdict[ch][iSCI[i]]))>mtch)[0]
                elif tmtch[0:2] == '%,': # Splitting a header keyword
                    splcom = tmtch.split(',')
                    try:
                        spltxt, argtxt, valtxt = splcom[1], np.int(splcom[2]), splcom[3]
                        tspl=[]
                        for sp in slf._fitsdict[ch]:
                            tspl.append(sp.split(spltxt)[argtxt])
                        tspl = np.array(tspl)
                        if valtxt == "''":
                            w = np.where(tspl==slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt])[0]
                        elif valtxt[0] == '=':
                            mtch = np.float64(slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt]) + np.float64(valtxt[1:])
                            w = np.where((tspl).astype(np.float64)==mtch)[0]
                        elif valtxt[0] == '<':
                            if valtxt[1] == '=':
                                mtch = np.float64(slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt]) + np.float64(valtxt[2:])
                                w = np.where((tspl).astype(np.float64)<=mtch)[0]
                            else:
                                mtch = np.float64(slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt]) + np.float64(valtxt[1:])
                                w = np.where((tspl).astype(np.float64)<mtch)[0]
                        elif valtxt[0] == '>':
                            if valtxt[1] == '=':
                                mtch = np.float64(slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt]) + np.float64(valtxt[2:])
                                w = np.where((tspl).astype(np.float64)>=mtch)[0]
                            else:
                                mtch = np.float64(slf._fitsdict[ch][iSCI[i]].split(spltxt)[argtxt]) + np.float64(valtxt[1:])
                                w = np.where((tspl).astype(np.float64)>mtch)[0]
                    except:
                        msgs.error("Bad form for matching criteria: {0:s}".format(tmtch))
                else:
                    msgs.bug("Matching criteria {0:s} is not supported".format(tmtch))
                n = np.intersect1d(n,w) # n corresponds to all frames with matching instrument setup to science frames
            # Find the time difference between the calibrations and science frames
            if slf._spect['fits']['calwin'] > 0.0:
                tdiff = np.abs(slf._fitsdict['time'][n].astype(np.float64)-np.float64(slf._fitsdict['time'][iSCI[i]]))
                w = np.where(tdiff<=slf._spect['fits']['calwin'])[0]
                n = np.intersect1d(n,w) # n corresponds to all frames with matching instrument setup
            # Now find which of the remaining n are the appropriate calibration frames
            n = np.intersect1d(n,iARR[ft])
            # How many frames are required
            numfr = slf._spect[ftag[ft]]['number']
            if slf._argflag['out']['verbose'] == 2:
                if numfr==1: areis = "is"
                else: areis = "are"
                if np.size(n) == 1:
                    msgs.info("  Found {0:d} {1:s} frame for {2:s} ({3:d} {4:s} required)".format(np.size(n),ftag[ft],slf._fitsdict['target'][iSCI[i]],numfr,areis))
                else:
                    msgs.info("  Found {0:d} {1:s} frames for {2:s} ({3:d} {4:s} required)".format(np.size(n),ftag[ft],slf._fitsdict['target'][iSCI[i]],numfr,areis))
            # Have we identified enough of these calibration frames to continue?
            if np.size(n) < numfr:
                msgs.warn("  Only {0:d}/{1:d} {2:s} frames for {3:s}".format(np.size(n),numfr,ftag[ft],slf._fitsdict['target'][iSCI[i]]))
                # Errors for insufficient BIAS frames
                if slf._argflag['reduce']['usebias'].lower() == ftag[ft]:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient PIXELFLAT frames
                if ftag[ft] == 'pixflat' and slf._argflag['reduce']['flatfield']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient BLAZEFLAT frames
                if ftag[ft] == 'blzflat' and slf._argflag['reduce']['flatfield']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient TRACE frames
                if ftag[ft] == 'trace':
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient ARC frames
                if ftag[ft] == 'arc' and slf._argflag['reduce']['calibrate']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient ARC frames
                if ftag[ft] == 'standard' and slf._argflag['reduce']['fluxcalibrate']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
            else:
                # Select the closest calibration frames to the science frame
                tdiff = np.abs(slf._fitsdict['time'][n].astype(np.float64)-np.float64(slf._fitsdict['time'][iSCI[i]]))
                wa = np.argsort(tdiff)
                slf._spect[ftag[ft]]['index'].append(n[wa[:numfr]])
        i+=1
    msgs.info("Science frames successfully matched to calibration frames")
    return


def match_frames(slf, frames, criteria, frametype='<None>', satlevel=None):
    """
    identify frames with a similar appearance (i.e. one frame appears to be a scaled version of another).
    """
    prob  = arutils.erf(criteria/np.sqrt(2.0))[0]
    frsh0, frsh1, frsh2 = frames.shape
    msgs.info("Matching {0:d} {1:s} frames with confidence interval {2:5.3f}%".format(frsh2,frametype,prob*100.0))
    srtframes = [np.zeros((frsh0,frsh1,1))]
    srtframes[0][:,:,0] = frames[:,:,0]
    tsrta = [frames[frsh0/2,:,0]]
    tsrtb = [frames[:,frsh1/2,0]]
    msgs.bug("Throughout this routine, you should probably search for the mean of the non-saturated pixels")
    tsrta[0] /= np.mean(tsrta[0])
    tsrtb[0] /= np.mean(tsrtb[0])
    for fr in range(1,frames.shape[2]):
        fm = None
        for st in range(len(srtframes)):
            tmata = frames[frsh0/2,:,fr]
            tmatb = frames[:,frsh1/2,fr]
            tmata /= np.mean(tmata)
            tmatb /= np.mean(tmatb)
            if satlevel is None:
                wa = np.where(tmata>0.0)
                wb = np.where(tmatb>0.0)
            else:
                wa = np.where((tmata>0.0)&(tmata<satlevel))
                wb = np.where((tmatb>0.0)&(tmatb<satlevel))
            testa, testb = np.mean(tsrta[st][wa]/tmata[wa]), np.mean(tsrtb[st][wb]/tmatb[wb])
            if np.size(wa[0])==0 or np.size(wb[0])==0:
                msgs.bug("I didn't expect to find a row of zeros in the middle of the chip!")
                sys.exit()
            if (testa>=prob) and (testa<=(2.0-prob)) and (testb>=prob) and (testb<=(2.0-prob)):
                fm = st
                break
        if fm is None:
            srtframes.append(np.zeros((frames.shape[0],frames.shape[1],1)))
            srtframes[-1][:,:,0] = frames[:,:,fr]
            tsrta.append(tmata)
            tsrtb.append(tmatb)
        else:
            srtframes[fm] = np.append(srtframes[fm],np.zeros((frames.shape[0],frames.shape[1],1)),axis=2)
            srtframes[fm][:,:,-1] = frames[:,:,fr]
    if len(srtframes) > 1:
        msgs.info("Found {0:d} different sets of {1:s} frames".format(len(srtframes),frametype))
    else:
        msgs.info("Found {0:d} set of {1:s} frames".format(len(srtframes),frametype))
    if frames.shape[2] > 1: del tsrta, tsrtb, tmata, tmatb, testa, testb
    return srtframes

def match_frames_old(slf, frames, frametype='<None>'):
    msgs.info("Matching {0:d} {1:s} frames".format(frames.shape[2],frametype))
    srtframes = [np.zeros((frames.shape[0],frames.shape[1],1))]
    srtframes[0][:,:,0] = frames[:,:,0]
    prob  = arutils.erf(slf._argflag['reduce']['flatmatch']/np.sqrt(2.0))[0]
#	chisqv = chisq.ppf(prob,frames.shape[0]*frames.shape[1])
    chisqv = frames.shape[0]*frames.shape[1]
    for fr in range(1,frames.shape[2]):
        fm = None
        for st in range(len(srtframes)):
            chisqc = arcyutils.checkmatch(srtframes[st][:,:,0],frames[:,:,fr],1048577.0)
            if chisqc < chisqv*10.0:
                fm = st
                break
        if fm is None:
            srtframes.append(np.zeros((frames.shape[0],frames.shape[1],1)))
            srtframes[-1][:,:,0] = frames[:,:,fr]
        else:
            srtframes[fm] = np.append(srtframes[fm],np.zeros((frames.shape[0],frames.shape[1],1)),axis=2)
            srtframes[fm][:,:,-1] = frames[:,:,fr]
    msgs.info("Found {0:d} different sets of {1:s} frames".format(len(srtframes),frametype))
    return srtframes

def make_dirs(slf):
    # First, get the current working directory
    currDIR=os.getcwd()
    msgs.info("Creating Science directory")
    newdir = "{0:s}/{1:s}".format(currDIR,slf._argflag['run']['scidir'])
    if os.path.exists(newdir):
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
        if not slf._argflag['out']['overwrite']:
            rmdir=''
            while os.path.exists(newdir):
                while rmdir != 'n' and rmdir != 'y' and rmdir != 'r':
                    rmdir=raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, [r]ename) - ")
                if rmdir == 'n':
                    msgs.warn("Any previous calibration files may be overwritten")
                    break
                elif rmdir == 'r':
                    newdir=raw_input(msgs.input()+"Enter a new directory name: ")
                elif rmdir == 'y':
                    shutil.rmtree(newdir)
                    os.mkdir(newdir)
                    break
            if rmdir == 'r': os.mkdir(newdir)
    else: os.mkdir(newdir)
    # Create a directory for each object in the Science directory
    msgs.info("Creating Object directories")
    #Go through objects creating directory tree structure
    w=slf._filesort['science']
    sci_targs = np.array(list(set(slf._fitsdict['target'][w])))
    # Loop through targets and replace spaces with underscores
    nored=np.array([])
    # Create directories
    rmalways = False
    for i in range(sci_targs.size):
        sci_targs[i] = sci_targs[i].replace(' ', '_')
        newdir = "{0:s}/{1:s}/{2:s}".format(currDIR,slf._argflag['run']['scidir'],sci_targs[i])
        if os.path.exists(newdir):
            if slf._argflag['out']['overwrite'] or rmalways:
                pass
#				shutil.rmtree(newdir)
#				os.mkdir(newdir)
            else:
                msgs.info("The following directory already exists:"+msgs.newline()+newdir)
                rmdir=''
                while rmdir != 'n' and rmdir != 'y' and rmdir != 'a':
                    rmdir=raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, or [a]lways) - ")
                if rmdir == 'n':
                    msgs.info("Not reducing {0:s}".format(sci_targs[i]))
                    nored = np.append(i)
                else:
                    shutil.rmtree(newdir)
                    os.mkdir(newdir)
                    if rmdir == 'a': rmalways=True
        else: os.mkdir(newdir)
    # Remove the entries from sci_targs which will not be reduced
    nored = nored.astype(np.int)
    while nored.size > 0:
        sci_targs = np.delete(sci_targs, nored[0])
        nored = np.delete(nored, 0)
    # Create a directory where all of the master calibration frames are stored.
    msgs.info("Creating Master Calibrations directory")
    newdir = "{0:s}/{1:s}".format(currDIR,slf._argflag['run']['masterdir'])
    if os.path.exists(newdir):
        if not slf._argflag['out']['overwrite']:
            msgs.info("The following directory already exists:"+msgs.newline()+newdir)
            rmdir=''
            while rmdir != 'n' and rmdir != 'y':
                rmdir=raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o) - ")
            if rmdir == 'n':
                msgs.warn("Any previous calibration files will be overwritten")
            else:
                shutil.rmtree(newdir)
                os.mkdir(newdir)
#		else:
#			shutil.rmtree(newdir)
#			os.mkdir(newdir)
    else: os.mkdir(newdir)
    # Create a directory where all of the master calibration frames are stored.
    msgs.info("Creating Plots directory")
    newdir = "{0:s}/{1:s}".format(currDIR,slf._argflag['run']['plotsdir'])
    if os.path.exists(newdir):
        if not slf._argflag['out']['overwrite']:
            msgs.info("The following directory already exists:"+msgs.newline()+newdir)
            rmdir=''
            while rmdir != 'n' and rmdir != 'y':
                rmdir=raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o) - ")
            if rmdir == 'n':
                msgs.warn("Any previously made plots will be overwritten")
            else:
                shutil.rmtree(newdir)
                os.mkdir(newdir)
        else:
            shutil.rmtree(newdir)
            os.mkdir(newdir)
    else: os.mkdir(newdir)
    # Return the name of the science targets
    return sci_targs

