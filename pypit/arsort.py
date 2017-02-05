from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import re
import sys
import shutil
import string
import numpy as np
import yaml

from pypit import armsgs
from pypit import arparse as settings
from pypit import arutils
from pypit.arflux import find_standard_file
from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
from astropy.table import Table as tTable, Column
from astropy import units as u

from linetools import utils as ltu

from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

# Logging
msgs = armsgs.get_logger()


def sort_data(fitsdict, flag_unknown=False):
    """ Generate a dict of filetypes from the input fitsdict object

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    flag_unknown : bool, optional
      Instead of crashing out if there are unidentified files,
      set to 'unknown' and continue

    Returns
    -------
    ftag : dict
      A dictionary of filetypes
    """
    msgs.bug("There appears to be a bug with the assignment of arc frames when only one science frame is supplied")
    msgs.info("Sorting files")
    numfiles = fitsdict['filename'].size
    # Set the filetype dictionary
    ftag = dict({'science': np.array([], dtype=np.int),
                 'standard': np.array([], dtype=np.int),
                 'bias': np.array([], dtype=np.int),
                 'dark': np.array([], dtype=np.int),
                 'pinhole': np.array([], dtype=np.int),
                 'pixelflat': np.array([], dtype=np.int),
                 'trace': np.array([], dtype=np.int),
                 'unknown': np.array([], dtype=np.int),
                 'arc': np.array([], dtype=np.int)})
    if len(settings.ftdict) > 0:
        for ifile,ftypes in settings.ftdict.items():
            idx = np.where(fitsdict['filename'] == ifile)[0]
            sptypes = ftypes.split(',')
            for iftype in sptypes:
                ftag[iftype] = np.concatenate([ftag[iftype], idx])
        # Sort
        for key in ftag.keys():
            ftag[key].sort()
        return ftag
    fkey = np.array(list(ftag.keys()))  # For Python 3 compatability
    # Create an array where 1 means it is a certain type of frame and 0 means it isn't.
    filarr = np.zeros((len(fkey), numfiles), dtype=np.int)
    setarr = np.zeros((len(fkey), numfiles), dtype=np.int)
    # Identify the frames:
    for i in range(len(fkey)):
        if fkey[i] == 'unknown':
            continue
        # Self identification
        if settings.argflag['run']['useIDname']:
            w = np.where(fitsdict['idname'] == settings.spect[fkey[i]]['idname'])[0]
            msgs.info("Sorting files")
        else:
            w = np.arange(numfiles)
        n = np.arange(numfiles)
        n = np.intersect1d(n, w)
        # Perform additional checks in order to make sure this identification is true
        if 'check' in settings.spect[fkey[i]].keys():
            chkk = settings.spect[fkey[i]]['check'].keys()
            for ch in chkk:
                if ch[0:9] == 'condition':
                    # Deal with a conditional argument
                    conds = re.split("(\||\&)", settings.spect[fkey[i]]['check'][ch])
                    ntmp = chk_condition(fitsdict, conds[0])
                    # And more
                    for cn in range((len(conds)-1)//2):
                        if conds[2*cn+1] == "|":
                            ntmp = ntmp | chk_condition(fitsdict, conds[2*cn+2])
                        elif conds[2*cn+1] == "&":
                            ntmp = ntmp & chk_condition(fitsdict, conds[2*cn+2])
                    w = np.where(ntmp)[0]
                else:
                    if fitsdict[ch].dtype.char == 'S':  # Numpy string array
                        # Strip numpy string array of all whitespace
                        w = np.where(np.char.strip(fitsdict[ch]) == settings.spect[fkey[i]]['check'][ch])[0]
                    else:
                        w = np.where(fitsdict[ch] == settings.spect[fkey[i]]['check'][ch])[0]
                n = np.intersect1d(n, w)
        # Assign these filetypes
        filarr[i, :][n] = 1
        # Check if these files can also be another type
        if settings.spect[fkey[i]]['canbe'] is not None:
            for cb in settings.spect[fkey[i]]['canbe']:
                # Assign these filetypes
                fa = np.where(fkey == cb)[0]
                if np.size(fa) == 1:
                    filarr[fa[0], :][n] = 1
                else:
                    msgs.error("Unknown type for argument 'canbe': {0:s}".format(cb))
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
#					if settings.spect[fkey[i]]['canbe'] is None:
#						print "{0:s} can only be of type: {1:s}, and is marked as {2:s}".format(fitsdict['filename'][n[bdf[b]]],fkey[i],fkey[a])
#						clashfound=True
#					elif fkey[a] not in settings.spect[fkey[i]]['canbe']:
#						print "{0:s}  current filetype: {1:s}".format(fitsdict['filename'][n[bdf[b]]],fkey[a])
#						clashfound=True
#			if clashfound: msgs.error("Check these files and your settings.{0:s} file before continuing.".format(settings.argflag['run']['spectrograph'])+msgs.newline()+"You can use the 'canbe' option to allow one frame to have multiple types.")
#			else: msgs.info("Clash permitted")
    # Identify the standard stars
    # Find the nearest standard star to each science frame
    wscistd = np.where(filarr[np.where(fkey == 'standard')[0], :].flatten() == 1)[0]
    for i in range(wscistd.size):
        radec = (fitsdict['ra'][wscistd[i]], fitsdict['dec'][wscistd[i]])
        if fitsdict['ra'][wscistd[i]] == 'None':
            msgs.warn("No RA and DEC information for file:" + msgs.newline() + fitsdict['filename'][wscistd[i]])
            msgs.warn("The above file could be a twilight flat frame that was" + msgs.newline() +
                      "missed by the automatic identification.")
            filarr[np.where(fkey == 'standard')[0], wscistd[i]] = 0
            continue
        # If an object exists within 20 arcmins of a listed standard, then it is probably a standard star
        foundstd = find_standard_file(radec, toler=20.*u.arcmin, check=True)
        if foundstd:
            filarr[np.where(fkey == 'science')[0], wscistd[i]] = 0
        else:
            filarr[np.where(fkey == 'standard')[0], wscistd[i]] = 0
    # Make any forced changes
    msgs.info("Making forced file identification changes")
    skeys = settings.spect['set'].keys()
    for sk in skeys:
        for j in settings.spect['set'][sk]:
            w = np.where(fitsdict['filename']==j)[0]
            filarr[:,w]=0
            setarr[np.where(fkey==sk)[0],w]=1
            del w
    filarr = filarr + setarr
    # Check that all files have an identification
    badfiles = np.where(np.sum(filarr, axis=0) == 0)[0]
    if np.size(badfiles) != 0:
        msgs.info("Couldn't identify the following files:")
        for i in range(np.size(badfiles)):
            msgs.info(fitsdict['filename'][badfiles[i]])
        if flag_unknown:
            filarr[np.where(fkey == 'unknown')[0],badfiles] = 1
        else:
            msgs.error("Check these files and your settings.{0:s} file before continuing".format(settings.argflag['run']['spectrograph']))
    # Now identify the dark frames
    wdark = np.where((filarr[np.where(fkey == 'bias')[0],:] == 1).flatten() &
        (fitsdict['exptime'].astype(np.float64) > settings.spect['mosaic']['minexp']))[0]
    ftag['dark'] = wdark
    # Store the frames in the ftag array
    for i in range(len(fkey)):
        ftag[fkey[i]] = np.where(filarr[i,:] == 1)[0]
    # Finally check there are no duplicates (the arrays will automatically sort with np.unique)
    msgs.info("Finalising frame sorting, and removing duplicates")
    for key in ftag.keys():
        ftag[key] = np.unique(ftag[key])
        if np.size(ftag[key]) == 1:
            msgs.info("Found {0:d} {1:s} frame".format(np.size(ftag[key]),key))
        else:
            msgs.info("Found {0:d} {1:s} frames".format(np.size(ftag[key]),key))
    # Return ftag!
    msgs.info("Sorting completed successfully")
    return ftag


def chk_condition(fitsdict, cond):
    """
    Code to perform condition.  A bit messy so a separate definition
    was generated.
    Create an exposure class for every science frame

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    cond : str
      A user-specified condition that is used to identify filetypes.
      This string is the fourth argument of the frame conditions that
      is specified in the settings file. For example, in the line:
      'bias check condition1 exptime=0'
      cond = 'exptime=0'

    Returns
    -------
    ntmp: bool array
      A boolean array of all frames that satisfy the input condition
    """
    if "<=" in cond:
        tcond = cond.split("<=")
        ntmp = fitsdict[tcond[0]] <= float(tcond[1])
    elif ">=" in cond:
        tcond = cond.split(">=")
        ntmp = fitsdict[tcond[0]] >= float(tcond[1])
    elif "!=" in cond:
        tcond = cond.split("!=")
        ntmp = fitsdict[tcond[0]] != tcond[1]
    elif "<" in cond:
        tcond = cond.split("<")
        ntmp = fitsdict[tcond[0]] < float(tcond[1])
    elif ">" in cond:
        tcond = cond.split(">")
        ntmp = fitsdict[tcond[0]] > float(tcond[1])
    elif "=" in cond:
        tcond = cond.split("=")
        ntmp = (fitsdict[tcond[0]] == tcond[1])
    else:
        ntmp = None
    return ntmp


def sort_write(fitsdict, filesort, space=3):
    """
    Write out an xml and ascii file that contains the details of the file sorting.
    By default, the filename is printed first, followed by the filetype.
    After these, all parameters listed in the 'keyword' item in the
    settings file will be printed

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    filesort : dict
      Details of the sorted files
    space : int
      Keyword to set how many blank spaces to place between keywords
    """
    msgs.info("Preparing to write out the data sorting details")
    nfiles = fitsdict['filename'].size
    # Specify which keywords to print after 'filename' and 'filetype'
    prord = ['filename', 'frametype', 'target', 'exptime', 'naxis0', 'naxis1', 'filter1', 'filter2']
    prdtp = ["char",     "char",      "char",   "double",  "int",    "int",    "char",     "char"]
    # Now insert the remaining keywords:
    fkey = settings.spect['keyword'].keys()
    for i in fkey:
        if i not in prord:
            prord.append(i)
            # Append the type of value this keyword holds
            typv = type(fitsdict[i][0])
            if typv is int or typv is np.int_:
                prdtp.append("int")
            elif isinstance(fitsdict[i][0], basestring) or typv is np.string_:
                prdtp.append("char")
            elif typv is float or typv is np.float_:
                prdtp.append("double")
            else:
                msgs.bug("I didn't expect useful headers to contain type {0:s}".format(typv).replace('<type ', '').replace('>', ''))

    '''
    # Open a VOTable for writing
    votable = VOTableFile()
    resource = Resource()
    votable.resources.append(resource)
    table = Table(votable)
    resource.tables.append(table)
    # Define VOTable fields
    tabarr=[]
    # Insert the filename and filetype first
    for i in range(len(prord)):
        tabarr.append(Field(votable, name=prord[i], datatype=prdtp[i], arraysize="*"))
    table.fields.extend(tabarr)
    table.create_arrays(nfiles)
    filtyp = filesort.keys()
    for i in range(nfiles):
        values = ()
        for pr in prord:
            if pr == 'frametype':
                addval = ""
                for ft in filtyp:
                    if i in filesort[ft]:
                        if len(addval) != 0: addval += ","
                        addval += ft
                addval = (addval,)
            else: addval = (fitsdict[pr][i],)
            values = values + addval
        table.array[i] = values
    #osspl = sortname.split('.')
    #if len(osspl) > 1:
    #    fname = sortname
    #else:
    fname = settings.argflag['output']['sorted']+'.xml'
    votable.to_xml(fname)
    msgs.info("Successfully written sorted data information file:"+msgs.newline() +
              "{0:s}".format(fname))
    '''

    # ASCII file
    asciiord = ['filename', 'date', 'frametype', 'frameno', 'target', 'exptime', 'binning',
        'dichroic', 'dispname', 'dispangle', 'decker']
    # Generate the columns except frametype
    ascii_tbl = tTable()
    badclms = []
    for pr in asciiord:
        if pr != 'frametype':
            try:  # No longer require that all of these be present
                ascii_tbl[pr] = fitsdict[pr]
            except KeyError:
                badclms.append(pr)
    # Remove
    for pr in badclms:
        asciiord.pop(asciiord.index(pr))
    # Now frame type
    ftypes = []
    filtyp = filesort.keys()
    for i in range(nfiles):
        addval = ""
        for ft in filtyp:
            if i in filesort[ft]:
                if len(addval) != 0: addval += ","
                addval += ft
        ftypes.append(addval)
    ascii_tbl['frametype'] = ftypes
    # Write
    ascii_name = settings.argflag['output']['sorted']+'.lst'
    ascii_tbl[asciiord].write(ascii_name, format='ascii.fixed_width')
    return ascii_tbl


def match_science(fitsdict, filesort):
    """
    For a given set of identified data, match calibration frames to science frames

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    filesort : dict
      Details of the sorted files
    """

    msgs.info("Matching calibrations to Science frames")
    ftag = ['standard', 'bias', 'dark', 'pixelflat', 'pinhole', 'trace', 'arc']
    setup_ftag = dict(standard=0, bias=0, dark=0, pixelflat=0, pinhole=0, trace=0, arc=1)
    nfiles = fitsdict['filename'].size
    iSCI = filesort['science']
    iSTD = filesort['standard']
    iBIA = filesort['bias']
    iDRK = filesort['dark']
    iPFL = filesort['pixelflat']
    iBFL = filesort['pinhole']
    iTRC = filesort['trace']
    iARC = filesort['arc']
    filesort['failures'] = []
    iARR = [iSTD, iBIA, iDRK, iPFL, iBFL, iTRC, iARC]
    nSCI = iSCI.size
    i = 0
    while i < nSCI:
        msgs.info("Matching calibrations to {:s}: {:s}".format(
                fitsdict['target'][iSCI[i]], fitsdict['filename'][iSCI[i]]))
        settings.spect['science']['index'].append(np.array([iSCI[i]]))
        # Find nearby calibration frames
        for ft in range(len(ftag)):
            # Some checks first to make sure we need to find matching frames
            if ftag[ft] == 'dark' and settings.argflag['bias']['useframe'] != 'dark':
                msgs.info("  Dark frames not required")
                continue
            if ftag[ft] == 'bias' and settings.argflag['bias']['useframe'] != 'bias' and not settings.argflag['reduce']['badpix']:
                msgs.info("  Bias frames not required")
                continue
            # Now go ahead and match the frames
            n = np.arange(nfiles)
            if 'match' not in settings.spect[ftag[ft]].keys():
                debugger.set_trace()
            chkk = settings.spect[ftag[ft]]['match'].keys()
            for ch in chkk:
                tmtch = settings.spect[ftag[ft]]['match'][ch]
                if tmtch == "''":
                    w = np.where(fitsdict[ch] == fitsdict[ch][iSCI[i]])[0]
                elif tmtch[0] == '=':
                    mtch = np.float64(fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                    w = np.where((fitsdict[ch]).astype(np.float64) == mtch)[0]
                elif tmtch[0] == '<':
                    if tmtch[1] == '=':
                        mtch = np.float64(fitsdict[ch][iSCI[i]]) + np.float64(tmtch[2:])
                        w = np.where((fitsdict[ch]).astype(np.float64) <= mtch)[0]
                    else:
                        mtch = np.float64(fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                        w = np.where((fitsdict[ch]).astype(np.float64) < mtch)[0]
                elif tmtch[0] == '>':
                    if tmtch[1] == '=':
                        mtch = np.float64(fitsdict[ch][iSCI[i]]) + np.float64(tmtch[2:])
                        w = np.where((fitsdict[ch]).astype(np.float64) >= mtch)[0]
                    else:
                        mtch = np.float64(fitsdict[ch][iSCI[i]]) + np.float64(tmtch[1:])
                        w = np.where((fitsdict[ch]).astype(np.float64) > mtch)[0]
                elif tmtch[0] == '|':
                    if tmtch[1] == '=':
                        mtch = np.float64(tmtch[2:])
                        w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][iSCI[i]])) == mtch)[0]
                    elif tmtch[1] == '<':
                        if tmtch[2] == '=':
                            mtch = np.float64(tmtch[3:])
                            w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][iSCI[i]])) <= mtch)[0]
                        else:
                            mtch = np.float64(tmtch[2:])
                            w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][iSCI[i]])) < mtch)[0]
                    elif tmtch[1] == '>':
                        if tmtch[2] == '=':
                            mtch = np.float64(tmtch[3:])
                            w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][iSCI[i]])) >= mtch)[0]
                        else:
                            mtch = np.float64(tmtch[2:])
                            w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][iSCI[i]])) > mtch)[0]
                elif tmtch[0:2] == '%,':  # Splitting a header keyword
                    splcom = tmtch.split(',')
                    try:
                        spltxt, argtxt, valtxt = splcom[1], np.int(splcom[2]), splcom[3]
                        tspl = []
                        for sp in fitsdict[ch]:
                            tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
                            tmpspl = re.split(tmpspl, sp)
                            if len(tmpspl) < argtxt+1:
                                tspl.append("-9999999")
                            else:
                                tspl.append(tmpspl[argtxt])
                        tspl = np.array(tspl)
#                        debugger.set_trace()
                        tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
                        tmpspl = re.split(tmpspl, fitsdict[ch][iSCI[i]])
                        if len(tmpspl) < argtxt + 1:
                            continue
                        else:
                            scispl = tmpspl[argtxt]
                        if valtxt == "''":
                            w = np.where(tspl == scispl)[0]
                        elif valtxt[0] == '=':
                            mtch = np.float64(scispl) + np.float64(valtxt[1:])
                            w = np.where(tspl.astype(np.float64) == mtch)[0]
                        elif valtxt[0] == '<':
                            if valtxt[1] == '=':
                                mtch = np.float64(scispl) + np.float64(valtxt[2:])
                                w = np.where(tspl.astype(np.float64) <= mtch)[0]
                            else:
                                mtch = np.float64(scispl) + np.float64(valtxt[1:])
                                w = np.where(tspl.astype(np.float64) < mtch)[0]
                        elif valtxt[0] == '>':
                            if valtxt[1] == '=':
                                mtch = np.float64(scispl) + np.float64(valtxt[2:])
                                w = np.where(tspl.astype(np.float64) >= mtch)[0]
                            else:
                                mtch = np.float64(scispl) + np.float64(valtxt[1:])
                                w = np.where(tspl.astype(np.float64) > mtch)[0]
                    except:
                        debugger.set_trace()
                        continue
                else:
                    msgs.bug("Matching criteria {0:s} is not supported".format(tmtch))
                n = np.intersect1d(n, w)  # n corresponds to all frames with matching instrument setup to science frames
            # Find the time difference between the calibrations and science frames
            if settings.spect['fits']['calwin'] > 0.0:
                tdiff = np.abs(fitsdict['time'][n].astype(np.float64)-np.float64(fitsdict['time'][iSCI[i]]))
                w = np.where(tdiff <= settings.spect['fits']['calwin'])[0]
                n = n[w] # n corresponds to all frames within a set time difference of the science target frame
            # Now find which of the remaining n are the appropriate calibration frames
            n = np.intersect1d(n, iARR[ft])
            # How many frames are required
            if settings.argflag['run']['setup']:
                numfr = setup_ftag[ftag[ft]]
            else:
                numfr = settings.spect[ftag[ft]]['number']
            if settings.argflag['output']['verbosity'] == 2:
                if numfr == 1: areis = "is"
                else: areis = "are"
                if np.size(n) == 1:
                    msgs.info("  Found {0:d} {1:s} frame for {2:s} ({3:d} {4:s} required)".format(
                        np.size(n), ftag[ft], fitsdict['target'][iSCI[i]], numfr, areis))
                else:
                    msgs.info("  Found {0:d} {1:s} frames for {2:s} ({3:d} {4:s} required)".format(np.size(n), ftag[ft], fitsdict['target'][iSCI[i]], numfr, areis))
            # Have we identified enough of these calibration frames to continue?
            if np.size(n) < np.abs(numfr):
                msgs.warn("  Only {0:d}/{1:d} {2:s} frames for {3:s}".format(np.size(n), numfr, ftag[ft],
                                                                             fitsdict['target'][iSCI[i]]))
                # Errors for insufficient BIAS frames
                if settings.argflag['bias']['useframe'].lower() == ftag[ft]:
                    msgs.warn("Expecting to use bias frames for bias subtraction. But insufficient frames found.")
                    msgs.warn("Either include more frames or modify bias method")
                    msgs.warn("  e.g.:   bias useframe overscan")
                    msgs.error("Unable to continue")
                # Errors for insufficient PIXELFLAT frames
                if ftag[ft] == 'pixelflat' and settings.argflag['reduce']['flatfield']['perform'] and (
                    settings.argflag['reduce']['flatfield']['useframe'] == 'pixelflat'):
                    msgs.warn("Either include more frames or reduce the required amount with:")
                    msgs.warn("     pixelflat number XX")
                    msgs.warn("in the spect read/end block")
                    msgs.warn("Or specify a pixelflat file with --  reduce flatfield useframe file_name")
                    msgs.error("Unable to continue")
                # Errors for insufficient PINHOLE frames
                if ftag[ft] == 'pinhole':
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient TRACE frames
                if ftag[ft] == 'trace' and settings.argflag['reduce']['flatfield']['perform']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient standard frames
                if ftag[ft] == 'standard' and settings.argflag['reduce']['calibrate']['flux']:
                    msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
                # Errors for insufficient ARC frames
                if ftag[ft] == 'arc' and settings.argflag['reduce']['calibrate']:
                    if settings.argflag['run']['setup']:
                        msgs.warn("No arc frame for {0:s}. Removing it from list of science frames".format(fitsdict['target'][iSCI[i]]))
                        msgs.warn("Add an arc and rerun one if you wish to reduce this with PYPIT!!")
                        # Remove
                        #tmp = list(filesort['science'])
                        #tmp.pop(tmp.index(iSCI[i]))
                        #filesort['science'] = np.array(tmp)
                        filesort['failures'].append(iSCI[i])
                        settings.spect['science']['index'].pop(-1)
                    else:
                        msgs.error("Unable to continue without more {0:s} frames".format(ftag[ft]))
            else:
                # Select the closest calibration frames to the science frame
                tdiff = np.abs(fitsdict['time'][n].astype(np.float64)-np.float64(fitsdict['time'][iSCI[i]]))
                wa = np.argsort(tdiff)
                if (settings.argflag['run']['setup']) or (numfr < 0):
                    settings.spect[ftag[ft]]['index'].append(n[wa])
                else:
                    settings.spect[ftag[ft]]['index'].append(n[wa[:numfr]])
        i += 1
    msgs.info("Science frames successfully matched to calibration frames")
    return


def match_frames(frames, criteria, frametype='<None>', satlevel=None):
    """
    identify frames with a similar appearance (i.e. one frame appears to be a scaled version of another).
    """

    prob = arutils.erf(criteria/np.sqrt(2.0))[0]
    frsh0, frsh1, frsh2 = frames.shape
    msgs.info("Matching {:d} {:s} frames with confidence interval {:5.3%}".format(frsh2, frametype, prob))
    srtframes = [np.zeros((frsh0, frsh1, 1))]
    srtframes[0][:,:,0] = frames[:,:,0]
    tsrta = [frames[frsh0/2,:,0]]
    tsrtb = [frames[:,frsh1/2,0]]
    msgs.bug("Throughout this routine, you should probably search for the mean of the non-saturated pixels")
    tsrta[0] /= np.mean(tsrta[0])
    tsrtb[0] /= np.mean(tsrtb[0])
    for fr in range(1, frames.shape[2]):
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
            if np.size(wa[0]) == 0 or np.size(wb[0]) == 0:
                msgs.bug("I didn't expect to find a row of zeros in the middle of the chip!")
                sys.exit()
            if (testa >= prob) and (testa <= (2.0-prob)) and (testb >= prob) and (testb <= (2.0-prob)):
                fm = st
                break
        if fm is None:
            srtframes.append(np.zeros((frames.shape[0], frames.shape[1], 1)))
            srtframes[-1][:,:,0] = frames[:,:,fr]
            tsrta.append(tmata)
            tsrtb.append(tmatb)
        else:
            srtframes[fm] = np.append(srtframes[fm],np.zeros((frames.shape[0], frames.shape[1], 1)), axis=2)
            srtframes[fm][:,:,-1] = frames[:,:,fr]
    if len(srtframes) > 1:
        msgs.info("Found {0:d} different sets of {1:s} frames".format(len(srtframes), frametype))
    else:
        msgs.info("Found {0:d} set of {1:s} frames".format(len(srtframes), frametype))
    if frames.shape[2] > 1:
        del tsrta, tsrtb, tmata, tmatb, testa, testb
    return srtframes


def make_dirs(fitsdict, filesort):
    """
    For a given set of identified data, match calibration frames to science frames

    Parameters
    ----------
    fitsdict : dict
      Contains relevant information from fits header files
    filesort : dict
      Details of the sorted files

    Returns
    -------
    sci_targs : str array
      Names of the science targets
    """

    # First, get the current working directory
    currDIR = os.getcwd()
    msgs.info("Creating Science directory")
    newdir = "{0:s}/{1:s}".format(currDIR, settings.argflag['run']['directory']['science'])
    if os.path.exists(newdir):
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
        if not settings.argflag['output']['overwrite']:
            rmdir = ''
            while os.path.exists(newdir):
                while rmdir != 'n' and rmdir != 'y' and rmdir != 'r':
                    rmdir = raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, [r]ename) - ")
                if rmdir == 'n':
                    msgs.warn("Any previous calibration files may be overwritten")
                    break
                elif rmdir == 'r':
                    newdir = raw_input(msgs.input()+"Enter a new directory name: ")
                elif rmdir == 'y':
                    shutil.rmtree(newdir)
                    os.mkdir(newdir)
                    break
            if rmdir == 'r': os.mkdir(newdir)
    else: os.mkdir(newdir)
    # Create a directory for each object in the Science directory
    msgs.info("Creating Object directories")
    #Go through objects creating directory tree structure
    w = filesort['science']
    sci_targs = np.array(list(set(fitsdict['target'][w])))
    '''
    # Loop through targets and replace spaces with underscores
    nored = np.array([])
    # Create directories
    rmalways = False
    for i in range(sci_targs.size):
        sci_targs[i] = sci_targs[i].replace(' ', '_')
        newdir = "{0:s}/{1:s}/{2:s}".format(currDIR, settings.argflag['run']['directory']['science'], sci_targs[i])
        if os.path.exists(newdir):
            if settings.argflag['output']['overwrite'] or rmalways:
                pass
#				shutil.rmtree(newdir)
#				os.mkdir(newdir)
            else:
                msgs.info("The following directory already exists:"+msgs.newline()+newdir)
                rmdir = ''
                while rmdir != 'n' and rmdir != 'y' and rmdir != 'a':
                    rmdir = raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, or [a]lways) - ")
                if rmdir == 'n':
                    msgs.info("Not reducing {0:s}".format(sci_targs[i]))
                    nored = np.append(i)
                else:
                    shutil.rmtree(newdir)
                    os.mkdir(newdir)
                    if rmdir == 'a': rmalways = True
        else: os.mkdir(newdir)
    # Remove the entries from sci_targs which will not be reduced
    nored = nored.astype(np.int)
    while nored.size > 0:
        sci_targs = np.delete(sci_targs, nored[0])
        nored = np.delete(nored, 0)
    '''
    # Create a directory where all of the master calibration frames are stored.
    msgs.info("Creating Master Calibrations directory")
    newdir = "{:s}/{:s}_{:s}".format(currDIR, settings.argflag['run']['directory']['master'],
                                     settings.argflag['run']['spectrograph'])
    if os.path.exists(newdir):
        if not settings.argflag['output']['overwrite']:
            msgs.info("The following directory already exists:"+msgs.newline()+newdir)
            rmdir = ''
            while rmdir != 'n' and rmdir != 'y':
                rmdir = raw_input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o) - ")
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
    newdir = "{0:s}/{1:s}".format(currDIR, settings.argflag['run']['directory']['qa'])
    if os.path.exists(newdir):
        if not settings.argflag['output']['overwrite']:
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


def calib_set(isetup_dict, sciexp, fitsdict):
    """
    Ordering is 'aa', 'ab', ...
    Parameters
    ----------
    isetup_dict
    cset

    Returns
    -------

    """
    cb_strs = []
    default = 'aa'
    for ll in ['a', 'b', 'c', 'd']:
        for lower in string.ascii_lowercase:
            cb_strs.append(ll+lower)
    # Build cbset from sciexp
    new_cbset = {}
    cbkeys = ['arcs', 'bias', 'trace', 'flat', 'sci', 'cent']
    for cbkey in cbkeys:
        nms = list(fitsdict['filename'][getattr(sciexp,'_idx_'+cbkey)])
        nms.sort()
        new_cbset[cbkey] = nms

    # Uninitialized?
    if default not in isetup_dict.keys():
        isetup_dict[default] = new_cbset
        return default

    # Not in default
    for cbkey in cbkeys:
        def_names = np.array(isetup_dict[default][cbkey])
        if np.array_equal(def_names, new_cbset[cbkey]):
            _ = new_cbset.pop(cbkey)
    # Science only or exactly the same?
    if len(new_cbset) == 0:
        return default
    elif len(new_cbset) == 1:
        assert list(new_cbset.keys()) == ['sci']
        if new_cbset['sci'][0] not in isetup_dict[default]['sci']:
            isetup_dict[default]['sci'] += new_cbset['sci']
        return default

    # New calib set?
    for cb_str in cb_strs[1:]:
        mtch = True
        if cb_str not in isetup_dict.keys():
            isetup_dict[cb_str] = new_cbset
            break
        for key in new_cbset:
            if key not in isetup_dict[cb_str]:
                mtch = False
                break
            if key in ['sci']:
                continue
            ar_names = np.array(isetup_dict[default][cbkey])
            mtch &= np.array_equal(ar_names, new_cbset[key])
        if mtch:
            # Add sci frames
            for sciframe in new_cbset['sci']:
                break
    # Return
    return cb_str


def det_setup(isetup_dict, ddict):
    """ Return detector setup on config dict or add to it if new

    Parameters
    ----------
    isetup_dict : dict
      Selected setup_dict
    ddict : dict
      detector dict

    Returns
    -------
    dkey : str
      Name of the detector string associated to the input ddict
      May be new or previously used

    """
    det_str = [str('{:02d}'.format(i+1)) for i in range(99)]
    # Init
    for dkey in det_str:
        mtch = True
        if dkey not in isetup_dict.keys():
            isetup_dict[dkey] = ddict
            break
        for key in isetup_dict[dkey].keys():
            mtch &= isetup_dict[dkey][key] == ddict[key]
        if mtch:
            break
    # Return
    return dkey


def instr_setup(sciexp, det, fitsdict, setup_dict, must_exist=False,
                skip_cset=False, config_name=None):
    """ Define instrument config
    Make calls to detector and calib set

    config: A, B, C
    detector setup: 01, 02, 03, ..
    calib set: aa, ab, ac, ..


    Parameters
    ----------
    sciexp : ScienceExposure
    det : int
      detector identifier
    fitsdict : dict
      contains header info
    setup_dict : dict
    skip_cset : bool, optional
      Skip calib_set;  only used when first generating instrument .setup file
    config_name : str, optional
      Can be used to choose the config value

    Returns
    -------
    setup : str
      Full setup ID, e.g. A_01_aa
    """
    dnum = settings.get_dnum(det)
    # Labels
    cfig_str = string.ascii_uppercase
    cstr = '--'
    # Arc
    idx = sciexp._idx_arcs
    disp_name = fitsdict["dispname"][idx[0]]
    disp_angle = fitsdict["dispangle"][idx[0]]
    # Common
    dichroic = fitsdict["dichroic"][idx[0]]
    decker = fitsdict["decker"][idx[0]]
    slitwid = fitsdict["slitwid"][idx[0]]
    slitlen = fitsdict["slitlen"][idx[0]]

    # Detector -- These may not be set properly from the header alone, e.g. LRIS
    binning = fitsdict["binning"][idx[0]]
    naxis0 = fitsdict["naxis0"][idx[0]]
    naxis1 = fitsdict["naxis1"][idx[0]]
    namp = settings.spect[dnum]["numamplifiers"]

    # Generate
    # Don't nest deeper than 1
    cdict = dict(disperser={'name': disp_name,
                            'angle': disp_angle},
                 dichroic=dichroic,
                 slit={'decker': decker,
                       'slitwid': slitwid,
                       'slitlen': slitlen},
                 )
    ddict = {'binning': binning,
              'det': det,
              'namp': namp}
              #'naxis0': naxis0,
              #'naxis1': naxis1}

    # Configuration
    setup = None
    if len(setup_dict) == 0: # Generate
        if config_name is None:
            setup = 'A'
        else:
            setup = config_name
        # Finish
        setup_dict[setup] = {}
        setup_dict[setup][cstr] = cdict
    else:
        for ckey in setup_dict.keys():
            mtch = True
            for key in setup_dict[ckey][cstr].keys():
                # Dict?
                if isinstance(setup_dict[ckey][cstr][key], dict):
                    for ikey in setup_dict[ckey][cstr][key].keys():
                        mtch &= setup_dict[ckey][cstr][key][ikey] == cdict[key][ikey]
                else:
                    mtch &= setup_dict[ckey][cstr][key] == cdict[key]
            if mtch:
                setup = ckey
                break
        # Augment setup_dict?
        if setup is None:
            if must_exist:
                msgs.error("This setup is not present in the setup_dict.  Something went wrong..")
            maxs = max(setup_dict.keys())
            setup = cfig_str[cfig_str.index(maxs)+1]
            setup_dict[setup] = {}
            setup_dict[setup][cstr] = cdict

    # Detector
    dkey = det_setup(setup_dict[setup], ddict)
    # Calib set
    if not skip_cset:
        calib_key = calib_set(setup_dict[setup], sciexp, fitsdict)
    else:
        calib_key = '--'

    return '{:s}_{:s}_{:s}'.format(setup, dkey, calib_key)


def get_setup_file(spectrograph=None):
    """ Passes back name of setup file
    Also checks for existing setup files

    Parameters
    ----------
    spectrograph : str, optional

    Returns
    -------
    setup_file : str
      Name for the setup file
    nexist : int
      Number of existing setup files (0 or 1)

    """
    import glob
    import datetime

    if spectrograph is None:
        spectrograph = settings.argflag['run']['spectrograph']
    setup_files = glob.glob('./{:s}*.setups'.format(spectrograph))
    nexist = len(setup_files)
    # Require 1 or 0
    if nexist == 1:
        return setup_files[0], nexist
    elif nexist == 0:
        setup_file = settings.argflag['run']['redname'].replace('.pypit', '.setups')
        if os.path.isfile(setup_file):
            nexist = 1
        #date = str(datetime.date.today().strftime('%Y-%b-%d'))
        #return '{:s}_{:s}.setups'.format(spectrograph,date), nexist
        return setup_file, nexist
    else:
        msgs.error("Found more than one .setup file in the working directory.  Limit to one.")


def compare_setup(s1, s2):
    """ Compare two setup dicts
    Parameters
    ----------
    s1
    s2

    Returns
    -------

    """
    pass
    #for key in s1.keys():
    #    for key2 in s1[key]

def load_setup(**kwargs):
    """ Load setup from the disk

    Returns
    -------
    setup_dict : dict
    setup_file : str

    """
    import yaml
    setup_file, nexist = get_setup_file(**kwargs)
    if nexist == 0:
        debugger.set_trace()
        msgs.error("No existing setup file.  Generate one first (e.g. pypit_setup)!")
    # YAML
    with open(setup_file, 'r') as infile:
        setup_dict = yaml.load(infile)
    # Return
    return setup_dict, setup_file


def write_calib(setup_dict):
    """ Output setup_dict with full calibrations to hard drive

    Parameters
    ----------
    setup_dict : dict
      setup dict
    """
    calib_file = settings.argflag['run']['redname'].replace('.pypit', '.calib')
    # Write
    ydict = arutils.yamlify(setup_dict)
    with open(calib_file, 'w') as yamlf:
        yamlf.write(yaml.dump(ydict))

def write_setup(setup_dict, use_json=False):
    """ Output setup_dict to hard drive

    Parameters
    ----------
    setup_dict : dict
      setup dict
    use_json : bool, optional
      Output with JSON instead of YAML (not recommended)

    Returns
    -------

    """
    import json, io
    # Write
    setup_file, nexist = get_setup_file()
    if nexist == 1:
        msgs.warn("Over-writing existing .setups file")
    if use_json:
        gddict = ltu.jsonify(setup_dict)
        with io.open(setup_file, 'w', encoding='utf-8') as f:
            f.write(unicode(json.dumps(gddict, sort_keys=True, indent=4,
                                       separators=(',', ': '))))
    else: # YAML
        ydict = arutils.yamlify(setup_dict)
        with open(setup_file, 'w') as yamlf:
            yamlf.write(yaml.dump(ydict))

def load_sorted(sorted_file):
    """ Load a .sorted file (mainly to generate a .pypit file)
    Parameters
    ----------
    sorted_file : str

    Returns
    -------
    all_setups : list
     list of all setups eg. ['A','B']
    all_setuplines : list
     list of lists of all setup lines which describe the setup
    all_setupfiles : list
     list of lists of all setup files including the header
    """
    all_setups, all_setuplines, all_setupfiles = [], [], []
    with open(sorted_file,'r') as ff:
        # Should begin with ####
        fline = ff.readline()
        if fline[0:4] != '####':
            msgs.error('Bad .sorted fomatting')
        # Loop on setups
        while fline[0:5] != '##end':
            # Setup lines
            setuplines = []
            while fline[0:2] != '#-':
                fline = ff.readline()
                # Setup name
                if 'Setup' in fline:
                    all_setups.append(fline[6:].strip())
                #
                setuplines.append(fline)
            all_setuplines.append(setuplines[:-1])
            # Data files
            datafiles = []
            while fline[0:2] != '##':
                fline = ff.readline()
                datafiles.append(fline)
            all_setupfiles.append(datafiles[:-1])
    # Return
    return all_setups, all_setuplines, all_setupfiles


def write_sorted(srt_tbl, group_dict, setup_dict):
    """ Write the .sorted file
    Parameters
    ----------
    group_dict
    setup_dict

    Returns
    -------

    """
    # Output file
    group_file = settings.argflag['run']['redname'].replace('.pypit', '.sorted')
    ff = open(group_file, 'w')
    # Keys
    setups = list(group_dict.keys())
    setups.sort()
    ftypes = list(group_dict[setups[0]].keys())
    ftypes.sort()
    # Loop on Setup
    asciiord = np.array(['filename', 'date', 'frameno', 'frametype', 'target', 'exptime', 'dispname', 'decker'])
    for setup in setups:
        ff.write('##########################################################\n')
        in_setup = []
        ff.write('Setup {:s}\n'.format(setup))
        ydict = arutils.yamlify(setup_dict[setup])
        ff.write(yaml.dump(ydict))
        ff.write('#---------------------------------------------------------\n')
        # ID files
        for key in ftypes:
            # Sort
            gfiles = group_dict[setup][key]
            gfiles.sort()
            for ifile in gfiles:
                mt = np.where(srt_tbl['filename'] == ifile)[0]
                if (len(mt) > 0) and (mt not in in_setup):
                    in_setup.append(mt[0])
        # Write overlapping keys
        gdkeys = np.in1d(asciiord, np.array(srt_tbl.keys()))
        subtbl = srt_tbl[asciiord[gdkeys].tolist()][np.array(in_setup)]
        subtbl.write(ff, format='ascii.fixed_width')
    ff.write('##end\n')
    ff.close()
