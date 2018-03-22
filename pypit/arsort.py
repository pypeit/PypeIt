""" Routines for sorting data to be reduced by PYPIT"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import re
import sys
import shutil
import string
import numpy as np
import yaml

from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
from astropy.table import Table as tTable, Column
from astropy import units as u

from linetools import utils as ltu

from pypit import armsgs
from pypit import arparse as settings
from pypit import arutils
from pypit.arflux import find_standard_file
from pypit import armeta
from pypit import ardebug as debugger

try:
    basestring
except NameError:
    basestring = str

try: input = raw_input
except NameError: pass

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
    filesort : dict
      A dictionary of filetypes
      Each key is a file type and contains an array of the file indices that qualify
    """
    msgs.bug("There appears to be a bug with the assignment of arc frames when only one science frame is supplied")
    msgs.info("Sorting files")
    numfiles = fitsdict['filename'].size
    # Set the filetype dictionary
    filesort = {}
    for ftype in armeta.allowed_file_types():
        filesort[ftype] = np.array([], dtype=np.int)
    # Set all filetypes by hand?
    if len(settings.ftdict) > 0:
        for ifile,ftypes in settings.ftdict.items():
            idx = np.where(fitsdict['filename'] == ifile)[0]
            sptypes = ftypes.split(',')
            for iftype in sptypes:
                filesort[iftype] = np.concatenate([filesort[iftype], idx])
        # Sort
        for key in filesort.keys():
            filesort[key].sort()
        return filesort
    #  Prepare to type
    fkeys = np.array(list(filesort.keys()))
    # Create an array where 1 means it is a certain type of frame and 0 means it isn't.
    filarr = np.zeros((fkeys.size, numfiles), dtype=np.int)
    setarr = np.zeros((fkeys.size, numfiles), dtype=np.int)

    # Identify the frames:
    # Loop on file type
    for i, fkey in enumerate(fkeys):
        if fkey == 'unknown':
            continue
        # Self identification (typically from Header; not recommended)
        if settings.argflag['run']['useIDname']:
            w = np.where(fitsdict['idname'] == settings.spect[fkey]['idname'])[0]
        else:
            w = np.arange(numfiles)
        n = np.arange(numfiles)
        n = np.intersect1d(n, w)
        # Perform additional checks in order to make sure this identification is true
        if 'check' in settings.spect[fkey].keys():
            n = chk_all_conditions(n, fkey, fitsdict)

        # Assign these images to the filetype
        filarr[i, :][n] = 1
        # Check if these files can also be another type
        #  e.g. some frames are used for pixelflat and slit tracing
        if settings.spect[fkey]['canbe'] is not None:
            for cb in settings.spect[fkey]['canbe']:
                # Assign these filetypes
                fa = np.where(fkeys == cb)[0]
                if np.size(fa) == 1:
                    filarr[fa[0], :][n] = 1
                else:
                    msgs.error("Unknown type for argument 'canbe': {0:s}".format(cb))

    # Identify the standard stars
    # Find the nearest standard star to each science frame
    wscistd = np.where(filarr[np.where(fkeys == 'standard')[0], :].flatten() == 1)[0]
    for i in range(wscistd.size):
        radec = (fitsdict['ra'][wscistd[i]], fitsdict['dec'][wscistd[i]])
        if fitsdict['ra'][wscistd[i]] == 'None':
            msgs.warn("No RA and DEC information for file:" + msgs.newline() + fitsdict['filename'][wscistd[i]])
            msgs.warn("The above file could be a twilight flat frame that was" + msgs.newline() +
                      "missed by the automatic identification.")
            filarr[np.where(fkeys == 'standard')[0], wscistd[i]] = 0
            continue
        # If an object exists within 20 arcmins of a listed standard, then it is probably a standard star
        foundstd = find_standard_file(radec, toler=20.*u.arcmin, check=True)
        if foundstd:
            filarr[np.where(fkeys == 'science')[0], wscistd[i]] = 0
        else:
            filarr[np.where(fkeys == 'standard')[0], wscistd[i]] = 0

    # Make any forced changes
    msgs.info("Making forced file identification changes")
    skeys = settings.spect['set'].keys()
    for sk in skeys:
        for j in settings.spect['set'][sk]:
            w = np.where(fitsdict['filename']==j)[0]
            filarr[:,w]=0
            setarr[np.where(fkeys==sk)[0],w]=1
    filarr = filarr + setarr

    # Check that all files have an identification
    badfiles = np.where(np.sum(filarr, axis=0) == 0)[0]
    if np.size(badfiles) != 0:
        msgs.info("Couldn't identify the following files:")
        for i in range(np.size(badfiles)):
            msgs.info(fitsdict['filename'][badfiles[i]])
        if flag_unknown:
            filarr[np.where(fkeys == 'unknown')[0],badfiles] = 1
        else:
            msgs.error("Check these files and your settings.{0:s} file before continuing".format(settings.argflag['run']['spectrograph']))

    # Now identify the dark frames
    wdark = np.where((filarr[np.where(fkeys == 'bias')[0], :] == 1).flatten() &
                     (fitsdict['exptime'].astype(np.float64) > settings.spect['mosaic']['minexp']))[0]
    filesort['dark'] = wdark

    # Store the frames in the filesort array
    for i,fkey in enumerate(fkeys):
        filesort[fkey] = np.where(filarr[i,:] == 1)[0]
    # Finally check there are no duplicates (the arrays will automatically sort with np.unique)
    msgs.info("Finalising frame sorting, and removing duplicates")
    for key in filesort.keys():
        filesort[key] = np.unique(filesort[key])
        if np.size(filesort[key]) == 1:
            msgs.info("Found {0:d} {1:s} frame".format(np.size(filesort[key]),key))
        else:
            msgs.info("Found {0:d} {1:s} frames".format(np.size(filesort[key]),key))
    # Return filesort!
    msgs.info("Sorting completed successfully")
    return filesort


def chk_all_conditions(n, fkey, fitsdict):
    """ Loop on the conditions for this given file type
    Parameters
    ----------
    n : ndarray
      Indices of images satisfying the filetype thus far
    fkey : str
      File type
    fitsdict : dict

    Returns
    -------
    n : ndarray
      Indices of images also satisfying the check conditions
    """
    chkk = settings.spect[fkey]['check'].keys()
    for ch in chkk:
        if ch[0:9] == 'condition':
            # Deal with a conditional argument
            conds = re.split("(\||\&)", settings.spect[fkey]['check'][ch])
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
                w = np.where(np.char.strip(fitsdict[ch]) == settings.spect[fkey]['check'][ch])[0]
            else:
                w = np.where(fitsdict[ch] == settings.spect[fkey]['check'][ch])[0]
        n = np.intersect1d(n, w)
    # Return
    return n


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
        if 'int' in fitsdict[tcond[0]].dtype.name:
            ntmp = fitsdict[tcond[0]] != int(tcond[1])
        elif 'float' in fitsdict[tcond[0]].dtype.name:
            ntmp = fitsdict[tcond[0]] != float(tcond[1])
        else:
            ntmp = fitsdict[tcond[0]] != tcond[1]
    elif "<" in cond:
        tcond = cond.split("<")
        ntmp = fitsdict[tcond[0]] < float(tcond[1])
    elif ">" in cond:
        tcond = cond.split(">")
        ntmp = fitsdict[tcond[0]] > float(tcond[1])
    elif "=" in cond:
        tcond = cond.split("=")
        if 'int' in fitsdict[tcond[0]].dtype.name:
            ntmp = fitsdict[tcond[0]] == int(tcond[1])
        elif 'float' in fitsdict[tcond[0]].dtype.name:
            ntmp = fitsdict[tcond[0]] == float(tcond[1])
        else:
            ntmp = fitsdict[tcond[0]] == tcond[1]
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
                msgs.bug("I didn't expect useful headers to contain type {!s:s}".format(typv).replace('<type ', '').replace('>', ''))

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
    if settings.argflag['run']['setup']:
        ascii_name = settings.argflag['run']['redname'].replace('.pypit', '.lst')
    else:
        ascii_name = settings.argflag['output']['sorted']+'.lst'
    ascii_tbl[asciiord].write(ascii_name, format='ascii.fixed_width')
    return ascii_tbl

def match_logic(ch, tmtch, fitsdict, idx):
    """ Perform logic on matching with fitsdict
    Parameters
    ----------
    ch : str
      Header card alias, eg. exptime
    tmtch : str
      Defines the logic
      any
      ''
      >, <, >=, <=, =, !=


    Returns
    -------
    w : ndarray, int
      indices that match the criterion
      None is returned if there is nothing to match
    """
    if tmtch == "any":
        w = np.arange(len(fitsdict['filename'])).astype(int)
    elif tmtch == "''":
        w = np.where(fitsdict[ch] == fitsdict[ch][idx])[0]
    elif tmtch[0] == '=':
        mtch = np.float64(fitsdict[ch][idx]) + np.float64(tmtch[1:])
        w = np.where((fitsdict[ch]).astype(np.float64) == mtch)[0]
    elif tmtch[0] == '<':
        if tmtch[1] == '=':
            mtch = np.float64(fitsdict[ch][idx]) + np.float64(tmtch[2:])
            w = np.where((fitsdict[ch]).astype(np.float64) <= mtch)[0]
        else:
            mtch = np.float64(fitsdict[ch][idx]) + np.float64(tmtch[1:])
            w = np.where((fitsdict[ch]).astype(np.float64) < mtch)[0]
    elif tmtch[0] == '>':
        if tmtch[1] == '=':
            mtch = np.float64(fitsdict[ch][idx]) + np.float64(tmtch[2:])
            w = np.where((fitsdict[ch]).astype(np.float64) >= mtch)[0]
        else:
            mtch = np.float64(fitsdict[ch][idx]) + np.float64(tmtch[1:])
            w = np.where((fitsdict[ch]).astype(np.float64) > mtch)[0]
    elif tmtch[0] == '|':
        if tmtch[1] == '=':
            mtch = np.float64(tmtch[2:])
            w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][idx])) == mtch)[0]
        elif tmtch[1] == '<':
            if tmtch[2] == '=':
                mtch = np.float64(tmtch[3:])
                w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][idx])) <= mtch)[0]
            else:
                mtch = np.float64(tmtch[2:])
                w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][idx])) < mtch)[0]
        elif tmtch[1] == '>':
            if tmtch[2] == '=':
                mtch = np.float64(tmtch[3:])
                w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][idx])) >= mtch)[0]
            else:
                mtch = np.float64(tmtch[2:])
                w = np.where(np.abs((fitsdict[ch]).astype(np.float64)-np.float64(fitsdict[ch][idx])) > mtch)[0]
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
            tmpspl = re.split(tmpspl, fitsdict[ch][idx])
            if len(tmpspl) < argtxt + 1:
                return None
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
            return None
    # Return
    return w


def match_science(fitsdict, filesort):
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
    cal_indx : list
      A list of dict's, one per science exposure, that contains the
      indices of the matched calibration files (and the science frame too)
      This is intended to replace settings.spect[ftag]['index']
    """

    msgs.info("Matching calibrations to Science frames")
    # Init
    ftag = armeta.allowed_file_types()
    assert ftag[0] == 'arc'
    for item in ['science', 'unknown']: # Remove undesired
        ftag.remove(item)
    setup_ftag = {}
    for item in ftag:
        setup_ftag[item] = 0
    setup_ftag['arc'] = 1
    # More setup
    nfiles = fitsdict['filename'].size
    iSCI = filesort['science']
    iARR = [filesort[itag] for itag in ftag]
    filesort['failures'] = []
    nSCI = iSCI.size
    tmp = {}
    for iftag in ftag:
        tmp[iftag] = []
    cal_index = []
    for kk in range(nSCI):
        cal_index.append(tmp.copy())
    # Loop on science frames
    i = 0
    while i < nSCI:
        msgs.info("Matching calibrations to {:s}: {:s}".format(
                fitsdict['target'][iSCI[i]], fitsdict['filename'][iSCI[i]]))
        # Science index (trivial)
        settings.spect['science']['index'].append(np.array([iSCI[i]]))
        cal_index[i]['science'] = np.array([iSCI[i]])
        # Find matching (and nearby) calibration frames
        for ft in range(len(ftag)):
            # Some checks first to make sure we need to find matching frames
            if ftag[ft] == 'dark' and settings.argflag['bias']['useframe'] != 'dark':
                msgs.info("  Dark frames not required.  Not matching..")
                continue
            if ftag[ft] == 'bias' and settings.argflag['bias']['useframe'] != 'bias' and not settings.argflag['reduce']['badpix']:
                msgs.info("  Bias frames not required.  Not matching..")
                continue
            # How many matching frames are required?  This is instrument specific
            if settings.argflag['run']['setup']:
                numfr = setup_ftag[ftag[ft]]
            else:
                numfr = settings.spect[ftag[ft]]['number']
            if (numfr == 0) and (not settings.argflag['run']['setup']):
                settings.spect[ftag[ft]]['index'].append(np.array([]))
                cal_index[i][ftag[ft]] = np.array([])
                msgs.info("No {0:s} frames are required.  Not matching..".format(ftag[ft]))
                continue
            # Now go ahead and match the frames
            n = np.arange(nfiles)
            if 'match' not in settings.spect[ftag[ft]].keys() and (not settings.argflag['run']['setup']):
                msgs.error("Need match criteria for {0:s}!!".format(ftag[ft]))
            elif 'match' not in settings.spect[ftag[ft]].keys():
                msgs.info("No matching criteria for {0:s} frames with this instrument".format(ftag[ft]))
            else:
                chkk = settings.spect[ftag[ft]]['match'].keys()
                for ch in chkk:
                    tmtch = settings.spect[ftag[ft]]['match'][ch]
                    w = match_logic(ch, tmtch, fitsdict, iSCI[i])
                    if w is not None:
                        n = np.intersect1d(n, w)  # n corresponds to all frames with matching instrument setup to science frames
            # Find the time difference between the calibrations and science frames
            if settings.spect['fits']['calwin'] > 0.0:
                tdiff = np.abs(fitsdict['time'][n].astype(np.float64)-np.float64(fitsdict['time'][iSCI[i]]))
                w = np.where(tdiff <= settings.spect['fits']['calwin'])[0]
                n = n[w] # n corresponds to all frames within a set time difference of the science target frame
            # Now find which of the remaining n are the appropriate calibration frames
            n = np.intersect1d(n, iARR[ft])
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
                code = match_warnings(ftag[ft], n.size, numfr, fitsdict, iSCI[i], filesort, cal_index[i])
                if code == 'break':
                    break
            else:
                # Select the closest calibration frames to the science frame
                tdiff = np.abs(fitsdict['time'][n].astype(np.float64)-np.float64(fitsdict['time'][iSCI[i]]))
                wa = np.argsort(tdiff)
                if (settings.argflag['run']['setup']) or (numfr < 0):
                    settings.spect[ftag[ft]]['index'].append(n[wa].copy())
                    cal_index[i][ftag[ft]] = n[wa].copy()
                else:
                    settings.spect[ftag[ft]]['index'].append(n[wa[:numfr]].copy())
                    cal_index[i][ftag[ft]] = n[wa[:numfr]].copy()
        i += 1
    msgs.info("Science frames successfully matched to calibration frames")
    return cal_index


def match_warnings(ftag, nmatch, numfr, fitsdict, idx, filesort, cindex):
    """ Give warnings related to matching calibration files to science
    Returns
    -------
    code : str or None
      None = no further action required
    """
    msgs.warn("  Only {0:d}/{1:d} {2:s} frames for {3:s}".format(nmatch, numfr, ftag,
                                                                 fitsdict['target'][idx]))
    # Errors for insufficient BIAS frames
    if settings.argflag['bias']['useframe'].lower() == ftag:
        msgs.warn("Expecting to use bias frames for bias subtraction. But insufficient frames found.")
        msgs.warn("Either include more frames or modify bias method" + msgs.newline() +
                  "  e.g.:   bias useframe overscan")
        msgs.error("Unable to continue")
    # Errors for insufficient PIXELFLAT frames
    if ftag == 'pixelflat' and settings.argflag['reduce']['flatfield']['perform'] and (
                settings.argflag['reduce']['flatfield']['useframe'] == 'pixelflat'):
        if settings.argflag['reduce']['masters']['force']:
            msgs.warn("Fewer pixelflat frames than expected for {0:s}, but will use MasterFrames".format(fitsdict['target'][idx]))
        else:
            msgs.warn("Either include more frames or reduce the required amount with:" + msgs.newline() +
                      "pixelflat number XX" + msgs.newline() +
                      "in the spect read/end block")
            msgs.warn("Or specify a pixelflat file with --  reduce flatfield useframe file_name")
            msgs.error("Unable to continue")
    # Errors for insufficient PINHOLE frames
    if ftag == 'pinhole':
        msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    # Errors for insufficient TRACE frames
    if ftag == 'trace' and settings.argflag['reduce']['flatfield']['perform']:
        if settings.argflag['reduce']['masters']['force']:
            msgs.warn("Fewer traceflat frames than expected for {0:s}, but will use MasterFrames".format(fitsdict['target'][idx]))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    # Errors for insufficient standard frames
    if ftag == 'standard' and settings.argflag['reduce']['calibrate']['flux']:
        if settings.argflag['reduce']['masters']['force']:
            msgs.warn("No standard star frames for {0:s}, but will use MasterFrames".format(fitsdict['target'][idx]))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    # Errors for insufficient ARC frames
    if ftag == 'arc' and settings.argflag['reduce']['calibrate']:
        if settings.argflag['run']['setup']:
            msgs.warn("No arc frame for {0:s}. Removing it from list of science frames".format(fitsdict['target'][idx]))
            msgs.warn("Add an arc and rerun one if you wish to reduce this with PYPIT!!")
            # Remove
            #tmp = list(filesort['science'])
            #tmp.pop(tmp.index(idx))
            #filesort['science'] = np.array(tmp)
            filesort['failures'].append(idx)
            settings.spect['science']['index'].pop(-1)
            cindex['science'].pop(-1)
            return 'break'
        elif settings.argflag['reduce']['masters']['force']:
            msgs.warn("No arc frame for {0:s}, but will use MasterFrames".format(fitsdict['target'][idx]))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))

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
                    rmdir = input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, [r]ename) - ")
                if rmdir == 'n':
                    msgs.warn("Any previous calibration files may be overwritten")
                    break
                elif rmdir == 'r':
                    newdir = input(msgs.input()+"Enter a new directory name: ")
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
                    rmdir = input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o, or [a]lways) - ")
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
                rmdir = input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o) - ")
            if rmdir == 'n':
                msgs.warn("Any previous calibration files will be overwritten")
            else:
                shutil.rmtree(newdir)
                os.mkdir(newdir)
#		else:
#			shutil.rmtree(newdir)
#			os.mkdir(newdir)
    else: os.mkdir(newdir)
    # Create a directory where all of the QA is stored
    msgs.info("Creating QA directory")
    newdir = "{0:s}/{1:s}".format(currDIR, settings.argflag['run']['directory']['qa'])
    if os.path.exists(newdir):
        msgs.warn("Pre-existing QA plots will be overwritten")
        '''
        if not settings.argflag['output']['overwrite']:
            msgs.info("The following directory already exists:"+msgs.newline()+newdir)
            rmdir=''
            while rmdir != 'n' and rmdir != 'y':
                rmdir=input(msgs.input()+"Remove this directory and it's contents? ([y]es, [n]o) - ")
            if rmdir == 'n':
                msgs.warn("Any previously made plots will be overwritten")
            else:
                shutil.rmtree(newdir)
                os.mkdir(newdir)
        else:
            shutil.rmtree(newdir)
            os.mkdir(newdir)
            os.mkdir(newdir+'/PNGs')
        '''
        if not os.path.exists(newdir+'/PNGs'):
            os.mkdir(newdir+'/PNGs')
    else:
        os.mkdir(newdir)
        os.mkdir(newdir+'/PNGs')
    # Return the name of the science targets
    return sci_targs

