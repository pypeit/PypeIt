""" Routines for sorting data to be reduced by PYPIT"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

try:
    basestring
except NameError:
    basestring = str

'''
try:
    input = raw_input
except NameError:
    pass
'''

import os
import re
import sys
import shutil

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units

from pypit import msgs
#from pypit import arparse as settings
from pypit import arutils
from pypit.arflux import find_standard_file
from pypit import ardebug as debugger

ftype_list = [     # NOTE:  arc must be listed first!
    'arc',         # Exposure of one or more arc calibration lamps for wavelength calibration
    'bias',        # Exposure for assessing detector bias (usually 0s)
    'dark',        # Exposure for assessing detector dark current
    'pinhole',     # Exposure for tracing the orders or slits
    'pixelflat',   # Exposure for assessing pixel-to-pixel variations
    'science',   # Exposure on one or more science targets
    'standard',    # Exposure on a 'standard star' used for flux calibration
    'trace',       # Exposure for tracing slits or echelle orders (usually twilight sky or flat lamp)
    'unknown',     # Unknown..
]


def ftype_indices(fitstbl, ftype, sci_ID):
    """

    Parameters
    ----------
    fitstbl : Table
    ftype : str
      e.g. arc, trace, science
    sci_ID : int
      ID value of the science exposure
      Binary, i.e.  1, 2, 4, 8, 16..

    Returns
    -------
    idx : ndarray (int)
      Indices of the rows of the Table matching the inputs

    """
    idx = np.where(fitstbl[ftype] & (fitstbl['sci_ID'] & sci_ID > 0))[0]
    return idx


def type_data(fitstbl, settings_spect, settings_argflag, flag_unknown=False, ftdict=None):
    """ Generate a table of filetypes from the input fitsdict object

    Parameters
    ----------
    fitstbl : dict
      Contains relevant information from fits header files
    flag_unknown : bool, optional
      Instead of crashing out if there are unidentified files,
      set to 'unknown' and continue

    Returns
    -------
    filetypes : Table
      A Table of filetypes
      Each key is a file type and contains False/True for each datafile
      This is stacked onto the fitstbl
    """
    msgs.info("Typing files")
    numfiles = fitstbl['filename'].size
    # Set the filetype dictionary
    filetypes = Table()
    for ftype in ftype_list:
        filetypes[ftype] = np.zeros(numfiles, dtype=bool)

    # Set all filetypes by hand?  (Typically read from PYPIT file)
    #if len(settings.ftdict) > 0:
    if ftdict is not None:
        for ifile,ftypes in ftdict.items():
            idx = fitstbl['filename'] == ifile
            sptypes = ftypes.split(',')
            for iftype in sptypes:
                filetypes[iftype][idx] = True
        # Sort
        #for key in filetypesfilesort.keys():
        #    filetypesfilesort[key].sort()
        return filetypes

    # Identify the frames:
    # Loop on file type
    for i, ftype in enumerate(ftype_list):
        if ftype == 'unknown':
            continue
        # Self identification (typically from Header; not recommended)
        if settings_argflag['run']['useIDname']:
            idx = fitstbl['idname'] == settings_spect[ftype]['idname']
            filetypes[ftype][idx] = True
            #w = np.where(fitsdict['idname'] == settings.spect[ftype]['idname'])[0]
        else: # Set all to True!
            filetypes[ftype] = True

        # Perform additional checks in order to make sure this identification is true
        if 'check' in settings_spect[ftype].keys():
            gd_chk = chk_all_conditions(ftype, fitstbl, settings_spect)
            filetypes[ftype] &= gd_chk

    # Identify the standard stars
    # Find the nearest standard star to each science frame
    wscistds = np.where(filetypes['standard'])[0]
    for wscistd in wscistds:
        radec = (fitstbl['ra'][wscistd], fitstbl['dec'][wscistd])
        if fitstbl['ra'][wscistd] == 'None':
            msgs.warn("No RA and DEC information for file:" + msgs.newline() + fitstbl['filename'][wscistd])
            msgs.warn("The above file could be a twilight flat frame that was" + msgs.newline() +
                      "missed by the automatic identification.")
            filetypes['standard'][wscistd] = False
            continue
        # If an object exists within 20 arcmins of a listed standard, then it is probably a standard star
        foundstd = find_standard_file(radec, toler=20.*units.arcmin, check=True)
        if foundstd:
            filetypes['science'][wscistd] = False
        else:
            filetypes['standard'][wscistd] = False

    # Make any forced changes
    skeys = settings_spect['set'].keys()
    if len(skeys) > 0:
        msgs.info("Making forced file identification changes")
        msgs.warn("Note that the image will have *only* the specified type")
        for sk in skeys:
            for jj in settings_spect['set'][sk]:
                idx = np.where(fitstbl['filename']==jj)[0]
                # Zero out the others
                for ftype in ftype_list:
                    filetypes[ftype][idx] = False
                # And set
                filetypes[sk][idx] = True

    # Check that all files have an identification
    chklist = []
    for ftype in ftype_list:
        if ftype == 'unknown':
            continue
        chklist.append(filetypes[ftype].data)
    badfiles = ~np.any(chklist,axis=0)
    if np.any(badfiles):
        msgs.info("Couldn't identify the following files:")
        for ifile in fitstbl['filename'][badfiles]:
            msgs.info(ifile)
        if flag_unknown:
            filetypes['unknown'][badfiles] = True
        else:
            msgs.error("Check these files and your settings.{0:s} file before continuing".format(settings_argflag['run']['spectrograph']))

    # Now identify the dark frames
    darks = filetypes['bias'] & (fitstbl['exptime'].data.astype(np.float64) >
                                   settings_spect['mosaic']['minexp'])
    filetypes['dark'] = darks

    # Return filesort!
    msgs.info("Typing completed!")
    return filetypes

'''
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
    for ftype in ftype_list:
        filesort[ftype] = np.array([], dtype=np.int)
    # Set all filetypes by hand (typically read from PYPIT file)?
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
    ftypes = np.array(list(filesort.keys()))
    # Create an array where 1 means it is a certain type of frame and 0 means it isn't.
    filarr = np.zeros((ftypes.size, numfiles), dtype=np.int)
    setarr = np.zeros((ftypes.size, numfiles), dtype=np.int)

    # Identify the frames:
    # Loop on file type
    for i, iftype in enumerate(ftypes):
        if iftype == 'unknown':
            continue
        # Self identification (typically from Header; not recommended)
        if settings.argflag['run']['useIDname']:
            w = np.where(fitsdict['idname'] == settings.spect[iftype]['idname'])[0]
        else:
            w = np.arange(numfiles)
        n = np.arange(numfiles)
        n = np.intersect1d(n, w)
        # Perform additional checks in order to make sure this identification is true
        if 'check' in settings.spect[iftype].keys():
            n = chk_all_conditions(n, iftype, fitsdict)

        # Assign these images to the filetype
        filarr[i, :][n] = 1
        # Check if these files can also be another type
        #  e.g. some frames are used for pixelflat and slit tracing
        if settings.spect[iftype]['canbe'] is not None:
            for cb in settings.spect[iftype]['canbe']:
                # Assign these filetypes
                fa = np.where(ftypes == cb)[0]
                if np.size(fa) == 1:
                    filarr[fa[0], :][n] = 1
                else:
                    msgs.error("Unknown type for argument 'canbe': {0:s}".format(cb))

    # Identify the standard stars
    # Find the nearest standard star to each science frame
    wscistd = np.where(filarr[np.where(iftype == 'standard')[0], :].flatten() == 1)[0]
    for i in range(wscistd.size):
        radec = (fitsdict['ra'][wscistd[i]], fitsdict['dec'][wscistd[i]])
        if fitsdict['ra'][wscistd[i]] == 'None':
            msgs.warn("No RA and DEC information for file:" + msgs.newline() + fitsdict['filename'][wscistd[i]])
            msgs.warn("The above file could be a twilight flat frame that was" + msgs.newline() +
                      "missed by the automatic identification.")
            filarr[np.where(ftypes == 'standard')[0], wscistd[i]] = 0
            continue
        # If an object exists within 20 arcmins of a listed standard, then it is probably a standard star
        foundstd = find_standard_file(radec, toler=20.*units.arcmin, check=True)
        if foundstd:
            filarr[np.where(ftypes == 'science')[0], wscistd[i]] = 0
        else:
            filarr[np.where(ftypes == 'standard')[0], wscistd[i]] = 0

    # Make any forced changes
    msgs.info("Making forced file identification changes")
    skeys = settings.spect['set'].keys()
    for sk in skeys:
        for j in settings.spect['set'][sk]:
            w = np.where(fitsdict['filename']==j)[0]
            filarr[:,w]=0
            setarr[np.where(ftypes==sk)[0],w]=1
    filarr = filarr + setarr

    # Check that all files have an identification
    badfiles = np.where(np.sum(filarr, axis=0) == 0)[0]
    if np.size(badfiles) != 0:
        msgs.info("Couldn't identify the following files:")
        for i in range(np.size(badfiles)):
            msgs.info(fitsdict['filename'][badfiles[i]])
        if flag_unknown:
            filarr[np.where(ftypes == 'unknown')[0],badfiles] = 1
        else:
            debugger.set_trace()
            msgs.error("Check these files and your settings.{0:s} file before continuing".format(settings.argflag['run']['spectrograph']))

    # Now identify the dark frames
    wdark = np.where((filarr[np.where(ftypes == 'bias')[0], :] == 1).flatten() &
                     (fitsdict['exptime'].astype(np.float64) > settings.spect['mosaic']['minexp']))[0]
    filesort['dark'] = wdark

    # Store the frames in the filesort array
    for i,iftype in enumerate(ftypes):
        filesort[iftype] = np.where(filarr[i,:] == 1)[0]
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
'''


def chk_all_conditions(fkey, fitstbl, settings_spect):
    """ Loop on the conditions for this given file type
    Parameters
    ----------
    fkey : str
      File type
    fitstbl : Table

    Returns
    -------
    gd_chk : ndarray
      Indices of images satisfying the check conditions
    """
    gd_chk = np.ones(len(fitstbl), dtype=bool)
    # Loop on the items to check
    chkk = settings_spect[fkey]['check'].keys()
    for ch in chkk:
        if ch[0:9] == 'condition':
            # Deal with a conditional argument
            conds = re.split("(\||\&)", settings_spect[fkey]['check'][ch])
            ntmp = chk_condition(fitstbl, conds[0])
            # And more
            for cn in range((len(conds)-1)//2):
                if conds[2*cn+1] == "|":
                    ntmp = ntmp | chk_condition(fitstbl, conds[2*cn+2])
                elif conds[2*cn+1] == "&":
                    ntmp = ntmp & chk_condition(fitstbl, conds[2*cn+2])
            gd_chk = gd_chk & ntmp
        else:
            if fitstbl[ch].dtype.char in ['S','U']:  # Numpy string array
                # Strip numpy string array of all whitespace
                gd_chk = gd_chk & (np.char.strip(fitstbl[ch]) == settings_spect[fkey]['check'][ch])
            else:
                gd_chk = gd_chk & (fitstbl[ch] == settings_spect[fkey]['check'][ch])
    # Return
    return gd_chk

'''
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
'''


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


def write_lst(fitstbl, settings_spect, settings_argflag):
    """
    Write out an ascii file that contains the details of the file sorting.
    By default, the filename is printed first, followed by the frametype.
    After these, all parameters listed in the 'keyword' item in the
    settings file will be printed

    Parameters
    ----------
    fitstbl : Table
      Contains relevant information from fits header files
    """
    msgs.info("Preparing to write out the data sorting details")
    nfiles = fitstbl['filename'].size
    # Specify which keywords to print after 'filename' and 'filetype'
    prord = ['filename', 'frametype', 'target', 'exptime', 'naxis0', 'naxis1', 'filter1', 'filter2']
    prdtp = ["char",     "char",      "char",   "double",  "int",    "int",    "char",     "char"]
    # Now insert the remaining keywords:
    fkey = settings_spect['keyword'].keys()
    for i in fkey:
        if i not in prord:
            prord.append(i)
            # Append the type of value this keyword holds
            typv = type(fitstbl[i][0])
            if typv is int or typv is np.int_:
                prdtp.append("int")
            elif isinstance(fitstbl[i][0], basestring) or typv is np.string_:
                prdtp.append("char")
            elif typv is float or typv is np.float_:
                prdtp.append("double")
            else:
                msgs.bug("I didn't expect useful headers to contain type {!s:s}".format(typv).replace('<type ', '').replace('>', ''))

    # ASCII file
    asciiord = ['filename', 'date', 'frametype', 'frameno', 'target', 'exptime', 'binning',
        'dichroic', 'dispname', 'dispangle', 'decker']
    # Generate the columns except frametype
    ascii_tbl = Table()
    badclms = []
    for pr in asciiord:
        if pr != 'frametype':
            try:  # No longer require that all of these be present
                ascii_tbl[pr] = fitstbl[pr]
            except KeyError:
                badclms.append(pr)
    # Remove
    for pr in badclms:
        asciiord.pop(asciiord.index(pr))
    # Frametype
    ascii_tbl['frametype'] = build_frametype_list(fitstbl)
    # Write
    if settings_argflag['run']['setup']:
        ascii_name = settings_argflag['run']['redname'].replace('.pypit', '.lst')
    else:
        ascii_name = settings_argflag['output']['sorted']+'.lst'
    ascii_tbl[asciiord].write(ascii_name, format='ascii.fixed_width')
    return ascii_tbl


def build_frametype_list(fitstbl):
    """

    Parameters
    ----------
    fitstbl : Table

    Returns
    -------
    ftypes : list
      List of frametype's for each frame, e.g.
        arc
        trace,pixelflat
    """
    # Now frame type
    ftypes = []
    for i in range(len(fitstbl)):
        addval = ""
        for ft in ftype_list:
            if fitstbl[ft][i]:
                if len(addval) != 0: addval += ","
                addval += ft
        ftypes.append(addval)
    # Return
    return ftypes

def match_ABBA(fitstbl, max_targ_sep=30, max_nod_sep=2):
    """

    Parameters
    ----------
    fitstbl : Table
        Contains relevant information from fits header files
    Returns
    -------
    fitstbl : Table (w/new column)

    """
    # Make coords for all observations, including calibration frames
    coords = SkyCoord(ra=fitstbl['ra'], dec=fitstbl['dec'], unit='deg')

    # Indices of files classified as science frames
    sci_idx = np.where(fitstbl['science'])[0]

    # Create a mask, currently all True
    mask = np.ones(np.sum(len(fitstbl)), dtype=bool)

    # Create empty dictionary of targets
    targets = {}

    # Set science frames to False in the mask
    mask[sci_idx] = False

    while np.any(~mask):
        idx = np.where(~mask)[0][0]
        # Find indices corresponding to any others files that match the coordinates of current science frame in consideration
        match = coords[idx].separation(coords[sci_idx]).arcsec < max_targ_sep # 30 arcsec separation is arbitrary (from 1e-2 deg) for target separation
        # Add this science target to the dictionary
        targ_name = fitstbl[idx]['target']  # Grab the target name from fitstbl
        # Create that target in dict and add corresponding indices that point to the relevant science frames
        targets[targ_name] = sci_idx[match]
        # Turn mask 'off' (i.e., to True) now that they have been grouped to this target
        mask[targets[targ_name]] = True

    # Create an empty list of filenames that will store each frame's buddy A/B frame
    AB_frame = [''] * len(fitstbl)

    for key, value in targets.items():

        files = fitstbl['filename'][value]

        # Check here that there are more than 1 files and that # of files is even
        if len(files) == 1:
            msgs.warn('Cannot perform NIR A-B reduction on targets with 1 file')
        elif len(files) % 2 != 0:
            msgs.warn('Expected an even number of files associated with target ' + key)
        #### Check for increasing time? Files are read in numerical sequential order -- should be in order of increasing time anyway..

        # Assume that the files are initially in ABBA order and proceed
        ABBA_coords = coords[value]

        # Break files into ABBA groups (includes 'remainder' if there are only 2 files)
        ABBA_groups = [ABBA_coords[i:i + 4] for i in range(0, len(ABBA_coords), 4)]
        value_groups = [value[i:i + 4] for i in range(0, len(ABBA_coords), 4)]

        for group in range(len(ABBA_groups)):
            # Warn user that if there are any groups of only 2 files, assuming they are in order of A and B
            if len(ABBA_groups[group]) == 2:
                msgs.info('Assuming these two frames are A and B frame')
            # Check that we have a 4-file, ABBA sequence
            elif len(ABBA_groups[group]) == 4:
                # Check that frames 1, 4 of an ABBA sequence are at the same nod position (A) based on their RA, DEC
                if ABBA_coords[0].separation(ABBA_coords[-1]).arcsec > max_nod_sep: # 5e-4 deg --> 1.8 arcsec --> ~2 arcsec set to be max nod separation
                    msgs.info('Are frames ' + str((group * 4) + 1) + ' and ' + str(
                        (group * 4) + 4) + ' for target ' + key + ' both A frames?')
                # Check that frames 2, 3 of an ABBA sequence are both at the same nod position (B)
                if ABBA_coords[1].separation(ABBA_coords[2]).arcsec > max_nod_sep:
                    msgs.info('Are frames ' + str((group * 4) + 2) + ' and ' + str(
                        (group * 4) + 3) + ' for target ' + key + ' both B frames?')
            else:
                msgs.error('Check number of frames for this target -- files are not grouped in ABBA or AB')

            # Create a copy of the array value_groups[group] (which gives the indices corresponding to the 4/2 ABBA/AB files in consideration)
            AB_idx_flip = np.copy(value_groups[group])
            # Flip such that, for example, (1, 2, 3, 4) --> (2, 1, 4, 3)
            AB_idx_flip[::2], AB_idx_flip[1::2] = value_groups[group][1::2], value_groups[group][::2]

            # Fill in AB_frame list
            for i in range(len(value_groups[group])):
                AB_frame[value_groups[group][i]] = fitstbl['filename'][AB_idx_flip[i]]

    fitstbl['AB_frame'] = AB_frame

    return fitstbl

'''
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

    """
    # Open a VOTable for writing
    votable = votable.tree.VOTableFile()
    resource = votable.tree.Resource()
    votable.resources.append(resource)
    table = votable.tree.Table(votable)
    resource.tables.append(table)
    # Define VOTable fields
    tabarr=[]
    # Insert the filename and filetype first
    for i in range(len(prord)):
        tabarr.append(votable.tree.Field(votable, name=prord[i], datatype=prdtp[i], arraysize="*"))
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
    """

    # ASCII file
    asciiord = ['filename', 'date', 'frametype', 'frameno', 'target', 'exptime', 'binning',
        'dichroic', 'dispname', 'dispangle', 'decker']
    # Generate the columns except frametype
    ascii_tbl = Table()
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
'''


def match_logic(ch, tmtch, fitstbl, idx):
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
      If tmtch begins with a "|", the match compares to the science frame
      else the value is added to the science frame
    fitstbl : Table
    idx : int
      Science index

    Returns
    -------
    w : ndarray, bool
      True/False for the rows in fitstbl satisfying the condition
    """
    if tmtch == "any":   # Anything goes
        w = np.ones_like(fitstbl, dtype=bool)
    elif tmtch == "''":  # Header value must match that of science
        w = fitstbl[ch] == fitstbl[ch][idx]
    elif tmtch[0] in ['=','<','>','|']: # Numerics
        mtch = np.float64(fitstbl[ch][idx]) + float(
            ''.join(c for c in tmtch if c not in ['=', '<', '>', '|']))
        operand = ''.join(c for c in tmtch if c in ['=', '<', '>'])
        if operand == '=':
            operand += '='
        #
        if tmtch[0] != '|':
            w = eval('fitstbl[ch].data.astype(np.float64) {:s} {:f}'.format(operand, mtch))
        else:
            w = eval('np.abs(fitstbl[ch].data.astype(np.float64) - np.float64(fitstbl[ch][idx])) {:s} {:f}'.format(operand, mtch))
    elif tmtch[0:2] == '%,':  # Splitting a header keyword
        splcom = tmtch.split(',')
        debugger.set_trace()
        spltxt, argtxt, valtxt = splcom[1], np.int(splcom[2]), splcom[3]
        tspl = []
        for sp in fitstbl[ch]:
            tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
            tmpspl = re.split(tmpspl, sp)
            if len(tmpspl) < argtxt+1:
                tspl.append("-9999999")
            else:
                tspl.append(tmpspl[argtxt])
        tspl = np.array(tspl)
        #                        debugger.set_trace()
        tmpspl = str(re.escape(spltxt)).replace("\\|", "|")
        tmpspl = re.split(tmpspl, fitstbl[ch][idx])
        msgs.warn("HAS NOT BEEN DEVELOPED SINCE THE SetupClass refactor;  no test case..")
        debugger.set_trace()  # HAS NOT BEEN DEVELOPED SINCE THE SetupClass refactor;  no test case..
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
    # Return
    return w


'''
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
    fitsdict
    idx : int
      Science index


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
'''


def match_to_science(fitstbl, settings_spect, settings_argflag):
    """
    For a given set of identified data, match calibration frames to science frames

    Parameters
    ----------
    fitstbl : Table
      Contains relevant information from fits header files
    settings_spect : dict
    settings_argflag : dict

    Returns
    -------
    fitstbl : Table
      Updated with failures and sci_ID columns
    """
    msgs.info("Matching calibrations to Science frames")

    # New columns
    fitstbl['failures'] = False
    fitstbl['sci_ID'] = 0

    # Loop on science frames
    #while i < nSCI:
    for ss, sci_idx in enumerate(np.where(fitstbl['science'])[0]):
        msgs.info("=================================================")
        msgs.info("Matching calibrations to {:s}: {:s}".format(
                fitstbl['target'][sci_idx], fitstbl['filename'][sci_idx]))

        # Science ID (trivial but key to the bit-wise that follows)
        fitstbl['sci_ID'][sci_idx] = 2**ss

        # Find matching (and nearby) calibration frames
        for ftag in ftype_list:
            gd_match = fitstbl[ftag].data.copy()  # All work for starters
            if ftag in ['science', 'unknown']:
                continue

            # bias/dark check to make sure we need to find matching frames
            if ftag == 'dark' and settings_argflag['bias']['useframe'] != 'dark':
                msgs.info("  Dark frames not required.  Not matching..")
                continue
            if ftag == 'bias' and settings_argflag['bias']['useframe'] != 'bias' and not settings_argflag['reduce']['badpix']:
                msgs.info("  Bias frames not required.  Not matching..")
                continue

            # How many matching frames are required?  This is instrument specific
            if settings_argflag['run']['setup']:
                numfr = 1 if ftag == 'arc' else 0
            else:
                numfr = settings_spect[ftag]['number']

            # If not required and not doing setup, continue
            if (numfr == 0) and (not settings_argflag['run']['setup']):
                msgs.info("   No {0:s} frames are required.  Not matching..".format(ftag))
                continue

            # Now go ahead and match the frames
            if 'match' not in settings_spect[ftag].keys() and (not settings_argflag['run']['setup']):
                msgs.error("Need match criteria for {0:s}!!".format(ftag))
            elif 'match' not in settings_spect[ftag].keys():
                msgs.info("No matching criteria for {0:s} frames with this instrument".format(ftag))
            else:
                chkk = settings_spect[ftag]['match'].keys()
                for ch in chkk:
                    tmtch = settings_spect[ftag]['match'][ch]
                    gd_match &= match_logic(ch, tmtch, fitstbl, sci_idx)

            # Find the time difference between the calibrations and science frames
            if settings_spect['fits']['calwin'] > 0.0:
                tdiff = np.abs(fitstbl['time']-fitstbl['time'][sci_idx])
                gd_match &= tdiff <= settings_spect['fits']['calwin']

            # Now find which of the remaining n are the appropriate calibration frames
            nmatch = np.sum(gd_match)
            if settings_argflag['output']['verbosity'] == 2:
                msgs.info("  Found {0:d} {1:s} frame for {2:s} ({3:d} required)".format(
                    nmatch, ftag, fitstbl['target'][sci_idx], numfr))

            # Have we identified enough of these calibration frames to continue?
            if nmatch < np.abs(numfr):
                code = match_warnings(ftag, nmatch, numfr, fitstbl['target'][sci_idx],
                                          settings_argflag)
                if code == 'break':
                    fitstbl['failure'][sci_idx] = True
                    fitstbl['sci_ID'][sci_idx] = -1  # This might break things but who knows..
                    break
            else:
                wa = np.where(gd_match)[0]
                # Select the closest calibration frames to the science frame
                tdiff = np.abs(fitstbl['time'][wa]-fitstbl['time'][sci_idx])
                wa = wa[np.argsort(tdiff)]
                #if ftag == 'bias':
                #    debugger.set_trace()
                if (settings_argflag['run']['setup']) or (numfr < 0):
                    fitstbl['sci_ID'][wa] |= 2**ss  # Flip the switch (if need be)
                else:
                    fitstbl['sci_ID'][wa[:numfr]] |= 2**ss  # Flip the switch (if need be)

    msgs.info("Science frames successfully matched to calibration frames")

    # How to do this if statement only if '--custom' is on?
    if settings_argflag['run']['spectrograph'] == 'keck_nirspec':
        fitstbl = match_ABBA(fitstbl)

    return fitstbl

'''
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
    #ftag = armeta.allowed_file_types()
    ftag = ftype_list
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
'''

def match_warnings(ftag, nmatch, numfr, target, settings_argflag):
    """
    Provide match warnings

    Parameters
    ----------
    ftag : str
      frametype, e.g. bias
    nmatch : int
    numfr : int
    target : str
      Name of the target
    settings_argflag : dict

    Returns
    -------
    code : str
      'None' = no further action required
    """
    code = 'None'
    msgs.warn("  Only {0:d}/{1:d} {2:s} frames for {3:s}".format(nmatch, numfr, ftag, target))
    # Errors for insufficient BIAS frames
    if settings_argflag['bias']['useframe'].lower() == ftag:
        msgs.warn("Expecting to use bias frames for bias subtraction. But insufficient frames found.")
        msgs.warn("Either include more frames or modify bias method" + msgs.newline() +
                  "  e.g.:   bias useframe overscan")
        msgs.error("Unable to continue")
    # Errors for insufficient PIXELFLAT frames
    if ftag == 'pixelflat' and settings_argflag['reduce']['flatfield']['perform'] and (
                settings_argflag['reduce']['flatfield']['useframe'] == 'pixelflat'):
        if settings_argflag['reduce']['masters']['force']:
            msgs.warn("Fewer pixelflat frames than expected for {0:s}, but will use MasterFrames".format(target))
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
    if ftag == 'trace' and settings_argflag['reduce']['flatfield']['perform']:
        if settings_argflag['reduce']['masters']['force']:
            msgs.warn("Fewer traceflat frames than expected for {0:s}, but will use MasterFrames".format(target))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    # Errors for insufficient standard frames
    if ftag == 'standard' and settings_argflag['reduce']['calibrate']['flux']:
        if settings_argflag['reduce']['masters']['force']:
            msgs.warn("No standard star frames for {0:s}, but will use MasterFrames".format(target))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    # Errors for insufficient ARC frames
    if ftag == 'arc' and settings_argflag['reduce']['calibrate']['wavelength']:
        if settings_argflag['run']['setup']:
            msgs.warn("No arc frame for {0:s}. Removing it from list of science frames".format(target))
            msgs.warn("Add an arc and rerun one if you wish to reduce this with PYPIT!!")
            return 'break'
        elif settings_argflag['reduce']['masters']['force']:
            msgs.warn("No arc frame for {0:s}, but will use MasterFrames".format(target))
        else:
            msgs.error("Unable to continue without more {0:s} frames".format(ftag))
    return code

'''
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
'''

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


def make_dirs(settings_argflag):
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
    newdir = "{0:s}/{1:s}".format(currDIR, settings_argflag['run']['directory']['science'])
    if os.path.exists(newdir):
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
        if not settings_argflag['output']['overwrite']:
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
    #w = filesort['science']
    #sci_targs = np.array(list(set(fitsdict['target'][w])))
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
    newdir = "{:s}/{:s}_{:s}".format(currDIR, settings_argflag['run']['directory']['master'],
                                     settings_argflag['run']['spectrograph'])
    if os.path.exists(newdir):
        if not settings_argflag['output']['overwrite']:
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
    newdir = "{0:s}/{1:s}".format(currDIR, settings_argflag['run']['directory']['qa'])
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
    return


def dummy_fitstbl(nfile=10, spectrograph='shane_kast_blue', directory='./', notype=False):
    """
    Generate a dummy fitstbl for testing

    Parameters
    ----------
    nfile : int, optional
      Number of files to mimic
    spectrograph : str, optional
      Name of spectrograph to mimic
    notype : bool (optional)
      If True, do not add image type info to the fitstbl

    Returns
    -------
    fitstbl : Table

    """
    fitsdict = dict({'directory': [], 'filename': [], 'utc': []})
    fitsdict['utc'] = ['2015-01-23']*nfile
    fitsdict['directory'] = [directory]*nfile
    fitsdict['filename'] = ['b{:03d}.fits'.format(i) for i in range(nfile)]
    fitsdict['date'] = ['2015-01-23T00:{:02d}:11.04'.format(i) for i in range(nfile)]  # Will fail at 60
    fitsdict['time'] = [(1432085758+i*60)/3600. for i in range(nfile)]
    fitsdict['target'] = ['Dummy']*nfile
    fitsdict['ra'] = ['00:00:00']*nfile
    fitsdict['dec'] = ['+00:00:00']*nfile
    fitsdict['exptime'] = [300.] * nfile
    fitsdict['naxis0'] = [2048] * nfile
    fitsdict['naxis1'] = [2048] * nfile
    fitsdict['dispname'] = ['600/4310'] * nfile
    fitsdict['dichroic'] = ['560'] * nfile
    fitsdict['dispangle'] = ['none'] * nfile
    fitsdict["binning"] = ['1x1']*nfile
    fitsdict["airmass"] = [1.0]*nfile
    #
    if spectrograph == 'shane_kast_blue':
        fitsdict['numamplifiers'] = [1] * nfile
        fitsdict['naxis0'] = [2112] * nfile
        fitsdict['naxis1'] = [2048] * nfile
        fitsdict['slitwid'] = [1.] * nfile
        fitsdict['slitlen'] = ['none'] * nfile
        # Lamps
        for i in range(1,17):
            fitsdict['lampstat{:02d}'.format(i)] = ['off'] * nfile
        fitsdict['exptime'][0] = 0        # Bias
        fitsdict['lampstat06'][1] = 'on'  # Arc
        fitsdict['exptime'][1] = 30       # Arc
        fitsdict['lampstat01'][2] = 'on'  # Trace, pixel, slit flat
        fitsdict['lampstat01'][3] = 'on'  # Trace, pixel, slit flat
        fitsdict['exptime'][2] = 30     # flat
        fitsdict['exptime'][3] = 30     # flat
        fitsdict['ra'][4] = '05:06:36.6'  # Standard
        fitsdict['dec'][4] = '52:52:01.0'
        fitsdict['airmass'][4] = 1.2
        fitsdict['ra'][5] = '07:06:23.45' # Random object
        fitsdict['dec'][5] = '+30:20:50.5'
        fitsdict['decker'] = ['0.5 arcsec'] * nfile
    elif spectrograph == 'none':
        pass
    # arrays
    for k in fitsdict.keys():
        fitsdict[k] = np.array(fitsdict[k])
    # Table me
    fitstbl = Table(fitsdict)
    # Image typing
    if not notype:
        for ftype in ftype_list:
            fitstbl[ftype] = np.zeros(len(fitstbl), dtype=bool)
        if spectrograph == 'shane_kast_blue':
            fitstbl['sci_ID'] = 1  # This links all the files to the science object
            fitstbl['bias'][0] = True
            fitstbl['arc'][1] = True
            fitstbl['trace'][2:4] = True
            fitstbl['pixelflat'][2:4] = True
            fitstbl['standard'][4] = True
            fitstbl['science'][5:] = True
    # Return
    return fitstbl
