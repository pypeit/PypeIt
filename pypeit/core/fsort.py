""" Routines for sorting data to be reduced by PYPIT"""
from __future__ import (print_function, absolute_import, division, unicode_literals)


import os
import re
import sys
import shutil

import numpy as np

from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units

from pypeit import msgs
#from pypeit import arparse as settings
from pypeit import utils
from pypeit.core.flux import find_standard_file
from pypeit import debugger

from pypeit.bitmask import BitMask

class FrameTypeBitMask(BitMask):
    """
    Define a bitmask to set the frame types.

    Frame types can be arc, bias, dark, pinhole, pixelflat, science,
    standard, or trace.
    """
    def __init__(self):
        frame_types = {         'arc': 'Arc lamp observation used for wavelength calibration',
                               'bias': 'Bias readout for detector bias subtraction',
                               'dark': 'Shuttered exposure to measure dark current',
                            'pinhole': 'Pinhole observation used for tracing slit centers',
                          'pixelflat': 'Flat-field exposure used for pixel-to-pixel response',
                            'science': 'On-sky observation of a primary target',
                           'standard': 'On-sky observation of a flux calibrator',
                              'trace': 'High-count exposure used to trace slit positions'
                      }
        super(FrameTypeBitMask, self).__init__(list(frame_types.keys()),
                                               descr=list(frame_types.values()))

## TODO: (KBW) You know my comment about this...
#ftype_list = [     # NOTE:  arc must be listed first!
#    'arc',         # Exposure of one or more 
#    'bias',        # Exposure for assessing detector bias (usually 0s)
#    'dark',        # Exposure for assessing detector dark current
#    'pinhole',     # Exposure for tracing the orders or slits
#    'pixelflat',   # Exposure for assessing pixel-to-pixel variations
#    'science',   # Exposure on one or more science targets
#    'standard',    # Exposure on a 'standard star' used for flux calibration
#    'trace',       # Exposure for tracing slits or echelle orders (usually twilight sky or flat lamp)
#    'unknown',     # Unknown..
#]

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


def list_of_files(fitstbl, ftype, sci_ID):
    """
    Generate a list of filenames with path for a given frametype and sci_ID

    Parameters
    ----------
    fitstbl : Table
    ftype : str
    sci_ID : int

    Returns
    -------
    file_list : list

    """
    file_list = []
    idx = ftype_indices(fitstbl, ftype, sci_ID)
    # Saving the match_frames algorithm (perhaps only for posterity!)
    #sframes = arsort.match_frames(frames, settings.argflag['trace']['combine']['match'], frametype='trace', satlevel=settings.spect[dnum]['saturation']*settings.spect['det'][det-1]['nonlinear'])
    for ii in idx:
        file_list.append(os.path.join(fitstbl['directory'][ii], fitstbl['filename'][ii]))
    # Return
    return file_list


def get_type_names(type_bits, bitmask=None, join=True):
    """
    Use the type bits to get the type names for each frame.

    Args:
        type_bits (int):
            The bit mask for each frame.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            The bit mask used to pull out the bit names.  Uses
            :class:`FrameTypeBitMask` by default.

    Returns:
        list: List of the frame types for each frame.  Each frame can
        have multiple types, meaning the 2nd axis is not necessarily the
        same length for all frames.
    """
    bm = FrameTypeBitMask() if bitmask is None else bitmask
    out = []
    for b in type_bits:
        n = bm.flagged_bits(b)
        if len(n) == 0:
            n = ['None']
        out += [','.join(n)] if join else [n]
    return out


def get_frame_types(spectrograph, fitstbl, flag_unknown=False, user=None, useIDname=False):
    """
    Generate a table of frame types from the input metadata object.

    Args:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The spectrograph used to collect the data save to each file.
            The class is used to provide the header keyword data to
            include in the table and specify any validation checks.
        fitstbl (:obj:`astropy.table.Table`):
            Table with the fits file metadata.
        flag_unknown (:obj:`bool`, optional):
            Instead of crashing out if there are unidentified files,
            leave without a type and continue.
        user (:obj:`dict`, optional):
            A dictionary with the types designated by the user.  The
            file name and type are expected to be the key and value of
            the dictionary, respectively.  The number of keys therefore
            *must* match the number of files in the provided `fitstbl`.
            For frames that have multiple types, the types should be
            provided as a string with comma-separated types.
        useIDname (:obj:`bool`, optional):
            Use ID name in the Header to image type

    Returns:
        :obj:`astropy.table.Table`: A Table with two columns, the type
        names and the type bits.  See :class:`FrameTypeBitMask` for the
        allowed frame types.
    """
    msgs.info("Typing files")

    # Prep the bits vector
    numfiles = fitstbl['filename'].size
    bm = FrameTypeBitMask()
    type_bits = np.zeros(numfiles, dtype=bm.minimum_dtype())

    # Use the user-defined frame types from the input dictionary
    if user is not None:
        if len(user.keys()) != numfiles:
            raise ValueError('The user-provided dictionary does not match the fitstbl.')
        msgs.info('Using user-provided frame types.')
        for ifile,ftypes in ftdict.items():
            indx = fitstbl['filename'] == ifile
            type_bits[indx] = bm.turn_on(type_bits[indx], flag=ftypes.split(','))
        return Table({'frametype': get_type_names(type_bits, bitmask=bm), 'framebit': type_bits})

    # Loop over the frame types
    for i, ftype in enumerate(bm.keys()):

        # Initialize: Flag frames with the correct ID name or start by
        # flagging all as true
        indx = fitstbl['idname'] == spectrograph.idname(ftype) if useIDname \
                    else np.ones(numfiles, dtype=bool)

        # Include a combination of instrument-specific checks using
        # combinations of the full set of metadata
        indx &= spectrograph.check_ftype(ftype, fitstbl)

        # Turn on the relevant bits
        type_bits[indx] = bm.turn_on(type_bits[indx], flag=ftype)

    # Find the nearest standard star to each science frame
    # TODO: Should this be 'standard' or 'science' or both?
    if 'ra' not in fitstbl.keys() or 'dec' not in fitstbl.keys():
        msgs.warn('Cannot associate standard with science frames without sky coordinates.')
    else:
        indx = bm.flagged(type_bits, flag='standard')
        for b, f, ra, dec in zip(type_bits[indx], fitstbl['filename'][indx], fitstbl['ra'][indx],
                                 fitstbl['dec'][indx]):
            if ra == 'None' or dec == 'None':
                msgs.warn('RA and DEC must not be None for file:' + msgs.newline() + f)
                msgs.warn('The above file could be a twilight flat frame that was' + msgs.newline()
                        + 'missed by the automatic identification.')
                b = bm.turn_off(b, flag='standard')
                continue

            # If an object exists within 20 arcmins of a listed standard,
            # then it is probably a standard star
            foundstd = find_standard_file(ra, dec, check=True)
            b = bm.turn_off(b, flag='science' if foundstd else 'standard')

#    # Make any forced changes
#    skeys = settings_spect['set'].keys()
#    if len(skeys) > 0:
#        msgs.info("Making forced file identification changes")
#        msgs.warn("Note that the image will have *only* the specified type")
#        for sk in skeys:
#            for jj in settings_spect['set'][sk]:
#                idx = np.where(fitstbl['filename']==jj)[0]
#                # Zero out the others
#                for ftype in ftype_list:
#                    filetypes[ftype][idx] = False
#                # And set
#                filetypes[sk][idx] = True

    # Find the files without any types
    indx = np.invert(bm.flagged(type_bits))
    if np.any(indx):
        msgs.info("Couldn't identify the following files:")
        for f in fitstbl['filename'][indx]:
            msgs.info(f)
        if not flag_unknown:
            msgs.error("Check these files before continuing")

    # Now identify the dark frames
    # TODO: Move this to Spectrograph.check_ftype?
    indx = bm.flagged(type_bits, flag='bias') \
                    & (fitstbl['exptime'].data.astype(float) > spectrograph.minexp)
    type_bits[indx] = bm.turn_on(type_bits[indx], 'dark')

    # Finish up (note that this is called above if user is not None!)
    msgs.info("Typing completed!")
    return Table({'frametype': get_type_names(type_bits, bitmask=bm), 'framebit': type_bits})


def chk_all_conditions(fitstbl, cond_dict):
    """ Loop on the conditions for this given file type

    Parameters
    ----------
    fitstbl : Table
    cond_dict : dict

    Returns
    -------
    gd_chk : ndarray (bool)
      True = Passes all checks
    """
    gd_chk = np.ones(len(fitstbl), dtype=bool)
    # Loop on the items to check
    chkk = cond_dict.keys()
    for ch in chkk:
        if ch[0:9] == 'condition':
            # Deal with a conditional argument
            conds = re.split("(\||\&)", cond_dict[ch])
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
                gd_chk = gd_chk & (np.char.strip(fitstbl[ch]) == cond_dict[ch])
            else:
                gd_chk = gd_chk & (fitstbl[ch] == cond_dict[ch])
    # Return
    return gd_chk


def chk_condition(fitstbl, cond):
    """
    Code to perform condition.  A bit messy so a separate definition
    was generated.
    Create an exposure class for every science frame

    Parameters
    ----------
    fitsdict : Table
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
        ntmp = fitstbl[tcond[0]] <= float(tcond[1])
    elif ">=" in cond:
        tcond = cond.split(">=")
        ntmp = fitstbl[tcond[0]] >= float(tcond[1])
    elif "!=" in cond:
        tcond = cond.split("!=")
        if 'int' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] != int(tcond[1])
        elif 'float' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] != float(tcond[1])
        else:
            ntmp = fitstbl[tcond[0]] != tcond[1]
    elif "<" in cond:
        tcond = cond.split("<")
        ntmp = fitstbl[tcond[0]] < float(tcond[1])
    elif ">" in cond:
        tcond = cond.split(">")
        ntmp = fitstbl[tcond[0]] > float(tcond[1])
    elif "=" in cond:
        tcond = cond.split("=")
        if 'int' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] == int(tcond[1])
        elif 'float' in fitstbl[tcond[0]].dtype.name:
            ntmp = fitstbl[tcond[0]] == float(tcond[1])
        else:
            ntmp = fitstbl[tcond[0]] == tcond[1]
    else:
        ntmp = None
    return ntmp


def write_fitstbl(fitstbl, ofile, columns=None):
    """
    Write the metadata for the files to reduce.

    The table is written with the filename and frametype columns first.
    All remaining columns, or a subset of them selected by `columns`,
    follow these first two.

    .. todo::
        Force a fault if 'filename' and 'frametype' are not columns in
        `fitstbl`?

    Args:
        fitstbl (:obj:`astropy.table.Table`):
            Table with the metadata for the fits files to reduce.
        ofile (:obj:`str`):
            Name for the output file.  Extension of the file should be
            '.lst'.
        columns (:obj:`list`, optional):
            A list of columns to include in the output file.  If None,
            all columns in `fitstbl` are written.
    
    Raises:
        ValueError:
            Raised if the columns to include are not unique.
    """
    msgs.info('Writing fits file metadata to {0}.'.format(ofile))

    # Set the columns to include and check that they are unique
    _columns = list(fitstbl.keys()) if columns is None else columns
    if len(numpy.unique(_columns)) != len(_columns):
        # TODO: A warning may suffice...
        raise ValueError('Column names must be unique!')

    # Force the filename and frametype columns to go first
    col_order = [ 'filename', 'frametype' ]
    col_order.append(list(set(_columns) - set(col_order)))

    # Remove any columns that don't exist
    for c in col_order:
        if c not in fitstbl.keys():
            msgs.warn('{0} is not a valid column!  Removing from output.'.format(c))
            col_order.remove(c)

    # Write the output
    fitstbl[col_order].write(ofile, format='ascii.fixed_width')


# TODO: Deprecated, but left here until I'm sure I've removed the calls
# to it.
def write_lst(fitstbl, skeys, pypeit_filename, setup=False,
              sort_dir=None):
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
    for i in skeys:
        if i not in prord:
            prord.append(i)
            # Append the type of value this keyword holds
            typv = type(fitstbl[i][0])
            if typv is int or typv is np.int_:
                prdtp.append("int")
            elif isinstance(fitstbl[i][0], str) or typv is np.string_:
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
    if setup:
        ascii_name = pypeit_filename.replace('.pypeit', '.lst')
    else:
        ascii_name = sort_dir+'.lst'
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
    elif tmtch == '':  # Header value must match that of science
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


def match_to_science(calib_par, match_dict, fitstbl, calwin, setup=False, verbose=True,
                     match_nods=False):
    """
    For a given set of identified data, match calibration frames to science frames

    Parameters
    ----------
    fitstbl : Table
      Contains relevant information from fits header files
    settings_spect : dict
    setup: bool, optional
      Running in setup mode?
    flux_calibrate : bool, optional
      Do checks related to flux calibration
    wave_calib : str
    calwin : float

    Returns
    -------
    fitstbl : Table
      Updated with failures and sci_ID columns
    """
    msgs.info("Matching calibrations to Science frames")

    # New columns

    fitstbl['failures'] = False
    fitstbl['sci_ID'] = 0

    # Frame type masking
    bm = FrameTypeBitMask()

    # Loop on science frames
    for ss, sci_idx in enumerate(np.where(fitstbl['science'])[0]):
        msgs.info("=================================================")
        msgs.info("Matching calibrations to {:s}: {:s}".format(
                fitstbl['target'][sci_idx], fitstbl['filename'][sci_idx]))

        # Science ID (trivial but key to the bit-wise that follows)
        fitstbl['sci_ID'][sci_idx] = 2**ss

        # Find matching (and nearby) calibration frames
        for ftag in bm.keys():
            # Skip science frames
            if ftag == 'science':
                continue

            # bias/dark check to make sure we need to find matching frames
            if ftag == 'dark' and calib_par['biasframe']['useframe'] != 'dark':
                msgs.info("  Dark frames not required.  Not matching..")
                continue
            if ftag == 'bias' and calib_par['biasframe']['useframe'] != 'bias' \
                        and not calib_par['badpix']:
                msgs.info("  Bias frames not required.  Not matching..")
                continue

            # Find the files of this type
            gd_match = bm.flagged(fitstbl['framebit'], flag=ftag)

            # How many matching frames are required?  This is instrument specific
            numfr = (1 if ftag == 'arc' else 0) if setup \
                        else calib_par['{0}frame'.format(ftag)]['number']

            # If not required and not doing setup, continue
            if numfr == 0 and not setup:
                msgs.info("   No {0:s} frames are required.  Not matching..".format(ftag))
                continue

            # Now go ahead and match the frames
            if 'match' not in match_dict[ftag].keys() and (not setup):
                msgs.error("Need match criteria for {0:s}!!".format(ftag))
            elif 'match' not in match_dict[ftag].keys():
                msgs.info("No matching criteria for {0:s} frames with this instrument".format(ftag))
            else:
                chkk = match_dict[ftag]['match'].keys()
                for ch in chkk:
                    tmtch = match_dict[ftag]['match'][ch]
                    gd_match &= match_logic(ch, tmtch, fitstbl, sci_idx)

            # Find the time difference between the calibrations and science frames
            if calwin > 0.0 and 'time' in fitstbl.keys():
                tdiff = np.abs(fitstbl['time']-fitstbl['time'][sci_idx])
                gd_match &= tdiff <= calwin

            # Now find which of the remaining n are the appropriate calibration frames
            nmatch = np.sum(gd_match)
            if verbose and 'target' in fitstbl.keys():
                msgs.info('  Found {0:d} {1:s} frame for {2:s} ({3:d} required)'.format(
                                nmatch, ftag, fitstbl['target'][sci_idx], numfr))

            # Have we identified enough of these calibration frames to continue?
            if nmatch < np.abs(numfr):
                code = match_warnings(calib_par, ftag, nmatch, numfr,
                                      fitstbl['target'][sci_idx] if 'target' in fitstbl.keys()
                                            else 'unknown')
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
                if setup or (numfr < 0):
                    fitstbl['sci_ID'][wa] |= 2**ss  # Flip the switch (if need be)
                else:
                    fitstbl['sci_ID'][wa[:numfr]] |= 2**ss  # Flip the switch (if need be)

    msgs.info("Science frames successfully matched to calibration frames")

    # Return with nods matched if requested
    return match_ABBA(fitstbl) if match_nods else fitstbl

#    # How to do this if statement only if '--custom' is on?
#    if spectrograph.spectrograph == 'keck_nirspec':
#        fitstbl = match_ABBA(fitstbl)


def insufficient_frame_error(frametype):
    msgs.error('Insufficient {0} frames found. Include more frames, '.format(frametype)
                + 'reduce the required amount by setting'
                + msgs.newline() + '[calibrations]'
                + msgs.newline() + '    [[{0}frame]]'.format(frametype)
                + msgs.newline() + '        number = XX'
                + msgs.newline() + 'in the pypeit file, or specify a specific'
                + 'pixelflat file by setting'
                + msgs.newline() + '[calibrations]'
                + msgs.newline() + '    [[{0}frame]]'.format(frametype)
                + msgs.newline() + '        useframe = XX'
                + msgs.newline() + 'in the pypeit file')


def match_warnings(calib_par, ftag, nmatch, numfr, target, setup=False):
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

    # TODO: Why does number of pixelflat, trace, and standard not matter
    # if you're not flat-fielding the data?  Particularly for trace...
    flatfield = calib_par['flatfield']['method'] is not None

    # Errors for insufficient BIAS frames
    if calib_par['biasframe']['useframe'].lower() == ftag:
        insufficient_frame_error(ftag)

    # Errors for insufficient PIXELFLAT frames
    if ftag == 'pixelflat' and flatfield and calib_par['flatfield']['frame'] == 'pixelflat':
        if calib_par['masters'] == 'force':
            msgs.warn('Fewer {0} frames than expected for {1}'.format(ftag, target)
                      +', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient PINHOLE frames
    if ftag == 'pinhole':
        insufficient_frame_error(ftag)

    # Errors for insufficient TRACE frames
    if ftag == 'trace' and flatfield:
        if calib_par['masters'] == 'force':
            msgs.warn('Fewer {0} frames than expected for {1}'.format(ftag, target)
                      +', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient standard frames
    if ftag == 'standard' and flatfield:
        if calib_par['masters'] == 'force':
            msgs.warn('No {0} frames for {1}'.format(ftag, target)
                      + ', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    # Errors for insufficient ARC frames
    if ftag == 'arc' and (calib_par['wavelengths']['reference'] not in ['pixel', 'sky']):
        if setup:
            msgs.warn('No {0} frame for {1}. '.format(ftag, target)
                      + 'Removing it from list of science frames.  Add an arc and rerun if '
                      + 'you wish to reduce this with PYPIT!!')
            return 'break'
        elif calib_par['masters'] == 'force':
            msgs.warn('No {0} frames for {1}'.format(ftag, target)
                      + ', but will use MasterFrames.')
        else:
            insufficient_frame_error(ftag)

    return code


def match_frames(frames, criteria, frametype='<None>', satlevel=None):
    """
    identify frames with a similar appearance (i.e. one frame appears to be a scaled version of another).
    """

    prob = utils.erf(criteria/np.sqrt(2.0))[0]
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


def make_dirs(spectrograph, caldir, scidir, qadir, redux_path=None, overwrite=False):
    """
    Make the directories for the pypeit output.

    .. todo::
        I think this should just fault if the directories exist and
        `overwrite` is False.

    Args:
        spectrograph (str):
            The name of the spectrograph that provided the data to be
            reduced.
        caldir (str):
            The directory to use for saving the master calibration
            frames.
        scidir (str):
            The directory to use for the main reduction output files.
        qadir (str):
            The directory to use for the quality assessment output.
        redux_dir (str):
            The directory for the reductions
        overwrite(:obj:`bool`, optional):
            Flag to overwrite any existing files/directories.
    """

    # First, get the current working directory
    if redux_path is None:
        redux_path = os.getcwd()
    msgs.info("Creating Science directory")
    newdir = os.path.join(redux_path, scidir)
    if os.path.exists(newdir):
        msgs.info("The following directory already exists:"+msgs.newline()+newdir)
        if not overwrite:
            rmdir = ''
            while os.path.exists(newdir):
                while rmdir != 'n' and rmdir != 'y' and rmdir != 'r':
                    rmdir = input(msgs.input() + 'Remove this directory and its contents?'
                                  + '([y]es, [n]o, [r]ename) - ')
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
    newdir = os.path.join(redux_path, caldir+'_'+spectrograph)
    if os.path.exists(newdir):
        if not overwrite:
            msgs.info("The following directory already exists:"+msgs.newline()+newdir)
            rmdir = ''
            while rmdir != 'n' and rmdir != 'y':
                rmdir = input(msgs.input() + 'Remove this directory and its contents?'
                              '([y]es, [n]o) - ')
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
    # TODO: I'd rather that this still consider overwrite and fault
    # instead of just proceeding
    msgs.info("Creating QA directory")
    newdir = os.path.join(redux_path, qadir)
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
    fitsdict['filename'] = ['b{:03d}.fits.gz'.format(i) for i in range(nfile)]
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
    fitstbl['instrume'] = spectrograph
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
