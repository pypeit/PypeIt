"""
Routines for matching frames to certain types or each other.
"""
from __future__ import print_function, absolute_import, division, unicode_literals

import os
import re

import numpy as np

from pypeit import msgs
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

    def type_names(self, type_bits, join=True):
        """
        Use the type bits to get the type names for each frame.

        .. todo::
            - This should probably be a general function in
              :class:`pypeit.bitmask.BitMask`
    
        Args:
            type_bits (int, list, numpy.ndarray):
                The bit mask for each frame.
            bitmask (:class:`pypeit.bitmask.BitMask`, optional):
                The bit mask used to pull out the bit names.  Uses
                :class:`FrameTypeBitMask` by default.
            join (:obj:`bool`, optional):
                Instead of providing a list of type names for items with
                multiple bits tripped, joint the list into a single,
                comma-separated string.
    
        Returns:
            list: List of the frame types for each frame.  Each frame can
            have multiple types, meaning the 2nd axis is not necessarily the
            same length for all frames.
        """
        _type_bits = np.atleast_1d(type_bits)
        out = []
        for b in _type_bits:
            n = self.flagged_bits(b)
            if len(n) == 0:
                n = ['None']
            out += [','.join(n)] if join else [n]
        return out[0] if isinstance(type_bits, np.integer) else out
    

def check_frame_exptime(exptime, exprng):
    """
    Check that the exposure time is within the provided range.
        
    Args:
        exptime (numpy.ndarray):
            Exposure times to check.
        exprng (array-like):
            An array with the minimum and maximum exposure.  The limits
            are *exclusive* and a limit of None means there is no limit.
        
    Returns:
        numpy.ndarray: A boolean array that is True for all times within
        the provided range.
        
    Raises:
        ValueError:
            Raised if the length of `exprng` is not 2.
    """
    if exprng is None:
        return np.ones(len(exptime), dtype=bool)
    if len(exprng) != 2:
        raise ValueError('exprng must have two elements.')
    indx = np.ones(len(exptime), dtype=bool) if exprng[0] is None \
                else exptime > exprng[0]
    if exprng[1] is not None:
        indx &= (exptime < exprng[1])
    return indx
    

def group_AB_frames(file_list, targets, coords, max_nod_sep=2):
    """
    Group files into a ABBA or AB sequences.

    Args:
        file_list (:obj:`list`):
            A list of file names.
        targets (:obj:`dict`):
            A dictionary that matches each file to a unique target name.
            The target name can be one of the files in the file list.
        coords (:class:`astropy.coordinates.SkyCoord`):
            The coordinates of all the exposures.  Number of coordinates
            should match the number of files.
        max_nod_sep (:obj:`int`, optional):
            The maximum separation (arcsec) between the 1st and 4th
            frame sky coordinates in the ABBA sequence that is allowed
            when identifying the sequence.  Note that the default (2
            arcsec) is arbitrary.
    
    Returns:
        list:
            A list that matches the length of the input list of files.
            Each file in an AB or ABBA sequence is identified with it's
            pair in the sequence.
    """

    AB_frame = [''] * len(file_list)

    for key, value in targets.items():
        files = file_list[value]

        # Check here that there are more than 1 files and that the
        # number of files is even
        if len(files) == 1:
            msgs.warn('Cannot perform ABBA reduction on targets with 1 file')
        elif len(files) % 2 != 0:
            msgs.warn('Expected an even number of files associated with target ' + key)

        # TODO: Check for increasing time? Files are read in numerical
        # sequential order -- should be in order of increasing time
        # anyway..

        # Assume that the files are initially in ABBA order and proceed
        ABBA_coords = coords[value]

        # Break files into ABBA groups (includes remainder if there are only 2 files)
        file_groups = [files[i:i+4] for i in range(0,len(files),4)]
        ABBA_groups = [ABBA_coords[i:i + 4] for i in range(0, len(ABBA_coords), 4)]
        value_groups = [value[i:i + 4] for i in range(0, len(ABBA_coords), 4)]

        for group in range(len(ABBA_groups)):
            if len(ABBA_groups[group]) == 2:
                # Warn user that if there are any groups of only 2
                # files, assuming they are in order of A and B
                msgs.info('Assuming these two frames are A and B frame:'
                          + msgs.newline() + file_groups[group][0]
                          + msgs.newline() + file_groups[group][1])
            elif len(ABBA_groups[group]) == 4:
                # Check that frames 1, 4 of an ABBA sequence are at the
                # same nod position (A) based on their RA, DEC
                AA_sep = ABBA_coords[0].separation(ABBA_coords[-1]).arcsec
                BB_sep = ABBA_coords[1].separation(ABBA_coords[2]).arcsec
                if AA_sep > max_nod_sep or BB_sep > max_nod_sep:
                    if AA_sep > max_nod_sep:
                        msgs.warn('Separation between 1st and 4th frame in presumed ABBA sequence '
                                  'have a large separation ({0}).'.format(AA_sep))
                    if BB_sep > max_nod_sep:
                        msgs.warn('Separation between 2nd and 3rd frame in presumed ABBA sequence '
                                  'have a large separation ({0}).'.format(BB_sep))
                    msgs.warn('Check ABBA identification for target {0} group {1}:'.format(
                                target, group) + msgs.newline() + 'A:' + file_groups[group][0]
                              + msgs.newline() + 'B:' + file_groups[group][1]
                              + msgs.newline() + 'B:' + file_groups[group][2]
                              + msgs.newline() + 'A:' + file_groups[group][3])
            else:
                msgs.error('BUG: This should never be reached.')

            # Flip group from ABBA to BABA, or AB to BA
            AB_idx_flip = np.copy(value_groups[group])
            AB_idx_flip[::2], AB_idx_flip[1::2] \
                    = value_groups[group][1::2], value_groups[group][::2]

            # Associate each file in the group with its AB pair
            for i,j in enumerate(value_groups[group]):
                AB_frame[j] = file_list[AB_idx_flip[i]]

    return AB_frame    


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
            import pdb; pdb.set_trace()
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




