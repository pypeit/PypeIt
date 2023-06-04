"""
Routines for matching frames to certain types or each other.

.. include:: ../include/links.rst
"""
# TODO -- Move this out of core?

from collections import OrderedDict

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
        # TODO JFH: We need a background image type
        frame_types = OrderedDict([
                       ('align', 'Trace constant spatial positions along the slit'),
                         ('arc', 'Arc lamp observation used for wavelength calibration'),
                        ('bias', 'Bias readout for detector bias subtraction'),
                        ('dark', 'Shuttered exposure to measure dark current'),
                     ('pinhole', 'Pinhole observation used for tracing slit centers'),
                   ('pixelflat', 'Flat-field exposure used for pixel-to-pixel response'),
                   ('illumflat', 'Flat-field exposure used for illumination flat'),
                   ('lampoffflats', 'Flat-field exposure with lamps off used to remove '
                                   'persistence from lamp on flat exposures and/or thermal emission '
                                   'from the telescope and dome'),
                     ('science', 'On-sky observation of a primary target'),
                    ('standard', 'On-sky observation of a flux calibrator'),
                       ('trace', 'High-count exposure used to trace slit positions'),
                        ('tilt', 'Exposure used to trace the tilt in the wavelength solution'),
                         ('sky', 'On-sky observation of the sky used for background subtraction'),
                                  ])
        super(FrameTypeBitMask, self).__init__(list(frame_types.keys()),
                                               descr=list(frame_types.values()))

    def type_names(self, type_bits, join=True):
        """
        Use the type bits to get the type names for each frame.

        .. todo::
            - This should probably be a general function in
              :class:`pypeit.bitmask.BitMask`
    
        Args:
            type_bits (int, list, `numpy.ndarray`_):
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


def valid_frametype(frametype, quiet=False, raise_error=False):
    """
    Confirm the provided frame type is known to ``PypeIt``.

    Args:
        frametype (:obj:`str`):
            The frame type name.
        quiet (:obj:`bool`, optional):
            Suppress output
        raise_error (:obj:`bool`, optional):
            Instead of issuing a warning, raise an exception.

    Returns:
        :obj:`bool`: Flag that the frametype name is valid.
    """
    good_frametype = frametype in FrameTypeBitMask().keys()
    if not good_frametype:
        _f = None
        if not quiet and not raise_error:
            _f = msgs.warn
        elif raise_error:
            _f = msgs.error
        if _f is not None:
            _f(f'{frametype} is not a valid PypeIt frame type.')
    return good_frametype
    

def check_frame_exptime(exptime, exprng):
    """
    Check that the exposure time is within the provided range.
        
    Args:
        exptime (`numpy.ndarray`_):
            Exposure times to check; allowed to be None.
        exprng (array-like):
            An array with the minimum and maximum exposure.  The limits
            are *exclusive* and a limit of None means there is no limit.
        
    Returns:
        `numpy.ndarray`_: A boolean array that is True for all times within
        the provided range. The value is False for any exposure time that is
        None or outside the provided range.
        
    Raises:
        ValueError:
            Raised if the length of `exprng` is not 2.
    """
    # Instantiate with all true
    indx = exptime != None
    if exprng is None:
        # No range specified
        return indx
    if len(exprng) != 2:
        # Range not correctly input
        raise ValueError('exprng must have two elements.')
    if exprng[0] is not None:
        indx[indx] &= (exptime[indx] > exprng[0])
    if exprng[1] is not None:
        indx[indx] &= (exptime[indx] <= exprng[1])
    return indx


# TODO: May want to keep this in case we ever try to bring it back....
#def group_AB_frames(file_list, targets, coords, max_nod_sep=2):
#    """
#    Group files into a ABBA or AB sequences.
#
#    Args:
#        file_list (:obj:`list`):
#            A list of file names.
#        targets (:obj:`dict`):
#            A dictionary that matches each file to a unique target name.
#            The target name can be one of the files in the file list.
#        coords (:class:`astropy.coordinates.SkyCoord`):
#            The coordinates of all the exposures.  Number of coordinates
#            should match the number of files.
#        max_nod_sep (:obj:`int`, optional):
#            The maximum separation (arcsec) between the 1st and 4th
#            frame sky coordinates in the ABBA sequence that is allowed
#            when identifying the sequence.  Note that the default (2
#            arcsec) is arbitrary.
#    
#    Returns:
#        list:
#            A list that matches the length of the input list of files.
#            Each file in an AB or ABBA sequence is identified with it's
#            pair in the sequence.
#    """
#
#    AB_frame = [''] * len(file_list)
#
#    for key, value in targets.items():
#        files = file_list[value]
#
#        # Check here that there are more than 1 files and that the
#        # number of files is even
#        if len(files) == 1:
#            msgs.warn('Cannot perform ABBA reduction on targets with 1 file')
#        elif len(files) % 2 != 0:
#            msgs.warn('Expected an even number of files associated with target ' + key)
#
#        # TODO: Check for increasing time? Files are read in numerical
#        # sequential order -- should be in order of increasing time
#        # anyway..
#
#        # Assume that the files are initially in ABBA order and proceed
#        ABBA_coords = coords[value]
#
#        # Break files into ABBA groups (includes remainder if there are only 2 files)
#        file_groups = [files[i:i+4] for i in range(0,len(files),4)]
#        ABBA_groups = [ABBA_coords[i:i + 4] for i in range(0, len(ABBA_coords), 4)]
#        value_groups = [value[i:i + 4] for i in range(0, len(ABBA_coords), 4)]
#
#        for group in range(len(ABBA_groups)):
#            if len(ABBA_groups[group]) == 2:
#                # Warn user that if there are any groups of only 2
#                # files, assuming they are in order of A and B
#                msgs.info('Assuming these two frames are A and B frame:'
#                          + msgs.newline() + file_groups[group][0]
#                          + msgs.newline() + file_groups[group][1])
#            elif len(ABBA_groups[group]) == 4:
#                # Check that frames 1, 4 of an ABBA sequence are at the
#                # same nod position (A) based on their RA, DEC
#                AA_sep = ABBA_coords[0].separation(ABBA_coords[-1]).arcsec
#                BB_sep = ABBA_coords[1].separation(ABBA_coords[2]).arcsec
#                if AA_sep > max_nod_sep or BB_sep > max_nod_sep:
#                    if AA_sep > max_nod_sep:
#                        msgs.warn('Separation between 1st and 4th frame in presumed ABBA sequence '
#                                  'have a large separation ({0}).'.format(AA_sep))
#                    if BB_sep > max_nod_sep:
#                        msgs.warn('Separation between 2nd and 3rd frame in presumed ABBA sequence '
#                                  'have a large separation ({0}).'.format(BB_sep))
#                    msgs.warn('Check ABBA identification for target {0} group {1}:'.format(
#                                target, group) + msgs.newline() + 'A:' + file_groups[group][0]
#                              + msgs.newline() + 'B:' + file_groups[group][1]
#                              + msgs.newline() + 'B:' + file_groups[group][2]
#                              + msgs.newline() + 'A:' + file_groups[group][3])
#            else:
#                msgs.error('BUG: This should never be reached.')
#
#            # Flip group from ABBA to BABA, or AB to BA
#            AB_idx_flip = np.copy(value_groups[group])
#            AB_idx_flip[::2], AB_idx_flip[1::2] \
#                    = value_groups[group][1::2], value_groups[group][::2]
#
#            # Associate each file in the group with its AB pair
#            for i,j in enumerate(value_groups[group]):
#                AB_frame[j] = file_list[AB_idx_flip[i]]
#
#    return AB_frame    


