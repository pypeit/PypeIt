"""
Routines for matching frames to certain types or each other.

.. include:: ../include/links.rst
"""
# TODO -- Move this out of core?

from collections import OrderedDict

import numpy as np

from pypeit import log
from pypeit import PypeItError
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
            ('slitless_pixflat', 'Flat-field exposure without slitmask used for pixel-to-pixel response'),
                  ('scattlight', 'Frame (ideally with lots of counts) used to determine the scattered light model'),
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
        message = f'{frametype} is not a valid PypeIt frame type.'
        if not quiet and not raise_error:
            log.warning(message)
        elif raise_error:
            raise PypeItError(message)
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
