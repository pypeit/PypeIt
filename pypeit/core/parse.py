"""
parse module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect

from IPython import embed

import numpy as np

# Logging
from pypeit import msgs


def load_sections(string, fmt_iraf=True):
    """
    From the input string, return the coordinate sections.
    In IRAF format (1 index) or Python

    Parameters
    ----------
    string : str
      character string of the form [x1:x2,y1:y2]
      x1 = left pixel
      x2 = right pixel
      y1 = bottom pixel
      y2 = top pixel
    fmt_iraf : bool
      Is the variable string in IRAF format (True) or
      python format (False)

    Returns
    -------
    sections : list or None
      the detector sections
    """
    xyrng = string.strip('[]()').split(',')
    if xyrng[0] == ":":
        xyarrx = [0, 0]
    else:
        xyarrx = xyrng[0].split(':')
        # If a lower/upper limit on the array slicing is not given (e.g. [:100] has no lower index specified),
        # set the lower/upper limit to be the first/last index.
        if len(xyarrx[0]) == 0: xyarrx[0] = 0
        if len(xyarrx[1]) == 0: xyarrx[1] = -1
    if xyrng[1] == ":":
        xyarry = [0, 0]
    else:
        xyarry = xyrng[1].split(':')
        # If a lower/upper limit on the array slicing is not given (e.g. [5:] has no upper index specified),
        # set the lower/upper limit to be the first/last index.
        if len(xyarry[0]) == 0: xyarry[0] = 0
        if len(xyarry[1]) == 0: xyarry[1] = -1
    if fmt_iraf:
        xmin = max(0, int(xyarry[0])-1)
        xmax = int(xyarry[1])
        ymin = max(0, int(xyarrx[0])-1)
        ymax = int(xyarrx[1])
    else:
        xmin = max(0, int(xyarrx[0]))
        xmax = int(xyarrx[1])
        ymin = max(0, int(xyarry[0]))
        ymax = int(xyarry[1])
    return [[xmin, xmax], [ymin, ymax]]


def get_dnum(det, caps=False, prefix=True):
    """ Convert a detector index into a string used by the settings dictionary
    or other bits of code.  Best to keep at two digits

    Parameters
    ----------
    det : int
      Detector index
    caps : bool, optional
      Return all caps?
    prefix : bool, optional
      Include the prefix?

    Returns
    -------
    dnum : str
      A string used by the settings dictionary
    """
    dnum = '{0:02d}'.format(det)
    if prefix:
        if caps:
            dnum = 'DET'+dnum
        else:
            dnum = 'det'+dnum
    # Return
    return dnum


def binning2string(binspectral, binspatial):
    """
    Convert the binning from integers to a string following the PypeIt
    convention order, spectral then spatial.

    Args:
        binspectral (:obj:`int`):
            Number of on-detector pixels binned in the spectral
            direction (along the first axis in the PypeIt convention).
        binspatial (int):
            Number of on-detector pixels binned in the spatial direction
            (along the second axis in the PypeIt convention).

    Returns:
        str: Comma-separated binning along the spectral and spatial
        directions; e.g., '2,1'
    """
    return '{0},{1}'.format(binspectral, binspatial)


def parse_binning(binning):
    """
    Parse input binning into binspectral, binspatial

    Note that for some instruments, the meaning will be swapped if
    parsed directly from the Header.  The developer needs to react accordingly..

    Args:
        binning (str, ndarray or tuple):

    Returns:
        tuple: binspectral, binspatial as integers

    """
    # comma separated format
    if isinstance(binning, str):
        if ',' in binning:
            binspectral, binspatial = [int(item) for item in binning.split(',')]  # Keck standard, I think
        elif 'x' in binning:
            binspectral, binspatial = [int(item) for item in binning.split('x')]  # LRIS
        elif binning == 'None':
            msgs.warn("Assuming unbinned, i.e.  1x1")
            binspectral, binspatial = 1,1
        else:
            binspectral, binspatial = [int(item) for item in binning.strip().split(' ')]  # Gemini
    elif isinstance(binning, tuple):
        binspectral, binspatial = binning
    elif isinstance(binning, np.ndarray):
        binspectral, binspatial = binning
    else:
        msgs.error("Unable to parse input binning: {}".format(binning))
    # Return
    return binspectral, binspatial

# TODO: Allow this to differentiate between detectors and mosaics.  Input syntax
# likely to become, e.g., DET01:175,DET01:205.
def parse_slitspatnum(slitspatnum):
    """
    Parse the ``slitspatnum`` into a list of detectors and SPAT_IDs.

    Args:
        slitspatnum (:obj:`str`, :obj:`list`):
            A single string or list of strings to parse.

    Returns:
        :obj:`tuple`:  Two arrays with the list of 1-indexed detectors (str)
        and spatial pixels coordinates for each slit.  The shape of each
        array is ``(nslits,)``, where ``nslits`` is the number of
        ``slitspatnum`` entries parsed (1 if only a single string is provided).
    """
    dets = []
    spat_ids = []
    if isinstance(slitspatnum,list):
        slitspatnum = ",".join(slitspatnum)
    for item in slitspatnum.split(','):
        spt = item.split(':')
        dets.append(spt[0])
        spat_ids.append(int(spt[1]))
    # Return
    return np.array(dets).astype(str), np.array(spat_ids).astype(int)


def sec2slice(subarray, one_indexed=False, include_end=False, require_dim=None, binning=None):
    """
    Convert a string representation of an array subsection (slice) into
    a list of slice objects.

    Args:
        subarray (str):
            The string to convert.  Should have the form of normal slice
            operation, 'start:stop:step'.  The parser ignores whether or
            not the string has the brackets '[]', but the string must
            contain the appropriate ':' and ',' characters.
        one_indexed (:obj:`bool`, optional):
            The string should be interpreted as 1-indexed.  Default
            is to assume python indexing.
        include_end (:obj:`bool`, optional):
            **If** the end is defined, adjust the slice such that
            the last element is included.  Default is to exclude the
            last element as with normal python slicing.
        require_dim (:obj:`int`, optional):
            Test if the string indicates the slice along the proper
            number of dimensions.
        binning (:obj:`str`, optional):
            Assume the slice is for an unbinned array and adjust the
            returned slice for this binning in each dimension.  If two
            dimensional, the format of this string must be, e.g., `1,2`
            for unbinned rows and a factor of 2 binning along columns.

    Returns:
        tuple: A tuple of slice objects, one per dimension of the
        prospective array.

    Raises:
        TypeError:
            Raised if the input `subarray` is not a string.
        ValueError:
            Raised if the string does not match the required
            dimensionality or if the string does not look like a
            slice.
    """
    # Check it's a string
    if not isinstance(subarray, str):
        raise TypeError('Can only parse string-based subarray sections.')
    # Remove brackets if they're included
    sections = subarray.strip('[]').split(',')
    # Check the dimensionality
    ndim = len(sections)
    _binning = [1]*ndim if binning is None else np.array(binning.split(',')).astype(int)
    if len(_binning) != ndim:
        raise ValueError('Incorrect binning dimensions (found {0}, expected {1}).'.format(
                            len(_binning), ndim))
    if require_dim is not None and ndim != require_dim:
        raise ValueError('Number of slices ({0}) in {1} does not match '.format(ndim, subarray) + 
                         'required dimensions ({0}).'.format(require_dim))
    # Convert the slice of each dimension from a string to a slice
    # object
    slices = []
    for s,b in zip(sections,_binning):
        flipped = False
        # Must be able to find the colon
        if ':' not in s:
            raise ValueError('Unrecognized slice string: {0}'.format(s))
        # Initial conversion
        _s = [ None if x == '' else int(x) for x in s.split(':') ]
        if len(_s) > 3:
            raise ValueError('String as too many sections.  Must have format \'start:stop:step\'.')
        if len(_s) < 3:
            # Include step
            _s += [ None ]
        # Must check order first so "include_last" and "one_indexed" are correctly applied
        # Check that the first two elements of the slice are ordered correctly
        if _s[0] is not None and _s[1] is not None:
            if _s[0] > _s[1]:
                flipped = True
                _s = [_s[1], _s[0], _s[2]]
        if one_indexed:
            # Decrement to convert from 1- to 0-indexing
            _s = [ None if x is None else x-1 for x in _s ]
        if include_end and _s[1] is not None:
            # Increment to include last 
            _s[1] += 1
        _s = [ None if ss is None else ss//b for ss in _s ]
        if flipped:
            if _s[0] == 0:
                _s = [_s[1]-1, None, -1]
            else:
                _s = [_s[1]-1, _s[0]-1, -1]
        # Append the new slice
        slices += [slice(*_s)]

    return tuple(slices)


def str2list(inp, length=None):
    """
    Expand a string with a comma-separated set of integers and slices
    into a list of the relevant integers.

    Setting a maximum length of the list to 10, examples of the allowed
    syntax and result are:

        - 'all': [0,1,2,3,4,5,6,7,8,9]
        - ':4': [0,1,2,3]
        - '3:5,8:': [3,4,8,9]
        - '3,1:5,6': [1,2,3,4,6]

    Note the function removes any non-unique integers (see the last
    example).

    Args:
        inp (:obj:`str`):
            String with a comma-separated set of integers and slices;
            can also be 'all'.
        length (:obj:`int`, optional):
            Maximum length of the list, which is needed to allow for open slices
            (e.g., '8:').  If None, open slices and setting ``inp='all'`` will
            raise a ValueError.
    
    Returns:
        list: List of parsed integers.
    """
    if inp == 'None':
        return None

    if inp == 'all':
        if length is None:
            raise ValueError('To use inp=all in str2list, must provide length.')
        return np.arange(length).tolist()

    if length is not None:
        gi = np.arange(length)

    groups = inp.split(',')
    indices = []
    for i, grp in enumerate(groups):
        if ':' in grp:
            if (grp[-1] == ':' or grp[0] == ':'):
                if length is None:
                    raise ValueError('To use open ended slices in str2lit, must provide length.')
                indices += gi[sec2slice(grp)].tolist()
            else:
                start, end = map(lambda x: int(x), grp.split(':'))
                indices += [i for i in range(start, end)]
        else:
            indices += [int(grp)]
    return np.unique(indices).tolist()

