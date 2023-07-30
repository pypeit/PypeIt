"""
parse module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import inspect
import numpy as np

# Logging
from pypeit import msgs


def get_current_name():
    """ Return the name of the function that called this function

    Returns
    -------
    name: str
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    return "'" + " ".join(ll) + "'"


def get_nmbr_name(anmbr=None, bnmbr=None, cnmbr=None):
    """ Return the name of the function that called this function,
    and append a two digit number to the first (when anmbr!=None),
    the second (when bnmbr!=None), or third (when cnmbr!=None) element
    of the function name.
    """
    cspl = inspect.currentframe().f_back.f_code.co_name.split('_')
    if anmbr is not None:
        cspl[0] += "{0:02d}".format(anmbr)
    if bnmbr is not None:
        cspl[1] += "{0:02d}".format(bnmbr)
    if cnmbr is not None:
        cspl[2] += "{0:02d}".format(cnmbr)
    return "_".join(cspl)


def load_sections(string, fmt_iraf=True):
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
    fmt_iraf : bool
      Is the variable string in IRAF format (True) or
      python format (False)

    Returns
    -------
    sections : list (or None)
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


def check_deprecated(v, deprecated, upper=False):
    """ Check if a keyword argument is deprecated.

    Parameters
    ----------
    v : str
      value of a keyword argument
    deprecated : list
      list of deprecated values that v might be
    upper : bool (optional)
      If True, the allowed list is expected to contain only
      uppercase strings. If False, the allowed list is expected
      to contain only lowercase strings.

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if upper:
        v = v.upper()
    else:
        v = v.lower()
    if v in deprecated:
        msgs.error("The argument of {0:s} is deprecated.".format(func_name) + msgs.newline() +
                   "Please choose one of the following:" + msgs.newline() +
                   ", ".join(deprecated))
    return


def key_allowed(v, allowed, upper=False):
    """ Check that a keyword argument is in an allowed list of parameters.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take
    upper : bool (optional)
      If True, the allowed list is expected to contain only
      uppercase strings. If False, the allowed list is expected
      to contain only lowercase strings.

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if upper:
        v = v.upper()
    else:
        v = v.lower()
    if v not in allowed:
        msgs.error("The argument of {0:s} must be one of".format(func_name) + msgs.newline() +
                   ", ".join(allowed))
    if v.lower() == "none":
        v = None
    return v


def key_allowed_filename(v, allowed):
    """ Check that a keyword argument is in an allowed list of parameters.
    If not, assume that it is a filename.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    vt = v.lower()
    if vt not in allowed:
        msgs.warn("Assuming the following is the name of a file:" + msgs.newline() + v)
    else:
        v = vt
    return v


def key_bool(v):
    """ Check that a keyword argument is a boolean variable.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : str
      A string used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "true":
        v = True
    elif v.lower() == "false":
        v = False
    else:
        msgs.error("The argument of {0:s} can only be 'True' or 'False'".format(func_name))
    return v


def key_check(v):
    """ Check that a keyword argument satisfies the form required of a
    keyword argument that checks frame types.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : str, list
      A value used by the settings dictionary
    """
    text = v.strip().replace('_', ' ')
    if ',' in text and text[0:2] != '%,':
        # There are multiple possibilities - split the text
        v = text.split(',')
    else:
        v = text
    return v


def key_float(v):
    """ Check that a keyword argument is a float.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : float
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = float(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type float".format(func_name))
    return v


def key_int(v):
    """ Check that a keyword argument is an int.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : int
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    try:
        v = int(v)
    except ValueError:
        msgs.error("The argument of {0:s} must be of type int".format(func_name))
    return v


def key_keyword(v, force_format=True):
    """ Check that a keyword argument satisfies the form required
    for specifying a header keyword.

    Parameters
    ----------
    v : str
      value of a keyword argument
    force_format : bool
      If True, v can only have the format of a header keyword.
      If False, v can take the form of a header keyword, or a
      user-specified value. In the latter case, the boolean
      variable 'valid' is returned so the parent function
      knows if the supplied value of v should be dealt with as
      a manual entry or as a header keyword.

    Returns
    -------
    v : str
      A value used by the settings dictionary
    valid : bool
      Is the input a valid header keyword format
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "none":
        v = None
    else:
        valid = is_keyword(v)
        if not valid and force_format:
            msgs.error("The argument of {0:s} must be of the form:".format(func_name) + msgs.newline() +
                       "##.NAME" + msgs.newline() +
                       "where ## is the fits extension (see command: fits headext##)," + msgs.newline() +
                       "and NAME is the header keyword name")
        elif valid and force_format:
            return v
        else:
            return v, valid
    return v


def key_list(strlist):
    """ Check that a keyword argument is a list. Set the
    appropriate type of the list based on the supplied values.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
    # Check if the input array is a null list
    if strlist == "[]" or strlist == "()":
        return []
    # Remove outer brackets and split by commas
    temp = strlist.lstrip('([').rstrip(')]').split(',')
    addarr = []
    # Find the type of the array elements
    for i in temp:
        if i.lower() == 'none':
            # None type
            addarr += [None]
        elif i.lower() == 'true' or i.lower() == 'false':
            # bool type
            addarr += [i.lower() in ['true']]
        elif ',' in i:
            # a list
            addarr += i.lstrip('([').rstrip('])').split(',')
            msgs.bug("nested lists could cause trouble if elements are not strings!")
        elif '.' in i:
            try:
                # Might be a float
                addarr += [float(i)]
            except ValueError:
                # Must be a string
                addarr += [i]
        else:
            try:
                # Could be an integer
                addarr += [int(i)]
            except ValueError:
                # Must be a string
                addarr += [i]
    return addarr


def key_list_allowed(v, allowed):
    """ Check that a keyword argument is a list. Set the
    appropriate type of the list based on the supplied values.
    Then, check that each value in the list is also in the
    supplied 'allowed' list.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    v = key_list(v)
    for ll in v:
        if ll not in allowed:
            msgs.error("The allowed list does not include: {0:s}".format(ll) + msgs.newline() +
                       "Please choose one of the following:" + msgs.newline() +
                       ", ".join(allowed) + msgs.newline() +
                       "for the argument of {0:s}".format(func_name))
    return v


def key_none_allowed(v, allowed):
    """ Check if a keyword argument is set to None. If not,
    check that its value is in the 'allowed' list.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
    """
    ll = inspect.currentframe().f_back.f_code.co_name.split('_')
    func_name = "'" + " ".join(ll) + "'"
    if v.lower() == "none":
        v = None
    elif v.lower() in allowed:
        for i in allowed:
            if v.lower() == i:
                v = i
                break
    else:
        msgs.error("The argument of {0:s} must be one of".format(func_name) + msgs.newline() +
                   ", ".join(allowed))
    return v


def key_none_allowed_filename(v, allowed):
    """ Check if a keyword argument is set to None. If not,
    check that its value is in the 'allowed' list. Finally,
    assume that the supplied value is the name of a filename.

    Parameters
    ----------
    v : str
      value of a keyword argument
    allowed : list
      list of allowed values that v can take

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    elif v.lower() in allowed:
        for i in allowed:
            if v.lower() == i:
                v = i
                break
    else:
        msgs.info("Assuming the following is the name of a file:" + msgs.newline() + v)
    return v

def key_none_int(v):
    """ Check if a keyword argument is set to None. If not,
    assume the supplied value is an int.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    else:
        v = key_int(v)
    return v

def key_none_list(v):
    """ Check if a keyword argument is set to None. If not,
    assume the supplied value is a list.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : list
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    else:
        if "," in v:
            v = v.strip("()[]").split(",")
        else:
            v = [v]
    return v


def key_none(v):
    """ Check if a keyword argument is set to None.

    Parameters
    ----------
    v : str
      value of a keyword argument

    Returns
    -------
    v : None, str
      A value used by the settings dictionary
    """
    if v.lower() == "none":
        v = None
    return v


def key_min_val(v, vmin):
    """ Check that the value exceeds a minimum
    Returns
    -------
    bool

    """
    if v < vmin:
        msgs.error("The argument of {0:s} must be >= -1".format(get_current_name()))
    else:
        return True

def combine_methods():
    """ The methods that can be used to combine a set of frames into a master frame
    """
    methods = ['mean', 'median', 'weightmean']
    return methods


def combine_replaces():
    """ The options that can be used to replace rejected pixels when combining a set of frames
    """
    methods = ['min', 'max', 'mean', 'median', 'weightmean', 'maxnonsat']
    return methods


def combine_satpixs():
    """ The options that can be used to replace saturated pixels when combining a set of frames
    """
    methods = ['reject', 'force', 'nothing']
    return methods


def is_keyword(v):
    """ Check if a value is of the format required to be a call to a header keyword

    Parameters
    ----------
    v : str
      Value to be tested

    Returns
    -------
    valid : bool
      True if 'v' has the correct format to be a header keyword, False otherwise.
    """
    valid = True
    if ("," in v) or ("." not in v):
        # Either an array or doesn't have the header keyword format (i.e. a fullstop)
        return False
    # Test if the first element is an integer
    vspl = v.split(".")
    try:
        int(vspl[0])
    except ValueError:
        valid = False
    # Test if there are two parts to the expression
    if len(vspl) != 2:
        return False
    # Test if the second element is a string
    try:
        if valid is True:
            int(vspl[1])
            # Input value must be a floating point number
            valid = False
    except ValueError:
        # Input value must be a string
        valid = True
    return valid


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
        int,int: binspectral, binspatial

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
        :obj:`tuple`:  Two integer arrays with the list of 1-indexed detector
        numbers and spatial pixels coordinates for each slit.  The shape of each
        array is ``(nslits,)``, where ``nslits`` is the number of
        ``slitspatnum`` entries parsed (1 if only a single string is provided).
    """
    dets = []
    spat_ids = []
    if isinstance(slitspatnum,list):
        slitspatnum = ",".join(slitspatnum)
    for item in slitspatnum.split(','):
        spt = item.split(':')
        dets.append(int(spt[0]))
        spat_ids.append(int(spt[1]))
    # Return
    return np.array(dets).astype(int), np.array(spat_ids).astype(int)


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


def str2list(inp, length):
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
        length (:obj:`int`):
            Maximum length of the list, which is needed to allow for
            open slices (e.g., '8:').
    
    Returns:
        list: List of parsed integers.
    """
    if inp == 'None':
        return None

    gi = np.arange(length)
    if inp == 'all':
        # Flag 'em all!
        return gi.tolist()

    # Parse the input string into a list of integers
    grp = np.concatenate([[int(g)] if g.find(':') == -1 else (gi[sec2slice(g)]).tolist()
                                for g in inp.split(',')])

    # Return the list, and ensure the integers are unique
    return np.unique(grp).tolist()

