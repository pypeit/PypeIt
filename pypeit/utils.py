"""
General utility functions.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import inspect
import pickle
import pathlib
import itertools
import glob
import collections.abc

from IPython import embed

import numpy as np
from numpy.lib.stride_tricks import as_strided

import scipy.ndimage

import matplotlib
import matplotlib.pyplot as plt

from astropy import units
from astropy import stats

from pypeit import msgs
from pypeit.move_median import move_median


def get_time_string(codetime):
    """
    Utility function that takes the codetime and
    converts this to a human readable String.

    Args:
        codetime (`float`):
            Code execution time in seconds (usually the difference of two time.time() calls)

    Returns:
        `str`: A string indicating the total execution time
    """
    if codetime < 60.0:
        retstr = 'Execution time: {0:.2f}s'.format(codetime)
    elif codetime / 60.0 < 60.0:
        mns = int(codetime / 60.0)
        scs = codetime - 60.0 * mns
        retstr = 'Execution time: {0:d}m {1:.2f}s'.format(mns, scs)
    else:
        hrs = int(codetime / 3600.0)
        mns = int(60.0 * (codetime / 3600.0 - hrs))
        scs = codetime - 60.0 * mns - 3600.0 * hrs
        retstr = 'Execution time: {0:d}h {1:d}m {2:.2f}s'.format(hrs, mns, scs)
    return retstr


def all_subclasses(cls):
    """
    Collect all the subclasses of the provided class.

    The search follows the inheritance to the highest-level class.  Intermediate
    base classes are included in the returned set, but not the base class itself.

    Thanks to:
    https://stackoverflow.com/questions/3862310/how-to-find-all-the-subclasses-of-a-class-given-its-name

    Args:
        cls (object):
            The base class

    Returns:
        :obj:`set`: The unique set of derived classes, including any
        intermediate base classes in the inheritance thread.
    """
    return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in all_subclasses(c)])


def embed_header():
    """
    Nominal header for an execution of `IPython.embed`_.

    Example:

        To include the returned string::

            from IPython import embed
            from pypeit.utils import embed_header

            embed(header=embed_header())

    Returns:
        :obj:`str`: String with the line in the calling module, the
        name of the calling function, and the name of the calling
        file.
    """
    info = inspect.getframeinfo(inspect.stack()[1][0])
    return '{0} {1} {2}'.format(info.lineno, info.function, os.path.split(info.filename)[1])


# Pulled from `pypeit.par.ParSet`. Maybe move these to
# doc/scripts/util.py?
def to_string(data, use_repr=True, verbatim=False):
    """
    Convert a single datum into a string

    Simply return strings, recursively convert the elements of any
    objects with a :attr:`__len__` attribute, and use the object's
    own :attr:`__repr__` attribute for all other objects.

    Args:
        data (object):
            The object to stringify.
        use_repr (:obj:`bool`, optional):
            Use the objects :attr:`__repr__` method; otherwise, use a
            direct string conversion.
        verbatim (:obj:`bool`, optional):
            Use quotes around the provided string to indicate that
            the string should be represented in a verbatim (fixed
            width) font.
        
    Returns:
        :obj:`str`: A string representation of the provided ``data``.
    """
    if isinstance(data, str):
        return data if not verbatim else '``' + data + '``'
    if hasattr(data, '__len__'):
        return '[]' if isinstance(data, list) and len(data) == 0 \
                    else ', '.join([to_string(d, use_repr=use_repr, verbatim=verbatim) 
                                        for d in data ])
    return data.__repr__() if use_repr else str(data)


def string_table(tbl, delimeter='print', has_header=True):
    """
    Provided the array of data, format it with equally spaced columns
    and add a header (first row) and contents delimeter.

    Args:
        tbl (`numpy.ndarray`_):
            Array of string representations of the data to print.
        delimeter (:obj:`str`, optional):
            If the first row in the table containts the column headers (see
            ``has_header``), this sets the delimeter between first table row and
            the column data. Use ``'print'`` for a simple line of hyphens,
            anything else results in an ``rst`` style table formatting.
        has_header (:obj:`bool`, optional):
            The first row in ``tbl`` contains the column headers.

    Returns:
        :obj:`str`: Single long string with the data table.
    """
    nrows, ncols = tbl.shape
    col_width = [np.amax([len(dij) for dij in dj]) for dj in tbl.T]

    _nrows = nrows
    start = 1
    if delimeter != 'print':
        _nrows += 2
        start += 1
    if has_header:
        _nrows += 1
        start += 1

    row_string = ['']*_nrows

    for i in range(start,nrows+start-1):
        row_string[i] = '  '.join([tbl[1+i-start,j].ljust(col_width[j]) for j in range(ncols)])
    if delimeter == 'print':
        # Heading row
        row_string[0] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
        # Delimiter
        if has_header:
            row_string[1] = '-'*len(row_string[0])
        return '\n'.join(row_string)+'\n'

    # For an rst table
    row_string[0] = '  '.join([ '='*col_width[j] for j in range(ncols)])
    row_string[1] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
    if has_header:
        row_string[2] = row_string[0]
    row_string[-1] = row_string[0]
    return '\n'.join(row_string)+'\n'


def spec_atleast_2d(wave, flux, ivar, gpm, copy=False):
    """
    Force spectral arrays to be 2D.

    Input and output spectra are ordered along columns; i.e., the flux vector
    for the first spectrum is in ``flux[:,0]``.
    
    Args:
        wave (`numpy.ndarray`_):
            Wavelength array. Must be 1D if the other arrays are 1D. If 1D
            and the other arrays are 2D, the wavelength vector is assumed to
            be the same for all spectra.
        flux (`numpy.ndarray`_):
            Flux array.  Can be 1D or 2D.
        ivar (`numpy.ndarray`_):
            Inverse variance array for the flux.  Shape must match ``flux``.
        gpm (`numpy.ndarray`_):
            Good pixel mask (i.e., True=Good). Shape must match ``flux``.
        copy (:obj:`bool`, optional):
            If the flux, inverse variance, and gpm arrays are already 2D on
            input, the function just returns the input arrays. This flag
            forces the returned arrays to be copies instead.

    Returns:
        :obj:`tuple`: Returns 6 objects. The first four are the reshaped
        wavelength, flux, inverse variance, and gpm arrays. The next two
        give the length of each spectrum and the total number of spectra;
        i.e., the last two elements are identical to the shape of the
        returned flux array.

    Raises:
        PypeItError:
            Raised if the shape of the input objects are not appropriately
            matched.
    """
    # Check the input
    if wave.shape[0] != flux.shape[0] or ivar.shape != flux.shape or gpm.shape != flux.shape \
            or wave.ndim == 2 and wave.shape != flux.shape:
        msgs.error('Input spectral arrays have mismatching shapes.')

    if flux.ndim == 1:
        # Input flux is 1D
        # NOTE: These reshape calls return copies of the arrays
        return wave.reshape(-1, 1), flux.reshape(-1, 1), ivar.reshape(-1, 1), \
                    gpm.reshape(-1, 1), flux.size, 1

    # Input is 2D
    nspec, norders = flux.shape
    _wave = np.tile(wave, (norders, 1)).T if wave.ndim == 1 else (wave.copy() if copy else wave)
    _flux = flux.copy() if copy else flux
    _ivar = ivar.copy() if copy else ivar
    _gpm = gpm.copy() if copy else gpm
    return _wave, _flux, _ivar, _gpm, nspec, norders


def nan_mad_std(data, axis=None, func=None):
    """

    Wrapper for astropy.stats.mad_std which ignores nans, so as to
    prevent bugs when using sigma_clipped_stats with the axis keyword
    and stdfunc=astropy.stats.mad_std

    Args:
        data (array-like):
            Data array or object that can be converted to an array.
        axis (int, sequence of int, None, optional):
            Axis along which the robust standard deviations are
            computed.  The default (`None`) is to compute the robust
            standard deviation of the flattened array.

    Returns:
        float, `numpy.ndarray`: The robust standard deviation of the
        input data.  If ``axis`` is `None` then a scalar will be
        returned, otherwise a `~numpy.ndarray` will be returned.
    """
    return stats.mad_std(data, axis=axis, func=func, ignore_nan=True)


def growth_lim(a, lim, fac=1.0, midpoint=None, default=[0., 1.]):
    """
    Calculate bounding limits for an array based on its growth.

    Args:
        a (array-like):
            Array for which to determine limits.
        lim (:obj:`float`):
            Percentage of the array values to cover. Set to 1 if
            provided value is greater than 1.
        fac (:obj:`float`, optional):
            Factor to increase the range based on the growth limits.
            Default is no increase.
        midpoint (:obj:`float`, optional):
            Force the midpoint of the range to be centered on this
            value. Default is the sample median.
        default (:obj:`list`, optional):
            Default limits to return if `a` has no data.

    Returns:
        :obj:`list`: Lower and upper boundaries for the data in `a`.
    """
    # Get the values to plot
    _a = a.compressed() if isinstance(a, np.ma.MaskedArray) else np.asarray(a).ravel()
    if len(_a) == 0:
        # No data so return the default range
        return default

    # Set the starting and ending values based on a fraction of the
    # growth
    _lim = 1.0 if lim > 1.0 else lim
    start, end = (len(_a)*(1.0+_lim*np.array([-1,1]))/2).astype(int)
    if end == len(_a):
        end -= 1

    # Set the full range and multiply it by the provided factor
    srt = np.ma.argsort(_a)
    Da = (_a[srt[end]] - _a[srt[start]])*fac

    # Set the midpoint
    mid = _a[srt[len(_a)//2]] if midpoint is None else midpoint

    # Return the range centered on the midpoint
    return [ mid - Da/2, mid + Da/2 ]


def nearest_unmasked(arr, use_indices=False):
    """
    Return the indices of the nearest unmasked element in a vector.

    .. warning::
        The function *uses the values of the masked data* for masked
        elements. This means that if you want to know the nearest
        unmasked element to one of the *masked* elements, the `data`
        attribute of the provided array should have meaningful values
        for these masked elements.

    Args:
        arr (`numpy.ma.MaskedArray`_):
            Array to analyze. Must be 1D.
        use_indices (:obj:`bool`, optional):
            The proximity of each element in the array is based on
            the difference in the array `data` values. Setting
            `use_indices` to `True` instead bases the calculation on
            the proximity of the element indices; i.e., find the
            index of the nearest unmasked element.

    Returns:
        `numpy.ndarray`_: Integer array with the indices of the
        nearest array elements, the definition of which depends on
        `use_indices`.
    """
    # Check the input
    if not isinstance(arr, np.ma.MaskedArray):
        raise TypeError('Must provide a numpy masked array.')
    if arr.ndim != 1:
        raise ValueError('Must be a 1D array.')
    if use_indices:
        return nearest_unmasked(np.ma.MaskedArray(np.arange(arr.size), mask=arr.mask.copy()))

    # Get the difference of each element with every other element
    nearest = np.absolute(arr[None,:]-arr.data[:,None])
    # Ignore the diagonal
    nearest[np.diag_indices(arr.size)] = np.ma.masked
    # Return the location of the minimum value ignoring the masked values
    return np.ma.argmin(nearest, axis=1)


def contiguous_true(m):
    """
    Find contiguous regions of True values in a boolean numpy array.

    This is identically what is done by `numpy.ma.flatnotmasked_contiguous`_,
    except the argument is the mask, not a masked array, and it selects
    contiguous True regions instead of contiguous False regions.

    Args:
        m (array-like):
            A boolean array.  Must be 1D.

    Returns:
        :obj:`list`: A list of slice objects that select contiguous regions of
        True values in the provided array.
    """
    _m = np.atleast_1d(m)
    if _m.ndim > 1:
        raise ValueError('contiguous_true only accepts 1D arrays.')
    if not np.any(_m):
        return [slice(0, _m.size)]
    i = 0
    result = []
    for (k, g) in itertools.groupby(_m.ravel()):
        n = len(list(g))
        if k:
            result.append(slice(i, i + n))
        i += n
    return result


def boxcar_smooth_rows(img, nave, wgt=None, mode='nearest', replace='original'):
    """
    Boxcar smooth an image along their first axis (rows).

    Constructs a boxcar kernel and uses `scipy.ndimage.convolve` to
    smooth the image.  Smoothing does not account for any masking.

    .. note::
        For images following the PypeIt convention, this smooths the
        data spectrally for each spatial position.

    Args:
        img (`numpy.ndarray`_):
            Image to convolve.
        nave (:obj:`int`):
            Number of pixels along rows for smoothing.
        wgt (`numpy.ndarray`_, optional):
            Image providing weights for each pixel in `img`.  Uniform
            weights are used if none are provided.
        mode (:obj:`str`, optional):
            See `scipy.ndimage.convolve`_.

    Returns:
        `numpy.ndarray`_: The smoothed image
    """
    if nave == 1:
        return img
    if img.ndim != 2:
        raise ValueError('Input image must be 2D.')
    if wgt is not None and img.shape != wgt.shape:
        raise ValueError('Input image to smooth and weights must have the same shape.')
    if nave > img.shape[0]:
        msgs.warn('Smoothing box is larger than the image size!')

    # Construct the kernel for mean calculation
    _nave = np.fmin(nave, img.shape[0])
    kernel = np.ones((_nave, 1))/float(_nave)

    if wgt is None:
        # No weights so just smooth
        return scipy.ndimage.convolve(img, kernel, mode='nearest')

    # Weighted smoothing
    cimg = scipy.ndimage.convolve(img*wgt, kernel, mode='nearest')
    wimg = scipy.ndimage.convolve(wgt, kernel, mode='nearest')
    smoothed_img = np.ma.divide(cimg, wimg)
    if replace == 'original':
        smoothed_img[smoothed_img.mask] = img[smoothed_img.mask]
    elif replace == 'zero':
        smoothed_img[smoothed_img.mask] = 0.0
    else:
        msgs.error('Unrecognized value of replace')
    return smoothed_img.data


# TODO: Could this use bisect?
def index_of_x_eq_y(x, y, strict=False):
    """
    Return an index array that maps the elements of `x` to those of
    `y`.

    This should return the index of the *first* element in array `x`
    equal to the associated value in array `y`. Inspired by:
    https://tinyurl.com/yyrx8acf

    Args:
        x (`numpy.ndarray`_):
            1D parent array
        y (`numpy.ndarray`_):
            1D reference array
        strict (:obj:`bool`, optional):
            Raise an exception unless every element of y is found in
            x. I.e., it must be true that::

                np.array_equal(x[index_of_x_eq_y(x,y)], y)

    Returns:
        `numpy.ndarray`_: An array with index of `x` that is equal to
        the given value of `y`.  Output shape is the same as `y`.
    """
    if y.ndim != 1 or y.ndim != 1:
        raise ValueError('Arrays must be 1D.')
    srt = np.argsort(x)
    indx = np.searchsorted(x[srt], y)
    x2y = np.take(srt, indx, mode='clip')
    if strict and not np.array_equal(x[x2y], y):
        raise ValueError('Not every element of y was found in x.')
    return x2y


def rebin(a, newshape):
    """

    Rebin an array to a new shape using slicing. This routine is taken
    from: https://scipy-cookbook.readthedocs.io/items/Rebinning.html.
    The image shapes need not be integer multiples of each other, but in
    this regime the transformation will not be reversible, i.e. if
    a_orig = rebin(rebin(a,newshape), a.shape) then a_orig will not be
    everywhere equal to a (but it will be equal in most places).

    Args:
        a (ndarray, any dtype):
            Image of any dimensionality and data type
        newshape (tuple):
            Shape of the new image desired. Dimensionality must be the
            same as a.

    Returns:
        ndarray: same dtype as input Image with same values as a
        rebinning to shape newshape
    """
    if not len(a.shape) == len(newshape):
        msgs.error('Dimension of a image does not match dimension of new requested image shape')

    slices = [slice(0, old, float(old) / new) for old, new in zip(a.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')  # choose the biggest smaller integer index
    return a[tuple(indices)]

# TODO This function is only used by procimg.lacosmic. Can it be replaced by above?
def rebin_evlist(frame, newshape):
    # This appears to be from
    # https://scipy-cookbook.readthedocs.io/items/Rebinning.html
    shape = frame.shape
    lenShape = len(shape)
    factor = (np.asarray(shape)/np.asarray(newshape)).astype(int)
    evList = ['frame.reshape('] + \
             ['int(newshape[%d]),factor[%d],'% (i, i) for i in range(lenShape)] + \
             [')'] + ['.sum(%d)' % (i+1) for i in range(lenShape)] + \
             ['/factor[%d]' % i for i in range(lenShape)]
    return eval(''.join(evList))




def pyplot_rcparams():
    """
    params for pretty matplotlib plots

    Returns:

    """
    # set some plotting parameters
    plt.rcParams["xtick.top"] = True
    plt.rcParams["ytick.right"] = True
    plt.rcParams["xtick.minor.visible"] = True
    plt.rcParams["ytick.minor.visible"] = True
    plt.rcParams["ytick.direction"] = 'in'
    plt.rcParams["xtick.direction"] = 'in'
    plt.rcParams["xtick.major.size"] = 6
    plt.rcParams["ytick.major.size"] = 6
    plt.rcParams["xtick.minor.size"] = 3
    plt.rcParams["ytick.minor.size"] = 3
    plt.rcParams["xtick.major.width"] = 1
    plt.rcParams["ytick.major.width"] = 1
    plt.rcParams["xtick.minor.width"] = 1
    plt.rcParams["ytick.minor.width"] = 1
    plt.rcParams["axes.linewidth"] = 1
    plt.rcParams["lines.linewidth"] = 3
    plt.rcParams["lines.markeredgewidth"] = 2
    plt.rcParams["patch.linewidth"] = 3
    plt.rcParams["hatch.linewidth"] = 3
    plt.rcParams["font.size"] = 13
    plt.rcParams["legend.frameon"] = False
    plt.rcParams["legend.handletextpad"] = 1


def pyplot_rcparams_default():
    """
    restore default rcparams

    Returns:

    """
    matplotlib.rcParams.update(matplotlib.rcParamsDefault)


def smooth(x, window_len, window='flat'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that edge effects are minimize at the beginning and end part of the signal.

     This code taken from this cookbook and slightly modified: https://scipy-cookbook.readthedocs.io/items/SignalSmooth.html

    .. todo::
        the window parameter could be the window itself if an array instead of a string

    Args:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing., default is 'flat'

    Returns:
        the smoothed signal, same shape as x

    Examples:

        >>> t=linspace(-2,2,0.1)
        >>> x=sin(t)+randn(len(t))*0.1
        >>> y=smooth(x)

    Notes:

        - See also: numpy.hanning, numpy.hamming, numpy.bartlett,
          numpy.blackman, numpy.convolve scipy.signal.lfilter

        - length(output) != length(input), to correct this, return
          y[(window_len/2-1):-(window_len/2)] instead of just y.

    """
    if x.ndim != 1:
        raise ValueError("smooth only accepts 1 dimension arrays.")

    if x.size < window_len:
        raise ValueError("Input vector needs to be bigger than window size.")

    if window_len < 3:
        return x

    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")

    s = np.r_[x[window_len - 1:0:-1], x, x[-2:-window_len - 1:-1]]
    # print(len(s))
    if window == 'flat':  # moving average
        w = np.ones(window_len, 'd')
    else:
        w = eval('np.' + window + '(window_len)')

    y = np.convolve(w / w.sum(), s, mode='same')

    return y[(window_len-1):(y.size-(window_len-1))]


def fast_running_median(seq, window_size):
    """

    Compute the median of sequence of numbers with a running window. The
    boundary conditions are identical to the scipy 'reflect' boundary
    codition:

    'reflect' (`d c b a | a b c d | d c b a`)

    The input is extended by reflecting about the edge of the last pixel.

    This code has been confirmed to produce identical results to
    scipy.ndimage.median_filter with the reflect boundary
    condition, but is ~ 100 times faster.

    Args:
        seq (list or 1-d numpy array of numbers):
        window_size (int): size of running window.

    Returns:
        ndarray: median filtered values

    Code originally contributed by Peter Otten, made to be consistent with
    scipy.ndimage.median_filter by Joe Hennawi.

    Now makes use of the Bottleneck library https://pypi.org/project/Bottleneck/.
    """
    # Enforce that the window_size needs to be smaller than the sequence, otherwise we get arrays of the wrong size
    # upon return (very bad). Added by JFH. Should we print out an error here?

    if (window_size > (len(seq)-1)):
        msgs.warn('window_size > len(seq)-1. Truncating window_size to len(seq)-1, but something is probably wrong....')
    if (window_size < 0):
        msgs.warn('window_size is negative. This does not make sense something is probably wrong. Setting window size to 1')

    window_size = int(np.fmax(np.fmin(int(window_size), len(seq)-1),1))
    # pad the array for the reflection
    seq_pad = np.concatenate((seq[0:window_size][::-1],seq,seq[-1:(-1-window_size):-1]))

    result = move_median.move_median(seq_pad, window_size)

    # This takes care of the offset produced by the original code deducec by trial and error comparison with
    # scipy.ndimage.medfilt

    result = np.roll(result, -window_size//2 + 1)
    return result[window_size:-window_size]


# Taken from stackoverflow
# https://stackoverflow.com/questions/30677241/how-to-limit-cross-correlation-window-width-in-numpy
# slightly modified to return lags

def cross_correlate(x, y, maxlag):
    """

    Cross correlation with a maximum number of lags. This computes the same result as::

        numpy.correlate(x, y, mode='full')[len(a)-maxlag-1:len(a)+maxlag]

    Edges are padded with zeros using ``np.pad(mode='constant')``.

    Args:
        x (ndarray):
            First vector of the cross-correlation.
        y (ndarray):
            Second vector of the cross-correlation. `x` and `y` must be
            one-dimensional numpy arrays with the same length.
        maxlag (int):
            The maximum lag for which to compute the cross-correlation.
            The cross correlation is computed at integer lags from
            (-maxlag, maxlag)

    Returns:
        tuple:  Returns are as follows:
            - lags (ndarray):  shape = (2*maxlag + 1); Lags for the
              cross-correlation. Integer spaced values from (-maxlag,
              maxlag)
            - xcorr (ndarray): shape = (2*maxlag + 1); Cross-correlation
              at the lags
    """

    x = np.asarray(x)
    y = np.asarray(y)
    if x.ndim != 1:
        msgs.error('x must be one-dimensional.')
    if y.ndim != 1:
        msgs.error('y must be one-dimensional.')


    #py = np.pad(y.conj(), 2*maxlag, mode=mode)
    py = np.pad(y, 2*maxlag, mode='constant')
    T = as_strided(py[2*maxlag:], shape=(2*maxlag+1, len(y) + 2*maxlag),
                   strides=(-py.strides[0], py.strides[0]))
    px = np.pad(x, maxlag, mode='constant')
    lags = np.arange(-maxlag, maxlag + 1,dtype=float)
    return lags, T.dot(px)


def clip_ivar(flux, ivar, sn_clip, gpm=None, verbose=False):
    """
    Add an error floor the the inverse variance array.

    This is primarily to prevent too much rejection at high-S/N (i.e.
    standard stars, bright objects).

    Args:
        flux (`numpy.ndarray`_):
            Flux array
        ivar (`numpy.ndarray`_):
            Inverse variance array
        sn_clip (:obj:`float`):
            This sets the small erorr that is added to the input ``ivar``
            such that the output inverse variance will never give S/N greater
            than ``sn_clip``. This prevents overly aggressive rejection in
            high S/N spectra, which nevertheless differ at a level greater
            than the formal S/N due to systematics. If None, the input
            inverse variance array is simply returned.
        gpm (`numpy.ndarray`_, optional):
            Good-pixel mask for the input fluxes.
        verbose (:obj:`bool`, optional):
            Write status messages to the terminal.

    Returns:
         `numpy.ndarray`_: The new inverse variance matrix that yields a S/N
         upper limit.
    """
    if sn_clip is None:
        return ivar

    if verbose:
        msgs.info('Inflating errors to keep S/N ratio below S/N_clip = {:5.3f}'.format(sn_clip))

    _gpm = ivar > 0.
    if gpm is not None:
        _gpm &= gpm
    adderr = 1.0/sn_clip
    ivar_cap = _gpm/(1.0/(ivar + np.logical_not(_gpm)) + adderr**2*(np.abs(flux))**2)
    return np.minimum(ivar, ivar_cap)


def inverse(array):
    """
    Calculate and return the inverse of the input array, enforcing
    positivity and setting values <= 0 to zero.  The input array should
    be a quantity expected to always be positive, like a variance or an
    inverse variance. The quantity::

        out = (array > 0.0)/(np.abs(array) + (array == 0.0))

    is returned.

    Args:
        array (`numpy.ndarray`_):
            Array to invert

    Returns:
        `numpy.ndarray`_: Result of controlled ``1/array`` calculation.
    """
    return (array > 0.0)/(np.abs(array) + (array == 0.0))


def calc_ivar(varframe):
    """

    Calculate the inverse variance based on the input array

    Wrapper to inverse()

    Args:
        varframe (ndarray):  Variance image

    Returns:
        ndarray:  Inverse variance image
    """
    # THIS WILL BE DEPRECATED!!
    return inverse(varframe)




def robust_meanstd(array):
    """
    Determine a robust measure of the mean and dispersion of array

    Args:
        array (ndarray): an array of values

    Returns:
        tuple: Median of the array and a robust estimate of the standand
        deviation (assuming a symmetric distribution).
    """
    med = np.median(array)
    mad = np.median(np.abs(array-med))
    return med, 1.4826*mad



def polyfitter2d(data, mask=None, order=2):
    """
    2D fitter

    Args:
        data:
        mask:
        order:

    Returns:

    """
    x, y = np.meshgrid(np.linspace(0.0, 1.0, data.shape[1]), np.linspace(0.0, 1.0, data.shape[0]))
    if isinstance(mask, (float, int)):
        # mask is the value that should be masked in data
        w = np.where(data != mask)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    elif mask is None or mask.size == 0:
        # There are no masks
        xf = x.flatten()
        yf = y.flatten()
        m = polyfit2d(xf, yf, data.T.flatten(), order)
    elif len(mask.shape) == 1:
        # mask is applied along one axis
        mskar = np.ones((data.shape[0], data.shape[1]))
        mskar[mask, :] = 0
        w = np.where(mskar == 1)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    elif mask.shape[0] == data.shape[0] and mask.shape[1] == data.shape[1]:
        # mask is an array that indicates the masked data
        w = np.where(mask == 0)
        xf = x[w].flatten()
        yf = y[w].flatten()
        m = polyfit2d(xf, yf, data[w].T.flatten(), order)
    # Return the best model
    return m, polyval2d(x, y, m).T


def polyfit2d(x, y, z, order=3):
    """
    Generate 2D polynomial

    Args:
        x:
        y:
        z:
        order:

    Returns:

    """
    ncols = (order + 1)**2
    G = np.zeros((x.size, ncols))
    ij = itertools.product(range(order+1), range(order+1))
    for k, (i,j) in enumerate(ij):
        G[:,k] = x**i * y**j
    m, null, null, null = np.linalg.lstsq(G, z)
    return m


def polyval2d(x, y, m):
    """
    Generate 2D polynomial

    Args:
        x:
        y:
        m:

    Returns:

    """
    order = int(np.sqrt(len(m))) - 1
    ij = itertools.product(range(order+1), range(order+1))
    z = np.zeros_like(x)
    for a, (i, j) in zip(m, ij):
        z += a * x**i * y**j
    return z


def subsample(frame):
    """
    Used by LACosmic

    Args:
        frame (ndarray):

    Returns:
        ndarray: Sliced image

    """
    newshape = (2*frame.shape[0], 2*frame.shape[1])
    slices = [slice(0, old, float(old)/new) for old, new in zip(frame.shape, newshape)]
    coordinates = np.mgrid[slices]
    indices = coordinates.astype('i')
    return frame[tuple(indices)]


def find_nearest(array, values):
    """For all elements of values, find the index of the nearest value in array

    Parameters
    ----------
    array : numpy.ndarray
        Array of values
    values : numpy.ndarray
        Values to be compared with the elements of `array`

    Return
    ------
    numpy.ndarray : indices of `array` that are closest to each element of value
    """
    # Make sure the input is a numpy array
    array = np.array(array)

    # get insert positions
    idxs = np.searchsorted(array, values, side="left")

    # find indexes where previous index is closer
    prev_idx_is_less = ((idxs == len(array)) | (np.fabs(values - array[np.maximum(idxs - 1, 0)]) <
                                                np.fabs(values - array[np.minimum(idxs, len(array) - 1)])))
    idxs[prev_idx_is_less] -= 1

    return idxs


def yamlify(obj, debug=False):
    """

    Recursively process an object so it can be serialised for yaml.

    Based on jsonify in `linetools`_.

    Also found in desiutils

    Note:
        All string-like keys in :class:`dict` s are converted to
        :class:`str`.

    Parameters
    ----------
    obj : :class:`object`
        Any object.
    debug : :class:`bool`, optional
        Print extra information if requested.

    Returns
    -------
    obj: :class:`object`
        An object suitable for yaml serialization.  For example
        :class:`numpy.ndarray` is converted to :class:`list`,
        :class:`numpy.int64` is converted to :class:`int`, etc.
    """
    # TODO: Change to np.floating?
    if isinstance(obj, (np.float64, np.float32)):
        obj = float(obj)
    # TODO: Change to np.integer?
    elif isinstance(obj, (np.int32, np.int64, np.int16)):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
#    elif isinstance(obj, bytes):
#        obj = obj.decode('utf-8')
    elif isinstance(obj, (np.string_, str)):
        # Worry about colons!
        if ':' in obj:
            obj = '"'+str(obj)+'"'
        else:
            obj = str(obj)
    elif isinstance(obj, units.Quantity):
        try:
            obj = obj.value.tolist()
        except AttributeError:
            obj = obj.value
    elif isinstance(obj, np.ndarray):  # Must come after Quantity
        obj = obj.tolist()
    elif isinstance(obj, dict):
        # First convert keys
        nobj = {}
        for key, value in obj.items():
            if isinstance(key, str):
                nobj[str(key)] = value
            else:
                nobj[key] = value
        # Now recursive
        obj = nobj
        for key, value in obj.items():
            obj[key] = yamlify(value, debug=debug)
    elif isinstance(obj, list):
        for i, item in enumerate(obj):
            obj[i] = yamlify(item, debug=debug)
    elif isinstance(obj, tuple):
        obj = list(obj)
        for i, item in enumerate(obj):
            obj[i] = yamlify(item, debug=debug)
        obj = tuple(obj)
    # elif isinstance(obj, Unit):
    #     obj = obj.name
    # elif obj is units.dimensionless_unscaled:
    #     obj = 'dimensionless_unit'
    if debug:
        print(type(obj))
    return obj


def add_sub_dict(d, key):
    """
    If a key is not present in the provided dictionary, add it as a new nested
    dictionary.

    Args:
        d (:obj:`dict`):
            Dictionary to alter
        key (:obj:`str`):
            Key to add

    Examples:
        >>> d = {}
        >>> add_sub_dict(d, 'test')
        >>> d
        {'test': {}}
        >>> d['test'] = 'this'
        >>> add_sub_dict(d, 'test')
        >>> d
        {'test': 'this'}
        >>> add_sub_dict(d, 'and')
        >>> d['and'] = 'that'
        >>> d
        {'test': 'this', 'and': 'that'}
    """
    if key not in d.keys():
        d[key] = {}


def recursive_update(d, u):
    """
    Update dictionary values with recursion to nested dictionaries.

    Thanks to:
    https://stackoverflow.com/questions/3232943/update-value-of-a-nested-dictionary-of-varying-depth

    Args:
        d (:obj:`dict`):
            Dictionary (potentially of other dictionaries) to be updated.  This
            is both edited in-place and returned.
        u (:obj:`dict`):
            Dictionary (potentially of other dictionaries) with the
            updated/additional values.

    Returns:
        :obj:`dict`: The updated dictionary.
    """
    for k, v in u.items():
        d[k] = recursive_update(d.get(k, {}), v) if isinstance(v, collections.abc.Mapping) else v
    return d


def save_pickle(fname, obj):
    """Save an object to a python pickle file

    Parameters
    ----------
    fname : :class:`str`
        Filename
    obj : :class:`object`
        An object suitable for pickle serialization.
    """
    if fname.split(".")[-1] != 'pkl':
        fname += '.pkl'
    with open(fname, 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)
        msgs.info('File saved: {0:s}'.format(fname))


def load_pickle(fname):
    """Load a python pickle file

    Parameters
    ----------
    fname : :class:`str`
        Filename

    Returns
    -------
    :class:`object`
        An object suitable for pickle serialization.
    """
    msgs.info('Loading file: {0:s}'.format(fname))
    with open(fname, 'rb') as f:
        return pickle.load(f)


##
##This code was originally published by the following individuals for use with
##Scilab:
##    Copyright (C) 2012 - 2013 - Michael Baudin
##    Copyright (C) 2012 - Maria Christopoulou
##    Copyright (C) 2010 - 2011 - INRIA - Michael Baudin
##    Copyright (C) 2009 - Yann Collette
##    Copyright (C) 2009 - CEA - Jean-Marc Martinez

##   website: forge.scilab.org/index.php/p/scidoe/sourcetree/master/macros
##Much thanks goes to these individuals. It has been converted to Python by
##Abraham Lee.
##"

## Python version taken from https://pythonhosted.org/pyDOE/randomized.html by JFH

def lhs(n, samples=None, criterion=None, iterations=None):
    """
    Generate a latin-hypercube design

    Parameters
    ----------
    n : int
        The number of factors to generate samples for

    Optional
    --------
    samples : int
        The number of samples to generate for each factor (Default: n)
    criterion : str
        Allowable values are "center" or "c", "maximin" or "m",
        "centermaximin" or "cm", and "correlation" or "corr". If no value
        given, the design is simply randomized.
    iterations : int
        The number of iterations in the maximin and correlations algorithms
        (Default: 5).

    Returns
    -------
    H : 2d-array
        An n-by-samples design matrix that has been normalized so factor values
        are uniformly spaced between zero and one.

    Example
    -------
    A 3-factor design (defaults to 3 samples)::

        >>> lhs(3)
        array([[ 0.40069325,  0.08118402,  0.69763298],
               [ 0.19524568,  0.41383587,  0.29947106],
               [ 0.85341601,  0.75460699,  0.360024  ]])

    A 4-factor design with 6 samples::

        >>> lhs(4, samples=6)
        array([[ 0.27226812,  0.02811327,  0.62792445,  0.91988196],
               [ 0.76945538,  0.43501682,  0.01107457,  0.09583358],
               [ 0.45702981,  0.76073773,  0.90245401,  0.18773015],
               [ 0.99342115,  0.85814198,  0.16996665,  0.65069309],
               [ 0.63092013,  0.22148567,  0.33616859,  0.36332478],
               [ 0.05276917,  0.5819198 ,  0.67194243,  0.78703262]])

    A 2-factor design with 5 centered samples::

        >>> lhs(2, samples=5, criterion='center')
        array([[ 0.3,  0.5],
               [ 0.7,  0.9],
               [ 0.1,  0.3],
               [ 0.9,  0.1],
               [ 0.5,  0.7]])

    A 3-factor design with 4 samples where the minimum distance between
    all samples has been maximized::

        >>> lhs(3, samples=4, criterion='maximin')
        array([[ 0.02642564,  0.55576963,  0.50261649],
               [ 0.51606589,  0.88933259,  0.34040838],
               [ 0.98431735,  0.0380364 ,  0.01621717],
               [ 0.40414671,  0.33339132,  0.84845707]])

    A 4-factor design with 5 samples where the samples are as uncorrelated
    as possible (within 10 iterations)::

        >>> lhs(4, samples=5, criterion='correlate', iterations=10)

    """
    H = None

    if samples is None:
        samples = n

    if criterion is not None:
        assert criterion.lower() in ('center', 'c', 'maximin', 'm',
                                     'centermaximin', 'cm', 'correlation',
                                     'corr'), 'Invalid value for "criterion": {}'.format(criterion)
    else:
        H = _lhsclassic(n, samples)

    if criterion is None:
        criterion = 'center'

    if iterations is None:
        iterations = 5

    if H is None:
        if criterion.lower() in ('center', 'c'):
            H = _lhscentered(n, samples)
        elif criterion.lower() in ('maximin', 'm'):
            H = _lhsmaximin(n, samples, iterations, 'maximin')
        elif criterion.lower() in ('centermaximin', 'cm'):
            H = _lhsmaximin(n, samples, iterations, 'centermaximin')
        elif criterion.lower() in ('correlate', 'corr'):
            H = _lhscorrelate(n, samples, iterations)

    return H

################################################################################

def _lhsclassic(n, samples):
    # Generate the intervals
    cut = np.linspace(0, 1, samples + 1)

    # Fill points uniformly in each interval
    u = np.random.rand(samples, n)
    a = cut[:samples]
    b = cut[1:samples + 1]
    rdpoints = np.zeros_like(u)
    for j in range(n):
        rdpoints[:, j] = u[:, j ] *( b -a) + a

    # Make the random pairings
    H = np.zeros_like(rdpoints)
    for j in range(n):
        order = np.random.permutation(range(samples))
        H[:, j] = rdpoints[order, j]

    return H

################################################################################

def _lhscentered(n, samples):
    # Generate the intervals
    cut = np.linspace(0, 1, samples + 1)

    # Fill points uniformly in each interval
    u = np.random.rand(samples, n)
    a = cut[:samples]
    b = cut[1:samples + 1]
    _center = (a + b ) /2

    # Make the random pairings
    H = np.zeros_like(u)
    for j in range(n):
        H[:, j] = np.random.permutation(_center)

    return H

################################################################################

def _lhsmaximin(n, samples, iterations, lhstype):
    maxdist = 0

    # Maximize the minimum distance between points
    for i in range(iterations):
        if lhstype=='maximin':
            Hcandidate = _lhsclassic(n, samples)
        else:
            Hcandidate = _lhscentered(n, samples)

        d = _pdist(Hcandidate)
        if maxdist <np.min(d):
            maxdist = np.min(d)
            H = Hcandidate.copy()

    return H

################################################################################

def _lhscorrelate(n, samples, iterations):
    mincorr = np.inf

    # Minimize the components correlation coefficients
    for i in range(iterations):
        # Generate a random LHS
        Hcandidate = _lhsclassic(n, samples)
        R = np.corrcoef(Hcandidate)
        if np.max(np.abs(R[ R!=1]) ) <mincorr:
            mincorr = np.max(np.abs( R -np.eye(R.shape[0])))
            print('new candidate solution found with max,abs corrcoef = {}'.format(mincorr))
            H = Hcandidate.copy()

    return H

################################################################################

def _pdist(x):
    """
    Calculate the pair-wise point distances of a matrix

    Parameters
    ----------
    x : 2d-array
        An m-by-n array of scalars, where there are m points in n dimensions.

    Returns
    -------
    d : array
        A 1-by-b array of scalars, where b = m*(m - 1)/2. This array contains
        all the pair-wise point distances, arranged in the order (1, 0),
        (2, 0), ..., (m-1, 0), (2, 1), ..., (m-1, 1), ..., (m-1, m-2).

    Examples
    --------
    ::

        >>> x = np.array([[0.1629447, 0.8616334],
        ...               [0.5811584, 0.3826752],
        ...               [0.2270954, 0.4442068],
        ...               [0.7670017, 0.7264718],
        ...               [0.8253975, 0.1937736]])
        >>> _pdist(x)
        array([ 0.6358488,  0.4223272,  0.6189940,  0.9406808,  0.3593699,
                0.3908118,  0.3087661,  0.6092392,  0.6486001,  0.5358894])

    """

    x = np.atleast_2d(x)
    assert len(x.shape )==2, 'Input array must be 2d-dimensional'

    m, n = x.shape
    if m< 2:
        return []

    d = []
    for i in range(m - 1):
        for j in range(i + 1, m):
            d.append((sum((x[j, :] - x[i, :]) ** 2)) ** 0.5)

    return np.array(d)


def is_float(s):
    """
    Detertmine if a string can be converted to a floating point number.
    """
    try:
        float(s)
    except:
        return False

    return True

def find_single_file(file_pattern) -> pathlib.Path:
    """Find a single file matching a wildcard pattern.

    Args:
        file_pattern (str): A filename pattern, see the python 'glob' module.

    Returns:
        :obj:`pathlib.Path`: A file name, or None if no filename was found. This will give a warning
             if multiple files are found and return the first one.
    """

    files = glob.glob(file_pattern)
    if len(files) == 1:
        return pathlib.Path(files[0])
    elif len(files) == 0:
        return None
    else:
        msgs.warn(f'Found multiple files matching {file_pattern}; using the first one.')
        return pathlib.Path(files[0])

def DFS(v: int, visited: list[bool], group: list[int], adj: np.ndarray):
    """
    Depth-First Search of graph given by matrix `adj` starting from `v`.
    Updates `visited` and `group`.

    Args:
        v (int): initial vertex
        visited (List[bool]): List keeping track of which vertices have been
            visited at any point in traversing the graph. `visited[i]` is True
            iff vertix `i` has been visited before.
        group (List[int]): List keeping track of which vertices have been
            visited in THIS CALL of DFS. After DFS returns, `group` contains
            all members of the connected component containing v. `i in group`
            is True iff vertex `i` has been visited in THIS CALL of DFS.
        adj (np.ndarray): Adjacency matrix description of the graph. `adj[i,j]`
            is True iff there is a vertex between `i` and `j`.
    """
    stack = []
    stack.append(v)
    while stack:
        u = stack.pop()
        if not visited[u]:
            visited[u] = True
            group.append(u)
            neighbors = [i for i in range(len(adj[u])) if adj[u,i]]
            for neighbor in neighbors:
                stack.append(neighbor)

def list_of_spectral_lines():
    """ Generate a list of spectral lines

    Returns:
        tuple: np.ndarray, np.ndarray
    """
    # spectral features
    CIVnam1, CIVwav1='CIV', 1548.
    CIVnam2, CIVwav2='CIV', 1550.

    HeIInam0, HeIIwav0='HeII', 1640.

    OIIInam01, OIIIwav01='OIII]', 1661.
    OIIInam02, OIIIwav02='OIII]', 1666.

    SiIIInam1, SiIIIwav1='SiIII]', 1882.
    SiIIInam2, SiIIIwav2='SiIII]', 1892.

    CIIInam1, CIIIwav1='CIII]', 1907.
    CIIInam2, CIIIwav2='CIII]', 1909.

    Lyalphanam, Lyalphawav='Lyalpha', 1215.7
    OIInam1, OIIwav1='[OII]', 3726
    OIInam2, OIIwav2='[OII]', 3729
    OIIInam1, OIIIwav1='[OIII]', 5007.
    OIIInam2, OIIIwav2='[OIII]', 4959.
    OIIInam3, OIIIwav3='[OIII]', 4363.
    Halphanam, Halphawav='Halpha', 6563.
    Hbetanam, Hbetawav='Hbeta', 4861.
    Hdeltanam, Hdeltawav='Hdelta', 4101.
    Hgammanam, Hgammawav='Hgamma', 4341.

    NeIIInam, NeIIIwav = '[NeIII]', 3869.
    NeVnam, NeVwav = '[NeV]', 3426.
    SIInam1, SIIwav1 = '[SII]', 6716.
    SIInam2, SIIwav2 = '[SII]', 6716.


    ##absorption
    H13nam, H13wav = 'H13', 3734.
    H12nam, H12wav = 'H12', 3750.
    H11nam, H11wav = 'H11', 3771.
    H10nam, H10wav = 'H10', 3798.
    H9nam, H9wav = 'H9', 3835.
    H8nam, H8wav = 'H8', 3889.
    HeInam, HeIwav = 'HeI', 3889.

    CAII_Knam, CaII_Kwav = 'CaK', 3934.
    CAII_Hnam, CaII_Hwav = 'CaH', 3968.

    Gbandnam, Gbandwav = 'Gband', 4305.

    line_names=np.array([CIVnam1, CIVnam2, HeIInam0, OIIInam01, OIIInam02, SiIIInam1, SiIIInam2,
                         CIIInam1, CIIInam2, Lyalphanam, OIInam1, OIInam2, OIIInam1, OIIInam2,
                         OIIInam3, Halphanam, Hbetanam, Hdeltanam, Hgammanam, NeIIInam, NeVnam,
                         SIInam1, SIInam2, H13nam, H12nam, H11nam, H10nam, H9nam, H8nam, HeInam,
                         CAII_Knam, CAII_Hnam, Gbandnam])

    line_wav = np.array([CIVwav2, CIVwav2, HeIIwav0, OIIIwav01, OIIIwav02, SiIIIwav1, SiIIIwav2,
                         CIIIwav1, CIIIwav2, Lyalphawav, OIIwav1, OIIwav2, OIIIwav1, OIIIwav2,
                         OIIIwav3, Halphawav, Hbetawav, Hdeltawav, Hgammawav, NeIIIwav, NeVwav,
                         SIIwav1,SIIwav2, H13wav, H12wav, H11wav, H10wav, H9wav, H8wav, HeIwav,
                         CaII_Kwav, CaII_Hwav, Gbandwav])

    return line_names, line_wav
