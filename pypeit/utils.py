"""
General utility functions.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import os
import inspect
import pickle
import warnings
import itertools
from collections import deque
from bisect import insort, bisect_left

from IPython import embed

import numpy as np
from numpy.lib.stride_tricks import as_strided

from scipy import interpolate, ndimage

import matplotlib
from matplotlib import pyplot as plt

from astropy import units
from astropy import stats

from pypeit.core import pydl
from pypeit import msgs

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


def string_table(tbl, delimeter='print'):
    """
    Provided the array of data, format it with equally spaced columns
    and add a header (first row) and contents delimeter.

    Args:
        tbl (`numpy.ndarray`_):
            Array of string representations of the data to print.
        delimeter (:obj:`str`, optional):
            Delimeter between first table row, which should contain
            the column headings, and the column data. Use ``'print'``
            for a simple line of hyphens, anything else results in an
            ``rst`` style table formatting.

    Returns:
        :obj:`str`: Single long string with the data table.
    """
    nrows, ncols = tbl.shape
    col_width = [np.amax([len(dij) for dij in dj]) for dj in tbl.T]
    row_string = ['']*(nrows+1) if delimeter == 'print' else ['']*(nrows+3)
    start = 2 if delimeter == 'print' else 3
    for i in range(start,nrows+start-1):
        row_string[i] = '  '.join([tbl[1+i-start,j].ljust(col_width[j]) for j in range(ncols)])
    if delimeter == 'print':
        # Heading row
        row_string[0] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
        # Delimiter
        row_string[1] = '-'*len(row_string[0])
        return '\n'.join(row_string)+'\n'

    # For an rst table
    row_string[0] = '  '.join([ '='*col_width[j] for j in range(ncols)])
    row_string[1] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
    row_string[2] = row_string[0]
    row_string[-1] = row_string[0]
    return '\n'.join(row_string)+'\n'


def spec_atleast_2d(wave, flux, ivar, mask):
    """
    Utility routine to repackage spectra to have shape (nspec, norders) or (nspec, ndetectors) or (nspec, nexp)

    Args:
        wave (`numpy.ndarray`_):
            Wavelength array
        flux (`numpy.ndarray`_):
            Flux array
        ivar (`numpy.ndarray`_):
            Inverse variance array
        mask (`numpy.ndarray`_, bool):
            Good pixel mask True=Good.

    Returns:
        wave_arr, flux_arr, ivar_arr, mask_arr, nspec, norders

            Reshaped arrays which all have shape (nspec, norders) or (nspec, ndetectors) or (nspec, nexp) along
            with nspec, and norders = total number of orders, detectors, or exposures


    """
    # Repackage the data into arrays of shape (nspec, norders)
    if flux.ndim == 1:
        nspec = flux.size
        norders = 1
        wave_arr = wave.reshape(nspec, 1)
        flux_arr = flux.reshape(nspec, 1)
        ivar_arr = ivar.reshape(nspec, 1)
        mask_arr = mask.reshape(nspec, 1)
    else:
        nspec, norders = flux.shape
        if wave.ndim == 1:
            wave_arr = np.tile(wave, (norders, 1)).T
        else:
            wave_arr = wave
        flux_arr = flux
        ivar_arr = ivar
        mask_arr = mask

    return wave_arr, flux_arr, ivar_arr, mask_arr, nspec, norders


def nan_mad_std(data, axis=None, func=None):
    """

    Wrapper for astropy.stats.mad_stad which ignores nans, so as to
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
        return ndimage.convolve(img, kernel, mode='nearest')

    # Weighted smoothing
    cimg = ndimage.convolve(img*wgt, kernel, mode='nearest')
    wimg = ndimage.convolve(wgt, kernel, mode='nearest')
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
    factor = np.asarray(shape)/np.asarray(newshape)
    evList = ['frame.reshape('] + \
             ['int(newshape[%d]),int(factor[%d]),'% (i, i) for i in range(lenShape)] + \
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
    scipy.ndimage.filters.median_filter with the reflect boundary
    condition, but is ~ 100 times faster.

    Args:
        seq (list or 1-d numpy array of numbers):
        window_size (int): size of running window.

    Returns:
        ndarray: median filtered values

    Code contributed by Peter Otten, made to be consistent with
    scipy.ndimage.filters.median_filter by Joe Hennawi.

    See discussion at:
    http://groups.google.com/group/comp.lang.python/browse_thread/thread/d0e011c87174c2d0
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

    seq_pad = iter(seq_pad)
    d = deque()
    s = []
    result = []
    for item in itertools.islice(seq_pad, window_size):
        d.append(item)
        insort(s, item)
        result.append(s[len(d)//2])
    m = window_size // 2
    for item in seq_pad:
        old = d.popleft()
        d.append(item)
        del s[bisect_left(s, old)]
        insort(s, item)
        result.append(s[m])

    # This takes care of the offset produced by the original code deducec by trial and error comparison with
    # scipy.ndimage.filters.medfilt

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



def clip_ivar(flux, ivar, sn_clip, mask=None, verbose=False):
    """

    This adds an error floor to the ivar, preventing too much rejection
    at high-S/N (i.e. standard stars, bright objects)

    Args:
        flux (ndarray):
            flux array
        ivar (ndarray):
            ivar array
        sn_clip (float):
            Small erorr is added to input ivar so that the output ivar_out will never give S/N greater than sn_clip.
            This prevents overly aggressive rejection in high S/N ratio spectra which neverthless differ at a
            level greater than the formal S/N due to systematics.

        mask (ndarray, bool): mask array, True=good

    Returns:
         ndarray: new ivar array
    """
    if sn_clip is None:
        return ivar
    else:
        if mask is None:
            mask = (ivar > 0.0)
        adderr = 1.0/sn_clip
        gmask = (ivar > 0) & mask
        ivar_cap = gmask/(1.0/(ivar + np.invert(gmask)) + adderr**2*(np.abs(flux))**2)
        ivar_out = np.minimum(ivar, ivar_cap)
        if verbose:
            msgs.info('Adding error to ivar to keep S/N ratio below S/N_clip = {:5.3f}'.format(sn_clip))
        return ivar_out


def inverse(array):
    """

    Calculate and return the inverse of the input array, enforcing
    positivity and setting values <= 0 to zero.  The input array should
    be a quantity expected to always be positive, like a variance or an
    inverse variance. The quantity::

        out = (array > 0.0)/(np.abs(array) + (array == 0.0))

    is returned.

    Args:
        a (np.ndarray):

    Returns:
        np.ndarray:

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








def robust_polyfit(xarray, yarray, order, weights=None, maxone=True, sigma=3.0,
                   function="polynomial", initialmask=None, forceimask=False,
                   minx=None, maxx=None, guesses=None, bspline_par=None, verbose=True):
    """
    A robust (equally weighted) polynomial fit is performed to the xarray, yarray pairs
    mask[i] = 1 are masked values

    :param xarray: independent variable values
    :param yarray: dependent variable values
    :param order: the order of the polynomial to be used in the fitting
    :param weights: weights to be used in the fitting (weights = 1/sigma)
    :param maxone: If True, only the most deviant point in a given iteration will be removed
    :param sigma: confidence interval for rejection
    :param function: which function should be used in the fitting (valid inputs: 'polynomial', 'legendre', 'chebyshev', 'bspline')
    :param initialmask: a mask can be supplied as input, these values will be masked for the first iteration. 1 = value masked
    :param forceimask: if True, the initialmask will be forced for all iterations
    :param minx: minimum value in the array (or the left limit for a legendre/chebyshev polynomial)
    :param maxx: maximum value in the array (or the right limit for a legendre/chebyshev polynomial)
    :return: mask, ct -- mask is an array of the masked values, ct is the coefficients of the robust polyfit.
    """
    # Setup the initial mask
    if initialmask is None:
        mask = np.zeros(xarray.size, dtype=np.int)
        if forceimask:
            msgs.warn("Initial mask cannot be enforced -- no initital mask supplied")
            forceimask = False
    else:
        mask = initialmask.copy()
    mskcnt = np.sum(mask)
    # Iterate, and mask out new values on each iteration
    ct = guesses

    while True:
        w = np.where(mask == 0)
        xfit = xarray[w]
        yfit = yarray[w]
        if weights is not None:
            wfit = weights[w]
        else:
            wfit = None
        ct = func_fit(xfit, yfit, function, order, w=wfit,
                      guesses=ct, minx=minx, maxx=maxx, bspline_par=bspline_par)
        yrng = func_val(ct, xarray, function, minx=minx, maxx=maxx)
        sigmed = 1.4826*np.median(np.abs(yfit-yrng[w]))
        #if xarray.size-np.sum(mask) <= order+2: JFH fixed this bug
        if xarray.size - np.sum(mask) <= order + 1:
            if verbose:
                msgs.warn("More parameters than data points - fit might be undesirable")
            break  # More data was masked than allowed by order
        if maxone:  # Only remove the most deviant point
            tst = np.abs(yarray[w]-yrng[w])
            m = np.argmax(tst)
            if tst[m] > sigma*sigmed:
                mask[w[0][m]] = 1
        else:
            if forceimask:
                w = np.where((np.abs(yarray-yrng) > sigma*sigmed) | (initialmask==1))
            else:
                w = np.where(np.abs(yarray-yrng) > sigma*sigmed)
            mask[w] = 1
        if mskcnt == np.sum(mask): break  # No new values have been included in the mask
        mskcnt = np.sum(mask)

    # Final fit
    w = np.where(mask == 0)
    xfit = xarray[w]
    yfit = yarray[w]
    if weights is not None:
        wfit = weights[w]
    else:
        wfit = None
    ct = func_fit(xfit, yfit, function, order, w=wfit, minx=minx, maxx=maxx, bspline_par=bspline_par)
    return mask, ct

def robust_optimize(ydata, fitfunc, arg_dict, maxiter=10, inmask=None, invvar=None,
                    lower=None, upper=None, maxdev=None, maxrej=None, groupdim=None,
                    groupsize=None, groupbadpix=False, grow=0, sticky=True, use_mad=False,
                    verbose=False,
                    **kwargs_optimizer):
    """
    A routine to perform robust optimization. It is completely analogous
    to :func:`robust_polyfit_djs`, but is more general in that it allows
    one to fit a more general model using the optimizer of the users
    choice. If you are fitting simple functions like Chebyshev or
    Legednre polynomials using a linear least-squares algorithm, you
    should use :func:robust_polyfit_djs` instead of this function.

    Args:
        ydata (`numpy.ndarray`_):
            Data to fit.
        fitfunc (callable):
            The callable object used to perform the fitting.  The
            calling sequence must be::

                ret_tuple = fitfunc(ydata, inmask, arg_dict, **kwargs_optimizer)

            See the descriptions of `ydata`, `inmask`, `arg_dict`, and
            `kwargs_optimizer`.  The returned object ret_tuple that can
            have two or three elements.  If it has two elements (result,
            ymodel):

                - `result`: Object returned by the specific
                  scipy.optimize method used to perform the fit.
                - `ymodel`: A `numpy.ndarray` with the model fit to
                  `ydata` and with the same shape.

            If it has three elements (result, ymodel, newivar):

                - `newivar`: new inverse variance for the ydata ymodel
                  comparison, in other words chi = (ydata -
                  ymodel)*np.sqrt(newivar). This functionality allows
                  one to deal with cases where the noise of the
                  data-model comaprison is model dependent.

        arg_dict (:obj:`dict`):
            Dictionary containing the other variables needed to evaluate
            the model fit.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations.  Set this to zero to
            disable rejection and simply do a fit.
        inmask (`numpy.ndarray`_, optional):
            Input mask.  Bad points are marked with a value that
            evaluates to `False`.  Must have the same number of
            dimensions as `ydata`.  Points masked as `False` in `inmask`
            will also always evaluate to `False` in the output mask.
        invvar (:obj:`float`, `numpy.ndarray`_, optional):
            Inverse variance of the data, used to reject points based on
            the values of `upper` and `lower`.  This can either be a
            single float for the entire yarray or a ndarray with the
            same shape as the yarray.
        lower (:obj:`int`, :obj:`float`, optional):
            If set, reject points with ``data < model - lower * sigma``, where
            ``sigma = 1/sqrt(invvar)``
        upper (:obj:`int`, :obj:`float`, optional):
            If set, reject points with ``data > model + upper * sigma``, where
            ``sigma = 1/sqrt(invvar)``.
        maxdev (:obj:`int` or :class:`float`, optional):
            If set, reject points with ``abs(data-model) > maxdev``.  It
            is permitted to set all three of `lower`, `upper` and
            `maxdev`.
        maxrej (:obj:`int`, `numpy.ndarray`_, optional):
            Maximum number of points to reject in this iteration.  If
            `groupsize` or `groupdim` are set to arrays, this should be
            an array, as well.
        groupdim (:obj:`int`, optional):
            Dimension along which to group the data. Set to 1 to group
            along the 1st dimension, 2 for the 2nd dimension, etc.  For
            example, if data has shape [100,200], then setting
            `groupdim=2` is equivalent to grouping the data with
            `groupsize=100`.  In either case, there are 200 groups,
            specified by `[*,i]`.  This functionality is **not well
            tested in python**!
        groupsize (:obj:`int`, optional):
            If this and `maxrej` are set, then reject a maximum of
            `maxrej` points per group of `groupsize` points.  If
            `groupdim` is also set, then this specifies sub-groups
            within that.  This functionality is **not well tested in
            python**!
        groupbadpix (:obj:`bool`, optional):
            If `True`, consecutive sets of bad pixels are considered
            groups, overriding the values of `groupsize`.
        grow (:obj:`int`, optional):
            If set to a non-zero integer, N, the N nearest neighbors of
            rejected pixels will also be rejected.
        sticky (:obj:`bool`, optional):
            If `True`, pixels rejected in one iteration remain rejected
            in subsequent iterations, even if the model changes.
        use_mad (:obj:`bool`, optional):
            It `True`, compute the median of the maximum absolute
            deviation between the data and use this for the rejection
            instead of the default, which is to compute the standard
            deviation of `ydata - modelfit`. Note that it is not
            possible to specify `use_mad=True` and also pass in a value for
            `invvar`, and the code will return an error if this is done.
        **kwargs_optimizer:
            Optional parameters passed to the optimizer.

    Returns:
        Three objects are returned:
            - The object returned by the `scipy.optimize` function used
              by the fitter.  See `fitfunc`.
            - A `numpy.ndarray`_ with the model value fit to `ydata` and
              has its same shape.
            - Boolean `numpy.ndarray`_ with the same shape as data
              indicating which pixels were masked in the final fit.
              Convention is that `True` are good values where `False`
              indicates bad values.

    """
    # Setup the initial mask
    if inmask is None:
        inmask = np.ones(ydata.size, dtype=bool)

    nin_good = np.sum(inmask)
    iter = 0
    qdone = False
    thismask = np.copy(inmask)

    while (not qdone) and (iter < maxiter):
        ret_tuple = fitfunc(ydata, thismask, arg_dict, **kwargs_optimizer)
        if (len(ret_tuple) == 2):
            result, ymodel = ret_tuple
            invvar_use = invvar
        elif (len(ret_tuple) == 3):
            result, ymodel, invvar_use = ret_tuple
        else:
            msgs.error('Invalid return value from fitfunc')

        thismask_iter = thismask.copy()
        thismask, qdone = pydl.djs_reject(ydata, ymodel, outmask=thismask, inmask=inmask, invvar=invvar_use,
                                          lower=lower, upper=upper, maxdev=maxdev, maxrej=maxrej,
                                          groupdim=groupdim, groupsize=groupsize, groupbadpix=groupbadpix, grow=grow,
                                          use_mad=use_mad, sticky=sticky)
        nrej = np.sum(thismask_iter & np.invert(thismask))
        nrej_tot = np.sum(inmask & np.invert(thismask))
        if verbose:
            msgs.info(
                'Iteration #{:d}: nrej={:d} new rejections, nrej_tot={:d} total rejections out of ntot={:d} '
                'total pixels'.format(iter, nrej, nrej_tot, nin_good))
        iter += 1

    if (iter == maxiter) & (maxiter != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in robust_optimize')
    outmask = np.copy(thismask)
    if np.sum(outmask) == 0:
        msgs.warn('All points were rejected!!! The fits will be zero everywhere.')

    # Perform a final fit using the final outmask if new pixels were rejected on the last iteration
    if qdone is False:
        ret_tuple = fitfunc(ydata, outmask, arg_dict, **kwargs_optimizer)

    return ret_tuple + (outmask,)

    #return result, ymodel, outmask

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

    Based on jsonify in `linetools
    <https://pypi.python.org/pypi/linetools>`_.

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
    if isinstance(obj, (np.float64, np.float32)):
        obj = float(obj)
    elif isinstance(obj, (np.int32, np.int64, np.int16)):
        obj = int(obj)
    elif isinstance(obj, np.bool_):
        obj = bool(obj)
#    elif isinstance(obj, bytes):
#        obj = obj.decode('utf-8')
    elif isinstance(obj, (np.string_, str)):
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
