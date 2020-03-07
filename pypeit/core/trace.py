"""
Module with generalized tracing routines.

Routines are primarily used for tracing slit edges.

TODO: Add object and wavelength tracing routines here?

TODO: Is there a way that we could define this link so that it's
accessible by the docstring of all modules?

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
from collections import Counter

from IPython import embed

import numpy as np
from scipy import ndimage, signal, interpolate
from matplotlib import pyplot as plt

from astropy.stats import sigma_clipped_stats, sigma_clip

from pypeit import msgs
from pypeit import utils
from pypeit import sampling
from pypeit.core import moment, pydl, arc

# TODO: Some of these functions could probably just live in pypeit.edges

def detect_slit_edges(flux, bpm=None, median_iterations=0, min_sqm=30., sobel_mode='nearest',
                      sigdetect=30., grow_bpm=5):
    r"""
    Find slit edges using the input image.

    The primary algorithm is to run a Sobel filter on the image and
    then trigger on all significant gradients. Positive gradients are
    left edges, negative gradients are right edges.

    Args:
        flux (`numpy.ndarray`_):
            Calibration frame used to identify slit edges. Likely a
            flat-field image that has been lightly smoothed in the
            spectral direction. The image should also have its bad
            pixels replaced (see
            :func:`pypeit.core.procimg.replace_columns`). The image
            *must* follow the pypeit convention, with shape
            :math:`(N_{\rm spec}, N_{\rm spat})`.
        bpm (`numpy.ndarray`_, optional):
            A boolean or integer bad-pixel mask.  If None, all pixels
            are assumed valid.  This is used to ignore features in the
            image that may be due to bad pixels.
        median_iterations (:obj:`int`, optional):
            Number of median smoothing iteration to perform on the trace
            image.  The size of the smoothing is always (7,3).  For
            long-slit data, we recommend `median_iterations=0`.
        min_sqm (:obj:`float`, optional):
            Minimum error used when detecting a slit edge.  TODO: This
            needs a better description.
        sobel_mode (:obj:`str`, optional):
            Mode to use with the Sobel filter.  See
            `scipy.ndimage.sobel`_.
        sigdetect (:obj:`float`, optional):
            Threshold for edge detection.
        grow_bpm (:int)
            The sobel_sig and edg_img are masked using the bpm. This is done by convolving the bpm with
            a spatial boxcar of width grow_bpm pixels, ensuring that pixels which touched bad pixels are also masked.

    Returns:
        Returns two `numpy.ndarray`_ objects: (1) The image of the
        significance of the edge detection in sigma and (2) the array
        isolating the slit edges. In the latter, left edges have a
        value of -1 and right edges have a value of 1.
    """
    # Checks
    if flux.ndim != 2:
        msgs.error('Trace image must be 2D.')
    if bpm is not None and bpm.shape != flux.shape:
        msgs.error('Mismatch in mask and trace image shapes.')

    # Specify how many times to repeat the median filter.  Even better
    # would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian
    # near the maximum in each column
    msgs.info("Detecting slit edges in the trace image")

    # Generate sqrt image
    sqmstrace = np.sqrt(np.abs(flux))

    # Median filter
    # TODO: Add size to parameter list
    for ii in range(median_iterations):
        sqmstrace = ndimage.median_filter(sqmstrace, size=(7, 3))

    # Replace pixel values near 0
    sqmstrace[(sqmstrace < 1.0) & (sqmstrace >= 0.0)] = 1.0
    sqmstrace[(sqmstrace > -1.0) & (sqmstrace <= 0.0)] = -1.0

    # Filter with a Sobel
    filt = ndimage.sobel(sqmstrace, axis=1, mode=sobel_mode)
    # Apply the bad-pixel mask
    if bpm is not None:
        # NOTE: Casts to float because filt is float
        filt *= (1.0 - bpm)
    # Significance of the edge detection
    sobel_sig = np.sign(filt)*np.power(filt,2)/np.maximum(sqmstrace, min_sqm)

    # First edges assigned according to S/N
    # TODO: why not match the sign of the Sobel image to the edge it
    # traces? I.e., why is the sign flipped?
    tedges = np.zeros(flux.shape, dtype=np.float)
    tedges[np.where(sobel_sig > sigdetect)] = -1.0  # A positive gradient is a left edge
    tedges[np.where(sobel_sig < -sigdetect)] = 1.0  # A negative gradient is a right edge
    
    # Clean the edges
    wcl = np.where((ndimage.maximum_filter1d(sobel_sig, 10, axis=1) == sobel_sig) & (tedges == -1))
    wcr = np.where((ndimage.minimum_filter1d(sobel_sig, 10, axis=1) == sobel_sig) & (tedges == 1))
    edge_img = np.zeros(sobel_sig.shape, dtype=np.int)
    edge_img[wcl] = -1
    edge_img[wcr] = 1

    if bpm is not None:
        msgs.info("Applying bad pixel mask")
        # JFH grow the bad pixel mask in the spatial direction
        _nave = np.fmin(grow_bpm, flux.shape[0])
        # Construct the kernel for mean calculation
        kernel = np.ones((1, _nave)) / float(_nave)
        bpm_grow = ndimage.convolve(bpm.astype(float), kernel, mode='nearest') > 0.0
        edge_img *= (1 - bpm_grow)
        sobel_sig *= (1 - bpm_grow)

    return sobel_sig, edge_img


def identify_traces(edge_img, max_spatial_separation=4, follow_span=10, minimum_spec_length=50):
    r"""
    Follow slit edges to identify unique slit traces.

    Args:
        edge_img (`numpy.ndarray`_):
            An array marked with -1 for left slit edges and +1 for
            right slit edges and 0 everywhere else. The image *must*
            follow the pypeit convention, with shape :math:`(N_{\rm
            spec},N_{\rm spat})`. See :func:`detect_slit_edges`.
        max_spatial_separation (:obj:`int`, optional):
            The maximum spatial separation between two edges in proximal
            spectral rows before they become separated into different
            slit traces.
        follow_span (:obj:`int`, optional):
            The number of previous spectral rows to consider when
            following slits forward.
        minimum_spec_length (:obj:`int`, optional):
            The minimum number of spectral rows in an edge trace.
            Traces that do not meet this criterion are ignored.

    Returns:
        `numpy.ndarray`_: An integer array with trace ID numbers at
        pixels locating the edge of that trace. Negative traces are
        for left edges, positive for right edges. The number of left
        and right traces can be determined using
        :func:`count_edge_traces`. Pixels not associated to any edge
        have a value of 0.
    """
    msgs.info('Finding unique traces among detected edges.')
    # Check the input
    if edge_img.ndim > 2:
        msgs.error('Provided edge image must be 2D.')
    if not np.all(np.isin(np.unique(edge_img), [-1,0,1])):
        msgs.error('Edge image must only have -1, 0, or 1 values.')

    # No edges were detected.
    if np.all(edge_img == 0):
        msgs.warn('No edges were found!')
        return np.zeros_like(edge_img, dtype=int)

    # Find the left and right coordinates
    lx, ly = np.where(edge_img == -1)
    rx, ry = np.where(edge_img == 1)
    x = np.concatenate((lx, rx))
    # Put left traces at negative y
    y = np.concatenate((-ly, ry))

    # The trace ID to associate with each coordinate
    traceid = np.full_like(x, -1)

    # Loop over spectral channels
    last = 0
    for row in range(np.amin(x), np.amax(x)+1):
        # Find the slit edges in this row
        indx = x == row
        in_row = np.sum(indx)
        if in_row == 0:
            # No slits found in this row
            continue

        # Find the unique edge y positions in the selected set of
        # previous rows and their trace IDs
        prev_indx = np.logical_and(x < row, x > row - follow_span)
        if not np.any(prev_indx):
            # This is likely the first row or the first row with any
            # slit edges
            traceid[indx] = np.arange(in_row)+last
            last += in_row
            continue
        uniq_y, uniq_i = np.unique(y[prev_indx], return_index=True)
        uniq_t = traceid[prev_indx][uniq_i]

        # Assign trace IDs to this row
        #   - First match to any previous IDs
        row_trace = np.full(in_row, -1)
        for i, _y in enumerate(y[indx]):
            dist = np.absolute(uniq_y-_y)
            mindist = np.argmin(dist)
            if dist[mindist] < max_spatial_separation:
                row_trace[i] = uniq_t[mindist]
        #   - Assign new trace IDs to unmatched edges
        unassigned = row_trace == -1
        n_unassigned = np.sum(unassigned)
        row_trace[unassigned] = np.arange(n_unassigned)+last
        last += n_unassigned
        #   - Assign all edges and continue
        traceid[indx] = row_trace

    # Reorder the traces and remove any that do not meet the specified
    # length.
    # TODO: This duplicates functionality in EdgeTraceSet. Rethink
    # this.
    #   - Left edges.  Given negative IDs starting with -1
    indx = y < 0
    left, reconstruct, counts = np.unique(traceid[indx], return_inverse=True, return_counts=True)
#    if np.any(counts > edge_img.shape[0]):
#        warnings.warn('Some traces have more pixels than allowed by the image.  The maximum '
#                      'spatial separation for the edges in a given trace may be too large.')
    good_trace = counts > minimum_spec_length
    left[:] = 0
    left[good_trace] = -1-np.arange(np.sum(good_trace))
    traceid[indx] = left[reconstruct]

    #   - Right edges.  Given positive IDs starting with 1
    indx = np.invert(indx)
    right, reconstruct, counts = np.unique(traceid[indx], return_inverse=True, return_counts=True)
#    if np.any(counts > edge_img.shape[0]):
#        warnings.warn('Some traces have more pixels than allowed by the image.  The maximum '
#                      'spatial separation for the edges in a given trace may be too large.')
    good_trace = counts > minimum_spec_length
    right[:] = 0
    right[good_trace] = 1+np.arange(np.sum(good_trace))
    traceid[indx] = right[reconstruct]

    # Construct the image with the trace IDs and return
    traceid_img = np.zeros_like(edge_img, dtype=int)
    traceid_img[x,np.absolute(y)] = traceid
    return traceid_img


def count_edge_traces(edge_img):
    """
    Count the number of left and right edges traced.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
    
    Returns:
        Two integers with the number of left and right edges,
        respectively.
    """
    # Avoid returning -0
    nleft = np.amin(edge_img)
    return 0 if nleft == 0 else -nleft, np.amax(edge_img)


def atleast_one_edge(edge_img, bpm=None, flux_valid=True, buffer=0, copy=False):
    """
    Ensure that there is at least one left and one right slit edge
    identified.

    This is especially useful for long slits that fill the full
    detector, e.g. Shane Kast.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
        bpm (`numpy.ndarray`_, optional):
            Integer (0 unmasked; 1 masked) or boolean array indicating
            bad pixels in the image.  If None, all pixels are considered
            good.
        flux_valid (:obj:`bool`, optional):
            The flux in the image used to construct the edge traces is
            valid meaning that any problems should not be an issue with
            the trace image itself.
        buffer (:obj:`int`, optional):
            If adding an edge, this is the minimum number of pixels
            near the detector edge at which to place the edge.
        copy (:obj:`bool`, optional):
            Copy `edge_img` to a new array before making any
            modifications.  Otherwise, `edge_img` is modified in-place.
   
    Returns:
        `numpy.ndarray`_: The modified trace image, which is either a
        new array or points to the in-place modification of `edge_img`
        according to the value of `copy`.  If no slit edges were found
        and the flux in the trace image is invalid (`flux_valid=False`),
        function returns `None`.
    """
    # Get the number of traces
    nleft, nright = count_edge_traces(edge_img)

    # Determine whether or not to edit the image in place
    _edge_img = edge_img.copy() if copy else edge_img

    if nleft != 0 and nright != 0:
        # Don't need to add anything
        return _edge_img

    if nleft == 0 and nright == 0 and not flux_valid:
        # No traces and fluxes are invalid.
        # TODO: This used to just be a warning, but I'm having it stop
        # the code if no traces are found and the flux is low.
        msgs.error('Unable to trace any edges!  Image flux is low; check trace image is correct.')

    # Use the mask to determine the first and last valid pixel column
    sum_bpm = np.zeros(edge_img.shape[1]) if bpm is None else np.sum(bpm, axis=0) 

    if nleft == 0:
        # Add a left edge trace at the first valid column
        msgs.warn('No left edge found. Adding one at the detector edge.')
        gdi0 = np.min(np.where(sum_bpm[buffer:] == 0)[0]) + buffer
        _edge_img[:,gdi0] = -1

    if nright == 0:
        # Add a right edge trace at the last valid column
        msgs.warn('No right edge found. Adding one at the detector edge.')
        gdi1 = np.max(np.where(sum_bpm[:-buffer] == 0)[0])
        _edge_img[:,gdi1] = 1

    return _edge_img


# TODO: This needs to be better tested
def handle_orphan_edges(edge_img, sobel_sig, bpm=None, flux_valid=True, buffer=0, copy=False):
    """
    In the case of single left/right traces and multiple matching
    traces, pick the most significant matching trace and remove the
    others.

    If *no* left and/or right edge is present, this will add one using
    :func:`atleast_one_edge`.

    Args:
        edge_img (`numpy.ndarray`_):
            Image with edge trace pixels numbered by their associated
            trace.  Pixels with positive numbers follow right slit edges
            and negative numbers follow left slit edges.
        sobel_sig (`numpy.ndarray`_):
            Image with the significance of the edge detection.  See
            :func:`detect_slit_edges`.
        bpm (`numpy.ndarray`_, optional):
            Integer (0 unmasked; 1 masked) or boolean array indicating
            bad pixels in the image.  If None, all pixels are considered
            good.
        flux_valid (:obj:`bool`, optional):
            The flux in the image used to construct the edge traces is
            valid meaning that any problems should not be an issue with
            the trace image itself.
        buffer (:obj:`int`, optional):
            If adding an edge, this is the minimum number of pixels
            near the detector edge at which to place the edge; see
            :func:`atleast_one_edge`.
        copy (:obj:`bool`, optional):
            Copy `edge_img` to a new array before making any
            modifications. Otherwise, `edge_img` is modified
            in-place.

    Returns:
        `numpy.ndarray`_: The modified trace image, which is either a
        new array or points to the in-place modification of
        `edge_img` according to the value of `copy`.
    """
    # Get the number of traces
    nleft, nright = count_edge_traces(edge_img)

#    if nleft == 0 or nright == 0:
#        # Deal with no left or right edges
#        # TODO: I think we should skip this and handle it in sync()
#        _edge_img = atleast_one_edge(edge_img, bpm=bpm, flux_valid=flux_valid, buffer=buffer,
#                                     copy=copy)
#        # Update the number of edges
#        nleft, nright = count_edge_traces(_edge_img)
#    else:
#        # Just do basic setup
#        _edge_img = edge_img.copy() if copy else edge_img

    _edge_img = edge_img.copy() if copy else edge_img
    
#    if nleft != 1 and nright != 1 or nleft == 1 and nright == 1:
    if nleft == 0 or nright == 0 or nleft != 1 and nright != 1 or nleft == 1 and nright == 1:
        # Nothing to do
        return _edge_img

    if nright > 1:
        # To get here, nleft must be 1.  This is mainly in here for
        # LRISb, which is a real pain..
        msgs.warn('Only one left edge, and multiple right edges.')
        msgs.info('Restricting right edge detection to the most significantly detected edge.')
        # Find the most significant right trace
        best_trace = np.argmin([-np.median(sobel_sig[_edge_img==t]) for t in range(nright)])+1
        # Remove the other right traces
        indx = _edge_img == best_trace
        _edge_img[(_edge_img > 0) & np.invert(indx)] = 0
        # Reset the number to a single right trace
        _edge_img[indx] = 1
        return _edge_img

    # To get here, nright must be 1.
    msgs.warn('Only one right edge, and multiple left edges.')
    msgs.info('Restricting left edge detection to the most significantly detected edge.')
    # Find the most significant left trace
    best_trace = np.argmax([np.median(sobel_sig[_edge_img == -t]) for t in range(nleft)])+1
    # Remove the other left traces
    indx = _edge_img == best_trace
    _edge_img[(_edge_img > 0) & np.invert(indx)] = 0
    # Reset the number to a single left trace
    _edge_img[indx] = 1

    return _edge_img


def most_common_trace_row(trace_bpm, valid_frac=1/3.):
    ## JFH DO not use row and column in the docs!!!! Change everywhere to spectral spatial. Traces always
    ## run spectral
    """
    Find the spectral position (row) that crosses the most traces.

    If provided the mask for a single trace, this just returns the
    median of the unmasked rows.

    Args:
        trace_bpm (`numpy.ndarray`_):
            Bad-pixel mask for the trace data (True is bad; False is
            good). Can be a 1D array for a single trace or a 2D array
            with shape (nspec, ntrace) for multiple traces.
        valid_frac (:obj:`float`, optional):
            The valid fraction of the detector from which to choose
            the row. For example, if 1/3, only choose the row from
            the central third of the rows.

    Returns:
        :obj:`int`: The row that crosses the most valid trace data.
    """
    if trace_bpm.ndim == 1 or trace_bpm.shape[1] == 1:
        # Only a single vector provided. Use the central valid pixel
        rows = np.where(np.invert(np.squeeze(trace_bpm)))[0]
        return rows[rows.size//2]

    s,e = ((0.5 + np.array([-1,1])*valid_frac/2)*trace_bpm.shape[0]).astype(int)
    gpm = np.invert(trace_bpm[s:e,:])
    n_good = np.sum(gpm, axis=0)
    if np.all(n_good == n_good[0]):
        # All the traces have the same number of unmasked pixels at any
        # spectral row, so use the central row
        return trace_bpm.shape[0]//2

    # Return the row with the most unmasked trace positions
    return Counter(np.where(gpm)[0]).most_common(1)[0][0] + s


def prepare_sobel_for_trace(sobel_sig, bpm=None, boxcar=5, side='left'):
    """
    Prepare the Sobel filtered image for tracing.

    This method:
        - Flips and/or truncates the pixel values based on the edge
          side to be traced (see `side`).
        - Smooths along rows (spatially)

    Args:           
        sobel_sig (`numpy.ndarray`_):
            Image with the significance of the edge detection; see
            :func:`detect_slit_edges`.
        boxcar (:obj:`int`, optional):
            Boxcar smooth the detection image along rows before
            recentering the edge centers; see
            :func:`pypeit.utils.boxcar_smooth_rows`. If `boxcar` is
            less than 1, no smoothing is performed.
        side (:obj:`str`, optional):
            The side that the image will be used to trace. In the
            Sobel image, positive values are for left traces,
            negative for right traces. If 'left', the image is
            clipped at a minimum value of -0.1. If 'right', the image
            sign is flipped and then clipped at a minimum of -0.1. If
            None, the image is not flipped or clipped, only smoothed.

    Returns:
        `numpy.ndarray`_: The smoothed image.
    """
    # NOTE: This performs the operations of what was previously
    # performed on the object passed to trace_crude_init, as well as
    # the smoothing done at the beginning of that function
    if side not in ['left', 'right', None]:
        raise ValueError('Side must be left, right, or None.')
    # Keep both sides if the side is undefined.
    # TODO: This 0.1 is drawn out of the ether and different from what is done in peak_trace
    img = sobel_sig if side is None else np.maximum((1 if side == 'left' else -1)*sobel_sig, 0.1)
    # Returned the smoothed image
    wgt = None if bpm is None else np.invert(bpm)
    return utils.boxcar_smooth_rows(img, boxcar, wgt=wgt, replace='zero')


def follow_centroid(flux, start_row, start_cen, ivar=None, bpm=None, fwgt=None, width=6.0,
                    maxshift_start=0.5, maxshift_follow=0.15, maxerror=0.2, continuous=True,
                    bitmask=None):
    """
    Follow the centroid of features in an image along the first axis.

    In the normal pypeit usage, this follows the spatial position
    (column) of a feature in the image as a function of spectral
    position (row).

    Starting from a specified row and input centers along each
    column, attempt to follow a set of features to both lower and
    higher rows in the provided image.

    Importantly, this function does not treat each row independently
    (as would be the case for a direct call to
    :func:`masked_centroid` for centroid positions along all rows),
    but treats the calculation of the centroids sequentially where
    the result for each row is dependent on and starts from the
    result from the previous row. The only independent measurement is
    the one performed at the input `start_row`. This function is much
    slower than :func:`masked_centroid` because of this introduced
    dependency.

    .. note::
        - This is an adaptation of ``trace_crude`` from ``idlspec2d``.
        - You should consider smoothing the input ``flux`` array
          before passing it to this function. See
          :func:`pypeit.utils.boxcar_smooth_rows`. For example::

            smimg = utils.boxcar_smooth_rows(img, nave, wgt=inmask)
            cen, cene, cenm = trace.follow_centroid(smimg, ...)

    Args:
        flux (`numpy.ndarray`_):
            Image used to weight the column coordinates when
            recentering. For example, when tracing slit edges, this
            should typically be the Sobel-filtered trace image after
            adjusting for the correct side and performing any
            smoothing; see :func:`prepare_sobel_for_trace`. In any
            case, consider that this image may need to be smoothed
            for robust output from this function. See
            :func:`pypeit.utils.boxcar_smooth_rows`.
        start_row (:obj:`int`):
            Row at which to start the calculation. The function
            begins with this row and then continues first to higher
            indices and then to lower indices.
        start_cen (:obj:`int`, `numpy.ndarray`_, optional):
            One or more coordinates to recenter. If an array, must be
            1D.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance in the image. If not provided, unity
            variance is assumed. Used for the calculation of the
            errors in the moment analysis. If this is not provided,
            be careful with the value set for `maxerror` (see below).
        bpm (`numpy.ndarray`_, optional):
            A boolean bad-pixel mask used to ignore pixels in the
            image (bad pixels are True).
        fwgt (`numpy.ndarray`_, optional):
            An additional weight to apply to each pixel in `flux`. If
            None, weights are uniform.
        width (:obj:`float`, `numpy.ndarray`_, optional):
            The size of the window about the provided starting center
            for the moment integration window. See
            :func:`pypeit.core.moment.moment1d`.
        maxshift_start (:obj:`float`, optional):
            Maximum shift in pixels allowed for the adjustment of the
            first row analyzed.
        maxshift_follow (:obj:`float`, optional):
            Maximum shift in pixels between centroids in adjacent
            rows as the routine follows the feature away from the
            first row analyzed.
        maxerror (:obj:`float`, optional):
             Maximum allowed error in the centroid returned by
            :func:`pypeit.core.moment.moment1d`.
        continuous (:obj:`bool`, optional):
            Keep only the continuous part of the feature trace from
            the starting row.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to flag the feature traces. If None,
            assessments use a boolean array to flag traces. If not
            None, errors will be raised if the object cannot
            interpret the correct flag names defined. In addition to
            flags used by :func:`_recenter_trace_row`, this function
            uses the DISCONTINUOUS flag.

    Returns:
        Three numpy arrays are returned: the optimized center, an
        estimate of the error, and a bad-pixel mask (masked values
        are True).
    """
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    # Shape of the image with pixel weights
    nr, nc = flux.shape

    # Instantiate theses supplementary arrays here to speed things up
    # in iterative calling of moment1d. moment1d will check the array
    # sizes. TODO: moment1d no longer instantiates these if they're not
    # provided, so this likely isn't necessary; testing required.
    _ivar = np.ones_like(flux, dtype=float) if ivar is None else ivar
    _bpm = np.zeros_like(flux, dtype=bool) if bpm is None else bpm
    _fwgt = np.ones_like(flux, dtype=float) if fwgt is None else fwgt

    # Number of starting coordinates
    _cen = np.atleast_1d(start_cen)
    nt = _cen.size
    # Check coordinates are within the image
    if np.any((_cen > nc-1) | (_cen < 0)) or start_row < 0 or start_row > nr-1:
        raise ValueError('Starting coordinates incompatible with input image!')
    # Check the dimensionality
    if _cen.ndim != 1:
        raise ValueError('Input coordinates to be at most 1D.')

    # Instantiate output; just repeat input for all image rows.
    xc = np.tile(_cen, (nr,1)).astype(float)
    xe = np.zeros_like(xc, dtype=float)
    xm = np.zeros_like(xc, dtype=bool) if bitmask is None \
                else np.zeros_like(xc, dtype=bitmask.minimum_dtype())

    # NOTE: This is effectively the old trace_crude_init

    # Recenter the starting row
    i = start_row
    xc[i,:], xe[i,:], xm[i,:] = masked_centroid(flux, xc[i,:], width, ivar=_ivar, bpm=_bpm,
                                                fwgt=_fwgt, row=i, maxshift=maxshift_start,
                                                maxerror=maxerror, bitmask=bitmask, fill='bound')

    # Go to higher indices using the result from the previous row
    for i in range(start_row+1,nr):
        xc[i,:], xe[i,:], xm[i,:] = masked_centroid(flux, xc[i-1,:], width, ivar=_ivar, bpm=_bpm,
                                                    fwgt=_fwgt, row=i, maxshift=maxshift_follow,
                                                    maxerror=maxerror, bitmask=bitmask,
                                                    fill='bound')

    # Go to lower indices using the result from the previous row
    for i in range(start_row-1,-1,-1):
        xc[i,:], xe[i,:], xm[i,:] = masked_centroid(flux, xc[i+1,:], width, ivar=_ivar, bpm=_bpm,
                                                    fwgt=_fwgt, row=i, maxshift=maxshift_follow,
                                                    maxerror=maxerror, bitmask=bitmask,
                                                    fill='bound')

    # NOTE: In edgearr_tcrude, skip_bad (roughly opposite of continuous
    # here) was True by default, meaning continuous would be False by
    # default.
    if not continuous:
        # Not removing discontinuous traces
        return xc, xe, xm

    # Keep only the continuous part of the trace starting from the
    # initial row
    # TODO: Instead keep the longest continuous segment?
    # TODO: Convert continuous from T/F to a tolerance that, e.g.,
    # allows for few-pixel discontinuities, but allows the trace to
    # pick back up again
    bad = xm > 0
    p = np.arange(nr)
    for i in range(nt):
        indx = bad[:,i] & (p > start_row)
        if np.any(indx):
            s = np.amin(p[indx])
            xm[s:,i] = True if bitmask is None else bitmask.turn_on(xm[s:,i], 'DISCONTINUOUS')
        indx = bad[:,i] & (p < start_row)
        if np.any(indx):
            e = np.amax(p[indx])-1
            xm[:e,i] = True if bitmask is None else bitmask.turn_on(xm[:e,i], 'DISCONTINUOUS')

    # Return centers, errors, and mask
    return xc, xe, xm


def masked_centroid(flux, cen, width, ivar=None, bpm=None, fwgt=None, row=None,
                    weighting='uniform', maxshift=None, maxerror=None, bitmask=None, fill='input',
                    fill_error=-1):
    """
    Measure the centroid within 1D apertures and flag and fill bad
    results.

    This is primarily a wrapper for
    :func:`pypeit.core.moment.moment1d` that flags the output.

    Data are flagged for the following reasons:
        - The centroids must always be within the aperture used by
          :func:`pypeit.core.moment.moment1d` (`OUTSIDEAPERTURE`).
        - The centroid cannot be within half an aperture of the image
          edge (`EDGEBUFFER`).
        - If `maxerror` is not None, the centroid error must be below
          the provided maximum (`MOMENTERROR`).
        - If `maxshift` is not None, the centroid must be different
          from the input value by less than the maximum
          (`LARGESHIFT`).

    The flags are provided as either a boolean array or an array of
    mask bits, depending on the value of `bitmask` (see below). The
    bitmask flags used for each criteria are given in the list above.
    If `fill` is 'input', the output centroids are replaced by the
    input value for any measurements captured by these flags. For
    centroids above `maxshift`, however, the replacement value can be
    at the maximum shift if `fill='bound'`.

    Args:
        flux (`numpy.ndarray`_):
            Array used for the centroid calculations.
        cen (`numpy.ndarray`_):
            Current estimate of the trace center.
        width (:obj:`float`):
            Passed directly to :func:`pypeit.core.moment.moment1d`;
            see the documentation there.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance in `flux`; passed directly to
            :func:`pypeit.core.moment.moment1d`.
        bpm (`numpy.ndarray`_, optional):
            Boolean bad-pixel mask for `flux`; passed directly to
            :func:`pypeit.core.moment.moment1d`.
        fwgt (`numpy.ndarray`_, optional):
            A weight to apply to each pixel in `flux`; passed
            directly to :func:`pypeit.core.moment.moment1d`.
        row (:obj:`int`, optional):
            Row (index along the first axis; spectral position) in
            `flux` at which to recenter the trace position. See
            `row` in :func:`pypeit.core.moment.moment1d`.
        weighting (:obj:`str`, optional):
            Passed directly to :func:`pypeit.core.moment.moment1d`;
            see the documentation there.
        maxshift (:obj:`float`, optional):
            Maximum shift allowed between the input and recalculated
            centroid.  If None, no limit is applied.
        maxerror (:obj:`float`, optional):
            Maximum error allowed in the calculated centroid.
            Measurements with errors larger than this value are
            returned as the input center value. If None, no limit is
            applied.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to toggle the returned bit masks. If
            provided, must be able to interpret MATHERROR,
            OUTSIDEAPERTURE, EDGEBUFFER, MOMENTERROR, and LARGESHIFT
            flags. If None, the function returns boolean flags set to
            True if there was an error in
            :func:`pypeit.core.moment.moment1d` or if the error is
            larger than `maxerror` (and `maxerror` is not None);
            centroids that have been altered by the maximum shift are
            *not* flagged.
        fill (:obj:`str`, optional):
            A string keyword specifying how flagged centroids should
            be replaced. Options are:

                - 'input': Replace the flagged centroids with their
                  input value.
                - 'bound': *Only* for the case where the
                  centroid is outside the provided `maxshift`,
                  replace the centroid by the value at either the
                  positive or negative boundary edge. In all other
                  cases, the centroid is still replaced with the
                  input value. See replacement cases above.

        fill_error (:obj:`float`, optional):
            For flagged centroids, this error is replaced with this
            dummy value.

    Returns:
        Returns three `numpy.ndarray`_ objects: the new centers, the
        center errors, and the measurement flags with a data type
        depending on `bitmask`.
    """
    # Calculate the moments
    radius = width/2
    xfit, xerr, matherr = moment.moment1d(flux, cen, width, ivar=ivar, bpm=bpm, fwgt=fwgt,
                                          row=row, weighting=weighting, order=1,
                                          fill_error=fill_error)

    # Flag centroids outide the aperture and too close to the image edge
    outside_ap = (np.absolute(xfit - cen) > radius + 0.5)
    edge_buffer = (xfit < radius - 0.5) | (xfit > flux.shape[1] - 0.5 - radius)
    indx = matherr | outside_ap | edge_buffer
    xfit[indx] = cen[indx]
    xerr[indx] = fill_error
    if maxshift is None and maxerror is None and bitmask is None:
        # Nothing else to do
        ## TODO: JFH It seems like this shold be returning indx here and not matherr if bitmask is None.
        return xfit, xerr, indx

    # Flag large shifts
    if maxshift is not None:
        large_shift = np.absolute(xfit-cen) > maxshift
        if fill == 'bound':
            xfit = np.clip(xfit - cen, -maxshift, maxshift)+cen
        else:
            indx |= large_shift

    # Flag large errors
    if maxerror is not None:
        large_error = xerr > maxerror
        indx |= large_error

    # Replace flagged data with the input and dummy error
    xfit[indx] = cen[indx]
    xerr[indx] = fill_error

    # Toggle the mask bits
    if bitmask is not None:
        xmsk = np.zeros_like(xfit, dtype=bitmask.minimum_dtype())
        xmsk[matherr] = bitmask.turn_on(xmsk[matherr], 'MATHERROR')
        xmsk[outside_ap] = bitmask.turn_on(xmsk[outside_ap], 'OUTSIDEAPERTURE')
        xmsk[edge_buffer] = bitmask.turn_on(xmsk[edge_buffer], 'EDGEBUFFER')
        if maxerror is not None:
            xmsk[large_error] = bitmask.turn_on(xmsk[large_error], 'MOMENTERROR')
        if maxshift is not None:
            xmsk[large_shift] = bitmask.turn_on(xmsk[large_shift], 'LARGESHIFT')

    # Return the new centers, errors, and flags
    return xfit, xerr, indx if bitmask is None else xmsk

# NOTE: keck_run_july changes: maxdev changed from 5.0 to 2.0
def fit_trace(flux, trace_cen, order, ivar=None, bpm=None, trace_bpm=None, weighting='uniform',
              fwhm=3.0, maxshift=None, maxerror=None, function='legendre', maxdev=2.0, maxiter=25,
              niter=9, bitmask=None, debug=False, idx=None, xmin=None, xmax=None):
    """
    Iteratively fit the trace of a feature in the provided image.

    Each iteration performs two steps:
        - Remeasure the trace centroid using :func:`masked_centroid`.
          The size of the integration window (see the definition of
          the `width` parameter for
          :func:`pypeit.core.moment.moment1d`) depends on the type of
          weighting: For *uniform weighting*, the code does a third
          of the iterations with `width = 2*1.3*fwhm`, a third with
          `width = 2*1.1*fhwm`, and a third with `width = 2*fwhm`.
          For *Gaussian weighting*, all iterations use `width =
          fwhm/2.3548`.
        - Fit the centroid measurements with a 1D function of the
          provided order. See :func:`pypeit.core.pydl.TraceSet`.

    The number of iterations performed is set by the keyword argument
    `niter`. There is no convergence test, meaning that this number
    of iterations is *always* performed.

    Notes:
        Revision History:
            - 23-June-2018  Written by J. Hennawi

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing. Must be 2D with shape (nspec,
            nspat).
        trace_cen (`numpy.ndarray`_):
            Initial guesses for trace centroids in the spatial
            dimension. This can either be an 2-d array with shape
            (nspec, nTrace) array, or a 1-d array with shape (nspec)
            for the case of a single trace.
        order (:obj:`int`):
            Order of function to fit to each trace.  See `function`.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance of the image intensity.  If not provided,
            unity variance is used.  If provided, must have the same
            shape as `flux`.
        bpm (`numpy.ndarray`_, optional):
            Boolean array with the input bad-pixel mask for the image
            (pixels to ignore are True). If not provided, all values
            in `flux` are considered valid. If provided, must have
            the same shape as `flux`.
        trace_bpm (`numpy.ndarray`_, optional):
            Boolean array with the trace bad-pixel mask; i.e., places
            where you know the trace is going to be bad that you
            always want to mask in the fits. Shape must match
            `trace_cen`.
        weighting (:obj:`str`, optional):
            The weighting to apply to the position within each
            integration window (see options in
            :func:`pypeit.core.moment.moment1d`).
        fwhm (:obj:`float`, optional):
            The expected width of the feature to trace, which is used
            to define the size of the integration window during the
            centroid calculation; see description above.
        maxshift (:obj:`float`, optional):
            Maximum shift allowed between the input and recalculated
            centroid (see :func:`masked_centroid`). If None, no limit
            is applied.
        maxerror (:obj:`float`, optional):
            Maximum error allowed in the calculated centroid (see
            :func:`masked_centroid`). Measurements with errors larger
            than this value are returned as the input center value.
            If None, no limit is applied.
        function (:obj:`str`, optional):
            The name of the function to fit. Must be a valid
            selection. See :class`pypeit.core.pydl.TraceSet`.
        maxdev (:obj:`float`, optional):
            If provided, reject points with `abs(data-model) >
            maxdev` during the fitting. If None, no points are
            rejected. See :func:`pypeit.utils.robust_polyfit_djs`.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allowed during the
            fitting. See :func:`pypeit.utils.robust_polyfit_djs`.
        niter (:obj:`int`, optional):
            The number of iterations for this method; i.e., the
            number of times the two-step fitting algorithm described
            above is performed.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to toggle the returned bit masks in edge
            centroid measurements; see :func:`masked_centroid`.
        debug (:obj:`bool`, optional):
            Plot the data and the fits.
        idx (`numpy.ndarray`_, optional):
            Array of strings with the IDs for each object. Used only
            if `debug` is true for the plotting. Default is just a
            running number.
        xmin (:obj:`float`, optional):
            Lower reference for robust_polyfit polynomial fitting.
            Default is to use zero
        xmax (:obj:`float`, optional):
            Upper reference for robust_polyfit polynomial fitting.
            Default is to use the image size in nspec direction

    Returns:
        tuple: Returns four `numpy.ndarray`_ objects all with the same
        shape as the input positions (`trace_cen`) and provide:

            - The best-fitting positions of each trace determined by the
              polynomial fit.
            - The centroids of the trace determined by either flux- or
              Gaussian-weighting, to which the polynomial is fit.
            - The errors in the centroids.
            - Boolean flags for each centroid measurement (see
              :func:`pypeit.core.moment.moment1d`).

    """
    # Ensure setup is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    if ivar is None:
        ivar = np.ones_like(flux, dtype=float)
    if ivar.shape != flux.shape:
        raise ValueError('Inverse variance array shape is incorrect.')
    if bpm is None:
        bpm = np.zeros_like(flux, dtype=bool)
    if bpm.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')
    fwgt = np.ones_like(flux, dtype=float)
    if trace_bpm is None:
        trace_bpm = np.zeros_like(trace_cen, dtype=bool)
    if trace_bpm.shape != trace_cen.shape:
        raise ValueError('Trace mask array shape is incorrect.')

    # Allow for single vectors as input
    _trace_cen = trace_cen.reshape(-1,1) if trace_cen.ndim == 1 else trace_cen
    _trace_bpm = trace_bpm.reshape(-1, 1) if trace_cen.ndim == 1 else trace_bpm
    nspec, ntrace = _trace_cen.shape
    if _trace_cen.shape != _trace_bpm.shape:
        raise ValueError('Trace data and its bad-pixel mask do not have the same shape.')

    # Define the fitting limits
    if xmin is None:
        xmin = 0.0
    if xmax is None:
        xmax = float(nspec-1)

    # Setup the width to use for each iteration depending on the weighting used
    width = np.full(niter, 2*fwhm if weighting == 'uniform' else fwhm/2.3548, dtype=float)
    if weighting == 'uniform':
        width[:niter//3] *= 1.3
        width[niter//3:2*niter//3] *= 1.1

    # Abscissa for fitting; needs to be float type when passed to
    # TraceSet
    trace_coo = np.tile(np.arange(nspec), (ntrace,1)).astype(float)
    # Values to fit
    trace_fit = np.copy(_trace_cen)

    # Uniform weighting during the fit
    trace_fit_ivar = np.ones_like(trace_fit)
    # NOTE: keck_run_july changes: Added down-weighting of masked parts
    # of the trace.
    # TODO: This feels arbitrary
    trace_fit_ivar[_trace_bpm] = 0.1

    for i in range(niter):
        # First recenter the trace using the previous trace fit/data.
        # See the replacement rules of masked_centroid for how
        # measurements that hit boundaries are treated. This replaces
        # any measurement that is outside the integration window,
        # within a buffer of the detector edge, above a maximum shift,
        # or is above the maximum error. In these cases, the output is
        # the same as the input and bad_trace has the relevant bit or
        # boolean.
        cen, err, msk = masked_centroid(flux, trace_fit, width[i], ivar=ivar, bpm=bpm, fwgt=fwgt,
                                        weighting=weighting, maxshift=maxshift, maxerror=maxerror,
                                        bitmask=bitmask)

        ################################################################
        # NOTE: keck_run_july changes: Now always replace the masked
        # data with the input trace data; downweight the masked points;
        # and include the masked points in the fit.

        # Always set the masked values to the initial input trace data.
        # The input trace data initially comes from either the standard
        # or the slit boundaries, and if we continually replace it for
        # all iterations we will naturally always extraplate the trace
        # to match the shape of a high S/N ratio fit (i.e. either the
        # standard or the flat which was used to determine the slit
        # edges.
        cen[_trace_bpm] = _trace_cen[_trace_bpm]

        # Do not do any kind of masking based on the trace recentering
        # errors. Trace fitting is much more robust when masked pixels
        # are simply replaced by the input trace values. Therefore,
        # masked pixels are not excluded from the fit, but we do give
        # them lower weight via that inverse variance. See the
        # instantation of trace_fit_ivar above.
        ################################################################

        # Fit the data
        traceset = pydl.TraceSet(trace_coo, cen.T,
                                 # Removed by keck_run_july:  inmask=np.invert(_trace_bpm.T),
                                 function=function, ncoeff=order, maxdev=maxdev, maxiter=maxiter,
                                 invvar=trace_fit_ivar.T, xmin=xmin, xmax=xmax)

        # TODO: Keep this around for now. I wanted to see how each
        # iteration affected the centroids and fit.
#        if debug:
#            bad = msk.astype(bool)
#            good = np.invert(bad)
#            for i in range(trace_fit.shape[1]):
#                plt.scatter(trace_coo[i,:], trace_fit[:,i], color='0.7', marker='.', s=50, lw=0,
#                            label='input')
#                plt.scatter(trace_coo[i,good[:,i]], trace_cen[good[:,i],i],
#                            color='k', marker='.', s=50, lw=0, label='output')
#                plt.scatter(trace_coo[i,bad[:,i]], trace_cen[bad[:,i],i],
#                            color='C1', marker='x', s=20, lw=0.5, label='bad output')
#                plt.scatter(trace_coo[i,:], traceset.yfit[i,:], color='r', marker='.', s=50, lw=0,
#                            label='fit')
#            plt.show()

        # TODO: Report iteration number and mean/stddev in difference
        # of coefficients with respect to previous iteration

        # TODO: Do this (as was done before)? This means in the next
        # iteration, the values being fit are based on the results of
        # this fit even for the bad traces instead of the original
        # input data.
        trace_fit = traceset.yfit.T

    # Plot the final fit if requested
    if debug:
        # Set the title based on the type of weighting used
        title_text = 'Flux Weighted' if weighting == 'uniform' else 'Gaussian Weighted'
        if idx is None:
            idx = np.arange(1,ntrace+1).astype(str)

        # Construct boolean flags
        inpgpm = np.invert(_trace_bpm)
        cengpm = np.invert(msk.astype(bool))
        fitgpm = traceset.outmask.T
        bpm_fit = _trace_bpm & fitgpm
        bpm_rej = _trace_bpm & np.invert(fitgpm)
        gpm_bdcen_fit = inpgpm & np.invert(cengpm) & fitgpm
        gpm_bdcen_rej = inpgpm & np.invert(cengpm) & np.invert(fitgpm)
        gpm_gdcen_fit = inpgpm & cengpm & fitgpm
        gpm_gdcen_rej = inpgpm & cengpm & np.invert(fitgpm)

        for i in range(ntrace):
            # Plot data masked on input and included in fit using input
            # locations and lower weight
            if np.any(bpm_fit[:,i]):
                plt.scatter(trace_coo[i,bpm_fit[:,i]], cen[bpm_fit[:,i],i], marker='o',
                            color='0.3', s=30, label='Input masked, fit')

            # Plot data masked on input and included in fit using input
            # locations and lower weight, but rejected by the fit
            if np.any(bpm_rej[:,i]):
                plt.scatter(trace_coo[i,bpm_rej[:,i]], cen[bpm_rej[:,i],i], marker='x',
                            color='C6', s=30, label='Input masked, fit, rejected')

            # Plot data with bad recentroid measurements, included in
            # fit using input locations and lower weight
            if np.any(gpm_bdcen_fit[:,i]):
                plt.scatter(trace_coo[i,gpm_bdcen_fit[:,i]], cen[gpm_bdcen_fit[:,i],i], marker='o',
                            color='0.7', s=30, label='Centroid masked, fit')

            # Plot data with bad recentroid measurements, included in
            # fit using input locations and lower weight, but rejected
            # by the fit
            if np.any(gpm_bdcen_rej[:,i]):
                plt.scatter(trace_coo[i,gpm_bdcen_rej[:,i]], cen[gpm_bdcen_rej[:,i],i], marker='x',
                            color='C1', s=30, label='Centroid masked, fit, rejected')

            # Plot data with good recentroid measurements and included
            # in fit
            if np.any(gpm_gdcen_fit[:,i]):
                plt.scatter(trace_coo[i,gpm_gdcen_fit[:,i]], cen[gpm_gdcen_fit[:,i],i], marker='o',
                            color='k', s=30, label='Remeasured and fit')

            # Plot data with good recentroid measurements and included
            # in fit but rejected
            if np.any(gpm_gdcen_rej[:,i]):
                plt.scatter(trace_coo[i,gpm_gdcen_rej[:,i]], cen[gpm_gdcen_rej[:,i],i], marker='x',
                            color='C3', s=30, label='Remeasured, fit, and rejected')

            # Plot all input trace locations as a line
            plt.plot(trace_coo[i,:], _trace_cen[:,i], color='C2', linewidth=1.5,
                     linestyle='--', label='Input Trace Data')

            # Plot all input trace locations as a line
            plt.plot(trace_coo[i,:], trace_fit[:,i], color='r', linewidth=2.0,
                     linestyle='--', label='Fit')

            plt.title(title_text + ' Centroid fit for trace {0}.'.format(idx[i]))
            plt.ylim((0.995*np.amin(trace_fit[:,i]), 1.005*np.amax(trace_fit[:,i])))
            plt.xlabel('Spectral Pixel')
            plt.ylabel('Spatial Pixel')
            plt.legend()
            plt.show()

    # Returns the fit, the actual weighted traces and errors, and
    # measurement flags for the last iteration
    return trace_fit, cen, err, msk, traceset


def build_trace_bpm(flux, trace_cen, bpm=None, boxcar=None, thresh=None, median_kernel=None):
    """
    Construct a bad-pixel mask for edge trace data.

    If no keyword arguments are provided, the traces are only masked
    when they land outside the bounds of the image.

    If both `boxcar` and `thresh` are provided, traces are also
    masked by extracting the provided image along the trace (see
    :func:`pypeit.core.moment.moment1d`) and flagging extracted
    values below the provided threshold.

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing. Shape is expected to be (nspec,
            nspat).
        trace_cen (`numpy.ndarray`_):
            Trace locations. Can be a 1D array for a single trace or a
            2D array with shape (nspec, ntrace) for multiple traces.
        bpm (`numpy.ndarray`_, optional):
            Boolean array with the input bad-pixel mask for the
            image. If not provided, all values in `flux` are
            considered valid. If provided, must have the same shape
            as `flux`.
        boxcar (:obj:`float`, optional):
            The width of the extraction window used for all traces
            and spectral rows. If None, the trace mask will not
            consider the extracted flux.
        thresh (:obj:`float`, optional):
            The minimum valid value of the extraced flux used to mask
            the traces. If None, the trace mask will not consider the
            extracted flux.
        median_kernel (:obj:`int`, optional):
            The spectral width of the kernel to use with
            `scipy.signal.medfilt` to filter the *extracted* data
            before setting the trace mask based on the provided
            threshold. If None, the extracted data are not filtered
            before flagging data below the threshold.

    Returns:
        `numpy.ndarray`_: The boolean mask for the traces.
    """
    # Setup and ensure input is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    nspec, nspat = flux.shape
    if bpm is None:
        bpm = np.zeros_like(flux, dtype=bool)
    if bpm.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')
    _trace_cen = trace_cen.reshape(-1,1) if trace_cen.ndim == 1 else trace_cen
    if _trace_cen.shape[0] != nspec:
        raise ValueError('Must provide trace position for each spectral pixel.')

    # Flag based on the trace positions
    trace_bpm = (_trace_cen < 0) | (_trace_cen > nspat - 1)

    if boxcar is None or thresh is None:
        # Only flagging based on the trace positions
        return trace_bpm

    # Get the extracted flux
    extract_flux = moment.moment1d(flux, _trace_cen, boxcar, bpm=bpm)[0]
    if median_kernel is not None:
        # Median filter the extracted data
        extract_flux = signal.medfilt(extract_flux, kernel_size=(median_kernel,1))
    return trace_bpm | (extract_flux < thresh)



# TODO: Add an option where the user specifies the number of slits, and
# so it takes only the highest peaks from detect_lines
def peak_trace(flux, ivar=None, bpm=None, trace_map=None, extract_width=None, smash_range=None,
               peak_thresh=100.0, peak_clip=None, trough=False, trace_median_frac=0.01,
               trace_thresh=10.0, fwhm_uniform=3.0, fwhm_gaussian=3.0, maxshift=None,
               maxerror=None, function='legendre', order=5, maxdev=5.0, maxiter=25,
               niter_uniform=9, niter_gaussian=6, bitmask=None, debug=False):
    """
    Find and trace features in an image by identifying peaks/troughs
    after collapsing along the spectral axis.

    The image is either compressed directly or after rectification
    using the supplied `trace_map`. The provided trace data *must*
    have the same shape as the input `flux` image and map each
    spatial position as a function of spectral position. This can be
    the output of :func:`pypeit.core.pca.pca_predict` where the
    provided coordinates are `np.arange(flux.shape[1])`; see also
    :func:`pypeit.edges.EdgeTracePCA.predict`. The rectification of
    the input `flux` is done using a boxcar extraction along the
    provided traces with a width of `extract_width`; see
    :func:`pypeit.core.moment.moment1d`.

    The (rectified) image is collapsed spectrally (see `smash_range`)
    giving the sigma-clipped mean flux as a function of spatial
    position. Peaks are then isolated in this vector (see
    :func:`pypeit.core.arc.detect_lines`).

    Traces that pass through these peak positions are then passed to
    two iterations of :func:`fit_trace`, which both remeasures the
    centroids of the trace and fits a polynomial to those trace data.
    The first iteration determines the centroids with uniform
    weighting, passing `fwhm=fwhm_uniform` to :func:`fit_trace`, and
    the second uses Gaussian weighting for the centroid measurements
    (passing `fwhm=fwhm_gaussian` to :func:`fit_trace`). The results
    of this second iteration of :func:`fit_trace` are the data
    returned.

    Troughs in the image can also be traced, which is done by
    flipping the sign of the image about its median and then
    repeating the "peak" finding and :func:`fit_trace` iterations. If
    troughs are fit, the traces are order with the set of peak traces
    first (the number of which is given by the last returned object
    of the function), followed by the trough traces.

    Args:
        flux (`numpy.ndarray`_):
            Image to use for tracing.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance of the image intensity.  If not provided,
            unity variance is used.  If provided, must have the same
            shape as `flux`.
        bpm (`numpy.ndarray`_, optional):
            Boolean array with the input bad-pixel mask for the
            image. If not provided, all values in `flux` are
            considered valid. If provided, must have the same shape
            as `flux`.
        trace_map (`numpy.ndarray`_, optional):
            Trace data that maps the spatial position of all spectra
            as a function of spectral row. For example, this can be
            the output of :func:`pypeit.core.pca.pca_predict` where
            the provided coordinates are `np.arange(flux.shape[1])`;
            see also :func:`pypeit.edges.EdgeTracePCA.predict`. This
            is used to rectify the input image so that spectra are
            identically organized along image rows (i.e., to select
            the `i`th spectrum, one would slice with `[:,i]`). Shape
            *must* be identical to `flux`. If None, `flux` is assumed
            to be rectified on input.
        extract_width (:obj:`float`, optional):
            The width of the extract aperture to use when rectifying
            the flux image. If None, set to `fwhm_gaussian`.
        smash_range (:obj:`tuple`, optional):
            Spectral range to over which to collapse the input image
            into a 1D flux as a function of spatial position. This 1D
            vector is used to detect features for tracing. This is
            useful (and recommended) for definining the relevant
            detector range for data with spectra that do not span the
            length of the detector. The tuple gives the minimum and
            maximum in the fraction of the full spectral length
            (nspec). If None, the full image is collapsed.
        peak_thresh (:obj:`float, optional):
            The threshold for detecting peaks in the image. See the
            `input_thresh` parameter for
            :func:`pypeit.core.arc.detect_lines`.
        peak_clip (:obj:`float, optional):
            Sigma-clipping threshold used to clip peaks with small
            values; no large values are clipped. If None, no clipping
            is performed. Generally, one should instead raise the
            value of `peak_thresh` instead, if the peak detection
            algorithm is finding insignificant peaks.
        trough (:obj:`bool`, optional):
            Trace both peaks **and** troughs in the input image. This
            is done by flipping the value of the smashed image about
            its median value, such that troughs can be identified as
            peaks.
        trace_median_frac (:obj:`float`, optional):
            After rectification of the image and before refitting the
            traces, the rectified image is median filtered with a
            kernel width of trace_median_frac*nspec along the
            spectral dimension.
        trace_thresh (:obj:`float`, optional):
            After rectification and median filtering of the image
            (see `trace_median_frac`), values in the resulting image
            that are *below* this threshold are masked in the
            refitting of the trace using :func:`fit_trace`.
        fwhm_uniform (:obj:`float`, optional):
            The `fwhm` parameter to use when using uniform weighting
            in the calls to :func:`fit_trace`. See description of the
            algorithm above.
        fwhm_gaussian (:obj:`float`, optional):
            The `fwhm` parameter to use when using Gaussian weighting
            in the calls to :func:`fit_trace`. See description of the
            algorithm above.
        maxshift (:obj:`float`, optional):
            Maximum shift allowed between the input and recalculated
            centroid (see :func:`fit_trace`).
        maxerror (:obj:`float`, optional):
            Maximum error allowed in the calculated centroid (see
            :func:`fit_trace`).
        function (:obj:`str`, optional):
            The type of polynomial to fit to the trace data. See
            :func:`fit_trace`.
        order (:obj:`int`, optional):
            Order of the polynomial to fit to each trace. See
            :func:`fit_trace`.
        maxdev (:obj:`float`, optional):
            See :func:`fit_trace`. If provided, reject points with
            `abs(data-model) > maxdev` when fitting the trace. If
            None, no points are rejected.
        maxiter (:obj:`int`, optional):
            Maximum number of rejection iterations allowed during the
            fitting. See :func:`fit_trace`.
        niter_uniform (:obj:`int`, optional):
            The number of iterations for :func:`fit_trace` when edge
            measurements are based on uniform weighting. See
            description above.
        niter_gaussian (:obj:`int`, optional):
            The number of iterations for :func:`fit_trace` when edge
            measurements are based on Gaussian weighting. See
            description above.
        bitmask (:class:`pypeit.bitmask.BitMask`, optional):
            Object used to toggle the returned bit masks in edge
            centroid measurements; see :func:`masked_centroid`.
        debug (:obj:`bool`, optional):
            Show plots useful for debugging.

    Returns:
        Returns four `numpy.ndarray`_ objects and the number of peak
        traces. The number of peak traces should be used to separate
        peak from trough traces; if `trough` is False, this will just be
        the total number of traces.  The four `numpy.ndarray`_ objects
        provide:

            - The best-fitting positions of each trace determined by the
              polynomial fit.
            - The centroids of the trace determined by the
              Gaussian-weighting iteration, to which the polynomial is
              fit.
            - The errors in the Gaussian-weighted centroids.
            - Boolean flags for each centroid measurement (see
              :func:`pypeit.core.moment.moment1d`).

    """
    # Setup and ensure input is correct
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    nspec, nspat = flux.shape
    if ivar is None:
        ivar = np.ones_like(flux, dtype=float)
    if ivar.shape != flux.shape:
        raise ValueError('Inverse variance array shape is incorrect.')
    if bpm is None:
        bpm = np.zeros_like(flux, dtype=bool)
    if bpm.shape != flux.shape:
        raise ValueError('Mask array shape is incorrect.')

    # Define the region to collapse
    if smash_range is None:
        smash_range = (0.,1.)

    # Set the image to collapse
    if trace_map is None:
        # Assume image already rectified
        flux_extract = flux
        # Just set the trace to the follow the spatial columns
        trace_map = np.tile(np.arange(nspat), (nspec,1))
    else:
        # Check there is a trace for each image pixel
        if trace_map.shape != flux.shape:
            raise ValueError('Provided trace data must match the image shape.')
        msgs.info('Rectifying image by extracting along trace for each spatial pixel')
        # TODO: JFH What should this aperture size be? I think fwhm=3.0
        # since that is the width of the sobel filter
        flux_extract = sampling.rectify_image(flux, trace_map, bpm=bpm, extract_width=fwhm_gaussian 
                                                if extract_width is None else extract_width)[0]
#        if debug:
#            ginga.show_image(flux_extract, chname ='rectified image')

    # Collapse the image along the spectral direction to isolate peaks/troughs
    start, end = np.clip(np.asarray(smash_range)*nspec, 0, nspec).astype(int)
    msgs.info('Collapsing image spectrally between pixels {0}:{1}'.format(start, end))
    flux_smash_mean, flux_smash_median, flux_smash_sig \
            = sigma_clipped_stats(flux_extract[start:end,:], axis=0, sigma=4.0)

    # Offset by the median
    # TODO: If tracing Sobel-filtered image, this should be close to,
    # or identically, 0
    flux_median = np.median(flux_smash_mean)
    flux_smash_mean -= flux_median

    # Trace peak or both peaks and troughs
    label = ['peak', 'trough'] if trough else ['peak']
    sign = [1, -1] if trough else [1]

    # Instantiate output
    npeak = 0
    fit = np.empty((nspec,0), dtype=float)
    cen = np.empty((nspec,0), dtype=float)
    err = np.empty((nspec,0), dtype=float)
    msk = np.empty((nspec,0), dtype=bool if bitmask is None else bitmask.minimum_dtype())

    # Get the smoothing kernel width and ensure it is odd
    median_kernel = None if trace_median_frac is None \
                        else int(np.ceil(nspec*trace_median_frac))//2 * 2 + 1

    # Identify and trace features in the image
    for i,(l,s) in enumerate(zip(label,sign)):

        # Identify the peaks
        msgs.info('Searching for peaks.')
        peak, _, _cen, _, _, best, _, _ \
                = arc.detect_lines(s*flux_smash_mean, cont_subtract=False, fwhm=fwhm_gaussian,
                                   input_thresh=peak_thresh, max_frac_fwhm=4.0,
                                   min_pkdist_frac_fwhm=5.0, debug=debug)

        if len(_cen) == 0 or not np.any(best):
            msgs.warn('No good {0}s found!'.format(l))
            continue
        msgs.info('Found {0} good {1}(s) in the rectified, collapsed image'.format(
                    len(_cen[best]),l))

        # Set the reference spatial locations to use for tracing the
        # detected peaks
        # TODO: Added this for a test case, not sure we should keep it
        # in the long run.
        _cen = _cen[best]
        loc = np.round(_cen).astype(int) 
        if peak_clip is not None:
            # Clip the peaks based on their amplitude as a stop-gap for
            # a detection threshold that may be too low. Only clip
            # aberrantly low values.
            clipped_peak = sigma_clip(peak[best], sigma_lower=peak_clip, sigma_higher=np.inf)
            peak_mask = np.ma.getmaskarray(clipped_peak)
            if np.any(peak_mask):
                msgs.warn('Clipping {0} detected peak(s) with aberrant amplitude(s).'.format(
                                np.sum(peak_mask)))
                loc = loc[np.invert(peak_mask)]
                _cen = _cen[np.invert(peak_mask)]

        # As the starting point for the iterative trace fitting, use
        # the input trace data at the positions of the detected peaks.
        # The trace at column `loc` is expected to pass through `loc`
        # at the reference spectral pixel. The offset of the trace is
        # to allow for non-integer measurements of the peak centroid.
        trace_peak = trace_map[:,loc] + (_cen-loc)[None,:]

        # Image to trace: flip when tracing the troughs and clip low
        # values
        # TODO: This 1 is drawn out of the ether; this is different
        # from what is done in prepare_sobel_for_trace
        _flux = np.clip(s*(flux - flux_median), 1, None)

        # Construct the trace mask
        trace_peak_bpm = np.zeros(trace_peak.shape, dtype=bool) if trace_thresh is None \
                            else build_trace_bpm(_flux, trace_peak, bpm=bpm, boxcar=fwhm_gaussian,
                                                 thresh=trace_thresh, median_kernel=median_kernel)

        # Remeasure and fit the trace using uniform weighting
        trace_peak, _cen, _err, _msk, _ \
                = fit_trace(_flux, trace_peak, order, ivar=ivar, bpm=bpm,
                            trace_bpm=trace_peak_bpm, fwhm=fwhm_uniform, maxshift=maxshift,
                            maxerror=maxerror, function=function, maxdev=maxdev, maxiter=maxiter,
                            niter=niter_uniform, bitmask=bitmask, debug=debug)

        # Reset the mask
        # TODO: Use or include `bad` resulting from fit_trace()?
        trace_peak_bpm = np.zeros(trace_peak.shape, dtype=bool) if trace_thresh is None \
                            else build_trace_bpm(_flux, trace_peak, bpm=bpm, boxcar=fwhm_gaussian,
                                                 thresh=trace_thresh, median_kernel=median_kernel)

        # Redo the measurements and trace fitting with Gaussian
        # weighting
        trace_peak, _cen, _err, _msk, _ \
                = fit_trace(_flux, trace_peak, order, ivar=ivar, bpm=bpm, trace_bpm=trace_peak_bpm,
                            weighting='gaussian', fwhm=fwhm_gaussian, maxshift=maxshift,
                            maxerror=maxerror, function=function, maxdev=maxdev, maxiter=maxiter,
                            niter=niter_gaussian, bitmask=bitmask, debug=debug)

        # Save the results
        fit = np.append(fit, trace_peak, axis=1)
        cen = np.append(cen, _cen, axis=1)
        err = np.append(err, _err, axis=1)
        msk = np.append(msk, _msk, axis=1)

        if i == 0:
            # Save the number of peaks (troughs are appended, if they're located)
            npeak = cen.shape[1]

    return fit, cen, err, msk, npeak


def parse_user_slits(add_slits, this_det, rm=False):
    """
    Parse the parset syntax for adding slits

    Args:
        add_slits (str, list):
          Taken from the parset
        this_det (int):
          current detector
        rm (bool, optional):
          Remove instead of add?

    Returns:
        list or None:
          if list,  [[x0,x1,yrow]] for add with one or more entries
          if list,  [[xcen,yrow]] for rm with one or more entries

    """
    # Might not be a list yet (only a str)
    if not isinstance(add_slits, list):
        add_slits = [add_slits]
    #
    user_slits = []
    for islit in add_slits:
        if not rm:
            det, x0, x1, yrow = [int(ii) for ii in islit.split(':')]
            if det == this_det:
                user_slits.append([x0,x1,yrow])
        else:
            det, xcen, yrow = [int(ii) for ii in islit.split(':')]
            if det == this_det:
                user_slits.append([xcen,yrow])
    # Finish
    if len(user_slits) == 0:
        return None
    else:
        return user_slits
