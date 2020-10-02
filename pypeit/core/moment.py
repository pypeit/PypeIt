"""
Module to compute moments.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

import numpy as np
from scipy import special


def moment1d(flux, col, width, ivar=None, bpm=None, fwgt=None, row=None, weighting='uniform',
             order=0, bounds=None, fill_error=-1., mesh=False):
    r"""
    Compute one-dimensional moments of the provided image within an
    aperture along its second axis (axis=1).

    This method allows for computations of the zeroth, first, and
    second moments (see `order`). The aperture used for the
    calculation is centered at the provided `col` pixel with a width
    defined by `width`; however, this definition depends on the type
    of weighting applied (see `weighting`). Formulae for each moment
    are as follows. The zeroth moment (`order=0`) computes the
    discrete weighted sum of the flux within the aperture:
    
    .. math::
        \mu_0 &= \sum_i w_i f_i \\
        \epsilon_{\mu_0}^2 &= \sum_i (w_i \epsilon_{f,i})^2,

    where :math:`f` is the flux in each pixel :math:`i`,
    :math:`\epsilon_f` is its error, and :math:`w` is the assigned
    weight (see `weighting`).  The first moment (`order=1`) computes the
    flux-weighted center of window:

    .. math::
        \mu_1 &= \frac{\sum_i x_i w_i f_i }{\mu_0} \\
        \epsilon_{\mu_1}^2 &= \mu_0^{-2}\ \sum_i [ w_i \epsilon_{f,i}\
                            (x_i - \mu_1)]^2,

    where :math:`x` is the pixel position along the 2nd axis (see
    `col`).  The second moment (`order=2`) computes the variance of the
    flux profile about its center within the window:
    
    .. math::
        \mu_2 &= \frac{\sum_i x^2_i w_i f_i }{\mu_0} - \mu_1^2 \\
        \epsilon_{\mu_2}^2 &= \mu_0^{-2}\ \sum_i w^2_i \epsilon^2_{f,i}\
                                [(x_i - \mu_1)^2 + \mu_2]^2.

    The values returned for the second-moment calculation are actually
    the standard deviation instead of the variance, where:

    .. math::
        \sigma = \sqrt{\mu_2} \\
        \epsilon_{\sigma} = \frac{\epsilon_{\mu_2}}{2 \sigma}.

    The method uses `numpy.ma.MaskedArray` objects to keep track of
    math errors, such as divisions by 0. The returned boolean array
    indicates when these errors occurred, and the method replaces
    these errors with the original centers.

    The shape of arrays depend on the shapes of the input image
    (`flux`), columns about which to perform the calculation (`col`),
    and rows from which to select the column data (`row`) as listed
    below.  In the following, the input flux array must be 2D and has
    shape :math:`(N_{\rm row}, N_{\rm col})`, whereas the number of rows
    and columns identified for a moment calculation are :math:`N_{\rm
    mom,row}` and :math:`N_{\rm mom,col}`, respectively.  The possible
    shapes of the output arrays are:

        - If `row` is None and `col` is 1D, the moments are calculated
          centered at each `col` for all rows in `flux`. The shape of
          the array for each moment is then :math:`(N_{\rm row}, N_{\rm
          mom,col})`
        - If `row` is None and `col` is 2D, `col` must have shape
          :math:`(N_{\rm row}, N_{\rm mom,col})` and the shape of the
          output array per moment calculated is the same.
        - If `row` is a single integer and `col` is 1D, the moments are
          calculated for each col only at the provided row in `flux`.
          The shape of the array for each moment is then :math:`(N_{\rm
          mom,col},)`.
        - If `row` is a 1D integer array and `col` is 1D, the default
          behavior is for the moments to be calculated at each col
          for each provided row. The shape of the array for each
          moment would then be :math:`(N_{\rm mom,row},N_{\rm
          mom,col})`. This will be true even if :math:`N_{\rm
          mom,row} = N_{\rm mom,col}`, as long as `mesh` is `True`.
          However, if :math:`N_{\rm mom,row} = N_{\rm mom,col}` and
          `mesh` is `False`, only the matched pairs of `col` and
          `row` will be used and the output shape per moment will be
          :math:`(N_{\rm mom,row},)`.
        - If `row` is a 1D integer array and `col` is 2D, the first axis
          of `col` must be the same length as `row`. The shape of the
          array for each moment is then :math:`(N_{\rm mom,row},N_{\rm
          mom,col})`, which is the same as the input `col`.
        - If `row` is 2D and `col` is 2D, they must have the same shape.
        - If `row` is a single integer and `col` is 2D, or if `row` is
          2D and `col` is 1D, the method raises an error.

    .. note::

        - This is an entirely general function, as reflected by the
          nomenclature used in the call. As used within PypeIt, the
          PypeIt image orientation convention means that moments are
          always taken along the spatial direction; i.e., `col` is
          the spatial coordinate and `row` is the spectral
          coordinate.

        - This function is a generalization of and builds on the
          heritage of functions in idlspec2d, specifically
          trace_fweight, trace_gweight, extrace_asymbox2,
          extract_boxcar.

    .. warning::

        The function has significant setup/input checking. Most of
        this introduces limited overhead with the exception of the
        handling of `ivar`, `bpm`, and `fwgt`. If any of these are
        provided as `None` on input, an array is constructed (unity
        for `ivar` and `fwgt` and all False for `bpm`) that serves as
        a place holder. If repetitive calls to the function are
        expected and any of these arrays are missing, significant
        efficiency gains can be made by providing pre-built values
        for these arrays so that time isn't lost in allocating the
        placeholder arrays in every call.

    .. todo::

        Optimize the code for efficiency, regardless of the input.

    Args:
        flux (`numpy.ndarray`_):
            Intensity image with shape :math:`(N_{\rm row}, N_{\rm
            col})`.
        col (`numpy.ndarray`_):
            Floating-point center along the 2nd axis for the integration
            window in pixel index (first value located at index 0).
            This can either be a 1D or 2D array.  See restrictions on
            the shape in the description above.
        width (:obj:`float`, `numpy.ndarray`_):
            The meaning of the parameter depends on the value of
            `weighting`.  If `weighting=='uniform'`, the width of the
            integration window in columns, centered at the input `col`.
            If `weighting=='gaussian'`, the :math:`\sigma` of a
            pixelated Gaussian weighting function.  The width of the
            integration window in columns, centered at the input `col`,
            is always `6*width` (i.e., the half-width is
            :math:`3\sigma`).  The provided value can be a scalar to use
            a constant window definition, or it can be an array variable
            integration window where the array must have the same shape
            as `col`.
        ivar (`numpy.ndarray`_, optional):
            Inverse variance of the image intensity.  If not provided,
            unity variance is used.  If provided, must have the same
            shape as `flux`.
        bpm (`numpy.ndarray`_, optional):
            Boolean bad-pixel mask for the input image. True values
            are ignored, False values are included. If not provided,
            all pixels are included. If provided, must have the same
            shape as `flux`.
        fwgt (`numpy.ndarray`_, optional):
            An additional weight to apply to each pixel in `flux`.  If
            None, weights are uniform.  Otherwise, the :math:`w_i` from
            above are the product of this weight and the result of the
            scheme set using the `weighting` argument.
        row (:obj:`int`, `numpy.ndarray`_, optional):
            Integer or integer array with the position along the first
            axis (axis=0) for the moment calculation.  This can either
            be None, an integer, or a 1D or 2D array.  See restrictions
            on the shape in the description above.
        weighting (:obj:`str`, optional):
            The weighting to apply to the position within each
            integration window (see `width` above). This must be
            (case-insensitive) either 'uniform' for uniform weighting or
            'gaussian' for weighting by a Gaussian centered at the input
            guess coordinates and integrated over the pixel width.
        order (:obj:`int`, array-like, optional):
            The order of the moment(s) to calculate. Can be a single
            integer or a list. Moments to calculate must be 0, 1, or
            2; at most order can be `[0,1,2]`. The shape of the
            output arrays depends on the number of moments
            calculated. Note that the calculation of the orders is
            necessarily sequential; i.e., setting `order=2` means
            that the zeroth and first moments have to be calculated
            anyway. The order must be provided in sorted order; i.e.,
            you cannot pass `order=[2,1]`.
        bounds (:obj:`tuple`, optional):
            A two-tuple with the lower and upper limit for each
            moment order. If None, no bounds are imposed. If not
            None, an upper and lower bound must be provided for each
            moment to compute; i.e., if more than one moment is
            computed, each element of the two-tuple must be an
            array-like object that matches the length of `order`. To
            set an upper or lower bound only, set the unbounded
            component to None. Bounds for the zeroth and second order
            moments are in an absolute sense, whereas first-order
            bounds are relative to the input `col`. Measurements that
            hit the bounds are masked; see the description of the
            returned objects. For example, to flag anything without a
            positive zeroth moment or a maximum shift from the input
            center of 1 pixel, call the method with arguments::

                order=[0,1], bounds=([0,-1], [None,1])
            
        fill_error (:obj:`float`, optional):
            Value to use as filler for undetermined moments,
            resulting from either the input bad-pixel mask or
            computational issues (division by zero, etc.; see return
            description below).
        mesh (:obj:`bool`, optional):
            If `col` and `row` are 1D vectors of the same length,
            this determines if each `col` and `row` should be paired
            (`mesh is False`) or used to construct a grid, i.e.,
            every `col` combined with every `row` (`mesh is True`).
            See the method description.

    Returns:
        Three `numpy.ndarray`_ objects are returned. If more than one
        moment order is requested, the moments are ordered along the
        first axis; e.g., if `order=[0,1]` the outputs `moment[0]` and
        `moment[1]` contain the zeroth and first moments, respectively.
        The subsequent dimensions of the output arrays are dictated by
        the input `row` and `col`; see the method description. The
        returned arrays are:

            - The moment calculated along the 2nd axis of the input
              image (axis=1).  Masked values (indicated by the third
              object returned) are 0 for the zeroth and second moments
              and equal to the input `col` value for the first moments.
            - The formal propagated error (see equations above) in the
              moments.  Errors are only meaningful if `ivar` is
              provided. Masked values (indicated by the third object
              returned) are set to `fill_error`.
            - A boolean bad-pixel mask for output data; True values
              should be ignored, False values are valid measurements.

    Raises:
        ValueError:
            Raised if input shapes are not correct or if the selected
            `weighting` is unknown.  See method description.

    Examples:

        First setup an image with some Gaussians:

        >>> from pypeit.core.moment import moment1d
        >>> import numpy as np
        >>> from scipy.special import erf
        >>> def gauss_comb():
        ...     c = [45,50,55]
        ...     img = np.zeros((len(c),100), dtype=float)
        ...     x = np.arange(100)
        ...     sig = 5.
        ...     for i,_c in enumerate(c):
        ...         img[i,:] = (erf((x-c[i]+0.5)/np.sqrt(2)/sig) 
        ...                         - erf((x-c[i]-0.5)/np.sqrt(2)/sig))/2.
        ...     return img
        ...
        >>> img = gauss_comb()

        Calculate all moments at one column and row:
        
        >>> mu, mue, flag = moment1d(img, 50, 40., row=0, order=[0,1,2])
        >>> print(mu)
        [ 0.99858297 45.02314924  4.97367636]
    
        Calculate all moments at one column for all rows:
        
        >>> moment1d(img, 50, 40., order=[0,1,2])[0]
        array([[ 0.99858297,  0.99993125,  0.99858297],
               [45.02314924, 50.        , 54.97685076],
               [ 4.97367636,  5.00545947,  4.97367636]])
    
        Calculate zeroth moments in all rows centered at column 50
        
        >>> moment1d(img, 50, 40., order=0)[0]
        array([0.99858297, 0.99993125, 0.99858297])

        Calculate zeroth moments in all rows for three column positions
        
        >>> moment1d(img, [45,50,55], 40., order=0, mesh=True)[0]
        array([[0.99993125, 0.99858297, 0.97670951],
               [0.99858297, 0.99993125, 0.99858297],
               [0.97670951, 0.99858297, 0.99993125]])

        Calculate first moments in all rows for the three column positions
        
        >>> moment1d(img, [45,50,55], 40., order=1, mesh=True)[0]
        array([[45.        , 45.02314924, 45.2814655 ],
               [49.97685076, 50.        , 50.02314924],
               [54.7185345 , 54.97685076, 55.        ]])

        Or pick the same column for all rows
        
        >>> moment1d(img, 50, 40., row=[0,1,2], order=1)[0]
        array([45.02314924, 50.        , 54.97685076])

        Or pick a column unique to each row
        
        >>> moment1d(img, [43,52,57], 40., row=[0,1,2], order=1)[0]
        array([44.99688181, 50.00311819, 55.00311819])

    """

    # TODO: Could be generalized further for higher dimensional
    # inputs...

    # TODO: Need to benchmark again after including the `bounds`
    # keyword.

    # Check mode input
    _weighting = weighting.lower()
    if _weighting not in ['uniform', 'gaussian']:
        raise ValueError('Weighting must be uniform or gaussian')

    # Check image input
    # TODO: For large images, instantiating ivar, bpm, and fwgt can be
    # a *significant* time sink if there are only a few moments
    # calculated
    if flux.ndim != 2:
        raise ValueError('Input image must be 2D.')
    nrow, ncol = flux.shape
    if ivar is not None and ivar.shape != flux.shape:
        raise ValueError('Inverse variance must have the same shape as the input image.')
    if bpm is not None and bpm.shape != flux.shape:
        raise ValueError('Pixel bad-pixel mask must have the same shape as the input image.')
    if fwgt is not None and fwgt.shape != flux.shape:
        raise ValueError('Pixel weights must have the same shape as the input image.')

    # Check moment order
    _order = np.atleast_1d(order)
    if not np.array_equal(_order, np.sort(_order)):
        # NOTE: This is rather strict, but it's needed to make sure
        # that the provided order and bounds make sense.
        raise ValueError('Order must be provided as a sorted array.')
    if _order.ndim != 1:
        raise ValueError('Order can be at most 1D.')
    if _order.size > 3:
        raise ValueError('Can return at most 3 moments.')
    if not np.array_equal(np.unique(_order), _order):
        raise ValueError('Moments must be unique!')
    if np.any((_order != 0) & (_order != 1) & (_order != 2)):
        raise ValueError('Selected moments must be either 0, 1, or 2.')
    norder = len(_order)

    # Check if the bounds are provided
    lower = np.array([None, None, None], dtype=object)
    upper = np.array([None, None, None], dtype=object)
    if bounds is not None:
        if len(bounds) != 2:
            raise ValueError('Bounds must be provided as a two-tuple.')
        _bounds = tuple([np.atleast_1d(b) for b in bounds])
        if np.any([len(b) != norder for b in _bounds]):
            raise ValueError('Number of bounds must match number of moments to calculate.')
        lower[_order] = _bounds[0]
        upper[_order] = _bounds[1]

    # Check coordinate input.  For both col and width, the atleast_1d
    # function does not result in a copy if the provided object is a
    # numpy.ndarray
    array_input = isinstance(col, np.ndarray)
    _col = np.atleast_1d(col)
    _width = np.atleast_1d(width)
    if _width.size == 1:
        _width = np.full(_col.shape, width, dtype=float)
    if _width.shape != _col.shape:
        raise ValueError('width must either be a single value or have the same shape as col.')
    if _col.ndim > 2:
        raise ValueError('Input columns for calculation must be a 1D or 2D array.')
    if row is None and _col.ndim == 2 and _col.shape[0] != nrow:
        raise ValueError('Length of first axis of col and flux must match.')
    if row is not None:
        _row = np.atleast_1d(row)
        if _row.ndim == 1 and _col.ndim == 2 and _col.shape[0] != _row.size:
            raise ValueError('First axis of col must match length of row.')
        if _row.ndim == 2 and _col.ndim == 2 and _col.shape != _row.shape:
            raise ValueError('col and row must have the same shape.')
        if len(_row) == 1 and _col.ndim == 2 or _row.ndim == 2 and _col.ndim == 1:
            raise ValueError('Invalid combination of row and col input shapes.')
    else:
        _row = np.arange(nrow)
    # Check column and row values are valid

    # TODO: This allows the window centers to be outside the image
    # range. This is okay because _col is never used when slicing the
    # image. However, _row values are, meaning that they have to be
    # within the image limits.
#    if np.any((_col < 0) | (_col >= ncol)):
#        raise ValueError('Column locations outside provided image.')
    if np.any((_row < 0) | (_row >= nrow)):
        raise ValueError('Row locations outside provided image.')

    # Fill the columns and rows so that they have matching shape
    rmrowdim = _row.size == 1
    rmcoldim = _col.size == 1 and not array_input
    if _col.ndim == 1 and _row.ndim == 1:
        if _row.size == 1 or _row.size == _col.size and not mesh:
            if _row.size == 1:
                # Repeat the row for each column
                _row = np.full(_col.size, _row[0], dtype=int)
            else:
                # mesh is False, so just match column and row numbers
                _row.reshape(1,-1)
                rmrowdim = True
            _col = _col.reshape(1,-1)
            _width = _width.reshape(1,-1)
        else:
            # Do the combinatorics for rows, columns, and integration
            # widths for the measurements
            nmomcol = _col.size
            _col = np.tile(_col, (_row.size,1))
            _width = np.tile(_width, (_row.size,1))
            _row = np.tile(_row, (nmomcol,1)).T
    elif _col.ndim == 2 and _row.ndim == 1:
        _row = np.tile(_row, (_col.shape[1],1)).T

    # Construct the shape for the output
    outshape = (norder,) + _col.shape
    if rmrowdim and rmcoldim:
        outshape = (outshape[0],)
    elif rmrowdim:
        outshape = (outshape[0],outshape[2])
    elif rmcoldim:
        outshape = (outshape[0],outshape[1])
    if norder == 1 and len(outshape) > 1:
        outshape = outshape[1:]

    # Flatten the coordinate arrays (creates a copy).
    _row = _row.flatten().astype(int)
    _col = _col.flatten().astype(float)
    _width = _width.flatten().astype(float)

    # The "radius" of the pixels to cover is either half of the
    # provided width for uniform weighting or 3*width for Gaussian
    # weighting, where width is the sigma of the Gaussian
    _radius = _width/2 if _weighting == 'uniform' else _width*3

    # Window for the integration for each coordinate. In the
    # calculation of `c`, the increase of the window size by 4 isn't
    # strictly necessary. At minimum it has to be 2, but the increase
    # by 4 ensures there are 0 pixels at either end of the integration
    # window. TODO: Should consider changing this to 2.
    i1 = np.floor(_col - _radius + 0.5).astype(int)
    i2 = np.floor(_col + _radius + 0.5).astype(int)
    c = i1[:,None]-1+np.arange(int(np.amax(np.amin(i2-i1)-1,0))+4)[None,:]
    ih = np.clip(c,0,ncol-1)

    # Set the weight over the window; masked pixels have 0 weight
    good = (c >= 0) & (c < ncol)
    if bpm is not None:
        # NOTE: `&=` doesn't work here because of the np.newaxis usage
        good = good & np.invert(bpm[_row[:,None],ih])
    if ivar is not None:
        good = good & (ivar[_row[:,None],ih] > 0)
    if _weighting == 'uniform':
        # Weight according to the fraction of each pixel within in the
        # integration window
        wt = good * np.clip(_radius[:,None] - np.abs(c - _col[:,None]) + 0.5,0,1)
    else:
        # Weight according to the integral of a Gaussian over the pixel
        coo = c - _col[:,None]
        wt = good * (special.erf((coo+0.5)/np.sqrt(2.)/_width[:,None])
                        - special.erf((coo-0.5)/np.sqrt(2.)/_width[:,None]))/2.

    # Construct the moment-independent component of the integrand and
    # the zeroth moment; the zeroth moment is always needed
    integ = flux[_row[:,None],ih] * wt
    if fwgt is not None:
        integ *= fwgt[_row[:,None],ih]
    mu = np.array([np.ma.sum(integ, axis=1), None, None], dtype=object)
    mue = np.array([None, None, None], dtype=object)
    mum = np.array([np.ma.getmaskarray(mu[0]).copy(), None, None], dtype=object)
    _var0 = np.square(wt) if ivar is None else np.ma.divide(np.square(wt), ivar[_row[:,None],ih])
    if 0 in _order:
        # Only calculate the error if the moment was requested
        mue[0] = np.ma.sqrt(np.ma.sum(_var0, axis=1))
        # Impose the boundary
        if lower[0] is not None:
            mum[0] |= mu[0] < lower[0]
        if upper[0] is not None:
            mum[0] |= mu[0] > upper[0]
    
    # Calculate the first moment if necessary
    if np.any(_order > 0):
        mu[1] = np.ma.divide(np.sum(integ*c, axis=1), mu[0])
        mum[1] = np.ma.getmaskarray(mu[1]).copy()
        if 1 in _order:
            # Only calculate the error if the moment was requested
            _var1 = np.square(wt * (c - mu[1][:,None]))
            if ivar is not None:
                _var1 = np.ma.divide(_var1, ivar[_row[:,None],ih])
            mue[1] = np.ma.divide(np.ma.sqrt(np.ma.sum(_var1, axis=1)), np.absolute(mu[0]))
            # Impose the boundary
            if lower[1] is not None:
                mum[1] |= mu[1] < _col - lower[1]
            if upper[1] is not None:
                mum[1] |= mu[1] > _col + upper[1]

    # Calculate the second moment if necessary
    if 2 in _order:
        mu[2] = np.ma.divide(np.sum(integ*np.square(c), axis=1), mu[0]) - np.square(mu[1])
        mue[2] = np.ma.divide(np.ma.sqrt(
                        np.ma.sum(_var0 * np.square(np.square(c - mu[1][:,None]) + mu[2][:,None]),
                                  axis=1)), np.absolute(mu[0]))
        mu[2] = np.ma.sqrt(mu[2])
        mue[2] = np.ma.divide(mue[2], 2*mu[2])
        mum[2] = np.ma.getmaskarray(mu[2]).copy()
        # Impose the boundary
        if lower[2] is not None:
            mum[2] |= mu[2] < lower[2]
        if upper[2] is not None:
            mum[2] |= mu[2] > upper[2]

    # Fill in the masked values
    for i in range(3):
        if mu[i] is None:
            continue
        mu[i][mum[i]] = _col[mum[i]] if i == 1 else 0.0
        mu[i] = mu[i].data
        if mue[i] is None:
            continue
        mue[i][mum[i]] = fill_error
        mue[i] = mue[i].filled(fill_error)

    # Return with the correct shape
    singlenum = outshape == (1,) and not array_input
    return (mu[_order][0][0], mue[_order][0][0], mum[_order][0][0]) if singlenum \
            else (np.concatenate(mu[_order]).reshape(outshape),
                  np.concatenate(mue[_order]).reshape(outshape),
                  np.concatenate(mum[_order]).reshape(outshape))

