# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""
Provides a set of functions to handle resampling.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import warnings
from IPython import embed
import numpy
from scipy import interpolate
import astropy.constants

from pypeit.core import moment

#from .covariance import Covariance


# TODO: spectral_coordinate_step, spectrum_velocity_scale, and
# angstroms_per_pixel should all use centers_to_borders

def spectral_coordinate_step(wave, log=False, base=10.0):
    """
    Return the *uniform* sampling step for the input wavelength
    vector.

    If the sampling is logarithmic, return the change in the
    logarithm of the wavelength; otherwise, return the linear step in
    angstroms.

    Args: 
        wave (`numpy.ndarray`_):
            Wavelength coordinates of each spectral channel in
            angstroms.
        log (:obj:`bool`, optional):
            Input spectrum has been sampled geometrically.
        base (:obj:`float`, optional):
            If sampled geometrically, the sampling is done using a
            logarithm with this base. For natural logarithm, use
            numpy.exp(1).

    Returns:
        :obj:`float`: Spectral sampling step in either angstroms
        (log=False) or the step in log(angstroms).

    Raises:
        ValueError:
            Raised if the wavelength vector is not linearly or
            log-linearly sampled.
    """
    dw = numpy.diff(numpy.log(wave))/numpy.log(base) if log else numpy.diff(wave)
    if numpy.any( numpy.absolute(numpy.diff(dw)) > 100*numpy.finfo(dw.dtype).eps):
        raise ValueError('Wavelength vector is not uniformly sampled to numerical accuracy.')
    return numpy.mean(dw)


def spectrum_velocity_scale(wave):
    """
    Determine the velocity sampling of an input wavelength vector when log sampled
    
    .. note::
        The wavelength vector is assumed to be geometrically sampled!
        However, the input units expected to be in angstroms, not, e.g.,
        log(angstrom).

    Args: 
        wave (numpy.ndarray): Wavelength coordinates of each spectral
            channel in angstroms.  It is expected that the spectrum has
            been sampled geometrically

    Returns:
        float: Velocity scale of the spectrum in km/s.

    """
    return astropy.constants.c.to('km/s').value*spectral_coordinate_step(wave, log=True,
                                                                         base=numpy.exp(1.))


def angstroms_per_pixel(wave, log=False, base=10.0, regular=True):
    """
    Return a vector with the angstroms per pixel at each channel.

    When `regular=True`, the function assumes that the wavelengths are
    either sampled linearly or geometrically.  Otherwise, it calculates
    the size of each pixel as the difference between the wavelength
    coordinates.  The first and last pixels are assumed to have a width
    as determined by assuming the coordinate is at its center.

    .. note::

        If the regular is False and log is True, the code does *not*
        assume the wavelength coordinates are at the geometric center of
        the pixel.

    Args:
        wave (`numpy.ndarray`_):
            (Geometric) centers of the spectrum pixels in angstroms.
        log (`numpy.ndarray`_, optional):
            The vector is geometrically sampled.
        base (:obj:`float`, optional):
            Base of the logarithm used in the geometric sampling.
        regular (:obj:`bool`, optional):
            The vector is regularly sampled.

    Returns:
        numpy.ndarray: The angstroms per pixel.
    """
    if regular:
        ang_per_pix = spectral_coordinate_step(wave, log=log, base=base)
        return ang_per_pix*wave*numpy.log(base) if log else numpy.repeat(ang_per_pix, len(wave))

    return numpy.diff([(3*wave[0]-wave[1])/2] + ((wave[1:] + wave[:-1])/2).tolist()
                      + [(3*wave[-1]-wave[-2])/2])


def grid_npix(rng=None, dx=None, log=False, base=10.0, default=None):
    """
    Determine the number of pixels needed for a given grid.

    Args:
        rng (array-like, optional):
            Two-element array with the
            starting and ending x coordinate of the pixel centers to
            divide into pixels of a given width.  If *log* is True, this
            must still be the linear value of the x coordinate, not
            log(x)!.
        dx (:obj:`float`, optional):
            Linear or logarithmic pixel width.
        log (:obj:`bool`, optional):
            Flag that the range should be logarithmically binned.
        base (:obj:`float`, optional):
            Base for the logarithm
        default (:obj:`int`, optional):
            Default number of pixels to use. The default is returned
            if either ``rng`` or ``dx`` are not provided.

    Returns:
        :obj:`tuple`: Returns the number of pixels to cover ``rng``
        with pixels of width ``dx`` and a two-element
        `numpy.ndarray`_ with the adjusted range such that number of
        pixels of size dx is the exact integer.

    Raises:
        ValueError:
            Raised if the range is not a two-element vector.
    """
    # If the range or sampling are not provided, the number of pixels is
    # already set
    if rng is None or dx is None:
        return default, rng
    if len(rng) != 2:
        raise ValueError('Range must be a 2-element vector.')

    _rng = numpy.atleast_1d(rng).copy()
    npix = int(numpy.floor(numpy.diff(numpy.log(_rng))[0]/numpy.log(base)/dx) + 1) if log else \
                    int(numpy.floor(numpy.diff(_rng)[0]/dx) + 1)
    _rng[1] = numpy.power(base, numpy.log(_rng[0])/numpy.log(base) + dx*(npix-1)) \
                            if log else _rng[0] + dx*(npix-1)

    # Fix for numerical precision
    if (not log and numpy.isclose(rng[1] - _rng[1], dx)) \
           or (log and numpy.isclose((numpy.log(rng[1]) - numpy.log(_rng[1]))/numpy.log(base), dx)):
        npix += 1
        _rng[1] = numpy.power(base, numpy.log(_rng[0])/numpy.log(base) + dx*(npix-1)) \
                                if log else _rng[0] + dx*(npix-1)

    return npix, _rng


def grid_borders(rng, npix, log=False, base=10.0):
    """
    Determine the borders of bin edges in a grid.

    Args:
        rng (array-like):
            Two-element array with the (geometric) centers of the
            first and last pixel in the grid.
        npix (:obj:`int`):
            Number of pixels in the grid.
        log (:obj:`bool`, optional):
            The input range is (to be) logarithmically sampled.
        base (:obj:`float`, optional):
            The base of the logarithmic sampling. Use
            ``numpy.exp(1.)`` for the natural logarithm.

    Returns:
        :obj:`tuple`: Returns a `numpy.ndarray`_ with the grid
        borders with shape ``(npix+1,)`` and the step size per grid
        point. If ``log=True``, the latter is the geometric step.
    """
    if log:
        _rng = numpy.log(rng)/numpy.log(base)
        dlogx = numpy.diff(_rng)[0]/(npix-1.)
        borders = numpy.power(base, numpy.linspace(*(_rng/dlogx + [-0.5, 0.5]), num=npix+1)*dlogx)
        return borders, dlogx
    dx = numpy.diff(rng)[0]/(npix-1.)
    borders = numpy.linspace(*(numpy.atleast_1d(rng)/dx + numpy.array([-0.5, 0.5])), num=npix+1)*dx
    return borders, dx


def grid_centers(rng, npix, log=False, base=10.0):
    """
    Determine the (geometric) center of pixels in a grid.

    Args:
        rng (array-like):
            Two-element array with the (geometric) centers of the
            first and last pixel in the grid.
        npix (:obj:`int`):
            Number of pixels in the grid.
        log (:obj:`bool`, optional):
            The input range is (to be) logarithmically sampled.
        base (:obj:`float`, optional):
            The base of the logarithmic sampling. Use
            ``numpy.exp(1.)`` for the natural logarithm.

    Returns:
        :obj:`tuple`: Returns a `numpy.ndarray`_ with the grid pixel
        (geometric) ceners with shape ``(npix,)`` and the step size
        per grid point. If ``log=True``, the latter is the geometric
        step.
    """
    if log:
        _rng = numpy.log(rng)/numpy.log(base)
        dlogx = numpy.diff(_rng)[0]/(npix-1.)
        centers = numpy.power(base, numpy.linspace(*(_rng/dlogx), num=npix)*dlogx)
        return centers, dlogx
    dx = numpy.diff(rng)[0]/(npix-1.)
    centers = numpy.linspace(*(numpy.atleast_1d(rng)/dx), num=npix)*dx
    return centers, dx


def borders_to_centers(borders, log=False):
    """
    Convert a set of bin borders to bin centers.

    Grid borders need not be regularly spaced.

    Args:
        borders (`numpy.ndarray`_):
            Borders for adjoining bins.
        log (:obj:`bool`, optional):
            Return the geometric center instead of the linear center
            of the bins.

    Returns:
        `numpy.ndarray`_: The vector of bin centers.
    """
    return numpy.sqrt(borders[:-1]*borders[1:]) if log else (borders[:-1]+borders[1:])/2.0


def centers_to_borders(x, log=False):
    """
    Convert a set of bin centers to bounding edges.

    Grid centers need not be regularly spaced. The first edge of the
    first bin and the last edge of the last bin are assumed to be
    equidistant from the center of the 2nd and penultimate bins,
    respectively.

    Args:
        x (`numpy.ndarray`_):
            Centers of adjoining bins.
        log (:obj:`bool`, optional):
            Adopt a geometric binning instead of a linear binning.

    Returns:
        `numpy.ndarray`_: The vector with the coordinates of
        adjoining bin edges.
    """
    if log:
        dx = numpy.diff(numpy.log(x))
        return numpy.exp(numpy.append(numpy.log(x[:-1]) - dx/2,
                                      numpy.log(x[-1]) + numpy.array([-1,1])*dx[-1]/2))
    dx = numpy.diff(x)
    return numpy.append(x[:-1] - dx/2, x[-1] + numpy.array([-1,1])*dx[-1]/2)


class Resample:
    r"""
    Resample regularly or irregularly sampled data to a new grid using
    integration.
    
    This is a generalization of the routine
    :func:`ppxf.ppxf_util.log_rebin` provided by Michele Cappellari in
    the pPXF package.

    The abscissa coordinates (``x``) or the pixel borders
    (``xBorders``) for the data (``y``) should be provided for
    irregularly sampled data. If the input data is linearly or
    geometrically sampled (``inLog=True``), the abscissa coordinates
    can be generated using the input range for the (geometric) center
    of each grid point. If ``x``, ``xBorders``, and ``xRange`` are
    all None, the function assumes grid coordinates of ``x =
    numpy.arange(y.shape[-1])``.

    The function resamples the data by constructing the borders of
    the output grid using the ``new*`` keywords and integrating the
    input function between those borders. The output data will be set
    to ``ext_value`` for any data beyond the abscissa limits of the
    input data.

    The data to resample (``y``) can be a 1D or 2D array; the
    abscissa coordinates must always be 1D. If ``y`` is 2D, the
    resampling is performed along the last axis (i.e., ``axis=-1``).

    The nominal assumption is that the provided function is a step
    function based on the provided input (i.e., ``step=True``). If
    the output grid is substantially finer than the input grid, the
    assumption of a step function will be very apparent. To assume
    the function is instead linearly interpolated between each
    provided point, choose ``step=False``; higher-order
    interpolations are not provided.

    If errors are provided, a nominal error propagation is performed
    to provide the errors in the resampled data.

    .. warning::

        Depending on the details of the resampling, the output errors
        are likely highly correlated.  Any later analysis of the
        resampled function should account for this.

    The covariance in the resampled pixels can be constructed by
    setting ``covar=True``; however, this is currently only supported
    when ``step=True``. If no errors are provided and ``covar=True``,
    the computed matrix is the *correlation* matrix instead of the
    *covariance* matrix. Given that the resampling is the same for all
    vectors, only one correlation matix will be calculated if no
    errors are provided, even if the input ``y`` is 2D. If the input
    data to be resampled is 2D and errors *are* provided, a
    covariance matrix is calculated for *each* vector in ``y``.
    Beware that this can be an expensive computation.

    The ``conserve`` keyword sets how the units of the input data
    should be treated. If ``conserve=False``, the input data are
    expected to be in density units (i.e., per ``x`` coordinate unit)
    such that the integral over :math:`dx` is independent of the
    units of :math:`x` (i.e., flux per unit angstrom or flux
    density). If ``conserve=True``, the value of the data is assumed
    to have been integrated over the size of each pixel (i.e., units
    of flux). If ``conserve=True``, :math:`y` is converted to units
    of per step in :math:`x` such that the integral before and after
    the resample is the same. For example, if :math:`y` is a spectrum
    in units of flux, the function first converts the units to flux
    density and then computes the integral over each new pixel to
    produce the new spectra with units of flux.

    .. todo::
        - Allow for higher order interpolations.
        - Enable covariance matrix calculations for ``step=False``.
        - Provide examples

    Args:
        y (`numpy.ndarray`_, `numpy.ma.MaskedArray`_):
            Data values to resample. The shape can be 1D or 2D. If
            1D, the shape must be :math:`(N_{\rm pix},)`; otherwise,
            it must be :math:`(N_y,N_{\rm pix})`. I.e., the length of
            the last axis must match the input coordinates.
        e (`numpy.ndarray`_, `numpy.ma.MaskedArray`_, optional):
            Errors in the data that should be resampled. The shape
            must match the input ``y`` array. These data are used to
            perform a nominal calculation of the error in the
            resampled array.
        mask (`numpy.ndarray`_, optional):
            A boolean array indicating values in ``y`` that should be
            ignored during the resampling (values to ignore have
            ``masked=True``, just like in a `numpy.ma.MaskedArray`_).
            The mask used during the resampling is the union of this
            object and the masks of ``y`` and ``e``, if either are
            provided as `numpy.ma.MaskedArray`_ objects.
        x (`numpy.ndarray`_, optional):
            Abscissa coordinates for the data, which do not need to
            be regularly sampled. If the pixel borders are not
            provided, they are assumed to be half-way between
            adjacent pixels, and the first and last borders are
            assumed to be equidistant about the provided value. If
            these coordinates are not provided, they are determined
            by the input borders, the input range, or just assumed to
            be the indices, :math:`0..N_{\rm pix}-1`.
        xRange (array-like, optional):
            A two-element array with the starting and ending value
            for the coordinates of the centers of the first and last
            pixels in ``y``. Default is :math:`[0,N_{\rm pix}-1]`.
        xBorders (`numpy.ndarray`_, optional):
            An array with the borders of each pixel that must have a
            length of :math:`N_{\rm pix}+1`.
        inLog (:obj:`bool`, optional):
            Flag that the input is logarithmically binned, primarily
            meaning that the coordinates are at the geometric center
            of each pixel and the centers are spaced logarithmically.
            If false, the sampling is expected to be linear.
        newx (array-like, optional):
            Abscissa coordinates for the *output* data, which do not
            need to be a regular grid. If this is provided, the pixel
            borders are assumed to be half-way between adjacent
            pixels, and the first and last borders are assumed to be
            equidistant about the provided value. If these
            coordinates are not provided, they are determined by the
            new range, the new number of pixels, and/or the new pixel
            width (and whether or not the new grid should be
            logarithmically binned). If this is provided,
            ``newRange``, ``newpix``, ``newLog``, and ``newdx`` are
            *all* ignored.
        newRange (array-like, optional):
            A two-element array with the (geometric) centers of the
            first and last pixel in the output vector. If not
            provided, assumed to be the same as the input range.
        newBorders (array-like, optional):
            An array with the borders of each pixel in the resampled
            vectors.
        newpix (:obj:`int`, optional): 
            Number of pixels for the output vector.  If not provided,
            assumed to be the same as the input vector.
        newLog (:obj:`bool`, optional):
            The output vector should be logarithmically binned.
        newdx (:obj:`float`, optional):
            The sampling step for the output vector. If
            `newLog=True`, this must be the change in the *logarithm*
            of :math:`x` for the output vector! If not provided, the
            sampling is set by the output range (see ``newRange``
            above) and number of pixels (see ``newpix`` above).
        base (:obj:`float`, optional):
            The base of the logarithm used for both input and output
            sampling, if specified. The default is 10; use
            ``numpy.exp(1)`` for natural logarithm.
        ext_value (:obj:`float`, optional):
            Set extrapolated values to the provided float. If set to
            None, values are just set to the linear extrapolation of
            the data beyond the provided limits; use `ext_value=None`
            with caution!
        conserve (:obj:`bool`, optional):
            Conserve the integral of the input vector.  For example, if
            the input vector is a spectrum in flux units, you should
            conserve the flux in the resampling; if the spectrum is in
            units of flux density, you do not want to conserve the
            integral.
        step (:obj:`bool`, optional):
            Treat the input function as a step function during the
            resampling integration.  If False, use a linear
            interpolation between pixel samples.
        covar (:obj:`bool`, optional):
            Calculate the covariance matrix between pixels in the
            resampled vector. Can only be used if ``step=True``. If
            no error vector is provided (``e``), the result is a
            *correlation* matrix.
    
    Attributes:
        x (`numpy.ndarray`_):
            The coordinates of the function on input.
        xborders (`numpy.ndarray`_):
            The borders of the input pixel samples.
        y (`numpy.ndarray`_):
            The function to resample.
        e (`numpy.ndarray`_):
            The 1-sigma errors in the function to resample.
        m (`numpy.ndarray`_):
            The boolean mask for the input function.
        outx (`numpy.ndarray`_):
            The coordinates of the function on output.
        outborders (`numpy.ndarray`_):
            The borders of the output pixel samples.
        outy (`numpy.ndarray`_):
            The resampled function.
        oute (`numpy.ndarray`_):
            The resampled 1-sigma errors.
        outf (`numpy.ndarray`_):
            The fraction of each output pixel that includes valid data
            from the input function.

    Raises:
        ValueError:
            Raised if more the one of ``x``, ``xRange``, or
            ``xBorders`` are provided, if more the one of ``newx``,
            ``newRange``, or ``newBorders`` are provided, if ``y`` is
            a `numpy.ndarray`_, if ``y`` is not 1D or 2D, if the
            covariance is requested but ``step`` is False, if the
            shapes of the provided errors or mask do not match ``y``,
            if there is insufficient information to construct the
            input or output grid, or if either ``xRange`` or
            ``newRange`` are not two-element arrays.
    """
#        covar (:class:`~mangadap.util.covariance.Covariance`):
#            The covariance or correlation matrices for the resampled
#            vectors.
    def __init__(self, y, e=None, mask=None, x=None, xRange=None, xBorders=None, inLog=False,
                 newx=None, newRange=None, newBorders=None, newpix=None, newLog=True, newdx=None,
                 base=10.0, ext_value=0.0, conserve=False, step=True): #, covar=False):

        # Check operation can be performed and is not ill-posed
        if numpy.sum([inp is not None for inp in [x, xRange, xBorders]]) != 1:
            raise ValueError('One and only one of the x, xRange, and xBorders arguments should be '
                             'provided.')
        if numpy.sum([inp is not None for inp in [newx, newRange, newBorders]]) != 1:
            raise ValueError('One and only one of the newx, newRange, and newBorders arguments '
                             'should be provided.')
        if not isinstance(y, numpy.ndarray):
            raise ValueError('Input vector must be a numpy.ndarray!')
        if y.ndim > 2:
            raise ValueError('Input must be a 1D or 2D array!')
#        if covar and not step:
#            raise ValueError('Covariance is currently only calculated for step resampling.')

        # Setup the data, errors, and mask
        self.y = y.filled(0.0) if isinstance(y, numpy.ma.MaskedArray) else y.copy()
        self.twod = self.y.ndim == 2
        self.e = None if e is None \
                    else e.filled(0.0) if isinstance(e, numpy.ma.MaskedArray) else e.copy()
        self.m = numpy.zeros(self.y.shape, dtype=bool) if mask is None else mask

        # Check the shapes
        if self.e is not None and self.e.shape != self.y.shape:
            raise ValueError('Error array shape mismatched!')
        if self.m.shape != self.y.shape:
            raise ValueError('Mask array shape mismatched!')

        # Get the union of all the relevant masks
        if isinstance(y, numpy.ma.MaskedArray):
            self.m |= y.mask
        if e is not None and isinstance(e, numpy.ma.MaskedArray):
            self.m |= e.mask

        # Get the input coordinates
        nx = self.y.shape[-1] if x is None and xBorders is None else None
        self.x, self.xborders = self._coordinate_grid(x=x, rng=xRange, nx=nx, borders=xBorders,
                                                      log=inLog, base=base)

        # If conserving integral, assume input is integrated over pixel
        # width and convert to a density function (divide by pixel size)
        if conserve:
            self.y /= (numpy.diff(self.xborders)[None,:] if self.twod \
                                else numpy.diff(self.xborders))

        # Get the output coordinates
        nx = self.x.size \
                if newx is None and newBorders is None and newpix is None and newdx is None \
                else newpix
        self.outx, self.outborders = self._coordinate_grid(x=newx, rng=newRange, nx=nx,
                                                           borders=newBorders, dx=newdx,
                                                           log=newLog, base=base)

#        if covar:
#            A = self._resample_step_matrix()
#            self.outy = numpy.dot(A, self.y.T).T
#            self.outf = numpy.dot(A, numpy.logical_not(self.m.T).astype(int)).T \
#                            / numpy.diff(self.outborders)[...,:]
#            if self.e is None:
#                # TODO: Repeat it N times if there are N y vectors?
#                self.covar = Covariance.from_matrix_multiplication(A, numpy.ones_like(self.x)
#                                                ).apply_new_variance(numpy.ones_like(self.outx))
#                self.oute = None
#            else:
#                if self.twod:
#                    self.covar = numpy.empty(self.y.shape[0], dtype=object)
#                    for i in range(self.y.shape[0]):
#                        self.covar[i] = Covariance.from_matrix_multiplication(A,
#                                            numpy.square(self.e[i])).full()
#                    self.covar = Covariance(self.covar, impose_triu=True)
#                    self.oute = numpy.sqrt(self.covar.variance().T)
#                else:
#                    self.covar = Covariance.from_matrix_multiplication(A, numpy.square(self.e))
#                    self.oute = numpy.sqrt(self.covar.variance())
#        else:
        # Perform the resampling
#        self.covar = None
        self.outy = self._resample_step(self.y) if step else self._resample_linear(self.y)
       
        # The mask and errors are always interpolated as a step function
        self.oute = None if self.e is None else self._resample_step(self.e, quad=True)
        
        self.outf = self._resample_step(numpy.logical_not(self.m).astype(int)) \
                        / numpy.diff(self.outborders)

        # Do not conserve the integral over the size of the pixel
        if not conserve:
            self.outy /= (numpy.diff(self.outborders)[None,:] if self.twod \
                            else numpy.diff(self.outborders))
            if self.oute is not None:
                self.oute /= (numpy.diff(self.outborders)[None,:] if self.twod \
                                    else numpy.diff(self.outborders))
#                if self.covar is not None:
#                    self.covar = self.covar.apply_new_variance(numpy.square(self.oute.T))

        # Set values for extrapolated regions
        if ext_value is not None:
            indx = (self.outborders[:-1] < self.xborders[0]) \
                        | (self.outborders[1:] > self.xborders[-1]) 
            if numpy.sum(indx) > 0:
                self.outy[...,indx] = ext_value
                self.outf[...,indx] = 0.
                if self.oute is not None:
                    self.oute[...,indx] = 0.

    @staticmethod
    def _coordinate_grid(x=None, rng=None, nx=None, dx=None, borders=None, log=False, base=10.0):
        """
        Use the provided information to construct the coordinate grid
        and the grid borders.
        """
        if x is not None and borders is not None:
            raise ValueError('Both x and borders provided.  Do not need to call _coordinate_grid, '
                             'but also _coordinate_grid does not check that x and borders are '
                             'consistenet with one another.')
        if (x is not None or borders is not None) and rng is not None:
            warnings.warn('Provided both x or borders and the range.  Ignoring range.')
        if x is None and borders is not None:
            # Use the borders to set the centers
            return borders_to_centers(borders, log=log), borders
        if x is not None and borders is None:
            # Use the centers to set the borders
            return x, centers_to_borders(x, log=log)

        # After this point, both x and borders should be None
        assert x is None and borders is None, 'Coding logic error'

        if rng is None and nx is None:
            raise ValueError('Insufficient input to construct coordinate grid.')

        if rng is None:
            # Just set the result to a uniform pixel grid
            return numpy.arange(nx, dtype=float) + 0.5, numpy.arange(nx+1, dtype=float)

        # After this point, rng cannot be None
        assert rng is not None, 'Coding logic error'

        if dx is not None and nx is not None:
            warnings.warn('Provided rng, dx, and nx, which over-specifies the grid; rng and nx '
                          'take precedence.')
        if nx is not None:
            borders = grid_borders(rng, nx, log=log, base=base)[0]
            return borders_to_centers(borders, log=log), borders

        nx, _rng = grid_npix(rng=rng, dx=dx, log=log, base=base)
        borders = grid_borders(_rng, nx, log=log, base=base)[0]
        return borders_to_centers(borders, log=log), borders

    def _resample_linear(self, v, quad=False):
        """Resample the vectors."""

        # Combine the input coordinates and the output borders
        combinedX = numpy.append(self.outborders, self.x)
        srt = numpy.argsort(combinedX)
        combinedX = combinedX[srt]

        # Get the indices where the data should be reduced
        border = numpy.ones(combinedX.size, dtype=bool)
        border[self.outborders.size:] = False
        k = numpy.arange(combinedX.size)[border[srt]]

        # Calculate the integrand
        if self.twod:
            # Linearly interpolate the input function at the output border positions
            interp = interpolate.interp1d(self.x, v, axis=-1, assume_sorted=True,
                                          fill_value='extrapolate')
            combinedY = numpy.append(interp(self.outborders), v, axis=-1)[:,srt]
            integrand = (combinedY[:,1:]+combinedY[:,:-1])*numpy.diff(combinedX)[None,:]/2.0
        else:
            # Linearly interpolate the input function at the output border positions
            interp = interpolate.interp1d(self.x, v, assume_sorted=True,
                                          fill_value='extrapolate')
            combinedY = numpy.append(interp(self.outborders), v)[srt]
            integrand = (combinedY[1:]+combinedY[:-1])*numpy.diff(combinedX)/2.0

        if quad:
            integrand = numpy.square(integrand)

        # Use reduceat to calculate the integral
        out = numpy.add.reduceat(integrand, k[:-1], axis=-1) if k[-1] == combinedX.size-1 \
                        else numpy.add.reduceat(integrand, k, axis=-1)[...,:-1]
    
        return numpy.sqrt(out) if quad else out

    def _resample_step(self, v, quad=False):
        """Resample the vectors."""

        # Convert y to a step function
        #  - repeat each element of the input vector twice
        _v = numpy.repeat(v, 2, axis=1) if self.twod else numpy.repeat(v, 2)
        #  - repeat each element of the border array twice, and remove
        #  the first and last elements
        _x = numpy.repeat(self.xborders, 2)[1:-1]

        # Combine the input coordinates and the output borders into a
        # single vector
        indx = numpy.searchsorted(_x, self.outborders)
        combinedX = numpy.insert(_x, indx, self.outborders)

        # Insert points at the borders of the output function
        v_indx = indx.copy()
        v_indx[indx >= _v.shape[-1]] = -1
        combinedY = numpy.array([ numpy.insert(__v, indx, __v[v_indx]) for __v in _v ]) \
                            if self.twod else numpy.insert(_v, indx, _v[v_indx])

        # Calculate the integrand
        integrand = combinedY[:,1:]*numpy.diff(combinedX)[None,:] if self.twod else \
                        combinedY[1:]*numpy.diff(combinedX)
        if quad:
            integrand = numpy.square(integrand)

        # Get the indices where the data should be reduced
        border = numpy.insert(numpy.zeros(_x.size, dtype=bool), indx,
                              numpy.ones(self.outborders.size, dtype=bool))
        k = numpy.arange(combinedX.size)[border]

        # Use reduceat to calculate the integral
        out = numpy.add.reduceat(integrand, k[:-1], axis=-1) if k[-1] == combinedX.size-1 \
                    else numpy.add.reduceat(integrand, k, axis=-1)[...,:-1]
        return numpy.sqrt(out) if quad else out

    def _resample_step_matrix(self):
        r"""
        Build a matrix such that

        .. math::
            y = \mathbf{A} x

        where :math:`x` is the input vector, :math:`y` is the resampled
        vector, and :math:`\mathbf{A}` is the matrix operations that
        resamples :math:`x`.
        """
        ny = self.outx.size
        nx = self.x.size
        dx = numpy.diff(self.xborders)

        # Repeat each element of the border array twice, and remove the
        # first and last elements
        _p = numpy.repeat(numpy.arange(self.x.size), 2)
        _x = numpy.repeat(self.xborders, 2)[1:-1]

        # Combine the input coordinates and the output borders into a
        # single vector
        indx = numpy.searchsorted(_x, self.outborders)
        combinedX = numpy.insert(_x, indx, self.outborders)

        # Insert points at the borders of the output function
        p_indx = indx.copy()
        p_indx[indx >= _p.shape[-1]] = -1
        combinedP = numpy.insert(_p, indx, _p[p_indx])

        # Get the indices where the data should be reduced
        border = numpy.insert(numpy.zeros(_x.size, dtype=bool), indx,
                              numpy.ones(self.outborders.size, dtype=bool))
        nn = numpy.where(numpy.logical_not(border))[0][::2]
        k = numpy.zeros(len(combinedX), dtype=int)
        k[border] = numpy.arange(numpy.sum(border))
        k[nn-1] = k[nn-2]
        k[nn] = k[nn-1]
        start,end = numpy.where(border)[0][[0,-1]]

        # Calculate the fraction of each pixel into each output pixel
        fraction = numpy.diff(combinedX[start:end+1])
        # Construct the output matrix
        indx = fraction > 0
        A = numpy.zeros((ny, nx), dtype=float)
        A[k[start:end][indx], combinedP[start:end][indx]] = fraction[indx]
        return A


def rectify_image(img, col, bpm=None, ocol=None, max_ocol=None, extract_width=None,
                  mask_threshold=0.5):
    r"""
    Rectify the image by shuffling flux along columns using the provided
    column mapping.

    The image recification is one dimensional, treating each image row
    independently. It can be done either by a direct resampling of the
    image columns using the provided mapping of output to input column
    location (see `col` and :class:`Resample`) or by an extraction along
    the provided column locations (see `extract_width`). The latter is
    generally faster; however, when resampling each row, the flux is
    explicitly conserved (see the `conserve` argument of
    :class:`Resample`).

    Args:
        img (`numpy.ndarray`_):
            The 2D image to rectify. Shape is :math:`(N_{\rm row},
            N_{\rm col})`.
        col (`numpy.ndarray`_):
            The array mapping each output column to its location in
            the input image. That is, e.g., `col[:,0]` provides the
            column coordinate in `img` that should be rectified to
            column 0 in the output image. Shape is :math:`(N_{\rm
            row}, N_{\rm map})`.
        bpm (`numpy.ndarray`_, optional):
            Boolean bad-pixel mask for pixels to ignore in input
            image. If None, no pixels are masked in the
            rectification. If provided, shape must match `img`.
        ocol (`numpy.ndarray`_, optional):
            The column in the output image for each column in `col`.
            If None, assume::

                ocol = numpy.arange(col.shape[1])

            These coordinates can fall off the output image (i.e.,
            :math:`<0` or :math:`\geq N_{\rm out,col}`), but those
            columns are removed from the output).
        max_ocol (:obj:`int`, optional):
            The last viable column *index* to include in the output
            image; ie., for an image with `ncol` columns, this should
            be `ncol-1`. If None, assume `max(ocol)`.
        extract_width (:obj:`float`, optional):
            The width of the extraction aperture to use for the image
            rectification. If None, the image recification is performed
            using :class:`Resample` along each row.
        mask_threshold (:obj:`float`, optional):
            Either due to `bpm` or the bounds of the provided `img`,
            pixels in the rectified image may not be fully covered by
            valid pixels in `img`. Pixels in the output image with
            less than this fractional coverage of an input pixel are
            flagged in the output.

    Returns:
        Two `numpy.ndarray`_ objects are returned both with shape
        `(nrow,max_ocol+1)`, the rectified image and its boolean
        bad-pixel mask.
    """
    # Check the input
    if img.ndim != 2:
        raise ValueError('Input image must be 2D.')
    if bpm is not None and bpm.shape != img.shape:
        raise ValueError('Image bad-pixel mask must match image shape.')
    _img = numpy.ma.MaskedArray(img, mask=bpm)
    nrow, ncol = _img.shape
    if col.ndim != 2:
        raise ValueError('Column mapping array must be 2D.')
    if col.shape[0] != nrow:
        raise ValueError('Number of rows in column mapping array must match image to rectify.')
    _ocol = numpy.arange(col.shape[1]) if ocol is None else numpy.atleast_1d(ocol)
    if _ocol.ndim != 1:
        raise ValueError('Output column indices must be provided as a vector.')
    if _ocol.size != col.shape[1]:
        raise ValueError('Output column indices must match columns in column mapping array.')
    _max_ocol = numpy.amax(_ocol) if max_ocol is None else max_ocol

    # Use an aperture extraction to rectify the image
    if extract_width is not None:
        # Select viable columns
        indx = (_ocol >= 0) & (_ocol <= _max_ocol)
        # Initialize the output image as all masked
        out_img = numpy.ma.masked_all((nrow,_max_ocol+1), dtype=float)
        # Perform the extraction
        out_img[:,_ocol[indx]] = moment.moment1d(_img.data, col[:,indx], extract_width,
                                                 bpm=_img.mask)[0]
        # Determine what fraction of the extraction fell off the image
        coo = col[:,indx,None] + numpy.arange(numpy.ceil(extract_width)).astype(int)[None,None,:] \
                    - extract_width/2
        in_image = (coo >= 0) | (coo < ncol)
        out_bpm = numpy.sum(in_image, axis=2)/extract_width < mask_threshold
        # Return the filled numpy.ndarray and boolean mask
        return out_img.filled(0.0)/extract_width, out_img.mask | out_bpm

    # Directly resample the image
    # 
    # `col` provides the column in the input image that should
    # resampled to a given position in the output image: the value of
    # the flux at img[:,col[:,0]] should be rectified to `outimg[:,0]`.
    # To run the resampling algorithm we need to invert this. That is,
    # instead of having a function that provides the output column as a
    # function of the input column, we want a function that provies the
    # input column as a function of the output column.

    # Instantiate the output image
    out_img = numpy.zeros((nrow,_max_ocol+1), dtype=float)
    out_bpm = numpy.zeros((nrow,_max_ocol+1), dtype=bool)
    icol = numpy.arange(ncol)
    for i in range(nrow):
        # Get the coordinate vector of the output column
        _icol = interpolate.interp1d(col[i,:], _ocol, copy=False, bounds_error=False,
                                     fill_value='extrapolate', assume_sorted=True)(icol)
        # Resample it
        r = Resample(_img[i,:], x=_icol, newRange=[0,_max_ocol], newpix=_max_ocol+1, newLog=False,
                     conserve=True)
        # Save the resampled data
        out_img[i,:] = r.outy
        # Flag pixels
        out_bpm[i,:] = r.outf < mask_threshold
    return out_img, out_bpm


