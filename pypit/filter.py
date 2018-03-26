# -*- coding: utf-8 -*-
"""
A set of functions used to filter arrays.

<<<<<<< HEAD
Originally pulled from SDSS-IV/MaNGA Data Analysis Pipeline, licensed
under BSD 3-clause license.
=======
*License*:
    Copyright (c) 2017, SDSS-IV/MaNGA Pipeline Group
        Licensed under BSD 3-clause license - see LICENSE.rst

*Source location*:
    $MANGADAP_DIR/python/mangadap/util/filter.py

*Imports and python version compliance*:
    ::

        from __future__ import division
        from __future__ import print_function
        from __future__ import absolute_import
        from __future__ import unicode_literals

        import sys
        if sys.version > '3':
            long = int

        import numpy

*Revision history*:
    | **26 Jan 2017**: Original implementation by K. Westfall (KBW)
    | **06 Apr 2017**: Interpolate sigma vectors as well as the smoothed
        vectors in :class:`BoxcarFilter`.
    | **21 Dec 2017**: Add weighting functionality to
        :class:`BoxcarFilter`.
>>>>>>> master
"""

from __future__ import division
from __future__ import print_function
from __future__ import absolute_import
from __future__ import unicode_literals

import sys
if sys.version > '3':
    long = int

import warnings
import numpy
from scipy import interpolate, sparse

<<<<<<< HEAD
class BoxcarFilter(object):
    """Boxcar filter a 2D image along columns."""
=======
# Debug
from matplotlib import pyplot
from matplotlib.ticker import NullFormatter

def high_pass_filter(flux, dw=0, k=None, Dw=None):
    """
    Pulled from FIREFLY's hpf() function and edited on 17 Nov 2016.

    Apply a high-pass filter to the input vector.

    """

    n = flux.size
    _dw = int(dw)

    if k is None and Dw is None:
        Dw = 100
        _k = n//Dw
    elif k is not None and Dw is None:
        _k = int(k)
        Dw = n//_k
    elif k is None and Dw is not None:
        _k = n//Dw
    else:
        raise ValueError('Cannot provide both k and Dw.')

    # PRINT A REPORT

    # Rita's typical inputs for SDSS:
    # w = 10 # 10
    # windowsize = 20 # 20

    # My MaNGA inputs:
    # if w == 0 and windowsize == 0:
    #     print "HPF parameters not specified: choosing default window (stepsize=40)"
    #     w = 40
    #     windowsize = 0

    h           = numpy.fft.fft(flux)
    h_filtered  = numpy.zeros(n, dtype=complex)
    window      = numpy.zeros(n)
    unwindow    = numpy.zeros(n)
    window[0]   = 1                             # keep F0 for normalisation
    unwindow[0] = 1

    for i in range(_dw):
        window[_k+i] = (i+1.0)/_dw
        window[n-1-(_k+_dw-i)] = (_dw-i)/_dw
    window[_k+_dw:n-(_k+_dw)] = 1

    unwindow        = 1 - window
    unwindow[0]     = 1
    
    h_filtered      = h * window
    un_h_filtered   = h * unwindow

    res     = numpy.real(numpy.fft.ifft(h_filtered))
    unres   = numpy.real(numpy.fft.ifft(un_h_filtered)) 
    res_out = (1.0+(res-numpy.median(res))/unres) * numpy.median(res) 

    return res_out, window, res, unres

#def off_diagonal_identity(size, win):
#    r"""
#    Construct a matrix with ones within a window along the diagonal.
#
#    Args:
#        size (int) : Size for the square matrix; i.e., :math:`N` for the
#            :math:`N\timesN` matrix.
#        win (int): Number of ones in each row along the diagonal.
#
#    Raises:
#        ValueError: Raised if the window is larger than 2*size-1.
#
#    """
#    if win > 2*size-1:
#        raise ValueError('Window too large for matrix size.')
#    if win == 2*size-1:
#        return numpy.ones((size,size), dtype=int)
#    x = numpy.zeros((size,size), dtype=int)#numpy.identity(size).astype(int)
#    for i in range(1,(win+1)//2):
#        x[:-i,i:] = x[:-i,i:] + numpy.identity(size-i).astype(int)
#    x += x.T
#    if win % 2 != 1:
#        x[:-i-1,i+1:] = x[:-i-1,i+1:] + numpy.identity(size-i-1).astype(int)
#    return x + numpy.identity(size).astype(int)

def off_diagonal_identity(size, win, return_sparse=False):
    r"""
    Construct a matrix with ones within a window along the diagonal.

    Args:
        size (int) : Size for the square matrix; i.e., :math:`N` for the
            :math:`N\timesN` matrix.
        win (int): Number of ones in each row along the diagonal.

    Raises:
        ValueError: Raised if the window is larger than 2*size-1.

    """
    if win > 2*size-1:
        raise ValueError('Window too large for matrix size.')
    if win == 2*size-1:
        return numpy.ones((size,size), dtype=int)

    # Indices of diagonal
    ii = numpy.arange(size).astype(int)

    # Build the upper triangle
    i = numpy.empty(0, dtype=int)
    j = numpy.empty(0, dtype=int)
    for k in range(1,(win+1)//2):
        i = numpy.append(i, ii[:size-k])
        j = numpy.append(j, ii[k:size])

    # Copy to the lower triangle
    _i = numpy.append(i,j)
    j = numpy.append(j,i)

    # Add the diagonal
    i = numpy.append(_i, ii)
    j = numpy.append(j, ii)

    # Accommodate an even window
    if win % 2 != 1:
        i = numpy.append(i, ii[:size-k-1])
        j = numpy.append(j, ii[k+1:size])

    # Construct and return the array
    if return_sparse:
        return sparse.coo_matrix((numpy.ones(len(i), dtype=int),(i,j)), shape=(size,size)).tocsr()

    a = numpy.zeros((size,size), dtype=int)
    a[i,j] = 1
    return a


def build_smoothing_mask(x, pix_buffer, default=None, mask_x=None):
    r"""
    Construct a mask for the provided vector that masks the first and
    last pix_buffer pixels in the coordinate vector.
    
    The mask is instantiated as fully unmasked, unless a default mask is
    provided.
    
    To mask specified ranges in x, provide mask_x with a shape
    :math:`N_{\rm mask}\times 2` where each mask is defined by the
    starting and ending value of x to exclude.
    """
    # Smooth the ratio of the data to the model
    if len(x.shape) != 1:
        raise ValueError('Input must be a vector.')
    npix = x.size
    mask = numpy.zeros(npix, dtype=bool) if default is None else default.copy()
    if mask.size != npix:
        raise ValueError('Provided default has an incorrect length.')
    mask[:pix_buffer] = True
    mask[-pix_buffer:] = True
    if mask_x is None:
        return mask

    for m in mask_x:
        mask |= numpy.logical_and(x > m[0], x < m[1])
    return mask


class BoxcarFilter():
>>>>>>> master
    def __init__(self, boxcar, lo=None, hi=None, niter=None, y=None, wgt=None, mask=None,
                 local_sigma=None):

        if boxcar <= 1:
            raise ValueError('Boxcar must be greater than 1!')

        self.y = None
        self.wgt = None
        self.input_mask = None
        self.output_mask = None
        self.boxcar = boxcar
        self.lo_rej = lo
        self.hi_rej = hi
        self.nrej = -1 if niter is None else niter
        self.nvec = None
        self.npix = None
        self.local_sigma = False if local_sigma is None else local_sigma

        self.smoothed_wgt = None
        self.smoothed_n = None
        self.smoothed_y = None
        self.sigma_y = None
        self.rdx_index = None

        if y is not None:
            self.smooth(y, wgt=wgt, mask=mask, boxcar=boxcar, lo=lo, hi=hi, niter=niter,
                        local_sigma=local_sigma)

<<<<<<< HEAD
=======

>>>>>>> master
    def _assign_par(self, boxcar, lo, hi, niter, local_sigma):
        if boxcar is not None:
            self.boxcar = boxcar
        if lo is not None:
            self.lo_rej = lo
        if hi is not None:
            self.hi_rej = hi
        if niter is not None:
            self.nrej = niter
        if local_sigma is not None:
            self.local_sigma = local_sigma

<<<<<<< HEAD
=======

>>>>>>> master
    def _check_par(self):
        if self.boxcar <= 1:
            raise ValueError('Boxcar must be greater than 1!')

        if self.nrej is not None and self.nrej > 0 and self.lo_rej is None and self.hi_rej is None:
            raise ValueError('Must provide lo or hi if niter provided')

<<<<<<< HEAD
=======

>>>>>>> master
    def _reduce_indices(self):
        indx = numpy.array([numpy.arange(self.npix)-self.boxcar//2,
                            numpy.arange(self.npix)+self.boxcar//2+1])
        if self.boxcar % 2 != 1:
            indx[:self.npix] += 1
        return indx.T.ravel().clip(0,self.npix-1)

<<<<<<< HEAD
    def _apply(self):
        self.smoothed_n = numpy.add.reduceat(numpy.invert(self.output_mask), self.rdx_index,
                                             axis=1, dtype=int)[:,::2]
=======

    def _apply(self):
        self.smoothed_n = numpy.add.reduceat(numpy.invert(self.output_mask), self.rdx_index, axis=1,
                                             dtype=int)[:,::2]
>>>>>>> master
        self.smoothed_wgt = numpy.add.reduceat(self.input_wgt, self.rdx_index, axis=1,
                                             dtype=float)[:,::2]
        self.smoothed_y = numpy.add.reduceat(self.input_wgt*self.y.filled(0.0), self.rdx_index,
                                             axis=1, dtype=float)[:,::2]
        smoothed_y2 = numpy.add.reduceat(self.input_wgt*numpy.square(self.y.filled(0.0)),
                                         self.rdx_index, axis=1, dtype=float)[:,::2]

        self.smoothed_y = numpy.ma.divide(self.smoothed_y, self.smoothed_wgt)
        self.smoothed_y[self.output_mask | (self.smoothed_n == 0)] = numpy.ma.masked

        self.sigma_y = numpy.ma.sqrt((numpy.ma.divide(smoothed_y2, self.smoothed_wgt)
                                            - numpy.square(self.smoothed_y))
                                        * numpy.ma.divide(self.smoothed_n, self.smoothed_n-1)) \
                        if self.local_sigma else \
                            numpy.ma.MaskedArray(numpy.array([numpy.std(self.y - self.smoothed_y,
                                                                        axis=1)]*self.npix).T)

        self.output_mask = numpy.ma.getmaskarray(self.smoothed_y) \
                                | numpy.ma.getmaskarray(self.sigma_y) 
        
        self.smoothed_y[self.output_mask | (self.smoothed_n == 0)] = numpy.ma.masked
        self.sigma_y[self.output_mask | (self.smoothed_n == 0)] = numpy.ma.masked

<<<<<<< HEAD
=======
        return

        w,h = pyplot.figaspect(1)
        fig = pyplot.figure(figsize=(1.5*w,1.5*h))

        minf = numpy.amin( numpy.ma.log10(numpy.append(self.y.compressed(), self.smoothed_y.compressed())) )
        maxf = numpy.amax( numpy.ma.log10(numpy.append(self.y.compressed(), self.smoothed_y.compressed())) )

        mins = numpy.amin(numpy.ma.log10(self.sigma_y).compressed())
        maxs = numpy.amax(numpy.ma.log10(self.sigma_y).compressed())

        minr = numpy.amin(numpy.ma.log10(numpy.ma.divide(self.y, self.smoothed_y)).compressed())
        maxr = numpy.amax(numpy.ma.log10(numpy.ma.divide(self.y, self.smoothed_y)).compressed())

        ax = fig.add_axes([0.05, 0.50, 0.45, 0.45])
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(self.y), origin='lower', interpolation='nearest', vmin=minf,
                  vmax=maxf, aspect='auto')
        ax = fig.add_axes([0.50, 0.50, 0.45, 0.45])
        ax.xaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(self.smoothed_y), origin='lower', interpolation='nearest',
                  vmin=minf, vmax=maxf, aspect='auto')
        ax = fig.add_axes([0.05, 0.05, 0.45, 0.45])
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(numpy.ma.divide(self.y,self.smoothed_y)), origin='lower',
                  interpolation='nearest', aspect='auto', vmin=minr, vmax=maxr)
        ax = fig.add_axes([0.50, 0.05, 0.45, 0.45])
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.imshow(numpy.ma.log10(self.sigma_y), origin='lower', interpolation='nearest',
                  aspect='auto', vmin=mins, vmax=maxs)

        pyplot.show()
        exit()

    
>>>>>>> master
    def _interpolate(self):
        """
        Interpolate the smoothed image across the masked regions.
        Leading and trailing masked regions of each row are set to the
        first and last unmasked value, respectively.

        interpolate both the smoothed data and the sigma
<<<<<<< HEAD
=======

>>>>>>> master
        """
        # Boolean arrays with bad pixels (should be the same for both
        # smoothed_y and sigma_y)
        badpix = numpy.ma.getmaskarray(self.smoothed_y)

        # Deal with any vectors that are fully masked
        goodvec = numpy.sum(numpy.invert(badpix), axis=1) > 0
        ngood = numpy.sum(goodvec)
        veci = numpy.arange(self.nvec)[goodvec]

        # Find the first and last unflagged pixels in the good vectors
        pixcoo = numpy.ma.MaskedArray(numpy.array([numpy.arange(self.npix)]*ngood),
                                      mask=badpix[goodvec,:]).astype(int)
        mini = numpy.ma.amin(pixcoo, axis=1)
        maxi = numpy.ma.amax(pixcoo, axis=1)

        # Copy those to the first and last pixels for each row
        self.smoothed_y[veci,0] = numpy.array([self.smoothed_y[i,m] for i,m in zip(veci, mini)])
        self.smoothed_y[veci,-1] = numpy.array([self.smoothed_y[i,m] for i,m in zip(veci, maxi)])

        self.sigma_y[veci,0] = numpy.array([self.sigma_y[i,m] for i,m in zip(veci, mini)])
        self.sigma_y[veci,-1] = numpy.array([self.sigma_y[i,m] for i,m in zip(veci, maxi)])

        # Detach badpix from self.smoothed_y
        badpix = numpy.ma.getmaskarray(self.smoothed_y).copy()

        # Interpolate the smoothed array
        x = numpy.ma.MaskedArray(numpy.arange(self.nvec*self.npix), mask=badpix)
        interpolator = interpolate.interp1d(x.compressed(), self.smoothed_y.compressed(),
                                            assume_sorted=True, fill_value='extrapolate')
        self.smoothed_y.ravel()[badpix.ravel()] = interpolator(x.data[badpix.ravel()])
        self.smoothed_y[numpy.invert(goodvec),:] = 0.0

        # Interpolate the sigma array
        interpolator = interpolate.interp1d(x.compressed(), self.sigma_y.compressed(),
                                            assume_sorted=True, fill_value='extrapolate')
        self.sigma_y.ravel()[badpix.ravel()] = interpolator(x.data[badpix.ravel()])
        self.sigma_y[numpy.invert(goodvec),:] = 0.0

<<<<<<< HEAD
=======

>>>>>>> master
    def smooth(self, y, wgt=None, mask=None, boxcar=None, lo=None, hi=None, niter=None,
               local_sigma=None):
        """
        Smooth a vector or array of vectors by a boxcar with specified
        rejection.
        """
        # Assign and check the binning parameters
        self._assign_par(boxcar, lo, hi, niter, local_sigma)
        self._check_par()

        # Check the input vector/array
        if len(y.shape) > 2:
            raise ValueError('Can only deal with vectors or matrices.')

        # Set the weighting
        self.input_wgt = numpy.ones(y.shape, dtype=bool) if wgt is None else wgt
        if self.input_wgt.shape != y.shape:
            raise ValueError('Weights shape does not match data.')
        if isinstance(y, numpy.ma.MaskedArray):
            self.input_wgt[numpy.ma.getmaskarray(y)] = 0.0

        # Set the mask
        self.input_mask = numpy.zeros(y.shape, dtype=bool) if mask is None else mask
        if self.input_mask.shape != y.shape:
            raise ValueError('Mask shape does not match data.')
        if isinstance(y, numpy.ma.MaskedArray):
            self.input_mask |= numpy.ma.getmaskarray(y)

        # Save the input
        provided_vector = len(y.shape) == 1
        self.y = numpy.ma.atleast_2d(numpy.ma.MaskedArray(y.copy(), mask=self.input_mask))
        self.input_wgt = numpy.atleast_2d(self.input_wgt)
        self.input_mask = numpy.atleast_2d(self.input_mask)
        self.output_mask = self.input_mask.copy()
        self.nvec, self.npix = self.y.shape

        # Create the list of reduceat indices
        self.rdx_index = self._reduce_indices()

        # Get the first iteration of the smoothed vectors
        self._apply()

        # No rejection iterations requested so return the result
        if self.lo_rej is None and self.hi_rej is None:
            self._interpolate()
            return self.smoothed_y[0,:] if provided_vector else self.smoothed_y

        # Iteratively reject outliers
        nmasked = numpy.sum(numpy.ma.getmaskarray(self.y))
        i=0
        while i > -1:
            nbad = numpy.sum(self.output_mask)
            if self.lo_rej is not None:
                self.output_mask = numpy.logical_or(self.output_mask, self.y.data - self.smoothed_y
                                                        < (-self.lo_rej*self.sigma_y))
            if self.hi_rej is not None:
                self.output_mask = numpy.logical_or(self.output_mask, self.y.data - self.smoothed_y
                                                        > (self.hi_rej*self.sigma_y))
            self.y[self.output_mask] = numpy.ma.masked
            self._apply()

            nmasked = numpy.sum(self.output_mask) - nbad
<<<<<<< HEAD
=======
#            print('Niter: {0}; Nmasked: {1}'.format(i+1, nmasked))
>>>>>>> master

            i += 1
            if i == self.nrej or nmasked == 0:
                break

        self._interpolate()
        return self.smoothed_y[0,:] if provided_vector else self.smoothed_y

<<<<<<< HEAD
=======

>>>>>>> master
    @property
    def mask(self):
        return self.output_mask.copy()

<<<<<<< HEAD
=======


def smooth_masked_vector(x, nsmooth):
    """
    Smooth a masked vector by a box of size nsmooth.
    """
    n = off_diagonal_identity(x.size, nsmooth)*numpy.invert(numpy.ma.getmaskarray(x))[None,:]
    nn = numpy.sum(n,axis=1)
    return numpy.ma.MaskedArray(numpy.dot(n,x.data), mask=numpy.invert(nn>0) | x.mask)/nn


def interpolate_masked_vector(y, quiet=True, extrap_with_median=False):
    """
    Interpolate over the masked pixels in an input vector using linear
    interpolation.
    """
    x = numpy.arange(y.size)
    indx = numpy.ma.getmaskarray(y)
    if numpy.sum(numpy.invert(indx)) < 2:
        if not quiet:
            warnings.warn('Input vector has fewer than 2 unmasked values!  Returning zero vector.')
        return numpy.zeros(y.size, dtype=y.dtype.name)
    interpolator = interpolate.interp1d(x[numpy.invert(indx)], y[numpy.invert(indx)],
                                        fill_value='extrapolate')
    _y = y.data.copy()
    _y[indx] = interpolator(x[indx])

    if extrap_with_median:
        med = numpy.median(interpolator.y)
        m_indx = (x[indx] < interpolator.x[0]) | (x[indx] > interpolator.x[-1])
        _y[indx][m_indx] = med

    return _y


def boxcar_smooth_vector(x, boxcar, mask=None, lo=None, hi=None, niter=None, return_mask=False):
    """
    Boxcar smooth an input vector.
        - Ignores masked pixels.
        - Allows for iterative positive and negative outliers.
        - Can return the mask separately from the smoothed vector.
    """
    if niter is not None and lo is None and hi is None:
        raise ValueError('Must provide lo or hi if niter provided')

    _mask = numpy.zeros(x.size, dtype=bool) if mask is None else mask
    _x = x if isinstance(x, numpy.ma.MaskedArray) else numpy.ma.MaskedArray(x)
    _x[_mask] = numpy.ma.masked

    sx = interpolate_masked_vector( smooth_masked_vector(_x, boxcar) )

#    pyplot.step(numpy.arange(_x.size), _x.data, where='mid', color='0.3', lw=0.5)
#    pyplot.step(numpy.arange(_x.size), _x, where='mid', color='k', lw=0.5)
#    pyplot.plot(numpy.arange(_x.size), sx, color='r', lw=1)
#    pyplot.show()

    if numpy.all([ x is None for x in [ lo, hi, niter ]]):
        if return_mask:
            return sx, numpy.ma.getmaskarray(_x).copy()
        return sx
    
    _niter = -1 if niter is None else niter
    nrej = numpy.sum(numpy.ma.getmaskarray(_x))
    i=0
    while i > -1:
        sig = numpy.std((_x - sx).compressed())
        mask = numpy.ma.getmaskarray(_x).copy()
        if lo is not None:
            mask = numpy.logical_or(mask, _x.data - sx < -lo*sig)
        if hi is not None:
            mask = numpy.logical_or(mask, _x.data - sx > hi*sig)
        nrej = numpy.sum(mask) - numpy.sum(_x.mask)
#        print('Niter: {0}; Nrej: {1}'.format(i+1, nrej))
        _x[mask] = numpy.ma.masked
        sx = interpolate_masked_vector( smooth_masked_vector(_x, boxcar) )

#        pyplot.step(numpy.arange(_x.size), _x.data, where='mid', color='cyan', lw=0.5)
#        pyplot.step(numpy.arange(_x.size), _x, where='mid', color='k', lw=0.5)
#        pyplot.plot(numpy.arange(_x.size), sx, color='r', lw=1)
#        pyplot.show()

        i += 1
        if i == _niter or nrej == 0:
            break

    pyplot.step(numpy.arange(_x.size), _x.data, where='mid', color='cyan', lw=0.5)
    pyplot.step(numpy.arange(_x.size), _x, where='mid', color='k', lw=0.5)
    pyplot.plot(numpy.arange(_x.size), sx, color='r', lw=1)
    pyplot.show()

    if return_mask:
        return sx, numpy.ma.getmaskarray(_x).copy()
    return sx


>>>>>>> master
