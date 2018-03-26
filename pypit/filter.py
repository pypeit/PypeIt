# -*- coding: utf-8 -*-
"""
A set of functions used to filter arrays.

Originally pulled from SDSS-IV/MaNGA Data Analysis Pipeline, licensed
under BSD 3-clause license.
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

class BoxcarFilter(object):
    """Boxcar filter a 2D image along columns."""
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

    def _check_par(self):
        if self.boxcar <= 1:
            raise ValueError('Boxcar must be greater than 1!')

        if self.nrej is not None and self.nrej > 0 and self.lo_rej is None and self.hi_rej is None:
            raise ValueError('Must provide lo or hi if niter provided')

    def _reduce_indices(self):
        indx = numpy.array([numpy.arange(self.npix)-self.boxcar//2,
                            numpy.arange(self.npix)+self.boxcar//2+1])
        if self.boxcar % 2 != 1:
            indx[:self.npix] += 1
        return indx.T.ravel().clip(0,self.npix-1)

    def _apply(self):
        self.smoothed_n = numpy.add.reduceat(numpy.invert(self.output_mask), self.rdx_index,
                                             axis=1, dtype=int)[:,::2]
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

    def _interpolate(self):
        """
        Interpolate the smoothed image across the masked regions.
        Leading and trailing masked regions of each row are set to the
        first and last unmasked value, respectively.

        interpolate both the smoothed data and the sigma
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

            i += 1
            if i == self.nrej or nmasked == 0:
                break

        self._interpolate()
        return self.smoothed_y[0,:] if provided_vector else self.smoothed_y

    @property
    def mask(self):
        return self.output_mask.copy()

