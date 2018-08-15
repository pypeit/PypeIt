# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
""" Methods taken from pydl
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from warnings import warn

from pypeit import msgs
from pypeit import debugger


# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
"""This module corresponds to the image directory in idlutils.
"""


def djs_maskinterp1(yval, mask, xval=None, const=False):
    """Interpolate over a masked, 1-d array.

    Parameters
    ----------
    yval : :class:`numpy.ndarray`
        The input values.
    mask : :class:`numpy.ndarray`
        The mask.
    xval : :class:`numpy.ndarray`, optional
        If set, use these x values, otherwise use an array.
    const : :class:`bool`, optional
        If set to ``True``, bad values around the edges of the array will be
        set to a constant value.  Because of the default behavior of
        :func:`numpy.interp`, this value actually makes no difference in
        the output.

    Returns
    -------
    :class:`numpy.ndarray`
        The `yval` array with masked values replaced by interpolated values.
    """
    import numpy as np
    good = mask == 0
    if good.all():
        return yval
    ngood = good.sum()
    igood = good.nonzero()[0]
    if ngood == 0:
        return yval
    if ngood == 1:
        return np.zeros(yval.shape, dtype=yval.dtype) + yval[igood[0]]
    ynew = yval.astype('d')
    ny = yval.size
    ibad = (mask != 0).nonzero()[0]
    if xval is None:
        ynew[ibad] = np.interp(ibad, igood, ynew[igood])
        if const:
            if igood[0] != 0:
                ynew[0:igood[0]] = ynew[igood[0]]
            if igood[ngood-1] != ny-1:
                ynew[igood[ngood-1]+1:ny] = ynew[igood[ngood-1]]
    else:
        ii = xval.argsort()
        ibad = (mask[ii] != 0).nonzero()[0]
        igood = (mask[ii] == 0).nonzero()[0]
        ynew[ii[ibad]] = np.interp(xval[ii[ibad]], xval[ii[igood]],
                                   ynew[ii[igood]])
        if const:
            if igood[0] != 0:
                ynew[ii[0:igood[0]]] = ynew[ii[igood[0]]]
            if igood[ngood-1] != ny-1:
                ynew[ii[igood[ngood-1]+1:ny]] = ynew[ii[igood[ngood-1]]]
    return ynew


def djs_maskinterp(yval, mask, xval=None, axis=None, const=False):
    """Interpolate over masked pixels in a vector, image or 3-D array.

    Parameters
    ----------
    yval : :class:`numpy.ndarray`
        The input values
    mask : :class:`numpy.ndarray`
        The mask
    xval : :class:`numpy.ndarray`, optional
        If set, use these x values, otherwise use an array
    axis : :class:`int`, optional
        Must be set if yval has more than one dimension. If set to zero,
        interpolate along the first axis of the array, if set to one,
        interpolate along the second axis of the array, and so on.
    const : :class:`bool`, optional
        This value is passed to a helper function, djs_maskinterp1.

    Returns
    -------
    :class:`numpy.ndarray`
        The interpolated array.
    """
    import numpy as np
    if mask.shape != yval.shape:
        raise ValueError('mask must have the same shape as yval.')
    if xval is not None:
        if xval.shape != yval.shape:
            raise ValueError('xval must have the same shape as yval.')
    ndim = yval.ndim
    if ndim == 1:
        ynew = djs_maskinterp1(yval, mask, xval=xval, const=const)
    else:
        if axis is None:
            raise ValueError('Must set axis if yval has more than one dimension.')
        if axis < 0 or axis > ndim-1 or axis - int(axis) != 0:
            raise ValueError('Invalid axis value.')
        ynew = np.zeros(yval.shape, dtype=yval.dtype)
        if ndim == 2:
            if xval is None:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        ynew[i, :] = djs_maskinterp1(yval[i, :], mask[i, :],
                                                     const=const)
                else:
                    for i in range(yval.shape[1]):
                        ynew[:, i] = djs_maskinterp1(yval[:, i], mask[:, i],
                                                     const=const)
            else:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        ynew[i, :] = djs_maskinterp1(yval[i, :], mask[i, :],
                                                     xval=xval[i, :],
                                                     const=const)
                else:
                    for i in range(yval.shape[1]):
                        ynew[:, i] = djs_maskinterp1(yval[:, i], mask[:, i],
                                                     xval=xval[:, i],
                                                     const=const)
        elif ndim == 3:
            if xval is None:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[1]):
                            ynew[i, j, :] = djs_maskinterp1(yval[i, j, :],
                                                            mask[i, j, :],
                                                            const=const)
                elif axis == 1:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[2]):
                            ynew[i, :, j] = djs_maskinterp1(yval[i, :, j],
                                                            mask[i, :, j],
                                                            const=const)
                else:
                    for i in range(yval.shape[1]):
                        for j in range(yval.shape[2]):
                            ynew[:, i, j] = djs_maskinterp1(yval[:, i, j],
                                                            mask[:, i, j],
                                                            const=const)
            else:
                if axis == 0:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[1]):
                            ynew[i, j, :] = djs_maskinterp1(yval[i, j, :],
                                                            mask[i, j, :],
                                                            xval=xval[i, j, :],
                                                            const=const)
                elif axis == 1:
                    for i in range(yval.shape[0]):
                        for j in range(yval.shape[2]):
                            ynew[i, :, j] = djs_maskinterp1(yval[i, :, j],
                                                            mask[i, :, j],
                                                            xval=xval[i, :, j],
                                                            const=const)
                else:
                    for i in range(yval.shape[1]):
                        for j in range(yval.shape[2]):
                            ynew[:, i, j] = djs_maskinterp1(yval[:, i, j],
                                                            mask[:, i, j],
                                                            xval=xval[:, i, j],
                                                            const=const)
        else:
            raise ValueError('Unsupported number of dimensions.')
    return ynew



class bspline(object):
    """Bspline class.

    Functions in the bspline library are implemented as methods on this
    class.

    Parameters
    ----------
    x : :class:`numpy.ndarray`
        The data.
    nord : :class:`int`, optional
        To be documented.
    npoly : :class:`int`, optional
        To be documented.
    bkpt : :class:`numpy.ndarray`, optional
        To be documented.
    bkspread : :class:`float`, optional
        To be documented.
    verbose : :class:`bool`, optional.
        If ``True`` print extra information.

    Attributes
    ----------
    breakpoints
        Breakpoints for bspline, spacing for these breakpoints are determinated by keywords inputs;
    nord
        Order of bspline; [default=4]
    npoly
        Polynomial order to fit over 2nd variable (when specified as x2): [default=1]
    mask
        Output mask, set =1 for good points, =0 for bad points;
    coeff
        Output coefficient of the bspline;
    icoeff
        Cholesky band matrix used to solve for the bspline coefficients;
    xmin
        Normalization minimum for x2; [default max(xdata)]
    xmax
        Normalization maximum for x2; [default min(xdata)]
    funcname
        Function for the second variable; [default 'legendre']
    from_dict
        If not None, create a bspline from a dictionary created by to_dict(). [default 'None']
        It is possible to instantiate a bspline from a dict without the x data:
        new_bspline = bspline(None, from_dict=dictionary)
    """

    def __init__(self, x, nord=4, npoly=1, bkpt=None, fullbkpt = None, bkspread=1.0,
                 verbose=False, from_dict=None, **kwargs):
        """Init creates an object whose attributes are similar to the
        structure returned by the create_bspline function.
        """
        #ToDO Consider refactoring the argument list so that there are no kwargs

        if from_dict is not None:
            self.nord=from_dict['nord'],
            self.npoly=from_dict['npoly'],
            self.breakpoints=np.array(from_dict['breakpoints']),
            self.mask=np.array(from_dict['mask']),
            self.coeff=np.array(from_dict['coeff']),
            self.icoeff=np.array(from_dict['icoeff']),
            self.xmin=from_dict['xmin'],
            self.xmax=from_dict['xmax'],
            self.funcname=from_dict['funcname']
        else:
            #
            # Set the breakpoints.
            #
            if fullbkpt is None:
                if bkpt is None:
                    startx = x.min()
                    rangex = x.max() - startx
                    if 'placed' in kwargs:
                        w = ((kwargs['placed'] >= startx) &
                             (kwargs['placed'] <= startx+rangex))
                        if w.sum() < 2:
                            bkpt = np.arange(2, dtype='f') * rangex + startx
                        else:
                            bkpt = kwargs['placed'][w]
                    elif 'bkspace' in kwargs:
                        nbkpts = int(rangex/kwargs['bkspace']) + 1
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt = np.arange(nbkpts, dtype='f')*tempbkspace + startx
                    elif 'nbkpts' in kwargs:
                        nbkpts = kwargs['nbkpts']
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt = np.arange(nbkpts, dtype='f') * tempbkspace + startx
                    elif 'everyn' in kwargs:
                        nx = x.size
                        nbkpts = max(nx/kwargs['everyn'], 1)
                        if nbkpts == 1:
                            xspot = [0]
                        else:
                            xspot = (nx/nbkpts)*np.arange(nbkpts)
                            # JFH This was a bug. Made fixes
                            #xspot = int(nx/(nbkpts-1)) * np.arange(nbkpts, dtype='i4')
                        #bkpt = x[xspot].astype('f')
                        bkpt = np.interp(xspot,np.arange(nx),x)
                    else:
                        raise ValueError('No information for bkpts.')
                imin = bkpt.argmin()
                imax = bkpt.argmax()
                if x.min() < bkpt[imin]:
                    if verbose:
                        print('Lowest breakpoint does not cover lowest x value: changing.')
                    bkpt[imin] = x.min()
                if x.max() > bkpt[imax]:
                    if verbose:
                        print('Highest breakpoint does not cover highest x value: changing.')
                    bkpt[imax] = x.max()
                nshortbkpt = bkpt.size
                fullbkpt = bkpt.copy()
                if nshortbkpt == 1:
                    bkspace = np.float32(bkspread)
                else:
                    bkspace = (bkpt[1] - bkpt[0]) * np.float32(bkspread)
                for i in np.arange(1, nord, dtype=np.float32):
                    fullbkpt = np.insert(fullbkpt, 0, bkpt[0]-bkspace*i)
                    fullbkpt = np.insert(fullbkpt, fullbkpt.shape[0],
                                         bkpt[nshortbkpt-1] + bkspace*i)

            #
            # Set the attributes
            #
            nc = fullbkpt.size - nord
            self.breakpoints = fullbkpt
            self.nord = nord
            self.npoly = npoly
            self.mask = np.ones((fullbkpt.size,), dtype='bool')
            if npoly > 1:
                self.coeff = np.zeros((npoly, nc), dtype='d')
                self.icoeff = np.zeros((npoly, nc), dtype='d')
            else:
                self.coeff = np.zeros((nc,), dtype='d')
                self.icoeff = np.zeros((nc,), dtype='d')
            self.xmin = 0.0
            self.xmax = 1.0
            if 'funcname' in kwargs:
                self.funcname = kwargs['funcname']
            else:
                self.funcname = 'legendre'

    def to_dict(self):
        """Write bspline parameters to a dict.

        Parameters
        ----------

        Returns
        -------
            A dict containing the relevant bspline paramenters.

        Notes
        -----
        The dictionary is JSON compatible.
        """

        # needs to move np.arrays to lists for JSON files
        return (dict(breakpoints=self.breakpoints.tolist(),
                     nord=self.nord,
                     npoly=self.npoly,
                     mask=self.mask.tolist(),
                     coeff=self.coeff.tolist(),
                     icoeff=self.icoeff.tolist(),
                     xmin=self.xmin,
                     xmax=self.xmax,
                     funcname=self.funcname))

    def fit(self, xdata, ydata, invvar, x2=None):
        """Calculate a B-spline in the least-squares sense.

        Fit is based on two variables: x which is sorted and spans a large range
        where bkpts are required y which can be described with a low order
        polynomial.

        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable.
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of `ydata`.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing an integer error code, and the evaluation of the
            b-spline at the input values.  An error code of -2 is a failure,
            -1 indicates dropped breakpoints, 0 is success, and positive
            integers indicate ill-conditioned breakpoints.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype='f')
            return (-2, yfit)
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        a1, lower, upper = self.action(xdata, x2=x2)
        foo = np.tile(invvar, bw).reshape(bw, invvar.size).transpose()
        a2 = a1 * foo
        alpha = np.zeros((bw, nfull+bw), dtype='d')
        beta = np.zeros((nfull+bw,), dtype='d')
        bi = np.arange(bw, dtype='i4')
        bo = np.arange(bw, dtype='i4')
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw-k, dtype='i4')+(bw+1)*k)
            bo = np.append(bo, np.arange(bw-k, dtype='i4')+bw*k)
        for k in range(nn-self.nord+1):
            itop = k*self.npoly
            ibottom = min(itop, nfull) + bw - 1
            ict = upper[k] - lower[k] + 1
            if ict > 0:
                work = np.dot(a1[lower[k]:upper[k]+1, :].T, a2[lower[k]:upper[k]+1, :])
                wb = np.dot(ydata[lower[k]:upper[k]+1], a2[lower[k]:upper[k]+1, :])
                alpha.T.flat[bo+itop*bw] += work.flat[bi]
                beta[itop:ibottom+1] += wb
        min_influence = 1.0e-10 * invvar.sum() / nfull
        errb = cholesky_band(alpha, mininf=min_influence)  # ,verbose=True)
        if isinstance(errb[0], int) and errb[0] == -1:
            a = errb[1]
        else:
            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
            return (self.maskpoints(errb[0]), yfit)
        errs = cholesky_solve(a, beta)
        if isinstance(errs[0], int) and errs[0] == -1:
            sol = errs[1]
        else:
            #
            # It is not possible for this to get called, because cholesky_solve
            # has only one return statement, & that statement guarantees that
            # errs[0] == -1
            #
            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
            return (self.maskpoints(errs[0]), yfit)
        if self.coeff.ndim == 2:
            # JFH made major bug fix here.
            self.icoeff[:, goodbk] = np.array(a[0, 0:nfull].T.reshape(self.npoly, nn,order='F'), dtype=a.dtype)
            self.coeff[:, goodbk] = np.array(sol[0:nfull].T.reshape(self.npoly, nn, order='F'), dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)
        yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
        return (0, yfit)

    def action(self, x, x2=None):
        """Construct banded bspline matrix, with dimensions [ndata, bandwidth].

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :func:`tuple`
            A tuple containing the b-spline action matrix; the 'lower' parameter,
            a list of pixel positions, each corresponding to the first
            occurence of position greater than breakpoint indx; and 'upper',
            Same as lower, but denotes the upper pixel positions.
        """
        nx = x.size
        nbkpt = self.mask.sum()
        if nbkpt < 2*self.nord:
            return (-2, 0, 0)
        n = nbkpt - self.nord
        gb = self.breakpoints[self.mask]
        bw = self.npoly*self.nord
        lower = np.zeros((n - self.nord + 1,), dtype=int)
        upper = np.zeros((n - self.nord + 1,), dtype=int) - 1
        indx = self.intrv(x)
        bf1 = self.bsplvn(x, indx)
        action = bf1
        aa = uniq(indx, np.arange(indx.size, dtype=int))
        upper[indx[aa]-self.nord+1] = aa
        rindx = indx[::-1]
        bb = uniq(rindx, np.arange(rindx.size, dtype=int))
        lower[rindx[bb]-self.nord+1] = nx - bb - 1
        if x2 is not None:
            if x2.size != nx:
                raise ValueError('Dimensions of x and x2 do not match.')
            x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
            if self.funcname == 'poly':
                temppoly = np.ones((nx, self.npoly), dtype='f')
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'poly1':
                temppoly = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
                for i in range(1, self.npoly):
                    temppoly[:, i] = temppoly[:, i-1] * x2norm
            elif self.funcname == 'chebyshev':
                # JFH fixed bug here where temppoly needed to be transposed because of different IDL and python array conventions
                temppoly = fchebyshev(x2norm, self.npoly).T
            elif self.funcname == 'legendre':
                temppoly = flegendre(x2norm, self.npoly).T
            else:
                raise ValueError('Unknown value of funcname.')
            action = np.zeros((nx, bw), dtype='d')
            counter = -1
            for ii in range(self.nord):
                for jj in range(self.npoly):
                    counter += 1
                    action[:, counter] = bf1[:, ii]*temppoly[:, jj]
        return (action, lower, upper)

    def intrv(self, x):
        """Find the segment between breakpoints which contain each value in the array x.

        The minimum breakpoint is nbkptord -1, and the maximum
        is nbkpt - nbkptord - 1.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Data values, assumed to be monotonically increasing.

        Returns
        -------
        :class:`numpy.ndarray`
            Position of array elements with respect to breakpoints.
        """
        gb = self.breakpoints[self.mask]
        n = gb.size - self.nord
        indx = np.zeros((x.size,), dtype=int)
        ileft = self.nord - 1
        for i in range(x.size):
            while x[i] > gb[ileft+1] and ileft < n - 1:
                ileft += 1
            indx[i] = ileft
        return indx

    def bsplvn(self, x, ileft):
        """To be documented.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            To be documented.
        ileft : :class:`int`
            To be documented

        Returns
        -------
        :class:`numpy.ndarray`
            To be documented.
        """
        bkpt = self.breakpoints[self.mask]
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype)
        deltap = vnikx.copy()
        deltam = vnikx.copy()
        j = 0
        vnikx[:, 0] = 1.0
        while j < self.nord - 1:
            ipj = ileft+j+1
            deltap[:, j] = bkpt[ipj] - x
            imj = ileft-j
            deltam[:, j] = x - bkpt[imj]
            vmprev = 0.0
            for l in range(j+1):
                vm = vnikx[:, l]/(deltap[:, l] + deltam[:, j-l])
                vnikx[:, l] = vm*deltap[:, l] + vmprev
                vmprev = vm*deltam[:, j-l]
            j += 1
            vnikx[:, j] = vmprev
        return vnikx

    def value(self, x, x2=None, action=None, lower=None, upper=None):
        """Evaluate a bspline at specified values.

        Parameters
        ----------
        x : :class:`numpy.ndarray`
            Independent variable.
        x2 : :class:`numpy.ndarray`, optional
            Orthogonal dependent variable for 2d fits.
        action : :class:`numpy.ndarray`, optional
            Action matrix to use.  If not supplied it is calculated.
        lower : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.
        upper : :class:`numpy.ndarray`, optional
            If the action parameter is supplied, this parameter must also
            be supplied.

        Returns
        -------
        :func:`tuple`
            A tuple containing the results of the bspline evaluation and a
            mask indicating where the evaluation was good.
        """

        xsort = x.argsort()
        xwork = x[xsort]
        if x2 is not None:
            x2work = x2[xsort]
        else:
            x2work = None
        if action is not None:
            if lower is None or upper is None:
                raise ValueError('Must specify lower and upper if action is set.')
        else:
            action, lower, upper = self.action(xwork, x2=x2work)
        yfit = np.zeros(x.shape, dtype=x.dtype)
        bw = self.npoly * self.nord
        spot = np.arange(bw, dtype='i4')
        goodbk = self.mask.nonzero()[0]
        coeffbk = self.mask[self.nord:].nonzero()[0]
        n = self.mask.sum() - self.nord
        if self.npoly > 1:
            goodcoeff = self.coeff[:, coeffbk]
        else:
            goodcoeff = self.coeff[coeffbk]
        # maskthis = np.zeros(xwork.shape,dtype=xwork.dtype)
        for i in range(n-self.nord+1):
            ict = upper[i] - lower[i] + 1
            if ict > 0:
                yfit[lower[i]:upper[i]+1] = np.dot(action[lower[i]:upper[i]+1, :], (goodcoeff.flatten('F'))[i*self.npoly+spot])
        yy = yfit.copy()
        yy[xsort] = yfit
        mask = np.ones(x.shape, dtype='bool')
        gb = self.breakpoints[goodbk]
        outside = ((x < gb[self.nord-1]) | (x > gb[n]))
        if outside.any():
            mask[outside] = False
        hmm = ((np.diff(goodbk) > 2).nonzero())[0]
        for jj in range(hmm.size):
            inside = ((x >= self.breakpoints[goodbk[hmm[jj]]]) &
                      (x <= self.breakpoints[goodbk[hmm[jj]+1]-1]))
            if inside.any():
                mask[inside] = False
        return (yy, mask)

    def maskpoints(self, err):
        """Perform simple logic of which breakpoints to mask.


        Parameters
        ----------
        err : :class:`numpy.ndarray` or int
            The list of indexes returned by the cholesky routines.
            This is indexed to the set of currently *good* breakpoints (i.e. self.mask=True)
            And the first nord are skipped

        Returns
        -------
        :class:`int`
            An integer indicating the results of the masking.  -1 indicates
            that the error points were successfully masked.  -2 indicates
            failure; the calculation should be aborted.

        Notes
        -----
        The mask attribute is modified, assuming it is possible to create the
        mask.
        """
        # Recast err as an array if a single value int was passed in (occasional)
        if not isinstance(err, np.ndarray):
            err = np.array([err])
        # Currently good points
        goodbkpt = np.where(self.mask)[0]
        nbkpt = len(goodbkpt)
        if nbkpt <= 2*self.nord:
            return -2
        # Find the unique ones for the polynomial
        hmm = err[uniq(err//self.npoly)]//self.npoly

        n = nbkpt - self.nord
        if np.any(hmm >= n):
            return -2
        test = np.zeros(nbkpt, dtype='bool')
        for jj in range(-int(np.ceil(self.nord/2)), int(self.nord/2.)):
            foo = np.where((hmm+jj) > 0, hmm+jj, np.zeros(hmm.shape, dtype=hmm.dtype))
            inside = np.where((foo+self.nord) < n-1, foo+self.nord, np.zeros(hmm.shape, dtype=hmm.dtype)+n-1)
            if len(inside)>0:
                test[inside] = True
        if test.any():
            reality = goodbkpt[test]
            if self.mask[reality].any():
                self.mask[reality] = False
                return -1
            else:
                return -2
        else:
            return -2

    def workit(self, xdata, ydata, invvar, action,lower,upper):
        """An internal routine for bspline_extract and bspline_radial which solve a general
        banded correlation matrix which is represented by the variable "action".  This routine
        only solves the linear system once, and stores the coefficients in sset. A non-zero return value
        signifies a failed inversion


        Parameters
        ----------
        xdata : :class:`numpy.ndarray`
            Independent variable.
        ydata : :class:`numpy.ndarray`
            Dependent variable.
        invvar : :class:`numpy.ndarray`
            Inverse variance of `ydata`.
        action : :class:`numpy.ndarray`
            Banded correlation matrix
        lower  : :class:`numpy.ndarray`
            A list of pixel positions, each corresponding to the first occurence of position greater than breakpoint indx
        upper  : :class:`numpy.ndarray`
            Same as lower, but denotes the upper pixel positions

        Returns
        -------
        :func:`tuple` (success, yfit)
            A tuple containing an boolean error code, and the evaluation of the b-spline yfit at the input values.  The error codes are as follows:

                 0 is good
                -1 is dropped breakpoints, try again
                -2 is failure, should abort

        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype='f')
            return (-2, yfit)
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        foo = np.sqrt(np.tile(invvar, bw).reshape(bw, invvar.size).transpose())
        a2 = action * foo
        #a2 = action*np.sqrt(np.outer(invvar,np.ones(bw)))

        alpha = np.zeros((bw, nfull+bw), dtype='d')
        beta = np.zeros((nfull+bw,), dtype='d')
        bi = np.arange(bw, dtype='i4')
        bo = np.arange(bw, dtype='i4')
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw-k, dtype='i4')+(bw+1)*k)
            bo = np.append(bo, np.arange(bw-k, dtype='i4')+bw*k)
        for k in range(nn-self.nord+1):
            itop = k*self.npoly
            ibottom = min(itop, nfull) + bw - 1
            try:
                ict = upper[k] - lower[k] + 1
            except:
                debugger.set_trace()
            if ict > 0:
                work = np.dot(a2[lower[k]:upper[k]+1, :].T, a2[lower[k]:upper[k]+1, :])
                wb = np.dot(ydata[lower[k]:upper[k]+1]*np.sqrt(invvar[lower[k]:upper[k]+1]), a2[lower[k]:upper[k]+1, :])
                alpha.T.flat[bo+itop*bw] += work.flat[bi]
                beta[itop:ibottom+1] += wb
        min_influence = 1.0e-10 * invvar.sum() / nfull
        # Right now we are not returning the covariance, although it may arise that we should
        covariance = alpha
        errb = cholesky_band(alpha, mininf=min_influence)  # ,verbose=True)
        if isinstance(errb[0], int) and errb[0] == -1: # successful cholseky_band returns -1
            a = errb[1]
        else:
            yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)
            return (self.maskpoints(errb[0]), yfit)
        errs = cholesky_solve(a, beta)
        if isinstance(errs[0], int) and errs[0] == -1:
            sol = errs[1]
        else:
            #
            # It is not possible for this to get called, because cholesky_solve
            # has only one return statement, & that statement guarantees that
            # errs[0] == -1
            #
            yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)
            return (self.maskpoints(errs[0]), yfit)

        if self.coeff.ndim == 2:
            self.icoeff[:, goodbk] = np.array(a[0, 0:nfull].T.reshape(self.npoly, nn,order='F'), dtype=a.dtype)
            self.coeff[:, goodbk] = np.array(sol[0:nfull].T.reshape(self.npoly, nn, order='F'), dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)

        yfit, foo = self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)
        return (0, yfit)



def cholesky_band(l, mininf=0.0):
    """Compute Cholesky decomposition of banded matrix.

    Parameters
    ----------
    l : :class:`numpy.ndarray`
        A matrix on which to perform the Cholesky decomposition.
    mininf : :class:`float`, optional
        Entries in the `l` matrix are considered negative if they are less
        than this value (default 0.0).

    Returns
    -------
    :func:`tuple`
        If problems were detected, the first item will be the index or
        indexes where the problem was detected, and the second item will simply
        be the input matrix.  If no problems were detected, the first item
        will be -1, and the second item will be the Cholesky decomposition.
    """
    #from . import PydlutilsUserWarning
    lower = l.copy()
    bw, nn = lower.shape
    n = nn - bw
    negative = lower[0, 0:n] <= mininf
    if negative.any() or not np.all(np.isfinite(lower)):
        msgs.warn('Bad entries: ' + str(negative.nonzero()[0]))
        return (negative.nonzero()[0], l)
    kn = bw - 1
    spot = np.arange(kn, dtype='i4') + 1
    bi = np.arange(kn, dtype='i4')
    for i in range(1, kn):
        bi = np.append(bi, np.arange(kn-i, dtype='i4') + (kn+1)*i)
    for j in range(n):
        lower[0, j] = np.sqrt(lower[0, j])
        lower[spot, j] /= lower[0, j]
        x = lower[spot, j]
        if not np.all(np.isfinite(x)):
            msgs.warn('NaN found in cholesky_band.')
            return (j, l)
        hmm = np.outer(x, x)
        here = bi+(j+1)*bw
        lower.T.flat[here] -= hmm.flat[bi]
    return (-1, lower)


def cholesky_solve(a, bb):
    """Solve the equation Ax=b where A is a Cholesky-banded matrix.

    Parameters
    ----------
    a : :class:`numpy.ndarray`
        :math:`A` in :math:`A x = b`.
    bb : :class:`numpy.ndarray`
        :math:`b` in :math:`A x = b`.

    Returns
    -------
    :func:`tuple`
        A tuple containing the status and the result of the solution.  The
        status is always -1.
    """
    b = bb.copy()
    bw = a.shape[0]
    n = b.shape[0] - bw
    kn = bw - 1
    spot = np.arange(kn, dtype='i4') + 1
    for j in range(n):
        b[j] /= a[0, j]
        b[j+spot] -= b[j]*a[spot, j]
    spot = kn - np.arange(kn, dtype='i4')
    for j in range(n-1, -1, -1):
        b[j] = (b[j] - np.sum(a[spot, j] * b[j+spot]))/a[0, j]
    return (-1, b)


def iterfit(xdata, ydata, invvar=None, upper=5, lower=5, x2=None,
            maxiter=10, nord = 4, bkpt = None, fullbkpt = None, kwargs_bspline={}, kwargs_reject={}):
    """Iteratively fit a b-spline set to data, with rejection.

    Parameters
    ----------
    xdata : :class:`numpy.ndarray`
        Independent variable.
    ydata : :class:`numpy.ndarray`
        Dependent variable.
    invvar : :class:`numpy.ndarray`
        Inverse variance of `ydata`.  If not set, it will be calculated based
        on the standard deviation.
    upper : :class:`int` or :class:`float`
        Upper rejection threshold in units of sigma, defaults to 5 sigma.
    lower : :class:`int` or :class:`float`
        Lower rejection threshold in units of sigma, defaults to 5 sigma.
    x2 : :class:`numpy.ndarray`, optional
        Orthogonal dependent variable for 2d fits.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations, default 10.  Set this to
        zero to disable rejection.

    Returns
    -------
    :func:`tuple`
        A tuple containing the fitted bspline object and an output mask.
    """
    #from .math import djs_reject
    nx = xdata.size
    if ydata.size != nx:
        raise ValueError('Dimensions of xdata and ydata do not agree.')
    if invvar is not None:
        if invvar.size != nx:
            raise ValueError('Dimensions of xdata and invvar do not agree.')
    else:
        #
        # This correction to the variance makes it the same
        # as IDL's variance()
        #
        var = ydata.var()*(float(nx)/float(nx-1))
        if var == 0:
            var = 1.0
        invvar = np.ones(ydata.shape, dtype=ydata.dtype)/var
    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & (invvar > 0))[xsort]
    if 'oldset' in kwargs_bspline:
        sset = kwargs_bspline['oldset']
        sset.mask = True
        sset.coeff = 0
    else:
        if not maskwork.any():
            raise ValueError('No valid data points.')
            # return (None,None)
        # JFH comment this out for now
        #        if 'fullbkpt' in kwargs:
        #            fullbkpt = kwargs['fullbkpt']
        else:
            sset = bspline(xdata[xsort[maskwork]], nord = nord, bkpt = bkpt, fullbkpt = fullbkpt, **kwargs_bspline)
            if maskwork.sum() < sset.nord:
                print('Number of good data points fewer than nord.')
                return (sset, outmask)
            if x2 is not None:
                if 'xmin' in kwargs_bspline:
                    xmin = kwargs_bspline['xmin']
                else:
                    xmin = x2.min()
                if 'xmax' in kwargs_bspline:
                    xmax = kwargs_bspline['xmax']
                else:
                    xmax = x2.max()
                if xmin == xmax:
                    xmax = xmin + 1
                sset.xmin = xmin
                sset.xmax = xmax
                if 'funcname' in kwargs_bspline:
                    sset.funcname = kwargs_bspline['funcname']
    xwork = xdata[xsort]
    ywork = ydata[xsort]
    invwork = invvar[xsort]
    if x2 is not None:
        x2work = x2[xsort]
    else:
        x2work = None
    iiter = 0
    error = -1
    # JFH fixed major bug here. Codes were not iterating
    qdone = False
    while (error != 0 or qdone is False) and iiter <= maxiter:
        goodbk = sset.mask.nonzero()[0]
        if maskwork.sum() <= 1 or not sset.mask.any():
            sset.coeff = 0
            iiter = maxiter + 1 # End iterations
        else:
            if 'requiren' in kwargs_bspline:
                i = 0
                while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx-1:
                    i += 1
                ct = 0
                for ileft in range(sset.nord, sset.mask.sum()-sset.nord+1):
                    while (xwork[i] >= sset.breakpoints[goodbk[ileft]] and
                           xwork[i] < sset.breakpoints[goodbk[ileft+1]] and
                           i < nx-1):
                        ct += invwork[i]*maskwork[i] > 0
                        i += 1
                    if ct >= kwargs_bspline['requiren']:
                        ct = 0
                    else:
                        sset.mask[goodbk[ileft]] = False
            error, yfit = sset.fit(xwork, ywork, invwork*maskwork,
                                   x2=x2work)
        iiter += 1
        inmask = maskwork
        if error == -2:

            return (sset, outmask)
        elif error == 0:
            maskwork, qdone = djs_reject(ywork, yfit, invvar=invwork,
                                         inmask=inmask, outmask=maskwork,
                                         upper=upper, lower=lower,**kwargs_reject)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset, outmask)

def uniq(x, index=None):
    """Replicates the IDL ``UNIQ()`` function.

    Returns the *subscripts* of the unique elements of an array.  The elements
    must actually be *sorted* before being passed to this function.  This can
    be done by sorting `x` explicitly or by passing the array subscripts that
    sort `x` as a second parameter.

    Parameters
    ----------
    x : array-like
        Search this array for unique items.
    index : array-like, optional
        This array provides the array subscripts that sort `x`.

    Returns
    -------
    array-like
        The subscripts of `x` that are the unique elements of `x`.

    Notes
    -----
    Given a sorted array, and assuming that there is a set of
    adjacent identical items, ``uniq()`` will return the subscript of the
    *last* unique item.  This charming feature is retained for
    reproducibility.

    References
    ----------
    http://www.harrisgeospatial.com/docs/uniq.html

    Examples
    --------
    >>> import numpy as np
    >>> from pydl import uniq
    >>> data = np.array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
    >>> print(uniq(np.sort(data)))
    [ 3  5  7  9 10 11 12 13]
    """
    from numpy import array, roll
    if index is None:
        indicies = (x != roll(x, -1)).nonzero()[0]
        if indicies.size > 0:
            return indicies
        else:
            return array([x.size - 1, ])
    else:
        q = x[index]
        indicies = (q != roll(q, -1)).nonzero()[0]
        if indicies.size > 0:
            return index[indicies]
        else:
            return array([q.size - 1, ], dtype=index.dtype)




def flegendre(x, m):
    """Compute the first `m` Legendre polynomials.

    Parameters
    ----------
    x : array-like
        Compute the Legendre polynomials at these abscissa values.
    m : :class:`int`
        The number of Legendre polynomials to compute.  For example, if
        :math:`m = 3`, :math:`P_0 (x)`, :math:`P_1 (x)` and :math:`P_2 (x)`
        will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    import numpy as np
    from scipy.special import legendre
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Number of Legendre polynomials must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = np.polyval(legendre(k), x)
    return leg


def fchebyshev(x, m):
    """Compute the first `m` Chebyshev polynomials.

    Parameters
    ----------
    x : array-like
        Compute the Chebyshev polynomials at these abscissa values.
    m : :class:`int`
        The number of Chebyshev polynomials to compute.  For example, if
        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and
        :math:`T_2 (x)` will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    from scipy.special import chebyt
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of Chebyshev polynomial must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = np.polyval(chebyt(k), x)
    return leg

def djs_reject(data, model, outmask=None, inmask=None, sigma=None,
               invvar=None, lower=None, upper=None, maxdev=None,
               maxrej=None, groupdim=None, groupsize=None, groupbadpix=False,
               grow=0, sticky=False):
    """Routine to reject points when doing an iterative fit to data.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The data
    model : :class:`numpy.ndarray`
        The model, must have the same number of dimensions as `data`.
    outmask : :class:`numpy.ndarray`, optional
        Output mask, generated by a previous call to `djs_reject`.  If not supplied,
        this mask will be initialized to a mask that masks nothing.  Although
        this parameter is technically optional, it will almost always be set.
    inmask : :class:`numpy.ndarray`, optional
        Input mask.  Bad points are marked with a value that evaluates to ``False``.
        Must have the same number of dimensions as `data`.
    sigma : :class:`numpy.ndarray`, optional
        Standard deviation of the data, used to reject points based on the values
        of `upper` and `lower`.
    invvar : :class:`numpy.ndarray`, optional
        Inverse variance of the data, used to reject points based on the values
        of `upper` and `lower`.  If both `sigma` and `invvar` are set, `invvar`
        will be ignored.
    lower : :class:`int` or :class:`float`, optional
        If set, reject points with data < model - lower * sigma.
    upper : :class:`int` or :class:`float`, optional
        If set, reject points with data > model + upper * sigma.
    maxdev : :class:`int` or :class:`float`, optional
        If set, reject points with abs(data-model) > maxdev.  It is permitted to
        set all three of `lower`, `upper` and `maxdev`.
    maxrej : :class:`int` or :class:`numpy.ndarray`, optional
        Maximum number of points to reject in this iteration.  If `groupsize` or
        `groupdim` are set to arrays, this should be an array as well.
    groupdim
        To be documented.
    groupsize
        To be documented.
    groupbadpix : :class:`bool`, optional
        If set to ``True``, consecutive sets of bad pixels are considered groups,
        overriding the values of `groupsize`.
    grow : :class:`int`, optional
        If set to a non-zero integer, N, the N nearest neighbors of rejected
        pixels will also be rejected.
    sticky : :class:`bool`, optional
        If set to ``True``, pixels rejected in one iteration remain rejected in
        subsequent iterations, even if the model changes.

    Returns
    -------
    :func:`tuple`
        A tuple containing a mask where rejected data values are ``False`` and
        a boolean value set to ``True`` if `djs_reject` believes there is no
        further rejection to be done.

    Raises
    ------
    ValueError
        If dimensions of various inputs do not match.
    """
    #from .misc import djs_laxisnum
    #
    # Create outmask setting = 1 for good data.
    #
    if outmask is None:
        outmask = np.ones(data.shape, dtype='bool')
    else:
        if data.shape != outmask.shape:
            raise ValueError('Dimensions of data and outmask do not agree.')
    #
    # Check other inputs.
    #
    if model is None:
        if inmask is not None:
            outmask = inmask
        return (outmask, False)
    else:
        if data.shape != model.shape:
            raise ValueError('Dimensions of data and model do not agree.')
    if inmask is not None:
        if data.shape != inmask.shape:
            raise ValueError('Dimensions of data and inmask do not agree.')
    if maxrej is not None:
        if groupdim is not None:
            if len(maxrej) != len(groupdim):
                raise ValueError('maxrej and groupdim must have the same number of elements.')
        else:
            groupdim = []
        if groupsize is not None:
            if len(maxrej) != len(groupsize):
                raise ValueError('maxrej and groupsize must have the same number of elements.')
        else:
            groupsize = len(data)
    if sigma is None and invvar is None:
        if inmask is not None:
            igood = (inmask & outmask).nonzero()[0]
        else:
            igood = outmask.nonzero()[0]
        if len(igood > 1):
            sigma = np.std(data[igood] - model[igood])
        else:
            sigma = 0
    diff = data - model
    #
    # The working array is badness, which is set to zero for good points
    # (or points already rejected), and positive values for bad points.
    # The values determine just how bad a point is, either corresponding
    # to the number of sigma above or below the fit, or to the number
    # of multiples of maxdev away from the fit.
    #
    badness = np.zeros(outmask.shape, dtype=data.dtype)
    #
    # Decide how bad a point is according to lower.
    #
    if lower is not None:
        if sigma is not None:
            qbad = diff < (-lower * sigma)
            badness += ((-diff/(sigma + (sigma == 0))) > 0) * qbad
        else:
            qbad = (diff * np.sqrt(invvar)) < -lower
            badness += ((-diff * np.sqrt(invvar)) > 0) * qbad
    #
    # Decide how bad a point is according to upper.
    #
    if upper is not None:
        if sigma is not None:
            qbad = diff > (upper * sigma)
            badness += ((diff/(sigma + (sigma == 0))) > 0) * qbad
        else:
            qbad = (diff * np.sqrt(invvar)) > upper
            badness += ((diff * np.sqrt(invvar)) > 0) * qbad
    #
    # Decide how bad a point is according to maxdev.
    #
    if maxdev is not None:
        qbad = np.absolute(diff) > maxdev
        badness += np.absolute(diff) / maxdev * qbad
    #
    # Do not consider rejecting points that are already rejected by inmask.
    # Do not consider rejecting points that are already rejected by outmask,
    # if sticky is set.
    #
    if inmask is not None:
        badness *= inmask
    if sticky:
        badness *= outmask
    #
    # Reject a maximum of maxrej (additional) points in all the data, or
    # in each group as specified by groupsize, and optionally along each
    # dimension specified by groupdim.
    #
    if maxrej is not None:
        #
        # Loop over each dimension of groupdim or loop once if not set.
        #
        for iloop in range(max(len(groupdim), 1)):
            #
            # Assign an index number in this dimension to each data point.
            #
            if len(groupdim) > 0:
                yndim = len(data.shape)
                if groupdim[iloop] > yndim:
                    raise ValueError('groupdim is larger than the number of dimensions for ydata.')
                dimnum = djs_laxisnum(data.shape, iaxis=groupdim[iloop]-1)
            else:
                dimnum = np.asarray([0])
            #
            # Loop over each vector specified by groupdim. For example, if
            # this is a 2-D array with groupdim=1, then loop over each
            # column of the data.  If groupdim=2, then loop over each row.
            # If groupdim is not set, then use the whole image.
            #
            for ivec in range(max(dimnum)):
                #
                # At this point it is not possible that dimnum is not set.
                #
                indx = (dimnum == ivec).nonzero()[0]
                #
                # Within this group of points, break it down into groups
                # of points specified by groupsize, if set.
                #
                nin = len(indx)
                if groupbadpix:
                    goodtemp = badness == 0
                    groups_lower = (-1*np.diff(np.insert(goodtemp, 0, 1)) == 1).nonzero()[0]
                    groups_upper = (np.diff(np.append(goodtemp, 1)) == 1).nonzero()[0]
                    ngroups = len(groups_lower)
                else:
                    #
                    # The IDL version of this test makes no sense because
                    # groupsize will always be set.
                    #
                    if 'groupsize' is not None:
                        ngroups = nin/groupsize + 1
                        groups_lower = np.arange(ngroups, dtype='i4')*groupsize
                        foo = (np.arange(ngroups, dtype='i4')+1)*groupsize
                        groups_upper = np.where(foo < nin, foo, nin) - 1
                    else:
                        ngroups = 1
                        groups_lower = [0, ]
                        groups_upper = [nin - 1, ]
                for igroup in range(ngroups):
                    i1 = groups_lower[igroup]
                    i2 = groups_upper[igroup]
                    nii = i2 - i1 + 1
                    #
                    # Need the test that i1 != -1 below to prevent a crash
                    # condition, but why is it that we ever get groups
                    # without any points?  Because this is badly-written,
                    # that's why.
                    #
                    if nii > 0 and i1 != -1:
                        jj = indx[i1:i2+1]
                        #
                        # Test if too many points rejected in this group.
                        #
                        if np.sum(badness[jj] != 0) > maxrej[iloop]:
                            isort = badness[jj].argsort()
                            #
                            # Make the following points good again.
                            #
                            badness[jj[isort[0:nii-maxrej[iloop]]]] = 0
                        i1 += groupsize[iloop]
    #
    # Now modify outmask, rejecting points specified by inmask=0, outmask=0
    # if sticky is set, or badness > 0.
    #
    # print(badness)
    newmask = badness == 0
    # print(newmask)
    if grow > 0:
        rejects = newmask == 0
        if rejects.any():
            irejects = rejects.nonzero()[0]
            for k in range(1, grow):
                newmask[(irejects - k) > 0] = 0
                newmask[(irejects + k) < (data.shape[0]-1)] = 0
    if inmask is not None:
        newmask = newmask & inmask
    if sticky:
        newmask = newmask & outmask
    #
    # Set qdone if the input outmask is identical to the output outmask.
    #
    qdone = bool(np.all(newmask == outmask)) # This needs to be a python (rather than a numpy) boolean to avoid painful problems with comparing
                                             # to python True and False cause problems
    outmask = newmask
    return (outmask, qdone)


def djs_laxisnum(dims, iaxis=0):
    """Returns an integer array where each element of the array is set equal
    to its index number in the specified axis.

    Parameters
    ----------
    dims : :class:`list`
        Dimensions of the array to return.
    iaxis : :class:`int`, optional
        Index along this dimension.

    Returns
    -------
    :class:`numpy.ndarray`
        An array of indexes with ``dtype=int32``.

    Raises
    ------
    ValueError
        If `iaxis` is greater than or equal to the number of dimensions, or
        if number of dimensions is greater than three.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine
    and :func:`~pydl.pydlutils.misc.djs_laxisgen`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisnum
    >>> print(djs_laxisnum([4,4]))
    [[0 0 0 0]
     [1 1 1 1]
     [2 2 2 2]
     [3 3 3 3]]
    """
    ndimen = len(dims)
    result = np.zeros(dims, dtype='i4')
    if ndimen == 1:
        pass
    elif ndimen == 2:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k, :] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:, k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    elif ndimen == 3:
        if iaxis == 0:
            for k in range(dims[0]):
                result[k, :, :] = k
        elif iaxis == 1:
            for k in range(dims[1]):
                result[:, k, :] = k
        elif iaxis == 2:
            for k in range(dims[2]):
                result[:, :, k] = k
        else:
            raise ValueError("Bad value for iaxis: {0:d}".format(iaxis))
    else:
        raise ValueError("{0:d} dimensions not supported.".format(ndimen))
    return result

