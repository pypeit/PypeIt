# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit import bspline
from pypeit.core import basis

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
    :class:`numpy.ndarray`:
        The `yval` array with masked values replaced by interpolated values.
    """
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



#class bspline(object):
#    """Bspline class.
#
#    Functions in the bspline library are implemented as methods on this
#    class.
#
#    Parameters
#    ----------
#    x : :class:`numpy.ndarray`
#        The data.
#    nord : :class:`int`, optional
#        To be documented.
#    npoly : :class:`int`, optional
#        To be documented.
#    bkpt : :class:`numpy.ndarray`, optional
#        To be documented.
#    bkspread : :class:`float`, optional
#        To be documented.
#    verbose : :class:`bool`, optional.
#        If ``True`` print extra information.
#
#    Attributes
#    ----------
#    breakpoints
#        Breakpoints for bspline, spacing for these breakpoints are determinated by keywords inputs;
#    nord
#        Order of bspline; [default=4]
#    npoly
#        Polynomial order to fit over 2nd variable (when specified as x2): [default=1]
#    mask
#        Output mask, set =1 for good points, =0 for bad points;
#    coeff
#        Output coefficient of the bspline;
#    icoeff
#        Cholesky band matrix used to solve for the bspline coefficients;
#    xmin
#        Normalization minimum for x2; [default max(xdata)]
#    xmax
#        Normalization maximum for x2; [default min(xdata)]
#    funcname
#        Function for the second variable; [default 'legendre']
#    from_dict
#        If not None, create a bspline from a dictionary created by to_dict(). [default 'None']
#        It is possible to instantiate a bspline from a dict without the x data:
#        new_bspline = bspline(None, from_dict=dictionary)
#    """
#
#    # ToDO Consider refactoring the argument list so that there are no kwargs
#    def __init__(self, x, fullbkpt = None, nord=4, npoly=1, bkpt=None, bkspread=1.0,
#                 verbose=False, from_dict=None, **kwargs):
#        """Init creates an object whose attributes are similar to the
#        structure returned by the create_bspline function.
#        """
#        # JFH added this to enforce immutability of these input arguments, as this code modifies bkpt and fullbkpt
#        # as it goes
#        fullbkpt1 = copy.copy(fullbkpt)
#        bkpt1 = copy.copy(bkpt)
#        if from_dict is not None:
#            self.nord=from_dict['nord']
#            self.npoly=from_dict['npoly']
#            self.breakpoints=np.array(from_dict['breakpoints'])
#            self.mask=np.array(from_dict['mask'])
#            self.coeff=np.array(from_dict['coeff'])
#            self.icoeff=np.array(from_dict['icoeff'])
#            self.xmin=from_dict['xmin']
#            self.xmax=from_dict['xmax']
#            self.funcname=from_dict['funcname']
#            return
#        # Instantiate empty if neither fullbkpt or x is set
#        elif x is None and fullbkpt is None:
#            self.nord = None
#            self.npoly = None
#            self.breakpoints= None
#            self.mask= None
#            self.coeff= None
#            self.icoeff= None
#            self.xmin= None
#            self.xmax= None
#            self.funcname= None
#            return
#        else:
#            #
#            # Set the breakpoints.
#            #
#            if fullbkpt1 is None:
#                if bkpt1 is None:
#                    startx = x.min()
#                    rangex = x.max() - startx
#                    if 'placed' in kwargs:
#                        w = ((kwargs['placed'] >= startx) &
#                             (kwargs['placed'] <= startx+rangex))
#                        if w.sum() < 2:
#                            bkpt1 = np.arange(2, dtype='f') * rangex + startx
#                        else:
#                            bkpt1 = kwargs['placed'][w]
#                    elif 'bkspace' in kwargs:
#                        nbkpts = int(rangex/kwargs['bkspace']) + 1
#                        if nbkpts < 2:
#                            nbkpts = 2
#                        tempbkspace = rangex/float(nbkpts-1)
#                        bkpt1 = np.arange(nbkpts, dtype='f')*tempbkspace + startx
#                    elif 'nbkpts' in kwargs:
#                        nbkpts = kwargs['nbkpts']
#                        if nbkpts < 2:
#                            nbkpts = 2
#                        tempbkspace = rangex/float(nbkpts-1)
#                        bkpt1 = np.arange(nbkpts, dtype='f') * tempbkspace + startx
#                    elif 'everyn' in kwargs:
#                        nx = x.size
#                        nbkpts = max(nx/kwargs['everyn'], 1)
#                        if nbkpts == 1:
#                            xspot = [0]
#                        else:
#                            xspot = (nx/nbkpts)*np.arange(nbkpts)
#                            # JFH This was a bug. Made fixes
#                            #xspot = int(nx/(nbkpts-1)) * np.arange(nbkpts, dtype='i4')
#                        #bkpt = x[xspot].astype('f')
#                        bkpt1 = np.interp(xspot,np.arange(nx),x)
#                    else:
#                        raise ValueError('No information for bkpts.')
#                # JFH added this new code, because bkpt.size = 1 implies fullbkpt has only 2*(nord-1) + 1 elements.
#                # This will cause a crash in action because nbkpt < 2*nord, i.e. for bkpt = 1, nord = 4 fullbkpt has
#                # seven elements which is less than 2*nord = 8. The codes above seem to require nbkpt >=2, so I'm implementing
#                # this requirement. Note that the previous code before this fix simply sets bkpt to bkpt[imax] =x.max()
#                # which is equally arbitrary, but still results in a crash. By requiring at least 2 bkpt, fullbkpt will
#                # have 8 elements preventing action from crashing
#                if (bkpt1.size < 2):
#                    bkpt1 = np.zeros(2,dtype=float)
#                    bkpt1[0] = x.min()
#                    bkpt1[1] = x.max()
#                else:
#                    imin = bkpt1.argmin()
#                    imax = bkpt1.argmax()
#                    if x.min() < bkpt1[imin]:
#                        if verbose:
#                            print('Lowest breakpoint does not cover lowest x value: changing.')
#                        bkpt1[imin] = x.min()
#                    if x.max() > bkpt1[imax]:
#                        if verbose:
#                            print('Highest breakpoint does not cover highest x value: changing.')
#                        bkpt1[imax] = x.max()
#
#                nshortbkpt = bkpt1.size
#                fullbkpt1 = bkpt1.copy()
#                # Note that with the JFH change above, this nshortbkpt ==1 is never realized beacause above I forced
#                # bkpt to have at least two elements. Not sure why this was even allowed, since bkpt.size = 1
#                #  basically results in action crashing as described above.
#                if nshortbkpt == 1:
#                    bkspace = bkspread
#                else:
#                    bkspace = (bkpt1[1] - bkpt1[0])*bkspread
#                for i in np.arange(1, nord):
#                    fullbkpt1 = np.insert(fullbkpt1, 0, bkpt1[0]-bkspace*i)
#                    fullbkpt1 = np.insert(fullbkpt1, fullbkpt1.shape[0],
#                                         bkpt1[nshortbkpt-1] + bkspace*i)
#
#
#            # JFH added this to fix bug in cases where fullbkpt is passed in but has < 2*nord elements
#            if fullbkpt1.size < 2*nord:
#                fullbkpt_init = fullbkpt1.copy()
#                nshortbkpt = fullbkpt_init.size
#                bkspace = (fullbkpt_init[1] - fullbkpt_init[0])*bkspread
#                for i in np.arange(1, nord):
#                    fullbkpt1 = np.insert(fullbkpt1, 0, fullbkpt_init[0] - bkspace * i)
#                    fullbkpt1 = np.insert(fullbkpt1, fullbkpt1.shape[0],
#                                          fullbkpt_init[nshortbkpt - 1] + bkspace * i)
#
#            nc = fullbkpt1.size - nord
#            self.breakpoints = fullbkpt1
#            self.nord = nord
#            self.npoly = npoly
#            self.mask = np.ones((fullbkpt1.size,), dtype='bool')
#            if npoly > 1:
#                self.coeff = np.zeros((npoly, nc), dtype='d')
#                self.icoeff = np.zeros((npoly, nc), dtype='d')
#            else:
#                self.coeff = np.zeros((nc,), dtype='d')
#                self.icoeff = np.zeros((nc,), dtype='d')
#            self.xmin = 0.0
#            self.xmax = 1.0
#            if 'funcname' in kwargs:
#                self.funcname = kwargs['funcname']
#            else:
#                self.funcname = 'legendre'
#
#        return
#
#    def copy(self):
#
#        bsp_copy = bspline(None)
#        bsp_copy.nord = self.nord
#        bsp_copy.npoly = self.npoly
#        bsp_copy.breakpoints = np.copy(self.breakpoints)
#        bsp_copy.mask = np.copy(self.mask)
#        bsp_copy.coeff = np.copy(self.coeff)
#        bsp_copy.icoeff = np.copy(self.icoeff)
#        bsp_copy.xmin = self.xmin
#        bsp_copy.xmax = self.xmax
#        bsp_copy.funcname = self.funcname
#        return bsp_copy
#
#    def to_dict(self):
#        """Write bspline parameters to a dict.
#
#        Parameters
#        ----------
#
#        Returns
#        -------
#            A dict containing the relevant bspline paramenters.
#
#        Notes
#        -----
#        The dictionary is JSON compatible.
#        """
#
#        # needs to move np.arrays to lists for JSON files
#        return (dict(breakpoints=self.breakpoints.tolist(),
#                     nord=self.nord,
#                     npoly=self.npoly,
#                     mask=self.mask.tolist(),
#                     coeff=self.coeff.tolist(),
#                     icoeff=self.icoeff.tolist(),
#                     xmin=self.xmin,
#                     xmax=self.xmax,
#                     funcname=self.funcname))
#
#    def fit(self, xdata, ydata, invvar, x2=None):
#        """Calculate a B-spline in the least-squares sense.
#
#        Fit is based on two variables: x which is sorted and spans a large range
#        where bkpts are required y which can be described with a low order
#        polynomial.
#
#        Parameters
#        ----------
#        xdata : :class:`numpy.ndarray`
#            Independent variable.
#        ydata : :class:`numpy.ndarray`
#            Dependent variable.
#        invvar : :class:`numpy.ndarray`
#            Inverse variance of `ydata`.
#        x2 : :class:`numpy.ndarray`, optional
#            Orthogonal dependent variable for 2d fits.
#
#        Returns
#        -------
#        :func:`tuple`
#            A tuple containing an integer error code, and the evaluation of the
#            b-spline at the input values.  An error code of -2 is a failure,
#            -1 indicates dropped breakpoints, 0 is success, and positive
#            integers indicate ill-conditioned breakpoints.
#        """
#        goodbk = self.mask[self.nord:]
#        nn = goodbk.sum()
#        if nn < self.nord:
#            yfit = np.zeros(ydata.shape, dtype='f')
#            return (-2, yfit)
#        nfull = nn * self.npoly
#        bw = self.npoly * self.nord
#        a1, lower, upper = self.action(xdata, x2=x2)
#        foo = np.tile(invvar, bw).reshape(bw, invvar.size).transpose()
#        a2 = a1 * foo
#        alpha = np.zeros((bw, nfull+bw), dtype='d')
#        beta = np.zeros((nfull+bw,), dtype='d')
#        bi = np.arange(bw, dtype='i4')
#        bo = np.arange(bw, dtype='i4')
#        for k in range(1, bw):
#            bi = np.append(bi, np.arange(bw-k, dtype='i4')+(bw+1)*k)
#            bo = np.append(bo, np.arange(bw-k, dtype='i4')+bw*k)
#        for k in range(nn-self.nord+1):
#            itop = k*self.npoly
#            ibottom = min(itop, nfull) + bw - 1
#            ict = upper[k] - lower[k] + 1
#            if ict > 0:
#                work = np.dot(a1[lower[k]:upper[k]+1, :].T, a2[lower[k]:upper[k]+1, :])
#                wb = np.dot(ydata[lower[k]:upper[k]+1], a2[lower[k]:upper[k]+1, :])
#                alpha.T.flat[bo+itop*bw] += work.flat[bi]
#                beta[itop:ibottom+1] += wb
#        min_influence = 1.0e-10 * invvar.sum() / nfull
#        errb = cholesky_band(alpha, mininf=min_influence)  # ,verbose=True)
#        if isinstance(errb[0], int) and errb[0] == -1:
#            a = errb[1]
#        else:
#            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
#            return (self.maskpoints(errb[0]), yfit)
#        errs = cholesky_solve(a, beta)
#        if isinstance(errs[0], int) and errs[0] == -1:
#            sol = errs[1]
#        else:
#            #
#            # It is not possible for this to get called, because cholesky_solve
#            # has only one return statement, & that statement guarantees that
#            # errs[0] == -1
#            #
#            yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
#            return (self.maskpoints(errs[0]), yfit)
#        if self.coeff.ndim == 2:
#            # JFH made major bug fix here.
#            self.icoeff[:, goodbk] = np.array(a[0, 0:nfull].T.reshape(self.npoly, nn,order='F'), dtype=a.dtype)
#            self.coeff[:, goodbk] = np.array(sol[0:nfull].T.reshape(self.npoly, nn, order='F'), dtype=sol.dtype)
#        else:
#            self.icoeff[goodbk] = np.array(a[0, 0:nfull], dtype=a.dtype)
#            self.coeff[goodbk] = np.array(sol[0:nfull], dtype=sol.dtype)
#        yfit, foo = self.value(xdata, x2=x2, action=a1, upper=upper, lower=lower)
#        return (0, yfit)
#
#    def action(self, x, x2=None):
#        """Construct banded bspline matrix, with dimensions [ndata, bandwidth].
#
#        Parameters
#        ----------
#        x : :class:`numpy.ndarray`
#            Independent variable.
#        x2 : :class:`numpy.ndarray`, optional
#            Orthogonal dependent variable for 2d fits.
#
#        Returns
#        -------
#        :func:`tuple`
#            A tuple containing the b-spline action matrix; the 'lower' parameter,
#            a list of pixel positions, each corresponding to the first
#            occurence of position greater than breakpoint indx; and 'upper',
#            Same as lower, but denotes the upper pixel positions.
#        """
#        nbkpt = self.mask.sum()
#        if nbkpt < 2*self.nord:
#            msgs.warn('Order ({0}) too low for {1} breakpoints.'.format(self.nord, nbkpt))
#            return -2, 0, 0
#        nx = x.size
#        n = nbkpt - self.nord
#        lower = np.zeros((n - self.nord + 1,), dtype=int)
#        upper = np.zeros((n - self.nord + 1,), dtype=int) - 1
#        indx = self.intrv(x)
#        bf1 = self.bsplvn(x, indx)
#        aa = uniq(indx)
#        upper[indx[aa]-self.nord+1] = aa
#        rindx = indx[::-1]
#        bb = uniq(rindx)
#        lower[rindx[bb]-self.nord+1] = nx - bb - 1
#        if x2 is None:
#            return bf1, lower, upper
#
#        if x2.size != nx:
#            raise ValueError('Dimensions of x and x2 do not match.')
#
#        # TODO: Below is unchanged.
#        x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
#        # TODO: Should consider faster ways of generating the temppoly arrays for poly and poly1
#        if self.funcname == 'poly':
#            temppoly = np.ones((nx, self.npoly), dtype='f')
#            for i in range(1, self.npoly):
#                temppoly[:, i] = temppoly[:, i-1] * x2norm
#        elif self.funcname == 'poly1':
#            temppoly = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
#            for i in range(1, self.npoly):
#                temppoly[:, i] = temppoly[:, i-1] * x2norm
#        elif self.funcname == 'chebyshev':
#            # JFH fixed bug here where temppoly needed to be transposed because of different IDL and python array conventions
#            temppoly = fchebyshev(x2norm, self.npoly).T
#        elif self.funcname == 'legendre':
#            temppoly = flegendre(x2norm, self.npoly).T
#        else:
#            raise ValueError('Unknown value of funcname.')
#
#        # TODO: Should consider faster way of calculating action that doesn't require a nested loop.
#        bw = self.npoly*self.nord
#        action = np.zeros((nx, bw), dtype='d')
#        counter = -1
#        for ii in range(self.nord):
#            for jj in range(self.npoly):
#                counter += 1
#                action[:, counter] = bf1[:, ii]*temppoly[:, jj]
#        return action, lower, upper
#
#    def intrv(self, x):
#        """Find the segment between breakpoints which contain each value in the array x.
#
#        The minimum breakpoint is nbkptord -1, and the maximum
#        is nbkpt - nbkptord - 1.
#
#        Parameters
#        ----------
#        x : :class:`numpy.ndarray`
#            Data values, assumed to be monotonically increasing.
#
#        Returns
#        -------
#        :class:`numpy.ndarray`
#            Position of array elements with respect to breakpoints.
#        """
#        gb = self.breakpoints[self.mask]
#        n = gb.size - self.nord
#        indx = np.zeros(x.size, dtype=int)
#        ileft = self.nord - 1
#        for i in range(x.size):
#            while x[i] > gb[ileft+1] and ileft < n - 1:
#                ileft += 1
#            indx[i] = ileft
#        return indx
#
#    def bsplvn(self, x, ileft):
#        """To be documented.
#
#        Parameters
#        ----------
#        x : :class:`numpy.ndarray`
#            To be documented.
#        ileft : :class:`int`
#            To be documented
#
#        Returns
#        -------
#        :class:`numpy.ndarray`
#            To be documented.
#        """
#        bkpt = self.breakpoints[self.mask]
#        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype)
#        deltap = vnikx.copy()
#        deltam = vnikx.copy()
#        j = 0
#        vnikx[:, 0] = 1.0
#        while j < self.nord - 1:
#            ipj = ileft+j+1
#            deltap[:, j] = bkpt[ipj] - x
#            imj = ileft-j
#            deltam[:, j] = x - bkpt[imj]
#            vmprev = 0.0
#            for l in range(j+1):
#                vm = vnikx[:, l]/(deltap[:, l] + deltam[:, j-l])
#                vnikx[:, l] = vm*deltap[:, l] + vmprev
#                vmprev = vm*deltam[:, j-l]
#            j += 1
#            vnikx[:, j] = vmprev
#        return vnikx
#
#    def value(self, x, x2=None, action=None, lower=None, upper=None):
#        """Evaluate a bspline at specified values.
#
#        Parameters
#        ----------
#        x : :class:`numpy.ndarray`
#            Independent variable.
#        x2 : :class:`numpy.ndarray`, optional
#            Orthogonal dependent variable for 2d fits.
#        action : :class:`numpy.ndarray`, optional
#            Action matrix to use.  If not supplied it is calculated.
#        lower : :class:`numpy.ndarray`, optional
#            If the action parameter is supplied, this parameter must also
#            be supplied.
#        upper : :class:`numpy.ndarray`, optional
#            If the action parameter is supplied, this parameter must also
#            be supplied.
#
#        Returns
#        -------
#        :func:`tuple`
#            A tuple containing the results of the bspline evaluation and a
#            mask indicating where the evaluation was good.
#        """
#        # TODO: Is the sorting necessary?
#        xsort = x.argsort()
#        if action is None:
#            action, lower, upper = self.action(x[xsort], x2=None if x2 is None else x2[xsort])
#        else:
#            if lower is None or upper is None:
#                raise ValueError('Must specify lower and upper if action is set.')
#
#        # TODO: Can we save some of these objects to self so that we
#        # don't have to recreate them?
#        yfit = np.zeros(x.shape, dtype=x.dtype)
#        bw = self.npoly * self.nord
#        spot = np.arange(bw, dtype=int)
#        goodbk = self.mask.nonzero()[0]
#        coeffbk = self.mask[self.nord:].nonzero()[0]
#        goodcoeff = self.coeff[...,coeffbk]
#
#        nowidth = np.invert(upper+1 > lower)
#        n = self.mask.sum() - self.nord
#        for i in range(n-self.nord+1):
#            if nowidth[i]:
#                continue
#            yfit[lower[i]:upper[i]+1] = np.dot(action[lower[i]:upper[i]+1,:],
#                                               goodcoeff.flatten('F')[i*self.npoly+spot])
#
#        mask = np.ones(x.shape, dtype='bool')
#        gb = self.breakpoints[goodbk]
#        mask[(x < gb[self.nord-1]) | (x > gb[n])] = False
#        hmm = (np.diff(goodbk) > 2).nonzero()[0]
#        if hmm.size == 0:
#            return yfit[np.argsort(xsort)], mask
#
#        for jj in range(hmm.size):
#            mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
#                    & (x <= self.breakpoints[goodbk[hmm[jj]+1]-1])] = False
#        return yfit[np.argsort(xsort)], mask
#
##    def value(self, x, x2=None, action=None, lower=None, upper=None):
##        """Evaluate a bspline at specified values.
##
##        Parameters
##        ----------
##        x : :class:`numpy.ndarray`
##            Independent variable.
##        x2 : :class:`numpy.ndarray`, optional
##            Orthogonal dependent variable for 2d fits.
##        action : :class:`numpy.ndarray`, optional
##            Action matrix to use.  If not supplied it is calculated.
##        lower : :class:`numpy.ndarray`, optional
##            If the action parameter is supplied, this parameter must also
##            be supplied.
##        upper : :class:`numpy.ndarray`, optional
##            If the action parameter is supplied, this parameter must also
##            be supplied.
##
##        Returns
##        -------
##        :func:`tuple`
##            A tuple containing the results of the bspline evaluation and a
##            mask indicating where the evaluation was good.
##        """
##        # TODO: Is the sorting necessary?
##        xsort = x.argsort()
##        if action is None:
##            action, lower, upper = self.action(x[xsort], x2=None if x2 is None else x2[xsort])
##        else:
##            if lower is None or upper is None:
##                raise ValueError('Must specify lower and upper if action is set.')
##
##        # TODO: Can we save some of these objects to self so that we
##        # don't have to recreate them?
##        yfit = np.zeros(x.shape, dtype=x.dtype)
##        bw = self.npoly * self.nord
##        spot = np.arange(bw, dtype=int)
##        goodbk = self.mask.nonzero()[0]
##        coeffbk = self.mask[self.nord:].nonzero()[0]
##        goodcoeff = self.coeff[...,coeffbk]
##
##        # maskthis = np.zeros(xwork.shape,dtype=xwork.dtype)
##        nowidth = np.invert(upper > lower)
##        n = self.mask.sum() - self.nord
##        for i in range(n-self.nord+1):
##            if nowidth[i]:
##                continue
##            yfit[lower[i]:upper[i]+1] = np.dot(action[lower[i]:upper[i]+1,:],
##                                               goodcoeff.flatten('F')[i*self.npoly+spot])
##
##        mask = np.ones(x.shape, dtype='bool')
##        gb = self.breakpoints[goodbk]
##        mask[(x < gb[self.nord-1]) | (x > gb[n])] = False
##        hmm = (np.diff(goodbk) > 2).nonzero()[0]
##        if hmm.size == 0:
##            return yfit[np.argsort(xsort)], mask
##
##        for jj in range(hmm.size):
##            mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
##                    & (x <= self.breakpoints[goodbk[hmm[jj]+1]-1])] = False
##        return yfit[np.argsort(xsort)], mask
#
#    def maskpoints(self, err):
#        """Perform simple logic of which breakpoints to mask.
#
#
#        Parameters
#        ----------
#        err : :class:`numpy.ndarray` or int
#            The list of indexes returned by the cholesky routines.
#            This is indexed to the set of currently *good* breakpoints (i.e. self.mask=True)
#            And the first nord are skipped
#
#        Returns
#        -------
#        :class:`int`
#            An integer indicating the results of the masking.  -1 indicates
#            that the error points were successfully masked.  -2 indicates
#            failure; the calculation should be aborted.
#
#        Notes
#        -----
#        The mask attribute is modified, assuming it is possible to create the
#        mask.
#        """
#        # Recast err as an array if a single value int was passed in (occasional)
#        if not isinstance(err, np.ndarray):
#            err = np.array([err])
#        # Currently good points
#        goodbkpt = np.where(self.mask)[0]
#        nbkpt = len(goodbkpt)
#        if nbkpt <= 2*self.nord:
#            msgs.warn('Fewer good break points than order of b-spline. Returning...')
#            return -2
#        # Find the unique ones for the polynomial
#        hmm = err[uniq(err//self.npoly)]//self.npoly
#
#        n = nbkpt - self.nord
#        if np.any(hmm >= n):
#            msgs.warn('Note enough unique points in cholesky_band decomposition of b-spline matrix. Returning...')
#            return -2
#        test = np.zeros(nbkpt, dtype='bool')
#        for jj in range(-int(np.ceil(self.nord/2)), int(self.nord/2.)):
#            foo = np.where((hmm+jj) > 0, hmm+jj, np.zeros(hmm.shape, dtype=hmm.dtype))
#            inside = np.where((foo+self.nord) < n-1, foo+self.nord, np.zeros(hmm.shape, dtype=hmm.dtype)+n-1)
#            if len(inside)>0:
#                test[inside] = True
#        if test.any():
#            reality = goodbkpt[test]
#            if self.mask[reality].any():
#                self.mask[reality] = False
#                return -1
#            else:
#                return -2
#        else:
#            return -2
#
#    def workit(self, xdata, ydata, invvar, action, lower, upper):
#        """An internal routine for bspline_extract and bspline_radial which solve a general
#        banded correlation matrix which is represented by the variable "action".  This routine
#        only solves the linear system once, and stores the coefficients in sset. A non-zero return value
#        signifies a failed inversion
#
#
#        Parameters
#        ----------
#        xdata : :class:`numpy.ndarray`
#            Independent variable.
#        ydata : :class:`numpy.ndarray`
#            Dependent variable.
#        invvar : :class:`numpy.ndarray`
#            Inverse variance of `ydata`.
#        action : :class:`numpy.ndarray`
#            Banded correlation matrix
#        lower  : :class:`numpy.ndarray`
#            A list of pixel positions, each corresponding to the first occurence of position greater than breakpoint indx
#        upper  : :class:`numpy.ndarray`
#            Same as lower, but denotes the upper pixel positions
#
#        Returns
#        -------
#        :func:`tuple` (success, yfit):
#            A tuple containing an boolean error code, and the evaluation
#            of the b-spline yfit at the input values.  The error codes
#            are as follows: 0 is good; -1 is dropped breakpoints, try
#            again; -2 is failure, should abort.
#
#        """
#        goodbk = self.mask[self.nord:]
#        # KBW: Interesting: x.sum() is actually a bit faster than np.sum(x)
#        nn = goodbk.sum()
#        if nn < self.nord:
#            msgs.warn('Fewer good break points than order of b-spline. Returning...')
#            # KBW: Why is the dtype set to 'f' = np.float32?
#            return -2, np.zeros(ydata.shape, dtype=np.float32)
#
#        nfull = nn * self.npoly
#        bw = self.npoly * self.nord
#        a2 = action * np.sqrt(invvar)[:,None]
#
#        alpha = np.zeros((bw, nfull+bw), dtype=float)
#        beta = np.zeros((nfull+bw,), dtype=float)
#        bi = np.concatenate([np.arange(i)+(bw-i)*(bw+1) for i in range(bw,0,-1)])
#        bo = np.concatenate([np.arange(i)+(bw-i)*bw for i in range(bw,0,-1)])
#        upper += 1
#        nowidth = np.invert(upper > lower)
#        for k in range(nn-self.nord+1):
#            if nowidth[k]:
#                continue
#            itop = k*self.npoly
#            alpha.T.flat[bo+itop*bw] \
#                    += np.dot(a2[lower[k]:upper[k],:].T, a2[lower[k]:upper[k],:]).flat[bi]
#            beta[itop:min(itop,nfull)+bw] \
#                    += np.dot(ydata[lower[k]:upper[k]] * np.sqrt(invvar[lower[k]:upper[k]]),
#                              a2[lower[k]:upper[k],:])
#        upper -= 1
#
#        # Right now we are not returning the covariance, although it may arise that we should
##        covariance = alpha
#        err, a = cholesky_band(alpha, mininf=1.0e-10 * invvar.sum() / nfull)
#
#        # successful cholseky_band returns -1
#        if not isinstance(err, int) or err != -1:
#            return self.maskpoints(err), \
#                        self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)[0]
#
#        # NOTE: cholesky_solve ALWAYS returns err == -1; don't even catch it.
#        sol = cholesky_solve(a, beta)[1]
#
#        if self.coeff.ndim == 2:
#            self.icoeff[:,goodbk] = np.array(a[0,:nfull].T.reshape(self.npoly, nn, order='F'), dtype=a.dtype)
#            self.coeff[:,goodbk] = np.array(sol[:nfull].T.reshape(self.npoly, nn, order='F'), dtype=sol.dtype)
#        else:
#            self.icoeff[goodbk] = np.array(a[0,:nfull], dtype=a.dtype)
#            self.coeff[goodbk] = np.array(sol[:nfull], dtype=sol.dtype)
#
#        return 0, self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)[0]
#
#
#
#def cholesky_band(l, mininf=0.0):
#    """Compute Cholesky decomposition of banded matrix.
#
#    Parameters
#    ----------
#    l : :class:`numpy.ndarray`
#        A matrix on which to perform the Cholesky decomposition.
#    mininf : :class:`float`, optional
#        Entries in the `l` matrix are considered negative if they are less
#        than this value (default 0.0).
#
#    Returns
#    -------
#    :func:`tuple`
#        If problems were detected, the first item will be the index or
#        indexes where the problem was detected, and the second item will simply
#        be the input matrix.  If no problems were detected, the first item
#        will be -1, and the second item will be the Cholesky decomposition.
#    """
#
#    # KBW: added isfinite check so that it doesn't need to be done in
#    # the for loop. This should be okay because sqrt and division of
#    # positive numbers by other positive numbers should never lead to
#    # non-finite numbers. However, need to watch out for numerical
#    # issues.
#    bw, nn = l.shape
#    n = nn - bw
#    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
#    # JFH changed this below to make it more consistent with IDL version. Not sure
#    # why the np.all(np.isfinite(lower)) was added. The code could return an empty
#    # list for negative.nonzero() and crash if all elements in lower are NaN.
#    if negative.any():
#        nz = negative.nonzero()[0]
#        msgs.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
#        return nz, l
#
##    negative = (lower[0, 0:n] <= mininf)
##    if negative.any() or not np.all(np.isfinite(lower)):
##        msgs.warn('Found {:d}'.format(len(negative.nonzero()[0])) +
##                  ' bad entries: ' + str(negative.nonzero()[0]))
##        return (negative.nonzero()[0], l)
#
#    # KBW: Moved the copy to after the initial check of the input
#    lower = l.copy()
#    kn = bw - 1
#    spot = np.arange(kn, dtype=int) + 1
#
#    # KBW: Faster by about factor of ~2 compared to previous version
#    #bi = np.arange(kn*kn).reshape(kn,kn)[np.triu_indices(kn)]
#    bi = np.concatenate([np.arange(i)+(kn-i)*(kn+1) for i in range(kn,0,-1)])
#    here = bi[:,None] + (np.arange(n)[None,:] + 1)*bw
#    for j in range(n):
#        lower[0,j] = np.sqrt(lower[0,j])
#        lower[spot, j] /= lower[0,j]
#        hmm = lower[spot,j,None] * lower[None,spot,j]
#        lower.T.flat[here[:,j]] -= hmm.flat[bi]
#
#    return -1, lower
#
#
#def cholesky_solve(a, bb):
#    """Solve the equation Ax=b where A is a Cholesky-banded matrix.
#
#    Parameters
#    ----------
#    a : :class:`numpy.ndarray`
#        :math:`A` in :math:`A x = b`.
#    bb : :class:`numpy.ndarray`
#        :math:`b` in :math:`A x = b`.
#
#    Returns
#    -------
#    :func:`tuple`
#        A tuple containing the status and the result of the solution.  The
#        status is always -1.
#    """
#    b = bb.copy()
#    n = b.shape[0] - a.shape[0]
#    kn = a.shape[0] - 1
#
#    spot = np.arange(kn, dtype=int) + 1
#    for j in range(n):
#        b[j] /= a[0,j]
#        b[j+spot] -= b[j]*a[spot,j]
#
#    spot = spot[::-1]
#    for j in range(n-1, -1, -1):
#        b[j] = (b[j] - np.sum(a[spot,j] * b[j+spot]))/a[0,j]
#
#    return -1, b


def iterfit(xdata, ydata, invvar=None, inmask = None, upper=5, lower=5, x2=None,
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

    if inmask is None:
        inmask = invvar > 0.0

    if x2 is not None:
        if x2.size != nx:
            raise ValueError('Dimensions of xdata and x2 do not agree.')
    yfit = np.zeros(ydata.shape)
    if invvar.size == 1:
        outmask = True
    else:
        outmask = np.ones(invvar.shape, dtype='bool')
    xsort = xdata.argsort()
    maskwork = (outmask & inmask & (invvar > 0.0))[xsort]
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
            sset = bspline.bspline(xdata[xsort[maskwork]], nord=nord, bkpt=bkpt, fullbkpt=fullbkpt,
                                   **kwargs_bspline)
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
        inmask_rej = maskwork
        if error == -2:

            return (sset, outmask)
        elif error == 0:
            # ToDO JFH by setting inmask to be tempin which is maskwork, we are basically implicitly enforcing sticky rejection
            # here. See djs_reject.py. I'm leaving this as is for consistency with the IDL version, but this may require
            # further consideration. I think requiring stick to be set is the more transparent behavior.
            maskwork, qdone = djs_reject(ywork, yfit, invvar=invwork,
                                         inmask=inmask_rej, outmask=maskwork,
                                         upper=upper, lower=lower,**kwargs_reject)
        else:
            pass
    outmask[xsort] = maskwork
    temp = yfit
    yfit[xsort] = temp
    return (sset, outmask)

## TODO: How important is it that the function return the last
## occurrence of the unique values in the sorted array than the first
## value? If it doesn't `numpy.unique(numpy.sort(x),
## return_index=True)[1]` is about twice as fast as
## `pydl.uniq(numpy.sort(x))`.
#def old_uniq(x, index=None):
#    """Replicates the IDL ``UNIQ()`` function.
#
#    Returns the *subscripts* of the unique elements of an array.  The elements
#    must actually be *sorted* before being passed to this function.  This can
#    be done by sorting `x` explicitly or by passing the array subscripts that
#    sort `x` as a second parameter.
#
#    Parameters
#    ----------
#    x : array-like
#        Search this array for unique items.
#    index : array-like, optional
#        This array provides the array subscripts that sort `x`.
#
#    Returns
#    -------
#    array-like
#        The subscripts of `x` that are the unique elements of `x`.
#
#    Notes
#    -----
#    Given a sorted array, and assuming that there is a set of
#    adjacent identical items, ``uniq()`` will return the subscript of the
#    *last* unique item.  This charming feature is retained for
#    reproducibility.
#
#    References
#    ----------
#    http://www.harrisgeospatial.com/docs/uniq.html
#
#    Examples
#    --------
#    >>> import numpy as np
#    >>> from pydl import uniq
#    >>> data = np.array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
#    >>> print(uniq(np.sort(data)))
#    [ 3  5  7  9 10 11 12 13]
#    """
#    if index is None:
#        indicies = (x != np.roll(x, -1)).nonzero()[0]
#        if indicies.size > 0:
#            return indicies
#        else:
#            return np.array([x.size - 1, ])
#    else:
#        q = x[index]
#        indicies = (q != np.roll(q, -1)).nonzero()[0]
#        if indicies.size > 0:
#            return index[indicies]
#        else:
#            return np.array([q.size - 1, ], dtype=index.dtype)
#
## Faster than previous version but not as fast as if we could switch to
## np.unique.
#def uniq(x, index=None):
#    """
#    Return the indices of the *last* occurrence of the unique
#    elements in a sorted array.
#
#    The input vector must be sorted before being passed to this
#    function. This can be done by sorting ``x`` directly or by
#    passing the array that sorts ``x`` (``index``).
#
#    Replicates the IDL ``UNIQ()`` function.
#
#    Parameters
#    ----------
#    x : array-like
#        Search this array for unique items.
#    index : array-like, optional
#        This array provides the array subscripts that sort `x`. That
#        is::
#
#            index = np.argsort(x)
#
#    Returns
#    -------
#    `np.ndarray`
#        The indices of the last occurence in `x` of its unique
#        values.
#
#    Notes
#    -----
#    Given a sorted array, and assuming that there is a set of
#    adjacent identical items, ``uniq()`` will return the subscript of
#    the *last* unique item. This charming feature is retained for
#    reproducibility.
#
#    References
#    ----------
#    http://www.harrisgeospatial.com/docs/uniq.html
#
#    Speed improvement thanks to discussion here:
#    https://stackoverflow.com/questions/47495510/numpy-in-a-sorted-list-find-the-first-and-the-last-index-for-each-unique-value
#
#    Examples
#    --------
#    >>> import numpy as np
#    >>> from pydl import uniq
#    >>> data = np.array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
#    >>> print(uniq(np.sort(data)))
#    [ 3  5  7  9 10 11 12 13]
#    """
#    if len(x) == 0:
#        raise ValueError('No unique elements in an empty array!')
#    if index is None:
#        return np.flatnonzero(np.concatenate(([True], x[1:] != x[:-1], [True])))[1:]-1
#    _x = x[index]
#    return np.flatnonzero(np.concatenate(([True], _x[1:] != _x[:-1], [True])))[1:]-1
#
#
#def flegendre(x, m):
#    """Compute the first `m` Legendre polynomials.
#
#    Parameters
#    ----------
#    x : array-like
#        Compute the Legendre polynomials at these abscissa values.
#    m : :class:`int`
#        The number of Legendre polynomials to compute.  For example, if
#        :math:`m = 3`, :math:`P_0 (x)`, :math:`P_1 (x)` and :math:`P_2 (x)`
#        will be computed.
#
#    Returns
#    -------
#    :class:`numpy.ndarray`
#    """
#    from scipy.special import legendre
#    if isinstance(x, np.ndarray):
#        n = x.size
#    else:
#        n = 1
#    if m < 1:
#        raise ValueError('Number of Legendre polynomials must be at least 1.')
#    try:
#        dt = x.dtype
#    except AttributeError:
#        dt = np.float64
#    leg = np.ones((m, n), dtype=dt)
#    if m >= 2:
#        leg[1, :] = x
#    if m >= 3:
#        for k in range(2, m):
#            leg[k, :] = np.polyval(legendre(k), x)
#    return leg
#
#
#def fchebyshev(x, m):
#    """Compute the first `m` Chebyshev polynomials.
#
#    Parameters
#    ----------
#    x : array-like
#        Compute the Chebyshev polynomials at these abscissa values.
#    m : :class:`int`
#        The number of Chebyshev polynomials to compute.  For example, if
#        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and
#        :math:`T_2 (x)` will be computed.
#
#    Returns
#    -------
#    :class:`numpy.ndarray`
#    """
#    from scipy.special import chebyt
#    if isinstance(x, np.ndarray):
#        n = x.size
#    else:
#        n = 1
#    if m < 1:
#        raise ValueError('Order of Chebyshev polynomial must be at least 1.')
#    try:
#        dt = x.dtype
#    except AttributeError:
#        dt = np.float64
#    leg = np.ones((m, n), dtype=dt)
#    if m >= 2:
#        leg[1, :] = x
#    if m >= 3:
#        for k in range(2, m):
#            leg[k, :] = np.polyval(chebyt(k), x)
#    return leg


#def fchebyshev_split(x, m):
#    """Compute the first `m` Chebyshev polynomials, but modified to allow a
#    split in the baseline at :math:`x=0`.  The intent is to allow a model fit
#    where a constant term is different for positive and negative `x`.
#
#    Parameters
#    ----------
#    x : array-like
#        Compute the Chebyshev polynomials at these abscissa values.
#    m : :class:`int`
#        The number of Chebyshev polynomials to compute.  For example, if
#        :math:`m = 3`, :math:`T_0 (x)`, :math:`T_1 (x)` and
#        :math:`T_2 (x)` will be computed.
#
#    Returns
#    -------
#    :class:`numpy.ndarray`
#    """
#    if isinstance(x, np.ndarray):
#        n = x.size
#    else:
#        n = 1
#    if m < 2:
#        raise ValueError('Order of polynomial must be at least 2.')
#    try:
#        dt = x.dtype
#    except AttributeError:
#        dt = np.float64
#    leg = np.ones((m, n), dtype=dt)
#    try:
#        leg[0, :] = (x >= 0).astype(x.dtype)
#    except AttributeError:
#        leg[0, :] = np.double(x >= 0)
#    if m > 2:
#        leg[2, :] = x
#    if m > 3:
#        for k in range(3, m):
#            leg[k, :] = 2.0 * x * leg[k-1, :] - leg[k-2, :]
#    return leg
#
#

#def fpoly(x, m):
#    """Compute the first `m` simple polynomials.
#
#    Parameters
#    ----------
#    x : array-like
#        Compute the simple polynomials at these abscissa values.
#    m : :class:`int`
#        The number of simple polynomials to compute.  For example, if
#        :math:`m = 3`, :math:`x^0`, :math:`x^1` and
#        :math:`x^2` will be computed.
#
#    Returns
#    -------
#    :class:`numpy.ndarray`
#    """
#    if isinstance(x, np.ndarray):
#        n = x.size
#    else:
#        n = 1
#    if m < 1:
#        raise ValueError('Order of polynomial must be at least 1.')
#    try:
#        dt = x.dtype
#    except AttributeError:
#        dt = np.float64
#    leg = np.ones((m, n), dtype=dt)
#    if m >= 2:
#        leg[1, :] = x
#    if m >= 3:
#        for k in range(2, m):
#            leg[k, :] = leg[k-1, :] * x
#    return leg


def func_fit(x, y, ncoeff, invvar=None, function_name='legendre', ia=None,
            inputans=None, inputfunc=None):
    """Fit `x`, `y` positions to a functional form.

    Parameters
    ----------
    x : array-like
        X values (independent variable).
    y : array-like
        Y values (dependent variable).
    ncoeff : :class:`int`
        Number of coefficients to fit.
    invvar : array-like, optional
        Weight values; inverse variance.
    function_name : :class:`str`, optional
        Function name, default 'legendre'.
    ia : array-like, optional
        An array of bool of length `ncoeff` specifying free (``True``) and
        fixed (``False``) parameters.
    inputans : array-like, optional
        An array of values of length `ncoeff` specifying the values of
        the fixed parameters.
    inputfunc : array-like, optional
        Multiply the function fit by these values.

    Returns
    -------
    :func:`tuple` of array-like
        Fit coefficients, length `ncoeff`; fitted values.

    Raises
    ------
    KeyError
        If an invalid function type is selected.
    ValueError
        If input dimensions do not agree.
    """
    if x.shape != y.shape:
        raise ValueError('Dimensions of X and Y do not agree!')
    if invvar is None:
        invvar = np.ones(x.shape, dtype=x.dtype)
    else:
        if invvar.shape != x.shape:
            raise ValueError('Dimensions of X and invvar do not agree!')
    if ia is None:
        ia = np.ones((ncoeff,), dtype=np.bool)
    if not ia.all():
        if inputans is None:
            inputans = np.zeros((ncoeff,), dtype=x.dtype)
    #
    # Select unmasked points
    #
    igood = (invvar > 0).nonzero()[0]
    ngood = len(igood)
    res = np.zeros((ncoeff,), dtype=x.dtype)
    yfit = np.zeros(x.shape, dtype=x.dtype)
    if ngood == 0:
        pass
    elif ngood == 1:
        res[0] = y[igood[0]]
        yfit += y[igood[0]]
    else:
        ncfit = min(ngood, ncoeff)
        function_map = {
            'legendre': basis.flegendre,
            'flegendre': basis.flegendre,
            'chebyshev': basis.fchebyshev,
            'fchebyshev': basis.fchebyshev,
            'chebyshev_split': basis.fchebyshev_split,
            'fchebyshev_split': basis.fchebyshev_split,
            'poly': basis.fpoly,
            'fpoly': basis.fpoly
            }
        try:
            legarr = function_map[function_name](x, ncfit).T
        except KeyError:
            raise KeyError('Unknown function type: {0}'.format(function_name))
        if inputfunc is not None:
            if inputfunc.shape != x.shape:
                raise ValueError('Dimensions of X and inputfunc do not agree!')
            legarr *= np.tile(inputfunc, ncfit).reshape(ncfit, x.shape[0])
        yfix = np.zeros(x.shape, dtype=x.dtype)
        nonfix = ia[0:ncfit].nonzero()[0]
        nparams = len(nonfix)
        fixed = (~ia[0:ncfit]).nonzero()[0]
        if len(fixed) > 0:
            yfix = np.dot(legarr.T, inputans * (1 - ia))
            ysub = y - yfix
            finalarr = legarr[nonfix, :]
        else:
            finalarr = legarr
            ysub = y
        # extra2 = finalarr * np.outer(np.ones((nparams,), dtype=x.dtype),
        #                             (invvar > 0))
        extra2 = finalarr * np.outer(np.ones((nparams,), dtype=x.dtype),
                                    invvar)
        alpha = np.dot(finalarr, extra2.T)
        # assert alpha.dtype == x.dtype
        if nparams > 1:
            # beta = np.dot(ysub * (invvar > 0), finalarr.T)
            beta = np.dot(ysub * invvar, finalarr.T)
            assert beta.dtype == x.dtype
            # uu,ww,vv = np.linalg.svd(alpha, full_matrices=False)
            res[nonfix] = np.linalg.solve(alpha, beta)
        else:
            # res[nonfix] = (ysub * (invvar > 0) * finalarr).sum()/alpha
            res[nonfix] = (ysub * invvar * finalarr).sum()/alpha
        if len(fixed) > 0:
            res[fixed] = inputans[fixed]
        yfit = np.dot(legarr.T, res[0:ncfit])
    return res, yfit


class TraceSet(object):
    """Implements the idea of a trace set.

    Attributes
    ----------
    func : :class:`str`
        Name of function type used to fit the trace set.
    xmin : float-like
        Minimum x value.
    xmax : float-like
        Maximum x value.
    coeff : array-like
        Coefficients of the trace set fit.
    nTrace : :class:`int`
        Number of traces in the object.
    ncoeff : :class:`int`
        Number of coefficients of the trace set fit.
    xjumplo : float-like
        Jump value, for BOSS readouts.
    xjumphi : float-like
        Jump value, for BOSS readouts.
    xjumpval : float-like
        Jump value, for BOSS readouts.
    lower : `int` or `float`, optional
        used for robust_polyfit_djs
        If set, reject points with data < model - lower * sigma.
    upper : `int` or `float`, optional
        used for robust_polyfit_djs
        If set, reject points with data > model + upper * sigma.
    outmask : array-like
        When initialized with x,y positions, this contains the rejected
        points.
    yfit : array-like
        When initialized with x,y positions, this contains the fitted y
        values.
    """
    #_func_map = {'poly': fpoly, 'legendre': flegendre,
    #                'chebyshev': fchebyshev}

    # ToDO Remove the kwargs and put in all the djs_reject parameters here
    def __init__(self, *args, **kwargs):
        """This class can be initialized either with a set of xy positions,
        or with a trace set HDU from a FITS file.
        """
        from astropy.io.fits.fitsrec import FITS_rec
        #from .math import djs_reject
        allowed_functions = ['polynomial', 'legendre', 'chebyshev']
        if len(args) == 1 and isinstance(args[0], FITS_rec):
            #
            # Initialize with FITS data
            #
            self.func = args[0]['FUNC'][0]
            self.xmin = args[0]['XMIN'][0]
            self.xmax = args[0]['XMAX'][0]
            self.coeff = args[0]['COEFF'][0]
            self.nTrace = self.coeff.shape[0]
            self.ncoeff = self.coeff.shape[1]
            if 'XJUMPLO' in args[0].dtype.names:
                self.xjumplo = args[0]['XJUMPLO'][0]
                self.xjumphi = args[0]['XJUMPHI'][0]
                self.xjumpval = args[0]['XJUMPVAL'][0]
            else:
                self.xjumplo = None
                self.xjumphi = None
                self.xjumpval = None
            self.lower = None
            self.upper = None
            self.outmask = None
            self.yfit = None
        elif len(args) == 2:
            #
            # Initialize with x, y positions.
            #
            xpos = args[0]
            ypos = args[1]



            self.nTrace = xpos.shape[0]
            if 'invvar' in kwargs:
                invvar = kwargs['invvar']
            else:
                invvar = None
                #invvar = np.ones(xpos.shape, dtype=xpos.dtype)
            if 'func' in kwargs:
                if kwargs['func'] not in allowed_functions:
                    msgs.error('Unrecognized function.')
                self.func = kwargs['func']
            else:
                self.func = 'legendre'
            if 'ncoeff' in kwargs:
                self.ncoeff = int(kwargs['ncoeff'])
            else:
                self.ncoeff = 3
            if 'xmin' in kwargs:
                self.xmin = np.float64(kwargs['xmin'])
            else:
                self.xmin = xpos.min()
            if 'xmax' in kwargs:
                self.xmax = np.float64(kwargs['xmax'])
            else:
                self.xmax = xpos.max()
            if 'maxiter' in kwargs:
                self.maxiter = int(kwargs['maxiter'])
            else:
                self.maxiter = 10
            if 'maxdev' in kwargs:
                self.maxdev = kwargs['maxdev']
            else:
                self.maxdev = None
            if 'inmask' in kwargs:
                inmask = kwargs['inmask']
            elif invvar is not None:
                inmask = (invvar > 0.0)
            else:
                inmask = np.ones(xpos.shape, dtype=np.bool)
            do_jump = False
            if 'xjumplo' in kwargs:
                do_jump = True
                self.xjumplo = np.float64(kwargs['xjumplo'])
            else:
                self.xjumplo = None
            if 'xjumphi' in kwargs:
                self.xjumphi = np.float64(kwargs['xjumphi'])
            else:
                self.xjumphi = None
            if 'xjumpval' in kwargs:
                self.xjumpval = np.float64(kwargs['xjumpval'])
            else:
                self.xjumpval = None
            if 'lower' in kwargs:
                self.lower = np.float64(kwargs['lower'])
            else:
                self.lower = None
            if 'upper' in kwargs:
                self.upper = np.float64(kwargs['upper'])
            else:
                self.upper = None

            self.coeff = np.zeros((self.nTrace, self.ncoeff+1), dtype=xpos.dtype)
            self.outmask = np.zeros(xpos.shape, dtype=np.bool)
            self.yfit = np.zeros(xpos.shape, dtype=xpos.dtype)
            for iTrace in range(self.nTrace):
                xvec = self.xnorm(xpos[iTrace, :], do_jump)
                if invvar is None:
                    thisinvvar = None
                else:
                    thisinvvar = invvar[iTrace, :]

                mask_djs, poly_coeff = utils.robust_polyfit_djs(xvec, ypos[iTrace, :], self.ncoeff,
                                                                function=self.func, maxiter = self.maxiter,
                                                                inmask = inmask[iTrace, :], invvar = thisinvvar,
                                                                lower = self.lower, upper = self.upper,
                                                                minx = self.xmin, maxx = self.xmax,
                                                                maxdev=self.maxdev,maxrej=None,groupdim=None,
                                                                groupsize=None,groupbadpix=None,grow=0,use_mad=False,sticky=False)
                ycurfit_djs = utils.func_val(poly_coeff, xvec, self.func, minx=self.xmin, maxx=self.xmax)

                ##Using robust_polyfit_djs to do the fitting and the following part are commented out by Feige
                #while (not qdone) and (iIter <= maxiter):
                #    res, ycurfit = func_fit(xvec, ypos[iTrace, :], self.ncoeff,
                #        invvar=tempivar*thismask, function_name=self.func)#function_name='poly')
                #    #ToDo: is this doing rejection? I think not???? THIS IS A MASSIVE BUG!!!! See IDL code.
                #    # ADd kwargs_reject in here like with iterfit
                #    thismask, qdone = djs_reject(ypos[iTrace, :], ycurfit,
                #                                invvar=tempivar)
                #    #thismask, qdone = djs_reject(ypos[iTrace, :], ycurfit,lower=3,upper=3)
                #    iIter += 1
                #self.yfit[iTrace, :] = ycurfit
                #self.coeff[iTrace, :] = res
                #self.outmask[iTrace, :] = thismask
                self.yfit[iTrace, :] = ycurfit_djs #ycurfit
                self.coeff[iTrace, :] = poly_coeff#[:-1] #res
                self.outmask[iTrace, :] = mask_djs #thismask

        else:
            msgs.error('Wrong number of arguments to TraceSet!')
            #raise PydlutilsException("Wrong number of arguments to TraceSet!")

    def xy(self, xpos=None, ignore_jump=False):
        """Convert from a trace set to an array of x,y positions.

        Parameters
        ----------
        xpos : array-like, optional
            If provided, evaluate the trace set at these positions.  Otherwise
            the positions will be constructed from the trace set object iself.
        ignore_jump : :class:`bool`, optional
            If ``True``, ignore any jump information in the `tset` object

        Returns
        -------
        :func:`tuple` of array-like
            The x, y positions.
        """
        #from .misc import djs_laxisgen
        do_jump = self.has_jump and (not ignore_jump)
        if xpos is None:
            xpos = djs_laxisgen([self.nTrace, self.nx], iaxis=1) + self.xmin
        ypos = np.zeros(xpos.shape, dtype=xpos.dtype)
        for iTrace in range(self.nTrace):
            xvec = self.xnorm(xpos[iTrace, :], do_jump)
            #legarr = self._func_map[self.func](xvec, self.ncoeff+1) #need to be norder+1 for utils functions
            ypos[iTrace, :] =  utils.func_val(self.coeff[iTrace, :], xvec, self.func, minx=self.xmin, maxx=self.xmax)
#            ypos[iTrace, :] = np.dot(legarr.T, self.coeff[iTrace, :])
        return (xpos, ypos)

    @property
    def has_jump(self):
        """``True`` if jump conditions are set.
        """
        return self.xjumplo is not None

    @property
    def xRange(self):
        """Range of x values.
        """
        return self.xmax - self.xmin

    @property
    def nx(self):
        """Number of x values.
        """
        return int(self.xRange + 1)

    # JFH This is no longer ndeeded. Below we changed to robust_polyfit_djs convention for writing this. Same value just one less parameter.
    @property
    def xmid(self):
        """Midpoint of x values.
        """
        return 0.5 * (self.xmin + self.xmax)

    def xnorm(self, xinput, jump):
        """Convert input x coordinates to normalized coordinates suitable
        for input to special polynomials.

        Parameters
        ----------
        xinput : array-like
            Input coordinates.
        jump : :class:`bool`
            Set to ``True`` if there is a jump.

        Returns
        -------
        array-like
            Normalized coordinates.
        """
        if jump:
            # Vector specifying what fraction of the jump has passed:
            jfrac = np.minimum(np.maximum(((xinput - self.xjumplo) /
                                (self.xjumphi - self.xjumplo)), 0.), 1.)
            # Conversion to "natural" x baseline:
            xnatural = xinput + jfrac * self.xjumpval
        else:
            xnatural = xinput
        return 2.0 * (xnatural - self.xmin)/self.xRange - 1.0


def traceset2xy(tset, xpos=None, ignore_jump=False):
    """Convert from a trace set to an array of x,y positions.

    Parameters
    ----------
    tset : :class:`TraceSet`
        A :class:`TraceSet` object.
    xpos : array-like, optional
        If provided, evaluate the trace set at these positions.  Otherwise
        the positions will be constructed from the trace set object iself.
    ignore_jump : bool, optional
        If ``True``, ignore any jump information in the `tset` object

    Returns
    -------
    :func:`tuple` of array-like
        The x, y positions.
    """
    return tset.xy(xpos, ignore_jump)


def xy2traceset(xpos, ypos, **kwargs):
    """Convert from x,y positions to a trace set.

    Parameters
    ----------
    xpos, ypos : array-like
        X,Y positions corresponding as [nx,Ntrace] arrays.
    invvar : array-like, optional
        Inverse variances for fitting.
    func : :class:`str`, optional
        Function type for fitting; defaults to 'legendre'.
    ncoeff : :class:`int`, optional
        Number of coefficients to fit.  Defaults to 3.
    xmin, xmax : :class:`float`, optional
        Explicitly set minimum and maximum values, instead of computing
        them from `xpos`.
    maxiter : :class:`int`, optional
        Maximum number of rejection iterations; set to 0 for no rejection;
        default to 10.
    inmask : array-like, optional
        Mask set to 1 for good points and 0 for rejected points;
        same dimensions as `xpos`, `ypos`.  Points rejected by `inmask`
        are always rejected from the fits (the rejection is "sticky"),
        and will also be marked as rejected in the outmask attribute.
    ia, inputans, inputfunc : array-like, optional
        These arguments will be passed to :func:`func_fit`.
    xjumplo : :class:`float`, optional
        x position locating start of an x discontinuity
    xjumphi : :class:`float`, optional
        x position locating end of that x discontinuity
    xjumpval : :class:`float`, optional
        magnitude of the discontinuity "jump" between those bounds
        (previous 3 keywords motivated by BOSS 2-phase readout)

    Returns
    -------
    :class:`TraceSet`
        A :class:`TraceSet` object.
    """
    return TraceSet(xpos, ypos, **kwargs)




def djs_reject(data, model, outmask=None, inmask=None,
               invvar=None, lower=None, upper=None, percentile=False, maxdev=None,
               maxrej=None, groupdim=None, groupsize=None, groupbadpix=False,
               grow=0, sticky=False, use_mad=False):
    """Routine to reject points when doing an iterative fit to data.

    Parameters
    ----------
    data : :class:`numpy.ndarray`
        The data
    model : :class:`numpy.ndarray`
        The model, must have the same number of dimensions as `data`.
    outmask : :class:`numpy.ndarray`, optional
        Output mask, generated by a previous call to `djs_reject`.  If sticky=True, then bad points accumulate in this mask
        between calls. Otherwise, this mask is only  used to determine if the rejection iterations are complete (e.g. to set qdone).
        Although this parameter is technically optional, it will almost always be set. If not supplied,
        this mask will be initialized to a mask that masks nothing and qdone will always be returned as True.
    inmask : :class:`numpy.ndarray`, optional
        Input mask.  Bad points are marked with a value that evaluates to ``False``.
        Must have the same number of dimensions as `data`. Points masked as bad "False" in the inmask
        will also always evaluate to "False" in the outmask
    invvar : :class:`numpy.ndarray`, optional
        Inverse variance of the data, used to reject points based on the values
        of `upper` and `lower`.
    lower : :class:`int` or :class:`float`, optional
        If set, reject points with data < model - lower * sigm, where sigma = 1.0/sqrt(invvar)
    upper : :class:`int` or :class:`float`, optional
        If set, reject points with data > model + upper * sigma, where sigma = 1.0/sqrt(invvar)
    maxdev : :class:`int` or :class:`float`, optional
        If set, reject points with abs(data-model) > maxdev.  It is permitted to
        set all three of `lower`, `upper` and `maxdev`.
    maxrej: :class:`int` or :class:`numpy.ndarray`, optional
        Maximum number of points to reject in this iteration.  If `groupsize` or
        `groupdim` are set to arrays, this should be an array as well.
    groupdim: class: `int`
        Dimension along which to group the data; set to 1 to group along the 1st dimension, 2 for the 2nd dimension, etc.
        If data has shape [100,200], then setting GROUPDIM=2 is equivalent to grouping the data with groupsize=100.
        In either case, there are 200 groups, specified by ``[*,i]``. NOT WELL TESTED IN PYTHON!
    groupsize: class: `int`
        If this and maxrej are set, then reject a maximum of maxrej points per group of groupsize points, where the grouping is performed in the
        along the dimension of the data vector. (For use in curve fitting, one probably wants to make sure that data is sorted according to the indpeendent
        variable. For multi-dimensional arrays where one desires this grouping along each dimension, then groupdim should be set.
        If groupdim is also set, then this specifies sub-groups within that.
    groupbadpix : :class:`bool`, optional
        If set to ``True``, consecutive sets of bad pixels are considered groups,
        overriding the values of `groupsize`.
    grow : :class:`int`, optional, default = 0
        If set to a non-zero integer, N, the N nearest neighbors of rejected
        pixels will also be rejected.
    sticky : :class:`bool`, optional
        If set to True then points rejected in outmask from a previous call to djs_reject are kept rejected. If
        set to False, if a fit (model) changes between iterations, points can alternate from being rejected to not rejected.
    use_mad : :class: `bool`, optional, defaul = False
        It set to ``True``, compute the median of the maximum absolute deviation between the data and use this for the rejection instead of
        the default which is to compute the standard deviation of the yarray - modelfit. Note that it is not possible to specify use_mad=True
        and also pass in values invvar, and the code will return an error if this is done.

    Returns
    -------
    outmask (np.ndarray, boolean):
        mask where rejected data values are ``False``
    qdone (boolean):
        a value set to "True" if  `djs_reject` believes there is no
        further rejection to be done. This will be set to "False" if the
        points marked as rejected in the outmask have changed. It will
        be set to "True" when the same points are rejected in outmask as
        from a previous call.  It will also be set to "False" if model
        is set to None. Recall that outmask is also an optional input
        parameter. If it is not set, then qdone will simply return true,
        so outmask needs to be input from the previous iteration for the
        routine to do something meaningful.

    Raises
    ------
    ValueError
        If dimensions of various inputs do not match.
    """
    #from .misc import djs_laxisnum
    #

    # ToDO It would be nice to come up with a way to use MAD but also use the errors in the rejection, i.e. compute the rejection threhsold using the mad.

    if upper is None and lower is None and maxdev is None:
        msgs.warn('upper, lower, and maxdev are all set to None. No rejection performed since no rejection criteria were specified.')

    if (use_mad and (invvar is not None)):
        raise ValueError('use_mad can only be set to True innvar = None. This code only computes a mad'
                         ' if errors are not input')

    # Create outmask setting = True for good data.
    #
    # ToDo JFH: I think it would actually make more sense for outmask be a required input parameter (named lastmask or something like that).
    if outmask is None:
        outmask = np.ones(data.shape, dtype='bool')
        msgs.warn('outmask was not specified as an input parameter. Cannot asess convergence of rejection -- qdone is automatically True')
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
            if isinstance(maxrej, (int,float)) | isinstance(groupsize, (int,float)):
                groupsize1=np.asarray([groupsize])
            else:
                if len(maxrej) != len(groupsize):
                    raise ValueError('maxrej and groupsize must have the same number of elements.')
                groupsize1=groupsize
        else:
            groupsize1 = np.asarray([len(data)])
        if isinstance(maxrej,(int,float)):
            maxrej1 = np.asarray([maxrej])
        else:
            maxrej1 = maxrej
    if invvar is None:
        if inmask is not None:
            igood = (inmask & outmask)
        else:
            igood = outmask
        if (np.sum(igood) > 1):
            if use_mad is True:
                sigma = 1.4826*np.median(np.abs(data[igood] - model[igood]))
            else:
                sigma = np.std(data[igood] - model[igood])
            invvar = utils.inverse(sigma**2)
        else:
            invvar = 0.0


    diff = data - model
    chi = diff * np.sqrt(invvar)

    #
    # The working array is badness, which is set to zero for good points
    # (or points already rejected), and positive values for bad points.
    # The values determine just how bad a point is, either corresponding
    # to the number of sigma above or below the fit, or to the number
    # of multiples of maxdev away from the fit.
    #
    badness = np.zeros(outmask.shape, dtype=data.dtype)

    if percentile:
        if inmask is not None:
            igood = (inmask & outmask)
        else:
            igood = outmask
        if (np.sum(igood)> 1):
            if lower is not None:
                lower_chi = np.percentile(chi[igood],lower)
            else:
                lower_chi = -np.inf
            if upper is not None:
                upper_chi = np.percentile(chi[igood], upper)
            else:
                upper_chi = np.inf
    #
    # Decide how bad a point is according to lower.
    #
    if lower is not None:
        if percentile:
            qbad = chi < lower_chi
        else:
            qbad = chi < -lower
        badness += np.fmax(-chi,0.0)*qbad
    #
    # Decide how bad a point is according to upper.
    #
    if upper is not None:
        if percentile:
            qbad = chi > upper_chi
        else:
            qbad = chi > upper
        badness += np.fmax(chi,0.0)*qbad
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
            for ivec in range(np.fmax(dimnum.max(),1)):
                #
                # At this point it is not possible that dimnum is not set.
                #
                if len(groupdim) == 0:
                    indx = np.arange(data.size)
                else:
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
                    # The IDL version of this test makes no sense because
                    # groupsize will always be set.
                    #
                    if False:
                        ngroups = 1
                        groups_lower = [0, ]
                        groups_upper = [nin - 1, ]
                    else:
                        ngroups = nin//groupsize1[iloop] + 1
                        groups_lower = np.arange(ngroups, dtype='i4')*groupsize1[iloop]
                        foo = (np.arange(ngroups, dtype='i4')+1)*groupsize1[iloop]
                        groups_upper = np.where(foo < nin, foo, nin) - 1

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
                        if np.sum(badness[jj] != 0) > maxrej1[iloop]:
                            isort = badness[jj].argsort()
                            #
                            # Make the following points good again.
                            #
                            badness[jj[isort[0:nii-maxrej1[iloop]]]] = 0
                        i1 += groupsize1[iloop]
    #
    # Now modify outmask, rejecting points specified by inmask=0, outmask=0
    # if sticky is set, or badness > 0.
    #
    # print(badness)
    newmask = badness == 0.0

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
    qdone = bool(np.all(newmask == outmask))
    # JFH This needs to be a python (rather than a numpy) boolean to avoid painful problems when comparing
    # to python True and False booleans

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



def djs_laxisgen(dims, iaxis=0):
    """Returns an integer array where each element of the array is set
    equal to its index number along the specified axis.

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
        If `iaxis` is greater than or equal to the number of dimensions.

    Notes
    -----
    For two or more dimensions, there is no difference between this routine
    and :func:`~pydl.pydlutils.misc.djs_laxisnum`.

    Examples
    --------
    >>> from pydl.pydlutils.misc import djs_laxisgen
    >>> print(djs_laxisgen([4,4]))
    [[0 0 0 0]
     [1 1 1 1]
     [2 2 2 2]
     [3 3 3 3]]
    """
    ndimen = len(dims)
    if ndimen == 1:
        return np.arange(dims[0], dtype='i4')
    return djs_laxisnum(dims, iaxis)



### Following part are imported from pydl spheregroup
class chunks(object):
    """chunks class

    Functions for creating and manipulating spherical chunks are implemented
    as methods on this class.
    """

    def __init__(self, ra, dec, minSize):
        """Init creates an object whose attributes are similar those created
        by the setchunks() function in the spheregroup library.
        """
        #
        # Save the value of minSize
        #
        self.minSize = minSize
        #
        # Find maximum and minimum dec (in degrees)
        #
        decMin = dec.min()
        decMax = dec.max()
        decRange = decMax - decMin
        #
        # Find the declination boundaries; make them an integer multiple of
        # minSize, with extra room (one cell) on the edges.
        #
        self.nDec = 3 + int(np.floor(decRange/minSize))
        decRange = minSize*float(self.nDec)
        decMin = decMin - 0.5*(decRange - decMax + decMin)
        decMax = decMin + decRange
        if decMin < -90.0 + 3.0*minSize:
            decMin = -90.0
        if decMax > 90.0 - 3.0*minSize:
            decMax = 90.0
        self.decBounds = decMin + ((decMax - decMin) * np.arange(self.nDec + 1,
                                    dtype='d'))/float(self.nDec)
        #
        # Find ra offset which minimizes the range in ra (this should take care
        # of the case that ra crosses zero in some parts
        #
        if abs(self.decBounds[self.nDec]) > abs(self.decBounds[0]):
            cosDecMin = np.cos(np.deg2rad(self.decBounds[self.nDec]))
        else:
            cosDecMin = np.cos(np.deg2rad(self.decBounds[0]))
        if cosDecMin <= 0.0:
            msgs.error("cosDecMin={0:f} not positive in setchunks().".format(cosDecMin))
        self.raRange, self.raOffset = self.rarange(ra, minSize/cosDecMin)
        self.raMin, self.raMax = self.getraminmax(ra, self.raOffset)
        #
        # Isn't this redundant?
        #
        self.raRange = self.raMax - self.raMin
        #
        # For each declination slice, find the number of ra divisions
        # necessary and set them
        #
        self.raBounds = list()
        self.nRa = list()
        for i in range(self.nDec):
            #
            # Get maximum declination and its cosine
            #
            if abs(self.decBounds[i]) > abs(self.decBounds[i+1]):
                cosDecMin = np.cos(np.deg2rad(self.decBounds[i]))
            else:
                cosDecMin = np.cos(np.deg2rad(self.decBounds[i+1]))
            if cosDecMin <= 0.0:
                msgs.error("cosDecMin={0:f} not positive in setchunks().".format(cosDecMin))
            #
            # Get raBounds array for this declination array, leave an extra
            # cell on each end
            #
            self.nRa.append(3 + int(np.floor(cosDecMin*self.raRange/minSize)))
            raRangeTmp = minSize*float(self.nRa[i])/cosDecMin
            raMinTmp = self.raMin - 0.5*(raRangeTmp-self.raMax+self.raMin)
            raMaxTmp = raMinTmp + raRangeTmp
            #
            # If we cannot avoid the 0/360 point, embrace it
            #
            if (raRangeTmp >= 360.0 or
                    raMinTmp <= minSize/cosDecMin or
                    raMaxTmp >= 360.0 - minSize/cosDecMin or
                    abs(self.decBounds[i]) == 90.0):
                raMinTmp = 0.0
                raMaxTmp = 360.0
                raRangeTmp = 360.0
            if self.decBounds[i] == -90.0 or self.decBounds[i+1] == 90.0:
                self.nRa[i] = 1
            self.raBounds.append(raMinTmp +
                (raMaxTmp - raMinTmp) * np.arange(self.nRa[i] + 1, dtype='d') /
                float(self.nRa[i]))
        #
        # Create an empty set of lists to hold the output of self.assign()
        #
        self.chunkList = [[list() for j in range(self.nRa[i])] for i in range(self.nDec)]
        #
        # nChunkMax will be the length of the largest list in chunkList
        # it is computed by chunks.assign()
        #
        self.nChunkMax = 0
        return

    def rarange(self, ra, minSize):
        """Finds the offset which yields the smallest raRange & returns both.

        Notes
        -----

        .. warning:: This is not (yet) well-defined for the case of only one point.
        """
        NRA = 6
        raRangeMin = 361.0
        raOffset = 0.0
        EPS = 1.0e-5
        for j in range(NRA):
            raMin, raMax = self.getraminmax(ra, 360.0*float(j)/float(NRA))
            raRange = raMax-raMin
            if (2.0*(raRange-raRangeMin)/(raRange+raRangeMin) < -EPS and
                    raMin > minSize and raMax < 360.0 - minSize):
                raRangeMin = raRange
                raOffset = 360.0*float(j)/float(NRA)
        return (raRangeMin, raOffset)

    def getraminmax(self, ra, raOffset):
        """Utility function used by rarange.
        """
        currRa = np.fmod(ra + raOffset, 360.0)
        return (currRa.min(), currRa.max())

    def cosDecMin(self, i):
        """Frequently used utility function.
        """
        if abs(self.decBounds[i]) > abs(self.decBounds[i+1]):
            return np.cos(np.deg2rad(self.decBounds[i]))
        else:
            return np.cos(np.deg2rad(self.decBounds[i+1]))

    def assign(self, ra, dec, marginSize):
        """Take the objects and the chunks (already defined in the constructor)
        and assign the objects to the appropriate chunks, with some leeway
        given by the parameter marginSize.  Basically, at the end, each
        chunk should be associated with a list of the objects that belong
        to it.
        """
        if marginSize >= self.minSize:
            msgs.error("marginSize>=minSize ({0:f}={1:f}) in chunks.assign().".format(marginSize, self.minSize))
        chunkDone = [[False for j in range(self.nRa[i])] for i in range(self.nDec)]
        for i in range(ra.size):
            currRa = np.fmod(ra[i] + self.raOffset, 360.0)
            try:
                raChunkMin, raChunkMax, decChunkMin, decChunkMax = self.getbounds(currRa, dec[i], marginSize)
            except:
                continue
            #
            # Reset chunkDone.  This is silly, but is necessary to
            # reproduce the logic.
            #
            for decChunk in range(decChunkMin, decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin]-1, raChunkMax[decChunk-decChunkMin]+2):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        chunkDone[decChunk][currRaChunk] = False
            for decChunk in range(decChunkMin, decChunkMax+1):
                for raChunk in range(raChunkMin[decChunk-decChunkMin], raChunkMax[decChunk-decChunkMin]+1):
                    if raChunk < 0:
                        currRaChunk = (raChunk+self.nRa[decChunk]) % self.nRa[decChunk]
                    elif raChunk > self.nRa[decChunk]-1:
                        currRaChunk = (raChunk-self.nRa[decChunk]) % self.nRa[decChunk]
                    else:
                        currRaChunk = raChunk
                    if currRaChunk >= 0 and currRaChunk <= self.nRa[decChunk]-1:
                        if not chunkDone[decChunk][currRaChunk]:
                            self.chunkList[decChunk][currRaChunk].append(i)
                            #
                            # Update nChunkMax
                            #
                            if len(self.chunkList[decChunk][currRaChunk]) > self.nChunkMax:
                                self.nChunkMax = len(self.chunkList[decChunk][currRaChunk])
                            chunkDone[decChunk][currRaChunk] = True
        return

    def getbounds(self, ra, dec, marginSize):
        """Find the set of chunks a point (with margin) belongs to.
        """
        #
        # Find the declination slice without regard to marginSize
        #
        decChunkMin = int(np.floor((dec - self.decBounds[0]) *
            float(self.nDec) /
            (self.decBounds[self.nDec]-self.decBounds[0])))
        decChunkMax = decChunkMin
        if decChunkMin < 0 or decChunkMin > self.nDec - 1:
            msgs.error("decChunkMin out of range in chunks.getbounds().")
        #
        # Set minimum and maximum bounds of dec
        #
        while dec - self.decBounds[decChunkMin] < marginSize and decChunkMin > 0:
            decChunkMin -= 1
        while self.decBounds[decChunkMax+1] - dec < marginSize and decChunkMax < self.nDec - 1:
            decChunkMax += 1
        #
        # Find ra chunk bounds for each dec chunk
        #
        raChunkMin = np.zeros(decChunkMax-decChunkMin+1, dtype='i4')
        raChunkMax = np.zeros(decChunkMax-decChunkMin+1, dtype='i4')
        for i in range(decChunkMin, decChunkMax+1):
            cosDecMin = self.cosDecMin(i)
            raChunkMin[i-decChunkMin] = int(np.floor((ra - self.raBounds[i][0]) *
                float(self.nRa[i]) /
                (self.raBounds[i][self.nRa[i]] - self.raBounds[i][0])))
            raChunkMax[i-decChunkMin] = raChunkMin[i-decChunkMin]
            if raChunkMin[i-decChunkMin] < 0 or raChunkMin[i-decChunkMin] > self.nRa[i]-1:
                msgs.error("raChunkMin out of range in chunks.getbounds().")
            #
            # Set minimum and maximum bounds of ra
            #
            raCheck = raChunkMin[i-decChunkMin]
            keepGoing = True
            while keepGoing and raCheck > -1:
                if raCheck >= 0 and raCheck < self.nRa[i]:
                    keepGoing = (ra - self.raBounds[i][raCheck])*cosDecMin < marginSize
                else:
                    keepGoing = False
                if keepGoing:
                    raCheck -= 1
            raChunkMin[i-decChunkMin] = raCheck
            raCheck = raChunkMax[i-decChunkMin]
            keepGoing = True
            while keepGoing and raCheck < self.nRa[i]:
                if raCheck >= 0 and raCheck < self.nRa[i]:
                    keepGoing = (self.raBounds[i][raCheck+1]-ra)*cosDecMin < marginSize
                else:
                    keepGoing = False
                if keepGoing:
                    raCheck += 1
            raChunkMax[i-decChunkMin] = raCheck
        return (raChunkMin, raChunkMax, decChunkMin, decChunkMax)

    def get(self, ra, dec):
        """Find the chunk to which a given point belongs.
        """
        #
        # Find dec chunk
        #
        decChunk = int(np.floor((dec - self.decBounds[0]) *
            float(self.nDec) /
            (self.decBounds[self.nDec]-self.decBounds[0])))
        #
        # Find ra chunk
        #
        if decChunk < self.nDec and decChunk >= 0:
            raChunk = int(np.floor((ra - self.raBounds[decChunk][0]) *
                float(self.nRa[decChunk]) /
                (self.raBounds[decChunk][self.nRa[decChunk]] - self.raBounds[decChunk][0])))
            if raChunk < 0 or raChunk > self.nRa[decChunk]-1:
                msgs.error("raChunk out of range in chunks.get()")
        else:
            raChunk = -1
        return (raChunk, decChunk)

    def friendsoffriends(self, ra, dec, linkSep):
        """Friends-of-friends using chunked data.
        """
        nPoints = ra.size
        inGroup = np.zeros(nPoints, dtype='i4') - 1
        #
        # mapGroups contains an equivalency mapping of groups.  mapGroup[i]=j
        # means i and j are actually the same group.  j<=i always, by design.
        # The largest number of groups you can get
        # (assuming linkSep < marginSize < minSize) is 9 times the number of
        # targets
        #
        mapGroups = np.zeros(9*nPoints, dtype='i4') - 1
        nMapGroups = 0
        for i in range(self.nDec):
            for j in range(self.nRa[i]):
                if len(self.chunkList[i][j]) > 0:
                    chunkGroup = self.chunkfriendsoffriends(ra, dec, self.chunkList[i][j], linkSep)
                    for k in range(chunkGroup.nGroups):
                        minEarly = 9*nPoints
                        l = chunkGroup.firstGroup[k]
                        while l != -1:
                            if inGroup[self.chunkList[i][j][l]] != -1:
                                checkEarly = inGroup[self.chunkList[i][j][l]]
                                while mapGroups[checkEarly] != checkEarly:
                                    checkEarly = mapGroups[checkEarly]
                                minEarly = min(minEarly, checkEarly)
                            else:
                                inGroup[self.chunkList[i][j][l]] = nMapGroups
                            l = chunkGroup.nextGroup[l]
                        if minEarly == 9*nPoints:
                            mapGroups[nMapGroups] = nMapGroups
                        else:
                            mapGroups[nMapGroups] = minEarly
                            l = chunkGroup.firstGroup[k]
                            while l != -1:
                                checkEarly = inGroup[self.chunkList[i][j][l]]
                                while mapGroups[checkEarly] != checkEarly:
                                    tmpEarly = mapGroups[checkEarly]
                                    mapGroups[checkEarly] = minEarly
                                    checkEarly = tmpEarly
                                mapGroups[checkEarly] = minEarly
                                l = chunkGroup.nextGroup[l]
                        nMapGroups += 1
        #
        # Now all groups which are mapped to themselves are the real groups
        # Make sure the mappings are set up to go all the way down.
        #
        nGroups = 0
        for i in range(nMapGroups):
            if mapGroups[i] != -1:
                if mapGroups[i] == i:
                    mapGroups[i] = nGroups
                    nGroups += 1
                else:
                    mapGroups[i] = mapGroups[mapGroups[i]]
            else:
                msgs.error("MapGroups[{0:d}]={1:d} in chunks.friendsoffriends().".format(i, mapGroups[i]))
        for i in range(nPoints):
            inGroup[i] = mapGroups[inGroup[i]]
        firstGroup = np.zeros(nPoints, dtype='i4') - 1
        nextGroup = np.zeros(nPoints, dtype='i4') - 1
        multGroup = np.zeros(nPoints, dtype='i4')
        for i in range(nPoints-1, -1, -1):
            nextGroup[i] = firstGroup[inGroup[i]]
            firstGroup[inGroup[i]] = i
        for i in range(nGroups):
            j = firstGroup[i]
            while j != -1:
                multGroup[i] += 1
                j = nextGroup[j]
        return (inGroup, multGroup, firstGroup, nextGroup, nGroups)

    def chunkfriendsoffriends(self, ra, dec, chunkList, linkSep):
        """Does friends-of-friends on the ra, dec that are defined by
        chunkList.
        """
        #
        # Convert ra, dec into something that can be digested by the
        # groups object.
        #
        x = np.deg2rad(np.vstack((ra[chunkList], dec[chunkList])))
        radLinkSep = np.deg2rad(linkSep)
        group = groups(x, radLinkSep, 'sphereradec')
        return group


class groups(object):
    """Group a set of objects (a list of coordinates in some space) based on
    a friends-of-friends algorithm
    """

    @staticmethod
    def euclid(x1, x2):
        """Pythagorean theorem in Euclidean space with arbitrary number
        of dimensions.
        """
        return np.sqrt(((x1-x2)**2).sum())

    @staticmethod
    def sphereradec(x1, x2):
        """Separation of two points on a 2D-sphere, assuming they are in
        longitude-latitude or right ascension-declination form.  Assumes
        everything is already in radians.
        """
        return gcirc(x1[0], x1[1], x2[0], x2[1], units=0)

    def __init__(self, coordinates, distance, separation='euclid'):
        """Init creates an object and performs the friends-of-friends
        algorithm.  The coordinates can have arbitrary dimensions, with each
        column representing one of the dimensions.  Each row defines an object.
        If separation is not defined it defaults to Euclidean space.
        """
        #
        # Find a separation function
        #
        if callable(separation):
            self.separation = separation
        elif isinstance(separation, str):
            if separation == 'euclid':
                self.separation = self.euclid
            elif separation == 'sphereradec':
                self.separation = self.sphereradec
            else:
                msgs.error("Unknown separation function: {0}.".format(separation))
        else:
            msgs.error("Improper type for separation!")
        #
        # Save information about the coordinates.
        #
        nGroups = 0
        nTargets = coordinates.shape[1]
        multGroup = np.zeros(nTargets, dtype='i4')
        firstGroup = np.zeros(nTargets, dtype='i4') - 1
        nextGroup = np.zeros(nTargets, dtype='i4') - 1
        inGroup = np.arange(nTargets, dtype='i4')
        #
        # Find all the other targets associated with each target
        #
        for i in range(nTargets):
            nTmp = 0
            minGroup = nGroups
            for j in range(nTargets):
                sep = self.separation(coordinates[:, i], coordinates[:, j])
                if sep <= distance:
                    multGroup[nTmp] = j
                    minGroup = min(minGroup, inGroup[j])
                    nTmp += 1
            #
            # Use this minimum for all
            #
            for j in range(nTmp):
                if inGroup[multGroup[j]] < nTargets:
                    k = firstGroup[inGroup[multGroup[j]]]
                    while k != -1:
                        inGroup[k] = minGroup
                        k = nextGroup[k]
                inGroup[multGroup[j]] = minGroup
            #
            # If it is a new group (no earlier groups), increment nGroups
            #
            if minGroup == nGroups:
                nGroups += 1
            for j in range(i+1):
                firstGroup[j] = -1
            for j in range(i, -1, -1):
                nextGroup[j] = firstGroup[inGroup[j]]
                firstGroup[inGroup[j]] = j
        #
        # Renumber to get rid of the numbers which were skipped
        #
        renumbered = np.zeros(nTargets, dtype='bool')
        nTmp = nGroups
        nGroups = 0
        for i in range(nTargets):
            if not renumbered[i]:
                j = firstGroup[inGroup[i]]
                while j != -1:
                    inGroup[j] = nGroups
                    renumbered[j] = True
                    j = nextGroup[j]
                nGroups += 1
        #
        # Reset the values of firstGroup and inGroup
        #
        firstGroup[:] = -1
        for i in range(nTargets-1, -1, -1):
            nextGroup[i] = firstGroup[inGroup[i]]
            firstGroup[inGroup[i]] = i
        #
        # Get the multiplicity
        #
        for i in range(nGroups):
            multGroup[i] = 0
            j = firstGroup[i]
            while j != -1:
                multGroup[i] += 1
                j = nextGroup[j]
        #
        # Set attributes
        #
        self.nGroups = nGroups
        self.nTargets = nTargets
        self.inGroup = inGroup
        self.multGroup = multGroup
        self.firstGroup = firstGroup
        self.nextGroup = nextGroup
        return


def spheregroup(ra, dec, linklength, chunksize=None):
    """Perform friends-of-friends grouping given ra/dec coordinates.

    Parameters
    ----------
    ra, dec : :class:`numpy.ndarray`
        Arrays of coordinates to group in decimal degrees.
    linklength : :class:`float`
        Linking length for the groups in decimal degrees.
    chunksize : :class:`float`, optional
        Break up the sphere into chunks of this size in decimal degrees.

    Returns
    -------
    :func:`tuple`
        A tuple containing the group number of each object, the multiplicity
        of each group, the first member of each group, and the next
        member of the group for each object.

    Raises
    ------
    msgs.error
        If the array of coordinates only contains one point.

    Notes
    -----
    It is important that `chunksize` >= 4 * `linklength`.  This is enforced.

    .. warning:: Behavior at the poles is not well tested.
    """
    npoints = ra.size
    if npoints == 1:
        msgs.error("Cannot group only one point!")
    #
    # Define the chunksize
    #
    if chunksize is not None:
        if chunksize < 4.0*linklength:
            chunksize = 4.0*linklength
            msgs.warn("chunksize changed to {0:.2f}.".format(chunksize))
    else:
        chunksize = max(4.0*linklength, 0.1)
    #
    # Initialize chunks
    #
    chunk = chunks(ra, dec, chunksize)
    chunk.assign(ra, dec, linklength)
    #
    # Run friends-of-friends
    #
    ingroup, multgroup, firstgroup, nextgroup, ngroups = chunk.friendsoffriends(ra, dec, linklength)
    #
    # Renumber the groups in order of appearance
    #
    renumbered = np.zeros(npoints, dtype='bool')
    iclump = 0
    for i in range(npoints):
        if not renumbered[i]:
            j = firstgroup[ingroup[i]]
            while j != -1:
                ingroup[j] = iclump
                renumbered[j] = True
                j = nextgroup[j]
            iclump += 1
    #
    # Reset the index lists
    #
    firstgroup[:] = -1
    for i in range(npoints-1, -1, -1):
        nextgroup[i] = firstgroup[ingroup[i]]
        firstgroup[ingroup[i]] = i
    #
    # Reset the multiplicities
    #
    multgroup[:] = 0
    for i in range(ngroups):
        j = firstgroup[i]
        while j != -1:
            multgroup[i] += 1
            j = nextgroup[j]
    return (ingroup, multgroup, firstgroup, nextgroup)


def spherematch(ra1, dec1, ra2, dec2, matchlength, chunksize=None,
                maxmatch=1):
    """Match points on a sphere.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : :class:`numpy.ndarray`
        The sets of coordinates to match.  Assumed to be in decimal degrees
    matchlength : :class:`float`
        Two points closer than this separation are matched. Assumed to be in decimal degrees.
    chunksize : :class:`float`, optional
        Value to pass to chunk assignment.
    maxmatch : :class:`int`, optional
        Allow up to `maxmatch` matches per coordinate.  Default 1. If set to zero,
        All possible matches will be returned.

    Returns
    -------
    :func:`tuple`
        A tuple containing the indices into the first set of points, the
        indices into the second set of points and the match distance in
        decimal degrees.

    Notes
    -----
    If you have sets of coordinates that differ in size, call this function
    with the larger list first.  This exploits the inherent asymmetry in the
    underlying code to reduce memory use.

    .. warning:: Behavior at the poles is not well tested.
    """
    #
    # Set default values
    #
    if chunksize is None:
        chunksize = max(4.0*matchlength, 0.1)
    #
    # Check input size
    #
    if ra1.size == 1:
        msgs.error("Change the order of the sets of coordinates!")
    #
    # Initialize chunks
    #
    chunk = chunks(ra1, dec1, chunksize)
    chunk.assign(ra2, dec2, matchlength)
    #
    # Create return arrays
    #
    match1 = list()
    match2 = list()
    distance12 = list()
    for i in range(ra1.size):
        currra = np.fmod(ra1[i]+chunk.raOffset, 360.0)
        rachunk, decchunk = chunk.get(currra, dec1[i])
        jmax = len(chunk.chunkList[decchunk][rachunk])
        if jmax > 0:
            for j in range(jmax):
                k = chunk.chunkList[decchunk][rachunk][j]
                sep = gcirc(ra1[i], dec1[i], ra2[k], dec2[k], units=2)/3600.0
                if sep < matchlength:
                    match1.append(i)
                    match2.append(k)
                    distance12.append(sep)
    #
    # Sort distances
    #
    omatch1 = np.array(match1)
    omatch2 = np.array(match2)
    odistance12 = np.array(distance12)
    s = odistance12.argsort()
    #
    # Retain only desired matches
    #
    if maxmatch > 0:
        gotten1 = np.zeros(ra1.size, dtype='i4')
        gotten2 = np.zeros(ra2.size, dtype='i4')
        nmatch = 0
        for i in range(omatch1.size):
            if (gotten1[omatch1[s[i]]] < maxmatch and
                    gotten2[omatch2[s[i]]] < maxmatch):
                gotten1[omatch1[s[i]]] += 1
                gotten2[omatch2[s[i]]] += 1
                nmatch += 1
        match1 = np.zeros(nmatch, dtype='i4')
        match2 = np.zeros(nmatch, dtype='i4')
        distance12 = np.zeros(nmatch, dtype='d')
        gotten1[:] = 0
        gotten2[:] = 0
        nmatch = 0
        for i in range(omatch1.size):
            if (gotten1[omatch1[s[i]]] < maxmatch and
                    gotten2[omatch2[s[i]]] < maxmatch):
                gotten1[omatch1[s[i]]] += 1
                gotten2[omatch2[s[i]]] += 1
                match1[nmatch] = omatch1[s[i]]
                match2[nmatch] = omatch2[s[i]]
                distance12[nmatch] = odistance12[s[i]]
                nmatch += 1
    else:
        match1 = omatch1[s]
        match2 = omatch2[s]
        distance12 = odistance12[s]
    return (match1, match2, distance12)


def gcirc(ra1, dec1, ra2, dec2, units=2):
    """Computes rigorous great circle arc distances.

    Parameters
    ----------
    ra1, dec1, ra2, dec2 : :class:`float` or array-like
        RA and Dec of two points.
    units : { 0, 1, 2 }, optional
        * units = 0: everything is already in radians
        * units = 1: RA in hours, dec in degrees, distance in arcsec.
        * units = 2: RA, dec in degrees, distance in arcsec (default)

    Returns
    -------
    :class:`float` or array-like
        The angular distance.  Units of the value returned depend on the
        input value of `units`.

    Notes
    -----
    The formula below is the one best suited to handling small angular
    separations.  See:
    http://en.wikipedia.org/wiki/Great-circle_distance
    """
    from numpy import arcsin, cos, deg2rad, rad2deg, sin, sqrt
    if units == 0:
        rarad1 = ra1
        dcrad1 = dec1
        rarad2 = ra2
        dcrad2 = dec2
    elif units == 1:
        rarad1 = deg2rad(15.0*ra1)
        dcrad1 = deg2rad(dec1)
        rarad2 = deg2rad(15.0*ra2)
        dcrad2 = deg2rad(dec2)
    elif units == 2:
        rarad1 = deg2rad(ra1)
        dcrad1 = deg2rad(dec1)
        rarad2 = deg2rad(ra2)
        dcrad2 = deg2rad(dec2)
    else:
        raise ValueError('units must be 0, 1 or 2!')
    deldec2 = (dcrad2-dcrad1)/2.0
    delra2 = (rarad2-rarad1)/2.0
    sindis = sqrt(sin(deldec2)*sin(deldec2) +
                  cos(dcrad1)*cos(dcrad2)*sin(delra2)*sin(delra2))
    dis = 2.0*arcsin(sindis)
    if units == 0:
        return dis
    else:
        return rad2deg(dis)*3600.0
### Above part are imported from pydl spheregroup


