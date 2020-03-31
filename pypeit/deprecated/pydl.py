# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit import bspline
from pypeit.core import basis


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

    # ToDO Consider refactoring the argument list so that there are no kwargs
    def __init__(self, x, fullbkpt = None, nord=4, npoly=1, bkpt=None, bkspread=1.0,
                 verbose=False, from_dict=None, **kwargs):
        """Init creates an object whose attributes are similar to the
        structure returned by the create_bspline function.
        """
        # JFH added this to enforce immutability of these input arguments, as this code modifies bkpt and fullbkpt
        # as it goes
        fullbkpt1 = copy.copy(fullbkpt)
        bkpt1 = copy.copy(bkpt)
        if from_dict is not None:
            self.nord=from_dict['nord']
            self.npoly=from_dict['npoly']
            self.breakpoints=np.array(from_dict['breakpoints'])
            self.mask=np.array(from_dict['mask'])
            self.coeff=np.array(from_dict['coeff'])
            self.icoeff=np.array(from_dict['icoeff'])
            self.xmin=from_dict['xmin']
            self.xmax=from_dict['xmax']
            self.funcname=from_dict['funcname']
            return
        # Instantiate empty if neither fullbkpt or x is set
        elif x is None and fullbkpt is None:
            self.nord = None
            self.npoly = None
            self.breakpoints= None
            self.mask= None
            self.coeff= None
            self.icoeff= None
            self.xmin= None
            self.xmax= None
            self.funcname= None
            return
        else:
            #
            # Set the breakpoints.
            #
            if fullbkpt1 is None:
                if bkpt1 is None:
                    startx = x.min()
                    rangex = x.max() - startx
                    if 'placed' in kwargs:
                        w = ((kwargs['placed'] >= startx) &
                             (kwargs['placed'] <= startx+rangex))
                        if w.sum() < 2:
                            bkpt1 = np.arange(2, dtype='f') * rangex + startx
                        else:
                            bkpt1 = kwargs['placed'][w]
                    elif 'bkspace' in kwargs:
                        nbkpts = int(rangex/kwargs['bkspace']) + 1
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt1 = np.arange(nbkpts, dtype='f')*tempbkspace + startx
                    elif 'nbkpts' in kwargs:
                        nbkpts = kwargs['nbkpts']
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt1 = np.arange(nbkpts, dtype='f') * tempbkspace + startx
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
                        bkpt1 = np.interp(xspot,np.arange(nx),x)
                    else:
                        raise ValueError('No information for bkpts.')
                # JFH added this new code, because bkpt.size = 1 implies fullbkpt has only 2*(nord-1) + 1 elements.
                # This will cause a crash in action because nbkpt < 2*nord, i.e. for bkpt = 1, nord = 4 fullbkpt has
                # seven elements which is less than 2*nord = 8. The codes above seem to require nbkpt >=2, so I'm implementing
                # this requirement. Note that the previous code before this fix simply sets bkpt to bkpt[imax] =x.max()
                # which is equally arbitrary, but still results in a crash. By requiring at least 2 bkpt, fullbkpt will
                # have 8 elements preventing action from crashing
                if (bkpt1.size < 2):
                    bkpt1 = np.zeros(2,dtype=float)
                    bkpt1[0] = x.min()
                    bkpt1[1] = x.max()
                else:
                    imin = bkpt1.argmin()
                    imax = bkpt1.argmax()
                    if x.min() < bkpt1[imin]:
                        if verbose:
                            print('Lowest breakpoint does not cover lowest x value: changing.')
                        bkpt1[imin] = x.min()
                    if x.max() > bkpt1[imax]:
                        if verbose:
                            print('Highest breakpoint does not cover highest x value: changing.')
                        bkpt1[imax] = x.max()

                nshortbkpt = bkpt1.size
                fullbkpt1 = bkpt1.copy()
                # Note that with the JFH change above, this nshortbkpt ==1 is never realized beacause above I forced
                # bkpt to have at least two elements. Not sure why this was even allowed, since bkpt.size = 1
                #  basically results in action crashing as described above.
                if nshortbkpt == 1:
                    bkspace = bkspread
                else:
                    bkspace = (bkpt1[1] - bkpt1[0])*bkspread
                for i in np.arange(1, nord):
                    fullbkpt1 = np.insert(fullbkpt1, 0, bkpt1[0]-bkspace*i)
                    fullbkpt1 = np.insert(fullbkpt1, fullbkpt1.shape[0],
                                         bkpt1[nshortbkpt-1] + bkspace*i)


            # JFH added this to fix bug in cases where fullbkpt is passed in but has < 2*nord elements
            if fullbkpt1.size < 2*nord:
                fullbkpt_init = fullbkpt1.copy()
                nshortbkpt = fullbkpt_init.size
                bkspace = (fullbkpt_init[1] - fullbkpt_init[0])*bkspread
                for i in np.arange(1, nord):
                    fullbkpt1 = np.insert(fullbkpt1, 0, fullbkpt_init[0] - bkspace * i)
                    fullbkpt1 = np.insert(fullbkpt1, fullbkpt1.shape[0],
                                          fullbkpt_init[nshortbkpt - 1] + bkspace * i)

            nc = fullbkpt1.size - nord
            self.breakpoints = fullbkpt1
            self.nord = nord
            self.npoly = npoly
            self.mask = np.ones((fullbkpt1.size,), dtype='bool')
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

        return

    def copy(self):

        bsp_copy = bspline(None)
        bsp_copy.nord = self.nord
        bsp_copy.npoly = self.npoly
        bsp_copy.breakpoints = np.copy(self.breakpoints)
        bsp_copy.mask = np.copy(self.mask)
        bsp_copy.coeff = np.copy(self.coeff)
        bsp_copy.icoeff = np.copy(self.icoeff)
        bsp_copy.xmin = self.xmin
        bsp_copy.xmax = self.xmax
        bsp_copy.funcname = self.funcname
        return bsp_copy

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
        nbkpt = self.mask.sum()
        if nbkpt < 2*self.nord:
            msgs.warn('Order ({0}) too low for {1} breakpoints.'.format(self.nord, nbkpt))
            return -2, 0, 0
        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros((n - self.nord + 1,), dtype=int)
        upper = np.zeros((n - self.nord + 1,), dtype=int) - 1
        indx = self.intrv(x)
        bf1 = self.bsplvn(x, indx)
        aa = uniq(indx)
        upper[indx[aa]-self.nord+1] = aa
        rindx = indx[::-1]
        bb = uniq(rindx)
        lower[rindx[bb]-self.nord+1] = nx - bb - 1
        if x2 is None:
            return bf1, lower, upper

        if x2.size != nx:
            raise ValueError('Dimensions of x and x2 do not match.')

        # TODO: Below is unchanged.
        x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
        # TODO: Should consider faster ways of generating the temppoly arrays for poly and poly1
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

        # TODO: Should consider faster way of calculating action that doesn't require a nested loop.
        bw = self.npoly*self.nord
        action = np.zeros((nx, bw), dtype='d')
        counter = -1
        for ii in range(self.nord):
            for jj in range(self.npoly):
                counter += 1
                action[:, counter] = bf1[:, ii]*temppoly[:, jj]
        return action, lower, upper

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
        indx = np.zeros(x.size, dtype=int)
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
        # TODO: Is the sorting necessary?
        xsort = x.argsort()
        if action is None:
            action, lower, upper = self.action(x[xsort], x2=None if x2 is None else x2[xsort])
        else:
            if lower is None or upper is None:
                raise ValueError('Must specify lower and upper if action is set.')

        # TODO: Can we save some of these objects to self so that we
        # don't have to recreate them?
        yfit = np.zeros(x.shape, dtype=x.dtype)
        bw = self.npoly * self.nord
        spot = np.arange(bw, dtype=int)
        goodbk = self.mask.nonzero()[0]
        coeffbk = self.mask[self.nord:].nonzero()[0]
        goodcoeff = self.coeff[...,coeffbk]

        nowidth = np.invert(upper+1 > lower)
        n = self.mask.sum() - self.nord
        for i in range(n-self.nord+1):
            if nowidth[i]:
                continue
            yfit[lower[i]:upper[i]+1] = np.dot(action[lower[i]:upper[i]+1,:],
                                               goodcoeff.flatten('F')[i*self.npoly+spot])

        mask = np.ones(x.shape, dtype='bool')
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord-1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size == 0:
            return yfit[np.argsort(xsort)], mask

        for jj in range(hmm.size):
            mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
                    & (x <= self.breakpoints[goodbk[hmm[jj]+1]-1])] = False
        return yfit[np.argsort(xsort)], mask

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
#        # maskthis = np.zeros(xwork.shape,dtype=xwork.dtype)
#        nowidth = np.invert(upper > lower)
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
            msgs.warn('Fewer good break points than order of b-spline. Returning...')
            return -2
        # Find the unique ones for the polynomial
        hmm = err[uniq(err//self.npoly)]//self.npoly

        n = nbkpt - self.nord
        if np.any(hmm >= n):
            msgs.warn('Note enough unique points in cholesky_band decomposition of b-spline matrix. Returning...')
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

    def workit(self, xdata, ydata, invvar, action, lower, upper):
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
        :func:`tuple` (success, yfit):
            A tuple containing an boolean error code, and the evaluation
            of the b-spline yfit at the input values.  The error codes
            are as follows: 0 is good; -1 is dropped breakpoints, try
            again; -2 is failure, should abort.

        """
        goodbk = self.mask[self.nord:]
        # KBW: Interesting: x.sum() is actually a bit faster than np.sum(x)
        nn = goodbk.sum()
        if nn < self.nord:
            msgs.warn('Fewer good break points than order of b-spline. Returning...')
            # KBW: Why is the dtype set to 'f' = np.float32?
            return -2, np.zeros(ydata.shape, dtype=np.float32)

        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        a2 = action * np.sqrt(invvar)[:,None]

        alpha = np.zeros((bw, nfull+bw), dtype=float)
        beta = np.zeros((nfull+bw,), dtype=float)
        bi = np.concatenate([np.arange(i)+(bw-i)*(bw+1) for i in range(bw,0,-1)])
        bo = np.concatenate([np.arange(i)+(bw-i)*bw for i in range(bw,0,-1)])
        upper += 1
        nowidth = np.invert(upper > lower)
        for k in range(nn-self.nord+1):
            if nowidth[k]:
                continue
            itop = k*self.npoly
            alpha.T.flat[bo+itop*bw] \
                    += np.dot(a2[lower[k]:upper[k],:].T, a2[lower[k]:upper[k],:]).flat[bi]
            beta[itop:min(itop,nfull)+bw] \
                    += np.dot(ydata[lower[k]:upper[k]] * np.sqrt(invvar[lower[k]:upper[k]]),
                              a2[lower[k]:upper[k],:])
        upper -= 1

        # Right now we are not returning the covariance, although it may arise that we should
#        covariance = alpha
        err, a = cholesky_band(alpha, mininf=1.0e-10 * invvar.sum() / nfull)

        # successful cholseky_band returns -1
        if not isinstance(err, int) or err != -1:
            return self.maskpoints(err), \
                        self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)[0]

        # NOTE: cholesky_solve ALWAYS returns err == -1; don't even catch it.
        sol = cholesky_solve(a, beta)[1]

        if self.coeff.ndim == 2:
            self.icoeff[:,goodbk] = np.array(a[0,:nfull].T.reshape(self.npoly, nn, order='F'), dtype=a.dtype)
            self.coeff[:,goodbk] = np.array(sol[:nfull].T.reshape(self.npoly, nn, order='F'), dtype=sol.dtype)
        else:
            self.icoeff[goodbk] = np.array(a[0,:nfull], dtype=a.dtype)
            self.coeff[goodbk] = np.array(sol[:nfull], dtype=sol.dtype)

        return 0, self.value(xdata, x2=xdata, action=action, upper=upper, lower=lower)[0]



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

    # KBW: added isfinite check so that it doesn't need to be done in
    # the for loop. This should be okay because sqrt and division of
    # positive numbers by other positive numbers should never lead to
    # non-finite numbers. However, need to watch out for numerical
    # issues.
    bw, nn = l.shape
    n = nn - bw
    negative = (l[0,:n] <= mininf) | np.invert(np.isfinite(l[0,:n]))
    # JFH changed this below to make it more consistent with IDL version. Not sure
    # why the np.all(np.isfinite(lower)) was added. The code could return an empty
    # list for negative.nonzero() and crash if all elements in lower are NaN.
    if negative.any():
        nz = negative.nonzero()[0]
        msgs.warn('Found {0} bad entries: {1}'.format(nz.size, nz))
        return nz, l

#    negative = (lower[0, 0:n] <= mininf)
#    if negative.any() or not np.all(np.isfinite(lower)):
#        msgs.warn('Found {:d}'.format(len(negative.nonzero()[0])) +
#                  ' bad entries: ' + str(negative.nonzero()[0]))
#        return (negative.nonzero()[0], l)

    # KBW: Moved the copy to after the initial check of the input
    lower = l.copy()
    kn = bw - 1
    spot = np.arange(kn, dtype=int) + 1

    # KBW: Faster by about factor of ~2 compared to previous version
    #bi = np.arange(kn*kn).reshape(kn,kn)[np.triu_indices(kn)]
    bi = np.concatenate([np.arange(i)+(kn-i)*(kn+1) for i in range(kn,0,-1)])
    here = bi[:,None] + (np.arange(n)[None,:] + 1)*bw
    for j in range(n):
        lower[0,j] = np.sqrt(lower[0,j])
        lower[spot, j] /= lower[0,j]
        hmm = lower[spot,j,None] * lower[None,spot,j]
        lower.T.flat[here[:,j]] -= hmm.flat[bi]

    return -1, lower


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
    n = b.shape[0] - a.shape[0]
    kn = a.shape[0] - 1

    spot = np.arange(kn, dtype=int) + 1
    for j in range(n):
        b[j] /= a[0,j]
        b[j+spot] -= b[j]*a[spot,j]

    spot = spot[::-1]
    for j in range(n-1, -1, -1):
        b[j] = (b[j] - np.sum(a[spot,j] * b[j+spot]))/a[0,j]

    return -1, b


# TODO: How important is it that the function return the last
# occurrence of the unique values in the sorted array than the first
# value? If it doesn't `numpy.unique(numpy.sort(x),
# return_index=True)[1]` is about twice as fast as
# `pydl.uniq(numpy.sort(x))`.
def old_uniq(x, index=None):
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
    if index is None:
        indicies = (x != np.roll(x, -1)).nonzero()[0]
        if indicies.size > 0:
            return indicies
        else:
            return np.array([x.size - 1, ])
    else:
        q = x[index]
        indicies = (q != np.roll(q, -1)).nonzero()[0]
        if indicies.size > 0:
            return index[indicies]
        else:
            return np.array([q.size - 1, ], dtype=index.dtype)

# Faster than previous version but not as fast as if we could switch to
# np.unique.
def uniq(x, index=None):
    """
    Return the indices of the *last* occurrence of the unique
    elements in a sorted array.

    The input vector must be sorted before being passed to this
    function. This can be done by sorting ``x`` directly or by
    passing the array that sorts ``x`` (``index``).

    Replicates the IDL ``UNIQ()`` function.

    Parameters
    ----------
    x : array-like
        Search this array for unique items.
    index : array-like, optional
        This array provides the array subscripts that sort `x`. That
        is::

            index = np.argsort(x)

    Returns
    -------
    `np.ndarray`
        The indices of the last occurence in `x` of its unique
        values.

    Notes
    -----
    Given a sorted array, and assuming that there is a set of
    adjacent identical items, ``uniq()`` will return the subscript of
    the *last* unique item. This charming feature is retained for
    reproducibility.

    References
    ----------
    http://www.harrisgeospatial.com/docs/uniq.html

    Speed improvement thanks to discussion here:
    https://stackoverflow.com/questions/47495510/numpy-in-a-sorted-list-find-the-first-and-the-last-index-for-each-unique-value

    Examples
    --------
    >>> import numpy as np
    >>> from pydl import uniq
    >>> data = np.array([ 1, 2, 3, 1, 5, 6, 1, 7, 3, 2, 5, 9, 11, 1 ])
    >>> print(uniq(np.sort(data)))
    [ 3  5  7  9 10 11 12 13]
    """
    if len(x) == 0:
        raise ValueError('No unique elements in an empty array!')
    if index is None:
        return np.flatnonzero(np.concatenate(([True], x[1:] != x[:-1], [True])))[1:]-1
    _x = x[index]
    return np.flatnonzero(np.concatenate(([True], _x[1:] != _x[:-1], [True])))[1:]-1


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


def fchebyshev_split(x, m):
    """Compute the first `m` Chebyshev polynomials, but modified to allow a
    split in the baseline at :math:`x=0`.  The intent is to allow a model fit
    where a constant term is different for positive and negative `x`.

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
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 2:
        raise ValueError('Order of polynomial must be at least 2.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    try:
        leg[0, :] = (x >= 0).astype(x.dtype)
    except AttributeError:
        leg[0, :] = np.double(x >= 0)
    if m > 2:
        leg[2, :] = x
    if m > 3:
        for k in range(3, m):
            leg[k, :] = 2.0 * x * leg[k-1, :] - leg[k-2, :]
    return leg



def fpoly(x, m):
    """Compute the first `m` simple polynomials.

    Parameters
    ----------
    x : array-like
        Compute the simple polynomials at these abscissa values.
    m : :class:`int`
        The number of simple polynomials to compute.  For example, if
        :math:`m = 3`, :math:`x^0`, :math:`x^1` and
        :math:`x^2` will be computed.

    Returns
    -------
    :class:`numpy.ndarray`
    """
    if isinstance(x, np.ndarray):
        n = x.size
    else:
        n = 1
    if m < 1:
        raise ValueError('Order of polynomial must be at least 1.')
    try:
        dt = x.dtype
    except AttributeError:
        dt = np.float64
    leg = np.ones((m, n), dtype=dt)
    if m >= 2:
        leg[1, :] = x
    if m >= 3:
        for k in range(2, m):
            leg[k, :] = leg[k-1, :] * x
    return leg


