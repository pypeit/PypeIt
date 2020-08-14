# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL

"""
Implements the bspline class

.. include:: ../links.rst
"""

import copy
import warnings

from IPython import embed

import numpy as np
from matplotlib import pyplot as plt

from pypeit.core import basis
from pypeit.core import pydl
from pypeit import datamodel
from pypeit import msgs

try:
    from pypeit.bspline.utilc import cholesky_band, cholesky_solve, solution_arrays, intrv, \
                                     bspline_model
except:
    warnings.warn('Unable to load bspline C extension.  Try rebuilding pypeit.  In the '
                  'meantime, falling back to pure python code.')
    from pypeit.bspline.utilpy import cholesky_band, cholesky_solve, solution_arrays, intrv, \
                                        bspline_model

# TODO: Used for testing.  Keep around for now.
#from pypeit.bspline.utilpy import bspline_model
#from pypeit.bspline.utilpy import cholesky_band, cholesky_solve, solution_arrays, intrv, \
#                                    bspline_model

# TODO: Types are important for the C extension. Types should be
# limited to int, float, bool!
# TODO: May need to add hooks to utilc.py that do the type conversion,
# but that should be a last resort for stability.
# TODO: This whole module needs to be cleaned up.


class bspline(datamodel.DataContainer):
    """Bspline class.

    Functions in the bspline library are implemented as methods on this
    class.

    Parameters
    ----------
    x : `numpy.ndarray`_
        The data.
    nord : :class:`int`, optional
        To be documented.
    npoly : :class:`int`, optional
        To be documented.
    bkpt : `numpy.ndarray`_, optional
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
    version = '1.0.0'

    datamodel = {
        'breakpoints':  dict(otype=np.ndarray, atype=np.floating, desc='Breakpoint locations'),
        'nord': dict(otype=int, desc='Order of the bspline fit'),
        'npoly': dict(otype=int, desc='Order of the bspline polynomial'),
        'mask': dict(otype=np.ndarray, atype=np.bool_, desc='Mask'),
        'coeff': dict(otype=np.ndarray, atype=np.floating, desc='Fit coefficients'),
        'icoeff': dict(otype=np.ndarray, atype=np.floating, desc='??'),
        'xmin': dict(otype=float, desc='Normalization for input data'),
        'xmax': dict(otype=float, desc='Normalization for input data'),
        'funcname': dict(otype=str, desc='Function of fit'),
    }

    # ToDO Consider refactoring the argument list so that there are no kwargs
    def __init__(self, x, fullbkpt=None, nord=4, npoly=1, bkpt=None, bkspread=1.0, verbose=False,
                 from_dict=None, **kwargs):
        """Init creates an object whose attributes are similar to the
        structure returned by the create_bspline function.
        """
        # Setup the DataContainer with everything None
        datamodel.DataContainer.__init__(self)
        # JFH added this to enforce immutability of these input arguments, as this code modifies bkpt and fullbkpt
        # as it goes
        fullbkpt1 = copy.copy(fullbkpt)
        bkpt1 = copy.copy(bkpt)
        if from_dict is not None:
            self.nord=from_dict['nord']
            self.npoly=from_dict['npoly']
            self.breakpoints=np.array(from_dict['breakpoints']).astype(float)   # Force type
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
                            bkpt1 = np.arange(2, dtype=float) * rangex + startx
                        else:
                            bkpt1 = kwargs['placed'][w]
                    elif 'bkspace' in kwargs:
                        nbkpts = int(rangex/kwargs['bkspace']) + 1
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt1 = np.arange(nbkpts, dtype=float)*tempbkspace + startx
                    elif 'nbkpts' in kwargs:
                        nbkpts = kwargs['nbkpts']
                        if nbkpts < 2:
                            nbkpts = 2
                        tempbkspace = rangex/float(nbkpts-1)
                        bkpt1 = np.arange(nbkpts, dtype=float) * tempbkspace + startx
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
                    bkpt1 = np.zeros(2, dtype=float)
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
            self.breakpoints = fullbkpt1.astype(float)      # Ensure type is float for C extension
            self.nord = nord
            self.npoly = npoly
            self.mask = np.ones((fullbkpt1.size,), dtype=bool)
            if npoly > 1:
                self.coeff = np.zeros((npoly, nc), dtype=float)
                self.icoeff = np.zeros((npoly, nc), dtype=float)
            else:
                self.coeff = np.zeros((nc,), dtype=float)
                self.icoeff = np.zeros((nc,), dtype=float)
            self.xmin = 0.0
            self.xmax = 1.0
            self.funcname = kwargs['funcname'] if 'funcname' in kwargs else 'legendre'

    def reinit_coeff(self):
        nc = self.breakpoints.size - self.nord
        self.coeff = np.zeros((self.npoly, nc), dtype=float) if self.npoly > 1 \
                        else np.zeros(nc, dtype=float)

    def _init_internals(self):
        self.hdu_prefix = None

    def _bundle(self):
        """
        Overload for the HDU name

        Returns:
            list:

        """
        return super(bspline, self)._bundle(ext='BSPLINE')


    def copy(self):
        """
        Return a copied instance of the object.
        """
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
        """
        Write bspline attributes to a dict.

        Attributes returned are: :attr:`breakpoints`, :attr:`nord`,
        :attr:`npoly`, :attr:`mask`, :attr:`coeff`, :attr:`icoeff`,
        :attr:`xmin`, :attr:`xmax`, and :attr:`funcname`.

        .. note::

            `numpy.ndarray`_ objects are converted to lists in the
            dictionary to make it JSON compatible.

        Returns:
            :obj:`dict`: A dictionary with the above keys and items.

        """
        return dict(breakpoints=self.breakpoints.tolist(),
                    nord=self.nord,
                    npoly=self.npoly,
                    mask=self.mask.tolist(),
                    coeff=self.coeff.tolist(),
                    icoeff=self.icoeff.tolist(),
                    xmin=self.xmin,
                    xmax=self.xmax,
                    funcname=self.funcname)

    # TODO: C this
    # TODO: Should this be used, or should we effectively replace it
    # with the content of utils.bspline_profile
    def fit(self, xdata, ydata, invvar, x2=None):
        """Calculate a B-spline in the least-squares sense.

        Fit is based on two variables: x which is sorted and spans a large range
        where bkpts are required y which can be described with a low order
        polynomial.

        Parameters
        ----------
        xdata : `numpy.ndarray`_
            Independent variable.
        ydata : `numpy.ndarray`_
            Dependent variable.
        invvar : `numpy.ndarray`_
            Inverse variance of `ydata`.
        x2 : `numpy.ndarray`_, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :obj:`tuple`
            A tuple containing an integer error code, and the evaluation of the
            b-spline at the input values.  An error code of -2 is a failure,
            -1 indicates dropped breakpoints, 0 is success, and positive
            integers indicate ill-conditioned breakpoints.
        """
        goodbk = self.mask[self.nord:]
        nn = goodbk.sum()
        if nn < self.nord:
            yfit = np.zeros(ydata.shape, dtype=float)
            return (-2, yfit)
        nfull = nn * self.npoly
        bw = self.npoly * self.nord
        a1, lower, upper = self.action(xdata, x2=x2)
        foo = np.tile(invvar, bw).reshape(bw, invvar.size).transpose()
        a2 = a1 * foo
        alpha = np.zeros((bw, nfull+bw), dtype=float)
        beta = np.zeros((nfull+bw,), dtype=float)
        bi = np.arange(bw, dtype=int)
        bo = np.arange(bw, dtype=int)
        for k in range(1, bw):
            bi = np.append(bi, np.arange(bw-k, dtype=int)+(bw+1)*k)
            bo = np.append(bo, np.arange(bw-k, dtype=int)+bw*k)
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
        x : `numpy.ndarray`_
            Independent variable.
        x2 : `numpy.ndarray`_, optional
            Orthogonal dependent variable for 2d fits.

        Returns
        -------
        :obj:`tuple`
            A tuple containing the b-spline action matrix; the 'lower' parameter,
            a list of pixel positions, each corresponding to the first
            occurence of position greater than breakpoint indx; and 'upper',
            Same as lower, but denotes the upper pixel positions.
        """
        nbkpt = self.mask.sum()
        if nbkpt < 2*self.nord:
            warnings.warn('Order ({0}) too low for {1} breakpoints.'.format(self.nord, nbkpt))
            return -2, 0, 0
        nx = x.size
        n = nbkpt - self.nord
        lower = np.zeros((n - self.nord + 1,), dtype=int)
        upper = np.zeros((n - self.nord + 1,), dtype=int) - 1
        indx = intrv(self.nord, self.breakpoints[self.mask], x)
        bf1 = self.bsplvn(x, indx)
#        print('F_CONTIGUOUS after bsplvn: {0}'.format(bf1.flags['F_CONTIGUOUS']))
        aa = uniq(indx)
        upper[indx[aa]-self.nord+1] = aa
        rindx = indx[::-1]
        bb = uniq(rindx)
        lower[rindx[bb]-self.nord+1] = nx - bb - 1
        if x2 is None:
            return bf1, lower, upper

#        print('x2!')

        if x2.size != nx:
            raise ValueError('Dimensions of x and x2 do not match.')

        # TODO: Below is unchanged.
        x2norm = 2.0 * (x2 - self.xmin) / (self.xmax - self.xmin) - 1.0
        # TODO: Should consider faster ways of generating the temppoly arrays for poly and poly1
        if self.funcname == 'poly':
            temppoly = np.ones((nx, self.npoly), dtype=float)
            for i in range(1, self.npoly):
                temppoly[:, i] = temppoly[:, i-1] * x2norm
        elif self.funcname == 'poly1':
            temppoly = np.tile(x2norm, self.npoly).reshape(nx, self.npoly)
            for i in range(1, self.npoly):
                temppoly[:, i] = temppoly[:, i-1] * x2norm
        elif self.funcname == 'chebyshev':
            # JFH fixed bug here where temppoly needed to be transposed because of different IDL and python array conventions
            # NOTE: Transposed them in the functions themselves
            temppoly = basis.fchebyshev(x2norm, self.npoly)
        elif self.funcname == 'legendre':
            temppoly = basis.flegendre(x2norm, self.npoly)
        else:
            raise ValueError('Unknown value of funcname.')

        # TODO: Should consider faster way of calculating action that
        # doesn't require a nested loop. Below might work, but it needs
        # to be tested.
#        _action = (bf1[:,:,None] * temppoly[:,None,:]).reshape(nx,-1)
        bw = self.npoly*self.nord
        action = np.zeros((nx, bw), dtype=float, order='F')
        counter = -1
        for ii in range(self.nord):
            for jj in range(self.npoly):
                counter += 1
                action[:, counter] = bf1[:, ii]*temppoly[:, jj]
        return action, lower, upper


    # TODO: C this?
    def bsplvn(self, x, ileft):
        """To be documented.

        Parameters
        ----------
        x : `numpy.ndarray`_
            To be documented.
        ileft : :class:`int`
            To be documented

        Returns
        -------
        vnikx : `numpy.ndarray`_
            To be documented.
        """
        bkpt = self.breakpoints[self.mask]
        # TODO: Had to set the order here to keep it consistent with
        # utils.bspline_profile, but is this going to break things
        # elsewhere? Ideally, we wouldn't be setting the memory order
        # anywhere...
        vnikx = np.zeros((x.size, self.nord), dtype=x.dtype, order='F')
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
        x : `numpy.ndarray`_
            Independent variable.
        x2 : `numpy.ndarray`_, optional
            Orthogonal dependent variable for 2d fits.
        action : `numpy.ndarray`_, optional
            Action matrix to use.  If not supplied it is calculated.
        lower : `numpy.ndarray`_, optional
            If the action parameter is supplied, this parameter must also
            be supplied.
        upper : `numpy.ndarray`_, optional
            If the action parameter is supplied, this parameter must also
            be supplied.

        Returns
        -------
        yfit : `numpy.ndarray`_
            Results of the bspline evaluation
        mask : `numpy.ndarray`_
            Mask indicating where the evaluation was good (i.e., True
            is good).
        """
        # TODO: Is the sorting necessary?
        xsort = x.argsort()
        if action is None:
            action, lower, upper = self.action(x[xsort], x2=None if x2 is None else x2[xsort])
        else:
            if lower is None or upper is None:
                raise ValueError('Must specify lower and upper if action is set.')

        n = self.mask.sum() - self.nord
        coeffbk = self.mask[self.nord:].nonzero()[0]
        goodcoeff = self.coeff[...,coeffbk]
        yfit = bspline_model(x, action, lower, upper, goodcoeff, n, self.nord, self.npoly)

        mask = np.ones(x.shape, dtype=bool)
        goodbk = self.mask.nonzero()[0]
        gb = self.breakpoints[goodbk]
        mask[(x < gb[self.nord-1]) | (x > gb[n])] = False
        hmm = (np.diff(goodbk) > 2).nonzero()[0]
        if hmm.size == 0:
            return yfit[np.argsort(xsort)], mask

        for jj in range(hmm.size):
            mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
                    & (x <= self.breakpoints[goodbk[hmm[jj]+1]-1])] = False
        return yfit[np.argsort(xsort)], mask

    def maskpoints(self, err):
        """Perform simple logic of which breakpoints to mask.


        Parameters
        ----------
        err : `numpy.ndarray`_, :obj:`int`
            The list of indexes returned by the cholesky routines.
            This is indexed to the set of currently *good*
            breakpoints (i.e. self.mask=True) and the first nord are
            skipped.

        Returns
        -------
        :obj:`int`
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
            warnings.warn('Fewer good break points than order of b-spline. Returning...')
            return -2
        # Find the unique ones for the polynomial
        hmm = err[uniq(err//self.npoly)]//self.npoly

        n = nbkpt - self.nord
        if np.any(hmm >= n):
            warnings.warn('Note enough unique points in cholesky_band decomposition of b-spline matrix. Returning...')
            return -2
        test = np.zeros(nbkpt, dtype=bool)
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
            return -2
        return -2

    def workit(self, xdata, ydata, invvar, action, lower, upper):
        """An internal routine for bspline_extract and bspline_radial which solve a general
        banded correlation matrix which is represented by the variable "action".  This routine
        only solves the linear system once, and stores the coefficients in sset. A non-zero return value
        signifies a failed inversion


        Parameters
        ----------
        xdata : `numpy.ndarray`_
            Independent variable.
        ydata : `numpy.ndarray`_
            Dependent variable.
        invvar : `numpy.ndarray`_
            Inverse variance of `ydata`.
        action : `numpy.ndarray`_
            Banded correlation matrix
        lower  : `numpy.ndarray`_
            A list of pixel positions, each corresponding to the first occurence of position greater than breakpoint indx
        upper  : `numpy.ndarray`_
            Same as lower, but denotes the upper pixel positions

        Returns
        -------
        success : :obj:`int`
            Method error code: 0 is good; -1 is dropped breakpoints,
            try again; -2 is failure, should abort.
        yfit : `numpy.ndarray`_
            Evaluation of the b-spline yfit at the input values.
        """
        goodbk = self.mask[self.nord:]
        # KBW: Interesting: x.sum() is actually a bit faster than np.sum(x)
        nn = goodbk.sum()
        if nn < self.nord:
            warnings.warn('Fewer good break points than order of b-spline. Returning...')
            # KBW: Why is the dtype set to 'f' = np.float32?
            return -2, np.zeros(ydata.shape, dtype=float)

        alpha, beta = solution_arrays(nn, self.npoly, self.nord, ydata, action, invvar, upper,
                                      lower)
        nfull = nn * self.npoly

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

# TODO: I don't think we need to make this reproducible with the IDL version anymore, and can opt for speed instead.
# TODO: Move this somewhere for more common access?
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
    >>> from pypeit.core.pydl import uniq
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


def iterfit(xdata, ydata, invvar=None, inmask=None, upper=5, lower=5, x2=None,
            maxiter=10, nord=4, bkpt=None, fullbkpt=None, kwargs_bspline={}, kwargs_reject={}):
        """Iteratively fit a b-spline set to data, with rejection. This is a utility function that allows
        the bspline to use via a direct function call.

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
        # from .math import djs_reject
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
            var = ydata.var() * (float(nx) / float(nx - 1))
            if var == 0:
                var = 1.0
            invvar = np.ones(ydata.shape, dtype=ydata.dtype) / var

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
                sset = bspline.bspline(xdata[xsort[maskwork]], nord=nord, bkpt=bkpt, fullbkpt=fullbkpt, **kwargs_bspline)
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
                iiter = maxiter + 1  # End iterations
            else:
                if 'requiren' in kwargs_bspline:
                    i = 0
                    while xwork[i] < sset.breakpoints[goodbk[sset.nord]] and i < nx - 1:
                        i += 1
                    ct = 0
                    for ileft in range(sset.nord, sset.mask.sum() - sset.nord + 1):
                        while (xwork[i] >= sset.breakpoints[goodbk[ileft]] and
                               xwork[i] < sset.breakpoints[goodbk[ileft + 1]] and
                               i < nx - 1):
                            ct += invwork[i] * maskwork[i] > 0
                            i += 1
                        if ct >= kwargs_bspline['requiren']:
                            ct = 0
                        else:
                            sset.mask[goodbk[ileft]] = False
                error, yfit = sset.fit(xwork, ywork, invwork * maskwork,
                                       x2=x2work)
            iiter += 1
            inmask_rej = maskwork
            if error == -2:

                return (sset, outmask)
            elif error == 0:
                # ToDO JFH by setting inmask to be tempin which is maskwork, we are basically implicitly enforcing sticky rejection
                # here. See djs_reject.py. I'm leaving this as is for consistency with the IDL version, but this may require
                # further consideration. I think requiring stick to be set is the more transparent behavior.
                maskwork, qdone = pydl.djs_reject(ywork, yfit, invvar=invwork, inmask=inmask_rej, outmask=maskwork,
                                             upper=upper, lower=lower, **kwargs_reject)
            else:
                pass
        outmask[xsort] = maskwork
        temp = yfit
        yfit[xsort] = temp
        return (sset, outmask)


def bspline_profile(xdata, ydata, invvar, profile_basis, ingpm=None, upper=5, lower=5, maxiter=25,
                    nord=4, bkpt=None, fullbkpt=None, relative=None, kwargs_bspline={},
                    kwargs_reject={}, quiet=False):
    """
    Fit a B-spline in the least squares sense with rejection to the
    provided data and model profiles.

    .. todo::
        Fully describe procedure.

    Parameters
    ----------
    xdata : `numpy.ndarray`_
        Independent variable.
    ydata : `numpy.ndarray`_
        Dependent variable.
    invvar : `numpy.ndarray`_
        Inverse variance of `ydata`.
    profile_basis : `numpy.ndarray`_
        Model profiles.
    ingpm : `numpy.ndarray`_, optional
        Input good-pixel mask. Values to fit in ``ydata`` should be
        True.
    upper : :obj:`int`, :obj:`float`, optional
        Upper rejection threshold in units of sigma, defaults to 5
        sigma.
    lower : :obj:`int`, :obj:`float`, optional
        Lower rejection threshold in units of sigma, defaults to 5
        sigma.
    maxiter : :obj:`int`, optional
        Maximum number of rejection iterations, default 10. Set this
        to zero to disable rejection.
    nord : :obj:`int`, optional
        Order of B-spline fit
    bkpt : `numpy.ndarray`_, optional
        Array of breakpoints to be used for the b-spline
    fullbkpt : `numpy.ndarray`_, optional
        Full array of breakpoints to be used for the b-spline,
        without letting the b-spline class append on any extra bkpts
    relative : `numpy.ndarray`_, optional
        Array of integer indices to be used for computing the reduced
        chi^2 of the fits, which then is used as a scale factor for
        the upper,lower rejection thresholds
    kwargs_bspline : :obj:`dict`, optional
        Keyword arguments used to instantiate
        :class:`pypeit.bspline.bspline`
    kwargs_reject : :obj:`dict`, optional
        Keyword arguments passed to :func:`pypeit.core.pydl.djs_reject`
    quiet : :obj:`bool`, optional
        Suppress output to the screen

    Returns
    -------
    sset : :class:`pypeit.bspline.bspline`
        Result of the fit.
    gpm : `numpy.ndarray`_
        Output good-pixel mask which the same size as ``xdata``. The
        values in this array for the corresponding data are not used in
        the fit, either because the input data was masked or the data
        were rejected during the fit, if they are False. Data
        rejected during the fit (if rejection is performed) are::

            rejected = ingpm & np.logical_not(gpm)

    yfit : `numpy.ndarray`_
        The best-fitting model; shape is the same as ``xdata``.
    reduced_chi : :obj:`float`
        Reduced chi-square of the best-fitting model.
    exit_status : :obj:`int`
        Indication of the success/failure of the fit.  Values are:

            - 0 = fit exited cleanly
            - 1 = maximum iterations were reached
            - 2 = all points were masked
            - 3 = all break points were dropped
            - 4 = Number of good data points fewer than nord

    """
    # Checks
    nx = xdata.size
    if ydata.size != nx:
        msgs.error('Dimensions of xdata and ydata do not agree.')

    # TODO: invvar and profile_basis should be optional

    # ToDO at the moment invvar is a required variable input
    #    if invvar is not None:
    #        if invvar.size != nx:
    #            raise ValueError('Dimensions of xdata and invvar do not agree.')
    #        else:
    #            #
    #            # This correction to the variance makes it the same
    #            # as IDL's variance()
    #            #
    #            var = ydata.var()*(float(nx)/float(nx-1))
    #            if var == 0:
    #                var = 1.0
    #            invvar = np.ones(ydata.shape, dtype=ydata.dtype)/var

    npoly = int(profile_basis.size/nx)
    if profile_basis.size != nx*npoly:
        msgs.error('Profile basis is not a multiple of the number of data points.')

    # Init
    yfit = np.zeros(ydata.shape)
    reduced_chi = 0.

    # TODO: Instanting these place-holder arrays can be expensive.  Can we avoid doing this?
    outmask = True if invvar.size == 1 else np.ones(invvar.shape, dtype=bool)

    if ingpm is None:
        ingpm = invvar > 0

    if not quiet:
        termwidth = 80-13
        msgs.info('B-spline fit:')
        msgs.info('    npoly = {0} profile basis functions'.format(npoly))
        msgs.info('    ngood = {0}/{1} measurements'.format(np.sum(ingpm), ingpm.size))
        msgs.info(' {0:>4}  {1:>8}  {2:>7}  {3:>6} '.format(
                    'Iter', 'Chi^2', 'N Rej', 'R. Fac').center(termwidth, '*'))
        hlinestr = ' {0}  {1}  {2}  {3} '.format('-'*4, '-'*8, '-'*7, '-'*6)
        nullval = '  {0:>8}  {1:>7}  {2:>6} '.format('-'*2, '-'*2, '-'*2)
        msgs.info(hlinestr.center(termwidth))

    maskwork = outmask & ingpm & (invvar > 0)
    if not maskwork.any():
        msgs.error('No valid data points in bspline_profile!.')

    # Init bspline class
    sset = bspline(xdata[maskwork], nord=nord, npoly=npoly, bkpt=bkpt, fullbkpt=fullbkpt,
                           funcname='Bspline longslit special', **kwargs_bspline)
    if maskwork.sum() < sset.nord:
        if not quiet:
            msgs.warn('Number of good data points fewer than nord.')
        # TODO: Why isn't maskwork returned?
        return sset, outmask, yfit, reduced_chi, 4

    # This was checked in detail against IDL for identical inputs
    # KBW: Tried a few things and this was about as fast as you can get.
    outer = np.outer(np.ones(nord, dtype=float), profile_basis.flatten('F')).T
    action_multiple = outer.reshape((nx, npoly * nord), order='F')
    #--------------------
    # Iterate spline fit
    iiter = 0
    error = -1                  # Indicates that the fit should be done
    qdone = False               # True if rejection iterations are done
    exit_status = 0
    relative_factor = 1.0
    nrel = 0 if relative is None else len(relative)
    # TODO: Why do we need both maskwork and tempin?
    tempin = np.copy(ingpm)
    while (error != 0 or qdone is False) and iiter <= maxiter and exit_status == 0:
        ngood = maskwork.sum()
        goodbk = sset.mask.nonzero()[0]
        if ngood <= 1 or not sset.mask.any():
            sset.coeff = 0
            exit_status = 2 # This will end iterations
        else:
            # Do the fit. Return values from workit for error are as follows:
            #    0 if fit is good
            #   -1 if some breakpoints are masked, so try the fit again
            #   -2 if everything is screwed

            # we'll do the fit right here..............
            if error != 0:
                bf1, laction, uaction = sset.action(xdata)
                if np.any(bf1 == -2) or bf1.size !=nx*nord:
                    msgs.error("BSPLINE_ACTION failed!")
                action = np.copy(action_multiple)
                for ipoly in range(npoly):
                    action[:, np.arange(nord)*npoly + ipoly] *= bf1
                del bf1 # Clear the memory

            if np.any(np.logical_not(np.isfinite(action))):
                msgs.error('Infinities in action matrix.  B-spline fit faults.')

            error, yfit = sset.workit(xdata, ydata, invvar*maskwork, action, laction, uaction)

        iiter += 1

        if error == -2:
            if not quiet:
                msgs.warn('All break points lost!!  Bspline fit failed.')
            exit_status = 3
            return sset, np.zeros(xdata.shape, dtype=bool), np.zeros(xdata.shape), reduced_chi, \
                        exit_status

        if error != 0:
            if not quiet:
                msgs.info((' {0:4d}'.format(iiter) + nullval).center(termwidth))
            continue

        # Iterate the fit -- next rejection iteration
        chi_array = (ydata - yfit)*np.sqrt(invvar * maskwork)
        reduced_chi = np.sum(np.square(chi_array)) / (ngood - npoly*(len(goodbk) + nord)-1)

        relative_factor = 1.0
        if relative is not None:
            this_chi2 = reduced_chi if nrel == 1 \
                            else np.sum(np.square(chi_array[relative])) \
                                    / (nrel - (len(goodbk) + nord) - 1)
            relative_factor = max(np.sqrt(this_chi2), 1.0)

        # Rejection

        # TODO: JFH by setting ingpm to be tempin which is maskwork, we
        #  are basically implicitly enforcing sticky rejection here. See
        #  djs_reject.py. I'm leaving this as is for consistency with
        #  the IDL version, but this may require further consideration.
        #  I think requiring sticky to be set is the more transparent
        #  behavior.
        maskwork, qdone = pydl.djs_reject(ydata, yfit, invvar=invvar, inmask=tempin,
                                          outmask=maskwork, upper=upper*relative_factor,
                                          lower=lower*relative_factor, **kwargs_reject)
        tempin = np.copy(maskwork)
        if not quiet:
            msgs.info(' {0:4d}  {1:8.3f}  {2:7d}  {3:6.2f} '.format(iiter,
                        reduced_chi, np.sum(maskwork == 0), relative_factor).center(termwidth))

    if iiter == (maxiter + 1):
        exit_status = 1

    # Exit status:
    #    0 = fit exited cleanly
    #    1 = maximum iterations were reached
    #    2 = all points were masked
    #    3 = all break points were dropped
    #    4 = Number of good data points fewer than nord

    if not quiet:
        msgs.info(' {0:>4}  {1:8.3f}  {2:7d}  {3:6.2f} '.format('DONE',
                    reduced_chi, np.sum(maskwork == 0), relative_factor).center(termwidth))
        msgs.info('*'*termwidth)

    # Finish
    # TODO: Why not return maskwork directly
    outmask = np.copy(maskwork)
    # Return
    return sset, outmask, yfit, reduced_chi, exit_status



def bspline_qa(xdata, ydata, sset, gpm, yfit, xlabel=None, ylabel=None, title=None, show=True):
    """
    Construct a QA plot of the bspline fit.

    Args:
        xdata (`numpy.ndarray`_):
            Array with the independent variable. Regardless of shape,
            data is treated as one-dimensional.
        ydata (`numpy.ndarray`_):
            Array with the dependent variable. Regardless of shape,
            data is treated as one-dimensional.
        sset (:class:`pypeit.bspline.bspline`):
            Object with the results of the fit. (First object
            returned by :func:`bspline_profile`).
        gpm (`numpy.ndarray`_):
            Boolean array with the same size as ``xdata``.
            Measurements rejected during the fit have ``gpm=False``.
            (Second object returned by :func:`bspline_profile`).
        yfit (`numpy.ndarray`_):
            Best-fitting model sampled at ``xdata``. (Third object
            returned by :func:`bspline_profile`).
        xlabel (:obj:`str`, optional):
            Label for the ordinate.  If None, none given.
        ylabel (:obj:`str`, optional):
            Label for the abcissa.  If None, none given.
        title (:obj:`str`, optional):
            Label for the plot.  If None, none given.
        show (:obj:`bool`, optional):
            Plot the result. If False, the axis instance is returned.
            This is done before any labels or legends are added to
            the plot.

    Returns:
        `matplotlib.axes.Axes`_: Axes instance with the data, model,
        and breakpoints.  Only returned if ``show`` is False.
    """
    goodbk = sset.mask
    bkpt, _ = sset.value(sset.breakpoints[goodbk])
    was_fit_and_masked = np.logical_not(gpm)

    plt.clf()
    ax = plt.gca()
    ax.plot(xdata, ydata, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
            linestyle='None', label='data')
    ax.plot(xdata[was_fit_and_masked], ydata[was_fit_and_masked], color='red', marker='+',
            markersize=1.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
    ax.plot(xdata, yfit, color='cornflowerblue', label='fit')
    ax.plot(sset.breakpoints[goodbk], bkpt, color='lawngreen', marker='o', markersize=2.0,
            mfc='lawngreen', fillstyle='full', linestyle='None', label='bspline breakpoints')
    ax.set_ylim(0.99*np.amin(yfit), 1.01*np.amax(yfit))
    if not show:
        return ax

    plt.legend()
    if xlabel is not None:
        plt.xlabel(xlabel)
    if ylabel is not None:
        plt.ylabel(ylabel)
    if title is not None:
        plt.title(title)
    plt.show()

