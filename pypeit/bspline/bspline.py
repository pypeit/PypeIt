# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL

"""
Implements the bspline class

.. include:: ../include/links.rst

"""

import copy
import warnings

from IPython import embed

import numpy as np

from pypeit.core import basis
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

    The datamodel attributes are:

    .. include:: ../include/class_datamodel_bspline.rst

    When written to an output-file HDU, all `numpy.ndarray`_ elements are
    bundled into an `astropy.io.fits.BinTableHDU`_, and the other elements are
    written as header keywords.  Any datamodel elements that are None are *not*
    included in the output.

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

    # TODO: Fix the description of icoeff
    datamodel = {'breakpoints':  dict(otype=np.ndarray, atype=np.floating,
                                      descr='Breakpoint locations'),
                 'nord': dict(otype=int, descr='Order of the bspline fit'),
                 'npoly': dict(otype=int, descr='Order of the bspline polynomial'),
                 'mask': dict(otype=np.ndarray, atype=np.bool_, descr='Mask'),
                 'coeff': dict(otype=np.ndarray, atype=np.floating, descr='Fit coefficients'),
                 'icoeff': dict(otype=np.ndarray, atype=np.floating, descr='??'),
                 'xmin': dict(otype=float, descr='Normalization for input data'),
                 'xmax': dict(otype=float, descr='Normalization for input data'),
                 'funcname': dict(otype=str, descr='Function of fit')}

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


