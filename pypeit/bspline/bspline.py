# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL

"""
Implements the bspline class

.. include:: ../include/links.rst

"""

import warnings

from IPython import embed

import numpy as np

from pypeit.core import basis
from pypeit import datamodel

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
    x : `numpy.ndarray`_, optional
        Independent variable for the definition of the b-spline.  If None,
        ``fullbkpt`` must be provided.
    fullbkpt : `numpy.ndarray`_, optional
        The full set of breakpoints.  The input vector is sorted and cast as a
        float.  If the length of the vector is less than twice ``nord``, it is
        also padded with ``nord-1`` values, as needed for the construction of
        the b-spline.  If None, ``x`` must be provided.
    nord : int, optional
        Order of the b-spline.
    npoly : int, optional
        Polynomial order to fit over 2nd variable (when specified using ``x2``;
        see :func:`~pypeit.bspline.bspline.bspline.fit`).
    bkpt : `numpy.ndarray`_, optional
        A precalculated set of breakpoints within the range of ``x``.  The input
        vector is sorted and any points beyond the range of x are omitted.  If
        only one breakpoint is provided (or not omitted), one breakpoint is set
        at each end of ``x``.  If the breakpoints do not cover the full range,
        the first and last breakpoints are moved such that they do.  If None,
        the breakpoints are determined using the keywords below.
    bkspread : float, optional
        Scale factor for the separation between breakpoints.
    bkspace : float, optional
        Defines the separation between breakpoints.  If provided, ``nbkpts`` and
        ``everyn`` are ignored.
    nbkpts : int, optional
        Defines the number of breakpoints used to span the full range of ``x``.
        Only used if ``bkspace`` is None.  If provided, ``everyn`` is ignored.
    everyn : int, float, optional
        Places a breakpoint at every Nth value of ``x``.  Only used if
        ``bkspace`` and ``nbkpts`` are both None.
    funcname : str, optional
        Function for the second variable (when specified using ``x2``; see
        :func:`~pypeit.bspline.bspline.bspline.fit`).
    """

    version = '1.0.0'
    """
    Datamodel version
    """

    datamodel = {
        'breakpoints': dict(otype=np.ndarray, atype=np.floating, descr='Breakpoint locations'),
        'nord': dict(otype=int, descr='Order of the bspline fit'),
        'npoly': dict(otype=int,
                      descr='Order of polynomial to fit over 2nd variable (when x2 is specified)'),
        'mask': dict(otype=np.ndarray, atype=np.bool_, descr='Output mask'),
        'coeff': dict(otype=np.ndarray, atype=np.floating, descr='Output fit coefficients'),
        'icoeff': dict(otype=np.ndarray, atype=np.floating,
                       descr='Cholesky band matrix used to solve for the bspline coefficients'),
        'xmin': dict(otype=float, descr='Normalization minimum for x2'),
        'xmax': dict(otype=float, descr='Normalization maximum for x2'),
        'funcname': dict(otype=str,
                         descr='Function type for the 2nd variable (when x2 is specified)'),
    }
    """
    Datamodel components.
    """

    def __init__(self, x=None, fullbkpt=None, nord=4, npoly=1, bkpt=None, bkspread=1.0,
                 bkspace=None, nbkpts=None, everyn=None, funcname='legendre'):

        # Instantiate the base class
        datamodel.DataContainer.__init__(self)

        # Instantiate empty if neither fullbkpt or x is set
        if x is None and fullbkpt is None:
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

        # Get the breakpoints
        self.breakpoints = bspline.get_breakpoints(
            x=x, bkpt=bkpt, fullbkpt=fullbkpt, nord=nord, bkspread=bkspread, bkspace=bkspace,
            nbkpts=nbkpts, everyn=everyn
        )

        # Finalize the setup
        nc = self.breakpoints.size - nord
        self.nord = nord
        self.npoly = npoly
        self.mask = np.ones((self.breakpoints.size,), dtype=bool)
        if npoly > 1:
            self.coeff = np.zeros((npoly, nc), dtype=float)
            self.icoeff = np.zeros((npoly, nc), dtype=float)
        else:
            self.coeff = np.zeros((nc,), dtype=float)
            self.icoeff = np.zeros((nc,), dtype=float)
        self.xmin = 0.0
        self.xmax = 1.0
        self.funcname = funcname

    @staticmethod
    def _fill_bkpt(bkpt, nord, bkspread):
        """
        Helper function used to pad the breakpoint vector according to the order
        of the b-spline.

        Parameters
        ----------
        bkpt : `numpy.ndarray`_
            The current set of breakpoints.
        nord : int
            Order of the b-spline.
        bkspread : float
            Scale factor for the separation between breakpoints.

        Returns
        -------
        `numpy.ndarray`_
            The padded set of breakpoints (typically ``fullbkpt`` as used by the class).
        """
        bkspace = (bkpt[1] - bkpt[0])*bkspread
        indx = np.arange(1, nord)
        return np.concatenate([bkpt[0] - bkspace*indx[::-1], bkpt, bkpt[-1] + bkspace*indx])

    @staticmethod
    def get_breakpoints(x=None, bkpt=None, fullbkpt=None, nord=4, bkspread=1.0, bkspace=None,
                        nbkpts=None, everyn=None):
        """
        Generate the set of breakpoints for the b-spline.

        Parameters
        ----------
        x : `numpy.ndarray`_, optional
            Independent variable for the definition of the b-spline.  If None,
            ``fullbkpt`` must be provided.
        bkpt : `numpy.ndarray`_, optional
            A precalculated set of breakpoints within the range of ``x``.  The
            input vector is sorted and any points beyond the range of x are
            omitted.  If only one breakpoint is provided (or not omitted), one
            breakpoint is set at each end of ``x``.  If the breakpoints do not
            cover the full range, the first and last breakpoints are moved such
            that they do.  If None, the breakpoints are determined using the
            keywords below.
        fullbkpt : `numpy.ndarray`_, optional
            The full set of breakpoints.  The input vector is sorted and cast as
            a float.  If the length of the vector is less than twice ``nord``,
            it is also padded with ``nord-1`` values, as needed for the
            construction of the b-spline.  If None, ``x`` must be provided.
        nord : int, optional
            Order of the b-spline.
        bkspread : float, optional
            Scale factor for the separation between breakpoints.
        bkspace : float, optional
            Defines the separation between breakpoints.  If provided, ``nbkpts``
            and ``everyn`` are ignored.
        nbkpts : int, optional
            Defines the number of breakpoints used to span the full range of
            ``x``.  Only used if ``bkspace`` is None.  If provided, ``everyn``
            is ignored.
        everyn : int, float, optional
            Places a breakpoint at every Nth value of ``x``.  Only used if
            ``bkspace`` and ``nbkpts`` are both None.

        Returns
        -------
        `numpy.ndarray`_
            Vector with the breakpoints

        Raises
        ------
        ValueError
            Raised if neither ``fullbkpt`` nor ``x`` are provided.
        """
        if fullbkpt is not None:
            _fullbkpt = np.sort(fullbkpt, kind='heapsort').astype(float)
            # JFH added this to fix bug in cases where fullbkpt is passed in but has
            # < 2*nord elements
            if _fullbkpt.size < 2*nord:
                _fullbkpt = bspline._fill_bkpt(_fullbkpt, nord, bkspread)
            return _fullbkpt

        if x is None:
            raise ValueError('Must provide `x` to determine breakpoints')

        sx = np.amin(x)
        ex = np.amax(x)
        if bkpt is None:
            if bkspace is not None:
                if bkspace >= ex - sx:
                    _bkpt = np.array([sx, ex])
                else:
                    _nbkpts = int((ex-sx)/bkspace) + 1
                    _bkpt = np.linspace(sx, ex, _nbkpts)
            elif nbkpts is not None:
                _bkpt = np.linspace(sx, ex, max(nbkpts,2))
            elif everyn is not None:
                # NOTE: There are places in the code where everyn is a float.
                # Need to continue to allow this.
                if everyn < x.size:
                    _nbkpts = max(x.size/everyn, 2.)
                    indx = (x.size/_nbkpts) * np.arange(_nbkpts)
                    _bkpt = np.interp(indx, np.arange(x.size, dtype=float), x)
                else:
                    _bkpt = np.array([sx, ex])
            else:
                raise ValueError('Insufficient information to set bkpts.')
        else:
            _bkpt = np.sort(bkpt, kind='heapsort')
            w = (_bkpt >= sx) & (_bkpt <= ex)
            _bkpt = np.array([sx, ex]) if np.sum(w) < 2 else _bkpt[w]

        # JFH added this new code, because bkpt.size = 1 implies fullbkpt
        # has only 2*(nord-1) + 1 elements.  This will cause a crash in
        # action because nbkpt < 2*nord, i.e. for bkpt = 1, nord = 4
        # fullbkpt has seven elements which is less than 2*nord = 8. The
        # codes above seem to require nbkpt >=2, so I'm implementing this
        # requirement. Note that the previous code before this fix simply
        # sets bkpt to bkpt[imax] =x.max() which is equally arbitrary, but
        # still results in a crash. By requiring at least 2 bkpt, fullbkpt
        # will have 8 elements preventing action from crashing.

        if _bkpt.size < 2:
            _bkpt = np.array([sx, ex])
        if _bkpt[0] > sx:
            _bkpt[0] = sx
        if _bkpt[-1] < ex:
            _bkpt[-1] = ex

        return bspline._fill_bkpt(_bkpt, nord, bkspread).astype(float)
    
    def reinit_coeff(self):
        nc = self.breakpoints.size - self.nord
        self.coeff = np.zeros((self.npoly, nc), dtype=float) if self.npoly > 1 \
                        else np.zeros(nc, dtype=float)

    def _bundle(self):
        """
        Overload the base class method (see
        :func:`~pypeit.datamodel.DataContainer._bundle`) to set the HDU name
        explicitly to BSPLINE.
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

    # TODO: C this
    # TODO: Should this be used, or should we effectively replace it
    # with the content of utils.bspline_profile
    def fit(self, xdata, ydata, invvar, x2=None):
        """
        Calculate a B-spline in the least-squares sense.

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
        xsort = x.argsort(kind='stable')
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
            return yfit[np.argsort(xsort, kind='stable')], mask

        for jj in range(hmm.size):
            mask[(x >= self.breakpoints[goodbk[hmm[jj]]])
                    & (x <= self.breakpoints[goodbk[hmm[jj]+1]-1])] = False
        return yfit[np.argsort(xsort, kind='stable')], mask

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
    result : `numpy.ndarray`_
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


