# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
from IPython import embed

import numpy as np

from pypeit import msgs
from pypeit import utils
from pypeit.core import basis
from pypeit.core import fitting

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
        ia = np.ones((ncoeff,), dtype=bool)
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


# TODO -- This class needs to become a DataContainer
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
    pypeitFits : list
        Holds a list of :class:`pypeit.fitting.PypeItFit` fits
    """
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
                inmask = np.ones(xpos.shape, dtype=bool)
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
            self.outmask = np.zeros(xpos.shape, dtype=bool)
            self.yfit = np.zeros(xpos.shape, dtype=xpos.dtype)
            self.pypeitFits = []
            for iTrace in range(self.nTrace):
                xvec = self.xnorm(xpos[iTrace, :], do_jump)
                if invvar is None:
                    thisinvvar = None
                else:
                    thisinvvar = invvar[iTrace, :]

                pypeitFit = fitting.robust_fit(xvec, ypos[iTrace, :], self.ncoeff,
                                                                function=self.func, maxiter = self.maxiter,
                                                                in_gpm = inmask[iTrace, :], invvar = thisinvvar,
                                                                lower = self.lower, upper = self.upper,
                                                                minx = self.xmin, maxx = self.xmax,
                                                                maxdev=self.maxdev,
                                                                grow=0,use_mad=False,sticky=False)
                ycurfit_djs = pypeitFit.eval(xvec)
                self.pypeitFits.append(pypeitFit)

                # Load
                self.yfit[iTrace, :] = ycurfit_djs
                self.coeff[iTrace, :] = pypeitFit.fitc
                self.outmask[iTrace, :] = pypeitFit.gpm

        else:
            msgs.error('Wrong number of arguments to TraceSet!')

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
            #ypos[iTrace, :] =  utils.func_val(self.coeff[iTrace, :], xvec, self.func, minx=self.xmin, maxx=self.xmax)
            ypos[iTrace, :] =  self.pypeitFits[iTrace].eval(xvec)#, self.func, minx=self.xmin, maxx=self.xmax)
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


