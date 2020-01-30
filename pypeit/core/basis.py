# Licensed under a 3-clause BSD style license - see PYDL_LICENSE.rst
# -*- coding: utf-8 -*-
# Also cite https://doi.org/10.5281/zenodo.1095150 when referencing PYDL
"""
Module implements a set of basis functions for fitting.
"""
from IPython import embed

import numpy as np
from scipy import special

def _init_basis(x, m):
    r"""
    Initialize the basis functions.

    Args:
        x (array-like):
            Compute the basis polynomials at these abscissa values.
        m (:obj:`int`):
            The number of polynomials to compute. For example, if
            :math:`m = 3`, :math:`P_0 (x)`, :math:`P_1 (x)` and
            :math:`P_2 (x)` will be computed. Must be :math:`\geq1`.

    Returns:
        tuple: Returns (1) the input :math:`x` array explicitly
        converted to a `numpy.ndarray`_ and (2) a unity array with
        the same data type as :math:`x`. The returned shape of the
        latter is the :math:`(N_x, m)`.

    Raises:
        ValueError:
            Raised if the input order is not at least 1.
    """
    if m < 1:
        raise ValueError('Order must be at least 1.')
    _x = np.atleast_1d(x)
    return _x, np.ones((_x.size, m), dtype=_x.dtype)


def _build_basis(x, m, func):
    r"""
    Perform initial checks of the basis function inputs.

    Args:
        x (array-like):
            Compute the basis polynomials at these abscissa values.
        m (:obj:`int`):
            The number of polynomials to compute. For example, if
            :math:`m = 3`, :math:`P_0 (x)`, :math:`P_1 (x)` and
            :math:`P_2 (x)` will be computed. Must be :math:`\geq1`.
        func (callable):
            Callable function that generates the basis polynomials.
            E.g., `scipy.special.legendre` for Legendre polynomials.

    Returns:
        `numpy.ndarray`_: An array of shape :math:`(N_x, m)` with
        the basis polynomials.

    Raises:
        TypeError:
            Raised if the provided ``func`` is not callable.
    """
    if not callable(func):
        raise TypeError('Must provide a callable function that constructs the basis polynomials.')
    _x, basis = _init_basis(x, m)
    if m >= 2:
        basis[:,1] = _x
    if m >= 3:
        for k in range(2, m):
            # TODO: Is there a faster way to set this up?
            basis[:,k] = np.polyval(func(k), _x)
    return basis


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
    return _build_basis(x, m, special.legendre)


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
    return _build_basis(x, m, special.chebyt)


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
    _x, basis = _init_basis(x, m)
    basis[:,0] = (_x >= 0).astype(_x.dtype)
    if m > 2:
        basis[:,2] = _x
    if m > 3:
        for k in range(3, m):
            basis[:,k] = 2.0 * _x * basis[:,k-1] - basis[:,k-2]
    return basis


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
    _x, basis = _init_basis(x, m)
    if m >= 2:
        basis[:,1] = _x
    if m >= 3:
        for k in range(2, m):
            basis[:,k] = basis[:,k-1] * _x
    return basis

