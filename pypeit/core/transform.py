"""
Provide basic coordinate tranformation functions.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""

from IPython import embed

import numpy as np

from pypeit import msgs


def affine_transform_matrix(scale=None, rotation=None, translation=None):
    r"""
    Construct an affine transformation matrix for a two-dimensional image in
    homologous coordinates.

    Method currently does *not* allow for a shear term.  The transformation is
    returned as a `numpy.ndarray`_ in homologous coordinates.
    
    When applied to an image, this follows the convention that coordinates are
    ordered as Cartesian :math:`x` along the first axis and Cartesian :math:`y`
    along the second axis.

    If no arguments are provided, this simply returns an identity matrix.

    Otherwise, the order of operations is to scale, rotate, and then translate
    (shift).  I.e., the following should pass:

    .. code-block:: python

        t = affine_transform_matrix(translation=[2,1])
        r = affine_transform_matrix(rotation=np.pi/3)
        a1 = affine_transform_matrix(rotation=np.pi/3, translation=[2,1])
        assert np.array_equal(t @ r, a1)

        s = affine_transform_matrix(scale=(1.5, 2.))
        a2 = affine_transform_matrix(scale=(1.5, 2.), rotation=np.pi/3, translation=[2,1])
        assert np.array_equal(t @ r @ s, a2)

    Args:
        scale (:obj:`float`, array-like, optional):
            Scale factor to apply to each axis.  Can be a single float to apply
            to both axes, or a two-element array-like with the scaling for each
            axis.
        rotation (:obj:`float`, optional):
            Counter-clockwise rotation in radians.  If None, rotation is 0.
        translation (array-like, optional):
            An array with two elements that provide the shift to apply in both
            image dimensions.

    Returns:
        `numpy.ndarray`_: A :math:`3\times3` affine-transformation matrix.
    """
    tform = np.eye(3)
    if scale is not None:
        _s = np.atleast_1d(scale)
        if _s.size == 1:
            sx = sy = scale
        elif _s.size == 2:
            sx, sy = scale
        else:
            msgs.error('Scale must be a single scalar or a two-element array.')
        tform[0,0] = float(sx)
        tform[1,1] = float(sy)
    if rotation is not None:
        tform[:2,:2] = np.array([[np.cos(rotation), -np.sin(rotation)],
                                 [np.sin(rotation), np.cos(rotation)]]) @ tform[:2,:2]
    if translation is not None:
        _t = np.atleast_1d(translation)
        if _t.size != 2:
            msgs.error('Translation must be a two-element array.')
        tform[0:2,2] = translation
    return tform


def affine_transform_series(steps):
    r"""
    Construct an affine transform from a set of transformation steps executed in
    series.

    Each step in the transformation is provided by a call to
    :func:`pypeit.core.transform.affine_transform_matrix`.  The order of the
    steps should provided in the order they should be preformed.  The following
    should pass:

    .. code-block:: python

        steps = [dict(scale=(1.5,2)), dict(rotation=np.pi/3), dict(translation=[2,1])]
        a1 = affine_transform_series(steps)
        a2 = affine_transform_matrix(scale=steps[0]['scale'], rotation=steps[1]['rotation'],
                                     translation=steps[2]['translation'])
        assert np.array_equal(a1, a2)

    Args:
        steps (array-like):
            A list of dictionaries with each dictionary containing the keyword
            arguments for a single call to 
            :func:`pypeit.core.transform.affine_transform_matrix`.
    
    Returns:
        `numpy.ndarray`_: A :math:`3\times3` affine-transformation matrix.
    """
    tform = np.eye(3)
    for kwargs in steps:
        tform = affine_transform_matrix(**kwargs) @ tform
    return tform


def coordinate_transform_2d(coo, matrix, inverse=True):
    r"""
    Apply a 2D coordinate transformation using an affine-transformation matrix
    in homologous coordinates.

    Args:
        coo (array-like):
            Coordinates to transform.  Shape must be :math:`(N,2)`, where
            :math:`N` is the number of coordinates.
        matrix (`numpy.ndarray`_):
            The affine-transformation matrix.  See
            :func:`pypeit.core.mosaic.affine_transform_matrix`.
        inverse (:obj:`bool`, optional):
            By default, the coordinate transformation 
    
    Returns:
        `numpy.ndarray`_: Transformed coordinates with shape :math:`(N,2)`.
    """
    _coo = np.atleast_2d(coo)
    if _coo.ndim != 2:
        msgs.error('Coordinate array must be 2D.')
    if _coo.shape[1] != 2:
        msgs.error('Coordinate array must have 2D coordinates along the last axis.')
    ncoo = _coo.shape[0]
    _m = np.linalg.inv(matrix) if inverse else matrix
    return (np.column_stack((_coo, np.ones(ncoo, dtype=_coo.dtype))) @ _m.T)[:,:2]



