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
    Construct a two-dimensional affine transformation matrix in homogeneous
    coordinates.

    The coordinate convention is to assume the coordinates are ordered in an
    array with the Cartesian :math:`x` in the first column and Cartesian
    :math:`y` in the second column.  See the examples below.
    
    This function currently does *not* allow for a shear term.  The
    transformation is returned as a `numpy.ndarray`_ in homogeneous coordinates.
    
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
            Counter-clockwise rotation in radians.  The ordinate is aligned with
            the first axis, and the abcissa is aligned with the second axis; see
            :func:`~pypeit.core.transform.coordinate_transform_2d`.  If None,
            the rotation is 0.
        translation (array-like, optional):
            An array with two elements that provide the shift to apply in each
            dimension.

    Returns:
        `numpy.ndarray`_: A :math:`3\times3` affine-transformation matrix.

    Examples:

        Rotate the unit box by 45 degrees:

        >>> import numpy as np
        >>> from pypeit.core import transform
        >>> # Cartesian:     x   y
        >>> coo = np.array([[0., 0.],
                            [1., 0.],
                            [1., 1.],
                            [0., 1.]])
        >>> tform = transform.affine_transform_matrix(rotation=np.radians(45.))
        >>> tform
        array([[ 0.70710678, -0.70710678,  0.        ],
               [ 0.70710678,  0.70710678,  0.        ],
               [ 0.        ,  0.        ,  1.        ]])
        >>> transform.coordinate_transform_2d(coo, tform)
        array([[ 0.00000000e+00,  0.00000000e+00],
               [ 7.07106781e-01,  7.07106781e-01],
               [ 1.11022302e-16,  1.41421356e+00],
               [-7.07106781e-01,  7.07106781e-01]])

        Rotation about the center by first shifting, then rotating, then shifting back.

        >>> shift = affine_transform_matrix(translation=[-0.5,-0.5])
        >>> _tform = np.linalg.inv(shift) @ tform @ shift
        >>> transform.coordinate_transform_2d(coo, _tform)
        array([[ 0.5       , -0.20710678],
               [ 1.20710678,  0.5       ],
               [ 0.5       ,  1.20710678],
               [-0.20710678,  0.5       ]])

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
    Construct an affine transform matrix from a set of transformation steps
    executed in series.

    Each step in the transformation is provided by a call to
    :func:`pypeit.core.transform.affine_transform_matrix`.  The order of the
    steps should be provided in the order they should be preformed.  The
    following should pass:

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


def coordinate_transform_2d(coo, matrix, inverse=False):
    r"""
    Apply a 2D coordinate transformation using an affine-transformation matrix
    in homologous coordinates.

    Args:
        coo (array-like):
            Coordinates to transform.  Shape must be :math:`(N,2)`, where
            :math:`N` is the number of coordinates, Cartesian :math:`x` is in
            the first column (``coo[:,0]``), and Cartesian :math:`y` is in the
            second column (``coo[:,1]``).
        matrix (`numpy.ndarray`_):
            The :math:3\times 3` affine-transformation matrix.  See
            :func:`~pypeit.core.mosaic.affine_transform_matrix`.
        inverse (:obj:`bool`, optional):
            By default, the function performs the *active* transformation; i.e.,
            applying the transformation to the coordinates, moving them within
            the existing coordinate frame.  Set ``inverse`` to true to instead
            perform the *passive* transformation; i.e., transforming the
            coordinate axes and providing the new coordinates in the transformed
            (e.g., rotated) frame.  
    
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



