"""
Module to run tests on affine transform
"""
from IPython import embed
import numpy as np

from pypeit.core import transform

def test_affine_transform_order():
    """
    Construct affine transforms with each of the three possible matrix elements
    --- scale, rotation, translation --- and check that the order is as
    expected.
    """

    # Translation matrix
    t = transform.affine_transform_matrix(translation=[2,1])
    # Rotation matrix
    r = transform.affine_transform_matrix(rotation=np.pi/3)
    # Composite
    a1 = transform.affine_transform_matrix(rotation=np.pi/3, translation=[2,1])
    # Order should be rotation then translation
    assert np.array_equal(t @ r, a1), 'Bad rotation+translation order'

    # Scale matrix
    s = transform.affine_transform_matrix(scale=(1.5, 2.))
    # Composite
    a2 = transform.affine_transform_matrix(scale=(1.5, 2.), rotation=np.pi/3, translation=[2,1])
    # Order should be scale, rotate, translate
    assert np.array_equal(t @ r @ s, a2), 'Bad scale+rotation+translation order'


def test_transform_series():
    """
    Test construction of composite from a series.
    """
    steps = [dict(scale=2), dict(rotation=np.pi/3), dict(translation=[2,1])]
    a1 = transform.affine_transform_series(steps)
    a2 = transform.affine_transform_matrix(scale=steps[0]['scale'], rotation=steps[1]['rotation'],
                                           translation=steps[2]['translation'])
    assert np.array_equal(a1, a2), 'Bad series order'


def test_coordinate_transform():
    tform = transform.affine_transform_matrix(rotation=np.pi/4)
    coo = np.array([[1, 1], [1, -1], [-1, -1], [-1, 1]])
    _coo = transform.coordinate_transform_2d(coo, tform, inverse=True)
    sqrt2 = np.sqrt(2.)
    assert np.allclose(np.array([[sqrt2, 0.], [0., -sqrt2], [-sqrt2, 0.], [0., sqrt2]]), _coo), \
            'Bad rotation'


