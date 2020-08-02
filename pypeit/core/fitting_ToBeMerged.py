""" Module for fitting codes

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst

"""

import numpy as np


def twoD_Gaussian(tup, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    """ A 2D Gaussian to be used to fit the cross-correlation

    Args:
        tup (tuple):
            A two element tuple containing the (x,y) coordinates where the 2D Gaussian will be evaluated

    Returns:
        model (`numpy.ndarray`_)
    """
    (x, y) = tup
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))
    return g.ravel()

