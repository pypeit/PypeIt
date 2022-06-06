

import os
import sys
import numpy as np
from pypeit.core import fitting
from astropy import constants as const
c_kms = const.c.to('km/s').value



def predict_order_coverage(arxiv_params, arxiv, xdisp, xdangle, norders, pad=0):
    """
    Predict the coverage of orders in the echelle spectrum using the disperser dependent
    fits of the reddest order as a function of xdangle.

    Args:
        arxiv_params (astropy.table.Table):
            Table holding the arxiv parameters
        arxiv (astropy.table.Table): \
            Table holding the arxiv data
        xdisp (str):
            Corss disperser. For HIRES this is either 'UV' or 'RED'
        xdangle (float):\
            Cross-disperser angle.
        norders (int):
            Number of orders identified on the detector
        pad (int):
            Number of orders to pad the coverage by on the blue and red side.

    Returns:
        order_vec (numpy.ndarray):
            Array of order numbers for the predicted coverage.
    """

    xd_min, xd_max = arxiv_params['xd_min'][0], arxiv_params['xd_max'][0]
    idisp = arxiv_params['xdisp_vec'] == xdisp
    reddest_order_fit = int(np.round(fitting.evaluate_fit(
        arxiv['xd_coeffs'][idisp, :], arxiv_params['xd_func'], xdangle, minx=xd_min, maxx=xd_max)))
    order_vec = reddest_order_fit + (np.arange(norders + 2*pad) - pad)[::-1]

    return order_vec
