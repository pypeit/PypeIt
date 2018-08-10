""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np


def arc_lines_from_spec(spec, min_ampl=300.):
    """
    Parameters
    ----------
    spec
    siglev
    min_ampl

    Returns
    -------

    """

    # imports
    from arclines.pypit_utils import find_peaks
    # Find peaks
    tampl, tcent, twid, w, yprep = find_peaks(spec)
    all_tcent = tcent[w]
    all_tampl = tampl[w]

    # Cut on Amplitude
    cut_amp = all_tampl > min_ampl
    cut_tcent = all_tcent[cut_amp]
    icut = np.where(cut_amp)[0]

    # Return
    return all_tcent, cut_tcent, icut
