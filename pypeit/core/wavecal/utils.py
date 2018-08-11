""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypeit.core import arc


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

    # Find peaks
    tampl, tcent, twid, w, yprep = arc.detect_lines(spec, nfitpix=7)
    all_tcent = tcent[w]
    all_tampl = tampl[w]

    # Cut on Amplitude
    cut_amp = all_tampl > min_ampl
    cut_tcent = all_tcent[cut_amp]
    icut = np.where(cut_amp)[0]

    # Return
    return all_tcent, cut_tcent, icut
