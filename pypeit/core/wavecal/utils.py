""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypeit.core import arc
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit


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


def match_peaks(inspec1, inspec2, smooth=5.0, debug=False):
    """ Stretch and shift inspec2 until it matches inspec1
    """

    # Initial estimate
    p0 = np.array([0.0])
    specs = (inspec1, inspec2, smooth,)
    res = curve_fit(shift_stretch, specs, np.array([0.0]), p0, epsfcn=1.0)
    bestpar = res[0][0]

    inspec2_adj = resample(inspec2, int(inspec1.size + bestpar))
    _, shift = shift_stretch(specs, bestpar, retshift=True)
    xval1 = np.arange(inspec2_adj.shape[0]) + shift
    xval2 = np.linspace(0.0, inspec2.size-1, inspec2_adj.shape[0])

    if debug:
        x1 = np.arange(inspec1.shape[0])
        x2 = np.arange(inspec2_adj.shape[0]) + shift
        from matplotlib import pyplot as plt
        plt.plot(x1, inspec1, 'k-', drawstyle='steps')
        plt.plot(x2, inspec2_adj, 'r-', drawstyle='steps')
        plt.show()

    return inspec2_adj, xval1, xval2, shift


def shift_stretch(specs, p, retshift=False):
    inspec1, inspec2, smooth = specs
    y1 = gaussian_filter(inspec1, smooth)
    y2 = gaussian_filter(inspec2, smooth)
    y1size = y1.size
    y2size = int(y1size + p)
    y2 = resample(y2, y2size)
    df = np.min([y1size // 2 - 1, y2size // 2 - 1])
    size = y1size + y2size - 1
    fsize = 2 ** np.int(np.ceil(np.log2(size)))  # Use this size for a more efficient computation
    conv = np.fft.fft(y1, fsize)
    conv *= scipy.conj(np.fft.fft(y2, fsize))
    cc = scipy.ifft(conv)  # [df:df+y1size]
    shift = np.argmax(np.abs(cc))
    stretch = 1.0 / np.max(np.abs(cc))
    if retshift:
        return np.array([stretch]), shift
    else:
        return np.array([stretch])
