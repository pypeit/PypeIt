""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypeit.core import arc
import numba as nb


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


@nb.jit(nopython=True, cache=True)
def hist_wavedisp(waves, disps, dispbin=None, wavebin=None, scale=1.0, debug=True):
    """ Generates a flexible 2D histogram of central wavelength and
    dispersion, where the wavelength grid spacing depends on the
    dispersion, so that each wavelength bin width roughly corresponds
    to the sample dispersion, scaled by a factor.

    Parameters
    ----------
    waves : ndarray
      array of central wavelengths (A). Must be the same size as disps.
    disps : ndarray
      logarithmic array of dispersions (A/pix). Must be the same size as waves.
    dispbin : ndarray
      bin widths for the dispersion dimension
    wavebin : two-element list
      minimum and maximum wavelength of the histogram bins
    scale : float
      Scale the sampling of the wavelength bin (see description above). Must be >= 1.0

    Returns
    -------
    hist_wd : ndarray
      An array of bin counts
    cent_w : ndarray
      The value of the wavelength at the centre of each bin
    cent_d : ndarray
      The value of the dispersion at the centre of each bin
    """
    if dispbin is None:
        dispbin = np.linspace(-3.0, 1.0, 1000)
    if wavebin is None:
        wavebin = [np.min(waves), np.max(waves)]

    lin_dispbin = 10.0**dispbin

    # Determine how many elements will be used for the histogram
    nelem = np.zeros(dispbin.size-1, dtype=nb.types.uint64)
    for dd in range(dispbin.size-1):
        dispval = 0.5*(lin_dispbin[dd] + lin_dispbin[dd+1])
        nelem[dd] = np.int(0.5 + scale*(wavebin[1]-wavebin[2])/dispval)

    # Generate a histogram
    nhistv = np.sum(nelem[dd])
    hist_wd = np.zeros(nhistv, dtype=nb.types.uint64)
    cent_w = np.zeros(nhistv, dtype=nb.types.uint64)
    cent_d = np.zeros(nhistv, dtype=nb.types.uint64)
    cntr = 0
    for dd in range(dispbin.size-1):
        wbin = np.linspace(wavebin[0], wavebin[1], nelem[dd])
        wdsp = np.where((disps > lin_dispbin[dd]) & (disps <= lin_dispbin[dd+1]))
        hist_wd[cntr:cntr+1+nelem[dd]], _ = np.histogram(waves[wdsp], bins=wbin)
        cent_d[cntr:cntr+1+nelem[dd]] = 0.5*(lin_dispbin[dd] + lin_dispbin[dd+1])
        cent_w[cntr:cntr+1+nelem[dd]] = wbin

    if debug:
        from matplotlib import pyplot as plt
        import pdb
        # Create and plot up the 2D plot
        nelem_mx = np.max(nelem)
        hist_wd_plt = np.zeros((nelem_mx, dispbin.size), dtype=nb.types.uint64)
        wbin_mx = np.linspace(wavebin[0], wavebin[1], nelem_mx)
        cntr = 0
        for dd in range(dispbin.size-1):
            wbin = np.linspace(wavebin[0], wavebin[1], nelem[dd])
            wdsp = np.where((disps > lin_dispbin[dd]) & (disps <= lin_dispbin[dd+1]))
            fval, _ = np.histogram(waves[wdsp], bins=wbin)
            hist_wd_plt[:, dd], _ = np.interp(wbin_mx, wbin, fval)
        plt.clf()
        plt.imshow(np.log10(np.abs(hist_wd_plt[:, ::-1].T)), extent=[wavebin[0], wavebin[1], dispbin[0], dispbin[-1]], aspect='auto')
        plt.show()
        pdb.set_trace()

    return hist_wd, cent_w, cent_d