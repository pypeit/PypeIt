""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import numba as nb
from pypeit.core import arc

from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit


def arc_lines_from_spec(spec, min_nsig =10.0, nonlinear_counts = 1e10):
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
    tampl, tcent, twid, centerr, w, yprep, nsig = arc.detect_lines(spec, nfitpix=7, nonlinear_counts = nonlinear_counts)
    all_tcent = tcent[w]
    all_tampl = tampl[w]
    all_ecent = centerr[w]
    all_nsig = nsig[w]

    # Old code cut on amplitude
    # Cut on Amplitude
    #cut_amp = all_tampl > min_ampl

    # Cut on significance
    cut_sig = all_nsig > min_nsig
    cut_tcent = all_tcent[cut_sig]
    icut = np.where(cut_sig)[0]
     # Return
    return all_tcent, all_ecent, cut_tcent, icut

# JFH ToDo This algorithm for computing the shift and stretch is unstable and often hangs. It needs to be replaced.
def match_peaks(inspec1, inspec2, smooth=5.0, debug=False):
    """ Stretch and shift inspec2 until it matches inspec1
    """

    # Initial estimate
    p0 = np.array([0.0])
    nspec = inspec1.size
    specs = (inspec1, inspec2, smooth,)

    try:
        res = curve_fit(shift_stretch, specs, np.array([0.0]), p0,  bounds = (-nspec, nspec))
        #res = curve_fit(shift_stretch, specs, np.array([0.0]), p0, epsfcn=1.0)
    except ValueError:
        # Probably no overlap of the two spectra
        return None, None
    stretch = res[0][0]
    _, shift = shift_stretch(specs, stretch, retshift=True)

    if debug:
        inspec2_adj = resample(inspec2, int(inspec1.size + stretch))
        x1 = np.arange(inspec1.shape[0])
        x2 = np.arange(inspec2_adj.shape[0]) + shift
        from matplotlib import pyplot as plt
        plt.plot(x1, inspec1, 'k-', drawstyle='steps')
        plt.plot(x2, inspec2_adj, 'r-', drawstyle='steps')
        plt.show()

    return stretch, shift


# JFH I think this should be done with scipy.optimize to find the maximum value of the cc correlation as a function of
# shift and stretch, rather than with curve_fit
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

# JFH These codes are under development by JFH to try to improve the problems with match_peaks.
def match_peaks2(inspec1, inspec2, smooth=5.0, debug=False):

    # Initial estimate
    p0 = np.array([0.0])
    nspec = inspec1.size
    # differential evolution optimizer appears to behave better
    print('Optimizing with IGM')
    guess = np.array([0.0,0.0])
    bounds = [(-nspec, nspec), (-nspec, nspec)]
    result = op.differential_evolution(xcorr_stretch, args=(inspec1, inspec2,smooth), bounds=bounds, popsize=25, recombination=0.7,
                                         disp=True, polish=True) # ToDO Should we polish?




# JFH I think this should be done with scipy.optimize to find the maximum value of the cc correlation as a function of
# shift and stretch, rather than with curve_fit
def xcorr_stretch(theta, inspec1, inspec2, smooth):

    y1 = gaussian_filter(inspec1, smooth)
    y2 = gaussian_filter(inspec2, smooth)
    y1size = y1.size
    y2size = int(y1size + theta[0])
    y2 = resample(y2, y2size)
    df = np.min([y1size // 2 - 1, y2size // 2 - 1])
    size = y1size + y2size - 1
    fsize = 2 ** np.int(np.ceil(np.log2(size)))  # Use this size for a more efficient computation
    conv = np.fft.fft(y1, fsize)
    conv *= scipy.conj(np.fft.fft(y2, fsize))
    cc = scipy.ifft(conv)  # [df:df+y1size]
    shift = np.argmax(np.abs(cc))
    return np.max(np.abs(cc))


@nb.jit(nopython=True, cache=True)
def hist_wavedisp(waves, disps, dispbin=None, wavebin=None, scale=1.0, debug=False):
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

    # Convert to linear
    lin_dispbin = np.power(10.0, dispbin)
    lin_disps = np.power(10.0, disps)

    # Determine how many elements will be used for the histogram
    nelem = np.zeros(dispbin.size-1, dtype=nb.types.uint64)
    for dd in range(dispbin.size-1):
        dispval = 0.5*(lin_dispbin[dd] + lin_dispbin[dd+1])
        nelem[dd] = np.int(0.5 + scale*(wavebin[1]-wavebin[0])/dispval)

    # Generate a histogram
    nhistv = np.sum(nelem)
    hist_wd = np.zeros(nhistv, dtype=nb.types.uint64)
    cent_w = np.zeros(nhistv, dtype=nb.types.uint64)
    cent_d = np.zeros(nhistv, dtype=nb.types.uint64)
    cntr = 0
    for dd in range(dispbin.size-1):
        wbin = np.linspace(wavebin[0], wavebin[1], nelem[dd]+1)
        wdsp = np.where((lin_disps > lin_dispbin[dd]) & (lin_disps <= lin_dispbin[dd+1]))
        if wdsp[0].size != 0:
            hist_wd[cntr:cntr+nelem[dd]], _ = np.histogram(waves[wdsp], bins=wbin)
            cent_d[cntr:cntr+nelem[dd]] = 0.5 * (lin_dispbin[dd] + lin_dispbin[dd + 1])
            cent_w[cntr:cntr+nelem[dd]] = 0.5 * (wbin[1:] + wbin[:-1])
        cntr += nelem[dd]

    """
    if debug:
        # Create and plot up the 2D plot
        nelem_mx = np.max(nelem)
        hist_wd_plt = np.zeros((nelem_mx, dispbin.size), dtype=nb.types.uint64)
        wbin_mx = np.linspace(wavebin[0], wavebin[1], nelem_mx)
        cntr = 0
        for dd in range(dispbin.size-1):
            wbin = np.linspace(wavebin[0], wavebin[1], nelem[dd]+1)
            wdsp = np.where((lin_disps > lin_dispbin[dd]) & (lin_disps <= lin_dispbin[dd+1]))
            if wdsp[0].size != 0:
                fval, _ = np.histogram(waves[wdsp], bins=wbin)
                fspl = interpolate.interp1d(cent_w[cntr:cntr+nelem[dd]], fval, bounds_error=False, fill_value="extrapolate")
                val = fspl(wbin_mx).astype(nb.types.uint64)
                hist_wd_plt[:, dd] = val
            cntr += nelem[dd]
        plt.clf()
        plt.imshow(np.log10(np.abs(hist_wd_plt[:, ::-1].T)), extent=[wavebin[0], wavebin[1], dispbin[0], dispbin[-1]], aspect='auto')
        plt.show()
        pdb.set_trace()
    """

    return hist_wd, cent_w, cent_d
