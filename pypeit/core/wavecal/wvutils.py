""" Module for basic utilties with holy grail
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import numba as nb

from matplotlib import pyplot as plt
from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit
from pypeit import msgs


from pypeit.core import arc
from pypeit import debugger

def arc_lines_from_spec(spec, sigdetect=10.0, fwhm=4.0,fit_frac_fwhm = 1.25, mask_frac_fwhm=1.0,max_frac_fwhm=2.0,
                        cont_samp=30, niter_cont=3,nonlinear_counts=1e10, debug=False):
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
    tampl, tampl_cont, tcent, twid, centerr, w, yprep, nsig = arc.detect_lines(spec, sigdetect = sigdetect, fwhm=fwhm,
                                                                               fit_frac_fwhm=fit_frac_fwhm,
                                                                               mask_frac_fwhm=mask_frac_fwhm,
                                                                               max_frac_fwhm=max_frac_fwhm,
                                                                               cont_samp=cont_samp,niter_cont = niter_cont,
                                                                               nonlinear_counts = nonlinear_counts, debug=debug)
    all_tcent = tcent[w]
    all_tampl = tampl[w]
    all_ecent = centerr[w]
    all_nsig = nsig[w]

    # Old code cut on amplitude
    # Cut on Amplitude
    #cut_amp = all_tampl > min_ampl

    # Cut on significance
    cut_sig = all_nsig > sigdetect
    cut_tcent = all_tcent[cut_sig]
    icut = np.where(cut_sig)[0]

    #debugger.set_trace()
    # Return
    return all_tcent, all_ecent, cut_tcent, icut


def shift_and_stretch(spec, shift, stretch):

    """ Utility function to shift and stretch a spectrum. This operation is being implemented in many steps and
    could be significantly optimized. But it works for now. Note that the stretch is applied *first* and then the
    shift is applied in stretch coordinates.

    Parameters
    ----------
    spec : ndarray
      spectrum to be shited and stretch
    shift: float
      shift to be applied
    stretch: float
      stretch to be applied

    Returns
    -------
    spec_out: ndarray
      shifted and stretch spectrum. Regions where there is no information are set to zero.
    """

    # Positive value of shift means features shift to larger pixel values

    nspec = spec.shape[0]
    # pad the spectrum on both sizes
    x1 = np.arange(nspec)/float(nspec)
    nspec_stretch = int(nspec*stretch)
    x2 = np.arange(nspec_stretch)/float(nspec_stretch)
    spec_str = (scipy.interpolate.interp1d(x1, spec, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(x2)
    # Now create a shifted version
    ind_shift = np.arange(nspec_stretch) - shift
    spec_str_shf = (scipy.interpolate.interp1d(np.arange(nspec_stretch), spec_str, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(ind_shift)
    # Now interpolate onto the original grid
    spec_out = (scipy.interpolate.interp1d(np.arange(nspec_stretch), spec_str_shf, kind = 'quadratic', bounds_error = False, fill_value = 0.0))(np.arange(nspec))

    return spec_out


def zerolag_shift_stretch(theta, y1, y2):

    """ Utility function which is run by the differential evolution optimizer in scipy. These is the fucntion we optimize.
    It is the zero lag cross-correlation coefficient of spectrum with a shift and stretch applied.

     Parameters
     ----------
     theta: ndarray
       Function parameters to optmize over. theta[0] = shift, theta[1] = stretch
     y1: ndarray, shape = (nspec,)
       First spectrum which acts as the refrence
     y2: ndarray,  shape = (nspec,)
       Second spectrum which will be transformed by a shift and stretch to match y1

     Returns
     -------
     corr_norm: float
       Negative of the zero lag cross-correlation coefficient (since we are miniziming with scipy.optimize). scipy.optimize will
       thus determine the shift,stretch that maximize the cross-correlation.
     """


    shift, stretch = theta
    y2_corr = shift_and_stretch(y2, shift, stretch)
    # Zero lag correlation
    corr_zero = np.sum(y1*y2_corr)
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr_zero/corr_denom
    return -corr_norm

def smooth_ceil_cont(inspec1, smooth, percent_ceil = None, use_raw_arc=False):
    """ Utility routine to smooth and apply a ceiling to spectra """

    # ToDO can we improve the logic here. Technically if use_raw_arc = True and perecent_ceil=None
    # we don't need to peak find or continuum subtract, but this makes the code pretty uggly.

    # Run line detection to get the continuum subtracted arc
    tampl1, tampl1_cont, tcent1, twid1, centerr1, w1, arc1, nsig1 = arc.detect_lines(inspec1, sigdetect=10.0)
    if use_raw_arc == True:
        ampl = tampl1
        use_arc = inspec1
    else:
        ampl = tampl1_cont
        use_arc = arc1

    if percent_ceil is not None:
        # If this is set, set a ceiling on the greater > 10sigma peaks
        ceil1 = np.percentile(ampl, percent_ceil)
        spec1 = np.fmin(arc1, ceil1)
    else:
        spec1 = np.copy(arc1)

    if smooth is not None:
        y1 = scipy.ndimage.filters.gaussian_filter(spec1, smooth)
    else:
        y1 = np.copy(spec1)

    return y1

# This code does not currently work yet and has bugs in the lag computation.
#def cross_correlate(y1,y2):

#    if y1.shape != y2.shape:
#        msgs.error('cross_correlate only works for equal sized arrays')
#    nspec =y1.shape
#    next2 = 2**(nspec-1).bit_length()
#    f1 = np.fft.fft(y1,n=next2)
#    f2 = np.fft.fft(np.flipud(y2),n=next2)
#    cc_raw = np.real(np.fft.ifft(f1 * f2))
#    cc = np.fft.fftshift(cc_raw)
#    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
#    cc_norm = cc/corr_denom
#    zero_index = int(next2/2) - 1
#    lags = zero_index - np.arange(next2)
#    return lags, cc_norm

# ToDO can we speed this code up? I've heard numpy.correlate is faster. Someone should investigate optimization. Also we don't need to compute
# all these lags.
def xcorr_shift(inspec1,inspec2,smooth=5.0,percent_ceil=90.0, use_raw_arc=False, debug=False):

    """ Determine the shift inspec2 relative to inspec1.  This routine computes the shift by finding the maximum of the
    the cross-correlation coefficient. The convention for the shift is that positive shift means inspec2 is shifted to the right
    (higher pixel values) relative to inspec1.

    Parameters
    ----------
    inspec1 : ndarray
      Reference spectrum
    inspec2 : ndarray
      Spectrum for which the shift and stretch are computed such that it will match inspec1

    Optional Parameters
    -------------------
    smooth: float, default
      Gaussian smoothing in pixels applied to both spectra for the computations. Default is 5.0
    percent_ceil: float, default=90.0
      Appply a ceiling to the input spectra at the percent_ceil percentile level of the distribution of peak amplitudes.
      This prevents extremely strong lines from completely dominating the cross-correlation, which can
      causes the cross-correlation to have spurious noise spikes that are not the real maximum.
    use_raw_arc: bool, default = False
      If this parameter is True the raw arc will be used rather than the continuum subtracted arc
    debug: boolean, default = False

    Returns
    -------
    shift: float
      the shift which was determined
    cross_corr: float
      the maximum of the cross-correlation coefficient at this shift
    """

    y1 = smooth_ceil_cont(inspec1,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc)
    y2 = smooth_ceil_cont(inspec2,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc)


    nspec = y1.shape[0]
    lags = np.arange(-nspec + 1, nspec)
    corr = scipy.signal.correlate(y1, y2, mode='full')
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr/corr_denom
    output = arc.detect_lines(corr_norm, sigdetect=5.0, fit_frac_fwhm=1.25, fwhm=20.0, mask_frac_fwhm=1.0, cont_samp=30, nfind = 1)
    pix_max = output[1]
    corr_max = np.interp(pix_max, np.arange(lags.shape[0]),corr_norm)
    lag_max  = np.interp(pix_max, np.arange(lags.shape[0]),lags)
    if debug:
        # Interpolate for bad lines since the fitting code often returns nan
        plt.figure(figsize=(14, 6))
        plt.plot(lags, corr_norm, color='black', drawstyle = 'steps-mid', lw=3, label = 'x-corr', linewidth = 1.0)
        plt.plot(lag_max[0], corr_max[0],'g+', markersize =6.0, label = 'peak')
        plt.title('Best shift = {:5.3f}'.format(lag_max[0]) + ',  corr_max = {:5.3f}'.format(corr_max[0]))
        plt.legend()
        plt.show()

    return lag_max[0], corr_max[0]


def xcorr_shift_stretch(inspec1, inspec2, cc_thresh=-1.0, smooth=5.0, percent_ceil=90.0, use_raw_arc=False,
                        shift_mnmx=(-0.05,0.05), stretch_mnmx=(0.95,1.05), debug=False, seed = None):

    """ Determine the shift and stretch of inspec2 relative to inspec1.  This routine computes an initial
    guess for the shift via maximimizing the cross-correlation. It then performs a two parameter search for the shift and stretch
    by optimizing the zero lag cross-correlation between the inspec1 and the transformed inspec2 (shifted and stretched via
    wvutils.shift_and_stretch()) in a narrow window about the initial estimated shift. The convention for the shift is that
    positive shift means inspec2 is shifted to the right (higher pixel values) relative to inspec1. The convention for the stretch is
    that it is float near unity that increases the size of the inspec2 relative to the original size (which is the size of inspec1)

    Parameters
    ----------
    inspec1 : ndarray
      Reference spectrum
    inspec2 : ndarray
      Spectrum for which the shift and stretch are computed such that it will match inspec1

    Optional Parameters
    -------------------
    cc_thresh: float, default = -1.0
      A number in the range [-1.0,1.0] which is the threshold on the initial cross-correlation coefficient for the shift/stretch.
      If the value of the initial cross-correlation is < cc_thresh the code will just exit and return this value and the best shift.
      This is desirable behavior since the shif/stretch optimization is slow and this allows one to test how correlated the spectra are
      before attempting it, since there is little value in that expensive computation for spectra with little overlap. The default cc_thresh =-1.0
      means shift/stretch is always attempted since the cross correlation coeficcient cannot be less than -1.0.
    smooth: float, default
      Gaussian smoothing in pixels applied to both spectra for the computations. Default is 5.0
    percent_ceil: float, default=90.0
      Appply a ceiling to the input spectra at the percent_ceil percentile level of the distribution of peak amplitudes.
      This prevents extremely strong lines from completely dominating the cross-correlation, which can
      causes the cross-correlation to have spurious noise spikes that are not the real maximum.
    use_raw_arc: bool, default = False
      If this parameter is True the raw arc will be used rather than the continuum subtracted arc
    shift_mnmx: tuple of floats, default = (-0.05,0.05)
      Range to search for the shift in the optimization about the initial cross-correlation based estimate of the shift.
      The optimization will search the window (shift_cc + nspec*shift_mnmx[0],shift_cc + nspec*shift_mnmx[1]) where nspec
      is the number of pixels in the spectrum
    stretch_mnmx: tuple of floats, default = (0.97,1.03)
      Range to search for the stretch in the optimization. The code may not work well if this range is significantly expanded
      because the linear approximation used to transform the arc starts to break down.
    seed: int or np.random.RandomState, optional, default = None
       Seed for scipy.optimize.differential_evolution optimizer. If not specified, the calculation will not be repeatable
    debug = False
       Show plots to the screen useful for debugging.

    Returns
    -------
    success: int
      A flag indicating the exist status.
          success  = 1, shift and stretch performed via sucessful optimization
          success  = 0, shift and stretch optimization failed
          success  = -1, initial x-correlation is below cc_thresh (see above), so shift/stretch optimization was not attempted
    shift: float
      the optimal shift which was determined.
      If cc_thresh is set, and the initial cross-correlation is < cc_thresh,  then this will be just the cross-correlation shift
    stretch: float
      the optimal stretch which was determined.
      If cc_thresh is set, and the initial cross-correlation is < cc_thresh,  then this will be just be 1.0
    cross_corr: float
      the value of the cross-correlation coefficient at the optimal shift and stretch. This is a number between zero and unity,
      which unity indicating a perfect match between the two spectra.
      If cc_thresh is set, and the initial cross-correlation is < cc_thresh, this will be just the initial cross-correlation
    shift_init:
      The initial shift determined by maximizing the cross-correlation coefficient without allowing for a stretch.
      If cc_thresh is set, and the initial cross-correlation is < cc_thresh, this will be just the shift from the initial cross-correlation
    cross_corr_init:
      The maximum of the initial cross-correlation coefficient determined without allowing for a stretch.
      If cc_thresh is set, and the initial cross-correlation is < cc_thresh, this will be just the initial cross-correlation
    """

    nspec = inspec1.size

    y1 = smooth_ceil_cont(inspec1,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc)
    y2 = smooth_ceil_cont(inspec2,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc)

    # Do the cross-correlation first and determine the initial shift
    shift_cc, cc_val = xcorr_shift(y1, y2, smooth = None, percent_ceil = None, debug = debug)

    if cc_val < cc_thresh:
        return -1.0, shift_cc, 1.0, cc_val, shift_cc, cc_val
    else:
        bounds = [(shift_cc + nspec*shift_mnmx[0],shift_cc + nspec*shift_mnmx[1]), stretch_mnmx]
        # TODO Can we make the differential evolution run faster?
        result = scipy.optimize.differential_evolution(zerolag_shift_stretch, args=(y1,y2), tol=1e-4,
                                                       bounds=bounds, disp=False, polish=True, seed=seed)
        if not result.success:
            msgs.warn('Fit for shift and stretch did not converge!')

        if debug:
            x1 = np.arange(nspec)
            inspec2_trans = shift_and_stretch(inspec2, result.x[0], result.x[1])
            plt.plot(x1,inspec1, 'k-', drawstyle='steps', label ='inspec1')
            plt.plot(x1,inspec2_trans, 'r-', drawstyle='steps', label = 'inspec2')
            plt.title('shift= {:5.3f}'.format(result.x[0]) +
                      ',  stretch = {:7.5f}'.format(result.x[1]) + ', corr = {:5.3f}'.format(-result.fun))
            plt.legend()
            plt.show()

        return int(result.success), result.x[0], result.x[1], -result.fun, shift_cc, cc_val



# JFH ToDo This algorithm for computing the shift and stretch is unstable. It was hanging but that has been fixed
# by ading the bounds. However, I think it is producing bogus results in many cases.
def match_peaks_old(inspec1, inspec2, smooth=5.0, debug=False):
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
def shift_stretch_old(specs, p, retshift=False):
    inspec1, inspec2, smooth = specs
    y1 = scipy.ndimage.filters.gaussian_filter(inspec1, smooth)
    y2 = scipy.ndimage.filters.gaussian_filter(inspec2, smooth)
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
