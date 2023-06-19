
import numpy as np
import scipy.ndimage


@nb.jit(nopython=True, cache=True)
def hist_wavedisp(waves, disps, dispbin=None, wavebin=None, scale=1.0, debug=False):
    """

    Generates a flexible 2D histogram of central wavelength and
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

def smooth_ceil_cont(inspec1, smooth, percent_ceil = None, use_raw_arc=False,sigdetect = 10.0, fwhm = 4.0):
    """ Utility routine to smooth and apply a ceiling to spectra """

    # ToDO can we improve the logic here. Technically if use_raw_arc = True and perecent_ceil=None
    # we don't need to peak find or continuum subtract, but this makes the code pretty uggly.

    # Run line detection to get the continuum subtracted arc
    tampl1, tampl1_cont, tcent1, twid1, centerr1, w1, arc1, nsig1 = arc.detect_lines(inspec1, sigdetect=sigdetect, fwhm=fwhm)
    if use_raw_arc == True:
        ampl = tampl1
        use_arc = inspec1
    else:
        ampl = tampl1_cont
        use_arc = arc1

    if percent_ceil is not None and (ampl.size > 0):
        # If this is set, set a ceiling on the greater > 10sigma peaks
        ceil1 = np.percentile(ampl, percent_ceil)
        spec1 = np.fmin(use_arc, ceil1)
    else:
        spec1 = np.copy(use_arc)

    if smooth is not None:
        y1 = scipy.ndimage.filters.gaussian_filter(spec1, smooth)
    else:
        y1 = np.copy(spec1)

    return y1

#This was JFHs keck_hires_dev development version to try to deal with the strongly saturated line regions
def smooth_ceil_cont_hires_dev(inspec1, smooth, percent_ceil = None, use_raw_arc=False, sigdetect = 10.0, fwhm = 4.0,
                     large_scale_fwhm_fact=20.0, large_scale_sigrej=20.0):

    """  Utility routine to smooth and apply a ceiling to spectra

    Args:
        inspec1 (ndarray):
          Input spectrum, shape = (nspec,)
        smooth (int):
          Number of pixels to smooth by
        percent_ceil (float, optional):
          Upper percentile threshold for thresholding positive and negative values
        use_raw_arc (bool, optional):
          If True, use the raw arc and do not continuum subtract. Default = False
        sigdetect (float, optional):
          Peak finding threshold which is used if continuum subtraction will occur (i.e. if use_raw_arc = False)
        fwhm (float, optional):
          Fwhm of arc lines for peak finding if continuum subtraction will occur (i.e. if use_raw_arc = False)
        large_scale_fwhm_fact (float, optional):
          Number of fwhms to use as the width of the running median to take out large scale features from saturated
          lines.
        large_scale_percentile (float, optional):
          Percentile of large_scale_smoothed arc to set as the threshold above which the arc is masked to zero.

    Returns:
        y1_out (ndarray):
          Spectrum with smoothing and ceiling applied, shape = (nspec,)

    """

    # ToDO can we improve the logic here. Technically if use_raw_arc = True and perecent_ceil=None
    # we don't need to peak find or continuum subtract, but this makes the code pretty uggly.

    # Run line detection to get the continuum subtracted arc
    tampl1, tampl1_cont, tcent1, twid1, centerr1, w1, arc1, nsig1 = arc.detect_lines(inspec1, sigdetect=sigdetect, fwhm=fwhm)

    if use_raw_arc:
        ampl = tampl1
        use_arc = inspec1
    else:
        ampl = tampl1_cont
        use_arc = arc1

    if percent_ceil is not None and (ampl.size > 0):
        # If this is set, set a ceiling on the greater > 10sigma peaks
        ceil_upper = np.percentile(ampl[ampl >= 0.0], percent_ceil) if np.any(ampl >= 0.0) else 0.0
        # Set a lower ceiling on negative fluctuations
        ceil_lower = np.percentile(ampl[ampl < 0.0], percent_ceil) if np.any(ampl < 0.0) else 0.0
        spec1 = np.clip(use_arc, ceil_lower, ceil_upper)
    else:
        spec1 = np.copy(use_arc)

    if smooth is not None:
        y1 = scipy.ndimage.filters.gaussian_filter(spec1, smooth)
    else:
        y1 = np.copy(spec1)

    # Mask out large scale features
    y1_ls = utils.fast_running_median(y1, fwhm * large_scale_fwhm_fact)
    mean, med, stddev = astropy.stats.sigma_clipped_stats(
        y1_ls[y1_ls != 0.0], sigma_lower=3.0, sigma_upper=3.0, cenfunc='median', stdfunc=utils.nan_mad_std)
    y1_ls_hi = np.percentile(y1_ls[y1_ls != 0.0], 100.0*scipy.stats.norm.cdf(1.0))
    y1_ls_lo = np.percentile(y1_ls[y1_ls != 0.0], 100.0*scipy.stats.norm.cdf(-1.0))
    stddev_percentile = (y1_ls_hi - y1_ls_lo)/2.0
    if (mean != 0.0) & (med != 0.0) & (stddev_percentile != 0.0):
        y1_ls_bpm = (y1_ls > (med + large_scale_sigrej*stddev_percentile)) | ((y1_ls < (med - large_scale_sigrej*stddev_percentile)))
    else:
        y1_ls_bpm = np.zeros_like(y1_ls, dtype=bool)

    y1_out = y1*np.logical_not(y1_ls_bpm)


    return y1_out



