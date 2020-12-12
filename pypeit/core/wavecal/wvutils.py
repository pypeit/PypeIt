""" Module for basic utilties with holy grail


.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import numpy as np
import os
import numba as nb


from matplotlib import pyplot as plt

from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit

from astropy.table import Table

from pypeit import msgs
from pypeit.core import arc

from IPython import embed

def parse_param(par, key, slit):
    # Find good lines for the tilts
    param_in = par[key]
    if isinstance(param_in, (float, int)):
        param = param_in
    elif isinstance(param_in, (list, np.ndarray)):
        param = param_in[slit]
    else:
        raise ValueError('Invalid input for parameter {:s}'.format(key))

    return param

def get_sampling(waves, pix_per_R=3.0):
    """
    Computes the median wavelength sampling of wavelength vector(s)

    Args:
        waves (float ndarray): shape = (nspec,) or (nspec, nimgs)
            Array of wavelengths. Can be one or two dimensional where
            the nimgs dimension can represent the orders, exposures, or
            slits
        pix_per_R (float):  default=3.0
            Number of pixels per resolution element used for the
            resolution guess. The default of 3.0 assumes roughly Nyquist
            smampling

    Returns:
        tuple: Returns dlam, dloglam, resln_guess, pix_per_sigma.
        Computes the median wavelength sampling of the wavelength
        vector(s)

    """

    if waves.ndim == 1:
        norders = 1
        nspec = waves.shape[0]
        waves_stack = waves.reshape((nspec, norders))
    elif waves.ndim == 2:
        waves_stack = waves
    elif waves.ndim == 3:
        nspec, norder, nexp = waves.shape
        waves_stack = np.reshape(waves, (nspec, norder * nexp), order='F')
    else:
        msgs.error('The shape of your wavelength array does not make sense.')

    wave_mask = waves_stack > 1.0
    waves_ma = np.ma.array(waves_stack, mask=np.invert(wave_mask))
    loglam = np.ma.log10(waves_ma)
    wave_diff = np.diff(waves_ma, axis=0)
    loglam_diff = np.diff(loglam, axis=0)
    dwave = np.ma.median(wave_diff)
    dloglam = np.ma.median(loglam_diff)
    #dloglam_ord = np.ma.median(loglam_roll, axis=0)
    #dloglam = np.median(dloglam_ord)
    resln_guess = 1.0 / (pix_per_R* dloglam * np.log(10.0))
    pix_per_sigma = 1.0 / resln_guess / (dloglam * np.log(10.0)) / (2.0 * np.sqrt(2.0 * np.log(2)))
    return dwave, dloglam, resln_guess, pix_per_sigma


def arc_lines_from_spec(spec, sigdetect=10.0, fwhm=4.0,
                        fit_frac_fwhm = 1.25, cont_frac_fwhm=1.0,max_frac_fwhm=2.0,
                        cont_samp=30, niter_cont=3,nonlinear_counts=1e10, debug=False):
    """
    Simple wrapper to arc.detect_lines.
    See that code for docs

    Args:
        spec:
        sigdetect:
        fwhm:
        fit_frac_fwhm:
        cont_frac_fwhm:
        max_frac_fwhm:
        cont_samp:
        niter_cont:
        nonlinear_counts (float, optional):
            Counts where the arc is presumed to go non-linear
        debug:

    Returns:
        tuple: all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub.
        See arc.detect_lines for details

    """

    # Find peaks
    tampl, tampl_cont, tcent, twid, centerr, w, arc_cont_sub, nsig = arc.detect_lines(
        spec, sigdetect = sigdetect, fwhm=fwhm, fit_frac_fwhm=fit_frac_fwhm,
        cont_frac_fwhm=cont_frac_fwhm, max_frac_fwhm=max_frac_fwhm,
        cont_samp=cont_samp,niter_cont = niter_cont, nonlinear_counts = nonlinear_counts,
        debug=debug)
    all_tcent = tcent[w]
    all_ecent = centerr[w]
    all_nsig = nsig[w]

    # Cut on significance
    cut_sig = all_nsig > sigdetect
    cut_tcent = all_tcent[cut_sig]
    icut = np.where(cut_sig)[0]

    # Return
    return all_tcent, all_ecent, cut_tcent, icut, arc_cont_sub


def shift_and_stretch(spec, shift, stretch):

    """
    Utility function to shift and stretch a spectrum. This operation is
    being implemented in many steps and could be significantly
    optimized. But it works for now. Note that the stretch is applied
    *first* and then the shift is applied in stretch coordinates.

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

    """
    Utility function which is run by the differential evolution
    optimizer in scipy. These is the fucntion we optimize.  It is the
    zero lag cross-correlation coefficient of spectrum with a shift and
    stretch applied.

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
        Negative of the zero lag cross-correlation coefficient (since we
        are miniziming with scipy.optimize). scipy.optimize will thus
        determine the shift,stretch that maximize the cross-correlation.

    """


    shift, stretch = theta
    y2_corr = shift_and_stretch(y2, shift, stretch)
    # Zero lag correlation
    corr_zero = np.sum(y1*y2_corr)
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr_zero/corr_denom
    return -corr_norm

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



# ToDO can we speed this code up? I've heard numpy.correlate is faster. Someone should investigate optimization. Also we don't need to compute
# all these lags.
def xcorr_shift(inspec1,inspec2, smooth=1.0, percent_ceil=80.0, use_raw_arc=False, sigdetect=10.0, fwhm=4.0, debug=False):

    """ Determine the shift inspec2 relative to inspec1.  This routine computes the shift by finding the maximum of the
    the cross-correlation coefficient. The convention for the shift is that positive shift means inspec2 is shifted to the right
    (higher pixel values) relative to inspec1.

    Args:
        inspec1 : ndarray
            Reference spectrum
        inspec2 : ndarray
            Spectrum for which the shift and stretch are computed such
            that it will match inspec1
        smooth: float, default=1.0
            Gaussian smoothing in pixels applied to both spectra for the
            computations. Default is 5.0
        percent_ceil: float, default=90.0
            Apply a ceiling to the input spectra at the percent_ceil
            percentile level of the distribution of peak amplitudes.
            This prevents extremely strong lines from completely
            dominating the cross-correlation, which can causes the
            cross-correlation to have spurious noise spikes that are not
            the real maximum.
        use_raw_arc: bool, default = False
            If this parameter is True the raw arc will be used rather
            than the continuum subtracted arc
        debug: boolean, default = False

    Returns:
       tuple: Returns the following:

            - shift: float; the shift which was determined
            - cross_corr: float; the maximum of the cross-correlation
              coefficient at this shift

    """

    y1 = smooth_ceil_cont(inspec1,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc, sigdetect = sigdetect, fwhm = fwhm)
    y2 = smooth_ceil_cont(inspec2,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc, sigdetect = sigdetect, fwhm = fwhm)

    nspec = y1.shape[0]
    lags = np.arange(-nspec + 1, nspec)
    corr = scipy.signal.correlate(y1, y2, mode='full')
    corr_denom = np.sqrt(np.sum(y1*y1)*np.sum(y2*y2))
    corr_norm = corr/corr_denom
    tampl_true, tampl, pix_max, twid, centerr, ww, arc_cont, nsig = arc.detect_lines(corr_norm, sigdetect=3.0,
                                                                                     fit_frac_fwhm=1.5, fwhm=5.0,
                                                                                     cont_frac_fwhm=1.0, cont_samp=30, nfind=1)
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


def xcorr_shift_stretch(inspec1, inspec2, cc_thresh=-1.0, smooth=1.0, percent_ceil=80.0, use_raw_arc=False,
                        shift_mnmx=(-0.05,0.05), stretch_mnmx=(0.95,1.05), sigdetect = 10.0, fwhm = 4.0,debug=False, seed = None):

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
    cc_thresh: float, default = -1.0
        A number in the range [-1.0,1.0] which is the threshold on the
        initial cross-correlation coefficient for the shift/stretch.  If
        the value of the initial cross-correlation is < cc_thresh the
        code will just exit and return this value and the best shift.
        This is desirable behavior since the shif/stretch optimization
        is slow and this allows one to test how correlated the spectra
        are before attempting it, since there is little value in that
        expensive computation for spectra with little overlap. The
        default cc_thresh =-1.0 means shift/stretch is always attempted
        since the cross correlation coeficcient cannot be less than
        -1.0.
    smooth: float, default
        Gaussian smoothing in pixels applied to both spectra for the computations. Default is 5.0
    percent_ceil: float, default=90.0
        Apply a ceiling to the input spectra at the percent_ceil
        percentile level of the distribution of peak amplitudes.  This
        prevents extremely strong lines from completely dominating the
        cross-correlation, which can causes the cross-correlation to
        have spurious noise spikes that are not the real maximum.
    use_raw_arc: bool, default = False
        If this parameter is True the raw arc will be used rather than the continuum subtracted arc
    shift_mnmx: tuple of floats, default = (-0.05,0.05)
        Range to search for the shift in the optimization about the
        initial cross-correlation based estimate of the shift.  The
        optimization will search the window (shift_cc +
        nspec*shift_mnmx[0],shift_cc + nspec*shift_mnmx[1]) where nspec
        is the number of pixels in the spectrum
    stretch_mnmx: tuple of floats, default = (0.97,1.03)
        Range to search for the stretch in the optimization. The code
        may not work well if this range is significantly expanded
        because the linear approximation used to transform the arc
        starts to break down.
    seed: int or np.random.RandomState, optional, default = None
        Seed for scipy.optimize.differential_evolution optimizer. If not
        specified, the calculation will not be repeatable
    debug = False
       Show plots to the screen useful for debugging.

    Returns
    -------
    success: int
        A flag indicating the exist status.  Values are:

          - success = 1, shift and stretch performed via sucessful
            optimization
          - success = 0, shift and stretch optimization failed
          - success = -1, initial x-correlation is below cc_thresh (see
            above), so shift/stretch optimization was not attempted

    shift: float
        the optimal shift which was determined.  If cc_thresh is set,
        and the initial cross-correlation is < cc_thresh,  then this
        will be just the cross-correlation shift
    stretch: float
        the optimal stretch which was determined.  If cc_thresh is set,
        and the initial cross-correlation is < cc_thresh,  then this
        will be just be 1.0
    cross_corr: float
        the value of the cross-correlation coefficient at the optimal
        shift and stretch. This is a number between zero and unity,
        which unity indicating a perfect match between the two spectra.
        If cc_thresh is set, and the initial cross-correlation is <
        cc_thresh, this will be just the initial cross-correlation
    shift_init:
        The initial shift determined by maximizing the cross-correlation
        coefficient without allowing for a stretch.  If cc_thresh is
        set, and the initial cross-correlation is < cc_thresh, this will
        be just the shift from the initial cross-correlation
    cross_corr_init:
        The maximum of the initial cross-correlation coefficient
        determined without allowing for a stretch.  If cc_thresh is set,
        and the initial cross-correlation is < cc_thresh, this will be
        just the initial cross-correlation

    """

    nspec = inspec1.size

    y1 = smooth_ceil_cont(inspec1,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc, sigdetect = sigdetect, fwhm = fwhm)
    y2 = smooth_ceil_cont(inspec2,smooth,percent_ceil=percent_ceil,use_raw_arc=use_raw_arc, sigdetect = sigdetect, fwhm = fwhm)

    # Do the cross-correlation first and determine the initial shift
    shift_cc, corr_cc = xcorr_shift(y1, y2, smooth = None, percent_ceil = None, use_raw_arc = True, sigdetect = sigdetect, fwhm=fwhm, debug = debug)

    if corr_cc < cc_thresh:
        return -1, shift_cc, 1.0, corr_cc, shift_cc, corr_cc
    else:
        bounds = [(shift_cc + nspec*shift_mnmx[0],shift_cc + nspec*shift_mnmx[1]), stretch_mnmx]
        # TODO Can we make the differential evolution run faster?
        result = scipy.optimize.differential_evolution(zerolag_shift_stretch, args=(y1,y2), tol=1e-4,
                                                       bounds=bounds, disp=False, polish=True, seed=seed)
        corr_de = -result.fun
        shift_de = result.x[0]
        stretch_de = result.x[1]
        if not result.success:
            msgs.warn('Fit for shift and stretch did not converge!')

        if(corr_de < corr_cc):
            # Occasionally the differential evolution crapps out and returns a value worse that the CC value. In these cases just use the cc value
            msgs.warn('Shift/Stretch optimizer performed worse than simple x-correlation.' +
                      'Returning simple x-correlation shift and no stretch:' + msgs.newline() +
                      '   Optimizer: corr={:5.3f}, shift={:5.3f}, stretch={:7.5f}'.format(corr_de, shift_de,stretch_de) + msgs.newline() +
                      '     X-corr : corr={:5.3f}, shift={:5.3f}'.format(corr_cc,shift_cc))
            corr_out = corr_cc
            shift_out = shift_cc
            stretch_out = 1.0
            result_out = 1
        else:
            corr_out = corr_de
            shift_out = shift_de
            stretch_out = stretch_de
            result_out = int(result.success)

        if debug:
            x1 = np.arange(nspec)
            y2_trans = shift_and_stretch(y2, shift_out, stretch_out)
            plt.figure(figsize=(14, 6))
            plt.plot(x1,y1, 'k-', drawstyle='steps', label ='inspec1, input spectrum')
            plt.plot(x1,y2_trans, 'r-', drawstyle='steps', label = 'inspec2, reference shift & stretch')
            plt.title('shift= {:5.3f}'.format(shift_out) +
                      ',  stretch = {:7.5f}'.format(stretch_out) + ', corr = {:5.3f}'.format(corr_out))
            plt.legend()
            plt.show()

        return result_out, shift_out, stretch_out, corr_out, shift_cc, corr_cc



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


def wavegrid(wave_min, wave_max, dwave, samp_fact=1.0, log10=False):
    """

    Utility routine to generate a uniform grid of wavelengths

    Args:
        wave_min (float):
           Mininum wavelength. Must be linear even if log10 is requested
        wave_max (float):
           Maximum wavelength. Must be linear even if log10 is requested.
        dwave (float):
           Delta wavelength interval. Must be linear if log10=False, or log10 if log10=True
        samp_fact (float):
           sampling factor to make the wavelength grid finer or coarser.
           samp_fact > 1.0 oversamples (finer), samp_fact < 1.0
           undersamples (coarser)

    Returns:
        `numpy.ndarray`_: Wavelength grid in Angstroms (i.e. log10 even
        if log10 is requested)

    """

    dwave_eff = dwave/samp_fact
    if log10:
        ngrid = np.ceil((np.log10(wave_max) - np.log10(wave_min))/dwave_eff).astype(int)
        loglam_grid = np.log10(wave_min) + dwave_eff*np.arange(ngrid)
        return np.power(10.0,loglam_grid)
    else:
        ngrid = np.ceil((wave_max - wave_min)/dwave_eff).astype(int)
        return wave_min + dwave_eff*np.arange(ngrid)

    return wave_grid


def write_template(nwwv, nwspec, binspec, outpath, outroot, det_cut=None,
                   order=None, overwrite=True):
    """
    Write the template spectrum into a binary FITS table

    Args:
        nwwv (`numpy.ndarray`_):
            Wavelengths for the template
        nwspec (`numpy.ndarray`_):
            Flux of the template
        binspec (int):
            Binning of the template
        outpath (str):
        outroot (str):
        det_cut (bool, optional):
            Cuts in wavelength for detector snippets
            Used primarily for DEIMOS
        order (`numpy.ndarray`_, optional):
            Echelle order numbers
        overwrite (bool, optional):
            If True, overwrite any existing file
    """
    tbl = Table()
    tbl['wave'] = nwwv
    tbl['flux'] = nwspec
    if order is not None:
        tbl['order'] = order

    tbl.meta['BINSPEC'] = binspec
    # Detector snippets??
    if det_cut is not None:
        tbl['det'] = 0
        for dets, wcuts in zip(det_cut['dets'], det_cut['wcuts']):
            gdwv = (tbl['wave'] > wcuts[0]) & (tbl['wave'] < wcuts[1])
            deti = np.sum([2**ii for ii in dets])
            tbl['det'][gdwv] += deti
    # Write
    outfile = os.path.join(outpath, outroot)
    tbl.write(outfile, overwrite=overwrite)
    print("Wrote: {}".format(outfile))
