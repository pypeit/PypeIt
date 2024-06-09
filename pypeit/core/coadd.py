"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import os
import sys
import copy
import string

from IPython import embed

import numpy as np
import scipy

import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter, NullLocator, MaxNLocator

from astropy import stats
from astropy import convolution

from pypeit import utils
from pypeit.core import fitting
from pypeit import specobjs
from pypeit import msgs
from pypeit.core import combine
from pypeit.core.wavecal import wvutils
from pypeit.core import pydl
from pypeit import data


def renormalize_errors_qa(chi, maskchi, sigma_corr, sig_range = 6.0,
                          title:str='', qafile:str=None):
    '''
    Generate a histogram QA plot of the input chi distribution.

    Args:
        chi (`numpy.ndarray`_):
            your chi values
        maskchi (`numpy.ndarray`_):
            True = good, mask for your chi array of type bool
        sigma_corr (float):
            corrected sigma
        sig_range (float):
            used to set binsize, default +- 6-sigma
        title (str, optional):
            plot title
        qafile (str, optional):
            Write figure to this output QA file, if provided
    '''
    # Prep
    n_bins = 50
    binsize = 2.0*sig_range/n_bins
    bins_histo = -sig_range + np.arange(n_bins)*binsize+binsize/2.0

    xvals = np.arange(-10.0,10,0.02)
    gauss = scipy.stats.norm(loc=0.0,scale=1.0)
    gauss_corr = scipy.stats.norm(loc=0.0,scale=sigma_corr)

    # Plot
    plt.figure(figsize=(12, 8))
    plt.hist(chi[maskchi],bins=bins_histo,density=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,gauss.pdf(xvals),'c-',lw=3,label='sigma=1')
    plt.plot(xvals,gauss_corr.pdf(xvals),'m--',lw=2,label='new sigma={:4.2f}'.format(round(sigma_corr,2)))
    plt.ylabel('Residual distribution')
    plt.xlabel('chi')
    plt.xlim([-6.05,6.05])
    plt.legend(fontsize=13,loc=2)
    plt.title(title, fontsize=16, color='red')
    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote QA: {:s}".format(qafile))
    plt.show()
    plt.close()


def renormalize_errors(chi, mask, clip=6.0, max_corr=5.0, title = '', debug=False):
    """
    Function for renormalizing errors. The distribution of input chi (defined by chi = (data - model)/sigma) values is
    analyzed, and a correction factor to the standard deviation sigma_corr is returned. This should be multiplied into
    the errors. In this way, a rejection threshold of i.e. 3-sigma, will always correspond to roughly the same percentile.
    This renormalization guarantees that rejection is not too agressive in cases where the empirical errors determined
    from the chi-distribution differ significantly from the noise model which was used to determine chi.

    Args:
        chi (`numpy.ndarray`_):
            input chi values
        mask (`numpy.ndarray`_):
            True = good, mask for your chi array of type bool
        mask (`numpy.ndarray`_):
            True = good, mask for your chi array of type bool
        clip (float, optional):
            threshold for outliers which will be clipped for the purpose of computing the renormalization factor
        max_corr (float, optional):
            maximum corrected sigma allowed.
        title (str, optional):
            title for QA plot, passed to renormalize_errors_qa
        debug (bool, optional):
            If True, show the QA plot created by renormalize_errors_qa

    Returns:
        tuple: (1) sigma_corr (float), corrected new sigma; (2) maskchi
        (`numpy.ndarray`_, bool): new mask (True=good) which indicates the values
        used to compute the correction (i.e it includes clipping)

    """
    chi2 = chi**2
    maskchi = (chi2 < clip**2) & mask
    if (np.sum(maskchi) > 0):
        gauss_prob = 1.0 - 2.0 * scipy.stats.norm.cdf(-1.0)
        chi2_sigrej = np.percentile(chi2[maskchi], 100.0*gauss_prob)
        sigma_corr = np.sqrt(chi2_sigrej)
        if sigma_corr < 1.0:
            msgs.warn("Error renormalization found correction factor sigma_corr = {:f}".format(sigma_corr) +
                      " < 1." + msgs.newline() +
                      " Errors are overestimated so not applying correction")
            sigma_corr = 1.0
        if sigma_corr > max_corr:
            msgs.warn(("Error renormalization found sigma_corr/sigma = {:f} > {:f}." + msgs.newline() +
                      "Errors are severely underestimated." + msgs.newline() +
                      "Setting correction to sigma_corr = {:4.2f}").format(sigma_corr, max_corr, max_corr))
            sigma_corr = max_corr

        if debug:
            renormalize_errors_qa(chi, maskchi, sigma_corr, title=title)

    else:
        msgs.warn('No good pixels in error_renormalize. There are probably issues with your data')
        sigma_corr = 1.0

    return sigma_corr, maskchi

def poly_model_eval(theta, func, model, wave, wave_min, wave_max):
    """
    Routine to evaluate the polynomial fit

    Args:
        theta (`numpy.ndarray`_):
           coefficient parameter vector of type=float
        func (str):
           polynomial type
        model (str):
           model type, valid model types are 'poly', 'square', or 'exp', corresponding to normal polynomial,
           squared polynomial, or exponentiated polynomial
        wave (`numpy.ndarray`_):
           array of wavelength values of type=float
        wave_min (float):
           minimum wavelength for polynomial fit range
        wave_max (float):
           maximum wavelength for polynomial fit range

    Returns:
        `numpy.ndarray`_: Array of evaluated polynomial with same shape as wave
    """
    # Evaluate the polynomial for rescaling
    if 'poly' in model:
        ymult = fitting.evaluate_fit(theta, func, wave, minx=wave_min, maxx=wave_max)
    elif 'square' in model:
        ymult = (fitting.evaluate_fit(theta, func, wave, minx=wave_min, maxx=wave_max)) ** 2
    elif 'exp' in model:
        # Clipping to avoid overflow.
        ymult = np.exp(np.clip(fitting.evaluate_fit(theta, func, wave, minx=wave_min, maxx=wave_max)
                               , None, 0.8 * np.log(sys.float_info.max)))

    else:
        msgs.error('Unrecognized value of model requested')

    return ymult


def poly_ratio_fitfunc_chi2(theta, gpm, arg_dict):
    """
    Function for computing the chi^2 loss function for solving for the polynomial rescaling of one spectrum to another.
    There are two non-standard things implemented here which increase ther robustness. The first is a non-standard error used for the
    chi, which adds robustness and increases the stability of the optimization. This was taken from the idlutils
    solve_poly_ratio code. The second thing is that the chi is remapped using the scipy huber loss function to
    reduce sensitivity to outliers, based on the scipy cookbook on robust optimization.

    Args:
        theta (`numpy.ndarray`_): parameter vector for the polymomial fit
        gpm (`numpy.ndarray`_): boolean mask for the current iteration of the optimization, True=good
        arg_dict (dict): dictionary containing arguments

    Returns:
        float: this is effectively the chi^2, i.e. the quantity to be
        minimized by the optimizer. Note that this is not formally the
        chi^2 since the huber loss function re-maps the chi to be less
        sensitive to outliers.
    """

    # Unpack the data to be rescaled, the mask for the reference spectrum, and the wavelengths
    mask = arg_dict['mask']
    flux_med = arg_dict['flux_med']
    ivar_med = arg_dict['ivar_med']
    flux_ref_med = arg_dict['flux_ref_med']
    ivar_ref_med = arg_dict['ivar_ref_med']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    model = arg_dict['model']
    ymult = poly_model_eval(theta, func, model, wave, wave_min, wave_max)

    flux_scale = ymult*flux_med
    mask_both = mask & gpm
    # This is the formally correct ivar used for the rejection, but not used in the fitting. This appears to yield
    # unstable results
    #totvar = utils.inverse(ivar_ref, positive=True) + ymult**2*utils.inverse(ivar, positive=True)
    #ivartot = mask_both*utils.inverse(totvar, positive=True)

    # The errors are rescaled at every function evaluation, but we only allow the errors to get smaller by up to a
    # factor of 1e4, and we only allow them to get larger slowly (as the square root).  This should very strongly
    # constrain the flux-corrrection vectors from going too small (or negative), or too large.
    ## Schlegel's version here
    vmult = np.fmax(ymult,1e-4)*(ymult <= 1.0) + np.sqrt(ymult)*(ymult > 1.0)
    ivarfit = mask_both/(1.0/(ivar_med + np.logical_not(mask_both)) + np.square(vmult)/(ivar_ref_med + np.logical_not(mask_both)))
    chi_vec = mask_both * (flux_ref_med - flux_scale) * np.sqrt(ivarfit)
    # Changing the Huber loss parameter from step to step results in instability during optimization --MSR.
    # Robustly characterize the dispersion of this distribution
    #chi_mean, chi_median, chi_std = stats.sigma_clipped_stats(
    #    chi_vec, np.invert(mask_both), cenfunc='median', stdfunc=utils.nan_mad_std, maxiters=5, sigma=2.0)
    chi_std = np.std(chi_vec)
    # The Huber loss function smoothly interpolates between being chi^2/2 for standard chi^2 rejection and
    # a linear function of residual in the outlying tails for large residuals. This transition occurs at the
    # value of the first argument, which we have set to be 2.0*chi_std, which is 2-sigma given the modified
    # errors described above from Schlegel's code.
    robust_scale = 2.0
    huber_vec = scipy.special.huber(robust_scale*chi_std, chi_vec)
    loss_function = np.sum(huber_vec*mask_both)
    #chi2 = np.sum(np.square(chi_vec))
    return loss_function

def poly_ratio_fitfunc(flux_ref, gpm, arg_dict, init_from_last=None, **kwargs_opt):
    """
    Function to be optimized by robust_optimize for solve_poly_ratio
    polynomial rescaling of one spectrum to match a reference
    spectrum. This function has the correct format for running
    robust_optimize optimization. In addition to running the
    optimization, this function recomputes the error vector ivartot
    for the error rejection that takes place at each iteration of the
    robust_optimize optimization. The ivartot is also renormalized
    using the renormalize_errors function enabling rejection. A scale
    factor is multiplied into the true errors to allow one to reject
    based on the statistics of the actual error distribution.

    Args:
        flux_ref (`numpy.ndarray`_):
            Reference flux that we are trying to rescale our spectrum
            to match
        gpm (`numpy.ndarray`_):
            Boolean array with mask for the current iteration of the
            optimization. True=good
        arg_dict (:obj:`dict`):
            dictionary containing arguments for the optimizing
            function. See poly_ratio_fitfunc_chi2 for how arguments
            are used. They are mask, flux_med, flux_ref_med,
            ivar_ref_med, wave, wave_min, wave_max, func
        init_from_last (obj, optional):
            Use this scipy optimization object from a previous iteration as the guess
        kwargs_opt (:obj:`dict`):
            arguments to be passed to the optimizer, which in this
            case is just vanilla scipy.minimize with the default
            optimizer

   Returns:
       tuple:
         Three objects are returned. (1) scipy optimization object,
         (2) scale factor to be applied to the data to match the
         reference spectrum flux_ref, (3) error vector to be used for
         the rejection that takes place at each iteration of the
         robust_optimize optimization

    """

    # flux_ref, ivar_ref act like the 'data', the rescaled flux will be the 'model'

    guess = arg_dict['guess'] if init_from_last is None else init_from_last.x
    result = scipy.optimize.minimize(poly_ratio_fitfunc_chi2, guess, args=(gpm, arg_dict),  **kwargs_opt)
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    model = arg_dict['model']
    # Evaluate the polynomial for rescaling
    ymult = poly_model_eval(result.x, func, model, wave, wave_min, wave_max)
    flux_scale = ymult*flux
    mask_both = mask & gpm
    totvar = utils.inverse(ivar_ref) + ymult**2*utils.inverse(ivar)
    ivartot1 = mask_both*utils.inverse(totvar)
    # Now rescale the errors
    chi = (flux_scale - flux_ref)*np.sqrt(ivartot1)
    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    sigma_corr, maskchi = renormalize_errors(chi, mask=mask_both, title = 'poly_ratio_fitfunc', debug=debug)
    ivartot = ivartot1/sigma_corr**2

    return result, flux_scale, ivartot

def median_filt_spec(flux, ivar, gpm, med_width):
    '''
    Utility routine to median filter a spectrum using the mask and propagating the errors using the
    utils.fast_running_median function.

    Parameters
    ----------
    flux : `numpy.ndarray`_
            flux array with shape (nspec,)
    ivar : `numpy.ndarray`_
            inverse variance with shape (nspec,)
    gpm : `numpy.ndarray`_
            Boolean mask on the spectrum with shape (nspec,). True = good
    med_width : float
            width for median filter in pixels

    Returns
    -------
    flux_med : `numpy.ndarray`_
        Median filtered flux
    ivar_med : `numpy.ndarray`_
        corresponding propagated variance
    '''

    flux_med = np.zeros_like(flux)
    ivar_med = np.zeros_like(ivar)
    flux_med0 = utils.fast_running_median(flux[gpm], med_width)
    flux_med[gpm] = flux_med0
    var = utils.inverse(ivar)
    var_med0 =  utils.smooth(var[gpm], med_width)
    ivar_med[gpm] = utils.inverse(var_med0)
    return flux_med, ivar_med

def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None,
                     scale_min = 0.05, scale_max = 100.0, func='legendre', model ='square',
                     maxiter=3, sticky=True, lower=3.0, upper=3.0, median_frac=0.01,
                     ref_percentile=70.0, debug=False):
    """
    Routine for solving for the polynomial rescaling of an input
    spectrum flux to match a reference spectrum flux_ref. The two
    spectra need to be defined on the same wavelength grid. The code
    will work best if you choose the reference to be the higher S/N
    ratio spectrum. Note that the code multiplies in the square of a
    polnomial of order norder to ensure positivity of the scale
    factor. It also operates on median filtered spectra to be more
    robust against outliers

    Parameters
    ----------
    wave : `numpy.ndarray`_
        wavelength array of shape (nspec,). flux, ivar, flux_ref, and ivar_ref
        must all be on the same wavelength grid
    flux : `numpy.ndarray`_
        flux that you want to rescale to match flux_ref
    ivar : `numpy.ndarray`_
        inverse varaiance of the array that you want to rescale to match flux_ref
    flux_ref : `numpy.ndarray`_
        reference flux that you want to rescale flux to match.
    ivar_ref : `numpy.ndarray`_
        inverse variance for reference flux
    norder : int
        Order of polynomial rescaling; norder=1 is a linear fit and norder must
        be >= 1 otherwise the code will fault.
    mask : `numpy.ndarray`_, optional
        boolean mask for spectrum that you want to rescale, True=Good
    mask_ref : `numpy.ndarray`_, optional
        boolean mask for reference flux
    scale_min : float, optional
        minimum scaling factor allowed. default =0.05
    scale_max : float, optional
        maximum scaling factor allowed. default=100.0
    func : str, optional
        function you want to use. default='legendre'
    model : str, optional
        model type, valid model types are 'poly', 'square', or 'exp',
        corresponding to normal polynomial, squared polynomial, or exponentiated
        polynomial. default = 'square'
    maxiter : int, optional
        maximum number of iterations for robust_optimize. default=3
    sticky : bool, optional
        whether you want the rejection to be sticky or not with robust_optimize.
        See docs for djs_reject for definition of sticky.  default=True
    lower : float, optional
        lower sigrej rejection threshold for robust_optimize. default=3.0
    upper : float, optional
        upper sigrej rejection threshold for robust_optimize. default=3.0
    median_frac : float, optional
        the code rescales median filtered spectra with 'reflect' boundary
        conditions. The with of the median filter will be median_frac*nspec,
        where nspec is the number of spectral pixels.  default = 0.01,
    debug : bool, optional
        If True, show interactive QA plot. default=False

    Returns
    -------
    ymult : `numpy.ndarray`_, (nspec,)
        rescaling factor to be multiplied into flux to match flux_ref.
    fit_tuple : :obj:`tuple`
        Tuple with the polynomial coefficients, the minimum wavelength
        coordinate and maximum wavelength coordinate used in the fit.
    flux_rescale : `numpy.ndarray`_, (nspec,)
        rescaled flux, i.e. ymult multiplied into flux.
    ivar_rescale : `numpy.ndarray`_, (nspec,)
        rescaled inverse variance
    outmask : `numpy.ndarray`_, bool, (nspec,)
        output mask determined from the robust_optimize optimization/rejection
        iterations. True=Good
    """

    if norder < 1:
        msgs.error('You cannot solve for the polynomial ratio for norder < 1. For rescaling by a constant use robust_median_ratio')

    if mask is None:
        mask = (ivar > 0.0)
    if mask_ref is None:
        mask_ref = (ivar_ref > 0.0)

    #
    nspec = wave.size
    # Determine an initial guess
    ratio = robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask, mask_ref=mask_ref,
                                ref_percentile=ref_percentile, max_factor=scale_max)
    # guess = np.append(ratio, np.zeros(norder))
    wave_min = wave.min()
    wave_max = wave.max()

    # Now compute median filtered versions of the spectra which we will actually operate on for the fitting. Note
    # that rejection will however work on the non-filtered spectra.
    med_width = (2.0*np.ceil(median_frac/2.0*nspec) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    flux_ref_med, ivar_ref_med = median_filt_spec(flux_ref, ivar_ref, mask_ref, med_width)

    if 'poly' in model:
        guess = np.append(ratio, np.zeros(norder))
    elif 'square' in model:
        guess = np.append(np.sqrt(ratio), np.zeros(norder))
    elif 'exp' in model:
        guess = np.append(np.log(ratio), np.zeros(norder))
    else:
        msgs.error('Unrecognized model type')

    ## JFH I'm not convinced any of this below is right or necessary. Going back to previous logic but
    ## leaving this here for now

    # Use robust_fit to get a best-guess linear fit as the starting point. The logic below deals with whether
    # we re fitting a polynomial model to the data model='poly', to the square model='square', or taking the exponential
    # of a polynomial fit model='exp'
    #if 'poly' in model:
    #    #guess = np.append(ratio, np.zeros(norder))
    #    yval = flux_ref_med
    #    yval_ivar = ivar_ref_med
    #    scale_mask = np.ones_like(flux_ref_med, dtype=bool) & (wave > 1.0)
    #elif 'square' in model:
    #    #guess = np.append(np.sqrt(ratio), np.zeros(norder))
    #    yval = np.sqrt(flux_ref_med + (flux_ref_med < 0))
    #    yval_ivar = 4.0*flux_ref_med*ivar_ref_med
    #    scale_mask = (flux_ref_med >= 0) & (wave > 1.0)
    #elif 'exp' in model:
    #    #guess = np.append(np.log(ratio), np.zeros(norder))
    #    yval = np.log(flux_ref_med + (flux_ref_med <= 0))
    #    yval_ivar = flux_ref_med**2*ivar_ref_med
    #    scale_mask = (flux_ref_med > 0) & (wave > 1.0)
    #else:
    #    msgs.error('Unrecognized model type')

    #pypfit = fitting.robust_fit(wave, yval, 1, function=func, in_gpm=scale_mask, invvar=yval_ivar,
    #           sticky=False, use_mad=False, debug=debug, upper=3.0, lower=3.0)
    #guess = np.append(pypfit.fitc, np.zeros(norder - 2)) if norder > 1 else pypfit.fitc

    arg_dict = dict(flux = flux, ivar = ivar, mask = mask,
                    flux_med = flux_med, ivar_med = ivar_med,
                    flux_ref_med = flux_ref_med, ivar_ref_med = ivar_ref_med,
                    ivar_ref = ivar_ref, wave = wave, wave_min = wave_min,
                    wave_max = wave_max, func = func, model=model, norder = norder, guess = guess, debug=debug)

    result, ymodel, ivartot, outmask = fitting.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict, inmask=mask_ref,
                                                             maxiter=maxiter, lower=lower, upper=upper, sticky=sticky)
    ymult1 = poly_model_eval(result.x, func, model, wave, wave_min, wave_max)
    ymult = np.fmin(np.fmax(ymult1, scale_min), scale_max)
    flux_rescale = ymult*flux
    ivar_rescale = ivar/ymult**2
    if debug:
        # Determine the y-range for the QA plots
        scale_spec_qa(wave, flux_med, ivar_med, wave, flux_ref_med, ivar_ref_med, ymult, 'poly', mask = mask, mask_ref=mask_ref,
                      title='Median Filtered Spectra that were poly_ratio Fit')

    return ymult, (result.x, wave_min, wave_max), flux_rescale, ivar_rescale, outmask


def interp_oned(wave_new, wave_old, flux_old, ivar_old, gpm_old, log10_blaze_function=None, sensfunc=False, kind='cubic'):
    """
    Interpolate a 1D spectrum onto a new wavelength grid.

    Interpolation is done using `scipy.interpolate.interp1d` with ``cubic``
    interpolation. Any wavelengths in ``wave_new`` that are beyond the range
    of ``wave_old`` are set to ``np.nan`` and masked via the output
    good-pixel mask.

    .. warning::

        Any wavelength in ``wave_old`` that is less than 1 is assumed to
        indicate that the wavelength is invalid!

    Args:
        wave_new (`numpy.ndarray`_):
            New wavelength grid for the output spectra.  Must be 1D.
        wave_old (`numpy.ndarray`_):
            Old wavelength grid.  Must be 1D, need not have the same size as
            ``wave_new``.
        flux_old (`numpy.ndarray`_):
            Old flux on the wave_old grid.  Shape must match ``wave_old``.
        ivar_old (`numpy.ndarray`_):
            Old ivar on the wave_old grid.  Shape must match ``wave_old``.
        gpm_old (`numpy.ndarray`_):
            Old good-pixel mask (True=Good) on the wave_old grid.  Shape must
            match ``wave_old``.
        log10_blaze_function: `numpy.ndarray`_ or None, optional
            Log10 of the blaze function. Shape must match ``wave_old``. Default=None.

        sensfunc (:obj:`bool`, optional):
            If True, the quantities ``flux*delta_wave`` and the corresponding
            ``ivar/delta_wave**2`` will be interpolated and returned instead of
            ``flux`` and ``ivar``. This is useful for sensitivity function
            computation where we need flux*(wavelength bin width). Beacause
            delta_wave is a difference of the wavelength grid, interpolating
            in the presence of masked data requires special care.
        kind (:obj:`int`, :obj:`str`, optional):
            Specifies the kind of interpolation as a string or as an integer
            specifying the order of the spline interpolator to use following the convention of
            scipy.interpolate.interp1d. The string has to be one of 'linear', 'nearest',
            'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', or 'next'. 'zero',
            'slinear', 'quadratic' and 'cubic' refer to a spline interpolation of
            zeroth, first, second or third order; 'previous' and 'next' simply
            return the previous or next value of the point; 'nearest-up' and
            'nearest' differ when interpolating half-integers (e.g. 0.5, 1.5)
            in that 'nearest-up' rounds up and 'nearest' rounds down. Default is 'cubic'.
    Returns:
        :obj:`tuple`: Returns four objects flux_new, ivar_new, gpm_new, log10_blaze_new
        interpolated flux, inverse variance, and good-pixel mask arrays with
        the length matching the new wavelength grid. They are all
        `numpy.ndarray`_ objects with except if log10_blaze_function is None in which case log10_blaze_new is None

    """
    # Check input
    if wave_new.ndim != 1 or wave_old.ndim != 1:
        msgs.error('All input vectors must be 1D.')
    if flux_old.shape != wave_old.shape or ivar_old.shape != wave_old.shape \
            or gpm_old.shape != wave_old.shape:
        msgs.error('All vectors to interpolate must have the same size.')

    # Do not interpolate if the wavelength is exactly same with wave_new
    if np.array_equal(wave_new, wave_old) and not sensfunc:
        return flux_old, ivar_old, gpm_old

    wave_gpm = wave_old > 1.0 # Deal with the zero wavelengths
    if sensfunc:
        delta_wave_interp = wvutils.get_delta_wave(wave_old, wave_gpm)
        flux_interp = flux_old[wave_gpm]/delta_wave_interp[wave_gpm]
        ivar_interp = ivar_old[wave_gpm]*delta_wave_interp[wave_gpm]**2
        if log10_blaze_function is not None:
            log10_blaze_interp = log10_blaze_function[wave_gpm] - np.log10(delta_wave_interp[wave_gpm])
    else:
        flux_interp = flux_old[wave_gpm]
        ivar_interp = ivar_old[wave_gpm]
        if log10_blaze_function is not None:
            log10_blaze_interp = log10_blaze_function[wave_gpm]

    flux_new = scipy.interpolate.interp1d(wave_old[wave_gpm], flux_interp, kind=kind,
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    ivar_new = scipy.interpolate.interp1d(wave_old[wave_gpm], ivar_interp, kind=kind,
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    log10_blaze_new = scipy.interpolate.interp1d(wave_old[wave_gpm], log10_blaze_interp, kind=kind,
                                                 bounds_error=False, fill_value=np.nan)(wave_new) \
        if log10_blaze_function is not None else None

    # Interpolate a floating-point version of the mask. Use linear interpolation here
    gpm_new_tmp = scipy.interpolate.interp1d(wave_old[wave_gpm], gpm_old.astype(float)[wave_gpm],
                                             kind='linear', bounds_error=False,
                                             fill_value=np.nan)(wave_new)
    # Don't allow the ivar to be ever be less than zero
    ivar_new = (ivar_new > 0.0)*ivar_new
    gpm_new = (gpm_new_tmp > 0.8) & (ivar_new > 0.0) & np.isfinite(flux_new) & np.isfinite(ivar_new)
    return flux_new, ivar_new, gpm_new, log10_blaze_new



# TODO: ``sensfunc`` should be something like "conserve_flux". It would be
# useful to compare these resampling routines against
# `pypeit.sampling.Resample`.
def interp_spec(wave_new, waves, fluxes, ivars, gpms, log10_blaze_function=None, sensfunc=False, kind='cubic'):
    """
    Interpolate a set of spectra onto a new wavelength grid.

    The method can perform two types of interpolation, depending on the
    shapes of the input arrays.

        1. If the new wavelength grid (``wave_new``) is 1D, all input spectra
           are interpolated to this new grid. The input spectra can be
           provided as either 1D or 2D arrays.

        2. If the new wavelength grid (``wave_new``) is 2D, all input spectra
           *must* be 1D. The single spectrum is then interpolated onto each of
           the new wavelength grids.

    Parameters
    ----------
    wave_new : `numpy.ndarray`_, shape (nspec,) or (nspec, nimgs),
        New wavelength grid for output spectra. Shape can be 1D or 2D.  See the
        method description for how this affects the code flow above.
    waves : `numpy.ndarray`_, shape (nspec,) or (nspec, nexp)
        Wavelength vector for current spectra. Shape can be 1D or 2D, where
        nexp, need not equal nimgs.  See the method description for how this
        affects the code flow above.
    fluxes : `numpy.ndarray`_
        Flux vectors.  Shape must match ``waves``.
    ivars : `numpy.ndarray`_
        Inverse variance vectors.  Shape must match ``waves``.
    gpms : `numpy.ndarray`_
        Boolean good-pixel masks for each spectrum (True=Good). Shape must match
        ``waves``.
    sensfunc : :obj:`bool`, optional
        If True, the quantities ``flux*delta_wave`` and the corresponding
        ``ivar/delta_wave**2`` will be interpolated and returned instead of
        ``flux`` and ``ivar``. This is useful for sensitivity function
        computation where we need flux*(wavelength bin width). Beacause
        delta_wave is a difference of the wavelength grid, interpolating in the
        presence of masked data requires special care.
    kind : str or int, optional
        Specifies the kind of interpolation as a string or as an integer
        specifying the order of the spline interpolator to use following the convention of
        scipy.interpolate.interp1d. The string has to be one of 'linear', 'nearest',
        'nearest-up', 'zero', 'slinear', 'quadratic', 'cubic', 'previous', or 'next'. 'zero',
        'slinear', 'quadratic' and 'cubic' refer to a spline interpolation of
        zeroth, first, second or third order; 'previous' and 'next' simply
        return the previous or next value of the point; 'nearest-up' and
        'nearest' differ when interpolating half-integers (e.g. 0.5, 1.5)
        in that 'nearest-up' rounds up and 'nearest' rounds down. Default is 'cubic'.

    Returns
    -------
    fluxes_inter : `numpy.ndarray`_
        interpolated flux with size and shape matching the new wavelength grid.
    ivars_inter : `numpy.ndarray`_
        interpolated inverse variance with size and shape matching the new
        wavelength grid.
    gpms_inter : `numpy.ndarray`_
        interpolated good-pixel mask with size and shape matching the new
        wavelength grid.
    log10_blazes_inter : `numpy.ndarray`_ or None
        interpolated log10 blaze function with size and shape matching the new
        wavelength grid.  If ``log10_blaze_function`` is not provided as an input,
        this will be None.
    """
    # Check input
    if wave_new.ndim > 2:
        msgs.error('Invalid shape for wave_new; must be 1D or 2D')
    if wave_new.ndim == 2 and fluxes.ndim != 1:
        msgs.error('If new wavelength grid is 2D, all other input arrays must be 1D.')
    if fluxes.shape != waves.shape or ivars.shape != waves.shape or gpms.shape != waves.shape:
        msgs.error('Input spectral arrays must all have the same shape.')

    # First case: interpolate either an (nspec, nexp) array of spectra onto a
    # single wavelength grid
    if wave_new.ndim == 1:
        if fluxes.ndim == 1:
            return interp_oned(wave_new, waves, fluxes, ivars, gpms, log10_blaze_function=log10_blaze_function,
                               sensfunc=sensfunc, kind=kind)

        nexp = fluxes.shape[1]
        # Interpolate spectra to have the same wave grid with the iexp spectrum.
        # And scale spectra to the same flux level with the iexp spectrum.
        fluxes_inter = np.zeros((wave_new.size, nexp), dtype=float)
        ivars_inter = np.zeros((wave_new.size, nexp), dtype=float)
        gpms_inter = np.zeros((wave_new.size, nexp), dtype=bool)
        if log10_blaze_function is not None:
            log10_blazes_inter = np.zeros((wave_new.size, nexp), dtype=float)
            for ii in range(nexp):
                fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii], log10_blazes_inter[:,ii] \
                    = interp_oned(wave_new, waves[:,ii], fluxes[:,ii], ivars[:,ii], gpms[:,ii],
                                  log10_blaze_function = log10_blaze_function[:, ii], sensfunc=sensfunc, kind=kind)
        else:
            for ii in range(nexp):
                fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii], _ \
                    = interp_oned(wave_new, waves[:,ii], fluxes[:,ii], ivars[:,ii], gpms[:,ii], sensfunc=sensfunc, kind=kind)
            log10_blazes_inter=None

        return fluxes_inter, ivars_inter, gpms_inter, log10_blazes_inter


    # Second case: interpolate a single spectrum onto an (nspec, nexp) array of
    # wavelengths. To make it here, wave_new.ndim must be 2.
    fluxes_inter = np.zeros_like(wave_new, dtype=float)
    ivars_inter = np.zeros_like(wave_new, dtype=float)
    gpms_inter = np.zeros_like(wave_new, dtype=bool)
    if log10_blaze_function is not None:
        log10_blazes_inter = np.zeros_like(wave_new, dtype=float)
        for ii in range(wave_new.shape[1]):
            fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii], log10_blazes_inter[:, ii] \
                = interp_oned(wave_new[:,ii], waves, fluxes, ivars, gpms, log10_blaze_function=log10_blaze_function, sensfunc=sensfunc)
    else:
        for ii in range(wave_new.shape[1]):
            fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii], _ \
                = interp_oned(wave_new[:,ii], waves, fluxes, ivars, gpms, log10_blaze_function=log10_blaze_function, sensfunc=sensfunc)
        log10_blazes_inter=None

    return fluxes_inter, ivars_inter, gpms_inter, log10_blazes_inter


def smooth_weights(inarr, gdmsk, sn_smooth_npix):
    """Smooth the input weights with a Gaussian 1D kernel.

    Args:
        inarr (`numpy.ndarray`_):
            S/N spectrum to be smoothed. shape = (nspec,)
        gdmsk (`numpy.ndarray`_):
            Boolean mask of good pixels. shape = (nspec,)
        sn_smooth_npix (float):
            Number of pixels used for determining smoothly varying S/N ratio weights.
            The sigma of the kernel is set by
            sig_res = max(sn_smooth_npix / 10.0, 3.0)

    Returns:
        `numpy.ndarray`_: smoothed version of inarr.
    """
    spec_vec = np.arange(gdmsk.size)
    sn_med1 = np.zeros(inarr.size)
    sn_med1[gdmsk] = utils.fast_running_median(inarr[gdmsk], sn_smooth_npix)
    sn_med2 = scipy.interpolate.interp1d(spec_vec[gdmsk], sn_med1[gdmsk], kind='cubic',
                                         bounds_error=False, fill_value=-999)(spec_vec)
    # Fill the S/N weight to the left and right with the nearest value
    mask_good = np.where(sn_med2 != -999)[0]
    idx_mn, idx_mx = np.min(mask_good), np.max(mask_good)
    sn_med2[:idx_mn] = sn_med2[idx_mn]
    sn_med2[idx_mx:] = sn_med2[idx_mx]
    # Smooth with a Gaussian kernel
    sig_res = np.fmax(sn_smooth_npix / 10.0, 3.0)
    gauss_kernel = convolution.Gaussian1DKernel(sig_res)
    sn_conv = convolution.convolve(sn_med2, gauss_kernel, boundary='extend')
    # TODO Should we be setting a floor on the weights to prevent tiny numbers?
    return sn_conv

# TODO add a wave_min, wave_max option here to deal with objects like high-z QSOs etc. which only have flux in a given wavelength range.
def calc_snr(fluxes, ivars, gpms):

    """
    Calculate the rms S/N of each input spectrum.

    Parameters
    ----------
    fluxes : list
        List of length nexp containing the `numpy.ndarray`_ 1d float spectra. The
        shapes in the list can be different. nexp = len(fluxes)
    ivars : list
        List of length nexp containing the `numpy.ndarray`_ 1d float inverse
        variances of the spectra.
    gpms : list
        List of length nexp containing the `numpy.ndarray`_ 1d float boolean good
        pixel masks of the spectra.
    verbose : bool, optional
        Verbosity of print out.

    Returns
    -------
    rms_sn : `numpy.ndarray`_
        Array of shape (nexp,) of root-mean-square S/N value for each input spectra where nexp=len(fluxes).
    sn_val : list
        List of length nexp containing the wavelength dependent S/N arrays for each input spectrum, i.e.
        each element contains the array flux*sqrt(ivar)
    """

    nexp = len(fluxes)

    # Calculate S/N
    sn_val, rms_sn, sn2 = [], [], []
    for iexp, (flux, ivar, gpm) in enumerate(zip(fluxes, ivars, gpms)):
        sn_val_iexp = flux*np.sqrt(ivar)
        sn_val.append(sn_val_iexp)
        sn_val_ma = np.ma.array(sn_val_iexp, mask=np.logical_not(gpm))
        sn_sigclip = stats.sigma_clip(sn_val_ma, sigma=3, maxiters=5)
        sn2_iexp = sn_sigclip.mean()**2  # S/N^2 value for each spectrum
        if sn2_iexp is np.ma.masked:
            msgs.error(f'No unmasked value in iexp={iexp+1}/{nexp}. Check inputs.')
        else:
            sn2.append(sn2_iexp)
            rms_sn.append(np.sqrt(sn2_iexp))  # Root Mean S/N**2 value for all spectra

    return np.array(rms_sn), sn_val

# TODO add a wave_min, wave_max option here to deal with objects like high-z QSOs etc. which only have flux in a given wavelength range.
def sn_weights(fluxes, ivars, gpms, sn_smooth_npix=None, weight_method='auto', verbose=False):

    """
    Calculate the S/N of each input spectrum and create an array of
    (S/N)^2 weights to be used for coadding.

    Parameters
    ----------
    fluxes : list
        List of len(nexp) containing the `numpy.ndarray`_ 1d float spectra. The
        shapes in the list can be different.
    ivars : list
        List of len(nexp) containing the `numpy.ndarray`_ 1d float inverse
        variances of the spectra.
    gpms : list
        List of len(nexp) containing the `numpy.ndarray`_ 1d float boolean good
        pixel masks of the spectra.
    sn_smooth_npix : float, optional
        Number of pixels used for determining smoothly varying S/N ratio
        weights. This must be passed for all weight methods excpet for weight_method='constant' or 'uniform', since then
        wavelength dependent weights are not used.
    weight_method : str, optional

        The weighting method to be used. Options are ``'auto'``, ``'constant'``, ``'uniform'``, ``'wave_dependent'``,
        ``'relative'``, or ``'ivar'``. The default is ``'auto'``.  Behavior is
        as follows:

            - ``'auto'``: Use constant weights if rms_sn < 3.0, otherwise use
              wavelength dependent.

            - ``'constant'``: Constant weights based on rms_sn**2

            - ``'uniform'``: Uniform weighting.

            - ``'wave_dependent'``: Wavelength dependent weights will be used
              irrespective of the rms_sn ratio. This option will not work well
              at low S/N ratio although it is useful for objects where only a
              small fraction of the spectral coverage has high S/N ratio (like
              high-z quasars).

            - ``'relative'``: Calculate weights by fitting to the ratio of
              spectra? Note, relative weighting will only work well when there
              is at least one spectrum with a reasonable S/N, and a continuum.
              RJC note - This argument may only be better when the object being
              used has a strong continuum + emission lines. The reference
              spectrum is assigned a value of 1 for all wavelengths, and the
              weights of all other spectra will be determined relative to the
              reference spectrum.  This is particularly useful if you are
              dealing with highly variable spectra (e.g. emission lines) and
              require a precision better than ~1 per cent.

            - ``'ivar'``: Use inverse variance weighting. This is not well
              tested and should probably be deprecated.

    verbose : bool, optional
        Verbosity of print out.

    Returns
    -------
    rms_sn : `numpy.ndarray`_
        Array of root-mean-square S/N value for each input spectra. Shape = (nexp,)
    weights : list
        List of  len(nexp) containing the signal-to-noise squared weights to be
        applied to the spectra. This output is aligned with the vector (or
        vectors) provided in waves, i.e. it is a list of arrays of type
        `numpy.ndarray`_  with the same shape as those in waves.
    """
    if weight_method not in ['auto', 'constant', 'uniform', 'wave_dependent', 'relative', 'ivar']:
        msgs.error('Unrecognized option for weight_method=%s').format(weight_method)

    nexp = len(fluxes)
    # Check that all the input lists have the same length
    if len(ivars) != nexp or len(gpms) != nexp:
        msgs.error("Input lists of spectra must have the same length")

    # Check that sn_smooth_npix if weight_method = constant or uniform
    if sn_smooth_npix is None and weight_method not in ['constant', 'uniform']:
        msgs.error("sn_smooth_npix cannot be None unless the weight_method='constant' or weight_method='uniform'")

    rms_sn, sn_val = calc_snr(fluxes, ivars, gpms)
    sn2 = np.square(rms_sn)

    # Check if relative weights input
    if verbose:
        msgs.info('Computing weights with weight_method=%s'.format(weight_method))

    weights = []

    weight_method_used = [] if weight_method == 'auto' else nexp*[weight_method]
    if weight_method == 'relative':
        # Relative weights are requested, use the highest S/N spectrum as a reference
        ref_spec = np.argmax(sn2)
        if verbose:
            msgs.info(
                "The reference spectrum (ref_spec={0:d}) has a typical S/N = {1:.3f}".format(ref_spec, sn2[ref_spec]))
        # Adjust the arrays to be relative
        refscale = utils.inverse(sn_val[ref_spec])
        for iexp in range(nexp):
            # Compute the relative (S/N)^2 and update the mask
            sn2[iexp] /= sn2[ref_spec]
            gpm_rel = gpms[iexp] & ((gpms[ref_spec]) | (sn_val[ref_spec] != 0))
            sn_val_rescaled = sn_val[iexp]*refscale
            weights.append(smooth_weights(sn_val_rescaled** 2, gpm_rel, sn_smooth_npix))
    elif weight_method == 'ivar':
        # TODO: Should ivar weights be deprecated??
        for ivar, mask in zip(ivars, gpms):
            weights.append(smooth_weights(ivar, mask, sn_smooth_npix))
    elif weight_method == 'constant':
        for iexp in range(nexp):
            weights.append(np.full(fluxes[iexp].size, np.fmax(sn2[iexp], 1e-2)))  # set the minimum  to be 1e-2 to avoid zeros
    elif weight_method == 'uniform':
        for iexp in range(nexp):
            weights.append(
                np.full(fluxes[iexp].size, 1.0))
    elif weight_method == 'wave_dependent':
        for iexp in range(nexp):
            weights.append(smooth_weights(sn_val[iexp] ** 2, gpms[iexp], sn_smooth_npix))
    elif weight_method == 'auto':
        for iexp in range(nexp):
            if rms_sn[iexp] < 3.0:
                weights.append(np.full(fluxes[iexp].size, np.fmax(sn2[iexp], 1e-2)))  # set the minimum  to be 1e-2 to avoid zeros
                weight_method_used.append('constant')
            else:
                weights.append(smooth_weights(sn_val[iexp]**2, gpms[iexp], sn_smooth_npix))
                weight_method_used.append('wavelength dependent')

    if verbose:
        for iexp in range(nexp):
            msgs.info('Using {:s} weights for coadding, S/N '.format(weight_method_used[iexp]) +
                      '= {:4.2f}, weight = {:4.2f} for {:}th exposure'.format(rms_sn[iexp], np.mean(weights[iexp]), iexp))

    # Finish
    return np.array(rms_sn), weights


def robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None, ref_percentile=70.0, min_good=0.05,
                        maxiters=5, sigrej=3.0, max_factor=10.0, snr_do_not_rescale=1.0,
                        verbose=False):
    """
    Robustly determine the ratio between input spectrum flux and reference spectrum flux_ref. The code will perform
    best if the reference spectrum is chosen to be the higher S/N ratio spectrum, i.e. a preliminary stack that you want
    to scale each exposure to match. Note that the flux and flux_ref need to be on the same wavelength grid!!

    Parameters
    ----------
    flux : `numpy.ndarray`_
        spectrum that will be rescaled. shape=(nspec,)
    ivar : `numpy.ndarray`_
        inverse variance for the spectrum that will be rescaled.  Same shape as
        flux
    flux_ref : `numpy.ndarray`_
        reference spectrum. Same shape as flux
    mask : `numpy.ndarray`_, optional
        boolean mask for the spectrum that will be rescaled. True=Good.  If not
        input, computed from inverse variance
    ivar_ref : `numpy.ndarray`_, optional
        inverse variance of reference spectrum.
    mask_ref : `numpy.ndarray`_, optional
        Boolean mask for reference spectrum. True=Good. If not input, computed
        from inverse variance.
    ref_percentile : float, optional, default=70.0
        Percentile fraction used for selecting the minimum SNR cut from the
        reference spectrum. Pixels above this percentile cut are deemed the
        "good" pixels and are used to compute the ratio. This must be a number
        between 0 and 100.
    min_good : float, optional, default = 0.05
        Minimum fraction of good pixels determined as a fraction of the total
        pixels for estimating the median ratio
    maxiters : int, optional, default = 5
        Maximum number of iterations for astropy.stats.SigmaClip
    sigrej : float, optional, default = 3.0
        Rejection threshold for astropy.stats.SigmaClip
    max_factor : float, optional, default = 10.0,
        Maximum allowed value of the returned ratio
    snr_do_not_rescale : float, optional, default = 1.0
        If the S/N ratio of the set of pixels (defined by upper ref_percentile
        in the reference spectrum) in the input spectrum have a median value
        below snr_do_not_rescale, median rescaling will not be attempted and the
        code returns ratio = 1.0. We also use this parameter to define the set
        of pixels (determined from the reference spectrum) to compare for the
        rescaling.

    Returns
    -------
    ratio : float
        the number that must be multiplied into flux in order to get it to match
        up with flux_ref
    """

    ## Mask for reference spectrum and your spectrum
    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0

    nspec = flux.size
    snr_ref = flux_ref * np.sqrt(ivar_ref)
    
    snr_ref_best = np.fmax(np.percentile(snr_ref[mask_ref], ref_percentile),snr_do_not_rescale)
    calc_mask = (snr_ref > snr_ref_best) & mask_ref & mask

    snr_resc = flux*np.sqrt(ivar)
    snr_resc_med = np.median(snr_resc[calc_mask])

    if (np.sum(calc_mask) > min_good*nspec) & (snr_resc_med > snr_do_not_rescale):
        # Take the best part of the higher SNR reference spectrum
        sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median', stdfunc=utils.nan_mad_std)

        flux_ref_ma = np.ma.MaskedArray(flux_ref, np.logical_not(calc_mask))
        flux_ref_clipped, lower, upper = sigclip(flux_ref_ma, masked=True, return_bounds=True)
        mask_ref_clipped = np.logical_not(flux_ref_clipped.mask)  # mask_stack = True are good values

        flux_ma = np.ma.MaskedArray(flux_ref, np.logical_not(calc_mask))
        flux_clipped, lower, upper = sigclip(flux_ma, masked=True, return_bounds=True)
        mask_clipped = np.logical_not(flux_clipped.mask)  # mask_stack = True are good values

        new_mask = mask_ref_clipped & mask_clipped

        flux_ref_median = np.median(flux_ref[new_mask])
        flux_dat_median = np.median(flux[new_mask])

        if (flux_ref_median < 0.0) or (flux_dat_median < 0.0):
            msgs.warn('Negative median flux found. Not rescaling')
            ratio = 1.0
        else:
            if verbose:
                msgs.info('Used {:} good pixels for computing median flux ratio'.format(np.sum(new_mask)))
            ratio = np.fmax(np.fmin(flux_ref_median/flux_dat_median, max_factor), 1.0/max_factor)
    else:
        if (np.sum(calc_mask) <= min_good*nspec):
            msgs.warn('Found only {:} good pixels for computing median flux ratio.'.format(np.sum(calc_mask))
            + msgs.newline() + 'No median rescaling applied')
        if (snr_resc_med <= snr_do_not_rescale):
            msgs.warn('Median flux ratio of pixels in reference spectrum {:} <= snr_do_not_rescale = {:}.'.format(snr_resc_med, snr_do_not_rescale)
                      + msgs.newline() + 'No median rescaling applied')
        ratio = 1.0

    return ratio


def scale_spec(wave, flux, ivar, sn, wave_ref, flux_ref, ivar_ref, mask=None, mask_ref=None, scale_method='auto', min_good=0.05,
               ref_percentile=70.0, maxiters=5, sigrej=3, max_median_factor=10.0,
               npoly=None, hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5, debug=False, show=False):
    """
    Routine for solving for the best way to rescale an input spectrum
    flux to match a reference spectrum flux_ref. The code will work
    best if you choose the reference to be the highest S/N ratio
    spectrum. If the scale_method is not specified, the code will
    make a decision about which method to use based on the input S/N
    ratio.

    Parameters
    ----------
    wave: `numpy.ndarray`_
            wavelengths grid for the spectra of shape (nspec,)
    flux: `numpy.ndarray`_
            spectrum that will be rescaled.
    ivar: `numpy.ndarray`_
            inverse variance for the spectrum that will be rescaled.
    sn: float
        S/N of the spectrum that is being scaled used to make decisions about the scaling method.
        This can be computed by sn_weights and passed in.
    flux_ref: `numpy.ndarray`_, (nspec,)
            reference spectrum.
    ivar_ref: `numpy.ndarray`_, (nspec,)
            inverse variance of reference spectrum.
    mask: `numpy.ndarray`_
            Boolean mask for the spectrum that will be rescaled. True=Good. If not input, computed from inverse variance
    mask_ref: `numpy.ndarray`_
            Boolean mask for reference spectrum. True=Good. If not input, computed from inverse variance.
    min_good: float, optional, default = 0.05
            minimum fraction of the total number of good pixels needed for estimate the median ratio
    maxiters: int, optional
            maximum number of iterations for rejecting outliers used
            by the robust_median_ratio routine if median rescaling is
            the method used.
    max_median_factor: float, optional, default=10.0
            maximum scale factor for median rescaling for robust_median_ratio if median rescaling is the method used.
    sigrej: float, optional, default=3.0
            rejection threshold used for rejecting outliers by robsut_median_ratio
    ref_percentile: float, optional, default=70.0
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
    npoly: int, optional, default=None
            order for the poly ratio scaling if polynomial rescaling
            is the method used. Default is to automatically compute
            this based on S/N ratio of data.
    scale_method: str, optional
            scale method, str, default='auto'. Options are auto,
            poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default
            is to let the code determine this automitically which
            works well.
    hand_scale: `numpy.ndarray`_, optional
            array of hand scale factors, not well tested. shape=(nexp,)
    sn_min_polyscale: float, optional, default=2.0
            maximum SNR for perforing median scaling
    sn_min_medscale: float, optional, default=0.5
            minimum SNR for perforing median scaling
    debug: bool, optional, default=False
            show interactive QA plot

    Returns
    -------
    flux_scale : `numpy.ndarray`_, shape is (nspec,)
        Scaled spectrum
    ivar_scale : `numpy.ndarray`_, shape is (nspec,)
        Inverse variance of scaled spectrum 
    scale : `numpy.ndarray`_, shape is (nspec,)
        Scale factor applied to the spectrum and inverse variance
    method_used : :obj:`str`
        Method that was used to scale the spectra.
    """

    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0


    flux_ref_int, ivar_ref_int, mask_ref_int, _ = interp_spec(wave, wave_ref, flux_ref, ivar_ref, mask_ref)

    # estimates the SNR of each spectrum and the stacked mean SNR
    #rms_sn, weights = sn_weights(wave, flux, ivar, mask, sn_smooth_npix)
    #sn = np.sqrt(np.mean(rms_sn**2))

    if scale_method == 'auto':
        if sn > sn_min_polyscale:
            method_used = 'poly'
        elif ((sn <= sn_min_polyscale) and (sn > sn_min_medscale)):
            method_used = 'median'
        else:
            method_used = 'none'
    else:
        method_used = scale_method

    # Estimate the scale factor
    if method_used == 'poly':
        # Decide on the order of the polynomial rescaling
        if npoly is None:
            if sn > 25.0:
                npoly = 5 # quintic, Is this stable?
            elif sn > 8.0:
                npoly = 3  # cubic
            elif sn >= 5.0:
                npoly = 2  # quadratic
            else:
                npoly = 1  # linear
        scale, fit_tuple, flux_scale, ivar_scale, outmask = solve_poly_ratio(
            wave, flux, ivar, flux_ref_int, ivar_ref_int, npoly,mask=mask, mask_ref=mask_ref_int,
            ref_percentile=ref_percentile, debug=debug)
    elif method_used == 'median':
        # Median ratio (reference to spectrum)
        med_scale = robust_median_ratio(flux, ivar, flux_ref_int, ivar_ref_int,ref_percentile=ref_percentile,min_good=min_good,
                                        mask=mask, mask_ref=mask_ref_int, maxiters=maxiters,
                                        max_factor=max_median_factor,sigrej=sigrej)
        # Apply
        flux_scale = flux * med_scale
        ivar_scale = ivar * 1.0/med_scale**2
        scale = np.full_like(flux,med_scale)
    elif method_used == 'hand':
        # Input?
        if hand_scale is None:
            msgs.error("Need to provide hand_scale parameter, single value")
        flux_scale = flux * hand_scale
        ivar_scale = ivar * 1.0 / hand_scale ** 2
        scale = np.full(flux.size, hand_scale)
    elif method_used == 'none':
        flux_scale = flux.copy()
        ivar_scale = ivar.copy()
        scale = np.ones_like(flux)
    else:
        msgs.error("Scale method not recognized! Check documentation for available options")
    # Finish
    if show:
        scale_spec_qa(wave, flux*mask, ivar*mask, wave_ref, flux_ref*mask_ref, ivar_ref*mask_ref, scale, method_used, mask = mask, mask_ref=mask_ref,
                      title='Scaling Applied to the Data')

    return flux_scale, ivar_scale, scale, method_used

def compute_stack(wave_grid, waves, fluxes, ivars, gpms, weights, min_weight=1e-8):
    """
    Compute a stacked spectrum from a set of exposures on the specified
    wave_grid with proper treatment of weights and masking. This code uses
    np.histogram to combine the data using NGP and does not perform any
    interpolations and thus does not correlate errors. It uses wave_grid to
    determine the set of wavelength bins that the data are averaged on. The
    final spectrum will be on an ouptut wavelength grid which is not the same as
    wave_grid.  The ouput wavelength grid is the weighted average of the
    individual wavelengths used for each exposure that fell into a given
    wavelength bin in the input wave_grid. This 1d coadding routine thus
    maintains the independence of the errors for each pixel in the combined
    spectrum and computes the weighted averaged wavelengths of each pixel in an
    analogous way to the 2d extraction procedure which also never interpolates
    to avoid correlating erorrs.

    Parameters
    ----------
    wave_grid : `numpy.ndarray`_
        new wavelength grid desired. This will typically be a reguarly spaced
        grid created by the get_wave_grid routine.  The reason for the ngrid+1
        is that this is the general way to specify a set of  bins if you desire
        ngrid bin centers, i.e. the output stacked spectra have ngrid elements.
        The spacing of this grid can be regular in lambda (better for multislit)
        or log lambda (better for echelle). This new wavelength grid should be
        designed with the sampling of the data in mind. For example, the code
        will work fine if you choose the sampling to be too fine, but then the
        number of exposures contributing to any given wavelength bin will be one
        or zero in the limiting case of very small wavelength bins. For larger
        wavelength bins, the number of exposures contributing to a given bin
        will be larger.  shape=(ngrid +1,)
    waves : list
        List of length nexp `numpy.ndarray`_ float wavelength arrays for spectra
        to be stacked.  Note that the wavelength grids can in general be
        different for each exposure, irregularly spaced, and can have different
        sizes.
    fluxes : list
        List of length nexp `numpy.ndarray`_ float flux arrays for spectra to be
        stacked. These are aligned to the wavelength grids in waves.
    ivars : list
        List of length nexp `numpy.ndarray`_ float inverse variance arrays for
        spectra to be stacked. These are aligned to the wavelength grids in
        waves.
    gpms : list
        List of length nexp `numpy.ndarray`_ boolean good pixel mask arrays for
        spectra to be stacked. These are aligned to the wavelength grids in
        waves. True=Good.
    weights : list
        List of length nexp `numpy.ndarray`_ float weights to be used for
        combining the spectra. These are aligned to the wavelength grids in
        waves and were computed using sn_weights.
    min_weight : float, optional
        Minimum allowed weight for any individual spectrum

    Returns
    -------
    wave_stack : `numpy.ndarray`_
        Wavelength grid for stacked spectrum. As discussed above, this is the
        weighted average of the wavelengths of each spectrum that contriuted to
        a bin in the input wave_grid wavelength grid. It thus has ngrid
        elements, whereas wave_grid has ngrid+1 elements to specify the ngrid
        total number of bins. Note that wave_stack is NOT simply the wave_grid
        bin centers, since it computes the weighted average. shape=(ngrid,)
    flux_stack : `numpy.ndarray`_, shape=(ngrid,)
        Final stacked spectrum on wave_stack wavelength grid
    ivar_stack : `numpy.ndarray`_, shape=(ngrid,)
        Inverse variance spectrum on wave_stack wavelength grid.
        Errors are propagated according to weighting and masking.
    gpm_stack : `numpy.ndarray`_, shape=(ngrid,)
        Boolean Mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
    nused : `numpy.ndarray`_, shape=(ngrid,)
        Number of exposures which contributed to
        each pixel in the wave_stack. Note that this is in general
        different from nexp because of masking, but also becuse of
        the sampling specified by wave_grid. In other words,
        sometimes more spectral pixels in the irregularly gridded
        input wavelength array waves will land in one bin versus
        another depending on the sampling.
    """

    #mask bad values and extreme values (usually caused by extreme low sensitivity at the edge of detectors)
    #TODO cutting on the value of ivar is dicey for data in different units. This should be removed.
    uber_gpms = [gpm & (weight > 0.0) & (wave > 1.0) & (ivar > 0.0) & (utils.inverse(ivar)<1e10)
                 for gpm, weight, wave, ivar in zip(gpms, weights, waves, ivars)]
    waves_flat, fluxes_flat, ivars_flat, weights_flat = [], [], [], []
    for wave, flux, ivar, weight, ugpm in zip(waves, fluxes, ivars, weights, uber_gpms):
        waves_flat   += wave[ugpm].tolist()
        fluxes_flat  += flux[ugpm].tolist()
        ivars_flat   += ivar[ugpm].tolist()
        weights_flat += weight[ugpm].tolist()

    waves_flat = np.array(waves_flat)
    fluxes_flat = np.array(fluxes_flat)
    ivars_flat = np.array(ivars_flat)
    weights_flat = np.array(weights_flat)
    vars_flat = utils.inverse(np.array(ivars_flat))

    # Counts how many pixels in each wavelength bin
    nused, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False)

    # Calculate the summed weights for the denominator
    weights_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=weights_flat)

    # Calculate the stacked wavelength
    ## TODO: JFH Made the minimum weight 1e-8 from 1e-4. I'm not sure what this min_weight is necessary for, or
    # is achieving FW.
    wave_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=waves_flat*weights_flat)
    wave_stack = (weights_total > min_weight)*wave_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked flux
    flux_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=fluxes_flat*weights_flat)
    flux_stack = (weights_total > min_weight)*flux_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked ivar
    var_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=vars_flat*weights_flat**2)
    var_stack = (weights_total > min_weight)*var_stack_total/(weights_total+(weights_total==0.))**2
    ivar_stack = utils.inverse(var_stack)

    # New mask for the stack
    gpm_stack = (weights_total > min_weight) & (nused > 0.0)
    return wave_stack, flux_stack, ivar_stack, gpm_stack, nused

def get_ylim(flux, ivar, mask):
    """
    Utility routine for setting the plot limits for QA plots.

    Args:
        flux (`numpy.ndarray`_):
            (nspec,) flux array
        ivar (`numpy.ndarray`_):
            (nspec,) inverse variance array
        mask (`numpy.ndarray`_):
            bool, (nspec,) mask array. True=Good

    Returns:
        tuple: lower and upper limits for plotting.

    """

    med_width = (2.0 * np.ceil(0.1 / 2.0 * np.size(flux[mask])) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    mask_lim = ivar_med > np.percentile(ivar_med, 20)
    ymax = 2.5 * np.max(flux_med[mask_lim])
    ymin = -0.15 * ymax
    return ymin, ymax

def scale_spec_qa(wave, flux, ivar, wave_ref, flux_ref, ivar_ref, ymult,
                  scale_method, mask=None, mask_ref=None, ylim = None, title=''):
    '''
    QA plot for spectrum scaling.

    Parameters
    ----------
    wave: `numpy.ndarray`_
            wavelength array for spectrum to be scaled and reference spectrum.
            shape=(nspec,)
    flux: `numpy.ndarray`_
            flux for spectrum to be scaled; shape=(nspec,)
    ivar: `numpy.ndarray`_
            inverse variance for spectrum to be scaled.  shape=(nspec,)
    wave_ref: `numpy.ndarray`_
            reference wavelengths; shape=(nspec,)
    flux_ref: `numpy.ndarray`_
            reference flux; shape=(nspec,)
    ivar_ref: `numpy.ndarray`_
            inverse variance of reference flux; shape=(nspec,)
    ymult: `numpy.ndarray`_
            scale factor array; shape=(nspec,)
    scale_method: str
            label of method used for rescaling which will be shown on QA plot.
    mask: `numpy.ndarray`_, optional
            Boolean mask for spectrum to be scaled. True=Good. If not specified determined form inverse variance
            shape=(nspec,)
    mask_ref: `numpy.ndarray`_, optional
            Boolean mask for reference flux. True=Good.
            shape=(nspec,)
    ylim: tuple, optional
            tuple for limits of the QA plot. If None, will be determined automtically with get_ylim
    title: str, optional
        QA plot title
    '''

    if mask is None:
        mask = ivar > 0.0
    if mask_ref  is None:
        mask_ref = ivar_ref > 0.0

    # This deals with spectrographs that have zero wavelength values. They are masked in mask, but this impacts plotting
    wave_mask = wave > 1.0
    wave_mask_ref = wave_ref > 1.0
    #dwave = wave[wave_mask].max() - wave[wave_mask].min()
    #dwave_ref = wave_ref[wave_mask_ref].max() - wave_ref[wave_mask_ref].min()
    # Get limits
    if ylim is None:
        ylim = get_ylim(flux_ref, ivar_ref, mask_ref)

    nullfmt = NullFormatter()  # no labels
    fig = plt.figure(figsize=(12, 8))
    # [left, bottom, width, height]
    poly_plot = fig.add_axes([0.1, 0.75, 0.8, 0.20])
    spec_plot = fig.add_axes([0.1, 0.10, 0.8, 0.65])
    poly_plot.xaxis.set_major_formatter(nullfmt)  # no x-axis labels for polynomial plot
    poly_plot.plot(wave[wave_mask], ymult[wave_mask], color='black', linewidth=3.0, label=scale_method + ' scaling')
    poly_plot.legend()
    # This logic below allows more of the spectrum to be plotted if wave_ref is a multi-order stack which has broader
    # wavelength coverage. For the longslit or single order case, this will plot the correct range as well
    wave_min = np.fmax(0.8*wave[wave_mask].min(), wave_ref[wave_mask_ref].min())
    wave_max = np.fmin(1.2*wave[wave_mask].max(), wave_ref[wave_mask_ref].max())
    poly_plot.set_xlim((wave_min, wave_max))
    spec_plot.set_xlim((wave_min, wave_max))
    spec_plot.set_ylim(ylim)

    spec_plot.plot(wave[wave_mask], flux[wave_mask], color='red', zorder=10,
                   marker='o', markersize=1.0, mfc='k', fillstyle='full', linestyle='None', label='original spectrum')
    spec_plot.plot(wave[wave_mask], flux[wave_mask]*ymult[wave_mask], color='dodgerblue', drawstyle='steps-mid', alpha=0.5, zorder=5, linewidth=2,
                   label='rescaled spectrum')
    spec_plot.plot(wave_ref[wave_mask_ref], flux_ref[wave_mask_ref], color='black', drawstyle='steps-mid', zorder=7, alpha = 0.5, label='reference spectrum')

    spec_plot.legend()
    fig.suptitle(title)
    plt.show()

# TODO: Change mask to gpm
def coadd_iexp_qa(wave, flux, rejivar, mask, wave_stack, flux_stack, ivar_stack, mask_stack,
                  outmask, norder=None, title='', qafile=None, show_telluric=False):
    """

    Routine to creqate QA for showing the individual spectrum
    compared to the combined stacked spectrum indicating which pixels
    were rejected.

    Args:
        wave (`numpy.ndarray`_):
            Wavelength array for spectrum of the exposure in
            question. Shape is (nspec,).
        flux (`numpy.ndarray`_):
            Flux for the exposure in question. Shape is (nspec,).
        ivar (`numpy.ndarray`_):
             Inverse variance for the exposure in question. Shape is
             (nspec,).
        mask (`numpy.ndarray`_): 
             Boolean array with mask for the exposure in question
             True=Good. If not specified determined form inverse
             variance.  Shape is (nspec,).
        flux_stack (`numpy.ndarray`_):
             Stacked spectrum to be compared to the exposure in
             question. Shape is (nspec,).
        ivar_stack (`numpy.ndarray`_): 
            Inverse variance of the stacked spectrum. Shape is
            (nspec,).
        mask_stack (`numpy.ndarray`_):
            Boolean array with mask for stacked spectrum. Shape is
            (nspec,).
        norder (:obj:`int`, optional):
            Indicate the number of orders if this is an echelle
            stack. If None, ...
        title (:obj:`str`, optional):
            Plot title
        qafile (:obj:`str`, optional):
            QA file name
        show_telluric (:obj:`bool`, optional):
            Show the atmospheric absorption if wavelengths > 9000A are covered by the spectrum

    """


    fig = plt.figure(figsize=(14, 8))
    spec_plot = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Get limits
    ymin, ymax = get_ylim(flux_stack, ivar_stack, mask_stack)

    # Plot spectrum
    rejmask = mask & np.logical_not(outmask)
    wave_mask = wave > 1.0
    wave_stack_mask = wave_stack > 1.0
    spec_plot.plot(wave[rejmask], flux[rejmask],'s',zorder=10,mfc='None', mec='r', label='rejected pixels')
    spec_plot.plot(wave[np.logical_not(mask)], flux[np.logical_not(mask)],'v', zorder=10, mfc='None', mec='orange',
                   label='originally masked')

    if norder is None:
        spec_plot.plot(wave[wave_mask], flux[wave_mask], color='dodgerblue', drawstyle='steps-mid',
                       zorder=2, alpha=0.5,label='single exposure')
        spec_plot.plot(wave[wave_mask], np.sqrt(utils.inverse(rejivar[wave_mask])),zorder=3,
                       color='0.7', alpha=0.5, drawstyle='steps-mid')
        spec_plot.plot(wave_stack[wave_stack_mask],flux_stack[wave_stack_mask]*mask_stack[wave_stack_mask],color='k',
                       drawstyle='steps-mid',lw=2,zorder=3, alpha=0.5, label='coadd')

        # TODO Use one of our telluric models here instead
        # Plot transmission
        if (np.max(wave[mask]) > 9000.0) and show_telluric:
            skytrans_file = data.get_skisim_filepath('atm_transmission_secz1.5_1.6mm.dat')
            skycat = np.genfromtxt(skytrans_file, dtype='float')
            scale = 0.8 * ymax
            spec_plot.plot(skycat[:, 0] * 1e4, skycat[:, 1] * scale, 'm-', alpha=0.5, zorder=11)
    else:
        npix = np.size(flux)
        nspec = int(npix / norder)
        spec_plot.plot(wave_stack[wave_stack_mask], flux_stack[wave_stack_mask] * mask_stack[wave_stack_mask],
                       color='k', drawstyle='steps-mid', lw=1, zorder=3, alpha=0.5, label='coadd')
        for iord in range(norder):
            spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           flux[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           drawstyle='steps-mid', zorder=1, alpha=0.7, label='order {:d}'.format(iord))

    # This logic below allows more of the spectrum to be plotted if wave_ref is a multi-order stack which has broader
    # wavelength coverage. For the longslit or single order case, this will plot the correct range as well
    wave_min = np.fmax(0.8*wave[wave_mask].min(), wave_stack[wave_stack_mask].min())
    wave_max = np.fmin(1.2*wave[wave_mask].max(), wave_stack[wave_stack_mask].max())

    # properties
    spec_plot.legend(fontsize=13)
    spec_plot.set_ylim([ymin, ymax])
    spec_plot.set_xlim((wave_min, wave_max))
    spec_plot.set_xlabel('Wavelength ($\\rm\\AA$)')
    spec_plot.set_ylabel('Flux')
    spec_plot.set_title(title, fontsize=16, color='red')
    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote QA: {:s}".format(qafile))
    plt.show()

def weights_qa(waves, weights, gpms, title='', colors=None):
    """
    Routine to make a QA plot for the weights used to compute a stacked spectrum.

    Parameters
    ----------
    wave : list
        List of `numpy.ndarray`_ float 1d wavelength arrays for spectra that
        went into a stack.
    weights : list
        List of `numpy.ndarray`_ float 1d (S/N)^2 weight arrays for the
        exposures that went into a stack. This would have been computed by
        sn_weights
    gpm : list
        List of `numpy.ndarray`_ boolean 1d good-pixel mask arrays for the
        exposures that went into a stack.  Good=True.
    title : str, optional
        Title for the plot.
    """

    if colors is None:
        colors = utils.distinct_colors(len(waves))

    fig = plt.figure(figsize=(12, 8))
    wave_min, wave_max = [], []
    for ii, (wave, weight, gpm) in enumerate(zip(waves, weights, gpms)):
        wave_gpm = wave > 1.0
        plt.plot(wave[wave_gpm], weight[wave_gpm]*gpm[wave_gpm], color=colors[ii])
        wave_min.append(wave[wave_gpm].min())
        wave_max.append(wave[wave_gpm].max())

    plt.xlim(np.min(wave_min), np.max(wave_max))
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Weights')
    plt.title(title, fontsize=16, color='red')
    plt.show()

def coadd_qa(wave, flux, ivar, nused, gpm=None, tell=None,
             title=None, qafile=None, show_telluric=False):
    """
    Routine to make QA plot of the final stacked spectrum. It works for both
    longslit/mulitslit, coadded individual order spectrum of the Echelle data
    and the final coadd of the Echelle data.

    Parameters
    ----------
    wave : `numpy.ndarray`_, shape=(nspec,)
        one-d wavelength array of your spectrum
    flux : `numpy.ndarray`_, shape=(nspec,)
        one-d flux array of your spectrum
    ivar : `numpy.ndarray`_, shape=(nspec,)
        one-d ivar array of your spectrum
    nused : `numpy.ndarray`_, shape=(nspec,)
        how many exposures used in the stack for each pixel, the same size with
        flux.
    gpm : `numpy.ndarray`_, optional, shape=(nspec,)
        boolean good pixel mask array for your spectrum. Good=True.
    tell : `numpy.ndarray`_, optional, shape=(nspec,)
        one-d telluric array for your spectrum
    title : str, optional
        plot title
    qafile : str, optional
        QA file name
    show_telluric : bool, optional
        Show a telluric absorptoin model on top of the data if wavelengths cover > 9000A
    """
    #TODO: This routine should take a parset

    if gpm is None:
        gpm = ivar > 0.0

    wave_gpm = wave > 1.0
    wave_min = wave[wave_gpm].min()
    wave_max = wave[wave_gpm].max()
    fig = plt.figure(figsize=(12, 8))
    # plot how may exposures you used at each pixel
    # [left, bottom, width, height]
    num_plot =  fig.add_axes([0.10, 0.70, 0.80, 0.23])
    spec_plot = fig.add_axes([0.10, 0.10, 0.80, 0.60], sharex=num_plot)
    num_plot.plot(wave[wave_gpm],nused[wave_gpm],drawstyle='steps-mid',color='k',lw=2)
    num_plot.set_xlim([wave_min, wave_max])
    num_plot.set_ylim([0.0, np.fmax(1.1*nused.max(), nused.max()+1.0)])
    num_plot.set_ylabel('$\\rm N_{EXP}$')
    num_plot.yaxis.set_major_locator(MaxNLocator(integer=True))
    num_plot.yaxis.set_minor_locator(NullLocator())

    # Plot spectrum
    spec_plot.plot(wave[wave_gpm], flux[wave_gpm], color='black', drawstyle='steps-mid',zorder=1,alpha=0.8, label='Single exposure')
    spec_plot.plot(wave[wave_gpm], np.sqrt(utils.inverse(ivar[wave_gpm])),zorder=2, color='red', alpha=0.7,
                   drawstyle='steps-mid', linestyle=':')

    # Get limits
    ymin, ymax = get_ylim(flux, ivar, gpm)

    # Plot transmission
    if (np.max(wave[gpm])>9000.0) and (tell is None) and show_telluric:
        skytrans_file = data.get_skisim_filepath('atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ymax
        spec_plot.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5,zorder=11)
    elif tell is not None:
        scale = 0.8*ymax
        spec_plot.plot(wave[wave_gpm], tell[wave_gpm]*scale, drawstyle='steps-mid', color='m',alpha=0.5,zorder=11)

    spec_plot.set_ylim([ymin, ymax])
    spec_plot.set_xlim([wave_min, wave_max])
    spec_plot.set_xlabel('Wavelength ($\\rm\\AA$)')
    spec_plot.set_ylabel('Flux')

    if title is not None:
        num_plot.set_title(title,color='red',fontsize=16)

    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote QA: {:s}".format(qafile))
    plt.show()

def update_errors(fluxes, ivars, masks, fluxes_stack, ivars_stack, masks_stack,
                  sn_clip=30.0, title='', debug=False):
    '''
    Determine corrections to errors using the residuals of each exposure about a preliminary stack. This routine is
    used as part of the iterative masking/stacking loop to determine the corrections to the errors used to reject pixels
    for the next iteration of the stack. The routine returns a set of corrections for each of the exposures that is input.

    Args:
        fluxes (`numpy.ndarray`_):
            fluxes for each exposure on the native wavelength grids
            shape=(nspec, nexp)
        ivars (`numpy.ndarray`_):
            Inverse variances for each exposure on the native wavelength grids
            shape=(nspec, nexp)
        masks (`numpy.ndarray`_):
            Boolean masks for each exposure on the native wavelength grids. True=Good.
            shape=(nspec, nexp)
        fluxes_stack (`numpy.ndarray`_):
            Stacked spectrum for this iteration interpolated on the native wavelength grid of the fluxes exposures.
            shape=(nspec, nexp)
        ivars_stack (`numpy.ndarray`_):
            Inverse variances of stacked spectrum for this iteration interpolated on the native wavelength grid of the
            fluxes exposures.
            shape=(nspec, nexp)
        masks_stack (`numpy.ndarray`_):
            Boolean mask of stacked spectrum for this iteration interpolated on the native wavelength grid of the fluxes exposures.
            shape=(nspec, nexp)
        sn_clip (float, optional):
            Errors are capped in output rejivars so that the S/N is never greater than sn_clip. This prevents overly
            aggressive rejection in high S/N ratio spectra which neverthless differ at a level greater than the implied S/N due to
            systematics.
        title (str, optional):
            Title for QA plot
        debug (bool, optional):
            If True, show QA plots useful for debuggin.

    Returns:
        tuple:
            - rejivars: `numpy.ndarray`_, (nspec, nexp): Updated inverse
              variances to be used in rejection
            - sigma_corrs, `numpy.ndarray`_, (nexp): Array of correction factors
              applied to the original ivars to get the new rejivars
            - outchi: `numpy.ndarray`_, (nspec, nexp): The original
              chi=(fluxes-fluxes_stack)*np.sqrt(ivars) used to determine
              the correction factors. This quantity is useful for
              plotting. Note that the outchi is computed using the
              original non-corrected errors.
            - maskchi: `numpy.ndarray`_, bool, (nspec, nexp): Mask returned by
              renormalize_erorrs indicating which pixels were used in
              the computation of the correction factors. This is
              basically the union of the input masks but with chi > clip
              (clip=6.0 is the default) values clipped out.
    '''

    if fluxes.ndim == 1:
        nexp = 1
    else:
        nexp = np.shape(fluxes)[1]

    outchi = np.zeros_like(ivars)
    maskchi = np.zeros_like(outchi,dtype=bool)
    rejivars = np.zeros_like(outchi)
    sigma_corrs = np.zeros(nexp)
    outmasks = np.copy(masks)

    # Loop on images to update noise model for rejection
    for iexp in range(nexp):
        if fluxes.ndim>1:
            # Grab the spectrum
            thisflux = fluxes[:, iexp]
            thisivar = ivars[:, iexp]
            thismask = outmasks[:,iexp]
            # Grab the stack interpolated with the same grid as the current exposure
            thisflux_stack = fluxes_stack[:, iexp]
            thisvar_stack = utils.inverse(ivars_stack[:, iexp])
            thismask_stack = masks_stack[:, iexp]
        else:
            thisflux = fluxes
            thisivar = ivars
            thismask = outmasks
            # Grab the stack interpolated with the same grid as the current exposure
            thisflux_stack = fluxes_stack
            thisvar_stack = utils.inverse(ivars_stack)
            thismask_stack = masks_stack

        # var_tot = total variance of the quantity (fluxes - fluxes_stack), i.e. the quadrature sum of the two variances
        var_tot = thisvar_stack + utils.inverse(thisivar)
        mask_tot = thismask & thismask_stack
        ivar_tot = utils.inverse(var_tot)

        # Impose the S/N clipping threshold before computing chi and renormalizing the errors
        ivar_clip = mask_tot*utils.clip_ivar(thisflux_stack, ivar_tot, sn_clip, gpm=mask_tot)
        # TODO Do we need the offset code to re-center the chi? If so add it right here into the chi
        chi = np.sqrt(ivar_clip)*(thisflux - thisflux_stack)
        # Adjust errors to reflect the statistics of the distribution of errors. This fixes cases where the
        # the noise model is not quite right
        this_sigma_corr, igood = renormalize_errors(chi, mask_tot, clip=6.0, max_corr=5.0, title=title, debug=debug)
        ivar_tot_corr = ivar_clip/this_sigma_corr ** 2
        # TODO is this correct below? JFH Thinks no
        #ivar_cap = utils.clip_ivar(thisflux_stack, ivar_tot_corr, sn_clip, mask=mask_tot)
        #ivar_cap = np.minimum(ivar_tot_corr, (sn_clip/(thisflux_stack + (thisflux_stack <= 0.0))) ** 2)
        if fluxes.ndim>1:
            sigma_corrs[iexp] = this_sigma_corr
            rejivars[:, iexp] = ivar_tot_corr
            outchi[:, iexp] = chi
            maskchi[:, iexp] = igood
        else:
            sigma_corrs = np.array([this_sigma_corr])
            rejivars = ivar_tot_corr
            outchi = chi
            maskchi = igood


    return rejivars, sigma_corrs, outchi, maskchi


def spec_reject_comb(wave_grid, wave_grid_mid, waves_list, fluxes_list, ivars_list, gpms_list, weights_list, sn_clip=30.0, lower=3.0, upper=3.0,
                     maxrej=None, maxiter_reject=5, title='', debug=False,
                     verbose=False):
    """
    Routine for executing the iterative combine and rejection of a set of
    spectra to compute a final stacked spectrum.

    Parameters
    ----------
    wave_grid : `numpy.ndarray`_
        new wavelength grid desired. This will typically be a reguarly spaced
        grid created by the get_wave_grid routine.  The reason for the ngrid+1
        is that this is the general way to specify a set of  bins if you desire
        ngrid bin centers, i.e. the output stacked spectra have ngrid elements.
        The spacing of this grid can be regular in lambda (better for multislit)
        or log lambda (better for echelle). This new wavelength grid should be
        designed with the sampling of the data in mind. For example, the code
        will work fine if you choose the sampling to be too fine, but then the
        number of exposures contributing to any given wavelength bin will be one
        or zero in the limiting case of very small wavelength bins. For larger
        wavelength bins, the number of exposures contributing to a given bin
        will be larger.  shape=(ngrid +1,)
    wave_grid_mid : `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See
        core.wavecal.wvutils.py for more.  shape=(ngrid,)
    waves : list
        List of length nexp float `numpy.ndarray`_ wavelength arrays for spectra
        to be stacked. Note that the wavelength grids can in general be
        different for each exposure,  irregularly spaced, and have different
        sizes.
    fluxes : list
        List of length nexp float `numpy.ndarray`_ fluxes for each exposure
        aligned with the wavelength arrays in waves
    ivars : list
        List of length nexp float `numpy.ndarray`_ inverse variances for each
        exposure aligned with the wavelength arrays in waves
    gpms : list
        List of length nexp boolean `numpy.ndarray`_ good pixel masks for each
        exposure aligned with the wavelength arrays in waves.  True=Good.
    weights : list
        List of length nexp float `numpy.ndarray`_ weights for each exposure
        aligned with the wavelength arrays in waves.  These are computed using
        sn_weights
    sn_clip : float, optional, default=30.0
        Errors are capped during rejection so that the S/N is never greater than
        sn_clip. This prevents overly aggressive rejection in high S/N ratio
        spectrum which neverthless differ at a level greater than the implied
        S/N due to systematics.
    lower : float, optional, default=3.0
        lower rejection threshold for djs_reject
    upper : float, optional, default=3.0
        upper rejection threshold for djs_reject
    maxrej : int, optional, default=None
        maximum number of pixels to reject in each iteration for djs_reject.
    maxiter_reject : int, optional, default=5
        maximum number of iterations for stacking and rejection. The code stops
        iterating either when the output mask does not change betweeen
        successive iterations or when maxiter_reject is reached.
    title : str, optional
        Title for QA plot
    debug : bool, optional, default=False
        Show QA plots useful for debugging.
    verbose : bool, optional, default=False
        Level of verbosity.

    Returns
    -------
    wave_stack : `numpy.ndarray`_
        Wavelength grid for stacked
        spectrum. As discussed above, this is the weighted average
        of the wavelengths of each spectrum that contriuted to a
        bin in the input wave_grid wavelength grid. It thus has
        ngrid elements, whereas wave_grid has ngrid+1 elements to
        specify the ngrid total number of bins. Note that
        wave_stack is NOT simply the wave_grid bin centers, since
        it computes the weighted average.
        shape=(ngrid,)
    flux_stack : `numpy.ndarray`_
        Final stacked spectrum on
        wave_stack wavelength grid
        shape=(ngrid,)
    ivar_stack : `numpy.ndarray`_
        Inverse variance spectrum
        on wave_stack wavelength grid. Erors are propagated
        according to weighting and masking.
        shape=(ngrid,)
    gpm_stack : `numpy.ndarray`_
        Boolean mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
        shape=(ngrid,)
    nused : `numpy.ndarray`_
        Number of exposures which
        contributed to each pixel in the wave_stack. Note that
        this is in general different from nexp because of masking,
        but also becuse of the sampling specified by wave_grid. In
        other words, sometimes more spectral pixels in the
        irregularly gridded input wavelength array waves will land
        in one bin versus another depending on the sampling.
        shape=(ngrid,)
    out_gpms : list
        List of length nexp output `numpy.ndarray`_ good pixel masks for each
        exposure aligned with the wavelength arrays in waves indicating which
        pixels are rejected in each exposure of the original input spectra after
        performing all of the iterations of combine/rejection
    """
    # Conver the input lists to arrays. This is currently needed since the rejection routine pydl.djs_reject requires
    # arrays as input.
    waves, nspec_list = utils.explist_to_array(waves_list, pad_value=0.0)
    fluxes, _ = utils.explist_to_array(fluxes_list, pad_value=0.0)
    ivars, _ = utils.explist_to_array(ivars_list, pad_value=0.0)
    weights, _ = utils.explist_to_array(weights_list, pad_value=0.0)
    gpms, _ = utils.explist_to_array(gpms_list, pad_value=False)
    this_gpms = np.copy(gpms)
    iter = 0
    qdone = False
    while (not qdone) and (iter < maxiter_reject):
        # Compute the stack
        wave_stack, flux_stack, ivar_stack, gpm_stack, nused = compute_stack(
            wave_grid, waves_list, fluxes_list, ivars_list, utils.array_to_explist(this_gpms, nspec_list=nspec_list), weights_list)
        # Interpolate the individual spectra onto the wavelength grid of the stack. Use wave_grid_mid for this
        # since it has no masked values
        flux_stack_nat, ivar_stack_nat, gpm_stack_nat, _ = interp_spec(
            waves, wave_grid_mid, flux_stack, ivar_stack, gpm_stack)
        ## TESTING
        #nused_stack_nat, _, _ = interp_spec(
        #    waves, wave_grid_mid, nused, ivar_stack, mask_stack)
        #embed()
        rejivars, sigma_corrs, outchi, chigpm = update_errors(fluxes, ivars, this_gpms,
                                                               flux_stack_nat,  ivar_stack_nat, gpm_stack_nat, sn_clip=sn_clip)
        this_gpms, qdone = pydl.djs_reject(fluxes, flux_stack_nat, outmask=this_gpms,inmask=gpms, invvar=rejivars,
                                          lower=lower,upper=upper, maxrej=maxrej, sticky=False)
        iter += 1


    if (iter == maxiter_reject) & (maxiter_reject != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter_reject) + ' reached in spec_reject_comb')
    out_gpms = np.copy(this_gpms)
    out_gpms_list = utils.array_to_explist(out_gpms, nspec_list=nspec_list)

    # print out a summary of how many pixels were rejected
    nexp = waves.shape[1]
    nrej = np.sum(np.logical_not(out_gpms) & gpms, axis=0)
    norig = np.sum((waves > 1.0) & np.logical_not(gpms), axis=0)

    if verbose:
        for iexp in range(nexp):
            # nrej = pixels that are now masked that were previously good
            msgs.info("Rejected {:d} pixels in exposure {:d}/{:d}".format(nrej[iexp], iexp, nexp))

    # Compute the final stack using this outmask
    wave_stack, flux_stack, ivar_stack, gpm_stack, nused = compute_stack(
        wave_grid, waves_list, fluxes_list, ivars_list, out_gpms_list, weights_list)

    # Used only for plotting below
    if debug:
        # TODO Add a line here to optionally show the distribution of all pixels about the stack as we do for X-shooter.
        for iexp in range(nexp):
            # plot the residual distribution for each exposure
            title_renorm = title + ': Error distriution about stack for exposure {:d}/{:d}'.format(iexp,nexp)
            renormalize_errors_qa(outchi[:, iexp], chigpm[:, iexp], sigma_corrs[iexp], title=title_renorm)
            # plot the rejections for each exposures
            title_coadd_iexp = title + ': nrej={:d} pixels rejected,'.format(nrej[iexp]) + \
                               ' norig={:d} originally masked,'.format(norig[iexp]) + \
                               ' for exposure {:d}/{:d}'.format(iexp,nexp)
            # JFH: QA should use wave_grid_mid
            coadd_iexp_qa(waves[:, iexp], fluxes[:, iexp], rejivars[:, iexp], gpms[:, iexp], wave_grid_mid, flux_stack,
                          ivar_stack, gpm_stack, out_gpms[:, iexp], qafile=None, title=title_coadd_iexp)
        # weights qa
        title_weights = title + ': Weights Used -- nrej={:d} total pixels rejected,'.format(np.sum(nrej)) + \
                        ' norig={:d} originally masked'.format(np.sum(norig))
        weights_qa(waves_list, weights_list, out_gpms_list, title=title_weights)



    return wave_stack, flux_stack, ivar_stack, gpm_stack, nused, out_gpms_list

def scale_spec_stack(wave_grid, wave_grid_mid, waves, fluxes, ivars, gpms, sns, weights,
                     ref_percentile=70.0, maxiter_scale=5,
                     sigrej_scale=3.0, scale_method='auto',
                     hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                     debug=False, show=False):

    """
    Scales a set of spectra to a common flux scale. This is done by first
    computing a stack of the spectra and then scaling each spectrum to match the
    composite with the :func:`scale_spec` algorithm which performs increasingly
    sophisticated scaling depending on the S/N ratio.

    Parameters
    ----------
    wave_grid : `numpy.ndarray`_
        New wavelength grid desired. This will typically be a reguarly spaced
        grid created by the get_wave_grid routine.  The reason for the ngrid+1
        is that this is the general way to specify a set of  bins if you desire
        ngrid bin centers, i.e. the output stacked spectra have ngrid elements.
        The spacing of this grid can be regular in lambda (better for multislit)
        or log lambda (better for echelle). This new wavelength grid should be
        designed with the sampling of the data in mind. For example, the code
        will work fine if you choose the sampling to be too fine, but then the
        number of exposures contributing to any given wavelength bin will be one
        or zero in the limiting case of very small wavelength bins. For larger
        wavelength bins, the number of exposures contributing to a given bin
        will be larger.  shape=(ngrid +1,)
    wave_grid_mid : `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See
        core.wavecal.wvutils.py for more.  shape=(ngrid,)
    waves : list
        List of length nexp float `numpy.ndarray`_ of wavelengths of the spectra
        to be stacked. Note that wavelength arrays can in general be different,
        irregularly spaced, and have different sizes.
    fluxes : list
        List of length nexp float `numpy.ndarray`_ of fluxes for each exposure
        aligned with the wavelength grids in waves.
    ivars : list
        List of length nexp float `numpy.ndarray`_ of inverse variances for each
        exposure aligned with the wavelength grids in waves.
    gpms : list
        List of length nexp bool `numpy.ndarray`_ of good pixel masks for each
        exposure aligned with the wavelength grids in waves. True=Good.
    sns : `numpy.ndarray`_
        sn of each spectrum in the list of exposures used to determine which
        scaling method should be used. This can be computed using sn_weights.
        shape=(nexp,)
    sigrej_scale : float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra
        with scale_spec.
    ref_percentile : float, optional, default=70.0
        percentile fraction cut used for selecting minimum SNR cut for
        robust_median_ratio
    maxiter_scale : int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    scale_method : str, optional, default='auto'
        Options are auto, poly, median, none, or hand. Hand is not well tested.
        User can optionally specify the rescaling method. Default is to let the
        code determine this automitically which works well.
    hand_scale : `numpy.ndarray`_, optional
        array of hand scale factors, not well tested
    sn_min_polyscale : float, optional, default=2.0
        maximum SNR for perforing median scaling
    sn_min_medscale : float, optional, default=0.5
        minimum SNR for perforing median scaling
    debug : bool, optional, default=False
        show interactive QA plot

    Returns
    -------
    fluxes_scales : list
        List of length nexp `numpy.ndarray`_ rescaled fluxes
    ivars_scales : list
        List of length nexp `numpy.ndarray`_ rescaled ivars.  shape=(nspec,
        nexp)
    scales : list
        List of length nexp `numpy.ndarray`_ scale factors applied to each
        individual spectra and their inverse variances.
    scale_method_used : list
        List of methods used for rescaling spectra.
    """

    # Compute an initial stack as the reference, this has its own wave grid based on the weighted averages
    wave_stack, flux_stack, ivar_stack, gpm_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, gpms, weights)

    # Rescale spectra to line up with our preliminary stack so that we can sensibly reject outliers
    fluxes_scale, ivars_scale, scales, scale_method_used = [], [], [], []
    for iexp, (wave, flux, ivar, gpm, sn) in enumerate(zip(waves, fluxes, ivars, gpms, sns)):
        hand_scale_iexp = None if hand_scale is None else hand_scale[iexp]
        # JFH Changed to used wave_grid_mid in interpolation to avoid interpolating onto zeros
        fluxes_scale_iexp, ivars_scale_iexp, scales_iexp, scale_method_iexp = scale_spec(
            wave, flux, ivar, sn, wave_grid_mid, flux_stack, ivar_stack,
            mask=gpm, mask_ref=gpm_stack, ref_percentile=ref_percentile, maxiters=maxiter_scale,
            sigrej=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale_iexp, sn_min_polyscale=sn_min_polyscale,
            sn_min_medscale=sn_min_medscale, debug=debug, show=show)
        fluxes_scale.append(fluxes_scale_iexp)
        ivars_scale.append(ivars_scale_iexp)
        scales.append(scales_iexp)
        scale_method_used.append(scale_method_iexp)

    return fluxes_scale, ivars_scale, scales, scale_method_used


def combspec(waves, fluxes, ivars, gpms, sn_smooth_npix,
             wave_method='linear', dwave=None, dv=None, dloglam=None,
             spec_samp_fact=1.0, wave_grid_min=None, wave_grid_max=None,
             ref_percentile=70.0, maxiter_scale=5, wave_grid_input=None,
             sigrej_scale=3.0, scale_method='auto', hand_scale=None,
             sn_min_polyscale=2.0, sn_min_medscale=0.5,
             weight_method='auto', maxiter_reject=5, sn_clip=30.0,
             lower=3.0, upper=3.0, maxrej=None, qafile=None, title='', debug=False,
             debug_scale=False, show_scale=False, show=False, verbose=True):

    '''
    Main driver routine for coadding longslit/multi-slit spectra.
    NEED LOTS MORE HERE

    Parameters
    ----------
    waves : list
        List of length nexp containing `numpy.ndarray`_ float wavelength arrays
        for spectra to be stacked. These wavelength  arrays can in general be
        different, irregularly spaced, and have different sizes.
    fluxes : list
        List of length nexp containing `numpy.ndarray`_ float flux arrays for
        spectra to be stacked. These should be aligned with the wavelength
        arrays in waves.
    ivars : list
        List of length nexp containing `numpy.ndarray`_ float ivar arrays for
        spectra to be stacked. These should be aligned with the wavelength
        arrays in waves.
    gpms : list
        List of length nexp containing `numpy.ndarray`_ boolean good pixel mask
        arrays for spectra to be stacked. These should be aligned with the
        wavelength arrays in waves.
    sn_smooth_npix : int
        Number of pixels to median filter by when computing S/N used to decide
        how to scale and weight spectra
    wave_method : str, optional
        method for generating new wavelength grid with get_wave_grid. Deafult is
        'linear' which creates a uniformly space grid in lambda. See
        docuementation on get_wave_grid for description of the options.
    dwave : float, optional
        dispersion in units of A in case you want to specify it for
        get_wave_grid, otherwise the code computes the median spacing from the
        data.
    dv : float, optional
        Dispersion in units of km/s in case you want to specify it in the
        get_wave_grid  (for the 'velocity' option), otherwise a median value is
        computed from the data.
    spec_samp_fact : float, optional
        Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or
        coarser (spec_samp_fact > 1.0) by this sampling factor. This basically
        multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
        spec_samp_fact are pixels.
    wave_grid_min : float, optional
        In case you want to specify the minimum wavelength in your wavelength
        grid, default=None computes from data.
    wave_grid_max : float, optional
        In case you want to specify the maximum wavelength in your wavelength
        grid, default=None computes from data.
    ref_percentile : float, optional
        percentile fraction cut used for selecting minimum SNR cut for
        robust_median_ratio
    maxiter_scale : int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    wave_grid_input : `numpy.ndarray`_, optional
        User input wavelength grid to be used with the 'user_input' wave_method.
        Shape is (nspec_input,)
    maxiter_reject : int, optional, default=5
        maximum number of iterations for stacking and rejection. The code stops
        iterating either when the output mask does not change betweeen
        successive iterations or when maxiter_reject is reached.
    sigrej_scale : float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra
        with scale_spec.
    scale_method : str, optional
        Options are auto, poly, median, none, or hand. Hand is not well tested.
        User can optionally specify the rescaling method.
        Default is 'auto' will let the code determine this automitically which works well.
    hand_scale : `numpy.ndarray`_, optional
        Array of hand scale factors, not well tested
    sn_min_polyscale : float, optional, default = 2.0
        maximum SNR for perforing median scaling
    sn_min_medscale : float, optional, default = 0.5
        minimum SNR for perforing median scaling
    weight_method (str):
        Weight method to use for coadding spectra (see
            :func:`~pypeit.core.coadd.sn_weights`) for documentation. Default='auto'
    sn_clip: float, optional, default=30.0
        Errors are capped during rejection so that the S/N is never greater than
        sn_clip. This prevents overly aggressive rejection in high S/N ratio
        spectrum which neverthless differ at a level greater than the implied
        S/N due to systematics.
    lower : float, optional, default=3.0
        lower rejection threshold for djs_reject
    upper : float, optional, default=3.0
        upper rejection threshold for djs_reject
    maxrej : int, optional
        maximum number of pixels to reject in each iteration for djs_reject.
    qafile : str, default=None
        Root name for QA, if None, it will be determined from the outfile
    title : str, optional
        Title for QA plots
    debug: bool, default=False
        Show all QA plots useful for debugging. Note there are lots of QA plots,
        so only set this to True if you want to inspect them all.
    debug_scale : bool, default=False
        show interactive QA plots for the rescaling of the spectra
    show : bool, default=False
        If True, show key QA plots or not
    show_scale: bool, default=False
        If True, show interactive QA plots for the rescaling of the spectra

    Returns
    -------
    wave_grid_mid : `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See
        core.wavecal.wvutils.py for more.  shape=(ngrid,)
    wave_stack : `numpy.ndarray`_
        Wavelength grid for stacked
        spectrum. As discussed above, this is the weighted average
        of the wavelengths of each spectrum that contriuted to a
        bin in the input wave_grid wavelength grid. It thus has
        ngrid elements, whereas wave_grid has ngrid+1 elements to
        specify the ngrid total number of bins. Note that
        wave_stack is NOT simply the wave_grid bin centers, since
        it computes the weighted average.
        shape=(ngrid,)
    flux_stack : `numpy.ndarray`_
        Final stacked spectrum on
        wave_stack wavelength grid
        shape=(ngrid,)
    ivar_stack : `numpy.ndarray`_
        Inverse variance spectrum on wave_stack
        wavelength grid. Erors are propagated according to
        weighting and masking.
        shape=(ngrid,)
    mask_stack : `numpy.ndarray`_
        Boolean mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
        shape=(ngrid,)
    '''

    #debug_scale=True
    #show_scale=True
    #debug=True
    #show=True

    #from IPython import embed
    #embed()
    # We cast to float64 because of a bug in np.histogram
    _waves = [np.float64(wave) for wave in waves]
    _fluxes = [np.float64(flux) for flux in fluxes]
    _ivars = [np.float64(ivar) for ivar in ivars]

    # Generate a giant wave_grid
    wave_grid, wave_grid_mid, _ = wvutils.get_wave_grid(
        waves=_waves, gpms=gpms, wave_method=wave_method,
        wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
        wave_grid_input=wave_grid_input,
        dwave=dwave, dv=dv, dloglam=dloglam, spec_samp_fact=spec_samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(_fluxes, _ivars, gpms, sn_smooth_npix=sn_smooth_npix, weight_method=weight_method, verbose=verbose)
    fluxes_scale, ivars_scale, scales, scale_method_used = scale_spec_stack(
        wave_grid, wave_grid_mid, _waves, _fluxes, _ivars, gpms, rms_sn, weights, ref_percentile=ref_percentile, maxiter_scale=maxiter_scale,
        sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
        sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale, debug=debug_scale, show=show_scale)
    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, gpm_stack, nused, outmask = spec_reject_comb(
        wave_grid, wave_grid_mid, _waves, fluxes_scale, ivars_scale, gpms, weights, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug, title=title)

    if show:
        # JFH Use wave_grid_mid for QA plots
        coadd_qa(wave_grid_mid, flux_stack, ivar_stack, nused, gpm=gpm_stack, title='Stacked spectrum', qafile=qafile)

    return wave_grid_mid, wave_stack, flux_stack, ivar_stack, gpm_stack

def multi_combspec(waves, fluxes, ivars, masks, sn_smooth_npix=None,
                   wave_method='linear', dwave=None, dv=None, dloglam=None, spec_samp_fact=1.0, wave_grid_min=None,
                   wave_grid_max=None, ref_percentile=70.0, maxiter_scale=5,
                   sigrej_scale=3.0, scale_method='auto', hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                   weight_method='auto', maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0,
                   maxrej=None,
                   qafile=None, debug=False, debug_scale=False, show_scale=False, show=False):
    """
    Routine for coadding longslit/multi-slit spectra. Calls combspec which is
    the main stacking algorithm.

    Parameters
    ----------
    waves : list
        List of `numpy.ndarray`_ wavelength arrays with shape (nspec_i,) with
        the wavelength arrays of the spectra to be coadded.
    fluxes : list
        List of `numpy.ndarray`_ wavelength arrays with shape (nspec_i,) with
        the flux arrays of the spectra to be coadded.
    ivars : list
        List of `numpy.ndarray`_ wavelength arrays with shape (nspec_i,) with
        the ivar arrays of the spectra to be coadded.
    masks : list
        List of `numpy.ndarray`_ wavelength arrays with shape (nspec_i,) with
        the mask arrays of the spectra to be coadded.
    sn_smooth_npix : int, optional
        Number of pixels to median filter by when computing S/N used to decide
        how to scale and weight spectra. If set to None, the code will determine
        the effective number of good pixels per spectrum in the stack that is
        being co-added and use 10% of this neff.
    wave_method : str, optional
        method for generating new wavelength grid with get_wave_grid. Deafult is
        'linear' which creates a uniformly space grid in lambda. See
        docuementation on get_wave_grid for description of the options.
    dwave : float, optional
        dispersion in units of A in case you want to specify it for
        get_wave_grid, otherwise the code computes the median spacing from the
        data.
    dv : float, optional
        Dispersion in units of km/s in case you want to specify it in the
        get_wave_grid  (for the 'velocity' option), otherwise a median value is
        computed from the data.
    spec_samp_fact : float, optional
        Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or
        coarser (spec_samp_fact > 1.0) by this sampling factor. This basically
        multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
        spec_samp_fact are pixels.
    wave_grid_min : float, optional
        In case you want to specify the minimum wavelength in your wavelength
        grid, default=None computes from data.
    wave_grid_max : float, optional
        In case you want to specify the maximum wavelength in your wavelength
        grid, default=None computes from data.
    wave_grid_input : `numpy.ndarray`_
        User input wavelength grid to be used with the 'user_input' wave_method.
        Shape is (nspec_input,)
    maxiter_reject : int, optional
        maximum number of iterations for stacking and rejection. The code stops
        iterating either when the output mask does not change betweeen
        successive iterations or when maxiter_reject is reached. Default=5.
    ref_percentile : float, optional
        percentile fraction cut used for selecting minimum SNR cut for
        robust_median_ratio. Should be a number between 0 and 100, default =
        70.0
    maxiter_scale : int, optional
        Maximum number of iterations performed for rescaling spectra. Default=5.
    sigrej_scale : float, optional
        Rejection threshold used for rejecting pixels when rescaling spectra
        with scale_spec. Default=3.0
    scale_method : str, optional
        Options are auto, poly, median, none, or hand. Hand is not well tested.
        User can optionally specify the rescaling method. Default='auto' will
        let the code determine this automitically which works well.
    hand_scale : `numpy.ndarray`_, optional
        Array of hand scale factors, not well tested
    sn_min_polyscale : float, optional
        maximum SNR for perforing median scaling
    sn_min_medscale : float, optional
        minimum SNR for perforing median scaling
    weight_method (str):
        Weight method to use for coadding spectra (see
            :func:`~pypeit.core.coadd.sn_weights`) for documentation. Default='auto'
    maxiter_reject : int, optional
        maximum number of iterations for stacking and rejection. The code stops
        iterating either when the output mask does not change betweeen
        successive iterations or when maxiter_reject is reached.
    sn_clip : float, optional
        Errors are capped during rejection so that the S/N is never greater than
        sn_clip. This prevents overly aggressive rejection in high S/N ratio
        spectrum which neverthless differ at a level greater than the implied
        S/N due to systematics.
    lower : float, optional
        lower rejection threshold for djs_reject
    upper : float, optional
        upper rejection threshold for djs_reject
    maxrej : int, optional
        maximum number of pixels to reject in each iteration for djs_reject.
    nmaskedge : int, optional
        Number of edge pixels to mask. This should be removed/fixed.
    qafile : str, optional, default=None
        Root name for QA, if None, it will be determined from the outfile
    outfile : str, optional, default=None,
        Root name for QA, if None, it will come from the target name from the
        fits header.
    debug : bool, optional, default=False,
        Show all QA plots useful for debugging. Note there are lots of QA plots,
        so only set this to True if you want to inspect them all.
    debug_scale : bool, optional, default=False
        show interactive QA plots for the rescaling of the spectra
    show : bool, optional, default=False
        Show key QA plots or not

    Returns
    -------
    wave_grid_mid : `numpy.ndarray`_, (ngrid,)
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity.  See
        core.wavecal.wvutils.py for more.
    wave_stack : `numpy.ndarray`_, (ngrid,)
        Wavelength grid for stacked spectrum. As discussed above, this is the
        weighted average of the wavelengths of each spectrum that contriuted to
        a bin in the input wave_grid wavelength grid. It thus has ngrid
        elements, whereas wave_grid has ngrid+1 elements to specify the ngrid
        total number of bins. Note that wave_stack is NOT simply the wave_grid
        bin centers, since it computes the weighted average.
    flux_stack : `numpy.ndarray`_, (ngrid,)
        Final stacked spectrum on wave_stack wavelength grid
    ivar_stack : `numpy.ndarray`_, (ngrid,)
        Inverse variance spectrum on wave_stack wavelength grid. Erors are
        propagated according to weighting and masking.
    mask_stack : `numpy.ndarray`_, bool, (ngrid,)
        Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
    """
    # Decide how much to smooth the spectra by if this number was not passed in
    if sn_smooth_npix is None:
        nexp = len(waves)
        # This is the effective good number of spectral pixels in the stack
        nspec_eff = np.sum([np.sum(wave > 1.0) for wave in waves]) / nexp
        sn_smooth_npix = int(np.round(0.1*nspec_eff))
        msgs.info('Using a sn_smooth_npix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    wave_grid_mid, wave_stack, flux_stack, ivar_stack, mask_stack = combspec(
        waves, fluxes,ivars, masks, wave_method=wave_method, dwave=dwave, dv=dv, dloglam=dloglam,
        spec_samp_fact=spec_samp_fact, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max, ref_percentile=ref_percentile,
        maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
        sn_min_medscale=sn_min_medscale, sn_min_polyscale=sn_min_polyscale, sn_smooth_npix=sn_smooth_npix,
        weight_method=weight_method, maxiter_reject=maxiter_reject, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej,  qafile=qafile, title='multi_combspec', debug=debug, debug_scale=debug_scale, show_scale=show_scale,
        show=show)

    return wave_grid_mid, wave_stack, flux_stack, ivar_stack, mask_stack


# TODO: Describe the "b_tuple"
def ech_combspec(waves_arr_setup, fluxes_arr_setup, ivars_arr_setup, gpms_arr_setup, weights_sens_arr_setup,
                 setup_ids = None, nbests=None,
                 wave_method='log10', dwave=None, dv=None, dloglam=None,
                 spec_samp_fact=1.0, wave_grid_min=None, wave_grid_max=None,
                 ref_percentile=70.0, maxiter_scale=5, niter_order_scale=3,
                 sigrej_scale=3.0, scale_method='auto',
                 hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                 maxiter_reject=5,
                 sn_clip=30.0, lower=3.0, upper=3.0,
                 maxrej=None, qafile=None, debug_scale=False, debug_order_stack=False,
                 debug_global_stack=False, debug=False,
                 show_order_stacks=False, show_order_scale=False,
                 show_exp=False, show=False, verbose=False):
    """
    Driver routine for coadding Echelle spectra. Calls combspec which is the
    main stacking algorithm. It will deliver three fits files:
    spec1d_order_XX.fits (stacked individual orders, one order per extension),
    spec1d_merge_XX.fits (straight combine of stacked individual orders),
    spec1d_stack_XX.fits (a giant stack of all exposures and all orders).  In
    most cases, you should use spec1d_stack_XX.fits for your scientific analyses
    since it reject most outliers.

    Parameters
    ----------
    waves_arr_setup : list of `numpy.ndarray`_
        List of wavelength arrays for spectra to be stacked. The length of the
        list is nsetups.  Each element of the list corresponds to a distinct
        setup and each numpy array has shape=(nspec, norder, nexp)
    fluxes_arr_setup : list of `numpy.ndarray`_
        List of flux arrays for spectra to be stacked. The length of the list is
        nsetups.  Each element of the list corresponds to a dinstinct setup and
        each numpy array has shape=(nspec, norder, nexp)
    ivars_arr_setup : list of `numpy.ndarray`_
        List of ivar arrays for spectra to be stacked. The length of the list is
        nsetups.  Each element of the list corresponds to a dinstinct setup and
        each numpy array has shape=(nspec, norder, nexp)
    gpms_arr_setup : list of `numpy.ndarray`_
        List of good pixel mask arrays for spectra to be stacked. The length of
        the list is nsetups.  Each element of the list corresponds to a
        dinstinct setup and each numpy array has shape=(nspec, norder, nexp)
    weights_sens_arr_setup : list of `numpy.ndarray`_
        List of sensitivity function weights required for relative weighting of
        the orders.  The length of the list is nsetups.  Each element of the
        list corresponds to a dinstinct setup and each numpy array has
        shape=(nspec, norder, nexp)
    setup_ids : list, optional, default=None
        List of strings indicating the name of each setup. If None uppercase
        letters A, B, C, etc. will be used
    nbests : int or list or  `numpy.ndarray`_, optional
        Integer or list of integers indicating the number of orders to use for
        estimating the per exposure weights per echelle setup.  Default is
        nbests=None, which will just use one fourth of the orders for a given
        setup.
    wave_method : str, optional
        method for generating new wavelength grid with get_wave_grid. Deafult is
        'log10' which creates a uniformly space grid in log10(lambda), which is
        typically the best for echelle spectrographs
    dwave : float, optional
        dispersion in units of A in case you want to specify it for
        get_wave_grid, otherwise the code computes the median spacing from the
        data.
    dv : float, optional
        Dispersion in units of km/s in case you want to specify it in the
        get_wave_grid  (for the 'velocity' option), otherwise a median value is
        computed from the data.
    dloglam : float, optional
        Dispersion in dimensionless units
    spec_samp_fact : float, optional
        Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or
        coarser (spec_samp_fact > 1.0) by this sampling factor. This basically
        multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
        spec_samp_fact are pixels.
    wave_grid_min : float, optional
        In case you want to specify the minimum wavelength in your wavelength
        grid, default=None computes from data.
    wave_grid_max : float, optional
        In case you want to specify the maximum wavelength in your wavelength
        grid, default=None computes from data.
    wave_grid_input : `numpy.ndarray`_, optional
        User input wavelength grid to be used with the 'user_input' wave_method.
        Shape is (nspec_input,)
    ref_percentile : float, optional
        percentile fraction cut used for selecting minimum SNR cut for
        robust_median_ratio
    maxiter_scale : int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    sigrej_scale : float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra
        with scale_spec.
    hand_scale : `numpy.ndarray`_, optional
        Array of hand scale factors, not well tested
    sn_min_polyscale : float, optional, default = 2.0
        maximum SNR for perforing median scaling
    sn_min_medscale : float, optional, default = 0.5
        minimum SNR for perforing median scaling
    maxiter_reject : int, optional, default=5
        maximum number of iterations for stacking and rejection. The code stops
        iterating either when the output mask does not change betweeen
        successive iterations or when maxiter_reject is reached.
    sn_clip : float, optional, default=30.0
        Errors are capped during rejection so that the S/N is never greater than
        sn_clip. This prevents overly aggressive rejection in high S/N ratio
        spectrum which neverthless differ at a level greater than the implied
        S/N due to
    lower : float, optional, default=3.0
        lower rejection threshold for djs_reject
    upper : float, optional, default=3.0
        upper rejection threshold for djs_reject
    maxrej : int, optional
        maximum number of pixels to reject in each iteration for djs_reject.
    debug_order_stack : bool, default=False
        show interactive QA plots for the order stacking
    debug_global_stack : bool, default=False
        show interactive QA plots for the global stacking of all the
        setups/orders/exposures
    debug : bool, optional, default=False
        show all QA plots useful for debugging. Note there are lots of QA plots,
        so only set this to True if you want to inspect them all.
    debug_scale : bool, optional, default=False
        show interactive QA plots for the rescaling of the spectra
    show_order_scale : bool, optional, default=False
        show interactive QA plots for the order scaling
    show : bool, optional, default=False,
        show coadded spectra QA plot or not
    show_exp : bool, optional, default=False
        show individual exposure spectra QA plot or not

    Returns
    -------
    wave_grid_mid : `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See
        core.wavecal.wvutils.py for more.  shape=(ngrid,)
    a_tuple : tuple
        Contains:

            - ``wave_giant_stack``: ndarray, (ngrid,): Wavelength grid for stacked
              spectrum. As discussed above, this is the weighted average of the
              wavelengths of each spectrum that contriuted to a bin in the input
              wave_grid wavelength grid. It thus has ngrid elements, whereas
              wave_grid has ngrid+1 elements to specify the ngrid total number
              of bins. Note that wave_giant_stack is NOT simply the wave_grid
              bin centers, since it computes the weighted average;

            - ``flux_giant_stack``: ndarray, (ngrid,): Final stacked spectrum on
              wave_stack wavelength grid;

            - ``ivar_giant_stack``: ndarray, (ngrid,): Inverse variance spectrum on
              wave_stack wavelength grid. Erors are propagated according to
              weighting and masking.;

            - ``mask_giant_stack``: ndarray, bool, (ngrid,): Mask for stacked
              spectrum on wave_stack wavelength grid. True=Good.

    b_tuple : tuple
        Contains:

            - ``waves_stack_orders``

            - ``fluxes_stack_orders``

            - ``ivars_stack_orders``

            - ``masks_stack_orders``

    """


    # Notes on object/list formats

    # waves_arr_setup  -- is a list of length nsetups, one for each setup. Each element is a numpy
    #                     array with shape = (nspec, norder, nexp) which is the data model for echelle spectra
    #                     for an individual setup. The utiltities utils.arr_setup_to_setup_list and
    #                     utils.setup_list_to_arr convert between arr_setup and setup_list

    #
    # waves_setup_list -- is a list of length nsetups, one for each setup. Each element of it is a list of length
    #                     norder*nexp elements, each of which contains the shape = (nspec1,) wavelength arrays
    #                     for the order/exposure in setup1. The list is arranged such that the nexp1 spectra
    #                     for iorder=0 appear first, then com nexp1 spectra for iorder=1, i.e. the outer or
    #                     fastest varying dimension in python array ordering is the exposure number. The utility
    #                     functions utils.echarr_to_echlist and utils.echlist_to_echarr convert between
    #                     the multi-d numpy arrays in the waves_arr_setup and the lists of numpy arrays in
    #                     waves_setup_list
    #
    # waves_concat     -- is a list of length = \Sum_i norder_i*nexp_i where the index i runs over the setups. The
    #                     elements of the list contains a numpy array of wavelengths for the
    #                     setup, order, exposure in question. The utility routines utils.setup_list_to_concat and
    #                     utils.concat_to_setup_list convert between waves_setup_lists and waves_concat


    if debug:
        show=True
        show_exp=True
        show_order_scale=True
        debug_order_stack=True
        debug_global_stack=True
        debug_scale=True


    if qafile is not None:
        qafile_stack = qafile.replace('.pdf', '_stack.pdf')
        qafile_chi = qafile.replace('.pdf', '_chi.pdf')
    else:
        qafile_stack = None
        qafile_chi = None

    # data shape
    nsetups=len(waves_arr_setup)

    if setup_ids is None:
        setup_ids = list(string.ascii_uppercase[:nsetups])

    setup_colors = utils.distinct_colors(nsetups)

    norders = []
    nexps = []
    nspecs = []
    for wave in waves_arr_setup:
        nspecs.append(wave.shape[0])
        norders.append(wave.shape[1])
        nexps.append(wave.shape[2])

    if nbests is None:
        _nbests = [int(np.ceil(norder/4)) for norder in norders]
    else:
        _nbests = nbests if isinstance(nbests, (list, np.ndarray)) else [nbests]*nsetups

    #  nspec, norder, nexp = shape
    # Decide how much to smooth the spectra by if this number was not passed in
    #nspec_good = []
    #ngood = []
    #if sn_smooth_npix is None:
    #    # Loop over setups
    #    for wave, norder, nexp in zip(waves_arr_setup, norders, nexps):
    #        # This is the effective good number of spectral pixels in the stack
    #        nspec_good.append(np.sum(wave > 1.0))
    #        ngood.append(norder*nexp)
    #    nspec_eff = np.sum(nspec_good)/np.sum(ngood)
    #    sn_smooth_npix = int(np.round(0.1 * nspec_eff))
    #    msgs.info('Using a sn_smooth_pix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    # Create the setup lists
    waves_setup_list = [utils.echarr_to_echlist(wave)[0] for wave in waves_arr_setup]
    fluxes_setup_list = [utils.echarr_to_echlist(flux)[0] for flux in fluxes_arr_setup]
    ivars_setup_list = [utils.echarr_to_echlist(ivar)[0] for ivar in ivars_arr_setup]
    gpms_setup_list= [utils.echarr_to_echlist(gpm)[0] for gpm in gpms_arr_setup]
    #shapes_setup_list = [wave.shape for wave in waves_arr_setup]

    # Generate some concatenated lists for wavelength grid determination
    waves_concat = utils.setup_list_to_concat(waves_setup_list)
    gpms_concat = utils.setup_list_to_concat(gpms_setup_list)

    # Generate a giant wave_grid
    wave_grid, wave_grid_mid, _ = wvutils.get_wave_grid(waves_concat, gpms=gpms_concat, wave_method=wave_method,
                                            wave_grid_min=wave_grid_min,
                                            wave_grid_max=wave_grid_max, dwave=dwave, dv=dv,
                                            dloglam=dloglam, spec_samp_fact=spec_samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    weights = []
    rms_sn_setup_list = []
    colors = []
    for isetup in range(nsetups):
        ## Need to feed have optional wave_min, wave_max here to deal with sources like high-z quasars?
        rms_sn_vec, _ = calc_snr(fluxes_setup_list[isetup], ivars_setup_list[isetup], gpms_setup_list[isetup])
        rms_sn = rms_sn_vec.reshape(norders[isetup], nexps[isetup])
        mean_sn_ord = np.mean(rms_sn, axis=1)
        best_orders = np.argsort(mean_sn_ord)[::-1][0:_nbests[isetup]]
        rms_sn_per_exp = np.mean(rms_sn[best_orders, :], axis=0)
        weights_exp = np.tile(np.square(rms_sn_per_exp), (nspecs[isetup], norders[isetup], 1))
        weights_isetup = weights_exp * weights_sens_arr_setup[isetup]
        weights.append(weights_isetup)
        rms_sn_setup_list.append(rms_sn)
        colors.append([setup_colors[isetup]]*norders[isetup]*nexps[isetup])



    # Create the waves_setup_list
    weights_setup_list = [utils.echarr_to_echlist(weight)[0] for weight in weights]

    if debug:
        weights_qa(utils.setup_list_to_concat(waves_setup_list),utils.setup_list_to_concat(weights_setup_list),
                   utils.setup_list_to_concat(gpms_setup_list), colors=utils.setup_list_to_concat(colors),
                   title='ech_combspec')

    #######################
    # Inter-order rescaling
    #######################
    #debug_scale=True
    fluxes_scl_interord_setup_list, ivars_scl_interord_setup_list, scales_interord_setup_list = [], [], []
    for isetup in range(nsetups):
        fluxes_scl_interord_isetup, ivars_scl_interord_isetup, scales_interord_isetup = [], [], []
        for iord in range(norders[isetup]):
            ind_start = iord*nexps[isetup]
            ind_end = (iord+1)*nexps[isetup]
            fluxes_scl_interord_iord, ivars_scl_interord_iord, scales_interord_iord, scale_method_used = \
                scale_spec_stack(wave_grid, wave_grid_mid, waves_setup_list[isetup][ind_start:ind_end],
                                 fluxes_setup_list[isetup][ind_start:ind_end],ivars_setup_list[isetup][ind_start:ind_end],
                                 gpms_setup_list[isetup][ind_start:ind_end], rms_sn_setup_list[isetup][iord, :],
                                 weights_setup_list[isetup][ind_start:ind_end], ref_percentile=ref_percentile,
                                 maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method,
                                 hand_scale=hand_scale,
                                 sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale, debug=debug_scale)
            fluxes_scl_interord_isetup += fluxes_scl_interord_iord
            ivars_scl_interord_isetup += ivars_scl_interord_iord
            scales_interord_isetup += scales_interord_iord

        fluxes_scl_interord_setup_list.append(fluxes_scl_interord_isetup)
        ivars_scl_interord_setup_list.append(ivars_scl_interord_isetup)
        scales_interord_setup_list.append(scales_interord_isetup)

    # TODO Add checking above in inter-order scaling such that orders with low S/N ratio are instead scaled using
    # scale factors from higher S/N ratio. The point is it makes no sense to take 0.0/0.0. In the low S/N regime,
    # i.e. DLAs, GP troughs, we should be rescaling using scale factors from orders with signal. This also applies
    # to the echelle combine below.


    #######################
    # Global Rescaling Computation -- Scale each setup/order/exp to match a preliminary global stack
    #######################
    #show_order_scale=True
    fluxes_concat = utils.setup_list_to_concat(fluxes_scl_interord_setup_list)
    ivars_concat = utils.setup_list_to_concat(ivars_scl_interord_setup_list)
    scales_concat = utils.setup_list_to_concat(scales_interord_setup_list)
    weights_concat = utils.setup_list_to_concat(weights_setup_list)
    rms_sn_concat = []
    for rms_sn in rms_sn_setup_list:
        rms_sn_concat += rms_sn.flatten().tolist()
    rms_sn_concat = np.array(rms_sn_concat)
    fluxes_pre_scale_concat = copy.deepcopy(fluxes_concat)
    ivars_pre_scale_concat  = copy.deepcopy(ivars_concat)
    # For the first iteration use the scale_method input as an argument (default=None, which will allow
    # soly_poly_ratio scaling which is very slow). For all the other iterations simply use median rescaling since
    # we are then applying tiny corrections and median scaling is much faster
    # ASC: This has been updated from:
    # scale_method_iter = [scale_method]  + ['median']*(niter_order_scale - 1)
    # to implement 2 median scaling prior to attempting other methods (particularly solve poly ratio). 
    # This may be helpful to correcting offsets between adjacent orders in multi-setup echelle observations. 
    # It was initially implemented with NIRSPEC in mind, but was not very effective in that case. We have tested
    # it against the X-Shooter dev suite example and found it does not harm the reduction, and may prove helpful in future. 
    scale_method_iter = ['median']*(2) + [scale_method]  + ['median']*(niter_order_scale - 3)
    # Iteratively scale and stack the entire set of spectra arcoss all setups, orders, and exposures
    for iteration in range(niter_order_scale):
        # JFH This scale_spec_stack routine takes a list of [(nspec1,), (nspec2,), ...] arrays, so a loop needs to be
        # added here over the outer setup dimension of the lists
        fluxes_scale_concat, ivars_scale_concat, scales_iter_concat, scale_method_used = scale_spec_stack(
            wave_grid, wave_grid_mid, waves_concat, fluxes_pre_scale_concat, ivars_pre_scale_concat, gpms_concat, rms_sn_concat,
            weights_concat, ref_percentile=ref_percentile,
            maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method_iter[iteration], hand_scale=hand_scale,
            sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale,
            show=(show_order_scale & (iteration == (niter_order_scale-1))))
        scales_concat = [scales_orig*scales_new for scales_orig, scales_new in zip(scales_concat, scales_iter_concat)]
        fluxes_pre_scale_concat = copy.deepcopy(fluxes_scale_concat)
        ivars_pre_scale_concat = copy.deepcopy(ivars_scale_concat)


    ###########################
    # Giant stack computation -- Perform the final stack with rejection for the globally rescaled spectra
    ###########################
    #debug=True
    wave_final_stack, flux_final_stack, ivar_final_stack, gpm_final_stack, nused_final_stack, out_gpms_concat = \
        spec_reject_comb(wave_grid, wave_grid_mid, waves_concat, fluxes_scale_concat, ivars_scale_concat, gpms_concat,
                         weights_concat, sn_clip=sn_clip, lower=lower, upper=upper, maxrej=maxrej,
                         maxiter_reject=maxiter_reject, debug=debug_global_stack)

    # Generate setup_lists and  arr_setup for some downstream computations
    fluxes_scale_setup_list = utils.concat_to_setup_list(fluxes_scale_concat, norders, nexps)
    ivars_scale_setup_list = utils.concat_to_setup_list(ivars_scale_concat, norders, nexps)
    out_gpms_setup_list = utils.concat_to_setup_list(out_gpms_concat, norders, nexps)
    fluxes_scale_arr_setup = utils.setup_list_to_arr_setup(fluxes_scale_setup_list, norders, nexps)
    ivars_scale_arr_setup = utils.setup_list_to_arr_setup(ivars_scale_setup_list, norders, nexps)
    out_gpms_arr_setup = utils.setup_list_to_arr_setup(out_gpms_setup_list, norders, nexps)


    ############################
    # Order Stack computation -- These are returned. CURRENTLY NOT USED FOR ANYTHING
    ############################
    #show_order_stacks=True
    waves_order_stack_setup, fluxes_order_stack_setup, ivars_order_stack_setup, gpms_order_stack_setup = [], [], [], []
    out_gpms_order_stack_setup_list = []
    for isetup in range(nsetups):
        waves_order_stack, fluxes_order_stack, ivars_order_stack, gpms_order_stack, out_gpms_order_stack = [], [], [], [], []
        for iord in range(norders[isetup]):
            ind_start = iord*nexps[isetup]
            ind_end = (iord+1)*nexps[isetup]
            wave_order_stack_iord, flux_order_stack_iord, ivar_order_stack_iord, gpm_order_stack_iord, \
                nused_order_stack_iord, outgpms_order_stack_iord = spec_reject_comb(
                wave_grid, wave_grid_mid, waves_setup_list[isetup][ind_start:ind_end],
                fluxes_scale_setup_list[isetup][ind_start:ind_end], ivars_scale_setup_list[isetup][ind_start:ind_end],
                gpms_setup_list[isetup][ind_start:ind_end], weights_setup_list[isetup][ind_start:ind_end],
                sn_clip=sn_clip, lower=lower, upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug_order_stack,
                title='order_stacks')
            waves_order_stack.append(wave_order_stack_iord)
            fluxes_order_stack.append(flux_order_stack_iord)
            ivars_order_stack.append(ivar_order_stack_iord)
            gpms_order_stack.append(gpm_order_stack_iord)
            out_gpms_order_stack.append(outgpms_order_stack_iord)
            if show_order_stacks:
                coadd_qa(wave_order_stack_iord, flux_order_stack_iord, ivar_order_stack_iord, nused_order_stack_iord,
                         gpm=gpm_order_stack_iord,
                         title='Coadded spectrum of order {:d}/{:d} for setup={:s}'.format(
                             iord, norders[isetup], setup_ids[isetup]))
        waves_order_stack_setup.append(waves_order_stack)
        fluxes_order_stack_setup.append(fluxes_order_stack)
        ivars_order_stack_setup.append(ivars_order_stack)
        gpms_order_stack_setup.append(gpms_order_stack)
        out_gpms_order_stack_setup_list.append(out_gpms_order_stack)

    ############################
    # QA Generation
    ############################
    if debug or show:
        fluxes_exps, ivars_exps, out_gpms_exps = np.array([],dtype=float), np.array([],dtype=float), np.array([],dtype=bool)
        flux_stack_exps, ivar_stack_exps, gpm_stack_exps = np.array([],dtype=float), np.array([],dtype=float), np.array([], dtype=bool)
        nrej_setup, norig_setup = [], []
        for isetup in range(nsetups):
            # Reshape everything now exposure-wise
            new_shape = (nspecs[isetup] * norders[isetup], nexps[isetup])
            waves_2d_exps = waves_arr_setup[isetup].reshape(new_shape, order='F')
            fluxes_2d_exps = fluxes_scale_arr_setup[isetup].reshape(new_shape, order='F')
            ivars_2d_exps = ivars_scale_arr_setup[isetup].reshape(new_shape, order='F')
            gpms_2d_exps = gpms_arr_setup[isetup].reshape(new_shape, order='F')
            out_gpms_2d_exps = out_gpms_arr_setup[isetup].reshape(new_shape, order='F')

            nrej = np.sum(np.logical_not(out_gpms_2d_exps) & gpms_2d_exps, axis=0)  # rejected pixels per exposure
            norig = np.sum((waves_2d_exps > 1.0) & np.logical_not(gpms_2d_exps), axis=0)  # originally masked pixels per exposure
            nrej_setup.append(np.sum(nrej))
            norig_setup.append(np.sum(norig))
            # Interpolate stack onto native 2d wavelength grids reshaped exposure-wise
            # JFH changed to wave_grid_mid
            flux_stack_2d_exps, ivar_stack_2d_exps, gpm_stack_2d_exps, _ = interp_spec(
                waves_2d_exps, wave_grid_mid, flux_final_stack, ivar_final_stack, gpm_final_stack)
            if show_exp:
                # Show QA plots for each exposure
                rejivars_2d_exps, sigma_corrs_2d_exps, outchi_2d_exps, gpm_chi_2d_exps = update_errors(
                    fluxes_2d_exps, ivars_2d_exps, out_gpms_2d_exps, flux_stack_2d_exps, ivar_stack_2d_exps,
                    gpm_stack_2d_exps, sn_clip=sn_clip)
                # QA for individual exposures
                for iexp in range(nexps[isetup]):
                    # plot the residual distribution
                    msgs.info('QA plots for exposure {:} with new_sigma = {:}'.format(iexp, sigma_corrs_2d_exps[iexp]))
                    # plot the residual distribution for each exposure
                    title_renorm = 'ech_combspec: Error distribution about stack for exposure {:d}/{:d} for setup={:s}'.format(iexp, nexps[isetup], setup_ids[isetup])
                    renormalize_errors_qa(outchi_2d_exps[:, iexp], gpm_chi_2d_exps[:, iexp], sigma_corrs_2d_exps[iexp],
                                          title=title_renorm)
                    title_coadd_iexp = 'ech_combspec: nrej={:d} pixels rejected,'.format(nrej[iexp]) + \
                                       ' norig={:d} originally masked,'.format(norig[iexp]) + \
                                       ' for exposure {:d}/{:d}'.format(iexp, nexps[isetup])
                    # JFH changed QA to use wave_grid_mid
                    coadd_iexp_qa(waves_2d_exps[:,iexp], fluxes_2d_exps[:,iexp], rejivars_2d_exps[:,iexp], gpms_2d_exps[:,iexp],
                                  wave_grid_mid, flux_final_stack, ivar_final_stack, gpm_final_stack, out_gpms_2d_exps[:, iexp],
                                  norder=norders[isetup], qafile=None, title=title_coadd_iexp)
            # Global QA for this setup
            rejivars_1d_isetup, sigma_corrs_1d_isetup, outchi_1d_isetup, gpm_chi_1d_isetup = update_errors(
                fluxes_2d_exps.flatten(), ivars_2d_exps.flatten(), out_gpms_2d_exps.flatten(),
                flux_stack_2d_exps.flatten(), ivar_stack_2d_exps.flatten(), gpm_stack_2d_exps.flatten(), sn_clip=sn_clip)
            renormalize_errors_qa(outchi_1d_isetup, gpm_chi_1d_isetup, sigma_corrs_1d_isetup[0],
                                  qafile=qafile_chi, title='Global Chi distribution for setup={:s}'.format(setup_ids[isetup]))
            fluxes_exps = np.append(fluxes_exps, fluxes_2d_exps.flatten())
            ivars_exps = np.append(ivars_exps, ivars_2d_exps.flatten())
            out_gpms_exps = np.append(out_gpms_exps, out_gpms_2d_exps.flatten())
            flux_stack_exps = np.append(flux_stack_exps,  flux_stack_2d_exps.flatten())
            ivar_stack_exps = np.append(ivar_stack_exps, ivar_stack_2d_exps.flatten())
            gpm_stack_exps = np.append(gpm_stack_exps, gpm_stack_2d_exps.flatten())

        # TODO Print out rejection statistics per setup?

        # Show a global chi distribution now for all setups, but only if there are multiple setups otherwise it is
        # identical to the isetup=0 plot.
        if nsetups > 1:
            rejivars_1d, sigma_corrs_1d, outchi_1d, gpm_chi_1d = update_errors(
                fluxes_exps, ivars_exps, out_gpms_exps, flux_stack_exps, ivar_stack_exps, gpm_stack_exps, sn_clip=sn_clip)
            renormalize_errors_qa(outchi_1d, gpm_chi_1d, sigma_corrs_1d[0], qafile=qafile_chi,
                                  title='Global Chi distribution for all nsetups={:d} setups'.format(nsetups))
        # show the final coadded spectrum
        coadd_qa(wave_grid_mid, flux_final_stack, ivar_final_stack, nused_final_stack, gpm=gpm_final_stack,
                 title='Final stacked spectrum', qafile=qafile_stack)


    return wave_grid_mid, (wave_final_stack, flux_final_stack, ivar_final_stack, gpm_final_stack), \
           (waves_order_stack_setup, fluxes_order_stack_setup, ivars_order_stack_setup, gpms_order_stack_setup,)


# ####################################################################
# ####################################################################
# Coadd2d routines follow this point
# ####################################################################

def get_wave_ind(wave_grid, wave_min, wave_max):
    """
    Utility routine used by coadd2d to determine the starting and ending indices of a wavelength grid.

    Parameters
    ----------
    wave_grid: `numpy.ndarray`_
          Wavelength grid.
    wave_min: float
          Minimum wavelength covered by the data in question.
    wave_max: float
          Maximum wavelength covered by the data in question.

    Returns
    -------
    ind_lower: int
        Integer lower indice corresponding to wave_min
    ind_upper: int
        Integer upper indice corresponding to wave_max
    """

    diff = wave_grid - wave_min
    diff[diff > 0] = np.inf
    if not np.any(diff < 0):
        ind_lower = 0
        msgs.warn('Your wave grid does not extend blue enough. Taking bluest point')
    else:
        ind_lower = np.argmin(np.abs(diff))
    diff = wave_max - wave_grid
    diff[diff > 0] = np.inf
    if not np.any(diff < 0):
        ind_upper = wave_grid.size-1
        msgs.warn('Your wave grid does not extend red enough. Taking reddest point')
    else:
        ind_upper = np.argmin(np.abs(diff))

    return ind_lower, ind_upper



def get_wave_bins(thismask_stack, waveimg_stack, wave_grid):
    """
    Utility routine to get the wavelength bins for 2d coadds from a mask

    Parameters
    ----------
    thismask_stack : list
        List of boolean arrays containing the masks indicating which pixels are
        on the slit in question.  `True` values are on the slit; `False` values
        are off the slit.  Length of the list is nimgs.   Shapes of the
        individual elements in the list are (nspec, nspat),  but each image can
        have a different shape.

    waveimg_stack : list
        List of the wavelength images, each of which is a float
        `numpy.ndarray`_. Length of the list is nimgs.  Shapes of the individual
        elements in the list are (nspec, nspat),  but each image can have a
        different shape.
    wave_grid : `numpy.ndarray`_, shape=(ngrid,)
        The wavelength grid created for the 2d coadd

    Returns
    -------
    wave_bins : `numpy.ndarray`_, shape = (ind_upper-ind_lower + 1, )
        Wavelength bins that are relevant given the illuminated pixels
        (thismask_stack) and wavelength coverage (waveimg_stack) of the image
        stack
    """

    # Determine the wavelength grid that we will use for the current slit/order
    # TODO This cut on waveimg_stack should not be necessary
    wave_lower = np.inf
    wave_upper = -np.inf
    for thismask, waveimg in zip(thismask_stack, waveimg_stack):
        wavemask = thismask & (waveimg > 1.0)
        wave_lower = min(wave_lower, np.amin(waveimg[wavemask]))
        wave_upper = max(wave_upper, np.amax(waveimg[wavemask]))
#        wave_min = waveimg[wavemask].min()
#        wave_max = waveimg[wavemask].max()
#        wave_lower = wave_min if wave_min < wave_lower else wave_lower
#        wave_upper = wave_max if wave_max > wave_upper else wave_upper
    ind_lower, ind_upper = get_wave_ind(wave_grid, wave_lower, wave_upper)
    return wave_grid[ind_lower:ind_upper + 1]


def get_spat_bins(thismask_stack, trace_stack, spat_samp_fact=1.0):
    """
    Determine the spatial bins for a 2d coadd and relative pixel coordinate
    images. This routine loops over all the images being coadded and creates an
    image of spatial pixel positions relative to the reference trace for each
    image in units of the desired rebinned spatial pixel sampling
    spat_samp_fact.  The minimum and maximum relative pixel positions in this
    frame are then used to define a spatial position grid with whatever desired
    pixel spatial sampling.

    Parameters
    ----------
    thismask_stack : list
        List of boolean arrays containing the masks indicating which pixels are
        on the slit in question.  `True` values are on the slit; `False` values
        are off the slit.  Length of the list is nimgs.   Shapes of the
        individual elements in the list are (nspec, nspat),  but each image can
        have a different shape.
    ref_trace_stack : list
        List of reference traces about which the images are rectified and
        coadded.  If the images were not dithered then this reference trace can
        simply be the center of the slit:

        .. code-block:: python

            slitcen = (slit_left + slit_righ)/2

        If the images were dithered, then this object can either be the slitcen
        appropriately shifted with the dither pattern, or it could be the trace
        of the object of interest in each exposure determined by running PypeIt
        on the individual images.  The list has nimgs elements, each of which is
        a 1D `numpy.ndarray`_ of shape (nspec,).
    spat_samp_fact : float, optional
        Spatial sampling for 2d coadd spatial bins in pixels. A value > 1.0
        (i.e. bigger pixels) will downsample the images spatially, whereas < 1.0
        will oversample. Default = 1.0

    Returns
    -------
    dspat_bins : `numpy.ndarray`_
        shape (spat_max_int +1 - spat_min_int,)
        Array of spatial bins for rectifying the image.
    dspat_stack : `numpy.ndarray`_
        shape (nimgs, nspec, nspat)
        Image stack which has the spatial position of each exposure relative to the trace in the trace_stack for that
        image.
    """

    nimgs = len(thismask_stack)
    dspat_stack = []
    spat_min = np.inf
    spat_max = -np.inf
    for thismask, trace in zip(thismask_stack, trace_stack):
        nspec, nspat = thismask.shape
        dspat_iexp = (np.arange(nspat)[np.newaxis, :] - trace[:, np.newaxis]) / spat_samp_fact
        dspat_stack.append(dspat_iexp)
        spat_min = min(spat_min, np.amin(dspat_iexp[thismask]))
        spat_max = max(spat_max, np.amax(dspat_iexp[thismask]))
#        spat_min_iexp = dspat_iexp[thismask].min()
#        spat_max_iexp = dspat_iexp[thismask].max()
#        spat_min = spat_min_iexp if spat_min_iexp < spat_min else spat_min
#        spat_max = spat_max_iexp if spat_max_iexp > spat_max else spat_max

    spat_min_int = np.floor(spat_min)
    spat_max_int = np.ceil(spat_max)
    dspat_bins = np.arange(spat_min_int, spat_max_int + 1.0, 1.0,dtype=float)
    return dspat_bins, dspat_stack

# TODO JFH I would like to modify this to take a stakc or coordinate or spatial
#  position image instead of a stack of reference traces
def compute_coadd2d(ref_trace_stack, sciimg_stack, sciivar_stack, skymodel_stack,
                    inmask_stack, thismask_stack, waveimg_stack,
                    wave_grid, spat_samp_fact=1.0, maskdef_dict=None,
                    weights=None, interp_dspat=True):
    """
    Construct a 2d co-add of a stack of PypeIt spec2d reduction outputs.

    Slits are 'rectified' onto a spatial and spectral grid, which
    encompasses the spectral and spatial coverage of the image stacks.
    The rectification uses nearest grid point interpolation to avoid
    covariant errors.  Dithering is supported as all images are centered
    relative to a set of reference traces in trace_stack.

    .. todo::
        These docs appear out-of-date

    Args:
        ref_trace_stack (list):
            List of reference traces about which the images are
            rectified and coadded.  If the images were not dithered then
            this reference trace can simply be the center of the slit::

                slitcen = (slit_left + slit_righ)/2

            If the images were dithered, then this object can either be
            the slitcen appropriately shifted with the dither pattern,
            or it could be the trace of the object of interest in each
            exposure determined by running PypeIt on the individual
            images.  Shape is (nspec, nimgs).
        sciimg_stack (list):
            List of science images, each of which is a float `numpy.ndarray`_. Length of the list is nimgs.
            Shapes of the individual elements in the list are (nspec, nspat),  but each image can have a different shape.
        sciivar_stack (list):
            List of inverse variance images, each of which is a float `numpy.ndarray`_.  Length of the list is nimgs.
            Shapes of the individual elements in the list are (nspec, nspat),  but each image can have a different shape.
        skymodel_stack (list):
            List of the model sky images, , each of which is a float `numpy.ndarray`_. Length of the list is nimgs.
            Shapes of the individual elements in the list are (nspec, nspat),  but each image can have a different shape.
        inmask_stack (list):
            List of input good pixel masks (i.e. `True` values are *good*, `False` values are *bad*.), each of which
            is a boolean `numpy.ndarray`_. Length of the list is nimgs.  Shapes of the individual elements in the list
            are (nspec, nspat),  but each image can have a different shape.
        waveimg_stack (list):
            List of the wavelength images, , each of which is a float `numpy.ndarray`_. Length of the list is nimgs.
            Shapes of the individual elements in the list are (nspec, nspat),  but each image can have a different shape.
        thismask_stack (list):
            List of boolean arrays containing the masks indicating which pixels are on
            the slit in question.  `True` values are on the slit;
            `False` values are off the slit.  Length of the list is nimgs.   Shapes of the individual elements in the list
            are (nspec, nspat),  but each image can have a different shape.
        weights (list, optional):
            The weights used when combining the rectified images (see
            :func:`~pypeit.core.combine.weighted_combine`).  If weights is None a
            uniform weighting is used.  If weights is not None then it must be a list of length nimgs.
            The individual elements of the list can either be floats, indiciating the weight to be used for each image, or
            arrays with  shape = (nspec,) or shape = (nspec, nspat), indicating pixel weights
            for the individual images. Weights are broadast to the correct size
            of the image stacks (see :func:`~pypeit.core.combine.broadcast_weights`), as necessary.
        spat_samp_fact (float, optional):
            Spatial sampling for 2d coadd spatial bins in pixels. A value > 1.0
            (i.e. bigger pixels) will downsample the images spatially, whereas <
            1.0 will oversample. Default = 1.0
        loglam_grid (`numpy.ndarray`_, optional):
            Wavelength grid in log10(wave) onto which the image stacks
            will be rectified.  The code will automatically choose the
            subset of this grid encompassing the wavelength coverage of
            the image stacks provided (see ``waveimg_stack``).
            Either `loglam_grid` or `wave_grid` must be provided.
        wave_grid (`numpy.ndarray`_, optional):
            Same as `loglam_grid` but in angstroms instead of
            log(angstroms). (TODO: Check units...)
        maskdef_dict (:obj:`dict`, optional): Dictionary containing all the maskdef info. The quantities saved
            are: maskdef_id, maskdef_objpos, maskdef_slitcen, maskdef_designtab. To learn what
            they are see :class:`~pypeit.slittrace.SlitTraceSet` datamodel.
        interp_dspat (bool, optional):
           Interpolate in the spatial coordinate image to faciliate running
           through core.extract.local_skysub_extract. This can be slow.   Default=True.



    Returns:
        dict: Returns a dict with the following keys:

            - wave_bins:
            - dspat_bins:
            - wave_mid:
            - wave_min:
            - wave_max:
            - dspat_mid:
            - sciimg: float ndarray shape = (nspec_coadd, nspat_coadd):
              Rectified and coadded science image
            - sciivar: float ndarray shape = (nspec_coadd, nspat_coadd):
              Rectified and coadded inverse variance image with correct
              error propagation
            - imgminsky: float ndarray shape = (nspec_coadd,
              nspat_coadd): Rectified and coadded sky subtracted image
            - outmask: bool ndarray shape = (nspec_coadd, nspat_coadd):
              Output mask for rectified and coadded images. True = Good,
              False=Bad.
            - nused: int ndarray shape = (nspec_coadd, nspat_coadd):
              Image of integers indicating the number of images from the
              image stack that contributed to each pixel
            - waveimg: float ndarray shape = (nspec_coadd, nspat_coadd):
              The averaged wavelength image corresponding to the
              rectified and coadded data.
            - dspat: float ndarray shape = (nspec_coadd, nspat_coadd):
              The average spatial offsets in pixels from the reference
              trace trace_stack corresponding to the rectified and
              coadded data.
            - nspec: int
            - nspat: int
            - maskdef_id: int
            - maskdef_slitcen: int
            - maskdef_objpos: int
            - maskdef_designtab: int

    """
    #nimgs, nspec, nspat = sciimg_stack.shape

    # TODO -- If weights is a numpy.ndarray, how can this not crash?
    #   Maybe the doc string above is inaccurate?

    nimgs =len(sciimg_stack)
    if weights is None:
        msgs.info('No weights were provided. Using uniform weights.')
        weights = (np.ones(nimgs)/float(nimgs)).tolist()

    shape_list = [sciimg.shape for sciimg in sciimg_stack]
    weights_stack = combine.broadcast_lists_of_weights(weights, shape_list)

    # Determine the wavelength grid that we will use for the current slit/order
    wave_bins = get_wave_bins(thismask_stack, waveimg_stack, wave_grid)
    dspat_bins, dspat_stack = get_spat_bins(thismask_stack, ref_trace_stack, spat_samp_fact=spat_samp_fact)

    skysub_stack = [sciimg - skymodel for sciimg, skymodel in zip(sciimg_stack, skymodel_stack)]
    #sci_list = [weights_stack, sciimg_stack, skysub_stack, tilts_stack,
    #            waveimg_stack, dspat_stack]
    sci_list = [weights_stack, sciimg_stack, skysub_stack, waveimg_stack, dspat_stack]

    var_list = [[utils.inverse(sciivar) for sciivar in sciivar_stack]]


    sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack \
            = rebin2d(wave_bins, dspat_bins, waveimg_stack, dspat_stack, thismask_stack,
                      inmask_stack, sci_list, var_list)
    # Now compute the final stack with sigma clipping
    sigrej = 3.0
    maxiters = 10
    # sci_list_rebin[0] = rebinned weights image stack
    # sci_list_rebin[1:] = stacks of images that we want to weighted combine
    # sci_list_rebin[2] = rebinned sciimg-sky_model images that we used for the sigma clipping
    # NOTE: outmask is a gpm
    sci_list_out, var_list_out, outmask, nused \
            = combine.weighted_combine(sci_list_rebin[0], sci_list_rebin[1:], var_list_rebin,
                               norm_rebin_stack != 0, sigma_clip=True,
                               sigma_clip_stack=sci_list_rebin[2], sigrej=sigrej,
                               maxiters=maxiters)
    sciimg, imgminsky, waveimg, dspat = sci_list_out
    sciivar = utils.inverse(var_list_out[0])

    # Compute the midpoints vectors, and lower/upper bins of the rectified image in spectral and spatial directions
    wave_mid = ((wave_bins + np.roll(wave_bins,1))/2.0)[1:]
    wave_min = wave_bins[:-1]
    wave_max = wave_bins[1:]
    dspat_mid = ((dspat_bins + np.roll(dspat_bins,1))/2.0)[1:]

    # Interpolate the dspat images wherever the coadds are masked
    # because a given pixel was not sampled. This is done because the
    # dspat image is not allowed to have holes if it is going to work
    # with local_skysub_extract
    nspec_coadd, nspat_coadd = imgminsky.shape
    spat_img_coadd, spec_img_coadd = np.meshgrid(np.arange(nspat_coadd), np.arange(nspec_coadd))

    if np.any(np.logical_not(outmask)) and interp_dspat:
        points_good = np.stack((spec_img_coadd[outmask], spat_img_coadd[outmask]), axis=1)
        points_bad = np.stack((spec_img_coadd[np.logical_not(outmask)],
                                spat_img_coadd[np.logical_not(outmask)]), axis=1)
        values_dspat = dspat[outmask]
        # JFH Changed to nearest on 5-26-20 because cubic is incredibly slow
        dspat_bad = scipy.interpolate.griddata(points_good, values_dspat, points_bad,
                                               method='nearest')
        dspat[np.logical_not(outmask)] = dspat_bad
        # Points outside the convex hull of the data are set to nan. We
        # identify those and simply assume them values from the
        # dspat_img_fake, which is what dspat would be on a regular
        # perfectly rectified image grid.
        nanpix = np.isnan(dspat)
        if np.any(nanpix):
            dspat_img_fake = spat_img_coadd + dspat_mid[0]
            dspat[nanpix] = dspat_img_fake[nanpix]
    else:
        dspat_img_fake = spat_img_coadd + dspat_mid[0]
        dspat[np.logical_not(outmask)] = dspat_img_fake[np.logical_not(outmask)]

    # initiate maskdef parameters
    # TODO I don't think this maskdef code belongs here. It should be rather be moved to the coadd2d class. This
    # is a core method that coadds images without any references to mask design tables
    maskdef_id = None
    maskdef_designtab = None
    new_maskdef_objpos = None
    new_maskdef_slitcen = None
    if maskdef_dict is not None and maskdef_dict['maskdef_id'] is not None:
        maskdef_id = maskdef_dict['maskdef_id']
        # update maskdef_objpos and maskdef_slitcen with the new value in the new slit
        if maskdef_dict['maskdef_objpos'] is not None and maskdef_dict['maskdef_slitcen'] is not None:
            new_maskdef_objpos = np.searchsorted(dspat[nspec_coadd//2, :], maskdef_dict['maskdef_objpos'])
            # maskdef_slitcen is the old slit center
            new_maskdef_slitcen = np.searchsorted(dspat[nspec_coadd//2, :], maskdef_dict['maskdef_slitcen'])
        if maskdef_dict['maskdef_designtab'] is not None:
            maskdef_designtab = maskdef_dict['maskdef_designtab']

    # TODO The rebin_stacks are indluded now for debugging but keeping them all may blow up memory usage so consider
    # removing
    return dict(wave_bins=wave_bins, dspat_bins=dspat_bins, wave_mid=wave_mid, wave_min=wave_min,
                wave_max=wave_max, dspat_mid=dspat_mid, sciimg=sciimg, sciivar=sciivar,
                imgminsky=imgminsky, outmask=outmask, nused=nused, waveimg=waveimg, # tilts=tilts,
                dspat=dspat, nspec=imgminsky.shape[0], nspat=imgminsky.shape[1],
                maskdef_id=maskdef_id, maskdef_slitcen=new_maskdef_slitcen, maskdef_objpos=new_maskdef_objpos,
                maskdef_designtab=maskdef_designtab)
#                rebin_weights_stack=sci_list_rebin[0], rebin_sciimg_stack=sci_list_rebin[1],
#                rebin_imgminsky_stack=sci_list_rebin[2], rebin_tilts_stack=sci_list_rebin[3],
#                rebin_waveimg_stack=sci_list_rebin[4], rebin_dspat_stack=sci_list_rebin[5],
#                rebin_var_stack=var_list_rebin[0], rebin_nsmp_stack=nsmp_rebin_stack,




def rebin2d(spec_bins, spat_bins, waveimg_stack, spatimg_stack,
            thismask_stack, inmask_stack, sci_list, var_list):
    """
    Rebin a set of images and propagate variance onto a new spectral and spatial grid. This routine effectively
    "recitifies" images using np.histogram2d which is extremely fast and effectively performs
    nearest grid point interpolation.

    Parameters
    ----------
    spec_bins : `numpy.ndarray`_, float, shape = (nspec_rebin)
        Spectral bins to rebin to.

    spat_bins : `numpy.ndarray`_, float ndarray, shape = (nspat_rebin)
        Spatial bins to rebin to.
    waveimg_stack : `numpy.ndarray`_, float , shape = (nimgs, nspec, nspat)
        Stack of nimgs wavelength images with shape = (nspec, nspat) each
    spatimg_stack : `numpy.ndarray`_, float, shape = (nimgs, nspec, nspat)
        Stack of nimgs spatial position images with shape = (nspec, nspat) each
    thismask_stack : `numpy.ndarray`_, bool, shape = (nimgs, nspec, nspat)
        Stack of nimgs images with shape = (nspec, nspat) indicating the
        locatons on the pixels on an image that are on the slit in question.
    inmask_stack : `numpy.ndarray`_, bool ndarray, shape = (nimgs, nspec, nspat)
        Stack of nimgs images with shape = (nspec, nspat) indicating which
        pixels on an image are masked.  True = Good, False = Bad
    sci_list : list
        Nested list of images, i.e. list of lists of images, where
        sci_list[i][j] is a shape = (nspec, nspat) where the shape can be
        different for each image. The ith index is the image type, i.e. sciimg,
        skysub, tilts, waveimg, the jth index is the exposure or image number,
        i.e. nimgs. These images are to be rebinned onto the commong grid.
    var_list : list
        Nested list of variance images, i.e. list of lists of images. The format
        is the same as for sci_list, but note that sci_list and var_list can
        have different lengths. Since this routine performs a NGP rebinning, it
        effectively comptues the average of a science image landing on a pixel.
        This means that the science is weighted by the 1/norm_rebin_stack, and
        hence variances must be weighted by that factor squared, which his why
        they must be input here as a separate list.

    Returns
    -------
    sci_list_out: list
        The list of ndarray rebinned images with new shape (nimgs, nspec_rebin,
        nspat_rebin)
    var_list_out : list
        The list of ndarray rebinned variance images with correct error
        propagation with shape (nimgs, nspec_rebin, nspat_rebin).
    norm_rebin_stack : int ndarray, shape (nimgs, nspec_rebin, nspat_rebin)
        An image stack indicating the integer occupation number of a given
        pixel. In other words, this number would be zero for empty bins, one for
        bins that were populated by a single pixel, etc. This image takes the
        input inmask_stack into account. The output mask for each image can be
        formed via outmask_rebin_stack = (norm_rebin_stack > 0).
    nsmp_rebin_stack : int ndarray, shape (nimgs, nspec_rebin, nspat_rebin)
        An image stack indicating the integer occupation number of a given pixel
        taking only the thismask_stack into account, but taking the inmask_stack
        into account. This image is mainly constructed for bookeeping purposes,
        as it represents the number of times each pixel in the rebin image was
        populated taking only the "geometry" of the rebinning into account (i.e.
        the thismask_stack), but not the masking (inmask_stack).
    """

    # allocate the output mages
    nimgs = len(sci_list[0])
    nspec_rebin = spec_bins.size - 1
    nspat_rebin = spat_bins.size - 1
    shape_out = (nimgs, nspec_rebin, nspat_rebin)
    nsmp_rebin_stack = np.zeros(shape_out)
    norm_rebin_stack = np.zeros(shape_out)
    sci_list_out = []
    for ii in range(len(sci_list)):
        sci_list_out.append(np.zeros(shape_out))
    var_list_out = []
    for jj in range(len(var_list)):
        var_list_out.append(np.zeros(shape_out))

    for img, (waveimg, spatimg, thismask, inmask) in enumerate(zip(waveimg_stack, spatimg_stack, thismask_stack, inmask_stack)):

        spec_rebin_this = waveimg[thismask]
        spat_rebin_this = spatimg[thismask]

        # This first image is purely for bookkeeping purposes to determine the number of times each pixel
        # could have been sampled
        nsmp_rebin_stack[img, :, :], spec_edges, spat_edges = np.histogram2d(spec_rebin_this, spat_rebin_this,
                                                               bins=[spec_bins, spat_bins], density=False)

        finmask = thismask & inmask
        spec_rebin = waveimg[finmask]
        spat_rebin = spatimg[finmask]
        norm_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                          bins=[spec_bins, spat_bins], density=False)
        norm_rebin_stack[img, :, :] = norm_img

        # Rebin the science images
        for indx, sci in enumerate(sci_list):
            weigh_sci, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                               bins=[spec_bins, spat_bins], density=False,
                                                               weights=sci[img][finmask])
            sci_list_out[indx][img, :, :] = (norm_img > 0.0) * weigh_sci/(norm_img + (norm_img == 0.0))

        # Rebin the variance images, note the norm_img**2 factor for correct error propagation
        for indx, var in enumerate(var_list):
            weigh_var, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                               bins=[spec_bins, spat_bins], density=False,
                                                               weights=var[img][finmask])
            var_list_out[indx][img, :, :] = (norm_img > 0.0)*weigh_var/(norm_img + (norm_img == 0.0))**2

    return sci_list_out, var_list_out, norm_rebin_stack.astype(int), nsmp_rebin_stack.astype(int)


