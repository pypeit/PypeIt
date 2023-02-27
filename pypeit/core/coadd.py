"""
Coadding module.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import os
import sys

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
    reduce sensitivity to outliers, ased on the scipy cookbook on robust optimization.

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
    ivarfit = mask_both/(1.0/(ivar_med + np.invert(mask_both)) + np.square(vmult)/(ivar_ref_med + np.invert(mask_both)))
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

    sigma_corr, maskchi = renormalize_errors(chi, mask=gpm, title = 'poly_ratio_fitfunc', debug=debug)
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


def interp_oned(wave_new, wave_old, flux_old, ivar_old, gpm_old, sensfunc=False):
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
        sensfunc (:obj:`bool`, optional):
            If True, the quantities ``flux*delta_wave`` and the corresponding
            ``ivar/delta_wave**2`` will be interpolated and returned instead of
            ``flux`` and ``ivar``. This is useful for sensitivity function
            computation where we need flux*(wavelength bin width). Beacause
            delta_wave is a difference of the wavelength grid, interpolating
            in the presence of masked data requires special care.

    Returns:
        :obj:`tuple`: Returns three `numpy.ndarray`_ objects with the
        interpolated flux, inverse variance, and good-pixel mask arrays with
        the length matching the new wavelength grid.
    """
    # Check input
    if wave_new.ndim != 1 or wave_old.ndim != 1:
        msgs.error('All input vectors must be 1D.')
    if flux_old.shape != wave_old.shape or ivar_old.shape != wave_old.shape \
            or gpm_old.shape != wave_old.shape:
        msgs.error('All vectors to interpolate must have the same size.')

    # Do not interpolate if the wavelength is exactly same with wave_new
    if np.array_equal(wave_new, wave_old):
        return flux_old, ivar_old, gpm_old

    wave_gpm = wave_old > 1.0 # Deal with the zero wavelengths
    if sensfunc:
        delta_wave_interp = wvutils.get_delta_wave(wave_old, wave_gpm)
        flux_interp = flux_old[wave_gpm]/delta_wave_interp[wave_gpm]
        ivar_interp = ivar_old[wave_gpm]*delta_wave_interp[wave_gpm]**2
    else:
        flux_interp = flux_old[wave_gpm]
        ivar_interp = ivar_old[wave_gpm]

    flux_new = scipy.interpolate.interp1d(wave_old[wave_gpm], flux_interp, kind='cubic',
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    ivar_new = scipy.interpolate.interp1d(wave_old[wave_gpm], ivar_interp, kind='cubic',
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    # Interpolate a floating-point version of the mask
    gpm_new_tmp = scipy.interpolate.interp1d(wave_old[wave_gpm], gpm_old.astype(float)[wave_gpm],
                                             kind='cubic', bounds_error=False,
                                             fill_value=np.nan)(wave_new)
    # Don't allow the ivar to be ever be less than zero
    ivar_new = (ivar_new > 0.0)*ivar_new
    gpm_new = (gpm_new_tmp > 0.8) & (ivar_new > 0.0) & np.isfinite(flux_new) & np.isfinite(ivar_new)
    return flux_new, ivar_new, gpm_new


# TODO: ``sensfunc`` should be something like "conserve_flux". It would be
# useful to compare these resampling routines against
# `pypeit.sampling.Resample`.
def interp_spec(wave_new, waves, fluxes, ivars, gpms, sensfunc=False):
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
            return interp_oned(wave_new, waves, fluxes, ivars, gpms, sensfunc=sensfunc)

        nexp = fluxes.shape[1]
        # Interpolate spectra to have the same wave grid with the iexp spectrum.
        # And scale spectra to the same flux level with the iexp spectrum.
        fluxes_inter = np.zeros((wave_new.size, nexp), dtype=float)
        ivars_inter = np.zeros((wave_new.size, nexp), dtype=float)
        gpms_inter = np.zeros((wave_new.size, nexp), dtype=bool)
        for ii in range(nexp):
            fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii] \
                    = interp_oned(wave_new, waves[:,ii], fluxes[:,ii], ivars[:,ii], gpms[:,ii],
                                  sensfunc=sensfunc)
        return fluxes_inter, ivars_inter, gpms_inter

    # Second case: interpolate a single spectrum onto an (nspec, nexp) array of
    # wavelengths. To make it here, wave_new.ndim must be 2.
    fluxes_inter = np.zeros_like(wave_new, dtype=float)
    ivars_inter = np.zeros_like(wave_new, dtype=float)
    gpms_inter = np.zeros_like(wave_new, dtype=bool)
    for ii in range(wave_new.shape[1]):
        fluxes_inter[:,ii], ivars_inter[:,ii], gpms_inter[:,ii] \
                = interp_oned(wave_new[:,ii], waves, fluxes, ivars, gpms, sensfunc=sensfunc)
    return fluxes_inter, ivars_inter, gpms_inter


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
    return sn_conv


def sn_weights(waves, fluxes, ivars, masks, sn_smooth_npix, const_weights=False,
               ivar_weights=False, relative_weights=False, verbose=False):

    """
    Calculate the S/N of each input spectrum and create an array of
    (S/N)^2 weights to be used for coadding.

    Parameters
    ----------
    waves : `numpy.ndarray`_
            Reference wavelength grid for all the spectra. If wave is a
            1d array the routine will assume that all spectra are on the
            same wavelength grid. If wave is a 2-d array, it will use
            the individual. shape = (nspec,) or (nspec, nexp)
    fluxes : `numpy.ndarray`_
            Stack of (nspec, nexp) spectra where nexp = number of
            exposures, and nspec is the length of the spectrum.
    ivars : `numpy.ndarray`_
            Inverse variance noise vectors for the spectra; shape = (nspec, nexp)
    masks : `numpy.ndarray`_
            Mask for stack of spectra. True=Good, False=Bad; shape = (nspec, nexp)
    sn_smooth_npix : float
            Number of pixels used for determining smoothly varying S/N ratio weights.
    const_weights : bool, optional
            Use a constant weights for each spectrum?
    ivar_weights : bool, optional
            Use inverse variance weighted scheme?
    relative_weights : bool, optional
            Calculate weights by fitting to the ratio of spectra? Note, relative weighting will
            only work well when there is at least one spectrum with a reasonable S/N, and a continuum.
            RJC note - This argument may only be better when the object being used has a strong
            continuum + emission lines. The reference spectrum is assigned a value of 1 for all
            wavelengths, and the weights of all other spectra will be determined relative to the
            reference spectrum. This is particularly useful if you are dealing with highly variable
            spectra (e.g. emission lines) and require a precision better than ~1 per cent.
    verbose : bool, optional
            Verbosity of print out.

    Returns
    -------
    rms_sn : `numpy.ndarray`_
        Root mean square S/N value for each input spectra; shape (nexp,)
    weights : `numpy.ndarray`_
        Weights to be applied to the spectra. These are
        signal-to-noise squared weights.  shape = (nspec, nexp)
    """

    # Give preference to ivar_weights
    if ivar_weights and relative_weights:
        msgs.warn("Performing inverse variance weights instead of relative weighting")
        relative_weights = False

    # Check input
    if fluxes.ndim == 1:
        nstack = 1
        nspec = fluxes.shape[0]
        wave_stack = waves.reshape((nspec, nstack))
        flux_stack = fluxes.reshape((nspec, nstack))
        ivar_stack = ivars.reshape((nspec, nstack))
        mask_stack = masks.reshape((nspec, nstack))
    elif fluxes.ndim == 2:
        nspec, nstack = fluxes.shape
        wave_stack = waves
        flux_stack = fluxes
        ivar_stack = ivars
        mask_stack = masks
    elif fluxes.ndim == 3:
        nspec, norder, nexp = fluxes.shape
        wave_stack = np.reshape(waves, (nspec, norder * nexp), order='F')
        flux_stack = np.reshape(fluxes, (nspec, norder * nexp), order='F')
        ivar_stack = np.reshape(ivars, (nspec, norder * nexp), order='F')
        mask_stack = np.reshape(masks, (nspec, norder * nexp), order='F')
        nstack = norder*nexp
    else:
        msgs.error('Unrecognized dimensionality for flux')

    # Calculate S/N
    sn_val = flux_stack*np.sqrt(ivar_stack)
    sn_val_ma = np.ma.array(sn_val, mask=np.logical_not(mask_stack))
    sn_sigclip = stats.sigma_clip(sn_val_ma, sigma=3, maxiters=5)
    # TODO: Update with sigma_clipped stats with our new cenfunc and std_func = mad_std
    sn2 = (sn_sigclip.mean(axis=0).compressed())**2  #S/N^2 value for each spectrum
    if sn2.shape[0] != nstack:
        msgs.error('No unmasked value in one of the exposures. Check inputs.')

    rms_sn = np.sqrt(sn2)  # Root Mean S/N**2 value for all spectra

    # Check if relative weights input
    if relative_weights:
        # Relative weights are requested, use the highest S/N spectrum as a reference
        ref_spec = np.argmax(sn2)
        if verbose:
            msgs.info(
                "The reference spectrum (ref_spec={0:d}) has a typical S/N = {1:.3f}".format(ref_spec, sn2[ref_spec]))
        # Adjust the arrays to be relative
        refscale = (sn_val[:, ref_spec] > 0) / (sn_val[:, ref_spec] + (sn_val[:, ref_spec] == 0))
        for iexp in range(nstack):
            # Compute the relative (S/N)^2 and update the mask
            sn2[iexp] /= sn2[ref_spec]
            mask_stack[:, iexp] = mask_stack[:, iexp] & ((mask_stack[:, ref_spec]) | (sn_val[:, ref_spec] != 0))
            sn_val[:, iexp] *= refscale

    # TODO: ivar weights is better than SN**2 or const_weights for merging orders. Eventually, we will change it to
    # TODO: Should ivar weights be deprecated??
    # Initialise weights
    weights = np.zeros_like(flux_stack)
    if ivar_weights:
        if verbose:
            msgs.info("Using ivar weights for merging orders")
        for iexp in range(nstack):
            weights[:, iexp] = smooth_weights(ivar_stack[:, iexp], mask_stack[:, iexp], sn_smooth_npix)
    else:
        for iexp in range(nstack):
            # Now
            if (rms_sn[iexp] < 3.0) or const_weights:
                weight_method = 'constant'
                weights[:, iexp] = np.full(nspec, np.fmax(sn2[iexp], 1e-2))  # set the minimum  to be 1e-2 to avoid zeros
            else:
                weight_method = 'wavelength dependent'
                # JFH THis line is experimental but it deals with cases where the spectrum drops to zero. We thus
                # transition to using ivar_weights. This needs more work because the spectra are not rescaled at this point.
                # RJC - also note that nothing should be changed to sn_val is relative_weights=True
                #sn_val[sn_val[:, iexp] < 1.0, iexp] = ivar_stack[sn_val[:, iexp] < 1.0, iexp]
                weights[:, iexp] = smooth_weights(sn_val[:, iexp]**2, mask_stack[:, iexp], sn_smooth_npix)
            if verbose:
                msgs.info('Using {:s} weights for coadding, S/N '.format(weight_method) +
                          '= {:4.2f}, weight = {:4.2f} for {:}th exposure'.format(
                              rms_sn[iexp], np.mean(weights[:, iexp]), iexp))

    if fluxes.ndim == 3:
        rms_sn = np.reshape(rms_sn, (norder, nexp), order='F')
        weights = np.reshape(weights, (nspec, norder, nexp), order='F')

    # Finish
    return rms_sn, weights


# TODO: This was commented out and would need to be refactored if brought back
# because of changes to the SensFunc and Telluric datamodels.
## TODO Rename this function to something sensfunc related
#def get_tell_from_file(sensfile, waves, masks, iord=None):
#    '''
#    Get the telluric model from the sensfile.
#
#    Args:
#        sensfile (str): the name of your fits format sensfile
#        waves (ndarray): wavelength grid for your output telluric model
#        masks (ndarray, bool): mask for the wave
#        iord (int or None): if None returns telluric model for all orders, otherwise return the order you want
#
#    Returns:
#         ndarray: telluric model on your wavelength grid
#    '''
#
#
#    sens_param = Table.read(sensfile, 1)
#    sens_table = Table.read(sensfile, 2)
#    telluric = np.zeros_like(waves)
#
#    if (waves.ndim == 1) and (iord is None):
#        msgs.info('Loading Telluric from Longslit sensfiles.')
#        tell_interp = scipy.interpolate.interp1d(sens_table[0]['WAVE'], sens_table[0]['TELLURIC'], kind='cubic',
#                                        bounds_error=False, fill_value=np.nan)(waves[masks])
#        telluric[masks] = tell_interp
#    elif (waves.ndim == 1) and (iord is not None):
#        msgs.info('Loading order {:} Telluric from Echelle sensfiles.'.format(iord))
#        wave_tell_iord = sens_table[iord]['WAVE']
#        tell_mask = (wave_tell_iord > 1.0)
#        tell_iord = sens_table[iord]['TELLURIC']
#        tell_iord_interp = scipy.interpolate.interp1d(wave_tell_iord[tell_mask], tell_iord[tell_mask], kind='cubic',
#                                        bounds_error=False, fill_value=np.nan)(waves[masks])
#        telluric[masks] = tell_iord_interp
#    else:
#        norder = np.shape(waves)[1]
#        for iord in range(norder):
#            wave_iord = waves[:, iord]
#            mask_iord = masks[:, iord]
#
#            # Interpolate telluric to the same grid with waves
#            # Since it will be only used for plotting, I just simply interpolate it rather than evaluate it based on the model
#            wave_tell_iord = sens_table[iord]['WAVE']
#            tell_mask = (wave_tell_iord > 1.0)
#            tell_iord = sens_table[iord]['TELLURIC']
#            tell_iord_interp = scipy.interpolate.interp1d(wave_tell_iord[tell_mask], tell_iord[tell_mask], kind='cubic',
#                                                    bounds_error=False, fill_value=np.nan)(wave_iord[mask_iord])
#            telluric[mask_iord, iord] = tell_iord_interp
#
#    return telluric


def robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None, ref_percentile=70.0, min_good=0.05,
                        maxiters=5, sigrej=3.0, max_factor=10.0, snr_do_not_rescale=1.0,
                        verbose=False):
    """
    Robustly determine the ratio between input spectrum flux and reference spectrum flux_ref. The code will perform
    best if the reference spectrum is chosen to be the higher S/N ratio spectrum, i.e. a preliminary stack that you want
    to scale each exposure to match. Note that the flux and flux_ref need to be on the same wavelength grid!!

    Parameters
    ----------
    flux: `numpy.ndarray`_
            spectrum that will be rescaled. shape=(nspec,)
    ivar: `numpy.ndarray`_
            inverse variance for the spectrum that will be rescaled.
            Same shape as flux
    flux_ref: `numpy.ndarray`_
            reference spectrum. Same shape as flux
    mask: `numpy.ndarray`_, optional
            boolean mask for the spectrum that will be rescaled. True=Good.
            If not input, computed from inverse variance
    ivar_ref: `numpy.ndarray`_, optional
            inverse variance of reference spectrum.
    mask_ref: `numpy.ndarray`_, optional
            Boolean mask for reference spectrum. True=Good. If not input, computed from inverse variance.
    ref_percentile: float, optional, default=70.0
            Percentile fraction used for selecting the minimum SNR cut from the reference spectrum. Pixels above this
            percentile cut are deemed the "good" pixels and are used to compute the ratio. This must be a number
            between 0 and 100.
    min_good: float, optional, default = 0.05
            Minimum fraction of good pixels determined as a fraction of the total pixels for estimating the median ratio
    maxiters: int, optional, default = 5
            Maximum number of iterations for astropy.stats.SigmaClip
    sigrej: float, optional, default = 3.0
            Rejection threshold for astropy.stats.SigmaClip
    max_factor: float, optional, default = 10.0,
            Maximum allowed value of the returned ratio
    snr_do_not_rescale: float, optional default = 1.0
            If the S/N ratio of the set of pixels (defined by upper ref_percentile in the reference spectrum) in the
            input spectrum have a median value below snr_do_not_rescale, median rescaling will not be attempted
            and the code returns ratio = 1.0. We also use this parameter to define the set of pixels (determined from
            the reference spectrum) to compare for the rescaling.

    Returns
    -------
    ratio: float
        the number that must be multiplied into flux in order to get it to match up with flux_ref
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

        flux_ref_ma = np.ma.MaskedArray(flux_ref, np.invert(calc_mask))
        flux_ref_clipped, lower, upper = sigclip(flux_ref_ma, masked=True, return_bounds=True)
        mask_ref_clipped = np.invert(flux_ref_clipped.mask)  # mask_stack = True are good values

        flux_ma = np.ma.MaskedArray(flux_ref, np.invert(calc_mask))
        flux_clipped, lower, upper = sigclip(flux_ma, masked=True, return_bounds=True)
        mask_clipped = np.invert(flux_clipped.mask)  # mask_stack = True are good values

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

def order_median_scale(waves, fluxes, ivars, masks, min_good=0.05, maxiters=5,
                       max_factor=10., sigrej=3, debug=False, show=False):
    '''
    Function to scaling different orders by the median S/N


    Args:
        waves (`numpy.ndarray`_): wavelength array of your spectra with the shape of (nspec, norder)
        fluxes (`numpy.ndarray`_): flux array of your spectra with the shape of (nspec, norder)
        ivars (`numpy.ndarray`_): ivar array of your spectra with the shape of (nspec, norder)
        masks (`numpy.ndarray`_): mask for your spectra with the shape of (nspec, norder)
        min_good (float, optional): minimum fraction of the total number of good pixels needed for estimate the median ratio
        maxiters (int or float, optional): maximum iterations for rejecting outliers
        max_factor (float, optional): maximum scale factor
        sigrej (float, optional): sigma used for rejecting outliers
        debug (bool, optional): if True show intermediate QA
        show (bool, optional): if True show the final QA

    Returns:
        tuple: (1) fluxes_new (`numpy.ndarray`_): re-scaled fluxes with the shape
        of (nspec, norder).  (2) ivars_new (`numpy.ndarray`_): re-scaled ivars
        with the shape of (nspec, norder) (3) order_ratios (`numpy.ndarray`_): an
        array of scale factor with the length of norder
    '''

    norder = np.shape(waves)[1]
    order_ratios = np.ones(norder)

    ## re-scale bluer orders to match the reddest order.
    # scaling spectrum order by order. We use the reddest order as the reference since slit loss in redder is smaller
    for ii in range(norder - 1):
        iord = norder - ii - 1
        wave_blue, flux_blue, ivar_blue, mask_blue = waves[:, iord-1], fluxes[:, iord-1],\
                                                     ivars[:, iord-1], masks[:, iord-1]

        wave_red_tmp, flux_red_tmp = waves[:, iord], fluxes[:, iord]*order_ratios[iord]
        ivar_red_tmp, mask_red_tmp = ivars[:, iord]*1.0/order_ratios[iord]**2, masks[:, iord]
        wave_mask = wave_red_tmp>1.0
        wave_red, flux_red, ivar_red, mask_red = wave_red_tmp[wave_mask], flux_red_tmp[wave_mask], \
                                                 ivar_red_tmp[wave_mask], mask_red_tmp[wave_mask],

        # interpolate iord-1 (bluer) to iord-1 (redder)
        flux_blue_inter, ivar_blue_inter, mask_blue_inter = interp_spec(wave_red, wave_blue, flux_blue, ivar_blue, mask_blue)

        npix_overlap = np.sum(mask_blue_inter & mask_red)
        percentile_iord = np.fmax(100.0 * (npix_overlap / np.sum(mask_red)-0.05), 10)

        mask_both = mask_blue_inter & mask_red
        snr_median_red = np.median(flux_red[mask_both]*np.sqrt(ivar_red[mask_both]))
        snr_median_blue = np.median(flux_blue_inter[mask_both]*np.sqrt(ivar_blue_inter[mask_both]))

        ## TODO: we set the SNR to be minimum of 300 to turn off the scaling but we need the QA plot
        ##       need to think more about whether we need to scale different orders, it seems make the spectra
        ##       much bluer than what it should be.
        if (snr_median_blue>300.0) & (snr_median_red>300.0):
            order_ratio_iord = robust_median_ratio(flux_blue_inter, ivar_blue_inter, flux_red, ivar_red, mask=mask_blue_inter,
                                                   mask_ref=mask_red, ref_percentile=percentile_iord, min_good=min_good,
                                                   maxiters=maxiters, max_factor=max_factor, sigrej=sigrej)
            order_ratios[iord - 1] = np.fmax(np.fmin(order_ratio_iord, max_factor), 1.0/max_factor)
            msgs.info('Scaled {}th order to {}th order by {:}'.format(iord-1, iord, order_ratios[iord-1]))
        else:
            if ii>0:
                order_ratios[iord - 1] = order_ratios[iord]
                msgs.warn('Scaled {}th order to {}th order by {:} using the redder order scaling '
                          'factor'.format(iord-1, iord, order_ratios[iord-1]))
            else:
                msgs.warn('The SNR in the overlapped region is too low or there is not enough overlapped pixels.'+ msgs.newline() +
                          'Median scale between order {:} and order {:} was not attempted'.format(iord-1, iord))

        if debug:
            plt.figure(figsize=(12, 8))
            plt.plot(wave_red[mask_red], flux_red[mask_red], 'k-', label='reference spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue],color='dodgerblue', lw=3, label='raw spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue]*order_ratios[iord-1], color='r',
                     alpha=0.5, label='re-scaled spectrum')
            ymin, ymax = get_ylim(flux_blue, ivar_blue, mask_blue)
            plt.ylim([ymin, ymax])
            plt.xlim([np.min(wave_blue[mask_blue]), np.max(wave_red[mask_red])])
            plt.legend()
            plt.xlabel('wavelength')
            plt.ylabel('Flux')
            plt.show()

    # Update flux and ivar
    fluxes_new = np.zeros_like(fluxes)
    ivars_new = np.zeros_like(ivars)
    for ii in range(norder):
        fluxes_new[:, ii] *= order_ratios[ii]
        ivars_new[:, ii] *= 1.0/order_ratios[ii]**2

    if show:
        plt.figure(figsize=(12, 8))
        ymin = []
        ymax = []
        for ii in range(norder):
            wave_stack_iord = waves[:, ii]
            flux_stack_iord = fluxes_new[:, ii]
            ivar_stack_iord = ivars_new[:, ii]
            mask_stack_iord = masks[:, ii]
            med_width = (2.0 * np.ceil(0.1 / 10.0 * np.size(wave_stack_iord[mask_stack_iord])) + 1).astype(int)
            flux_med, ivar_med = median_filt_spec(flux_stack_iord, ivar_stack_iord, mask_stack_iord, med_width)
            plt.plot(wave_stack_iord[mask_stack_iord], flux_med[mask_stack_iord], alpha=0.7)
            #plt.plot(wave_stack_iord[mask_stack_iord], flux_stack_iord[mask_stack_iord], alpha=0.5)
            # plt.plot(wave_stack_iord[mask_stack_iord],1.0/np.sqrt(ivar_stack_iord[mask_stack_iord]))
            ymin_ii, ymax_ii = get_ylim(flux_stack_iord, ivar_stack_iord, mask_stack_iord)
            ymax.append(ymax_ii)
            ymin.append(ymin_ii)
        plt.xlim([np.min(waves[masks]), np.max(waves[masks])])
        plt.ylim([-0.15*np.median(ymax), 1.5*np.median(ymax)])
        plt.xlabel('Wavelength ($\\rm\\AA$)')
        plt.ylabel('Flux')
        plt.show()

    return fluxes_new, ivars_new, order_ratios


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
    Multiple items: tuple
        (1) flux_scale: ndarray (nspec,) scaled spectrum; (2)
        ivar_scale: ndarray (nspec,) inverse variance for scaled
        spectrum; (3) scale: `numpy.ndarray`_ (nspec,) scale factor applied to
        the spectrum and inverse variance; (4) scale_method: str, method
        that was used to scale the spectra.
    """

    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0


    # Interpolate the reference spectrum onto the wavelengths of the spectrum that will be rescaled
    flux_ref_int, ivar_ref_int, mask_ref_int = interp_spec(wave, wave_ref, flux_ref, ivar_ref, mask_ref)

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
        scale_spec_qa(wave, flux, ivar, wave_ref, flux_ref, ivar_ref, scale, method_used, mask = mask, mask_ref=mask_ref,
                      title='Scaling Applied to the Data')

    return flux_scale, ivar_scale, scale, method_used


def compute_stack(wave_grid, waves, fluxes, ivars, masks, weights, min_weight=1e-8):
    '''
    Compute a stacked spectrum from a set of exposures on the specified wave_grid with proper treatment of
    weights and masking. This code uses np.histogram to combine the data using NGP and does not perform any
    interpolations and thus does not correlate errors. It uses wave_grid to determine the set of wavelength bins that
    the data are averaged on. The final spectrum will be on an ouptut wavelength grid which is not the same as wave_grid.
    The ouput wavelength grid is the weighted average of the individual wavelengths used for each exposure that fell into
    a given wavelength bin in the input wave_grid. This 1d coadding routine thus maintains the independence of the
    errors for each pixel in the combined spectrum and computes the weighted averaged wavelengths of each pixel
    in an analogous way to the 2d extraction procedure which also never interpolates to avoid correlating erorrs.

    Parameters
    ----------
    wave_grid: `numpy.ndarray`_
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the get_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.  shape=(ngrid +1,)
    waves: `numpy.ndarray`_
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
            shape=(nspec, nexp)
    fluxes: `numpy.ndarray`_
            fluxes for each exposure on the waves grid
            shape=(nspec, nexp)
    ivars: `numpy.ndarray`_
            Inverse variances for each exposure on the waves grid
            shape=(nspec, nexp)
    masks: `numpy.ndarray`_
            Boolean masks for each exposure on the waves grid. True=Good.
            shape=(nspec, nexp)
    weights: `numpy.ndarray`_
            Weights to be used for combining your spectra. These are computed using sn_weights
            shape=(nspec, nexp)
    min_weight: float, optional
            Minimum allowed weight for any individual spectrum

    Returns
    -------
    wave_stack: `numpy.ndarray`_
        Wavelength grid for stacked
        spectrum. As discussed above, this is the weighted average
        of the wavelengths of each spectrum that contriuted to a
        bin in the input wave_grid wavelength grid. It thus has
        ngrid elements, whereas wave_grid has ngrid+1 elements to
        specify the ngrid total number of bins. Note that
        wave_stack is NOT simply the wave_grid bin centers, since
        it computes the weighted average. shape=(ngrid,)
    flux_stack: `numpy.ndarray`_
        Final stacked spectrum on wave_stack wavelength grid
        shape=(ngrid,)
    ivar_stack: `numpy.ndarray`_
        Inverse variance spectrum on wave_stack wavelength grid.
        Errors are propagated according to weighting and masking.
        shape=(ngrid,)
    mask_stack: `numpy.ndarray`_
        Boolean Mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
        shape=(ngrid,)
    nused: `numpy.ndarray`_
        Number of exposures which contributed to
        each pixel in the wave_stack. Note that this is in general
        different from nexp because of masking, but also becuse of
        the sampling specified by wave_grid. In other words,
        sometimes more spectral pixels in the irregularly gridded
        input wavelength array waves will land in one bin versus
        another depending on the sampling.
        shape=(ngrid,)
    '''

    #mask bad values and extreme values (usually caused by extreme low sensitivity at the edge of detectors)
    ubermask = masks & (weights > 0.0) & (waves > 1.0) & (ivars > 0.0) & (utils.inverse(ivars)<1e10)
    waves_flat = waves[ubermask].flatten()
    fluxes_flat = fluxes[ubermask].flatten()
    ivars_flat = ivars[ubermask].flatten()
    vars_flat = utils.inverse(ivars_flat)
    weights_flat = weights[ubermask].flatten()

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
    mask_stack = (weights_total > min_weight) & (nused > 0.0)

    return wave_stack, flux_stack, ivar_stack, mask_stack, nused

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
                  outmask, norder=None, title='', qafile=None):
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

    """


    fig = plt.figure(figsize=(14, 8))
    spec_plot = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Get limits
    ymin, ymax = get_ylim(flux_stack, ivar_stack, mask_stack)

    # Plot spectrum
    rejmask = mask & np.invert(outmask)
    wave_mask = wave > 1.0
    wave_stack_mask = wave_stack > 1.0
    spec_plot.plot(wave[rejmask], flux[rejmask],'s',zorder=10,mfc='None', mec='r', label='rejected pixels')
    spec_plot.plot(wave[np.invert(mask)], flux[np.invert(mask)],'v', zorder=10, mfc='None', mec='orange',
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
        if (np.max(wave[mask]) > 9000.0):
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


def weights_qa(waves, weights, masks, title=''):
    '''
    Routine to make a QA plot for the weights used to compute a stacked spectrum.

    Parameters
    ----------
    wave: `numpy.ndarray`_
            wavelength array for spectra that went into a stack;
            shape=(nspec, nexp)
    weights: `numpy.ndarray`_
            (S/N)^2 weights for the exposures that went into a stack. This would have been computed by sn_weights
            shape=(nspec, nexp)
    mask: `numpy.ndarray`_
            Boolean array indicating pixels which were masked
            in each individual exposure which go into the stack.
            shape=(nspec, nexp)
    title: str, optional
            Title for the plot.
    '''

    if waves.ndim == 1:
        nstack = 1
        nspec = waves.shape[0]
        waves_stack = waves.reshape((nspec, nstack))
        weights_stack = weights.reshape((nspec, nstack))
        masks_stack = masks.reshape((nspec, nstack))
    elif waves.ndim == 2:
        nspec, nstack = waves.shape
        waves_stack = waves
        weights_stack = weights
        masks_stack = masks
    elif weights.ndim == 3:
        nspec, norder, nexp = waves.shape
        waves_stack = np.reshape(waves, (nspec, norder * nexp), order='F')
        weights_stack = np.reshape(weights, (nspec, norder * nexp), order='F')
        masks_stack = np.reshape(masks, (nspec, norder * nexp), order='F')
        nstack = norder*nexp
    else:
        msgs.error('Unrecognized dimensionality for waves')


    fig = plt.figure(figsize=(12, 8))
    for iexp in range(nstack):
        wave_mask = waves_stack[:, iexp] > 1.0
        plt.plot(waves_stack[wave_mask,iexp], weights_stack[wave_mask,iexp]*masks_stack[wave_mask,iexp])

    plt.xlim(waves_stack[(waves_stack > 1.0)].min(), waves_stack[(waves_stack > 1.0)].max())
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Weights')
    plt.title(title, fontsize=16, color='red')
    plt.show()

def coadd_qa(wave, flux, ivar, nused, mask=None, tell=None,
             title=None, qafile=None):
    '''
    Routine to make QA plot of the final stacked spectrum. It works for both longslit/mulitslit, coadded individual
    order spectrum of the Echelle data and the final coadd of the Echelle data.

    Parameters
    ----------
    wave: `numpy.ndarray`_
            one-d wavelength array of your spectrum;
            shape=(nspec,)
    flux: `numpy.ndarray`_
            one-d flux array of your spectrum;
            shape=(nspec,)
    ivar: `numpy.ndarray`_
            one-d ivar array of your spectrum;
            shape=(nspec,)
    nused: `numpy.ndarray`_
            how many exposures used in the stack for each pixel, the same size with flux
            shape=(nspec,)
    mask: `numpy.ndarray`_, optional
            boolean mask array for your spectrum;
            shape=(nspec,)
    tell: `numpy.ndarray`_, optional
            one-d telluric array for your spectrum; shape=(nspec,)
    title: str, optional
            plot title
    qafile: str, optional
           QA file name
    '''
    #TODO: This routine should take a parset

    if mask is None:
        mask = ivar > 0.0

    wave_mask = wave > 1.0
    wave_min = wave[wave_mask].min()
    wave_max = wave[wave_mask].max()
    fig = plt.figure(figsize=(12, 8))
    # plot how may exposures you used at each pixel
    # [left, bottom, width, height]
    num_plot =  fig.add_axes([0.10, 0.70, 0.80, 0.23])
    spec_plot = fig.add_axes([0.10, 0.10, 0.80, 0.60])
    num_plot.plot(wave[wave_mask],nused[wave_mask],drawstyle='steps-mid',color='k',lw=2)
    num_plot.set_xlim([wave_min, wave_max])
    num_plot.set_ylim([0.0, np.fmax(1.1*nused.max(), nused.max()+1.0)])
    num_plot.set_ylabel('$\\rm N_{EXP}$')
    num_plot.yaxis.set_major_locator(MaxNLocator(integer=True))
    num_plot.yaxis.set_minor_locator(NullLocator())

    # Plot spectrum
    spec_plot.plot(wave[wave_mask], flux[wave_mask], color='black', drawstyle='steps-mid',zorder=1,alpha=0.8, label='Single exposure')
    spec_plot.plot(wave[wave_mask], np.sqrt(utils.inverse(ivar[wave_mask])),zorder=2, color='red', alpha=0.7,
                   drawstyle='steps-mid', linestyle=':')

    # Get limits
    ymin, ymax = get_ylim(flux, ivar, mask)

    # Plot transmission
    if (np.max(wave[mask])>9000.0) and (tell is None):
        skytrans_file = data.get_skisim_filepath('atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ymax
        spec_plot.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5,zorder=11)
    elif tell is not None:
        scale = 0.8*ymax
        spec_plot.plot(wave[wave_mask], tell[wave_mask]*scale, drawstyle='steps-mid', color='m',alpha=0.5,zorder=11)

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


def spec_reject_comb(wave_grid, waves, fluxes, ivars, masks, weights, sn_clip=30.0, lower=3.0, upper=3.0,
                     maxrej=None, maxiter_reject=5, title='', debug=False,
                     verbose=False):
    """
    Routine for executing the iterative combine and rejection of a set of spectra to compute a final stacked spectrum.

    Parameters
    ----------
    wave_grid: `numpy.ndarray`_
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the get_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.
            shape=(ngrid +1,)
    waves: `numpy.ndarray`_
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
            shape=(nspec, nexp)
    fluxes: `numpy.ndarray`_
            fluxes for each exposure on the waves grid
    ivars: `numpy.ndarray`_
            Inverse variances for each exposure on the waves grid
    masks: `numpy.ndarray`_
            Boolean masks for each exposure on the waves grid. True=Good.
    weights: `numpy.ndarray`_
            Weights to be used for combining your spectra. These are computed using sn_weights
    sn_clip: float, optional, default=30.0
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
    lower: float, optional, default=3.0
            lower rejection threshold for djs_reject
    upper: float, optional, default=3.0
            upper rejection threshold for djs_reject
    maxrej: int, optional, default=None
            maximum number of pixels to reject in each iteration for djs_reject.
    maxiter_reject: int, optional, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
    title: str, optional
             Title for QA plot
    debug: bool, optional, default=False
            Show QA plots useful for debugging.
    verbose: bool, optional, default=False
            Level

    Returns
    -------
    wave_stack: `numpy.ndarray`_
        Wavelength grid for stacked
        spectrum. As discussed above, this is the weighted average
        of the wavelengths of each spectrum that contriuted to a
        bin in the input wave_grid wavelength grid. It thus has
        ngrid elements, whereas wave_grid has ngrid+1 elements to
        specify the ngrid total number of bins. Note that
        wave_stack is NOT simply the wave_grid bin centers, since
        it computes the weighted average.
        shape=(ngrid,)
    flux_stack: `numpy.ndarray`_
        Final stacked spectrum on
        wave_stack wavelength grid
        shape=(ngrid,)
    ivar_stack: `numpy.ndarray`_
        Inverse variance spectrum
        on wave_stack wavelength grid. Erors are propagated
        according to weighting and masking.
        shape=(ngrid,)
    mask_stack: `numpy.ndarray`_
        Boolean mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
        shape=(ngrid,)
    outmask: `numpy.ndarray`_
        Output bool mask with shape=(nspec, nexp)
        indicating which pixels are rejected in each exposure of
        the original input spectra after performing all of the
        iterations of combine/rejection
    nused: `numpy.ndarray`_
        Number of exposures which
        contributed to each pixel in the wave_stack. Note that
        this is in general different from nexp because of masking,
        but also becuse of the sampling specified by wave_grid. In
        other words, sometimes more spectral pixels in the
        irregularly gridded input wavelength array waves will land
        in one bin versus another depending on the sampling.
        shape=(ngrid,)

    """
    thismask = np.copy(masks)
    iter = 0
    qdone = False
    while (not qdone) and (iter < maxiter_reject):
        wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(
            wave_grid, waves, fluxes, ivars, thismask, weights)
        flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(
            waves, wave_stack, flux_stack, ivar_stack, mask_stack)
        rejivars, sigma_corrs, outchi, maskchi = update_errors(fluxes, ivars, thismask,
                                                               flux_stack_nat, ivar_stack_nat, mask_stack_nat,
                                                               sn_clip=sn_clip)
        thismask, qdone = pydl.djs_reject(fluxes, flux_stack_nat, outmask=thismask,inmask=masks, invvar=rejivars,
                                          lower=lower,upper=upper, maxrej=maxrej, sticky=False)
        iter += 1


    if (iter == maxiter_reject) & (maxiter_reject != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter_reject) + ' reached in spec_reject_comb')
    outmask = np.copy(thismask)

    # print out a summary of how many pixels were rejected
    nexp = waves.shape[1]
    nrej = np.sum(np.invert(outmask) & masks, axis=0)
    norig = np.sum((waves > 1.0) & np.invert(masks), axis=0)

    if verbose:
        for iexp in range(nexp):
            # nrej = pixels that are now masked that were previously good
            msgs.info("Rejected {:d} pixels in exposure {:d}/{:d}".format(nrej[iexp], iexp, nexp))

    # Compute the final stack using this outmask
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, outmask, weights)

    # Used only for plotting below
    if debug:
        # TODO Add a line here to optionally show the distribution of all pixels about the stack as we do for X-shooter.
        #flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(waves, wave_stack, flux_stack, ivar_stack, mask_stack)
        for iexp in range(nexp):
            # plot the residual distribution for each exposure
            title_renorm = title + ': Error distriution about stack for exposure {:d}/{:d}'.format(iexp,nexp)
            renormalize_errors_qa(outchi[:, iexp], maskchi[:, iexp], sigma_corrs[iexp], title=title_renorm)
            # plot the rejections for each exposures
            title_coadd_iexp = title + ': nrej={:d} pixels rejected,'.format(nrej[iexp]) + \
                               ' norig={:d} originally masked,'.format(norig[iexp]) + \
                               ' for exposure {:d}/{:d}'.format(iexp,nexp)
            coadd_iexp_qa(waves[:, iexp], fluxes[:, iexp], rejivars[:, iexp], masks[:, iexp], wave_stack, flux_stack,
                          ivar_stack, mask_stack, outmask[:, iexp], qafile=None, title=title_coadd_iexp)
        # weights qa
        title_weights = title + ': Weights Used -- nrej={:d} total pixels rejected,'.format(np.sum(nrej)) + \
                        ' norig={:d} originally masked'.format(np.sum(norig))
        weights_qa(waves, weights, outmask, title=title_weights)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused


def scale_spec_stack(wave_grid, waves, fluxes, ivars, masks, sn, weights,
                     ref_percentile=70.0, maxiter_scale=5,
                     sigrej_scale=3.0, scale_method='auto',
                     hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                     debug=False, show=False):

    """
    THIS NEEDS A PROPER DESCRIPTION

    Parameters
    ----------
    wave_grid: `numpy.ndarray`_
        New wavelength grid desired. This will typically be a reguarly spaced grid created by the get_wave_grid routine.
        The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
        bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
        lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
        with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
        too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
        limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
        to a given bin will be larger.
        shape=(ngrid +1,)
    waves: `numpy.ndarray`_
        wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
        each exposure and irregularly spaced.
        shape=(nspec, nexp)
    fluxes: `numpy.ndarray`_
        fluxes for each exposure on the waves grid
        shape=(nspec, nexp)
    ivars: `numpy.ndarray`_
        Inverse variances for each exposure on the waves grid
        shape=(nspec, nexp)
    masks: `numpy.ndarray`_
        Bool masks for each exposure on the waves grid. True=Good.
        shape=(nspec, nexp)
    sn:  `numpy.ndarray`_
        sn of each spectrum in the stack used to determine which scaling method should be used. This can
        be computed using sn_weights. shape=(nexp,)
    sigrej_scale: float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
    ref_percentile: float, optional, default=70.0
        percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
    maxiter_scale: int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    scale_method: str, optional, default='auto'
        Options are auto, poly, median, none, or hand. Hand is not well tested.
        User can optionally specify the rescaling method. Default is to let the
        code determine this automitically which works well.
    hand_scale: `numpy.ndarray`_, optional
        array of hand scale factors, not well tested
    sn_min_polyscale: float, optional, default=2.0
        maximum SNR for perforing median scaling
    sn_min_medscale: float, optional default=0.5
        minimum SNR for perforing median scaling
    debug: bool, optional, default=False
        show interactive QA plot

    Returns
    -------
    fluxes_scales: `numpy.ndarray`_
        Scale factors applied to the fluxes
        shape=(nspec, nexp)
    ivars_scales: `numpy.ndarray`_
        Scale factors applied to the ivars;
        shape=(nspec, nexp)
    scales: `numpy.ndarray`_
        shape=(nspec, nexp); Scale factors applied to
        each individual spectrum before the combine computed by
        scale_spec
    scale_method_used: list
        List of methods used for rescaling spectra.
    """

    # Compute an initial stack as the reference, this has its own wave grid based on the weighted averages
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, masks, weights)

    # Rescale spectra to line up with our preliminary stack so that we can sensibly reject outliers
    nexp = np.shape(fluxes)[1]
    fluxes_scale = np.zeros_like(fluxes)
    ivars_scale = np.zeros_like(ivars)
    scales = np.zeros_like(fluxes)
    scale_method_used = []
    for iexp in range(nexp):
        hand_scale_iexp = None if hand_scale is None else hand_scale[iexp]
        fluxes_scale[:, iexp], ivars_scale[:, iexp], scales[:, iexp], scale_method_iexp = scale_spec(
            waves[:, iexp], fluxes[:, iexp], ivars[:, iexp], sn[iexp], wave_stack, flux_stack, ivar_stack,
            mask=masks[:, iexp], mask_ref=mask_stack, ref_percentile=ref_percentile, maxiters=maxiter_scale,
            sigrej=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale_iexp, sn_min_polyscale=sn_min_polyscale,
            sn_min_medscale=sn_min_medscale, debug=debug, show=show)
        scale_method_used.append(scale_method_iexp)

    return fluxes_scale, ivars_scale, scales, scale_method_used


def combspec(waves, fluxes, ivars, masks, sn_smooth_npix,
             wave_method='linear', dwave=None, dv=None, dloglam=None,
             spec_samp_fact=1.0, wave_grid_min=None, wave_grid_max=None,
             ref_percentile=70.0, maxiter_scale=5, wave_grid_input=None,
             sigrej_scale=3.0, scale_method='auto', hand_scale=None,
             sn_min_polyscale=2.0, sn_min_medscale=0.5,
             const_weights=False, maxiter_reject=5, sn_clip=30.0,
             lower=3.0, upper=3.0, maxrej=None, qafile=None, title='', debug=False,
             debug_scale=False, show_scale=False, show=False,
             verbose=True):

    '''
    Main driver routine for coadding longslit/multi-slit spectra.
    NEED LOTS MORE HERE

    Parameters
    ----------
    waves: `numpy.ndarray`_
        Wavelength arrays for spectra to be stacked.
        shape=(nspec, nexp)
    fluxes: `numpy.ndarray`_
        Flux arrays for spectra to be stacked.
        shape=(nspec, nexp)
    ivars: `numpy.ndarray`_
        ivar arrays for spectra to be stacked.
        shape=(nspec, nexp)
    sn_smooth_npix: int
        Number of pixels to median filter by when computing S/N used to decide how to scale and weight spectra
    wave_method: str, optional
        method for generating new wavelength grid with get_wave_grid. Deafult is 'linear' which creates a uniformly
        space grid in lambda. See docuementation on get_wave_grid for description of the options.
    dwave: float, optional
        dispersion in units of A in case you want to specify it for get_wave_grid, otherwise the code computes the
        median spacing from the data.
    dv: float, optional
        Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option),
        otherwise a median value is computed from the data.
    spec_samp_fact: float, optional
        Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser (spec_samp_fact > 1.0) by this
        sampling factor. This basically multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
        spec_samp_fact are pixels.
    wave_grid_min: float, optional
        In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
    wave_grid_max: float, optional
        In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
    ref_percentile:
        percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
    maxiter_scale: int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    wave_grid_input : `numpy.ndarray`_
        User input wavelength grid to be used with the 'user_input' wave_method. Shape is (nspec_input,)
    maxiter_reject: int, optional, default=5
        maximum number of iterations for stacking and rejection. The code stops iterating either when
        the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
    sigrej_scale: float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
    scale_method: str, optional
        Options are auto, poly, median, none, or hand. Hand is not well tested.
        User can optionally specify the rescaling method.
        Default is 'auto' will let the code determine this automitically which works well.
    hand_scale: `numpy.ndarray`_, optional
            Array of hand scale factors, not well tested
    sn_min_polyscale: float, optional, default = 2.0
            maximum SNR for perforing median scaling
    sn_min_medscale: float, optional, default = 0.5
            minimum SNR for perforing median scaling
    const_weights: bool, optional
            If True, apply constant weight
    sn_clip: float, optional, default=30.0
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
    lower: float, optional, default=3.0
            lower rejection threshold for djs_reject
    upper: float: optional, default=3.0
            upper rejection threshold for djs_reject
    maxrej: int, optional
            maximum number of pixels to reject in each iteration for djs_reject.
    qafile: str, default=None
            Root name for QA, if None, it will be determined from the outfile
    title: str, optional
        Title for QA plots
    debug: bool, default=False
            Show all QA plots useful for debugging. Note there are lots of QA plots, so only set this to True if you want to inspect them all.
    debug_scale: bool, default=False
            show interactive QA plots for the rescaling of the spectra
    show: bool, default=False
        If True, show key QA plots or not
    show_scale: bool, default=False
        If True, show interactive QA plots for the rescaling of the spectra

    Returns
    -------
    wave_grid_mid: `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See core.wavecal.wvutils.py for more.
        shape=(ngrid,)
    wave_stack: `numpy.ndarray`_
        Wavelength grid for stacked
        spectrum. As discussed above, this is the weighted average
        of the wavelengths of each spectrum that contriuted to a
        bin in the input wave_grid wavelength grid. It thus has
        ngrid elements, whereas wave_grid has ngrid+1 elements to
        specify the ngrid total number of bins. Note that
        wave_stack is NOT simply the wave_grid bin centers, since
        it computes the weighted average.
        shape=(ngrid,)
    flux_stack: `numpy.ndarray`_
        Final stacked spectrum on
        wave_stack wavelength grid
        shape=(ngrid,)
    ivar_stack: `numpy.ndarray`_
        Inverse variance spectrum on wave_stack
        wavelength grid. Erors are propagated according to
        weighting and masking.
        shape=(ngrid,)
    mask_stack: `numpy.ndarray`_
        Boolean mask for stacked
        spectrum on wave_stack wavelength grid. True=Good.
        shape=(ngrid,)
    '''

    # We cast to float64 because of a bug in np.histogram
    waves = np.float64(waves)
    fluxes = np.float64(fluxes)
    ivars = np.float64(ivars)

    # Generate a giant wave_grid
    wave_grid, wave_grid_mid, _ = wvutils.get_wave_grid(
        waves=waves, masks = masks, wave_method=wave_method, 
        wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max, 
        wave_grid_input=wave_grid_input, 
        dwave=dwave, dv=dv, dloglam=dloglam, spec_samp_fact=spec_samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves, fluxes, ivars, masks, sn_smooth_npix, const_weights=const_weights, verbose=verbose)

    fluxes_scale, ivars_scale, scales, scale_method_used = scale_spec_stack(
        wave_grid, waves, fluxes, ivars, masks, rms_sn, weights, ref_percentile=ref_percentile, maxiter_scale=maxiter_scale,
        sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
        sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale, debug=debug_scale, show=show_scale)

    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves, fluxes_scale, ivars_scale, masks, weights, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug, title=title)

    if show:
        coadd_qa(wave_stack, flux_stack, ivar_stack, nused, mask=mask_stack, title='Stacked spectrum', qafile=qafile)

    return wave_grid_mid, wave_stack, flux_stack, ivar_stack, mask_stack

def multi_combspec(waves, fluxes, ivars, masks, sn_smooth_npix=None,
                   wave_method='linear', dwave=None, dv=None, dloglam=None, spec_samp_fact=1.0, wave_grid_min=None,
                   wave_grid_max=None, ref_percentile=70.0, maxiter_scale=5,
                   sigrej_scale=3.0, scale_method='auto', hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                   const_weights=False, maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0,
                   maxrej=None,
                   qafile=None, debug=False, debug_scale=False, show_scale=False, show=False):
    """
    Routine for coadding longslit/multi-slit spectra. Calls combspec which is
    the main stacking algorithm.

    Args:
        waves (`numpy.ndarray`_):
            Wavelength array  with shape (nspec, nexp) containing the spectra to be coadded.
        fluxes (`numpy.ndarray`_):
            Flux array with shape (nspec, nexp) containing the spectra to be coadded.
        ivars (`numpy.ndarray`_):
            Ivar array with shape (nspec, nexp) containing the spectra to be coadded.
        masks (`numpy.ndarray`_):
            Maks array with shape (nspec, nexp) containing the spectra to be coadded.
        sn_smooth_npix (int): optional
           Number of pixels to median filter by when computing S/N used to decide how to scale and weight spectra. If
           set to None, the code will determine the effective number of good pixels per spectrum
           in the stack that is being co-added and use 10% of this neff.
        wave_method (str, optional):
           method for generating new wavelength grid with get_wave_grid. Deafult is 'linear' which creates a uniformly
           space grid in lambda. See docuementation on get_wave_grid for description of the options.
        dwave (float, optional):
           dispersion in units of A in case you want to specify it for get_wave_grid, otherwise the code computes the
           median spacing from the data.
        dv (float, optional):
           Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option),
           otherwise a median value is computed from the data.
        spec_samp_fact (float, optional):
            Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser (spec_samp_fact > 1.0) by this
            sampling factor. This basically multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
            spec_samp_fact are pixels.
        wave_grid_min (float, optional):
           In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
        wave_grid_max (float, optional):
           In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
        wave_grid_input (`numpy.ndarray`_):
            User input wavelength grid to be used with the 'user_input' wave_method. Shape is (nspec_input,)
        maxiter_reject (int): optional
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached. Default=5.
        ref_percentile (float, optional):
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio. Should be a number between
            0 and 100, default = 70.0
        maxiter_scale (int, optional):
            Maximum number of iterations performed for rescaling spectra. Default=5.
        sigrej_scale (float, optional):
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec. Default=3.0
        scale_method (str, optional):
            Options are auto, poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default='auto' will let the
            code determine this automitically which works well.
        hand_scale (`numpy.ndarray`_, optional):
            Array of hand scale factors, not well tested
        sn_min_polyscale (float, optional):
            maximum SNR for perforing median scaling
        sn_min_medscale (float, optional):
            minimum SNR for perforing median scaling
        const_weights (`numpy.ndarray`_, optional):
             Constant weight factors specified
        maxiter_reject (int, optional):
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        sn_clip (float, optional):
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
        lower (float, optional):
            lower rejection threshold for djs_reject
        upper (float, optional):
            upper rejection threshold for djs_reject
        maxrej (int, optional):
            maximum number of pixels to reject in each iteration for djs_reject.
        nmaskedge (int, optional):
            Number of edge pixels to mask. This should be removed/fixed.
        qafile (str, optional): optional, default=None
            Root name for QA, if None, it will be determined from the outfile
        outfile (str, optional): optional, default=None,
            Root name for QA, if None, it will come from the target name from the fits header.
        debug (bool, optional): optinoal, default=False,
            Show all QA plots useful for debugging. Note there are lots of QA plots, so only set this to True if you want to inspect them all.
        debug_scale (bool, optional): optional, default=False
            show interactive QA plots for the rescaling of the spectra
        show (bool, optional): optional, default=False,
             Show key QA plots or not

    Returns:
        :obj:`tuple`:

            - wave_grid_mid: `numpy.ndarray`_, (ngrid,): Wavelength grid (in Angstrom)
              evaluated at the bin centers, uniformly-spaced either in lambda or
              log10-lambda/velocity.  See core.wavecal.wvutils.py for more.
            - wave_stack: `numpy.ndarray`_, (ngrid,): Wavelength grid for stacked
              spectrum. As discussed above, this is the weighted average of the
              wavelengths of each spectrum that contriuted to a bin in the input
              wave_grid wavelength grid. It thus has ngrid elements, whereas
              wave_grid has ngrid+1 elements to specify the ngrid total number
              of bins. Note that wave_stack is NOT simply the wave_grid bin
              centers, since it computes the weighted average.
            - flux_stack: `numpy.ndarray`_, (ngrid,): Final stacked spectrum on
              wave_stack wavelength grid
            - ivar_stack: `numpy.ndarray`_, (ngrid,): Inverse variance spectrum on
              wave_stack wavelength grid. Erors are propagated according to
              weighting and masking.
            - mask_stack: `numpy.ndarray`_, bool, (ngrid,): Mask for stacked spectrum on
              wave_stack wavelength grid. True=Good.

    """
    # Decide how much to smooth the spectra by if this number was not passed in
    if sn_smooth_npix is None:
        nspec, nexp = waves.shape
        # This is the effective good number of spectral pixels in the stack
        nspec_eff = np.sum(waves > 1.0)/nexp
        sn_smooth_npix = int(np.round(0.1*nspec_eff))
        msgs.info('Using a sn_smooth_npix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    wave_grid_mid, wave_stack, flux_stack, ivar_stack, mask_stack = combspec(
        waves, fluxes,ivars, masks, wave_method=wave_method, dwave=dwave, dv=dv, dloglam=dloglam,
        spec_samp_fact=spec_samp_fact, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max, ref_percentile=ref_percentile,
        maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
        sn_min_medscale=sn_min_medscale, sn_min_polyscale=sn_min_polyscale, sn_smooth_npix=sn_smooth_npix,
        const_weights=const_weights, maxiter_reject=maxiter_reject, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej,  qafile=qafile, title='multi_combspec', debug=debug, debug_scale=debug_scale, show_scale=show_scale,
        show=show)

    # Write to disk?
    #if outfile is not None:
    #    save.save_coadd1d_to_fits(outfile, wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
    #                              ex_value=ex_value, overwrite=True)

    return wave_grid_mid, wave_stack, flux_stack, ivar_stack, mask_stack


def ech_combspec(waves, fluxes, ivars, masks, weights_sens, nbest=None,
                 wave_method='log10', dwave=None, dv=None, dloglam=None,
                 spec_samp_fact=1.0, wave_grid_min=None, wave_grid_max=None,
                 ref_percentile=70.0, maxiter_scale=5, niter_order_scale=3,
                 sigrej_scale=3.0, scale_method='auto',
                 hand_scale=None, sn_min_polyscale=2.0, sn_min_medscale=0.5,
                 sn_smooth_npix=None, const_weights=False, maxiter_reject=5,
                 sn_clip=30.0, lower=3.0, upper=3.0,
                 maxrej=None, qafile=None, debug_scale=False, debug=False,
                 show_order_stacks=False, show_order_scale=False,
                 show_exp=False, show=False, verbose=False):
    """
    Driver routine for coadding Echelle spectra. Calls combspec which is the main stacking algorithm. It will deliver
    three fits files: spec1d_order_XX.fits (stacked individual orders, one order per extension), spec1d_merge_XX.fits
    (straight combine of stacked individual orders), spec1d_stack_XX.fits (a giant stack of all exposures and all orders).
    In most cases, you should use spec1d_stack_XX.fits for your scientific analyses since it reject most outliers.

    ..todo.. -- Clean up the doc formatting

    Parameters
    ----------
    waves: `numpy.ndarray`_
        Wavelength arrays for spectra to be stacked.
        shape=(nspec, norder, nexp)
    fluxes: `numpy.ndarray`_
        Flux arrays for spectra to be stacked.
        shape=(nspec, norder, nexp)
    ivars: `numpy.ndarray`_
        ivar arrays for spectra to be stacked.
        shape=(nspec, norder, nexp)
    masks: `numpy.ndarray`_
        Mask array with shape (nspec, norders, nexp) containing the spectra to be coadded.
    weights_sens: `numpy.ndarray`_
        Sensitivity function weights required for relatively weighting of the
        orders.  Must have the same shape as waves, etc.
    nbest: int, optional
        Number of orders to use for estimating the per exposure weights.
        Default is nbest=None, which will just use one fourth of the orders.
    wave_method: str, optional
        method for generating new wavelength grid with get_wave_grid. Deafult is 'log10' which creates a uniformly
        space grid in log10(lambda), which is typically the best for echelle spectrographs
    dwave: float, optional
        dispersion in units of A in case you want to specify it for get_wave_grid, otherwise the code computes the
        median spacing from the data.
    dv: float, optional
        Dispersion in units of km/s in case you want to specify it in the get_wave_grid  (for the 'velocity' option),
        otherwise a median value is computed from the data.
    dloglam: float, optional
        Dispersion in dimensionless units
    spec_samp_fact: float, optional
            Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser (spec_samp_fact > 1.0) by this
            sampling factor. This basically multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
            spec_samp_fact are pixels.
    wave_grid_min: float, optional
           In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
    wave_grid_max: float, optional
           In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
    wave_grid_input : `numpy.ndarray`_
            User input wavelength grid to be used with the 'user_input' wave_method. Shape is (nspec_input,)
    ref_percentile:
        percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
    maxiter_scale: int, optional, default=5
        Maximum number of iterations performed for rescaling spectra.
    sigrej_scale: float, optional, default=3.0
        Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
    hand_scale: `numpy.ndarray`_, optional
            Array of hand scale factors, not well tested
    sn_min_polyscale: float, optional, default = 2.0
            maximum SNR for perforing median scaling
    sn_min_medscale: float, optional, default = 0.5
            minimum SNR for perforing median scaling
    const_weights: bool, optional
            If True, apply constant weight
    maxiter_reject: int, optional, default=5
        maximum number of iterations for stacking and rejection. The code stops iterating either when
        the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
    sn_clip: float, optional, default=30.0
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
    lower: float, optional, default=3.0
            lower rejection threshold for djs_reject
    upper: float: optional, default=3.0
            upper rejection threshold for djs_reject
    maxrej: int, optional
            maximum number of pixels to reject in each iteration for djs_reject.

    Returns
    -------
    wave_grid_mid: `numpy.ndarray`_
        Wavelength grid (in Angstrom) evaluated at the bin centers,
        uniformly-spaced either in lambda or log10-lambda/velocity. See core.wavecal.wvutils.py for more.
        shape=(ngrid,)
    a_tuple: tuple
        (wave_giant_stack: ndarray, (ngrid,): Wavelength grid for
        stacked spectrum. As discussed above, this is the weighted
        average of the wavelengths of each spectrum that
        contriuted to a bin in the input wave_grid wavelength
        grid. It thus has ngrid elements, whereas wave_grid has
        ngrid+1 elements to specify the ngrid total number of
        bins. Note that wave_giant_stack is NOT simply the
        wave_grid bin centers, since it computes the weighted
        average;
        flux_giant_stack: ndarray, (ngrid,): Final stacked
        spectrum on wave_stack wavelength grid;
        ivar_giant_stack: ndarray, (ngrid,): Inverse variance
        spectrum on wave_stack wavelength grid. Erors are
        propagated according to weighting and masking.;
        mask_giant_stack: ndarray, bool, (ngrid,): Mask for
        stacked spectrum on wave_stack wavelength grid. True=Good.
        )
    another_tuple: tuple
        (waves_stack_orders,
        fluxes_stack_orders,
        ivars_stack_orders,
        masks_stack_orders,
        None)
    """

# TODO: Please leave this commented docstring entry here for now.
#        merge_stack: bool, default=False,
#            Compute an experimental combine of the high S/N combined orders in addition to the default algorithm,
#            which is to compute one giant stack using all order overlaps

    # output filenams for fits and QA plots
    #outfile_order = outfile.replace('.fits', '_order.fits') if outfile is not None else None

    if qafile is not None:
        qafile_stack = qafile.replace('.pdf', '_stack.pdf')
        qafile_chi = qafile.replace('.pdf', '_chi.pdf')
    else:
        qafile_stack = None
        qafile_chi = None

    # data shape
    nspec, norder, nexp = waves.shape
    if nbest is None:
        # Decide how many orders to use for estimating the per exposure weights. If nbest was not passed in
        # default to using one fourth of the orders
        nbest = int(np.ceil(norder/4))


    # Decide how much to smooth the spectra by if this number was not passed in
    if sn_smooth_npix is None:
        # This is the effective good number of spectral pixels in the stack
        nspec_eff = np.sum(waves > 1.0)/(norder*nexp)
        sn_smooth_npix = int(np.round(0.1 * nspec_eff))
        msgs.info('Using a sn_smooth_pix={:d} to decide how to scale and weight your spectra'.format(sn_smooth_npix))

    # create some arrays
    scales = np.zeros_like(waves)

    # Generate a giant wave_grid
    wave_grid, wave_grid_mid, _ = wvutils.get_wave_grid(waves, masks=masks, wave_method=wave_method,
                                            wave_grid_min=wave_grid_min,
                                            wave_grid_max=wave_grid_max, dwave=dwave, dv=dv,
                                            dloglam=dloglam, spec_samp_fact=spec_samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights_sn = sn_weights(waves, fluxes, ivars, masks, sn_smooth_npix, const_weights=const_weights, verbose=verbose)
    # Isolate the nbest best orders, and then use the average S/N of these to determine the per exposure relative weights.
    mean_sn_ord = np.mean(rms_sn, axis=1)
    best_orders = np.argsort(mean_sn_ord)[::-1][0:nbest]
    rms_sn_per_exp = np.mean(rms_sn[best_orders, :], axis=0)
    weights_exp = np.tile(rms_sn_per_exp**2, (nspec, norder, 1))
    weights = weights_exp*weights_sens
    #
    # Old code below for ivar weights if the sensfile was not passed in
    #msgs.error('Using ivar weights is deprecated.')
    #msgs.warn('No sensfunc is available for weighting, using smoothed ivar weights which is not optimal!')
    #_, weights_ivar = sn_weights(waves, fluxes, ivars, masks, sn_smooth_npix, const_weights=const_weights,
    # ivar_weights=True, verbose=True)
    #weights = weights_exp*weights_ivar

    if debug:
        weights_qa(waves, weights, masks, title='ech_combspec')

    fluxes_scl_interord = np.zeros_like(fluxes)
    ivars_scl_interord = np.zeros_like(ivars)
    scales_interord = np.zeros_like(fluxes)
    # First perform inter-order scaling once
    for iord in range(norder):
        # TODO Add checking here such that orders with low S/N ratio are instead scaled using scale factors from
        # higher S/N ratio. The point is it makes no sense to take 0.0/0.0. In the low S/N regime, i.e. DLAs,
        # GP troughs, we should be rescaling using scale factors from orders with signal. This also applies
        # to the echelle combine below.
        fluxes_scl_interord[:, iord], ivars_scl_interord[:,iord], scales_interord[:,iord], scale_method_used = \
            scale_spec_stack(wave_grid, waves[:, iord, :], fluxes[:, iord, :], ivars[:, iord, :], masks[:, iord, :],
                             rms_sn[iord, :], weights[:, iord, :], ref_percentile=ref_percentile,
                             maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method,
                             hand_scale=hand_scale,
                             sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale, debug=debug_scale)

    # Arrays to store rescaled spectra. Need Fortran like order reshaping to create a (nspec, norder*nexp) stack of spectra.
    # The order of the reshaping in the second dimension is such that blocks norder long for each exposure are stacked
    # sequentially, i.e. for order number [:, 0:norder] would be the 1st exposure, [:,norder:2*norder] would be the
    # 2nd exposure, etc.
    shape_2d = (nspec, norder * nexp)
    waves_2d = np.reshape(waves, shape_2d, order='F')
    fluxes_2d = np.reshape(fluxes_scl_interord, shape_2d, order='F')
    ivars_2d = np.reshape(ivars_scl_interord, shape_2d, order='F')
    masks_2d = np.reshape(masks, shape_2d, order='F')
    scales_2d = np.reshape(scales_interord, shape_2d, order='F')
    weights_2d = np.reshape(weights, shape_2d, order='F')
    rms_sn_2d = np.reshape(rms_sn, (norder*nexp), order='F')
    # Iteratively scale and stack the spectra, this takes or the order re-scaling we were doing previously
    fluxes_pre_scale = fluxes_2d.copy()
    ivars_pre_scale = ivars_2d.copy()
    # For the first iteration use the scale_method input as an argument (default=None, which will allow
    # soly_poly_ratio scaling which is very slow). For all the other iterations simply use median rescaling since
    # we are then applying tiny corrections and median scaling is much faster
    scale_method_iter = [scale_method]  + ['median']*(niter_order_scale - 1)
    for iter in range(niter_order_scale):
        fluxes_scale_2d, ivars_scale_2d, scales_iter, scale_method_used = scale_spec_stack(
            wave_grid, waves_2d, fluxes_pre_scale, ivars_pre_scale, masks_2d, rms_sn_2d, weights_2d, ref_percentile=ref_percentile,
            maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method_iter[iter], hand_scale=hand_scale,
            sn_min_polyscale=sn_min_polyscale, sn_min_medscale=sn_min_medscale,
            show=(show_order_scale & (iter == (niter_order_scale-1))))
        scales_2d *= scales_iter
        fluxes_pre_scale = fluxes_scale_2d.copy()
        ivars_pre_scale = ivars_scale_2d.copy()

    # Reshape the outputs to be (nspec, norder, nexp)
    fluxes_scale = np.reshape(fluxes_scale_2d, (nspec, norder, nexp), order='F')
    ivars_scale = np.reshape(ivars_scale_2d, (nspec, norder, nexp), order='F')
    scales = np.reshape(scales_2d, (nspec, norder, nexp), order='F')

    # Arrays to store stacked individual order spectra.
    waves_stack_orders = np.zeros((np.size(wave_grid) - 1, norder))
    fluxes_stack_orders = np.zeros_like(waves_stack_orders)
    ivars_stack_orders = np.zeros_like(waves_stack_orders)
    masks_stack_orders = np.zeros_like(waves_stack_orders, dtype=bool)
    outmasks_orders = np.zeros_like(masks)
    # Now perform stacks order by order
    for iord in range(norder):
        # Rejecting and coadding
        waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], \
        masks_stack_orders[:, iord],  outmasks_orders[:,iord,:], nused_iord = spec_reject_comb(
            wave_grid, waves[:, iord, :], fluxes_scale[:, iord, :], ivars_scale[:, iord, :], masks[:, iord, :], weights[:, iord, :],
            sn_clip=sn_clip, lower=lower, upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug,
            title='order_stacks')
        if show_order_stacks:
            # TODO This will probably crash since sensfile is not guarnetted to have telluric.
            #if sensfile is not None:
            #    tell_iord = get_tell_from_file(sensfile, waves_stack_orders[:, iord], masks_stack_orders[:, iord], iord=iord)
            #else:
            #    tell_iord = None
            tell_iord=None
            coadd_qa(waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], nused_iord,
                     mask=masks_stack_orders[:, iord], tell=tell_iord,
                     title='Coadded spectrum of order {:d}/{:d}'.format(iord, norder))

    # Now compute the giant stack
    wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmask_giant_stack, nused_giant_stack = \
        spec_reject_comb(wave_grid, waves_2d, fluxes_2d, ivars_2d, masks_2d, weights_2d, sn_clip=sn_clip,
                         lower=lower, upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug)

    # Reshape everything now exposure-wise
    waves_2d_exps = waves_2d.reshape((nspec * norder, nexp), order='F')
    fluxes_2d_exps = fluxes_2d.reshape(np.shape(waves_2d_exps), order='F')
    ivars_2d_exps = ivars_2d.reshape(np.shape(waves_2d_exps), order='F')
    masks_2d_exps = masks_2d.reshape(np.shape(waves_2d_exps), order='F')
    outmasks_2d_exps = outmask_giant_stack.reshape(np.shape(waves_2d_exps), order='F')
    # rejection statistics, exposure by exposure
    nrej = np.sum(np.invert(outmasks_2d_exps) & masks_2d_exps, axis=0)  # rejected pixels
    norig = np.sum((waves_2d_exps > 1.0) & np.invert(masks_2d_exps), axis=0) # originally masked pixels
    if debug or show:
        # Interpolate stack onto native 2d wavelength grids reshaped exposure-wise
        flux_stack_2d_exps, ivar_stack_2d_exps, mask_stack_2d_exps = interp_spec(
            waves_2d_exps, wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack)
        if show_exp:
            # Show QA plots for each exposure
            rejivars_2d_exps, sigma_corrs_2d_exps, outchi_2d_exps, maskchi_2d_exps = update_errors(
                fluxes_2d_exps, ivars_2d_exps, outmasks_2d_exps, flux_stack_2d_exps, ivar_stack_2d_exps,
                mask_stack_2d_exps, sn_clip=sn_clip)
            # QA for individual exposures
            for iexp in range(nexp):
                # plot the residual distribution
                msgs.info('QA plots for exposure {:} with new_sigma = {:}'.format(iexp, sigma_corrs_2d_exps[iexp]))
                # plot the residual distribution for each exposure
                title_renorm = 'ech_combspec: Error distribution about stack for exposure {:d}/{:d}'.format(iexp, nexp)
                renormalize_errors_qa(outchi_2d_exps[:, iexp], maskchi_2d_exps[:, iexp], sigma_corrs_2d_exps[iexp],
                                      title=title_renorm)
                title_coadd_iexp = 'ech_combspec: nrej={:d} pixels rejected,'.format(nrej[iexp]) + \
                                   ' norig={:d} originally masked,'.format(norig[iexp]) + \
                                   ' for exposure {:d}/{:d}'.format(iexp, nexp)
                coadd_iexp_qa(waves_2d_exps[:,iexp], fluxes_2d_exps[:,iexp], rejivars_2d_exps[:,iexp], masks_2d_exps[:,iexp],
                              wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmasks_2d_exps[:, iexp],
                              norder=norder, qafile=None, title=title_coadd_iexp)
        # Global QA
        rejivars_1d, sigma_corrs_1d, outchi_1d, maskchi_1d = update_errors(
            fluxes_2d_exps.flatten(), ivars_2d_exps.flatten(), outmasks_2d_exps.flatten(),
            flux_stack_2d_exps.flatten(), ivar_stack_2d_exps.flatten(), mask_stack_2d_exps.flatten(), sn_clip=sn_clip)
        renormalize_errors_qa(outchi_1d, maskchi_1d, sigma_corrs_1d[0], qafile=qafile_chi, title='Global Chi distribution')
        # show the final coadded spectrum
        coadd_qa(wave_giant_stack, flux_giant_stack, ivar_giant_stack, nused_giant_stack, mask=mask_giant_stack,
                 title='Final stacked spectrum', qafile=qafile_stack)

# TODO: Please leave this commented code in for now.
#    ## Stack with an altnernative method: combine the stacked individual order spectra directly. This is deprecated
#    merge_stack = False
#    if merge_stack:
#        order_weights = sensfunc_weights(sensfile, waves_stack_orders, debug=debug)
#        wave_merge, flux_merge, ivar_merge, mask_merge, nused_merge = compute_stack(
#            wave_grid, waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders, order_weights)
#        if show_order_stacks:
#            qafile_merge = 'spec1d_merge_{:}'.format(qafile)
#            coadd_qa(wave_merge, flux_merge, ivar_merge, nused_merge, mask=mask_merge, tell = None,
#                     title='Straight combined spectrum of the stacked individual orders', qafile=qafile_merge)
#        #if outfile is not None:
#        #    outfile_merge = outfile.replace('.fits', '_merge.fits')
#        #    save.save_coadd1d_to_fits(outfile_merge, wave_merge, flux_merge, ivar_merge, mask_merge, header=header,
#        #                              ex_value=ex_value, overwrite=True)

    # Save stacked individual order spectra
    #save.save_coadd1d_to_fits(outfile_order, waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
    #                          header=header, ex_value = ex_value, overwrite=True)
    #save.save_coadd1d_to_fits(outfile, wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack,
    #                          header=header, ex_value=ex_value, overwrite=True)


    return wave_grid_mid, (wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack), \
           (waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,)


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
    thismask_stack: list
        List of boolean arrays containing the masks indicating which pixels are on
        the slit in question.  `True` values are on the slit;
        `False` values are off the slit.  Length of the list is nimgs.   Shapes of the individual elements in the list
        are (nspec, nspat),  but each image can have a different shape.
    waveimg_stack: list
        List of the wavelength images, each of which is a float `numpy.ndarray`_. Length of the list is nimgs.
        Shapes of the individual elements in the list are (nspec, nspat),  but each image can have a different shape.
    wave_grid: array  shape (ngrid)
        The wavelength grid created for the 2d coadd

    Returns
    -------
    wave_bins : `numpy.ndarray`_
        shape (ind_upper-ind_lower + 1, )
        Wavelength bins that are relevant given the illuminated pixels (thismask_stack) and
        wavelength coverage (waveimg_stack) of the image stack

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
                    weights='uniform', interp_dspat=True):
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
        weights (list or str, optional):
            The weights used when combining the rectified images (see
            :func:`weighted_combine`).  If weights is set to 'uniform' then a
            uniform weighting is used.  Weights are broadast to the correct size
            of the image stacks (see :func:`broadcast_weights`), as necessary.
            If a list is passed it must have a length of nimgs. The individual elements
            of the list can either be floats, indiciating the weight to be used for each image, or
            arrays with  shape = (nspec,) or shape = (nspec, nspat), indicating pixel weights
            for the individual images. Weights are broadast to the correct size
            of the image stacks (see :func:`broadcast_weights`), as necessary.
            (TODO: JFH I think the str option should be changed here, but am leaving it for now.)
        spat_samp_fact (float, optional):
            Spatial sampling for 2d coadd spatial bins in pixels. A value > 1.0
            (i.e. bigger pixels) will downsample the images spatially, whereas <
            1.0 will oversample. Default = 1.0
        loglam_grid (`numpy.ndarray`_, optional):
            Wavelength grid in log10(wave) onto which the image stacks
            will be rectified.  The code will automatically choose the
            subset of this grid encompassing the wavelength coverage of
            the image stacks provided (see :func:`waveimg_stack`).
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
    if isinstance(weights,str) and weights == 'uniform':
        msgs.info('No weights were provided. Using uniform weights.')
        weights = np.ones(nimgs)/float(nimgs)

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
    spec_bins: `numpy.ndarray`_
        Spectral bins to rebin to.
        float, shape = (nspec_rebin)
    spat_bins: `numpy.ndarray`_
        Spatial bins to rebin to.
        float ndarray, shape = (nspat_rebin)
    waveimg_stack: `numpy.ndarray`_
        Stack of nimgs wavelength images with shape = (nspec, nspat) each
        float , shape = (nimgs, nspec, nspat)
    spatimg_stack: `numpy.ndarray`_
        Stack of nimgs spatial position images with shape = (nspec, nspat) each
        float, shape = (nimgs, nspec, nspat)
    thismask_stack: `numpy.ndarray`_
        Stack of nimgs images with shape = (nspec, nspat) indicating the locatons on the pixels on an image that
        are on the slit in question.
        bool, shape = (nimgs, nspec, nspat)
    inmask_stack: `numpy.ndarray`_
        Stack of nimgs images with shape = (nspec, nspat) indicating which pixels on an image are masked.
        True = Good, False = Bad
        bool ndarray, shape = (nimgs, nspec, nspat)

    sci_list: list
        Nested list of images, i.e. list of lists of images, where sci_list[i][j] is a shape = (nspec, nspat) where
        the shape can be different for each image. The ith index is the image type, i.e. sciimg, skysub, tilts, waveimg,
        the jth index is the exposure or image number, i.e. nimgs. These images are to be rebinned onto
        the commong grid.
    var_list: list
        Nested list of variance images, i.e. list of lists of images. The format is the same as for sci_list, but
        note that sci_list and var_list can have different lengths. Since this routine performs a NGP rebinning,
        it effectively comptues the average of a science image landing on a pixel. This means that the science
        is weighted by the 1/norm_rebin_stack, and hence variances must be weighted by that factor squared,
        which his why they must be input here as a separate list.


    Returns
    -------
    sci_list_out: list. The list of ndarray rebinned images
              with new shape (nimgs, nspec_rebin, nspat_rebin)
    var_list_out: list. The list of ndarray rebinned variance
              images with correct error propagation with shape (nimgs,
              nspec_rebin, nspat_rebin)
    norm_rebin_stack: int ndarray, shape (nimgs, nspec_rebin,
              nspat_rebin). An image stack indicating the integer
              occupation number of a given pixel. In other words, this
              number would be zero for empty bins, one for bins that
              were populated by a single pixel, etc. This image takes
              the input inmask_stack into account. The output mask for
              each image can be formed via outmask_rebin_satck =
              (norm_rebin_stack > 0)
    nsmp_rebin_stack: int ndarray, shape (nimgs, nspec_rebin,
              nspat_rebin). An image stack indicating the integer
              occupation number of a given pixel taking only the
              thismask_stack into account, but taking the inmask_stack
              into account. This image is mainly constructed for
              bookeeping purposes, as it represents the number of times
              each pixel in the rebin image was populated taking only
              the "geometry" of the rebinning into account (i.e. the
              thismask_stack), but not the masking (inmask_stack).
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

        # This fist image is purely for bookeeping purposes to determine the number of times each pixel
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
