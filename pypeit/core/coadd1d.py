
import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
from astropy import stats
from astropy.io import fits
from astropy import convolution

from pkg_resources import resource_filename
from pypeit import utils
from pypeit import msgs
from pypeit.core import load, save
from pypeit.core.wavecal import wvutils
from pypeit.core import pydl
from astropy import constants as const
c_kms = const.c.to('km/s').value

from matplotlib.ticker import NullFormatter, NullLocator

## Plotting parameters
plt.rcdefaults()
plt.rcParams['font.family'] = 'times new roman'
plt.rcParams["xtick.top"] = True
plt.rcParams["ytick.right"] = True
plt.rcParams["xtick.minor.visible"] = True
plt.rcParams["ytick.minor.visible"] = True
plt.rcParams["ytick.direction"] = 'in'
plt.rcParams["xtick.direction"] = 'in'
plt.rcParams["xtick.labelsize"] = 15
plt.rcParams["ytick.labelsize"] = 15
plt.rcParams["axes.labelsize"] = 17

# ToDo: update README and descriptions

def new_wave_grid(waves, wave_method='iref',iref=0, wave_grid_min=None, wave_grid_max=None,
                  A_pix=None,v_pix=None,samp_fact=1.0, **kwargs):
    """ Create a new wavelength grid for the spectra to be rebinned and coadded on

    Parameters
    ----------
    waves : masked ndarray
        Set of N original wavelength arrays
        nspec, nexp
    wave_method : str, optional
        Desired method for creating new wavelength grid.
        'iref' -- Use the first wavelength array (default)
        'velocity' -- Constant velocity
        'pixel' -- Constant pixel grid
        'concatenate' -- Meld the input wavelength arrays
    iref : int, optional
      Reference spectrum
    wave_grid_min: float, optional
      min wavelength value for the final grid
    wave_grid_max: float, optional
      max wavelength value for the final grid
    A_pix : float
      Pixel size in same units as input wavelength array (e.g. Angstroms)
      If not input, the median pixel size is calculated and used
    v_pix : float
      Pixel size in km/s for velocity method
      If not input, the median km/s per pixel is calculated and used
    samp_fact: float
      sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
      samp_fact < 1.0 undersamples (coarser)

    Returns
    -------
    wave_grid : ndarray
        New wavelength grid, not masked
    """

    #if not isinstance(waves, np.ma.MaskedArray):
    #    waves = np.ma.array(waves,mask=waves<1.0)
    wave_mask = waves>1.0

    if wave_method == 'velocity':  # Constant km/s
        spl = 299792.458
        if v_pix is None:
            # Find the median velocity of a pixel in the input
            dv = spl * np.abs(waves - np.roll(waves,1,axis=0)) / waves   # km/s
            v_pix = np.median(dv)

        # to make the wavelength grid finer or coarser
        v_pix = v_pix/samp_fact

        # Generate wavelength array
        if wave_grid_min is None:
            wave_grid_min = np.min(waves[wave_mask])
        if wave_grid_max is None:
            wave_grid_max = np.max(waves[wave_mask])
        x = np.log10(v_pix/spl + 1)
        npix = int(np.log10(wave_grid_max/wave_grid_min) / x) + 1
        wave_grid = wave_grid_min * 10.0**(x*np.arange(npix))

    elif wave_method == 'pixel': # Constant Angstrom
        if A_pix is None:
            dA =  np.abs(waves - np.roll(waves,1,axis=0))
            A_pix = np.median(dA)

        # Generate wavelength array
        if wave_grid_min is None:
            wave_grid_min = np.min(waves[wave_mask])
        if wave_grid_max is None:
            wave_grid_max = np.max(waves[wave_mask])
        wave_grid = wvutils.wavegrid(wave_grid_min, wave_grid_max + A_pix, \
                                     A_pix,samp_fact=samp_fact)

    elif wave_method == 'loggrid':
        dloglam_n = np.log10(waves) - np.roll(np.log10(waves), 1,axis=0)
        logwave_mask = wave_mask & np.roll(wave_mask, 1, axis=0)
        dloglam = np.median(dloglam_n[logwave_mask])
        wave_grid_max = np.max(waves[wave_mask])
        wave_grid_min = np.min(waves[wave_mask])
        loglam_grid = wvutils.wavegrid(np.log10(wave_grid_min), np.log10(wave_grid_max)+dloglam, \
                                       dloglam,samp_fact=samp_fact)
        wave_grid = 10.0**loglam_grid

    elif wave_method == 'concatenate':  # Concatenate
        # Setup
        loglam = np.log10(waves) # This deals with padding (0's) just fine, i.e. they get masked..
        nexp = waves.shape[1]
        newloglam = loglam[:, iref]  # Deals with mask
        # Loop
        for j in range(nexp):
            if j == iref:
                continue
            #
            iloglam = loglam[:, j]
            dloglam_0 = (newloglam[1]-newloglam[0])
            dloglam_n =  (newloglam[-1] - newloglam[-2]) # Assumes sorted
            if (newloglam[0] - iloglam[0]) > dloglam_0:
                kmin = np.argmin(np.abs(iloglam - newloglam[0] - dloglam_0))
                newloglam = np.concatenate([iloglam[:kmin], newloglam])
            #
            if (iloglam[-1] - newloglam[-1]) > dloglam_n:
                kmin = np.argmin(np.abs(iloglam - newloglam[-1] - dloglam_n))
                newloglam = np.concatenate([newloglam, iloglam[kmin:]])
        # Finish
        wave_grid = 10**newloglam

    elif wave_method == 'iref':
        wave_tmp = waves[:, iref]
        wave_grid = wave_tmp[wave_tmp>1.0]

    else:
        msgs.error("Bad method for scaling: {:s}".format(wave_method))

    return wave_grid

def gauss1(x, mean, sigma, area):
    ygauss = np.exp(-np.power(x - mean, 2.) / (2 * np.power(sigma, 2.)))
    norm = area / (sigma * np.sqrt(2 * np.pi))
    return norm * ygauss

def renormalize_errors_qa(chi, maskchi, sigma_corr, sig_range = 6.0, title='', qafile=None):

    n_bins = 50
    binsize = 2.0*sig_range/n_bins
    bins_histo = -sig_range + np.arange(n_bins)*binsize+binsize/2.0

    xvals = np.arange(-10.0,10,0.02)
    ygauss = gauss1(xvals,0.0,1.0,1.0)
    ygauss_new = gauss1(xvals,0.0,sigma_corr,1.0)

    plt.figure(figsize=(12, 8))
    plt.hist(chi[maskchi],bins=bins_histo,normed=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,ygauss,'c-',lw=3,label='sigma=1')
    plt.plot(xvals,ygauss_new,'m--',lw=2,label='new sigma={:4.2f}'.format(round(sigma_corr,2)))
    plt.xlabel('Residual distribution')
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

    return

def renormalize_errors(chi, mask, clip = 6.0, max_corr = 5.0, title = '', debug=False):

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
            msgs.warn("Error renormalization found sigma_corr/sigma = {:f} > {:f}." + msgs.newline() +
                      "Errors are severely underestimated." + msgs.newline() +
                      "Setting correction to sigma_corr = {:4.2f}".format(sigma_corr, max_corr, max_corr))
            sigma_corr = max_corr

        if debug:
            renormalize_errors_qa(chi, maskchi, sigma_corr, title=title)

    else:
        msgs.warn('No good pixels in error_renormalize. There are probably issues with your data')
        sigma_corr = 1.0

    return sigma_corr, maskchi

def poly_ratio_fitfunc_chi2(theta, flux_ref, thismask, arg_dict):
    """
    Function to be optimized for poly_ratio rescaling

    Args:
        theta:
        flux_ref:
        ivar_ref:
        thismask:
        arg_dict:

    Returns:

    """

    # Unpack the data to be rescaled, the mask for the reference spectrum, and the wavelengths
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    flux_med = arg_dict['flux_med']
    ivar_med = arg_dict['ivar_med']
    flux_ref_med = arg_dict['flux_ref_med']
    ivar_ref_med = arg_dict['ivar_ref_med']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = (utils.func_val(theta, wave, func, minx=wave_min, maxx=wave_max))**2
    flux_scale = ymult*flux_med
    mask_both = mask & thismask
    # This is the formally correct ivar used for the rejection, but not used in the fitting. This appears to yield
    # unstable results
    #totvar = utils.calc_ivar(ivar_ref) + ymult**2*utils.calc_ivar(ivar)
    #ivartot = mask_both*utils.calc_ivar(totvar)

    # The errors are rescaled at every function evaluation, but we only allow the errors to get smaller by up to a
    # factor of 1e4, and we only allow them to get larger slowly (as the square root).  This should very strongly
    # constrain the flux-corrrection vectors from going too small (or negative), or too large.
    ## Schlegel's version here
    vmult = np.fmax(ymult,1e-4)*(ymult <= 1.0) + np.sqrt(ymult)*(ymult > 1.0)
    ivarfit = mask_both/(1.0/(ivar_med + np.invert(mask_both)) + np.square(vmult)/(ivar_ref_med + np.invert(mask_both)))
    chi_vec = mask_both * (flux_ref_med - flux_scale) * np.sqrt(ivarfit)
    # Robustly characterize the dispersion of this distribution
    chi_mean, chi_median, chi_std = \
        stats.sigma_clipped_stats(chi_vec, np.invert(mask_both), cenfunc='median', stdfunc=stats.mad_std,
                                  maxiters=5, sigma=2.0)
    # The Huber loss function smoothly interpolates between being chi^2/2 for standard chi^2 rejection and
    # a linear function of residual in the outlying tails for large residuals. This transition occurs at the
    # value of the first argument, which we have set to be 2.0*chi_std, which is 2-sigma given the modified
    # errors described above from Schlegel's code.
    robust_scale = 2.0
    huber_vec = scipy.special.huber(robust_scale*chi_std, chi_vec)
    loss_function = np.sum(np.square(huber_vec*mask_both))
    #chi2 = np.sum(np.square(chi_vec))
    return loss_function

def poly_ratio_fitfunc(flux_ref, thismask, arg_dict, **kwargs_opt):

    # flux_ref, ivar_ref act like the 'data', the rescaled flux will be the 'model'

    # Function that we are optimizing
    #result = scipy.optimize.differential_evolution(poly_ratio_fitfunc_chi2, args=(flux_ref, ivar_ref, thismask, arg_dict,), **kwargs_opt)
    guess = arg_dict['guess']
    result = scipy.optimize.minimize(poly_ratio_fitfunc_chi2, guess, args=(flux_ref, thismask, arg_dict),  **kwargs_opt)
    #result = scipy.optimize.least_squares(poly_ratio_fitfunc_chi, guess, args=(flux_ref, thismask, arg_dict),  **kwargs_opt)
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    flux_med = arg_dict['flux_med']
    ivar_med = arg_dict['ivar_med']
    flux_ref_med = arg_dict['flux_ref_med']
    ivar_ref_med = arg_dict['ivar_ref_med']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = (utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max))**2
    flux_scale = ymult*flux
    mask_both = mask & thismask
    totvar = utils.calc_ivar(ivar_ref) + ymult**2*utils.calc_ivar(ivar)
    ivartot1 = mask_both*utils.calc_ivar(totvar)
    # Now rescale the errors
    chi = (flux_scale - flux_ref)*np.sqrt(ivartot1)
    try:
        debug = arg_dict['debug']
    except KeyError:
        debug = False

    sigma_corr, maskchi = renormalize_errors(chi, mask=thismask, title = 'poly_ratio_fitfunc', debug=debug)
    ivartot = ivartot1/sigma_corr**2

    return result, flux_scale, ivartot

def median_filt_spec(flux, ivar, mask, med_width):

    flux_med = np.zeros_like(flux)
    ivar_med = np.zeros_like(ivar)
    flux_med0 = utils.fast_running_median(flux[mask], med_width)
    flux_med[mask] = flux_med0
    var = utils.calc_ivar(ivar)
    var_med0 =  utils.smooth(var[mask], med_width)
    ivar_med[mask] = utils.calc_ivar(var_med0)
    return flux_med, ivar_med

def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None,
                     scale_min = 0.05, scale_max = 100.0, func='legendre',
                     maxiter=3, sticky=True, lower=3.0, upper=3.0, median_frac=0.01, debug=False):

    if mask is None:
        mask = (ivar > 0.0)
    if mask_ref is None:
        mask_ref = (ivar_ref > 0.0)

    #
    nspec = wave.size
    # Determine an initial guess
    ratio = robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=mask, mask_ref=mask_ref)
    guess = np.append(np.sqrt(ratio), np.zeros(norder-1))
    wave_min = wave.min()
    wave_max = wave.max()

    # Now compute median filtered versions of the spectra which we will actually operate on for the fitting. Note
    # that rejection will however work on the non-filtered spectra.
    med_width = (2.0*np.ceil(median_frac/2.0*nspec) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    flux_ref_med, ivar_ref_med = median_filt_spec(flux_ref, ivar_ref, mask_ref, med_width)

    arg_dict = dict(flux = flux, ivar = ivar, mask = mask,
                    flux_med = flux_med, ivar_med = ivar_med,
                    flux_ref_med = flux_ref_med, ivar_ref_med = ivar_ref_med,
                    ivar_ref = ivar_ref, wave = wave, wave_min = wave_min,
                    wave_max = wave_max, func = func, norder = norder, guess = guess, debug=debug)

    result, ymodel, ivartot, outmask = utils.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict, inmask=mask_ref,
                                                             maxiter=maxiter, lower=lower, upper=upper, sticky=sticky)
    ymult1 = (utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max))**2
    ymult = np.fmin(np.fmax(ymult1, scale_min), scale_max)
    flux_rescale = ymult*flux
    ivar_rescale = ivar/ymult**2

    if debug:
        # Determine the y-range for the QA plots
        scale_spec_qa(wave, flux_med, ivar_med, flux_ref_med, ivar_ref_med, ymult, 'poly', mask = mask, mask_ref=mask_ref,
                      title='Median Filtered Spectra that were poly_ratio Fit')

    return ymult, flux_rescale, ivar_rescale, outmask


def interp_oned(wave_new, wave_old, flux_old, ivar_old, mask_old):
    '''
    Args:
       wave_new: (one-D array) New wavelength
       wave_old: (one-D array) Old wavelength
       flux_old: (one-D array) Old flux
       ivar_old: (one-D array) Old ivar
       mask_old: (one-D array) Old float mask
    Returns :
       flux_new, ivar_new, mask_new (bool)
    '''

    # make the mask array to be float, used for interpolation
    masks_float = np.zeros_like(flux_old)
    masks_float[mask_old] = 1.0

    flux_new = interpolate.interp1d(wave_old[mask_old], flux_old[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    ivar_new = interpolate.interp1d(wave_old[mask_old], ivar_old[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    mask_new_tmp = interpolate.interp1d(wave_old[mask_old], masks_float[mask_old], kind='cubic', \
                                         bounds_error=False, fill_value=0.)(wave_new)
    mask_new = (mask_new_tmp > 0.5) & (ivar_new > 0.) & (flux_new != 0.)
    return flux_new,ivar_new,mask_new

def interp_spec(wave_new, waves, fluxes, ivars, masks):
    '''
    Interpolate all spectra to the page of wave_new
    Args:
        waves:
        fluxes:
        sigs:
        masks:
        iref:
    Returns:
    '''

    if (fluxes.ndim==2) and (wave_new.ndim==1):
        nexp = np.shape(fluxes)[1]
        # Interpolate spectra to have the same wave grid with the iexp spectrum.
        # And scale spectra to the same flux level with the iexp spectrum.
        fluxes_inter = np.zeros((wave_new.size, nexp))
        ivars_inter  = np.zeros((wave_new.size, nexp))
        masks_inter  = np.zeros((wave_new.size, nexp), dtype=bool)
        for ii in range(nexp):
            mask_ii = masks[:, ii]
            if np.sum(wave_new == waves[:, ii]) == np.size(wave_new):
                # do not interpolate if the wavelength is exactly same with wave_new
                fluxes_inter[:, ii] = fluxes[:, ii].copy()
                ivars_inter[:, ii] = ivars[:, ii].copy()
                masks_inter[:, ii] = mask_ii.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = \
                    interp_oned(wave_new, waves[:, ii],fluxes[:, ii],ivars[:, ii], masks[:, ii])
                fluxes_inter[:, ii] = flux_inter_ii  # * ratio_ii
                ivars_inter[:, ii] = ivar_inter_ii  # * ratio_ii
                masks_inter[:, ii] = mask_inter_ii

    elif (fluxes.ndim==1) and (wave_new.ndim==1):
        fluxes_inter, ivars_inter, masks_inter = interp_oned(wave_new,waves,fluxes,ivars,masks)

    elif (fluxes.ndim==1) and (wave_new.ndim==2):
        nexp = np.shape(wave_new)[1]
        fluxes_inter = np.zeros_like(wave_new)
        ivars_inter = np.zeros_like(wave_new)
        masks_inter = np.zeros_like(wave_new, dtype=bool)

        for ii in range(nexp):
            if np.sum(wave_new[:, ii] == waves) == np.size(waves):
                # do not interpolate if the wavelength is exactly same with wave_new
                fluxes_inter[:, ii] = fluxes.copy()
                ivars_inter[:, ii] = ivars.copy()
                masks_inter[:, ii] = masks.copy()
            else:
                flux_inter_ii, ivar_inter_ii, mask_inter_ii = \
                    interp_oned(wave_new[:, ii], waves, fluxes, ivars, masks)
                fluxes_inter[:, ii] = flux_inter_ii  # * ratio_ii
                ivars_inter[:, ii] = ivar_inter_ii  # * ratio_ii
                masks_inter[:, ii] = mask_inter_ii

    return fluxes_inter, ivars_inter, masks_inter

def sn_weights(waves, fluxes, ivars, masks, dv_smooth=10000.0, const_weights=False, sens_weights=False, verbose=False):
    """ Calculate the S/N of each input spectrum and create an array of (S/N)^2 weights to be used
    for coadding.

    Parameters
    ----------
    fluxes: float ndarray, shape = (nspec, nexp)
        Stack of (nspec, nexp) spectra where nexp = number of exposures, and nspec is the length of the spectrum.
    sigs: float ndarray, shape = (nspec, nexp)
        1-sigm noise vectors for the spectra
    masks: bool ndarray, shape = (nspec, nexp)
        Mask for stack of spectra. True=Good, False=Bad.
    waves: flota ndarray, shape = (nspec,) or (nspec, nexp)
        Reference wavelength grid for all the spectra. If wave is a 1d array the routine will assume
        that all spectra are on the same wavelength grid. If wave is a 2-d array, it will use the individual

    Optional Parameters:
    --------------------
    dv_smooth: float, 10000.0
         Velocity smoothing used for determining smoothly varying S/N ratio weights.

    Returns
    -------
    rms_sn : array
        Root mean square S/N value for each input spectra
    weights : ndarray
        Weights to be applied to the spectra. These are signal-to-noise squared weights.
    """

    sigs = np.sqrt(utils.calc_ivar(ivars))

    if fluxes.ndim == 1:
        nstack = 1
        nspec = fluxes.shape[0]
        flux_stack = fluxes.reshape((nspec, nstack))
        sig_stack = sigs.reshape((nspec, nstack))
        ivar_stack = ivars.reshape((nspec, nstack))
        mask_stack = masks.reshape((nspec, nstack))
    elif fluxes.ndim == 2:
        nspec = fluxes.shape[0]
        nstack = fluxes.shape[1]
        flux_stack = fluxes
        sig_stack = sigs
        ivar_stack = ivars
        mask_stack = masks
    else:
        msgs.error('Unrecognized dimensionality for flux')

    # if the wave
    if waves.ndim == 1:
        wave_stack = np.outer(waves, np.ones(nstack))
    elif waves.ndim == 2:
        wave_stack = waves
    else:
        msgs.error('wavelength array has an invalid size')

    # Calculate S/N
    sn_val = flux_stack*np.sqrt(ivar_stack)
    sn_val_ma = np.ma.array(sn_val, mask = np.invert(mask_stack))
    sn_sigclip = stats.sigma_clip(sn_val_ma, sigma=3, maxiters=5)
    ## TODO Update with sigma_clipped stats with our new cenfunc and std_func = mad_std
    sn2 = (sn_sigclip.mean(axis=0).compressed())**2 #S/N^2 value for each spectrum
    rms_sn = np.sqrt(sn2) # Root Mean S/N**2 value for all spectra
    rms_sn_stack = np.sqrt(np.mean(sn2))

    if sens_weights:
        # TODO: ivar weights is better than SN**2 or const_weights for merging orders. Enventially, we will change it to
        if verbose:
            msgs.info("Using sensfunc weights for merging orders")
        #weights = ivar_stack
        weights = np.ones_like(flux_stack) #((fluxes.shape[0], fluxes.shape[1]))
        spec_vec = np.arange(nspec)
        for iexp in range(nstack):
            imask = mask_stack[:, iexp]
            wave_now = wave_stack[imask, iexp]
            spec_now = spec_vec[imask]
            dwave = (wave_now - np.roll(wave_now,1))[1:]
            dv = (dwave/wave_now[1:])*c_kms
            dv_pix = np.median(dv)
            med_width = int(np.round(dv_smooth/dv_pix))
            sn_med1 = utils.fast_running_median(ivar_stack[imask,iexp], med_width)
            sn_med2 = np.interp(spec_vec, spec_now, sn_med1)
            sig_res = np.fmax(med_width/10.0, 3.0)
            gauss_kernel = convolution.Gaussian1DKernel(sig_res)
            sn_conv = convolution.convolve(sn_med2, gauss_kernel)
            weights[:, iexp] = sn_conv
    # ToDo: For high-z quasar, one would never want to use const_weight !!!
    elif rms_sn_stack <= 3.0 or const_weights:
    #if const_weights:
        weights = np.outer(np.ones(nspec), np.fmax(sn2,1e-5)) # set the minimum value to be 1e-5 to avoid zeros
        if verbose:
            msgs.info("Using constant weights for coadding, RMS S/N = {:g}".format(rms_sn_stack))
            for iexp in np.arange(nstack):
                msgs.info('S/N = {:4.2f}, weight = {:4.2f} for {:}th exposure'.format(
                    rms_sn[iexp],np.mean(weights[:,iexp]), iexp))
    else:
        if verbose:
            msgs.info("Using wavelength dependent weights for coadding")
        weights = np.ones_like(flux_stack) #((fluxes.shape[0], fluxes.shape[1]))
        spec_vec = np.arange(nspec)
        for iexp in range(nstack):
            imask = mask_stack[:, iexp]
            wave_now = wave_stack[imask, iexp]
            spec_now = spec_vec[imask]
            dwave = (wave_now - np.roll(wave_now,1))[1:]
            dv = (dwave/wave_now[1:])*c_kms
            dv_pix = np.median(dv)
            med_width = int(np.round(dv_smooth/dv_pix))
            sn_med1 = utils.fast_running_median(sn_val[imask,iexp]**2, med_width)
            #sn_med1 = utils.fast_running_median(ivar_stack[imask,iexp]**2, med_width)
            ##sn_med1 = scipy.ndimage.filters.median_filter(sn_val[imask,iexp]**2, size=med_width, mode='reflect')
            sn_med2 = np.interp(spec_vec, spec_now, sn_med1)
            ##sn_med2 = np.interp(wave_stack[iexp,:], wave_now,sn_med1)
            sig_res = np.fmax(med_width/10.0, 3.0)
            gauss_kernel = convolution.Gaussian1DKernel(sig_res)
            sn_conv = convolution.convolve(sn_med2, gauss_kernel)
            weights[:, iexp] = sn_conv
            if verbose:
                msgs.info('S/N = {:4.2f}, averaged weight = {:4.2f} for {:}th exposure'.format(
                    rms_sn[iexp],np.mean(weights[:, iexp]), iexp))
    # Finish
    return rms_sn, weights

def robust_median_ratio(flux, ivar, flux_ref, ivar_ref, ref_percentile=20.0, min_good=0.05, mask=None, mask_ref=None,
                        maxiters=5, max_factor = 10.0, sigrej = 3.0):
    '''
    Calculate the ratio between reference spectrum and your spectrum.
    Need to be in the same wave grid !!!
    Args:
        flux:
        sig:
        flux_ref:
        sig_ref:
        mask:
        mask_ref:
        snr_cut:
        maxiters:
        sigma:
    Returns:
        ratio_median: median ratio
    '''
    ## Mask for reference spectrum and your spectrum
    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0

    nspec = flux.size
    snr_ref = flux_ref * np.sqrt(ivar_ref)
    snr_ref_best = np.fmax(np.percentile(snr_ref[mask_ref], ref_percentile),0.5)
    calc_mask = (snr_ref > snr_ref_best) & mask_ref & mask

    if (np.sum(calc_mask) > min_good*nspec):
        # Take the best part of the higher SNR reference spectrum
        sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median', stdfunc=stats.mad_std)

        flux_ref_ma = np.ma.MaskedArray(flux_ref, np.invert(calc_mask))
        flux_ref_clipped, lower, upper = sigclip(flux_ref_ma, masked=True, return_bounds=True)
        mask_ref_clipped = np.invert(flux_ref_clipped.mask)  # mask_stack = True are good values

        flux_ma = np.ma.MaskedArray(flux_ref, np.invert(calc_mask))
        flux_clipped, lower, upper = sigclip(flux_ma, masked=True, return_bounds=True)
        mask_clipped = np.invert(flux_clipped.mask)  # mask_stack = True are good values

        new_mask = mask_ref_clipped & mask_clipped

        flux_ref_median = np.median(flux_ref[new_mask])
        flux_dat_median = np.median(flux[new_mask])

        #flux_ref_mean, flux_ref_median, flux_ref_std = \
        #    stats.sigma_clipped_stats(flux_ref,np.invert(calc_mask),cenfunc='median', stdfunc=stats.mad_std,
        #                              maxiters=maxiters,sigma=sigrej)
        #flux_dat_mean, flux_dat_median, flux_dat_std = \
        #    stats.sigma_clipped_stats(flux,np.invert(calc_mask),cenfunc='median', stdfunc=stats.mad_std,
        #                              maxiters=maxiters,sigma=sigrej)
        if (flux_ref_median < 0.0) or (flux_dat_median < 0.0):
            msgs.warn('Negative median flux found. Not rescaling')
            ratio = 1.0
        else:
            msgs.info('Used {:} good pixels for computing median flux ratio'.format(np.sum(new_mask)))
            ratio = np.fmax(np.fmin(flux_ref_median/flux_dat_median, max_factor), 1.0/max_factor)
    else:
        msgs.warn('Found only {:} good pixels for computing median flux ratio.'.format(np.sum(calc_mask))
                  + msgs.newline() + 'No median rescaling applied')
        ratio = 1.0

    return ratio

def order_median_scale(waves, fluxes, ivars, masks, min_good=0.05, maxiters=5, max_factor=10., sigrej=3,
                       debug=False, show=False):
    '''
    Args:
        waves: 2D array of wavelength
        fluxes: 2D array of flux
        ivars: 2D array of ivar
        masks: 2D array of mask
        min_good: minmum fraction of the total number of good pixels needed for estimate the median ratio
        maxiters: maximum iterations for rejecting outliers
        max_factor: maximum scale factor
        sigrej: sigma used for rejecting outliers
        debug: True or False
    Returns:
        re-scaled fluxes, ivars and an array of scale factor
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

        if (snr_median_blue>0.5) & (snr_median_red>0.5):
            order_ratio_iord = robust_median_ratio(flux_blue_inter, ivar_blue_inter, flux_red, ivar_red, mask=mask_blue_inter,
                                                   mask_ref=mask_red, ref_percentile=percentile_iord, min_good=min_good,
                                                   maxiters=maxiters, max_factor=max_factor, sigrej=sigrej)
            order_ratios[iord - 1] = np.fmax(np.fmin(order_ratio_iord, max_factor), 1.0/max_factor)
            msgs.info('Scaled {}th order to {}th order by {:}'.format(iord-1, iord, order_ratios[iord-1]))
        else:
            msgs.warn('The SNR in the overlapped region is too low or there is not enough overlapped pixels.'+ msgs.newline() +
                      'Median scale between order {:} and order {:} was not attempted'.format(iord-1, iord))

        if debug:
            plt.figure(figsize=(12, 8))
            plt.plot(wave_red[mask_red], flux_red[mask_red], 'k-', label='reference spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue],color='dodgerblue', lw=3, label='raw spectrum')
            plt.plot(wave_blue[mask_blue], flux_blue[mask_blue]*order_ratios[iord-1], color='r',
                     alpha=0.5, label='re-scaled spectrum')

            med_width = (2.0 * np.ceil(0.03 / 2.0 * wave_blue[mask_blue].size) + 1).astype(int)
            flux_med, ivar_med = median_filt_spec(flux_blue, ivar_blue, mask_blue, med_width)
            ymax = 1.5 * flux_med.max()
            ymin = -0.15 * ymax

            plt.ylim([ymin, ymax])
            plt.xlim([np.min(wave_blue[mask_blue]), np.max(wave_red[mask_red])])
            plt.legend()
            plt.xlabel('wavelength')
            plt.ylabel('Flux')
            plt.show()

    # Update flux and ivar
    fluxes_new = np.copy(fluxes)
    ivars_new = np.copy(ivars)
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
            plt.plot(wave_stack_iord[mask_stack_iord], flux_stack_iord[mask_stack_iord])
            # plt.plot(wave_stack_iord[mask_stack_iord],1.0/np.sqrt(ivar_stack_iord[mask_stack_iord]))

            med_width = (2.0 * np.ceil(0.1 / 2.0 * np.size(wave_stack_iord[mask_stack_iord])) + 1).astype(int)
            flux_med, ivar_med = median_filt_spec(flux_stack_iord, ivar_stack_iord, mask_stack_iord, med_width)
            ymax.append(1.5 * flux_med.max())
            ymin.append(-0.2 * flux_med.max())
        plt.xlim([np.min(waves[masks]), np.max(waves[masks])])
        plt.ylim([np.min(ymin), np.max(ymax)])
        plt.xlabel('Wavelength ($\\rm\\AA$)')
        plt.ylabel('Flux')
        plt.show()

    return fluxes_new, ivars_new, order_ratios

def scale_spec_qa(wave, flux, ivar, flux_ref, ivar_ref, ymult, scale_method,
                  mask=None, mask_ref=None, ylim = None, median_frac = 0.03, title=''):

    if mask is None:
        mask = ivar > 0.0
    if mask_ref  is None:
        mask_ref = ivar_ref > 0.0

    # This deals with spectrographs that have zero wavelength values. They are masked in mask, but this impacts plotting
    wave_mask = wave > 1.0
    # Get limits
    if ylim is None:
        nspec = flux.size
        med_width = (2.0 * np.ceil(median_frac / 2.0 * nspec) + 1).astype(int)
        flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
        flux_ref_med, ivar_ref_med = median_filt_spec(flux_ref, ivar_ref, mask_ref, med_width)
        flux_max = 1.5 * (np.fmax(flux_med.max(), flux_ref_med.max()))
        flux_min = np.fmin(-0.15 * flux_max, np.fmin(flux_med.min(), flux_ref_med.min()))
        ylim = (flux_min, flux_max)

    nullfmt = NullFormatter()  # no labels
    fig = plt.figure(figsize=(12, 8))
    # [left, bottom, width, height]
    poly_plot = fig.add_axes([0.1, 0.75, 0.8, 0.20])
    spec_plot = fig.add_axes([0.1, 0.1, 0.8, 0.65])
    poly_plot.xaxis.set_major_formatter(nullfmt)  # no x-axis labels for polynomial plot
    poly_plot.plot(wave[wave_mask], ymult[wave_mask], color='black', linewidth=3.0, label=scale_method + ' scaling')
    poly_plot.legend()
    spec_plot.set_ylim(ylim)
    spec_plot.plot(wave[wave_mask], flux_ref[wave_mask], color='black', drawstyle='steps-mid', zorder=3, label='reference spectrum')
    spec_plot.plot(wave[wave_mask], flux[wave_mask], color='dodgerblue', drawstyle='steps-mid', zorder=10, alpha=0.5,
                   label='original spectrum')
    spec_plot.plot(wave[wave_mask], flux[wave_mask]*ymult[wave_mask], color='red', drawstyle='steps-mid', alpha=0.7, zorder=1, linewidth=2,
                   label='rescaled spectrum')
    spec_plot.legend()
    fig.suptitle(title)
    plt.show()


def scale_spec(wave, flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None, min_good=0.05,
               ref_percentile=20.0, maxiters=5, sigrej=3, max_median_factor=10,
               npoly = None, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, debug=True):
    '''
    Scale the spectra into the same page with the reference spectrum.

    Args:
        waves:
        fluxes:
        ivars:
        masks:
        flux_iref:
        ivar_iref:
        mask_iref:
        iref:
        cenfunc:
        snr_cut:
        maxiters:
        sigma:
    Returns:
    '''

    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0

    # estimates the SNR of each spectrum and the stacked mean SNR
    rms_sn, weights = sn_weights(wave, flux, ivar, mask)
    rms_sn_stack = np.sqrt(np.mean(rms_sn**2))

    if scale_method is None:
        if rms_sn_stack > sn_max_medscale:
            scale_method = 'poly'
        elif ((rms_sn_stack <= sn_max_medscale) and (rms_sn_stack > sn_min_medscale)):
            scale_method = 'median'
        else:
            scale_method = 'none'

    # Estimate the scale factor
    if scale_method == 'poly':
        # Decide on the order of the polynomial rescaling
        if npoly is None:
            if rms_sn_stack > 25.0:
                npoly = 5 # Is this stable?
            elif rms_sn_stack > 8.0:
                npoly = 3
            elif rms_sn_stack >= 5.0:
                npoly = 2
            else:
                npoly = 1
        scale, flux_scale, ivar_scale, outmask = solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, npoly,
                                                                      mask=mask, mask_ref=mask_ref, debug=debug)
    elif scale_method == 'median':
        # Median ratio (reference to spectrum)
        med_scale = robust_median_ratio(flux,ivar,flux_ref,ivar_ref,ref_percentile=ref_percentile,min_good=min_good,\
                                        mask=mask, mask_ref=mask_ref, maxiters=maxiters,\
                                        max_factor=max_median_factor,sigrej=sigrej)
        # Apply
        flux_scale = flux * med_scale
        ivar_scale = ivar * 1.0/med_scale**2
        scale = np.full_like(flux,med_scale)
    elif scale_method == 'hand':
        # Input?
        if hand_scale is None:
            msgs.error("Need to provide hand_scale parameter, single value")
        flux_scale = flux * hand_scale
        ivar_scale = ivar * 1.0 / hand_scale ** 2
        scale = np.full(flux.size, hand_scale)
    elif scale_method == 'none':
        flux_scale = flux.copy()
        ivar_scale = ivar.copy()
        scale = np.ones_like(flux)
    else:
        msgs.error("Scale method not recognized! Check documentation for available options")
    # Finish
    if debug:
        scale_spec_qa(wave, flux, ivar, flux_ref, ivar_ref, scale, scale_method, mask = mask, mask_ref=mask_ref,
                      title='Scaling Applied to the Data')

    return flux_scale, ivar_scale, scale, scale_method


def compute_stack(wave_grid, waves, fluxes, ivars, masks, weights):
    '''
    Compute the stacked spectrum based on spectra and wave_grid with weights being taken into account.
    Args:
        waves:
        fluxes:
        ivars:
        masks:
        wave_grid:
        weights:
    Returns:
        weighted stacked wavelength, flux and ivar
    '''

    ubermask = masks & (weights > 0.0) & (waves > 1.0) & (ivars > 0.0)
    waves_flat = waves[ubermask].flatten()
    fluxes_flat = fluxes[ubermask].flatten()
    ivars_flat = ivars[ubermask].flatten()
    vars_flat = utils.calc_ivar(ivars_flat)
    weights_flat = weights[ubermask].flatten()

    # Counts how many pixels in each wavelength bin
    nused, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False)

    # Calculate the summed weights for the denominator
    weights_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=weights_flat)

    # Calculate the stacked wavelength
    wave_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=waves_flat*weights_flat)
    wave_stack = (weights_total > 0.0)*wave_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked flux
    flux_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=fluxes_flat*weights_flat)
    flux_stack = (weights_total > 0.0)*flux_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked ivar
    var_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=vars_flat*weights_flat**2)
    var_stack = (weights_total > 0.0)*var_stack_total/(weights_total+(weights_total==0.))**2
    ivar_stack = utils.calc_ivar(var_stack)

    # New mask for the stack
    mask_stack = (weights_total > 0.0) & (nused > 0.0)

    return wave_stack, flux_stack, ivar_stack, mask_stack, nused

def coadd_iexp_qa(wave, flux, ivar, flux_stack, ivar_stack, mask=None, mask_stack=None,
                  norder=None, qafile=None):

    if mask is None:
        mask = ivar > 0.0
    if mask_stack is None:
        mask_stack = ivar_stack > 0.0

    wave_mask = wave > 1.0
    fig = plt.figure(figsize=(12, 8))
    spec_plot = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Get limits
    npix = flux.size
    med_width = (2.0 * np.ceil(0.03 / 2.0 * npix) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    ymax = 2.0 * flux_med.max()
    ymin = -0.15 * ymax

    # Plot spectrum
    rejmask = (mask == False) & (wave_mask == True)
    spec_plot.plot(wave[rejmask], flux[rejmask],'s',zorder=10,mfc='None', mec='r', label='Rejected pixels')

    if norder is None:
        spec_plot.plot(wave[wave_mask], flux[wave_mask], color='dodgerblue', linestyle='steps-mid',
                       zorder=2, alpha=0.7,label='Single exposure')
        spec_plot.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(ivar[wave_mask])),zorder=3,
                       color='0.7', linestyle='steps-mid')
        spec_plot.plot(wave[wave_mask],flux_stack[wave_mask]*mask_stack[wave_mask],color='k',
                       linestyle='steps-mid',lw=2,zorder=1,label='Coadd')

        # Plot transmission
        if (np.max(wave[mask]) > 9000.0):
            skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
            skycat = np.genfromtxt(skytrans_file, dtype='float')
            scale = 0.8 * ymax
            spec_plot.plot(skycat[:, 0] * 1e4, skycat[:, 1] * scale, 'm-', alpha=0.5, zorder=11)
    else:
        nspec = int(npix / norder)
        for iord in range(norder):
            spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           flux[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           linestyle='steps-mid', zorder=1, alpha=0.7)
            #spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
            #               np.sqrt(utils.calc_ivar(ivar[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]])),
            #               zorder=3, color='0.7', linestyle='steps-mid')
            spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           flux_stack[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]]*
                           mask_stack[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           color='k', linestyle='steps-mid',lw=1,zorder=2)

    # properties
    spec_plot.legend(fontsize=13)
    spec_plot.set_ylim([ymin, ymax])
    spec_plot.set_xlim([wave[wave_mask].min(), wave[wave_mask].max()])
    spec_plot.set_xlabel('Wavelength (Angstrom)')
    spec_plot.set_ylabel('Flux')

    if qafile is not None:
        if len(qafile.split('.'))==1:
            msgs.info("No fomat given for the qafile, save to PDF format.")
            qafile = qafile+'.pdf'
        plt.savefig(qafile,dpi=300)
        msgs.info("Wrote QA: {:s}".format(qafile))
    plt.show()

    return

def weights_qa(waves, weights, masks):

    nexp = np.shape(waves)[1]
    fig = plt.figure(figsize=(12, 8))

    wave_mask_all = np.zeros_like(masks)
    for iexp in range(nexp):
        this_wave, this_weights, this_mask = waves[:, iexp], weights[:, iexp], masks[:, iexp]
        wave_mask = this_wave > 1.0
        wave_mask_all[wave_mask, iexp] = True
        plt.plot(this_wave[wave_mask], this_weights[wave_mask] * this_mask[wave_mask])
    plt.xlim([waves[wave_mask_all].min(), waves[wave_mask_all].max()])
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('(S/N)^2 weights')
    plt.show()

def coadd_qa(wave, flux, ivar, nused, mask=None, title=None, qafile=None):

    from matplotlib.ticker import MaxNLocator

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
    num_plot.plot(wave[wave_mask],nused[wave_mask],linestyle='steps-mid',color='k',lw=2)
    num_plot.set_xlim([wave_min, wave_max])
    num_plot.set_ylim([0.0, np.fmax(1.1*nused.max(), nused.max()+1.0)])
    num_plot.set_ylabel('$\\rm N_{EXP}$')
    num_plot.yaxis.set_major_locator(MaxNLocator(integer=True))
    num_plot.yaxis.set_minor_locator(NullLocator())

    # Plot spectrum
    spec_plot.plot(wave[wave_mask], flux[wave_mask], color='blue', linestyle='steps-mid',zorder=1,label='Single exposure')
    spec_plot.plot(wave[wave_mask], np.sqrt(utils.calc_ivar(ivar[wave_mask])),zorder=2, color='0.7', linestyle='steps-mid')

    # Get limits
    nspec = flux.size
    med_width = (2.0 * np.ceil(0.03 / 2.0 * nspec) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    ymax = 2.0 * flux_med.max()
    ymin = -0.15 * ymax

    # Plot transmission
    if (np.max(wave[mask])>9000.0):
        skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ymax
        spec_plot.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5,zorder=11)

    spec_plot.set_ylim([ymin, ymax])
    spec_plot.set_xlim([wave_min, wave_max])
    spec_plot.set_xlabel('Wavelength (Angstrom)')
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

    return

def update_errors(waves, fluxes, ivars, masks, fluxes_stack, ivars_stack, masks_stack, sn_cap=20.0, debug=False):

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
        #if nexp>1: # JXP TOUCHED THIS
        if fluxes.ndim>1:
            # Grab the spectrum
            thisflux = fluxes[:, iexp]
            thisivar = ivars[:, iexp]
            thismask = outmasks[:,iexp]
            # Grab the stack interpolated with the same grid as the current exposure
            thisflux_stack = fluxes_stack[:, iexp]
            thisvar_stack = utils.calc_ivar(ivars_stack[:, iexp])
            thismask_stack = masks_stack[:, iexp]
        else:
            thisflux = fluxes
            thisivar = ivars
            thismask = outmasks
            # Grab the stack interpolated with the same grid as the current exposure
            thisflux_stack = fluxes_stack
            thisvar_stack = utils.calc_ivar(ivars_stack)
            thismask_stack = masks_stack

        # var_tot
        var_tot = thisvar_stack + utils.calc_ivar(thisivar)
        ivar_tot = utils.calc_ivar(var_tot)
        mask_tot = thismask & thismask_stack
        # TODO Do we need the offset code? If so add it right here into the chi
        chi = np.sqrt(ivar_tot)*(thisflux - thisflux_stack)
        # Adjust errors to reflect the statistics of the distribution of errors. This fixes cases where the
        # the noise model is not quite right
        this_sigma_corr, igood = renormalize_errors(chi, mask_tot, clip=6.0, max_corr=5.0, title='spec_reject', debug=debug)
        ivar_tot_corr = ivar_tot/this_sigma_corr ** 2
        ivar_cap = np.minimum(ivar_tot_corr, (sn_cap/(thisflux_stack + (thisflux_stack <= 0.0))) ** 2)
        # if nexp>1:  #JXP TOUCHED THIS
        if fluxes.ndim>1:
            sigma_corrs[iexp] = this_sigma_corr
            rejivars[:, iexp] = ivar_cap
            outchi[:, iexp] = chi
            maskchi[:, iexp] = igood
        else:
            sigma_corrs = np.array([this_sigma_corr])
            rejivars = ivar_cap
            outchi = chi
            maskchi = igood

    return rejivars, sigma_corrs, outchi, maskchi



def spec_reject_comb(wave_grid, waves, fluxes, ivars, masks, weights, sn_cap=20.0, lower=3.0, upper=3.0,
                     maxrej=None, maxiter_reject=5, title=None, qafile=None, debug=False, show=False):

    nexp = np.shape(waves)[1]
    iIter = 0
    qdone = False
    thismask = np.copy(masks)
    while (not qdone) and (iIter < maxiter_reject):
        wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(
            wave_grid, waves, fluxes, ivars, thismask, weights)
        flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(
            waves, wave_stack, flux_stack, ivar_stack,mask_stack)
        rejivars, sigma_corrs, outchi, maskchi = update_errors(waves, fluxes, ivars, thismask,
                                                               flux_stack_nat, ivar_stack_nat, mask_stack_nat,
                                                               sn_cap=sn_cap)
        thismask, qdone = pydl.djs_reject(fluxes, flux_stack_nat, outmask=thismask,inmask=masks, invvar=rejivars,
                                          lower=lower,upper=upper, maxrej=maxrej, sticky=False)
        # print out how much was rejected
        for iexp in range(nexp):
            thisreject = thismask[:, iexp]
            nrej = np.sum(np.invert(thisreject))
            if nrej > 0:
                msgs.info("Rejecting {:d} pixels in exposure {:d}".format(nrej, iexp))

        iIter = iIter +1

    if (iIter == maxiter_reject) & (maxiter_reject != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter_reject) + ' reached in combspec')
    outmask = np.copy(thismask)

    # Compute the final stack using this outmask
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(
        wave_grid, waves, fluxes, ivars, outmask, weights)
    # Used only for plotting below
    flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(waves, wave_stack, flux_stack, ivar_stack, mask_stack)
    if debug:
        for iexp in range(nexp):
            # plot the residual distribution
            renormalize_errors_qa(outchi[:, iexp], maskchi[:, iexp], sigma_corrs[iexp])
            # plot the rejections for each exposures
            coadd_iexp_qa(waves[:, iexp], fluxes[:, iexp], ivars[:, iexp], flux_stack_nat[:, iexp],
                          ivar_stack_nat[:, iexp], mask=outmask[:, iexp], mask_stack=mask_stack_nat[:, iexp],
                          qafile=None)
        # weights qa
        weights_qa(waves, weights, outmask)

    # Plot the final coadded spectrum
    if show or debug:
        coadd_qa(wave_stack,flux_stack,ivar_stack, nused, mask=mask_stack, title= title, qafile=qafile)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused

def long_combspec(wave_grid, waves, fluxes, ivars, masks, ref_percentile=30.0, maxiter_scale=5, sigrej=3,
                  scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, dv_smooth=10000.0,
                  const_weights=False, maxiter_reject=5, sn_cap=20.0, lower=3.0, upper=3.0, maxrej=None,
                  header=None, ex_value='OPT', title=None, qafile=None, outfile=None, debug=False, show=False):

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves,fluxes,ivars,masks, dv_smooth=dv_smooth, const_weights=const_weights, verbose=True)

    # Compute an initial stack as the reference, this has its own wave grid based on the weighted averages
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, masks, weights)
    # Interpolate the stack onto each individual exposures native wavelength grid
    flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(waves, wave_stack, flux_stack, ivar_stack, mask_stack)

    # Rescale spectra to line up with our preliminary stack so that we can sensibly reject outliers
    nexp = np.shape(fluxes)[1]
    fluxes_scale = np.zeros_like(fluxes)
    ivars_scale = np.zeros_like(ivars)
    scales = np.zeros_like(fluxes)
    for iexp in range(nexp):
        # TODO Create a parset for the coadd parameters!!!
        fluxes_scale[:, iexp], ivars_scale[:, iexp], scales[:, iexp], omethod = scale_spec(
            waves[:, iexp],fluxes[:, iexp],ivars[:, iexp], flux_stack_nat[:, iexp], ivar_stack_nat[:, iexp],
            mask=masks[:, iexp], mask_ref=mask_stack_nat[:, iexp], ref_percentile=ref_percentile, maxiters=maxiter_scale,
            sigrej=sigrej, scale_method=scale_method, hand_scale=hand_scale, sn_max_medscale=sn_max_medscale,
            sn_min_medscale=sn_min_medscale, debug=debug)

    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves, fluxes_scale, ivars_scale, masks, weights, sn_cap=sn_cap, lower=lower, upper=upper,
        maxrej=maxrej, maxiter_reject=maxiter_reject, title=title, qafile=qafile, debug=debug, show=show)

    # Write to disk?
    if outfile is not None:
        if len(outfile.split('.'))==1:
            msgs.info("No fomat given for the outfile, save to fits format.")
            outfile = outfile+'.fits'
        #write_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, outfile, clobber=True)
        save.save_coadd1d_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
                                  outfile=outfile, ex_value=ex_value, overwrite=True)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales, rms_sn

#Todo: This should probaby take a parset?
#Todo: Make this works for multiple objects after the coadd script input file format is fixed.
def multi_combspec(fnames, objids, ex_value='OPT', flux_value=True, wave_method='pixel', A_pix=None, v_pix=None,
                   samp_fact = 1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
                   sigrej=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
                   dv_smooth=1000.0, const_weights=False, maxiter_reject=5, sn_cap=20.0, lower=3.0, upper=3.0,
                   maxrej=None, max_factor=10.0, maxiters=5, min_good=0.05, phot_scale_dicts=None, nmaskedge=2,
                   qafile=None, outfile = None, debug=False, show=False):

    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=objids, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value, nmaskedge=nmaskedge)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Coadd
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, weights, scales, rms_sn = \
        long_combspec(wave_grid, waves, fluxes, ivars, masks, ref_percentile=ref_percentile,
                      maxiter_scale=maxiter_scale, sigrej=sigrej, scale_method=scale_method, hand_scale=hand_scale,
                      sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
                      const_weights=const_weights, maxiter_reject=maxiter_reject, sn_cap=sn_cap, lower=lower,
                      upper=upper, maxrej=maxrej, title='Stacked spectrum', qafile=qafile,
                      outfile=outfile, debug=debug, show=show)

    return wave_stack, flux_stack, ivar_stack, mask_stack

def ech_combspec(fnames, objids, ex_value='OPT', flux_value=True, wave_method='loggrid', A_pix=None, v_pix=None,
                 samp_fact=1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
                 sigrej=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
                 dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_cap=20.0, lower=3.0, upper=3.0,
                 maxrej=None, max_factor=10.0, maxiters=5, min_good=0.05, phot_scale_dicts=None, nmaskedge=2,
                 qafile=None, outfile = None, debug=False, show=False):

    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=objids, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value, nmaskedge=nmaskedge)
    # data shape
    data_shape = np.shape(waves)
    npix = data_shape[0] # detector size in the wavelength direction
    norder = data_shape[1]
    nexp = data_shape[2]

    # create some arrays
    scales = np.zeros_like(waves)
    weights = np.zeros_like(waves)
    outmasks = np.zeros_like(waves,dtype=bool)

    # output name root for fits and QA plots
    if outfile is None:
        outfile = header['TARGET']
    outfile_order = 'spec1d_order_{:}'.format(outfile)
    outfile_stack = 'spec1d_stack_{:}'.format(outfile)
    outfile_giant_stack = 'spec1d_giant_stack_{:}'.format(outfile)

    if qafile is None:
        qafile = header['TARGET']
    qafile_stack = 'spec1d_stack_{:}'.format(outfile)
    qafile_giant_stack = 'spec1d_giant_stack_{:}'.format(outfile)
    qafile_chi = 'spec1d_chi_{:}'.format(outfile)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Arrays to store stacked individual order spectra.
    waves_stack_orders = np.zeros((np.size(wave_grid)-1, norder))
    fluxes_stack_orders = np.zeros_like(waves_stack_orders)
    ivars_stack_orders = np.zeros_like(waves_stack_orders)
    masks_stack_orders = np.zeros_like(waves_stack_orders,dtype=bool)

    # Loop over orders to get the initial stacks of each individual order
    for ii in range(norder):

        # get the slice of iord spectra of all exposures
        waves_iord, fluxes_iord, ivars_iord, masks_iord = waves[:,ii,:], fluxes[:,ii,:], ivars[:,ii,:], masks[:,ii,:]

        # Get the stacked spectrum for each order
        waves_stack_orders[:, ii], fluxes_stack_orders[:, ii], ivars_stack_orders[:, ii],  masks_stack_orders[:, ii], \
        outmask_iord, weights_iord, scales_iord, rms_sn_iord = \
                    long_combspec(wave_grid, waves_iord, fluxes_iord, ivars_iord, masks_iord, ref_percentile=ref_percentile,
                                  maxiter_scale=maxiter_scale, sigrej=sigrej, scale_method=scale_method, hand_scale=hand_scale,
                                  sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
                                  const_weights=const_weights, maxiter_reject=maxiter_reject, sn_cap=sn_cap, lower=lower,
                                  upper=upper, maxrej=maxrej, title='Stacked spectrum of order {:}'.format(ii+1), qafile=None,
                                  outfile=None, debug=debug, show=show)

        # store new masks, scales and weights, all of these arrays are in native wave grid
        scales[:,ii,:] = scales_iord
        weights[:,ii,:] = weights_iord
        outmasks[:,ii,:] = outmask_iord

    # Now that we have high S/N ratio individual order stacks, let's compute re-scaling fractors from the order
    # overlaps. We will work from red to blue.
    fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = \
        order_median_scale(waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
                           min_good=min_good, maxiters=maxiters, max_factor=max_factor, sigrej=sigrej,
                           debug=debug, show=show)

    ## Stack with the first method: combine the stacked individual order spectra directly
    # Get weights for individual order stacks
    # ToDo: The sensfunc need to be put in here for the order merge, i.e. the weights should be sensfunc**2
    rms_sn_stack, weights_stack = sn_weights(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale,
                                             masks_stack_orders, dv_smooth=dv_smooth, const_weights=const_weights,
                                             sens_weights=True, verbose=True)
    '''
    from IPython import embed
    embed()
    for i in range(norder):
        wavei = waves_stack_orders[:,i]
        fluxi = fluxes_stack_orders_scale[:,i]
        sigi = np.sqrt(utils.calc_ivar(ivars_stack_orders_scale[:,i]))
        maski = masks_stack_orders[:,i]
        weightsi = weights_stack[:,i]
        plt.plot(wavei[maski], fluxi[maski]/sigi[maski],zorder=1)
        plt.plot(wavei[maski], weightsi[maski],zorder=2)
    plt.show()
    '''

    # TODO Will we use this reject/stack below? It is the straight combine of the stacked individual orders.
    #  This does not take advantage
    #  of the fact that we have many samples in the order overlap regions allowing us to better reject. It does
    #  however have the advatnage that it operates on higher S/N ratio stacked spectra.
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
        weights_stack, sn_cap=sn_cap, lower=lower, upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject,
        title='Order merged stacked spectra', qafile=qafile_stack, debug=debug, show=show)

    # Save stacked individual order spectra
    save.save_coadd1d_to_fits(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
                              header=header, outfile=outfile_order, ex_value = ex_value, overwrite=True)
    save.save_coadd1d_to_fits(wave_stack, flux_stack, ivar_stack, mask_stack, header=header, outfile=outfile_stack,
                              ex_value = ex_value, overwrite=True)

    # apply order_ratios to the scales array: order_ratio*scale
    scales_new = np.copy(scales)
    for ii in range(norder):
        scales_new[:,ii,:] *= order_ratios[ii]

    fluxes_scale = fluxes * scales_new
    ivars_scale = ivars/scales_new**2

    # apply the sensfunc weights to the orginal weights: sensfunc_weights*weights
    # TODO: need to make sure whether we can times these two weights directly, seems working correctly
    weights_new = np.copy(weights)
    for iord in range(norder):
        weight_iord = weights_stack[:, iord]
        wave_weight_iord = waves_stack_orders[:, iord]
        mask_weight_iord = masks_stack_orders[:, iord] & (weight_iord>0) & (wave_weight_iord>1.0)
        for iexp in range(nexp):
            weight_iord_interp = np.interp(waves[:,iord,iexp],wave_weight_iord[mask_weight_iord],
                                           weight_iord[mask_weight_iord])
            weights_new[:,iord,iexp] = weights[:,iord,iexp] * weight_iord_interp

    # reshaping 3D arrays (npix, norder, nexp) to 2D arrays (npix, norder*nexp)
    # need Fortran like order reshaping to make sure you are getting the right spectrum for each expsoure
    waves_2d = np.reshape(waves,(npix, norder*nexp),order='F')
    fluxes_2d = np.reshape(fluxes_scale, np.shape(waves_2d),order='F')
    ivars_2d = np.reshape(ivars_scale, np.shape(waves_2d),order='F')
    masks_2d = np.reshape(outmasks, np.shape(waves_2d),order='F')
    weights_2d = np.reshape(weights_new, np.shape(waves_2d),order='F')

    wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmask_giant_stack, nused_giant_stack = \
        spec_reject_comb(wave_grid, waves_2d, fluxes_2d, ivars_2d, masks_2d, weights_2d, sn_cap=sn_cap, lower=lower,
                         upper=upper, maxrej=maxrej, maxiter_reject=maxiter_reject, title='Final stacked spectra',
                         qafile=qafile_giant_stack, debug=debug, show=show)

    if debug or show:
        # Estimate chi for each exposure after interpolating and reshaping
        flux_stack_2d, ivar_stack_2d, mask_stack_2d = interp_spec(waves_2d, wave_giant_stack, flux_giant_stack,
                                                                  ivar_giant_stack, mask_giant_stack)
        waves_2d_exps = waves_2d.reshape((npix*norder, nexp),order='F')
        fluxes_2d_exps = fluxes_2d.reshape(np.shape(waves_2d_exps),order='F')
        ivars_2d_exps = ivars_2d.reshape(np.shape(waves_2d_exps),order='F')
        masks_2d_exps = masks_2d.reshape(np.shape(waves_2d_exps),order='F')
        flux_stack_2d_exps = flux_stack_2d.reshape(np.shape(waves_2d_exps),order='F')
        ivar_stack_2d_exps = ivar_stack_2d.reshape(np.shape(waves_2d_exps),order='F')
        mask_stack_2d_exps = mask_stack_2d.reshape(np.shape(waves_2d_exps),order='F')

        rejivars_2d, sigma_corrs_2d, outchi_2d, maskchi_2d = update_errors(
            waves_2d_exps, fluxes_2d_exps, ivars_2d_exps, masks_2d_exps, flux_stack_2d_exps, ivar_stack_2d_exps,
            mask_stack_2d_exps, sn_cap=sn_cap)

        if debug:
            # QA for individual exposures
            for iexp in range(nexp):
                # plot the residual distribution
                msgs.info('QA plots for exposure {:} with new_sigma = {:}'.format(iexp+1, sigma_corrs_2d[iexp]))
                renormalize_errors_qa(outchi_2d[:, iexp], maskchi_2d[:, iexp], sigma_corrs_2d[iexp],
                                      title='Chi distribution of exposure {:}'.format(iexp+1))

                # plot coadd_qa
                coadd_iexp_qa(waves_2d_exps[:,iexp], fluxes_2d_exps[:,iexp], ivars_2d_exps[:,iexp],
                              flux_stack_2d_exps[:,iexp], ivar_stack_2d_exps[:,iexp], mask=masks_2d_exps[:,iexp],
                              mask_stack=mask_stack_2d_exps[:,iexp], norder=norder, qafile=None)
        # Global QA
        rejivars_1d, sigma_corrs_1d, outchi_1d, maskchi_1d = update_errors(
            waves_2d_exps.flatten(), fluxes_2d_exps.flatten(), ivars_2d_exps.flatten(), masks_2d_exps.flatten(),
            flux_stack_2d_exps.flatten(), ivar_stack_2d_exps.flatten(), mask_stack_2d_exps.flatten(), sn_cap=sn_cap)
        renormalize_errors_qa(outchi_1d, maskchi_1d, sigma_corrs_1d[0], qafile=qafile_chi, title='Global Chi distribution')

    save.save_coadd1d_to_fits(wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack,
                              header=header, outfile=outfile_giant_stack, ex_value=ex_value, overwrite=True)

    return wave_stack, flux_stack, ivar_stack, mask_stack