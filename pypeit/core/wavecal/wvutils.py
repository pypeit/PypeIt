""" Module for basic utilties with holy grail


.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""
import numpy as np
import os

from matplotlib import pyplot as plt

from scipy.ndimage.filters import gaussian_filter
from scipy.signal import resample
import scipy
from scipy.optimize import curve_fit

from astropy.table import Table
from astropy import convolution
from astropy import constants

from pypeit import msgs
from pypeit import utils
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

# TODO: Should this code allow the user to skip the smoothing steps and just
# provide the raw delta_wave vector? I would think there are cases where you
# want the *exact* pixel width, as opposed to the smoothed version.
def get_delta_wave(wave, wave_gpm, frac_spec_med_filter=0.03):
    r"""
    Compute the change in wavelength per pixel.

    Given an input wavelength vector and an input good pixel mask, the *raw*
    change in wavelength is defined to be ``delta_wave[i] =
    wave[i+1]-wave[i]``, with ``delta_wave[-1] = delta_wave[-2]`` to force
    ``wave`` and ``delta_wave`` to have the same length.

    The method imposes a smoothness on the change in wavelength by (1)
    running a median filter over the raw values and (2) smoothing
    ``delta_wave`` with a Gaussian kernel. The boxsize for the median filter
    is set by ``frac_spec_med_filter``, and the :math:`\sigma` for the
    Gaussian kernel is either a 10th of that boxsize or 3 pixels (whichever
    is larger).

    Parameters
    ---------- 
    wave : float `numpy.ndarray`_, shape = (nspec,)
        Array of input wavelengths. Must be 1D.
    wave_gpm : bool `numpy.ndarray`_, shape = (nspec)
        Boolean good-pixel mask defining where the ``wave`` values are
        good.
    frac_spec_med_filter : :obj:`float`, optional
        Fraction of the length of the wavelength vector to use to median
        filter the raw change in wavelength, used to impose a smoothness.
        Default is 0.03, which means the boxsize for the running median
        filter will be approximately ``0.03*nspec`` (forced to be an odd
        number of pixels).

    Returns
    -------
    delta_wave : `numpy.ndarray`_, float, shape = (nspec,)
        A smooth estimate for the change in wavelength for each pixel in the
        input wavelength vector.
    """
    # Check input
    if wave.ndim != 1:
        msgs.error('Input wavelength array must be 1D.')

    nspec = wave.size
    # This needs to be an odd number
    nspec_med_filter = 2*int(np.round(nspec*frac_spec_med_filter/2.0)) + 1
    delta_wave = np.zeros_like(wave)
    wave_diff = np.diff(wave[wave_gpm])
    wave_diff = np.append(wave_diff, wave_diff[-1])
    wave_diff_filt = utils.fast_running_median(wave_diff, nspec_med_filter)

    # Smooth with a Gaussian kernel
    sig_res = np.fmax(nspec_med_filter/10.0, 3.0)
    gauss_kernel = convolution.Gaussian1DKernel(sig_res)
    wave_diff_smooth = convolution.convolve(wave_diff_filt, gauss_kernel, boundary='extend')
    delta_wave[wave_gpm] = wave_diff_smooth
    return delta_wave


def get_sampling(waves, pix_per_R=3.0):
    """
    Computes the median wavelength sampling of wavelength vector(s)

    Args:
        waves (float `numpy.ndarray`_): shape = (nspec,) or (nspec, nimgs)
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


# TODO: the other methods iref should be deprecated or removed
def get_wave_grid(waves, masks=None, wave_method='linear', iref=0, wave_grid_min=None,
                  wave_grid_max=None, dwave=None, dv=None, dloglam=None, spec_samp_fact=1.0):
    """
    Create a new wavelength grid for spectra to be rebinned and coadded.

    Args:
        waves (`numpy.ndarray`_):
            Set of N original wavelength arrays.  Shape is (nspec, nexp).
        masks (`numpy.ndarray`_, optional):
            Good-pixel mask for wavelengths.  Shape must match waves.
        wave_method (:obj:`str`, optional):
            Desired method for creating new wavelength grid:

                * 'iref' -- Use the first wavelength array (default)
                * 'velocity' -- Grid is uniform in velocity
                * 'log10'  -- Grid is uniform in log10(wave). This is the same as velocity.
                * 'linear' -- Constant pixel grid
                * 'concatenate' -- Meld the input wavelength arrays

        iref (:obj:`int`, optional):
            Index in waves array for reference spectrum
        wave_grid_min (:obj:`float`, optional):
            min wavelength value for the final grid
        wave_grid_max (:obj:`float`, optional):
            max wavelength value for the final grid
        dwave (:obj:`float`, optional):
            Pixel size in same units as input wavelength array (e.g. Angstroms).
            If not input, the median pixel size is calculated and used.
        dv (:obj:`float`, optional):
            Pixel size in km/s for velocity method.  If not input, the median
            km/s per pixel is calculated and used
        dloglam (:obj:`float`, optional):
            Pixel size in log10(wave) for the log10 method.
        spec_samp_fact (:obj:`float`, optional):
            Make the wavelength grid sampling finer (spec_samp_fact < 1.0) or
            coarser (spec_samp_fact > 1.0) by this sampling factor. This
            basically multiples the 'native' spectral pixels by
            ``spec_samp_fact``, i.e. units ``spec_samp_fact`` are pixels.

    Returns:
        :obj:`tuple`: Returns two `numpy.ndarray`_ objects and a float:
            - ``wave_grid``: New wavelength grid, not masked
            - ``wave_grid_mid``: New wavelength grid evaluated at the centers of
              the wavelength bins, that is this grid is simply offset from
              ``wave_grid`` by ``dsamp/2.0``, in either linear space or log10
              depending on whether linear or (log10 or velocity) was requested.
              For iref or concatenate the linear wavelength sampling will be
              calculated.
            - ``dsamp``: The pixel sampling for wavelength grid created.
    """
    c_kms = constants.c.to('km/s').value

    if masks is None:
        masks = waves > 1.0

    if wave_grid_min is None:
        wave_grid_min = waves[masks].min()
    if wave_grid_max is None:
        wave_grid_max = waves[masks].max()

    dwave_data, dloglam_data, resln_guess, pix_per_sigma = get_sampling(waves)

    # TODO: These tests of the string value should not use 'in', they should use ==
    if ('velocity' in wave_method) or ('log10' in wave_method):
        if dv is not None and dloglam is not None:
            msgs.error('You can only specify dv or dloglam but not both')
        elif dv is not None:
            dloglam_pix = dv/c_kms/np.log(10.0)
        elif dloglam is not None:
            dloglam_pix = dloglam
        else:
            dloglam_pix = dloglam_data

        # Generate wavelength array
        wave_grid = wavegrid(wave_grid_min, wave_grid_max, dloglam_pix,
                             spec_samp_fact=spec_samp_fact, log10=True)
        loglam_grid_mid = np.log10(wave_grid) + dloglam_pix*spec_samp_fact/2.0
        wave_grid_mid = np.power(10.0, loglam_grid_mid)
        dsamp = dloglam_pix

    elif 'linear' in wave_method: # Constant Angstrom
        if dwave is not None:
            dwave_pix = dwave
        else:
            dwave_pix = dwave_data
        # Generate wavelength array
        wave_grid = wavegrid(wave_grid_min, wave_grid_max, dwave_pix, spec_samp_fact=spec_samp_fact)
        wave_grid_mid = wave_grid + dwave_pix*spec_samp_fact/2.0
        dsamp = dwave_pix

    elif 'concatenate' in wave_method:  # Concatenate
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
        wave_grid = np.power(10.0,newloglam)

    elif 'iref' in wave_method:
        wave_tmp = waves[:, iref]
        wave_grid = wave_tmp[ wave_tmp > 1.0]

    else:
        msgs.error("Bad method for wavelength grid: {:s}".format(wave_method))

    if ('iref' in wave_method) | ('concatenate' in wave_method):
        wave_grid_diff = np.diff(wave_grid)
        wave_grid_diff = np.append(wave_grid_diff, wave_grid_diff[-1])
        wave_grid_mid = wave_grid + wave_grid_diff / 2.0
        dsamp = np.median(wave_grid_diff)

    return wave_grid, wave_grid_mid, dsamp


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
    theta (float `numpy.ndarray`_):
        Function parameters to optmize over. theta[0] = shift, theta[1] = stretch
    y1 (float `numpy.ndarray`_):  shape = (nspec,)
        First spectrum which acts as the refrence
    y2 (float `numpy.ndarray`_):  shape = (nspec,)
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
        plt.plot(lags, corr_norm, color='black', drawstyle = 'steps-mid', lw=3, label = 'x-corr')
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



def wavegrid(wave_min, wave_max, dwave, spec_samp_fact=1.0, log10=False):
    """

    Utility routine to generate a uniform grid of wavelengths

    Args:
        wave_min (float):
           Mininum wavelength. Must be linear even if log10 is requested
        wave_max (float):
           Maximum wavelength. Must be linear even if log10 is requested.
        dwave (float):
           Delta wavelength interval. Must be linear if log10=False, or log10 if log10=True
        spec_samp_fact (float, optional):
            Make the wavelength grid  sampling finer (spec_samp_fact < 1.0) or coarser (spec_samp_fact > 1.0) by this
            sampling factor. This basically multiples the 'native' spectral pixels by spec_samp_fact, i.e. units
            spec_samp_fact are pixels.

    Returns:
        `numpy.ndarray`_: Wavelength grid in Angstroms (i.e. log10 even
        if log10 is requested)

    """

    dwave_eff = dwave*spec_samp_fact
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
