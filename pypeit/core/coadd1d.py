import os
import scipy
import numpy as np
import matplotlib.pyplot as plt
from astropy import stats
from astropy.io import fits
from astropy import convolution
from IPython import embed
from astropy.table import Table

from pkg_resources import resource_filename
from pypeit import utils
from pypeit import msgs
from pypeit.core import load, save
from pypeit.core.wavecal import wvutils
from pypeit.core import pydl
from astropy import constants


from matplotlib.ticker import NullFormatter, NullLocator, MaxNLocator

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

# Used it in several places, so make it global.
c_kms = constants.c.to('km/s').value


# TODO: merge with wavegrid routine in wvutils


def ech_combspec_old(fnames, objids, sensfile=None, ex_value='OPT', flux_value=True, wave_method='loggrid', A_pix=None, v_pix=None,
                 samp_fact=1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
                 sigrej_scale=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
                 dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0,
                 maxrej=None, max_factor=10.0, maxiters=5, min_good=0.05, phot_scale_dicts=None, nmaskedge=2,
                 qafile=None, outfile = None, order_scale=False,
                 merge_stack=False, debug_scale=False, debug_order_scale=False, debug=False, show=False):
    '''
    Driver routine for coadding Echelle spectra. Calls combspec which is the main stacking algorithm. It will deliver
    three fits files: spec1d_order_XX.fits (stacked individual orders, one order per extension), spec1d_merge_XX.fits
    (straight combine of stacked individual orders), spec1d_stack_XX.fits (a giant stack of all exposures and all orders).
    In most cases, you should use spec1d_stack_XX.fits for your scientific analyses since it reject most outliers.

    Args:
        fnames: list
           a list of spec1d fits file names
        objids: list
           objids (e.g. 'OBJ0001') you want to combine of that spectrum in the spec1d fits files
        sensfile: str, default = None for a smoothed ivar weighting when sticking different orders
        ex_value: str, default = 'OPT' for optimal extraction, 'BOX' for boxcar extraction.
        flux_value: bool, default=True
           if True coadd fluxed spectrum, if False coadd spectra in counts
        wave_method: str, default=pixel
           method for generating new wavelength grid with new_wave_grid. Deafult is 'pixel' which creates a uniformly
           space grid in lambda
        A_pix: float,
           dispersion in units of A in case you want to specify it for new_wave_grid, otherwise the code computes the
           median spacing from the data.
        v_pix: float,
           Dispersion in units of km/s in case you want to specify it in the new_wave_grid  (for the 'velocity' option),
           otherwise a median value is computed from the data.
        samp_fact: float, default=1.0
           sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
           samp_fact < 1.0 undersamples (coarser).
        wave_grid_min: float, default=None
           In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
        wave_grid_max: float, default=None
           In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
        ref_percentile:
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
        maxiter_scale: int, default=5
            Maximum number of iterations performed for rescaling spectra.
        max_median_factor: float, default=10.0
            maximum scale factor for median rescaling for robust_median_ratio if median rescaling is the method used.
        sigrej_scale: flaot, default=3.0
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
        scale_method: scale method, str, default=None.
            Options are poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default is to let the
            code determine this automitically which works well.
        hand_scale: ndarray,
            Array of hand scale factors, not well tested
        sn_max_medscale: float, default = 2.0,
            maximum SNR for perforing median scaling
        sn_min_medscale: float, default = 0.5
            minimum SNR for perforing median scaling
        dv_smooth: float, 10000.0
            Velocity smoothing used for determining smoothly varying S/N ratio weights by sn_weights
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        const_weights: ndarray, (nexp,)
             Constant weight factors specif
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        max_factor: float, default = 10.0,
            Maximum allowed value of the returned ratio
        maxiters: int, defrault = 5,
            Maximum number of iterations for astropy.stats.SigmaClip
        min_good: float, default = 0.05
            Minimum fraction of good pixels determined as a fraction of the total pixels for estimating the median ratio
        phot_scale_dicts: dict,
            Dictionary for rescaling spectra to match photometry. Not yet implemented.
        nmaskedge: int, default=2
            Number of edge pixels to mask. This should be removed/fixed.
        qafile: str, default=None
            Root name for QA, if None, it will be determined either the outfile
        outfile: str, default=None,
            Root name for QA, if None, it will come from the target name from the fits header.
        order_scale: bool, default=False,
            Re-scale the orders to match up in the overlap regions. This is currently producing weird results for IR spectra
        merge_stack: bool, default=False,
            Compute an experimental combine of the high S/N combined orders in addition to the default algorithm,
            which is to compute one giant stack using all order overlaps

        debug: bool, default=False,
            Show all QA plots useful for debugging. Note there are lots of QA plots, so only set this to True if you want to inspect them all.
        debug_scale (bool): default=False
            Show interactive QA plots for the rescaling of the spectra for each individua order
        debug_order_scale (bool): default=False
            Show interactive QA plots for the rescaling of the spectra so that the overlap regions match from order to order
        show: bool, default=False,
             Show key QA plots or not

    Returns:
        wave_giant_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_giant_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_giant_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_giant_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_giant_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
    '''

    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=objids, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value, nmaskedge=nmaskedge)
    # data shape
    nspec, norder, nexp = waves.shape

    # create some arrays
    scales = np.zeros_like(waves)
    weights = np.zeros_like(waves)
    outmasks = np.zeros_like(waves,dtype=bool)

    # output name root for fits and QA plots
    if outfile is None:
        outfile = header['TARGET']+'.fits'
    elif len(outfile.split('.'))==1:
        outfile = outfile+'.fits'

    outfile_order = 'spec1d_order_{:}'.format(outfile)
    outfile_stack = 'spec1d_stack_{:}'.format(outfile)

    if qafile is None:
        qafile = outfile.split('.')[0]+'.pdf'
    qafile_stack = 'spec1d_stack_{:}'.format(qafile)
    qafile_chi = 'spec1d_chi_{:}'.format(qafile)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Arrays to store stacked individual order spectra.
    waves_stack_orders = np.zeros((np.size(wave_grid)-1, norder))
    fluxes_stack_orders = np.zeros_like(waves_stack_orders)
    ivars_stack_orders = np.zeros_like(waves_stack_orders)
    masks_stack_orders = np.zeros_like(waves_stack_orders,dtype=bool)

    # Loop over orders to get the initial stacks of each individual order
    for iord in range(norder):

        # Get the stacked spectrum for each order
        waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], masks_stack_orders[:, iord], \
        outmasks[:,iord,:], nused_iord, weights[:,iord,:], scales[:,iord,:], rms_sn_iord = combspec_old(
            wave_grid, waves[:,iord,:], fluxes[:,iord,:], ivars[:,iord,:], masks[:,iord,:], ref_percentile=ref_percentile,
            maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
            sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
            const_weights=const_weights, maxiter_reject=maxiter_reject, sn_clip=sn_clip, lower=lower,
            upper=upper, maxrej=maxrej, debug_scale=debug_scale, title='Order-by-Order Combine', debug=debug)

        if show:
            # TODO can we make this bit below more modular for the telluric?
            if sensfile is not None:
                tell_iord = get_tell_from_file(sensfile, waves_stack_orders[:, iord], masks_stack_orders[:, iord], iord=iord)
            else:
                tell_iord = None
            coadd_qa(waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], nused_iord,
                     mask=masks_stack_orders[:, iord], tell=tell_iord,
                     title='Coadded spectrum of order {:}'.format(iord+1))

    # TODO Is this order rescaling currently taking place??
    # Now that we have high S/N ratio individual order stacks, let's compute re-scaling fractors from the order
    # overlaps. We will work from red to blue.
    if order_scale:
        fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = order_median_scale(
            waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
            min_good=min_good, maxiters=maxiters, max_factor=max_factor, sigrej=sigrej_scale,
            debug=debug_order_scale, show=show)
    else:
        fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = fluxes_stack_orders, ivars_stack_orders, np.ones(norder)

    # apply order_ratios to the scales array: order_ratio*scale
    scales_new = np.zeros_like(scales)
    for iord in range(norder):
        scales_new[:,iord,:] = order_ratios[iord]*scales[:,iord,:]
    fluxes_scale = fluxes * scales_new
    ivars_scale = ivars/scales_new**2


    # Get the new ech_weights for the stack which will merge all the orders
    if sensfile is None:
        rms_sn_stack, order_weights = sn_weights(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale,
                                                 masks_stack_orders, dv_smooth=dv_smooth, const_weights=const_weights,
                                                 ivar_weights=True, verbose=True)
    else:
        rms_sn_stack = None
        order_weights, masks_stack_orders = sensfunc_weights(sensfile, waves_stack_orders, masks_stack_orders, debug=debug)

    #TODO think through whether this is the correct approach of multiplying weights?
    # apply the sensfunc weights to the orginal weights: sensfunc_weights*weightsf
    ech_weights = np.zeros_like(weights)
    for iord in range(norder):
        mask_weight_iord = masks_stack_orders[:, iord] & (order_weights[:, iord] > 0.0) & (waves_stack_orders[:, iord] > 1.0)
        # Interpolate these order_weights onto the native wavelength grid of each exposure for this order
        for iexp in range(nexp):
            order_weight_interp = scipy.interpolate.interp1d(
                waves_stack_orders[mask_weight_iord, iord], order_weights[mask_weight_iord, iord],  kind = 'cubic',
                bounds_error = False, fill_value = np.nan)(waves[:,iord,iexp])
            ech_weights[:,iord,iexp] = weights[:,iord,iexp] * order_weight_interp


    # TODO Will we use this reject/stack below? It is the straight combine of the stacked individual orders.
    #  This does not take advantage
    #  of the fact that we have many samples in the order overlap regions allowing us to better reject. It does
    #  however have the advatnage that it operates on higher S/N ratio stacked spectra.
    #  should we compute the stack directly with compute_stack or do more rejections with spec_reject_comb?
    #  spec_reject_comb will reject tons of pixels for overlap in telluric region.
    if merge_stack:
        ## Stack with the first method: combine the stacked individual order spectra directly
        wave_merge, flux_merge, ivar_merge, mask_merge, nused = compute_stack(
            wave_grid, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
            order_weights)
        if debug or show:
            qafile_merge = 'spec1d_merge_{:}'.format(qafile)
            coadd_qa(wave_merge, flux_merge, ivar_merge, nused, mask=mask_merge, tell = None,
                     title='Straight combined spectrum of the stacked individual orders', qafile=qafile_merge)


    #TODO Add a note here clarifyng how these reshaped spectra are arranged, i.e. are they packed by the order or by
    # by exposure.

    # reshaping 3D arrays (nspec, norder, nexp) to 2D arrays (nspec, norder*nexp)
    # need Fortran like order reshaping to make sure you are getting the right spectrum for each exposure
    waves_2d = np.reshape(waves,(nspec, norder*nexp),order='F')
    fluxes_2d = np.reshape(fluxes_scale, np.shape(waves_2d),order='F')
    ivars_2d = np.reshape(ivars_scale, np.shape(waves_2d),order='F')
    masks_2d = np.reshape(masks, np.shape(waves_2d),order='F')
    outmasks_2d = np.reshape(outmasks, np.shape(waves_2d),order='F')
    ech_weights_2d = np.reshape(ech_weights, np.shape(waves_2d),order='F')

    wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmask_giant_stack, nused_giant_stack = \
        spec_reject_comb(wave_grid, waves_2d, fluxes_2d, ivars_2d, outmasks_2d, ech_weights_2d, sn_clip=sn_clip,
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
        if debug:
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
                # plot coadd_qa
                coadd_iexp_qa(waves_2d_exps[:,iexp], fluxes_2d_exps[:,iexp], masks_2d_exps[:,iexp],
                              flux_stack_2d_exps[:,iexp], mask_stack_2d_exps[:,iexp],
                              rejivars_2d_exps[:,iexp], outmasks_2d_exps[:,iexp], norder=norder, qafile=None,
                              title=title_coadd_iexp)
        # Global QA
        rejivars_1d, sigma_corrs_1d, outchi_1d, maskchi_1d = update_errors(
            fluxes_2d_exps.flatten(), ivars_2d_exps.flatten(), outmasks_2d_exps.flatten(),
            flux_stack_2d_exps.flatten(), ivar_stack_2d_exps.flatten(), mask_stack_2d_exps.flatten(), sn_clip=sn_clip)
        renormalize_errors_qa(outchi_1d, maskchi_1d, sigma_corrs_1d[0], qafile=qafile_chi, title='Global Chi distribution')
        # show the final coadded spectrum
        coadd_qa(wave_giant_stack, flux_giant_stack, ivar_giant_stack, nused_giant_stack, mask=mask_giant_stack,
                 title='Final stacked spectrum', qafile=qafile_stack)

    # Save stacked individual order spectra
    save.save_coadd1d_to_fits(outfile_order, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
                              header=header, ex_value = ex_value, overwrite=True)
    save.save_coadd1d_to_fits(outfile_stack, wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack,
                              header=header, ex_value=ex_value, overwrite=True)
    if merge_stack:
        outfile_merge = 'spec1d_merge_{:}'.format(outfile)
        save.save_coadd1d_to_fits(outfile_merge, wave_merge, flux_merge, ivar_merge, mask_merge, header=header,
                                  ex_value=ex_value, overwrite=True)

    return wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack


def combspec_old(wave_grid, waves, fluxes, ivars, masks, ref_percentile=30.0, maxiter_scale=5, sigrej_scale=3,
             scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, dv_smooth=10000.0,
             const_weights=False, maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0, maxrej=None, debug_scale=False,
             debug=False, title=''):

    '''
    Routine for optimally combining long or multi-slit spectra or echelle spectra of individual orders. It will
    compute a stacked spectrum from a set of exposures on the specified wave_grid with proper treatment of
    weights and masking. This code calls the stacking code compute_stack, which uses np.histogram to combine the data using
    NGP and does not perform any interpolations and thus does not correlate errors. It uses wave_grid to determine the set
    of wavelength bins that the data are averaged on. The final spectrum will be on an ouptut wavelength grid which is not
    the same as wave_grid. The ouput wavelength grid is the weighted average of the individual wavelengths used for each
    exposure that fell into a given wavelength bin in the input wave_grid. This 1d coadding routine thus maintains the
    independence of the errors for each pixel in the combined spectrum and computes the weighted averaged wavelengths of
    each pixel in an analogous way to the 2d extraction procedure which also never interpolates to avoid correlating
    erorrs. It performs a number of iterations where it combines the spectra and performs rejection of outlier pixels
    using the spec_reject_comb code. The outliers are rejected using the true noise of the individual exposures, but
    uses the distribution of the pixel values about the stack to apply correction factors to the errors before rejecting.
    These corrected errors are currently only used in rejection but are not applied to the data.  This code is based
    on the xidl long_combpsec.pro routine but with significant improvements.

    Args:
        wave_grid: ndarray, (ngrid +1,)
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the new_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.
        waves: ndarray, (nspec, nexp)
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
        fluxes: ndarray, (nspec, nexp)
            fluxes for each exposure on the waves grid
        ivars: ndarray, (nspec, nexp)
            Inverse variances for each exposure on the waves grid
        masks: ndarray, bool, (nspec, nexp)
            Masks for each exposure on the waves grid. True=Good.
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
                definition of sticky.
        sigrej_scale: float, default=3.0
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        ref_percentile: float, default=20.0
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
        maxiter_scale: int, default=5
            Maximum number of iterations performed for rescaling spectra.
        scale_method: scale method, str, default=None. Options are poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default is to let the
            code determine this automitically which works well.
        dv_smooth: float, 10000.0
            Velocity smoothing used for determining smoothly varying S/N ratio weights by sn_weights
        hand_scale:
            array of hand scale factors, not well tested
        sn_max_medscale (float): default=2.0
            maximum SNR for perforing median scaling
        sn_min_medscale (float): default=0.5
            minimum SNR for perforing median scaling
        debug_scale (bool): default=False
            show interactive QA plots for the rescaling of the spectra
        title (str):
            Title prefix for spec_reject_comb QA plots
        debug (bool): default=False
            show interactive QA plot

    Returns:
        wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused, weights, scales, rms_sn

        wave_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
        outmask: ndarray, bool, (nspec, nexp)
             Output mask indicating which pixels are rejected in each exposure of the original input spectra after
             performing all of the iterations of combine/rejection
        nused: ndarray, (ngrid,)
             Numer of exposures which contributed to each pixel in the wave_stack. Note that this is in general
             different from nexp because of masking, but also becuse of the sampling specified by wave_grid. In other
             words, sometimes more spectral pixels in the irregularly gridded input wavelength array waves will land in
             one bin versus another depending on the sampling.
        weights: ndarray, (nspec, nexp)
            Weights used for combining your spectra which are computed using sn_weights
        scales: ndarray, (nspec, nexp)
            Scale factors applied to each individual spectrum before the combine computed by scale_spec
        rms_sn: ndarray, (nexp,)
            Root mean square S/N ratio of each of your individual exposures computed by sn_weights
    '''

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves,fluxes,ivars, masks, dv_smooth=dv_smooth, const_weights=const_weights, verbose=True)

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
            sigrej=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale, sn_max_medscale=sn_max_medscale,
            sn_min_medscale=sn_min_medscale, debug=debug_scale)

    # TODO Move this out of this routine and into the routine that does the actual coadd?
    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves, fluxes_scale, ivars_scale, masks, weights, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug, title=title)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused, weights, scales, rms_sn


def new_wave_grid(waves, wave_method='iref',iref=0, wave_grid_min=None, wave_grid_max=None,
                  A_pix=None,v_pix=None,samp_fact=1.0):
    """
    Create a new wavelength grid for the spectra to be rebinned and coadded on

    Args:
        waves (ndarray): (nspec, nexp,)
            Set of N original wavelength arrays
        wave_method (str): optional
            Desired method for creating new wavelength grid.
            'iref' -- Use the first wavelength array (default)
            'velocity' -- Constant velocity
            'pixel' -- Constant pixel grid
            'concatenate' -- Meld the input wavelength arrays
        iref (int): optional
            Reference spectrum
        wave_grid_min (float): optional
            min wavelength value for the final grid
        wave_grid_max (float): optional
            max wavelength value for the final grid
        A_pix (float):
            Pixel size in same units as input wavelength array (e.g. Angstroms)
            If not input, the median pixel size is calculated and used
        v_pix (float):
            Pixel size in km/s for velocity method
            If not input, the median km/s per pixel is calculated and used
        samp_fact (float):
            sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
            samp_fact < 1.0 undersamples (coarser)

    Returns
        wave_grid (ndarray): New wavelength grid, not masked
    """

    wave_mask = waves>1.0

    if wave_grid_min is None:
        wave_grid_min = np.min(waves[wave_mask])
    if wave_grid_max is None:
        wave_grid_max = np.max(waves[wave_mask])

    if wave_method == 'velocity':  # Constant km/s
        if v_pix is None:
            # Find the median velocity of a pixel in the input
            dv = c_kms * np.diff(waves, axis=0)/waves[1:]   # km/s
            v_pix = np.median(dv)

        # to make the wavelength grid finer or coarser
        v_pix = v_pix/samp_fact

        # Generate wavelength array
        x = np.log10(v_pix/c_kms + 1.0)
        npix = int(np.log10(wave_grid_max/wave_grid_min)/x) + 1
        wave_grid = wave_grid_min * 10.0**(x*np.arange(npix))

    elif wave_method == 'pixel': # Constant Angstrom
        if A_pix is None:
            dA =  np.abs(waves - np.roll(waves,1,axis=0))
            A_pix = np.median(dA)

        # Generate wavelength array
        wave_grid = wvutils.wavegrid(wave_grid_min, wave_grid_max + A_pix, \
                                     A_pix,samp_fact=samp_fact)

    elif wave_method == 'loggrid':
        dloglam_n = np.log10(waves) - np.roll(np.log10(waves), 1,axis=0)
        logwave_mask = wave_mask & np.roll(wave_mask, 1, axis=0)
        dloglam = np.median(dloglam_n[logwave_mask])
        wave_grid_max = np.max(waves[wave_mask])
        wave_grid_min = np.min(waves[wave_mask])
        # TODO: merge  wvutils.wavegrid with this function
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

def renormalize_errors_qa(chi, maskchi, sigma_corr, sig_range = 6.0, title='', qafile=None):
    '''
    Histogram QA plot of your chi distribution.
    Args:
        chi (ndarray): your chi values
        maskchi (ndarray, bool): True = good, mask for your chi array
        sigma_corr (float): corrected sigma
        sig_range (float): used for set binsize, default 6-sigma
        title (str):  plot title
        qafile (str or None): output QA file name
    Return:
        None
    '''

    n_bins = 50
    binsize = 2.0*sig_range/n_bins
    bins_histo = -sig_range + np.arange(n_bins)*binsize+binsize/2.0

    xvals = np.arange(-10.0,10,0.02)
    gauss = scipy.stats.norm(loc=0.0,scale=1.0)
    gauss_corr = scipy.stats.norm(loc=0.0,scale=sigma_corr)

    plt.figure(figsize=(12, 8))
    plt.hist(chi[maskchi],bins=bins_histo,normed=True,histtype='step', align='mid',color='k',linewidth=3,label='Chi distribution')
    plt.plot(xvals,gauss.pdf(xvals),'c-',lw=3,label='sigma=1')
    plt.plot(xvals,gauss_corr.pdf(xvals),'m--',lw=2,label='new sigma={:4.2f}'.format(round(sigma_corr,2)))
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
    '''
    Function for renormalizing errors. The distribbution of input chi (defined by chi = (data - model)/sigma) values is
    analyzed, and a correction factor to the standard deviation sigma_corr is returned. This should be multiplied into
    the errors. In this way, a rejection threshold of i.e. 3-sigma, will always correspond to roughly the same percentile.
    This renormalization guarantees that rejection is not too agressive in cases where the empirical errors determined
    from the chi-distribution differ significantly from the noise model which was used to determine chi.

    Args:
        chi (ndarray):
            input chi values
        mask (ndarray, bool):
            True = good, mask for your chi array
        clip (float):
            threshold for outliers which will be clipped for the purpose of computing the renormalization factor
        max_corr (float):
            maximum corrected sigma allowed.
        title (str):
            title for QA plot, will parsed to renormalize_errors_qa
        debug (bool):
            whether or not show the QA plot created by renormalize_errors_qa
    Returns:
        sigma_corr (float): corrected new sigma
        maskchi (ndarray, bool): new mask (True=good) which indicates the values used to compute the correction (i.e it includes clipping)
    '''

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
    Function for computing the chi^2 loss function for solving for the polynomial rescaling of one spectrum to another.
    There are two non-standard things implemented here which increase ther robustness. The first is a non-standard error used for the
    chi, which adds robustness and increases the stability of the optimization. This was taken from the idlutils
    solve_poly_ratio code. The second thing is that the chi is remapped using the scipy huber loss function to
    reduce sensitivity to outliers, ased on the scipy cookbook on robust optimization.


    Args:
        theta (ndarray): parameter vector for the polymomial fit
        flux_ref (ndarray): reference flux that data will be rescaled to match
        thismask (ndarray, bool): mask for the current iteration of the optimization, True=good
        arg_dict (dict): dictionary containing arguments

    Returns:
        loss_function (float): this is effectively the chi^2, i.e. the quantity to be minimized by the optimizer. Note
                       that this is not formally the chi^2 since the huber loss function re-maps the chi to be less
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
    # Evaluate the polynomial for rescaling
    ymult = (utils.func_val(theta, wave, func, minx=wave_min, maxx=wave_max))**2
    flux_scale = ymult*flux_med
    mask_both = mask & thismask
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
    '''
    Function to be optimized by robust_optimize for solve_poly_ratio polynomial rescaling of one spectrum to
    match a reference spectrum. This function has the correct format for running robust_optimize optimization. In addition
    to running the optimization, this function recomputes the error vector ivartot for the error rejection that takes
    place at each iteration of the robust_optimize optimization. The ivartot is also renormalized using the
    renormalize_errors function enabling rejection. A scale factor is multiplied into the true errors to allow one
    to reject based on the statistics of the actual error distribution.

    Args:
        flux_ref: ndarray, reference flux that we are trying to rescale our spectrum to match
        thismask: ndarray, bool, mask for the current iteration of the optimization. True=good
        arg_dict: dictionary containing arguments for the optimizing function. See poly_ratio_fitfunc_chi2 for how arguments
                  are used. They are mask, flux_med, flux_ref_med, ivar_ref_med, wave, wave_min, wave_max, func
        kwargs_opt: arguments to be passed to the optimizer, which in this case is just vanilla scipy.minimize with
                    the default optimizer
    Return:
        result: scipy optimization object
        flux_scale: scale factor to be applied to the data to match the reference spectrum flux_ref
        ivartot: error vector to be used for the rejection that takes place at each iteration of the robust_optimize
                 optimization
    '''

    # flux_ref, ivar_ref act like the 'data', the rescaled flux will be the 'model'

    guess = arg_dict['guess']
    result = scipy.optimize.minimize(poly_ratio_fitfunc_chi2, guess, args=(flux_ref, thismask, arg_dict),  **kwargs_opt)
    flux = arg_dict['flux']
    ivar = arg_dict['ivar']
    mask = arg_dict['mask']
    ivar_ref = arg_dict['ivar_ref']
    wave = arg_dict['wave']
    wave_min = arg_dict['wave_min']
    wave_max = arg_dict['wave_max']
    func = arg_dict['func']
    # Evaluate the polynomial for rescaling
    ymult = (utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max))**2
    flux_scale = ymult*flux
    mask_both = mask & thismask
    totvar = utils.inverse(ivar_ref) + ymult**2*utils.inverse(ivar)
    ivartot1 = mask_both*utils.inverse(totvar)
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
    '''
    Utility routine to median filter a spectrum using the mask and propagating the errors using the
    utils.fast_running_median function.
    Args:
        flux: ndarray, (nspec,) flux
        ivar: ndarray, (nspec,) inverse variance
        mask: ndarray, bool, (nspec,) True = good
        med_width: width for median filter in pixels
    Return:
        flux_med, ivar_med: median filtered flux and corresponding propagated errors
    '''

    flux_med = np.zeros_like(flux)
    ivar_med = np.zeros_like(ivar)
    flux_med0 = utils.fast_running_median(flux[mask], med_width)
    flux_med[mask] = flux_med0
    var = utils.inverse(ivar)
    var_med0 =  utils.smooth(var[mask], med_width)
    ivar_med[mask] = utils.inverse(var_med0)
    return flux_med, ivar_med

def solve_poly_ratio(wave, flux, ivar, flux_ref, ivar_ref, norder, mask = None, mask_ref = None,
                     scale_min = 0.05, scale_max = 100.0, func='legendre',
                     maxiter=3, sticky=True, lower=3.0, upper=3.0, median_frac=0.01, debug=False):
    '''
    Routine for solving for the polynomial rescaling of an input spectrum flux to match a reference spectrum flux_ref.
    The two spectra need to be defined on the same wavelength grid. The code will work best if you choose the reference
    to be the higher S/N ratio spectrum. Note that the code
    multiplies in the square of a polnomial of order norder to ensure positivity of the scale factor.  It also
    operates on median filtered spectra to be more robust against outliers

    Args:
        wave: ndarray, (nspec,)
            wavelength. flux, ivar, flux_ref, and ivar_ref must all be on the same wavelength grid
        flux: ndarray, (nspec,)
             flux that you want to rescale to match flux_ref
        mask: ndarray, bool, (nspec,)
            mask for spectrum that you want to rescale, True=Good
        flux_ref: ndarray, (nspec,)
            reference flux that you want to rescale flux to match.
        ivar_ref: ndarray, (nspec,)
            inverse variance for reference flux
        mask_ref: ndarray, bool (nspec,)
            mask for reference flux
        norder: int
            order of polynomial rescaling.  Note that the code multiplies in by the square of a polynomail of order
            norder to ensure positivity of the scale factor.
        scale_min: float, default =0.05
             minimum scaling factor allowed
        scale_max: float, default=100.0
             maximum scaling factor allowed
        func: str, default='legendre'
             function you want to use,
        maxiter: int, default=3
             maximum number of iterations for robust_optimize
        sticky: bool, default=True
             whether you want the rejection to be sticky or not with robust_optimize. See docs for djs_reject for
             definition of sticky.
        lower: float, default=3.0
             lower sigrej rejection threshold for robust_optimize
        upper: float, default=3.0
             upper sigrej rejection threshold for robust_optimize
        median_frac: float default = 0.01,
             the code rescales median filtered spectra with 'reflect' boundary conditions. The
              with of the median filter will be median_frac*nspec, where nspec is the number of spectral pixels.
        debug: bool, default=False
             show interactive QA plot
    Return:
        ymult, flux_rescale, ivar_rescale, outmask

        ymult: ndarray, (nspec,)
            rescaling factor to be multiplied into flux to match flux_ref
        flux_rescale: ndarray, (nspec,)
            rescaled flux, i.e. ymult multiplied into flux
        ivar_rescale: ndarray, (nspec,)
            rescaled inverse variance
        outmask: ndarray, bool, (nspec,)
            output mask determined from the robust_optimize optimization/rejection iterations. True=Good
    '''

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
                    wave_max = wave_max, func = func, norder = norder, guess = guess, debug=False)

    result, ymodel, ivartot, outmask = utils.robust_optimize(flux_ref, poly_ratio_fitfunc, arg_dict, inmask=mask_ref,
                                                             maxiter=maxiter, lower=lower, upper=upper, sticky=sticky)
    ymult1 = (utils.func_val(result.x, wave, func, minx=wave_min, maxx=wave_max))**2
    ymult = np.fmin(np.fmax(ymult1, scale_min), scale_max)
    flux_rescale = ymult*flux
    ivar_rescale = ivar/ymult**2

    if debug:
        # Determine the y-range for the QA plots
        scale_spec_qa(wave, flux_med, ivar_med, wave, flux_ref_med, ivar_ref_med, ymult, 'poly', mask = mask, mask_ref=mask_ref,
                      title='Median Filtered Spectra that were poly_ratio Fit')

    return ymult, flux_rescale, ivar_rescale, outmask


def interp_oned(wave_new, wave_old, flux_old, ivar_old, mask_old):
    '''
    Utility routine to perform 1d linear nterpolation of spectra onto a new wavelength grid

    Args:
       wave_new: ndarray, (nspec_new)
            New wavelengths that you want to interpolate onto.
       wave_old: ndarray, (nspec_old)
            Old wavelength grid
       flux_old: ndarray, (nspec_old)
            Old flux on the wave_old grid
       ivar_old: ndarray, (nspec_old)
            Old ivar on the wave_old grid
       mask_old: ndarray, bool, (nspec_old),
            Old mask on the wave_old grid. True=Good
    Returns :
       flux_new, ivar_new, mask_new

       flux_new: ndarray, (nspec_new,)
            interpolated flux
       ivar_new: ndarray, (nspec_new,)
            interpolated ivar
       mask_new: ndarray, bool, (nspec_new,)
            interpolated mask. True=Good
    '''

    # Do not interpolate if the wavelength is exactly same with wave_new
    if np.array_equal(wave_new, wave_old):
        return flux_old, ivar_old, mask_old

    # make the mask array to be float, used for interpolation
    masks_float = mask_old.astype(float)
    wave_mask = wave_old > 1.0 # Deal with the zero wavelengths
    flux_new = scipy.interpolate.interp1d(wave_old[wave_mask], flux_old[wave_mask], kind='cubic',
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    ivar_new = scipy.interpolate.interp1d(wave_old[wave_mask], ivar_old[wave_mask], kind='cubic',
                                    bounds_error=False, fill_value=np.nan)(wave_new)
    mask_new_tmp = scipy.interpolate.interp1d(wave_old[wave_mask], masks_float[wave_mask], kind='cubic',
                                        bounds_error=False, fill_value=np.nan)(wave_new)
    # Don't allow the ivar to be every less than zero
    ivar_new = (ivar_new > 0.0)*ivar_new
    mask_new = (mask_new_tmp > 0.8) & (ivar_new > 0.0) & np.isfinite(flux_new) & np.isfinite(ivar_new)
    return flux_new, ivar_new, mask_new

def interp_spec(wave_new, waves, fluxes, ivars, masks):
    '''
    Utility routine to interpolate a set of spectra onto a new wavelength grid, wave_new
    Args:
        wave_new: ndarray, shape (nspec,) or (nspec, nimgs),
             new wavelength grid
        waves:  ndarray, shape (nspec,) or (nspec, nexp)
             where nexp, need not equal nimgs. Old wavelength grids
        fluxes: ndarray,
             same shape as waves, old flux
        ivars: ndarray,
             same shape as waves, old ivar
        masks: ndarray, bool,
             same shape as waves, old mask, True=Good
    Returns:
        fluxes_inter, ivars_inter, masks_inter

        Interpolated flux, ivar and mask with the size and shape matching wave_new. masks_inter is bool with True=Good
    '''

    # First case: interpolate either an (nspec, nexp) array of spectra onto a single wavelength grid
    if (wave_new.ndim == 1):
        if fluxes.ndim == 1:
            fluxes_inter, ivars_inter, masks_inter = interp_oned(wave_new, waves, fluxes, ivars, masks)
        else:
            nexp = fluxes.shape[1]
            # Interpolate spectra to have the same wave grid with the iexp spectrum.
            # And scale spectra to the same flux level with the iexp spectrum.
            fluxes_inter = np.zeros((wave_new.size, nexp))
            ivars_inter  = np.zeros((wave_new.size, nexp))
            masks_inter  = np.zeros((wave_new.size, nexp), dtype=bool)
            for ii in range(nexp):
                fluxes_inter[:, ii], ivars_inter[:, ii], masks_inter[:, ii] = interp_oned(
                    wave_new, waves[:, ii], fluxes[:, ii], ivars[:, ii], masks[:, ii])

        return fluxes_inter, ivars_inter, masks_inter

    # Second case: interpolate a single spectrum onto an (nspec, nexp) array of wavelengths
    elif (wave_new.ndim == 2):
        if fluxes.ndim != 1:
            msgs.error('If wave_new is two dimensional, all other input arrays must be one dimensional')
        nexp = wave_new.shape[1]
        fluxes_inter = np.zeros_like(wave_new)
        ivars_inter = np.zeros_like(wave_new)
        masks_inter = np.zeros_like(wave_new, dtype=bool)

        for ii in range(nexp):
            fluxes_inter[:, ii], ivars_inter[:, ii], masks_inter[:, ii] = interp_oned(
                wave_new[:, ii], waves, fluxes, ivars, masks)

        return fluxes_inter, ivars_inter, masks_inter

    else:
        msgs.error('Invalid size for wave_new')

def sn_weights(waves, fluxes, ivars, masks, dv_smooth=10000.0, const_weights=False, ivar_weights=False, verbose=False):

    """ Calculate the S/N of each input spectrum and create an array of (S/N)^2 weights to be used for coadding.

    Args:
    ----------
    fluxes: float ndarray, shape = (nspec, nexp)
        Stack of (nspec, nexp) spectra where nexp = number of exposures, and nspec is the length of the spectrum.
    ivars: float ndarray, shape = (nspec, nexp)
        Inverse variance noise vectors for the spectra
    masks: bool ndarray, shape = (nspec, nexp)
        Mask for stack of spectra. True=Good, False=Bad.
    waves: flota ndarray, shape = (nspec,) or (nspec, nexp)
        Reference wavelength grid for all the spectra. If wave is a 1d array the routine will assume
        that all spectra are on the same wavelength grid. If wave is a 2-d array, it will use the individual
    dv_smooth: float, optional, default = 10000.0
         Velocity smoothing used for determining smoothly varying S/N ratio weights.

    Returns
    -------
    rms_sn : ndarray, shape (nexp)
        Root mean square S/N value for each input spectra
    weights : ndarray, shape = (nspec, nexp)
        Weights to be applied to the spectra. These are signal-to-noise squared weights.
    """

    sigs = np.sqrt(utils.inverse(ivars))

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
    sn_val_ma = np.ma.array(sn_val, mask = np.invert(mask_stack))
    sn_sigclip = stats.sigma_clip(sn_val_ma, sigma=3, maxiters=5)
    ## TODO Update with sigma_clipped stats with our new cenfunc and std_func = mad_std
    sn2 = (sn_sigclip.mean(axis=0).compressed())**2 #S/N^2 value for each spectrum
    rms_sn = np.sqrt(sn2) # Root Mean S/N**2 value for all spectra
    rms_sn_stack = np.sqrt(np.mean(sn2))

    # TODO: ivar weights is better than SN**2 or const_weights for merging orders. Enventially, we will change it to
    # TODO Should ivar weights be deprecated??
    if ivar_weights:
        if verbose:
            msgs.info("Using sensfunc weights for merging orders")
        #weights = ivar_stack
        weights = np.ones_like(flux_stack) # Should this be zeros_like?
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
            # TODO Change to scipy.interpolate
            sn_med2 = np.interp(spec_vec, spec_now, sn_med1)
            sig_res = np.fmax(med_width/10.0, 3.0)
            gauss_kernel = convolution.Gaussian1DKernel(sig_res)
            sn_conv = convolution.convolve(sn_med2, gauss_kernel)
            weights[:, iexp] = sn_conv
    elif rms_sn_stack <= 3.0 or const_weights:
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
            sn_med2 = np.interp(spec_vec, spec_now, sn_med1)
            ##sn_med2 = np.interp(wave_stack[iexp,:], wave_now,sn_med1)
            sig_res = np.fmax(med_width/10.0, 3.0)
            gauss_kernel = convolution.Gaussian1DKernel(sig_res)
            sn_conv = convolution.convolve(sn_med2, gauss_kernel)
            weights[:, iexp] = sn_conv
            if verbose:
                msgs.info('S/N = {:4.2f}, averaged weight = {:4.2f} for {:}th exposure'.format(
                    rms_sn[iexp],np.mean(weights[:, iexp]), iexp))

    if fluxes.ndim == 3:
        rms_sn = np.reshape(rms_sn, (norder, nexp), order='F')
        weights = np.reshape(weights, (nspec, norder, nexp), order='F')

    # Finish
    return rms_sn, weights

# TODO Rename this function to something sensfunc related
def get_tell_from_file(sensfile, waves, masks, iord=None):
    '''
    Get the telluric model from the sensfile.
    Args:
        sensfile (str): the name of your fits format sensfile
        waves (ndarray): wavelength grid for your output telluric model
        masks (ndarray, bool): mask for the wave
        iord (int or None): if None returns telluric model for all orders, otherwise return the order you want
    Returns:
         telluric (ndarray): telluric model on your wavelength grid
    '''


    sens_param = Table.read(sensfile, 1)
    sens_table = Table.read(sensfile, 2)
    telluric = np.zeros_like(waves)

    if (waves.ndim == 1) and (iord is None):
        msgs.info('Loading Telluric from Longslit sensfiles.')
        tell_interp = scipy.interpolate.interp1d(sens_table[0]['WAVE'], sens_table[0]['TELLURIC'], kind='cubic',
                                        bounds_error=False, fill_value=np.nan)(waves[masks])
        telluric[masks] = tell_interp
    elif (waves.ndim == 1) and (iord is not None):
        msgs.info('Loading order {:} Telluric from Echelle sensfiles.'.format(iord))
        wave_tell_iord = sens_table[iord]['WAVE']
        tell_mask = (wave_tell_iord > 1.0)
        tell_iord = sens_table[iord]['TELLURIC']
        tell_iord_interp = scipy.interpolate.interp1d(wave_tell_iord[tell_mask], tell_iord[tell_mask], kind='cubic',
                                        bounds_error=False, fill_value=np.nan)(waves[masks])
        telluric[masks] = tell_iord_interp
    else:
        norder = np.shape(waves)[1]
        for iord in range(norder):
            wave_iord = waves[:, iord]
            mask_iord = masks[:, iord]

            # Interpolate telluric to the same grid with waves
            # Since it will be only used for plotting, I just simply interpolate it rather than evaluate it based on the model
            wave_tell_iord = sens_table[iord]['WAVE']
            tell_mask = (wave_tell_iord > 1.0)
            tell_iord = sens_table[iord]['TELLURIC']
            tell_iord_interp = scipy.interpolate.interp1d(wave_tell_iord[tell_mask], tell_iord[tell_mask], kind='cubic',
                                                    bounds_error=False, fill_value=np.nan)(wave_iord[mask_iord])
            telluric[mask_iord, iord] = tell_iord_interp

    return telluric


def sensfunc_weights(sensfile, waves, masks, debug=False):
    '''
    Get the weights based on the sensfunc
    Args:
        sensfile (str): the name of your fits format sensfile
        waves (ndarray): wavelength grid for your output weights
        masks (ndarray, bool): mask for the wave s
        debug (bool): whether you want show the weights QA
    Returns:
        weights (ndarray): weights on you wavelength grid
        masks (ndarray, bool): mask for your weights
    '''

    sens_meta= Table.read(sensfile, 1)
    sens_table = Table.read(sensfile, 2)
    func = sens_meta['FUNC'][0]
    polyorder_vec = sens_meta['POLYORDER_VEC'][0]

    weights = np.zeros_like(waves)
    norder = waves.shape[1]

    if norder != len(sens_table):
        msgs.error('The number of orders in {:} does not agree with your data. Wrong sensfile?'.format(sensfile))

    for iord in range(norder):
        wave_iord = waves[:, iord]
        mask_iord = masks[:, iord]
        wave_mask = wave_iord > 1.0

        # get sensfunc from the sens_table
        coeff = sens_table[iord]['OBJ_THETA'][0:polyorder_vec[iord] + 2]
        wave_min=sens_table[iord]['WAVE_MIN']
        wave_max = sens_table[iord]['WAVE_MAX']
        sensfunc_iord = np.exp(utils.func_val(coeff, wave_iord[wave_mask], func, minx=wave_min, maxx=wave_max))
        mask_sens_iord = sensfunc_iord > 0.0
        weights[mask_iord, iord] = 1.0/(sensfunc_iord+(sensfunc_iord==0.))
        masks[mask_iord, iord] = mask_sens_iord

    if debug:
        weights_qa(waves, weights, masks)


    return weights, masks

def robust_median_ratio(flux, ivar, flux_ref, ivar_ref, mask=None, mask_ref=None, ref_percentile=20.0, min_good=0.05,
                        maxiters=5, sigrej=3.0, max_factor=10.0):
    '''
    Robustly determine the ratio between input spectrum flux and reference spectrum flux_ref. The code will perform
    best if the reference spectrum is chosen to be the higher S/N ratio spectrum, i.e. a preliminary stack that you want
    to scale each exposure to match. Note that the flux and flux_ref need to be on the same wavelength grid!!

    Args:
        wave: ndarray, (nspec,)
            wavelengths grid for the spectra
        flux: ndarray, (nspec,)
            spectrum that will be rescaled.
        ivar: ndarray, (nspec,)
            inverse variance for the spectrum that will be rescaled.
        mask: ndarray, bool, (nspec,)
            mask for the spectrum that will be rescaled. True=Good. If not input, computed from inverse variance
        flux_ref: ndarray, (nspec,)
            reference spectrum.
        ivar_ref: ndarray, (nspec,)
            inverse variance of reference spectrum.
        mask_ref: ndarray, bool, (nspec,)
            mask for reference spectrum. True=Good. If not input, computed from inverse variance.
        ref_percentile: float, default=20.0
            Percentile fraction used for selecting the minimum SNR cut. Pixels above this cut are deemed the "good"
            pixels and are used to compute the ratio. This must be a number between 0 and 100.
        min_good: float, default = 0.05
            Minimum fraction of good pixels determined as a fraction of the total pixels for estimating the median ratio
        maxiters: int, defrault = 5,
            Maximum number of iterations for astropy.stats.SigmaClip
        sigrej: float, default = 3.0
            Rejection threshold for astropy.stats.SigmaClip
        max_factor: float, default = 10.0,
            Maximum allowed value of the returned ratio
    Returns:
        ratio: float, the number that must be multiplied into flux in order to get it to match up with flux_ref
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
    Function for scaling different orders
    Args:
        waves (ndarray): wavelength array of your spectra with the shape of (nspec, norder)
        fluxes (ndarray): flux array of your spectra with the shape of (nspec, norder)
        ivars (ndarray): ivar array of your spectra with the shape of (nspec, norder)
        masks (ndarray, bool): mask for your spectra with the shape of (nspec, norder)
        min_good (float): minmum fraction of the total number of good pixels needed for estimate the median ratio
        maxiters (int or float): maximum iterations for rejecting outliers
        max_factor (float): maximum scale factor
        sigrej (float): sigma used for rejecting outliers
        debug (bool): if True show the QA
    Returns:
        fluxes_new (ndarray): re-scaled fluxes with the shape of (nspec, norder)
        ivars_new (ndarray): re-scaled ivars with the shape of (nspec, norder)
        order_ratios (ndarray): an array of scale factor with the length of norder
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


def scale_spec(wave, flux, ivar, wave_ref, flux_ref, ivar_ref, mask=None, mask_ref=None, scale_method=None, min_good=0.05,
               ref_percentile=20.0, maxiters=5, sigrej=3, max_median_factor=10.0,
               npoly=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, debug=False, show=False):
    '''
    Routine for solving for the best way to rescale an input spectrum flux to match a reference spectrum flux_ref.
    The two spectra need to be defined on the same wavelength grid. The code will work best if you choose the reference
    to be the higher S/N ratio spectrum. If the scale_method is not specified, the code will make a decision about
    which method to use based on the S/N ratio of the input spectrum flux.

    Args:
    wave: ndarray, (nspec,)
       wavelengths grid for the spectra
    flux: ndarray, (nspec,)
       spectrum that will be rescaled.
    ivar: ndarray, (nspec,)
       inverse variance for the spectrum that will be rescaled.
    mask: ndarray, bool, (nspec,)
       mask for the spectrum that will be rescaled. True=Good. If not input, computed from inverse variance
    flux_ref: ndarray, (nspec,)
       reference spectrum.
    ivar_ref: ndarray, (nspec,)
       inverse variance of reference spectrum.
    mask_ref: ndarray, bool, (nspec,)
       mask for reference spectrum. True=Good. If not input, computed from inverse variance.
    min_good: float, default = 0.05
       minmum fraction of the total number of good pixels needed for estimate the median ratio
    maxiters: int,
       maximum number of iterations for rejecting outliers used by the robust_median_ratio routine if median
       rescaling is the  method used.
    max_median_factor: float, default=10.0
       maximum scale factor for median rescaling for robust_median_ratio if median rescaling is the method used.
    sigrej: float, default=3.0
       rejection threshold used for rejecting outliers by robsut_median_ratio
    ref_percentile: float, default=20.0
       percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
    npoly: int, default=None
       order for the poly ratio scaling if polynomial rescaling is the method used. Default is to automatically compute
       this based on S/N ratio of data.
    scale_method: scale method, str, default=None. Options are poly, median, none, or hand. Hand is not well tested.
                  User can optionally specify the rescaling method. Default is to let the
                  code determine this automitically which works well.
    hand_scale: ndarray, (nexp,)
        array of hand scale factors, not well tested
    sn_max_medscale: float, default=2.0
       maximum SNR for perforing median scaling
    sn_min_medscale: float, default=0.5
       minimum SNR for perforing median scaling
    debug: bool, default=False
       show interactive QA plot

    Returns:
        flux_scale: ndarray (nspec,) scaled spectrum
        ivar_scale: ndarray (nspec,) inverse variance for scaled spectrum
        scale: ndarray (nspec,) scale factor applied to the spectrum and inverse variance
        scale_method: str, method that was used to scale the spectra.
    '''

    if mask is None:
        mask = ivar > 0.0
    if mask_ref is None:
        mask_ref = ivar_ref > 0.0


    # Interpolate the reference spectrum onto the wavelengths of the spectrum that will be rescaled
    flux_ref_int, ivar_ref_int, mask_ref_int = interp_spec(wave, wave_ref, flux_ref, ivar_ref, mask_ref)

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
        scale, flux_scale, ivar_scale, outmask = solve_poly_ratio(wave, flux, ivar, flux_ref_int, ivar_ref_int, npoly,
                                                                      mask=mask, mask_ref=mask_ref_int, debug=debug)
    elif scale_method == 'median':
        # Median ratio (reference to spectrum)
        med_scale = robust_median_ratio(flux, ivar, flux_ref_int, ivar_ref_int,ref_percentile=ref_percentile,min_good=min_good,
                                        mask=mask, mask_ref=mask_ref_int, maxiters=maxiters,
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
    if show:
        scale_spec_qa(wave, flux, ivar, wave_ref, flux_ref, ivar_ref, scale, scale_method, mask = mask, mask_ref=mask_ref,
                      title='Scaling Applied to the Data')

    return flux_scale, ivar_scale, scale, scale_method


def compute_stack(wave_grid, waves, fluxes, ivars, masks, weights):
    '''
    Compute a stacked spectrum from a set of exposures on the specified wave_grid with proper treatment of
    weights and masking. This code uses np.histogram to combine the data using NGP and does not perform any
    interpolations and thus does not correlate errors. It uses wave_grid to determine the set of wavelength bins that
    the data are averaged on. The final spectrum will be on an ouptut wavelength grid which is not the same as wave_grid.
    The ouput wavelength grid is the weighted average of the individual wavelengths used for each exposure that fell into
    a given wavelength bin in the input wave_grid. This 1d coadding routine thus maintains the independence of the
    errors for each pixel in the combined spectrum and computes the weighted averaged wavelengths of each pixel
    in an analogous way to the 2d extraction procedure which also never interpolates to avoid correlating erorrs.

    Args:
        wave_grid: ndarray, (ngrid +1,)
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the new_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.
        waves: ndarray, (nspec, nexp)
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
        fluxes: ndarray, (nspec, nexp)
            fluxes for each exposure on the waves grid
        ivars: ndarray, (nspec, nexp)
            Inverse variances for each exposure on the waves grid
        masks: ndarray, bool, (nspec, nexp)
            Masks for each exposure on the waves grid. True=Good.
        weights: ndarray, (nspec, nexp)
            Weights to be used for combining your spectra. These are computed using sn_weights
    Returns:
        wave_stack, flux_stack, ivar_stack, mask_stack, nused

        wave_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
        nused: ndarray, (ngrid,)
             Numer of exposures which contributed to each pixel in the wave_stack. Note that this is in general
             different from nexp because of masking, but also becuse of the sampling specified by wave_grid. In other
             words, sometimes more spectral pixels in the irregularly gridded input wavelength array waves will land in
             one bin versus another depending on the sampling.
    '''

    ubermask = masks & (weights > 0.0) & (waves > 1.0) & (ivars > 0.0)
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
    wave_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=waves_flat*weights_flat)
    wave_stack = (weights_total > 0.0)*wave_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked flux
    flux_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=fluxes_flat*weights_flat)
    flux_stack = (weights_total > 0.0)*flux_stack_total/(weights_total+(weights_total==0.))

    # Calculate the stacked ivar
    var_stack_total, wave_edges = np.histogram(waves_flat,bins=wave_grid,density=False,weights=vars_flat*weights_flat**2)
    var_stack = (weights_total > 0.0)*var_stack_total/(weights_total+(weights_total==0.))**2
    ivar_stack = utils.inverse(var_stack)

    # New mask for the stack
    mask_stack = (weights_total > 0.0) & (nused > 0.0)

    return wave_stack, flux_stack, ivar_stack, mask_stack, nused

def get_ylim(flux, ivar, mask):
    """
    Utility routine for setting the plot limits for QA plots.
    Args:
        flux: ndarray, (nspec,) flux array
        ivar: ndarray, (nspec,) inverse variance array
        mask: ndarray, bool, (nspec,) mask array. True=Good

    Returns:
        ymin, ymax: limits for plotting.

    """

    med_width = (2.0 * np.ceil(0.1 / 2.0 * np.size(flux[mask])) + 1).astype(int)
    flux_med, ivar_med = median_filt_spec(flux, ivar, mask, med_width)
    mask_lim = ivar_med > np.percentile(ivar_med, 20)
    ymax = 2.5 * np.max(flux_med[mask_lim])
    ymin = -0.15 * ymax
    return ymin, ymax

def scale_spec_qa(wave, flux, ivar, wave_ref, flux_ref, ivar_ref, ymult, scale_method,
                  mask=None, mask_ref=None, ylim = None, title=''):
    '''
    QA plot for spectrum scaling.

    Args:
        wave: ndarray, (nspec,)
            wavelength array for spectrum to be scaled and reference spectrum.
        flux: ndarray, (nspec,)
            flux for spectrum to be scaled
        ivar: ndarray, (nspec,)
             inverse variance for spectrum to be scaled.
        mask: ndarray, bool, (nspec,) optional,
             mask for spectrum to be scaled. True=Good. If not specified determined form inverse variance
        flux_ref: ndarray (nspec,)
             reference flux
        ivar_ref: ndarray (nspec,)
            inverse variance of reference flux
        mask_ref: ndarray, bool, (nspec,)
            mask for reference flux. True=Good.
        ymult: ndarray (nspec,)
            scale factor array
        scale_method: str,
            method used for rescaling which will be shown on QA plot.
        ylim: ylim, default=None,
            tuple for limits of the QA plot. If None, will be determined automtically with get_ylim
        title: str, QA plot title
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
        ylim = get_ylim(flux, ivar, mask)

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
    wave_min = np.fmax(0.7*wave[wave_mask].min(), wave_ref[wave_mask_ref].min())
    wave_max = np.fmin(1.3*wave[wave_mask].max(), wave_ref[wave_mask_ref].max())
    poly_plot.set_xlim((wave_min, wave_max))
    spec_plot.set_xlim((wave_min, wave_max))
    spec_plot.set_ylim(ylim)

    spec_plot.set_xlim((wave_ref[wave_mask_ref].min(), wave_ref[wave_mask_ref].max()))
    spec_plot.plot(wave[wave_mask], flux[wave_mask], color='red', zorder=10,
                   marker='o', markersize=1.0, mfc='k', fillstyle='full', linestyle='None', label='original spectrum')
    spec_plot.plot(wave[wave_mask], flux[wave_mask]*ymult[wave_mask], color='dodgerblue', drawstyle='steps-mid', alpha=0.5, zorder=5, linewidth=2,
                   label='rescaled spectrum')
    spec_plot.plot(wave_ref[wave_mask_ref], flux_ref[wave_mask_ref], color='black', drawstyle='steps-mid', zorder=7, alpha = 0.5, label='reference spectrum')

    spec_plot.legend()
    fig.suptitle(title)
    plt.show()

def coadd_iexp_qa(wave, flux, mask, flux_stack, mask_stack, rejivar, outmask, norder=None, title='', qafile=None):
    '''
    Routine to creqate QA for showing the individual spectrum compared to the combined stacked spectrum indicating
     which pixels were rejected.
    Args:
        wave: ndarray, (nspec,)
            wavelength array for spectrum of the exposure in question.
        flux: ndarray, (nspec,)
            flux for the exposure in question
        ivar: ndarray, (nspec,)
             inverse variance for the exposure in question
        mask: ndarray, bool, (nspec,) optional,
             mask for the exposure in question True=Good. If not specified determined form inverse variance
        flux_stack: ndarray (nspec,)
             Stacked spectrum to be compared to the exposure in question.
        ivar_ref: ndarray (nspec,)
            inverse variance of the stacked spectrum
        mask_ref: ndarray, bool, (nspec,)
            mask for stacked spectrum
        norder: int, default=None, Indicate the number of orders if this is an echelle stack
        title (str):
            plot title
        qafile: QA file name

    '''


    fig = plt.figure(figsize=(12, 8))
    spec_plot = fig.add_axes([0.1, 0.1, 0.8, 0.85])

    # Get limits
    ymin, ymax = get_ylim(flux, rejivar, (mask & mask_stack & np.isfinite(rejivar)))

    # Plot spectrum
    rejmask = mask & np.invert(outmask)
    wave_mask = wave > 1.0
    spec_plot.plot(wave[rejmask], flux[rejmask],'s',zorder=10,mfc='None', mec='r', label='rejected pixels')
    spec_plot.plot(wave[np.invert(mask)], flux[np.invert(mask)],'v', zorder=10, mfc='None', mec='orange',
                   label='originally masked')

    if norder is None:
        spec_plot.plot(wave[wave_mask], flux[wave_mask], color='dodgerblue', linestyle='steps-mid',
                       zorder=2, alpha=0.5,label='single exposure')
        spec_plot.plot(wave[wave_mask], np.sqrt(utils.inverse(rejivar[wave_mask])),zorder=3,
                       color='0.7', alpha=0.5, linestyle='steps-mid')
        spec_plot.plot(wave[wave_mask],flux_stack[wave_mask]*mask_stack[wave_mask],color='k',
                       linestyle='steps-mid',lw=2,zorder=3, alpha=0.5, label='coadd')

        # TODO Use one of our telluric models here instead
        # Plot transmission
        if (np.max(wave[mask]) > 9000.0):
            skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
            skycat = np.genfromtxt(skytrans_file, dtype='float')
            scale = 0.8 * ymax
            spec_plot.plot(skycat[:, 0] * 1e4, skycat[:, 1] * scale, 'm-', alpha=0.5, zorder=11)
    else:
        npix = np.size(flux)
        nspec = int(npix / norder)
        for iord in range(norder):
            spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           flux[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           linestyle='steps-mid', zorder=1, alpha=0.7)
            #spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
            #               np.sqrt(utils.inverse(ivar[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]], positive=True)),
            #               zorder=3, color='0.7', linestyle='steps-mid')
            spec_plot.plot(wave[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           flux_stack[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]]*
                           mask_stack[nspec*iord:nspec*(iord+1)][wave_mask[nspec*iord:nspec*(iord+1)]],
                           color='k', linestyle='steps-mid',lw=1,zorder=2)

    # properties
    spec_plot.legend(fontsize=13)
    spec_plot.set_ylim([ymin, ymax])
    spec_plot.set_xlim([wave[wave_mask].min(), wave[wave_mask].max()])
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

    Args:
        wave: ndarray, (nspec, nexp)
            wavelength array for spectra that went into a stack
        weights: ndarray, (nspec, nexp,)
            (S/N)^2 weights for the exposures that went into a stack. This would have been computed by sn_weights
        mask: ndarray, bool, (nspec, nexp)
            Pixels which were masked in each individual exposure which go into the stack.
    '''

    nexp = np.shape(waves)[1]
    fig = plt.figure(figsize=(12, 8))
    wave_mask = waves > 1.0
    for iexp in range(nexp):
        plt.plot(waves[wave_mask[:,iexp],iexp], weights[wave_mask[:,iexp],iexp]*masks[wave_mask[:,iexp],iexp])
    plt.xlim((waves[wave_mask].min(), waves[wave_mask].max()))
    plt.xlabel('Wavelength (Angstrom)')
    plt.ylabel('Weights')
    plt.title(title, fontsize=16, color='red')
    plt.show()

def coadd_qa(wave, flux, ivar, nused, mask=None, tell=None, title=None, qafile=None):
    '''
    Routine to make QA plot of the final stacked spectrum. It works for both longslit/mulitslit, coadded individual
    order spectrum of the Echelle data and the final coadd of the Echelle data.

    Args:
        wave: ndarray, (nspec,)
            one-d wavelength array of your spectrum
        flux: ndarray, (nspec,)
            one-d flux array of your spectrum
        ivar: ndarray, (nspec,)
            one-d ivar array of your spectrum
        mask: ndarray, bool (nspec,)
            mask array for your spectrum
        nused: ndarray, (nspec,)
            how many exposures used in the stack for each pixel, the same size with flux
        title: str
            plot title
        qafile: str
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
    num_plot.plot(wave[wave_mask],nused[wave_mask],linestyle='steps-mid',color='k',lw=2)
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
        skytrans_file = resource_filename('pypeit', '/data/skisim/atm_transmission_secz1.5_1.6mm.dat')
        skycat = np.genfromtxt(skytrans_file,dtype='float')
        scale = 0.8*ymax
        spec_plot.plot(skycat[:,0]*1e4,skycat[:,1]*scale,'m-',alpha=0.5,zorder=11)
    elif (tell is not None):
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

    return

def update_errors(fluxes, ivars, masks, fluxes_stack, ivars_stack, masks_stack, sn_clip=30.0, title='', debug=False):
    '''
    Deterimine corrections to errors using the residuals of each exposure about a preliminary stack. This routine is
    used as part of the iterative masking/stacking loop to determine the corrections to the errors used to reject pixels
    for the next iteration of the stack. The routine returns a set of corrections for each of the exposores that is input.
    Args:
        fluxes (ndarray): (nspec, nexp)
            fluxes for each exposure on the native wavelength grids
        ivars (ndarray): (nspec, nexp)
            Inverse variances for each exposure on the native wavelength grids
        masks (ndarray): bool, (nspec, nexp)
            Masks for each exposure on the native wavelength grids. True=Good.
        fluxes_stack (ndarray): (nspec, nexp)
            Stacked spectrum for this iteration interpolated on the native wavelength grid of the fluxes exposures.
        ivars_stack (ndarray): (nspec, nexp)
            Inverse variances of stacked spectrum for this iteration interpolated on the native wavelength grid of the
            fluxes exposures.
        masks_stack (ndarray,=): bool, (nspec, nexp)
            Mask of stacked spectrum for this iteration interpolated on the native wavelength grid of the fluxes exposures.
        sn_clip (float): default=30.0,
            Errors are capped in output rejivars so that the S/N is never greater than sn_clip. This prevents overly
            aggressive rejection in high S/N ratio spectra which neverthless differ at a level greater than the implied S/N due to
            systematics.
        title (str):
            Title for QA plot
        debug (bool): default=False
            Show QA plots useful for debuggin.

    Returns:
        rejivars, sigma_corrs, outchi, maskchi

        rejivars: ndarray, (nspec, nexp)
             Updated inverse variances to be used in rejection
        sigma_corrs, ndarray, (nexp)
             Array of correction factors applied to the original ivars to get the new rejivars
        outchi: ndarray, (nspec, nexp)
             The original chi=(fluxes-fluxes_stack)*np.sqrt(ivars) used to determine the correction factors. This
             quantity is useful for plotting. Note that the outchi is computed using the original non-corrected errors.
        maskchi: ndarray, bool, (nspec, nexp)
             Mask returned by renormalize_erorrs indicating which pixels were used in the computation of the correction
             factors. This is basically the union of the input masks but with chi > clip (clip=6.0 is the default)
             values clipped out.
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
        ivar_clip = mask_tot*utils.clip_ivar(thisflux_stack, ivar_tot, sn_clip, mask=mask_tot)
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
                     maxrej=None, maxiter_reject=5, title='', debug=False):
    '''
    Routine for executing the iterative combine and rejection of a set of spectra to compute a final stacked spectrum.

    Args:
        wave_grid: ndarray, (ngrid +1,)
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the new_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.
        waves: ndarray, (nspec, nexp)
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
        fluxes: ndarray, (nspec, nexp)
            fluxes for each exposure on the waves grid
        ivars: ndarray, (nspec, nexp)
            Inverse variances for each exposure on the waves grid
        masks: ndarray, bool, (nspec, nexp)
            Masks for each exposure on the waves grid. True=Good.
        weights: ndarray, (nspec, nexp)
            Weights to be used for combining your spectra. These are computed using sn_weights
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
                definition of sticky.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        title (str):
             Title for QA plot
        debug: bool, default=False,
            Show QA plots useful for debugging.

    Return:

        wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused

        wave_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
        outmask: ndarray, bool, (nspec, nexp)
             Output mask indicating which pixels are rejected in each exposure of the original input spectra after
             performing all of the iterations of combine/rejection
        nused: ndarray, (ngrid,)
             Numer of exposures which contributed to each pixel in the wave_stack. Note that this is in general
             different from nexp because of masking, but also becuse of the sampling specified by wave_grid. In other
             words, sometimes more spectral pixels in the irregularly gridded input wavelength array waves will land in
             one bin versus another depending on the sampling.

    '''
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

    for iexp in range(nexp):
        # nrej = pixels that are now masked that were previously good
        msgs.info("Rejected {:d} pixels in exposure {:d}/{:d}".format(nrej[iexp], iexp, nexp))

    # Compute the final stack using this outmask
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, outmask, weights)

    # Used only for plotting below
    if debug:
        flux_stack_nat, ivar_stack_nat, mask_stack_nat = interp_spec(waves, wave_stack, flux_stack, ivar_stack, mask_stack)
        for iexp in range(nexp):
            # plot the residual distribution for each exposure
            title_renorm = title + ': Error distriution about stack for exposure {:d}/{:d}'.format(iexp,nexp)
            renormalize_errors_qa(outchi[:, iexp], maskchi[:, iexp], sigma_corrs[iexp], title=title_renorm)
            # plot the rejections for each exposures
            title_coadd_iexp = title + ': nrej={:d} pixels rejected,'.format(nrej[iexp]) + \
                               ' norig={:d} originally masked,'.format(norig[iexp]) + \
                               ' for exposure {:d}/{:d}'.format(iexp,nexp)
            coadd_iexp_qa(waves[:, iexp], fluxes[:, iexp], masks[:, iexp], flux_stack_nat[:, iexp], mask_stack_nat[:, iexp],
                          rejivars[:, iexp], outmask[:, iexp], qafile=None, title=title_coadd_iexp)
        # weights qa
        title_weights = title + ': Weights Used -- nrej={:d} total pixels rejected,'.format(np.sum(nrej)) + \
                        ' norig={:d} originally masked'.format(np.sum(norig))
        weights_qa(waves, weights, outmask, title=title_weights)

    return wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused


def scale_spec_stack(wave_grid, waves, fluxes, ivars, masks, weights, ref_percentile=30.0, maxiter_scale=5, sigrej_scale=3,
                     scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5, debug=False, show=False):

    '''
    Routine for optimally combining long or multi-slit spectra or echelle spectra of individual orders. It will
    compute a stacked spectrum from a set of exposures on the specified wave_grid with proper treatment of
    weights and masking. This code calls the stacking code compute_stack, which uses np.histogram to combine the data using
    NGP and does not perform any interpolations and thus does not correlate errors. It uses wave_grid to determine the set
    of wavelength bins that the data are averaged on. The final spectrum will be on an ouptut wavelength grid which is not
    the same as wave_grid. The ouput wavelength grid is the weighted average of the individual wavelengths used for each
    exposure that fell into a given wavelength bin in the input wave_grid. This 1d coadding routine thus maintains the
    independence of the errors for each pixel in the combined spectrum and computes the weighted averaged wavelengths of
    each pixel in an analogous way to the 2d extraction procedure which also never interpolates to avoid correlating
    erorrs. It performs a number of iterations where it combines the spectra and performs rejection of outlier pixels
    using the spec_reject_comb code. The outliers are rejected using the true noise of the individual exposures, but
    uses the distribution of the pixel values about the stack to apply correction factors to the errors before rejecting.
    These corrected errors are currently only used in rejection but are not applied to the data.  This code is based
    on the xidl long_combpsec.pro routine but with significant improvements.

    Args:
        wave_grid: ndarray, (ngrid +1,)
            new wavelength grid desired. This will typically be a reguarly spaced grid created by the new_wave_grid routine.
            The reason for the ngrid+1 is that this is the general way to specify a set of  bins if you desire ngrid
            bin centers, i.e. the output stacked spectra have ngrid elements.  The spacing of this grid can be regular in
            lambda (better for multislit) or log lambda (better for echelle). This new wavelength grid should be designed
            with the sampling of the data in mind. For example, the code will work fine if you choose the sampling to be
            too fine, but then the number of exposures contributing to any given wavelength bin will be one or zero in the
            limiting case of very small wavelength bins. For larger wavelength bins, the number of exposures contributing
            to a given bin will be larger.
        waves: ndarray, (nspec, nexp)
            wavelength arrays for spectra to be stacked. Note that the wavelength grids can in general be different for
            each exposure and irregularly spaced.
        fluxes: ndarray, (nspec, nexp)
            fluxes for each exposure on the waves grid
        ivars: ndarray, (nspec, nexp)
            Inverse variances for each exposure on the waves grid
        masks: ndarray, bool, (nspec, nexp)
            Masks for each exposure on the waves grid. True=Good.
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
                definition of sticky.
        sigrej_scale: float, default=3.0
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        ref_percentile: float, default=20.0
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
        maxiter_scale: int, default=5
            Maximum number of iterations performed for rescaling spectra.
        scale_method: scale method, str, default=None. Options are poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default is to let the
            code determine this automitically which works well.
        dv_smooth: float, 10000.0
            Velocity smoothing used for determining smoothly varying S/N ratio weights by sn_weights
        hand_scale:
            array of hand scale factors, not well tested
        sn_max_medscale (float): default=2.0
            maximum SNR for perforing median scaling
        sn_min_medscale (float): default=0.5
            minimum SNR for perforing median scaling
        debug_scale (bool): default=False
            show interactive QA plots for the rescaling of the spectra
        title (str):
            Title prefix for spec_reject_comb QA plots
        debug (bool): default=False
            show interactive QA plot

    Returns:
        wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused, weights, scales, rms_sn

        wave_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
        outmask: ndarray, bool, (nspec, nexp)
             Output mask indicating which pixels are rejected in each exposure of the original input spectra after
             performing all of the iterations of combine/rejection
        nused: ndarray, (ngrid,)
             Numer of exposures which contributed to each pixel in the wave_stack. Note that this is in general
             different from nexp because of masking, but also becuse of the sampling specified by wave_grid. In other
             words, sometimes more spectral pixels in the irregularly gridded input wavelength array waves will land in
             one bin versus another depending on the sampling.
        weights: ndarray, (nspec, nexp)
            Weights used for combining your spectra which are computed using sn_weights
        scales: ndarray, (nspec, nexp)
            Scale factors applied to each individual spectrum before the combine computed by scale_spec
        rms_sn: ndarray, (nexp,)
            Root mean square S/N ratio of each of your individual exposures computed by sn_weights
    '''

    # Compute an initial stack as the reference, this has its own wave grid based on the weighted averages
    wave_stack, flux_stack, ivar_stack, mask_stack, nused = compute_stack(wave_grid, waves, fluxes, ivars, masks, weights)

    # Rescale spectra to line up with our preliminary stack so that we can sensibly reject outliers
    nexp = np.shape(fluxes)[1]
    fluxes_scale = np.zeros_like(fluxes)
    ivars_scale = np.zeros_like(ivars)
    scales = np.zeros_like(fluxes)
    scale_method_used = []
    for iexp in range(nexp):
        # TODO Create a parset for the coadd parameters!!!
        fluxes_scale[:, iexp], ivars_scale[:, iexp], scales[:, iexp], scale_method_iexp = scale_spec(
            waves[:, iexp],fluxes[:, iexp],ivars[:, iexp], wave_stack, flux_stack, ivar_stack,
            mask=masks[:, iexp], mask_ref=mask_stack, ref_percentile=ref_percentile, maxiters=maxiter_scale,
            sigrej=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale, sn_max_medscale=sn_max_medscale,
            sn_min_medscale=sn_min_medscale, debug=debug, show=show)
        scale_method_used.append(scale_method_iexp)

    return fluxes_scale, ivars_scale, scales, scale_method_used


#Todo: This should probaby take a parset?
#Todo: Make this works for multiple objects after the coadd script input file format is fixed.
def combspec(fnames, objids, ex_value='OPT', flux_value=True, wave_method='pixel', A_pix=None, v_pix=None,
             samp_fact=1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
             sigrej_scale=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
             dv_smooth=1000.0, const_weights=False, maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0,
             maxrej=None, phot_scale_dicts=None, nmaskedge=2,
             qafile=None, outfile = None, debug=False, debug_scale=False, show_scale=False, show=False):

    '''
    Driver routine for coadding longslit/multi-slit spectra. Calls combspec which is the main stacking algorithm.

    Args:
        fnames: list
           a list of spec1d fits file names
        objids: list
           objids you want to combine, i.e the extension name (e.g. 'SPAT0764-SLIT0000-DET07') of
           that spectrum in the spec1d fits files
        ex_value: str, default = 'OPT' for optimal extraction, 'BOX' for boxcar extraction.
        flux_value: bool, default=True
           if True coadd fluxed spectrum, if False coadd spectra in counts
        wave_method: str, default=pixel
           method for generating new wavelength grid with new_wave_grid. Deafult is 'pixel' which creates a uniformly
           space grid in lambda
        A_pix: float,
           dispersion in units of A in case you want to specify it for new_wave_grid, otherwise the code computes the
           median spacing from the data.
        v_pix: float,
           Dispersion in units of km/s in case you want to specify it in the new_wave_grid  (for the 'velocity' option),
           otherwise a median value is computed from the data.
        samp_fact: float, default=1.0
           sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
           samp_fact < 1.0 undersamples (coarser).
        wave_grid_min: float, default=None
           In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
        wave_grid_max: float, default=None
           In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        ref_percentile:
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
        maxiter_scale: int, default=5
            Maximum number of iterations performed for rescaling spectra.
        sigrej_scale: flaot, default=3.0
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
        scale_method: scale method, str, default=None.
            Options are poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default is to let the
            code determine this automitically which works well.
        hand_scale: ndarray,
            Array of hand scale factors, not well tested
        sn_max_medscale: float, default = 2.0,
            maximum SNR for perforing median scaling
        sn_min_medscale: float, default = 0.5
            minimum SNR for perforing median scaling
        dv_smooth: float, 10000.0
            Velocity smoothing used for determining smoothly varying S/N ratio weights by sn_weights
        const_weights: ndarray, (nexp,)
             Constant weight factors specif
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        phot_scale_dicts: dict,
            Dictionary for rescaling spectra to match photometry. Not yet implemented.
        nmaskedge: int, default=2
            Number of edge pixels to mask. This should be removed/fixed.
        qafile: str, default=None
            Root name for QA, if None, it will be determined from the outfile
        outfile: str, default=None,
            Root name for QA, if None, it will come from the target name from the fits header.
        debug: bool, default=False,
            Show all QA plots useful for debugging. Note there are lots of QA plots, so only set this to True if you want to inspect them all.
        debug_scale (bool): default=False
            show interactive QA plots for the rescaling of the spectra
        show: bool, default=False,
             Show key QA plots or not

        Returns:
            wave_stack, flux_stack, ivar_stack, mask_stack

        wave_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
    '''
    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=objids, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value, nmaskedge=nmaskedge)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves, fluxes, ivars, masks, dv_smooth=dv_smooth, const_weights=const_weights, verbose=True)

    fluxes_scale, ivars_scale, scales, scale_method_used = scale_spec_stack(
        wave_grid, waves, fluxes, ivars, masks, weights, ref_percentile=ref_percentile, maxiter_scale=maxiter_scale,
        sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
        sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, debug=debug_scale, show=show_scale)

    # Coadd
    #wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused, weights, scales, rms_sn = \
    #    combspec(wave_grid, waves, fluxes, ivars, masks, ref_percentile=ref_percentile,
    #             maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
    #             sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
    #             const_weights=const_weights, maxiter_reject=maxiter_reject, sn_clip=sn_clip, lower=lower,
    #             upper=upper, maxrej=maxrej, debug_scale=debug_scale, debug=debug, title='multi_combpsec')

    # Rejecting and coadding
    wave_stack, flux_stack, ivar_stack, mask_stack, outmask, nused = spec_reject_comb(
        wave_grid, waves, fluxes_scale, ivars_scale, masks, weights, sn_clip=sn_clip, lower=lower, upper=upper,
        maxrej=maxrej, maxiter_reject=maxiter_reject, debug=debug, title='multi_combspec')

    if show:
        coadd_qa(wave_stack, flux_stack, ivar_stack, nused, mask=mask_stack, title='Stacked spectrum', qafile=qafile)

    # Write to disk?
    if outfile is not None:
        save.save_coadd1d_to_fits(outfile, wave_stack, flux_stack, ivar_stack, mask_stack, header=header,
                                  ex_value=ex_value, overwrite=True)


    return wave_stack, flux_stack, ivar_stack, mask_stack

def ech_combspec(fnames, objids, sensfile=None, ex_value='OPT', flux_value=True, wave_method='loggrid', A_pix=None, v_pix=None,
                 samp_fact=1.0, wave_grid_min=None, wave_grid_max=None, ref_percentile=20.0, maxiter_scale=5,
                 niter_order_scale=3, sigrej_scale=3, scale_method=None, hand_scale=None, sn_max_medscale=2.0, sn_min_medscale=0.5,
                 dv_smooth=10000.0, const_weights=False, maxiter_reject=5, sn_clip=30.0, lower=3.0, upper=3.0,
                 maxrej=None, max_factor=10.0, maxiters=5, min_good=0.05, phot_scale_dicts=None, nmaskedge=2,
                 qafile=None, outfile = None, order_scale=False,
                 merge_stack=False, debug_scale=False, debug=False, show_order_scale=False, show=False):
    '''
    Driver routine for coadding Echelle spectra. Calls combspec which is the main stacking algorithm. It will deliver
    three fits files: spec1d_order_XX.fits (stacked individual orders, one order per extension), spec1d_merge_XX.fits
    (straight combine of stacked individual orders), spec1d_stack_XX.fits (a giant stack of all exposures and all orders).
    In most cases, you should use spec1d_stack_XX.fits for your scientific analyses since it reject most outliers.

    Args:
        fnames: list
           a list of spec1d fits file names
        objids: list
           objids (e.g. 'OBJ0001') you want to combine of that spectrum in the spec1d fits files
        sensfile: str, default = None for a smoothed ivar weighting when sticking different orders
        ex_value: str, default = 'OPT' for optimal extraction, 'BOX' for boxcar extraction.
        flux_value: bool, default=True
           if True coadd fluxed spectrum, if False coadd spectra in counts
        wave_method: str, default=pixel
           method for generating new wavelength grid with new_wave_grid. Deafult is 'pixel' which creates a uniformly
           space grid in lambda
        A_pix: float,
           dispersion in units of A in case you want to specify it for new_wave_grid, otherwise the code computes the
           median spacing from the data.
        v_pix: float,
           Dispersion in units of km/s in case you want to specify it in the new_wave_grid  (for the 'velocity' option),
           otherwise a median value is computed from the data.
        samp_fact: float, default=1.0
           sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
           samp_fact < 1.0 undersamples (coarser).
        wave_grid_min: float, default=None
           In case you want to specify the minimum wavelength in your wavelength grid, default=None computes from data.
        wave_grid_max: float, default=None
           In case you want to specify the maximum wavelength in your wavelength grid, default=None computes from data.
        ref_percentile:
            percentile fraction cut used for selecting minimum SNR cut for robust_median_ratio
        maxiter_scale: int, default=5
            Maximum number of iterations performed for rescaling spectra.
        max_median_factor: float, default=10.0
            maximum scale factor for median rescaling for robust_median_ratio if median rescaling is the method used.
        sigrej_scale: flaot, default=3.0
            Rejection threshold used for rejecting pixels when rescaling spectra with scale_spec.
        scale_method: scale method, str, default=None.
            Options are poly, median, none, or hand. Hand is not well tested.
            User can optionally specify the rescaling method. Default is to let the
            code determine this automitically which works well.
        hand_scale: ndarray,
            Array of hand scale factors, not well tested
        sn_max_medscale: float, default = 2.0,
            maximum SNR for perforing median scaling
        sn_min_medscale: float, default = 0.5
            minimum SNR for perforing median scaling
        dv_smooth: float, 10000.0
            Velocity smoothing used for determining smoothly varying S/N ratio weights by sn_weights
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        const_weights: ndarray, (nexp,)
             Constant weight factors specif
        maxiter_reject: int, default=5
            maximum number of iterations for stacking and rejection. The code stops iterating either when
            the output mask does not change betweeen successive iterations or when maxiter_reject is reached.
        sn_clip: float, default=30.0,
            Errors are capped during rejection so that the S/N is never greater than sn_clip. This prevents overly aggressive rejection
            in high S/N ratio spectrum which neverthless differ at a level greater than the implied S/N due to
            systematics.
        lower: float, default=3.0,
            lower rejection threshold for djs_reject
        upper: float: default=3.0,
            upper rejection threshold for djs_reject
        maxrej: int, default=None,
            maximum number of pixels to reject in each iteration for djs_reject.
        max_factor: float, default = 10.0,
            Maximum allowed value of the returned ratio
        maxiters: int, defrault = 5,
            Maximum number of iterations for astropy.stats.SigmaClip
        min_good: float, default = 0.05
            Minimum fraction of good pixels determined as a fraction of the total pixels for estimating the median ratio
        phot_scale_dicts: dict,
            Dictionary for rescaling spectra to match photometry. Not yet implemented.
        nmaskedge: int, default=2
            Number of edge pixels to mask. This should be removed/fixed.
        qafile: str, default=None
            Root name for QA, if None, it will be determined either the outfile
        outfile: str, default=None,
            Root name for QA, if None, it will come from the target name from the fits header.
        order_scale: bool, default=False,
            Re-scale the orders to match up in the overlap regions. This is currently producing weird results for IR spectra
        merge_stack: bool, default=False,
            Compute an experimental combine of the high S/N combined orders in addition to the default algorithm,
            which is to compute one giant stack using all order overlaps

        debug: bool, default=False,
            Show all QA plots useful for debugging. Note there are lots of QA plots, so only set this to True if you want to inspect them all.
        debug_scale (bool): default=False
            Show interactive QA plots for the rescaling of the spectra for each individua order
        debug_order_scale (bool): default=False
            Show interactive QA plots for the rescaling of the spectra so that the overlap regions match from order to order
        show: bool, default=False,
             Show key QA plots or not

    Returns:
        wave_giant_stack: ndarray, (ngrid,)
             Wavelength grid for stacked spectrum. As discussed above, this is the weighted average of the wavelengths
             of each spectrum that contriuted to a bin in the input wave_grid wavelength grid. It thus has ngrid
             elements, whereas wave_grid has ngrid+1 elements to specify the ngrid total number of bins. Note that
             wave_giant_stack is NOT simply the wave_grid bin centers, since it computes the weighted average.
        flux_giant_stack: ndarray, (ngrid,)
             Final stacked spectrum on wave_stack wavelength grid
        ivar_giant_stack: ndarray, (ngrid,)
             Inverse variance spectrum on wave_stack wavelength grid. Erors are propagated according to weighting and
             masking.
        mask_giant_stack: ndarray, bool, (ngrid,)
             Mask for stacked spectrum on wave_stack wavelength grid. True=Good.
    '''

    # Loading Echelle data
    waves, fluxes, ivars, masks, header = load.load_1dspec_to_array(fnames, gdobj=objids, order=None, ex_value=ex_value,
                                                                    flux_value=flux_value, nmaskedge=nmaskedge)
    # data shape
    nspec, norder, nexp = waves.shape

    # create some arrays
    scales = np.zeros_like(waves)
    weights = np.zeros_like(waves)
    outmasks = np.zeros_like(waves,dtype=bool)

    # output name root for fits and QA plots
    if outfile is None:
        outfile = header['TARGET']+'.fits'
    elif len(outfile.split('.'))==1:
        outfile = outfile+'.fits'

    outfile_order = 'spec1d_order_{:}'.format(outfile)
    outfile_stack = 'spec1d_stack_{:}'.format(outfile)

    if qafile is None:
        qafile = outfile.split('.')[0]+'.pdf'
    qafile_stack = 'spec1d_stack_{:}'.format(qafile)
    qafile_chi = 'spec1d_chi_{:}'.format(qafile)

    # Generate a giant wave_grid
    wave_grid = new_wave_grid(waves, wave_method=wave_method, wave_grid_min=wave_grid_min, wave_grid_max=wave_grid_max,
                              A_pix=A_pix, v_pix=v_pix, samp_fact=samp_fact)

    # Evaluate the sn_weights. This is done once at the beginning
    rms_sn, weights = sn_weights(waves, fluxes, ivars, masks, dv_smooth=dv_smooth, const_weights=const_weights,
                                 verbose=True)

    fluxes_scl_interord = np.zeros_like(fluxes)
    ivars_scl_interord = np.zeros_like(ivars)
    scales_interord = np.zeros_like(fluxes)
    # First perform inter-order scaling once
    for iord in range(norder):
        fluxes_scl_interord[:, iord], ivars_scl_interord[:,iord], scales_interord[:,iord], scale_method_used = \
            scale_spec_stack(wave_grid, waves[:, iord, :], fluxes[:, iord, :], ivars[:, iord, :], masks[:, iord, :],
                             weights[:, iord, :], ref_percentile=ref_percentile,
                             maxiter_scale=maxiter_scale,sigrej_scale=sigrej_scale, scale_method=scale_method,
                             hand_scale=hand_scale,
                             sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, debug=debug_scale)

    # Arrays to store rescaled spectra. Need Fortran like order reshaping to create a (nspec, norder*nexp) stack of spectra
    shape_2d = (nspec, norder * nexp)
    waves_2d = np.reshape(waves, shape_2d, order='F')
    fluxes_2d = np.reshape(fluxes_scl_interord, shape_2d, order='F')
    ivars_2d = np.reshape(ivars_scl_interord, shape_2d, order='F')
    masks_2d = np.reshape(masks, shape_2d, order='F')
    scales_2d = np.reshape(scales_interord, shape_2d, order='F')
    weights_2d = np.reshape(weights, shape_2d, order='F')
    # Iteratively scale and stack the spectra
    fluxes_pre_scale = fluxes_2d.copy()
    ivars_pre_scale = ivars_2d.copy()
    for iter in range(niter_order_scale):
        fluxes_scale, ivars_scale, scales_iter, scale_method_used = scale_spec_stack(
            wave_grid, waves_2d, fluxes_pre_scale, ivars_pre_scale, masks_2d, weights_2d, ref_percentile=ref_percentile,
            maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
            sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale,
            show=(show_order_scale & (iter == (niter_order_scale-1))))
        scales_2d *= scales_iter
        fluxes_pre_scale = fluxes_scale.copy()
        ivars_pre_scale = ivars_scale.copy()

    embed()
    # Arrays to store stacked individual order spectra.
    waves_stack_orders = np.zeros((np.size(wave_grid)-1, norder))
    fluxes_stack_orders = np.zeros_like(waves_stack_orders)
    ivars_stack_orders = np.zeros_like(waves_stack_orders)
    masks_stack_orders = np.zeros_like(waves_stack_orders,dtype=bool)

    # Loop over orders to get the initial stacks of each individual order
    for iord in range(norder):
        # Get the stacked spectrum for each order
        waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], masks_stack_orders[:, iord], \
        outmasks[:,iord,:], nused_iord, weights[:,iord,:], scales[:,iord,:], rms_sn_iord = scale_spec_stack(
            wave_grid, waves[:,iord,:], fluxes[:,iord,:], ivars[:,iord,:], masks[:,iord,:], ref_percentile=ref_percentile,
            maxiter_scale=maxiter_scale, sigrej_scale=sigrej_scale, scale_method=scale_method, hand_scale=hand_scale,
            sn_max_medscale=sn_max_medscale, sn_min_medscale=sn_min_medscale, dv_smooth=dv_smooth,
            const_weights=const_weights, maxiter_reject=maxiter_reject, sn_clip=sn_clip, lower=lower,
            upper=upper, maxrej=maxrej, debug_scale=debug_scale, title='Order-by-Order Combine', debug=debug)

        if show:
            # TODO can we make this bit below more modular for the telluric?
            if sensfile is not None:
                tell_iord = get_tell_from_file(sensfile, waves_stack_orders[:, iord], masks_stack_orders[:, iord], iord=iord)
            else:
                tell_iord = None
            coadd_qa(waves_stack_orders[:, iord], fluxes_stack_orders[:, iord], ivars_stack_orders[:, iord], nused_iord,
                     mask=masks_stack_orders[:, iord], tell=tell_iord,
                     title='Coadded spectrum of order {:}'.format(iord+1))

    # TODO Is this order rescaling currently taking place??
    # Now that we have high S/N ratio individual order stacks, let's compute re-scaling fractors from the order
    # overlaps. We will work from red to blue.
    if order_scale:
        fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = order_median_scale(
            waves_stack_orders, fluxes_stack_orders, ivars_stack_orders, masks_stack_orders,
            min_good=min_good, maxiters=maxiters, max_factor=max_factor, sigrej=sigrej_scale,
            debug=debug_order_scale, show=show)
    else:
        fluxes_stack_orders_scale, ivars_stack_orders_scale, order_ratios = fluxes_stack_orders, ivars_stack_orders, np.ones(norder)

    # apply order_ratios to the scales array: order_ratio*scale
    scales_new = np.zeros_like(scales)
    for iord in range(norder):
        scales_new[:,iord,:] = order_ratios[iord]*scales[:,iord,:]
    fluxes_scale = fluxes * scales_new
    ivars_scale = ivars/scales_new**2


    # Get the new ech_weights for the stack which will merge all the orders
    if sensfile is None:
        rms_sn_stack, order_weights = sn_weights(waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale,
                                                 masks_stack_orders, dv_smooth=dv_smooth, const_weights=const_weights,
                                                 ivar_weights=True, verbose=True)
    else:
        rms_sn_stack = None
        order_weights, masks_stack_orders = sensfunc_weights(sensfile, waves_stack_orders, masks_stack_orders, debug=debug)

    #TODO think through whether this is the correct approach of multiplying weights?
    # apply the sensfunc weights to the orginal weights: sensfunc_weights*weightsf
    ech_weights = np.zeros_like(weights)
    for iord in range(norder):
        mask_weight_iord = masks_stack_orders[:, iord] & (order_weights[:, iord] > 0.0) & (waves_stack_orders[:, iord] > 1.0)
        # Interpolate these order_weights onto the native wavelength grid of each exposure for this order
        for iexp in range(nexp):
            order_weight_interp = scipy.interpolate.interp1d(
                waves_stack_orders[mask_weight_iord, iord], order_weights[mask_weight_iord, iord],  kind = 'cubic',
                bounds_error = False, fill_value = np.nan)(waves[:,iord,iexp])
            ech_weights[:,iord,iexp] = weights[:,iord,iexp] * order_weight_interp


    # TODO Will we use this reject/stack below? It is the straight combine of the stacked individual orders.
    #  This does not take advantage
    #  of the fact that we have many samples in the order overlap regions allowing us to better reject. It does
    #  however have the advatnage that it operates on higher S/N ratio stacked spectra.
    #  should we compute the stack directly with compute_stack or do more rejections with spec_reject_comb?
    #  spec_reject_comb will reject tons of pixels for overlap in telluric region.
    if merge_stack:
        ## Stack with the first method: combine the stacked individual order spectra directly
        wave_merge, flux_merge, ivar_merge, mask_merge, nused = compute_stack(
            wave_grid, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
            order_weights)
        if debug or show:
            qafile_merge = 'spec1d_merge_{:}'.format(qafile)
            coadd_qa(wave_merge, flux_merge, ivar_merge, nused, mask=mask_merge, tell = None,
                     title='Straight combined spectrum of the stacked individual orders', qafile=qafile_merge)


    #TODO Add a note here clarifyng how these reshaped spectra are arranged, i.e. are they packed by the order or by
    # by exposure.

    # reshaping 3D arrays (nspec, norder, nexp) to 2D arrays (nspec, norder*nexp)
    # need Fortran like order reshaping to make sure you are getting the right spectrum for each exposure
    waves_2d = np.reshape(waves,(nspec, norder*nexp),order='F')
    fluxes_2d = np.reshape(fluxes_scale, np.shape(waves_2d),order='F')
    ivars_2d = np.reshape(ivars_scale, np.shape(waves_2d),order='F')
    masks_2d = np.reshape(masks, np.shape(waves_2d),order='F')
    outmasks_2d = np.reshape(outmasks, np.shape(waves_2d),order='F')
    ech_weights_2d = np.reshape(ech_weights, np.shape(waves_2d),order='F')

    wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack, outmask_giant_stack, nused_giant_stack = \
        spec_reject_comb(wave_grid, waves_2d, fluxes_2d, ivars_2d, outmasks_2d, ech_weights_2d, sn_clip=sn_clip,
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
        if debug:
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
                # plot coadd_qa
                coadd_iexp_qa(waves_2d_exps[:,iexp], fluxes_2d_exps[:,iexp], masks_2d_exps[:,iexp],
                              flux_stack_2d_exps[:,iexp], mask_stack_2d_exps[:,iexp],
                              rejivars_2d_exps[:,iexp], outmasks_2d_exps[:,iexp], norder=norder, qafile=None,
                              title=title_coadd_iexp)
        # Global QA
        rejivars_1d, sigma_corrs_1d, outchi_1d, maskchi_1d = update_errors(
            fluxes_2d_exps.flatten(), ivars_2d_exps.flatten(), outmasks_2d_exps.flatten(),
            flux_stack_2d_exps.flatten(), ivar_stack_2d_exps.flatten(), mask_stack_2d_exps.flatten(), sn_clip=sn_clip)
        renormalize_errors_qa(outchi_1d, maskchi_1d, sigma_corrs_1d[0], qafile=qafile_chi, title='Global Chi distribution')
        # show the final coadded spectrum
        coadd_qa(wave_giant_stack, flux_giant_stack, ivar_giant_stack, nused_giant_stack, mask=mask_giant_stack,
                 title='Final stacked spectrum', qafile=qafile_stack)

    # Save stacked individual order spectra
    save.save_coadd1d_to_fits(outfile_order, waves_stack_orders, fluxes_stack_orders_scale, ivars_stack_orders_scale, masks_stack_orders,
                              header=header, ex_value = ex_value, overwrite=True)
    save.save_coadd1d_to_fits(outfile_stack, wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack,
                              header=header, ex_value=ex_value, overwrite=True)
    if merge_stack:
        outfile_merge = 'spec1d_merge_{:}'.format(outfile)
        save.save_coadd1d_to_fits(outfile_merge, wave_merge, flux_merge, ivar_merge, mask_merge, header=header,
                                  ex_value=ex_value, overwrite=True)

    return wave_giant_stack, flux_giant_stack, ivar_giant_stack, mask_giant_stack


