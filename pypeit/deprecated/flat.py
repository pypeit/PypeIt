""" Core module for methods related to flat fielding
"""
import inspect

import numpy as np
import os

from scipy import interpolate
from pypeit import msgs
from pypeit.core import parse
from pypeit.core import pixels
from pypeit.core import tracewave
from scipy.interpolate import interp1d

from pypeit import debugger
from pypeit import utils
from pypeit.core import pydl
from matplotlib import pyplot as plt
import copy
from IPython import embed

import scipy

# NOTE: This is replaced by pypeit.core.flat.tweak_slit_edges
def tweak_slit_edges(slit_left_in, slit_righ_in, ximg_fit, normimg, tweak_slits_thresh,
                     tweak_slits_maxfrac):
    """
    DOC THIS!
    """


    tweak_left = False
    tweak_righ = False
    slit_left_out = np.copy(slit_left_in)
    slit_righ_out = np.copy(slit_righ_in)
    # How many pixels wide is the slit at each Y?
    slitwidth = np.median(slit_righ_in - slit_left_in)
    # Determine the maximum at the left and right end of the slit
    ileft = (ximg_fit > 0.1) & (ximg_fit < 0.4)
    irigh = (ximg_fit > 0.6) & (ximg_fit < 0.9)
#    if (not np.any(ileft)) or (not np.any(irigh)):
#        msgs.error('Cannot tweak slits because much of the slit is masked. You probably have a bad slit')
#        tweak_dict = {'xleft': 0.0, 'xrigh': 0.0,
#                      'norm_max_left': 0.0, 'norm_max_righ': 0.0,
#                      'tweak_left': tweak_left, 'tweak_righ': tweak_righ}
#        return slit_left_out, slit_righ_out, tweak_dict

    #xleft = ximg_fit[ileft]
    #xrigh = ximg_fit[irigh]
    norm_max_left = normimg[ileft].max()
    norm_max_righ = normimg[irigh].max()

    msgs.info('Tweaking slit boundaries using slit illumination function')
    step = 0.001
    # march out from middle to find left edge
    msgs.info('Left threshold = {:5.3f}'.format(tweak_slits_thresh * norm_max_left) +
              ' --  or {:5.3f}'.format(
                  100.0 * tweak_slits_thresh) + ' % of left side max of illumination function = {:5.3f}'.format(norm_max_left))
    for xleft in np.arange(0.5, ximg_fit.min(), -step):
        norm_now = np.interp(xleft, ximg_fit, normimg)
        if (norm_now < tweak_slits_thresh * norm_max_left) & (xleft < tweak_slits_maxfrac):
            slit_left_out += xleft * slitwidth
            tweak_left = True
            msgs.info('Tweaking left slit boundary by {:5.3f}'.format(100 * xleft) +
                      ' %, or {:7.3f}'.format(xleft * slitwidth) + ' pixels')
            if np.abs(xleft - tweak_slits_maxfrac) < 0.01:
                msgs.warn(
                    'Left slit boundary tweak limited by maximum changed allowed by tweak_slits_maxfracn={:5.3f}'.format(
                        100.0 * tweak_slits_maxfrac) + ' %')
            break
    msgs.info('Right threshold = {:5.3f}'.format(tweak_slits_thresh * norm_max_righ) +
              ' --  or {:5.3f}'.format(
                  100.0 * tweak_slits_thresh) + ' % of right side max of illumination function = {:5.3f}'.format(norm_max_righ))
    # march out from middle  to find right edge
    for xrigh in np.arange(0.5, ximg_fit.max(), step):
        norm_now = np.interp(xrigh, ximg_fit, normimg)
        if (norm_now < tweak_slits_thresh * norm_max_righ) & ((1.0 - xrigh) < tweak_slits_maxfrac):
            slit_righ_out -= (1.0 - xrigh) * slitwidth
            tweak_righ = True
            msgs.info('Tweaking right slit boundary by {:5.3f}'.format(100 * (1.0 - xrigh)) +
                      ' %, or {:7.3f}'.format((1.0 - xrigh) * slitwidth) + ' pixels')
            if np.abs((1.0 - xrigh) - tweak_slits_maxfrac) < 0.01:
                msgs.warn(
                    'Right slit boundary tweak limited by maximum changed allowed by tweak_slits_maxfracn={:5.3f}'.format(
                        100.0 * tweak_slits_maxfrac) + ' %')
            break

    tweak_dict = {'xleft': xleft, 'xrigh': xrigh,
                  'norm_max_left': norm_max_left, 'norm_max_righ': norm_max_righ,
                  'tweak_left': tweak_left, 'tweak_righ': tweak_righ}

    return slit_left_out, slit_righ_out, tweak_dict

# NOTE: This was moved to pypeit.flatfield.FlatField.fit and refactored.
def fit_flat(flat, tilts_dict, tslits_dict_in, slit, inmask = None,
             spec_samp_fine = 1.2, spec_samp_coarse = 50.0, spat_samp = 5.0, npoly = None, trim_edg = (3.0,3.0), pad =5.0,
             tweak_slits = True, tweak_slits_thresh = 0.93, tweak_slits_maxfrac = 0.10, nonlinear_counts =1e10, debug=False):


    """ Compute pixelflat and illumination flat from a flat field image.

    Parameters
    ----------
    flat :  float ndarray, shape (nspec, nspat)
        Flat field image in units of electrons.


    tilts_dict: dict
          Dictionary containing wavelength tilts image and other information indicating how wavelengths move across the slit

    tslits_dict: dict
          Dictionary with information on the slit boundaries
    slit: int
          Slit currently being considered
    inmask: boolean ndarray, shape (nspec, nspat), default inmask = None, optional
      Input mask for pixels not to be included in sky subtraction fits. True = Good (not masked), False = Bad (masked)

    spec_samp_fine: float, default = 1.2, optional
      bspline break point spacing in units of pixels for spectral fit to flat field blaze function.

    spec_samp_coarse: float, default = 50.0, optional
      bspline break point spacing in units of pixels for 2-d bspline-polynomial fit to flat field image residuals.
      This should be a large number unless you are trying to fit a sky flat with lots of features.

    spat_samp: float, default = 5.0, optional
      Spatial sampling for spatial slit illumination function. This is the width of the median filter in pixels used to
      determine the slit illumination function, and thus sets the minimum scale on which the illumination function will
      have features.

    trim_edg: tuple of floats  (left_edge, right_edge), default (3,3), optional
      indicates how many pixels to trim from left and right slit edges for creating the edgemask, which is used to mask
      the edges from the initial (fine) spectroscopic fit to the blaze function.

    pad: int, default = 5, optional
      Padding window used to create expanded slitmask images used for re-determining slit boundaries. Tilts are also
      computed using this expanded slitmask in cases the slit boundaries need to be moved outward.

    npoly: int, default = None, optional
      Order of polynomial for 2-d bspline-polynomial fit to flat field image residuals. The code determines the order of
      these polynomials to each slit automatically depending on the slit width, which is why the default is None.
      Do not attempt to set this paramter unless you know what you are doing.


    tweak_slits: bool, default = True, optional
      Slit edges will be tweaked such the left and right bounadaries intersect the location where the illumination
      function falls below tweak_slits_thresh (see below) of its maximum value near the center (moving out from the center)

    tweak_slits_thresh: float, default = 0.93, optional
      If tweak_slits is True, this sets the illumination function threshold used to tweak the slits

    tweak_slits_maxfrac: float, default = 0.10, optional
      Maximum fractinoal amount (of slit width) allowed for each trimming the left and right slit boundaries, i.e. the
      default is 10% which means slits would shrink by at most 20% (10% on each side)

    debug: bool, default = False, optional
      Show plots useful for debugging. This will block further execution of the code until the plot windows are closed.

    Returns
    -------
    pixeflat:   ndarray with same shape as flat
      Pixelflat gives pixel-to-pixel variations of detector response. Values are centered about unity.

    illumflat:  ndarray with same shape as flat
      Illumination flat gives variations of the slit illumination function across the spatial direction of the detect.
      Values are centered about unity. The slit illumination function is computed by dividing out the spectral response and
      collapsing out the spectral direction.

    flat_model:  ndarray with same shape as flat
      Full 2-d model image of the input flat image in units of electrons.  The pixelflat is defined to be flat/flat_model.

    tilts: ndarray with same shape as flat
      Tilts image fit for this slit evaluated using the new slit boundaries

    thismask_out: ndarray with same shape as flat, bool
       Boolean mask indicating which pixels are on the slit now with the new slit boundaries

    slit_left_out: ndarray with shape (nspec,)
       Tweaked left slit bounadries

    slit_righ_out: ndarray with shape (nspec,)
       Tweaked right slit bounadries

    Notes
    -----
    
    Revision History
        - 11-Mar-2005  First version written by Scott Burles.
        - 2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
        - 3-Sep-2018 Ported to python by J. F. Hennawi and significantly improved
    """

    shape = flat.shape
    nspec = shape[0]
    nspat = shape[1]

    # Get the thismask_in and input slit bounadries from the tslits_dict
    slit_left_in = tslits_dict_in['slit_left'][:,slit]
    slit_righ_in = tslits_dict_in['slit_righ'][:,slit]
    thismask_in = pixels.tslits2mask(tslits_dict_in) == slit

    # Check for saturation of the flat. If there are not enough pixels do not attempt a fit
    good_frac = np.sum(thismask_in & (flat < nonlinear_counts))/np.sum(thismask_in)
    if good_frac < 0.5:
        msgs.warn(msgs.newline() + 'Only {:4.2f}'.format(100*good_frac) + '% of the pixels on this slit are not saturated.' +
                  msgs.newline() + 'Consider raising nonlinear_counts={:5.3f}'.format(nonlinear_counts) +
                  msgs.newline() + 'Not attempting to flat field slit# {:d}'.format(slit))
        return np.ones_like(flat), np.ones_like(flat), np.zeros_like(flat), tilts_dict['tilts'], thismask_in, slit_left_in, slit_righ_in

    # Approximate number of pixels sampling each spatial pixel for this (original) slit.
    npercol = np.fmax(np.floor(np.sum(thismask_in)/nspec),1.0)
    # Demand at least 10 pixels per row (on average) per degree of the polynomial
    if npoly is None:
        npoly_in = 7
        npoly  = np.clip(npoly_in, 1, np.ceil(npercol/10.).astype(int))
        #npoly = np.fmax(np.fmin(npoly_in, (np.ceil(npercol/10.)).astype(int)),1)


    ximg_in, edgmask_in = pixels.ximg_and_edgemask(slit_left_in, slit_righ_in, thismask_in, trim_edg=trim_edg)
    # Create a fractional position image ximg that encompasses the whole image, rather than just the thismask_in slit pixels
    spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
    slit_left_img = np.outer(slit_left_in, np.ones(nspat))   # left slit boundary replicated spatially
    slitwidth_img = np.outer(slit_righ_in - slit_left_in, np.ones(nspat)) # slit width replicated spatially
    ximg = (spat_img - slit_left_img)/slitwidth_img

    # Create a wider slitmask image with shift pixels padded on each side
    slitmask_pad = pixels.tslits2mask(tslits_dict_in, pad = pad)
    thismask = (slitmask_pad == slit) # mask enclosing the wider slit bounadries
    # Create a tilts image using this padded thismask, rather than using the original thismask_in slit pixels
    tilts = tracewave.fit2tilts(shape, tilts_dict['coeffs'], tilts_dict['func2d'])
    piximg = tilts * (nspec-1)
    pixvec = np.arange(nspec)

    if inmask is None:
        inmask = np.copy(thismask)

    # Fit the spectral direction of the blaze. We do this in the log
    log_flat = np.log(np.fmax(flat, 1.0))
    inmask_log = ((flat > 1.0) & inmask)
    log_ivar = inmask_log.astype(float)/0.5**2 # set errors to just be 0.5 in the log

    # Flat field pixels for fitting spectral direction. Restrict to original slit pixels
    fit_spec = thismask_in & inmask & np.invert(edgmask_in) #& (flat < nonlinear_counts)
    nfit_spec = np.sum(fit_spec)
    spec_frac = nfit_spec/np.sum(thismask_in)
    msgs.info('Spectral fit of flatfield for {:}'.format(nfit_spec) + ' pixels')
    if spec_frac < 0.5:
        msgs.warn('Spectral flatfield fit is to only {:4.2f}'.format(100*spec_frac) + '% of the pixels on this slit.' +
                  msgs.newline() + '          Something appears to be wrong here')

    isrt_spec = np.argsort(piximg[fit_spec])
    pix_fit = piximg[fit_spec][isrt_spec]
    log_flat_fit = log_flat[fit_spec][isrt_spec]
    log_ivar_fit = log_ivar[fit_spec][isrt_spec]
    inmask_log_fit = inmask_log[fit_spec][isrt_spec]
    logrej = 0.5 # rejectino threshold for spectral fit in log(image)

    # ToDo Figure out how to deal with the fits going crazy at the edges of the chip in spec direction
    spec_set_fine, outmask_spec, specfit, _, exit_status = \
        utils.bspline_profile(pix_fit, log_flat_fit, log_ivar_fit,np.ones_like(pix_fit), inmask = inmask_log_fit,
        nord = 4, upper=logrej, lower=logrej,
        kwargs_bspline = {'bkspace':spec_samp_fine},kwargs_reject={'groupbadpix':True, 'maxrej': 5})


    # Debugging/checking spectral fit
    if debug:
        goodbk = spec_set_fine.mask
        specfit_bkpt, _ = spec_set_fine.value(spec_set_fine.breakpoints[goodbk])
        was_fit_and_masked = (outmask_spec == False)
        plt.clf()
        ax = plt.gca()
        ax.plot(pix_fit,log_flat_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',
                linestyle='None', label = 'all pixels')
        ax.plot(pix_fit[was_fit_and_masked],log_flat_fit[was_fit_and_masked], color='red', marker='+',
                markersize=1.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
        ax.plot(pix_fit, specfit, color='cornflowerblue', label = 'fit to blaze')
        ax.plot(spec_set_fine.breakpoints[goodbk], specfit_bkpt, color='lawngreen', marker='o', markersize=2.0,
                mfc='lawngreen', fillstyle='full', linestyle='None', label='bspline breakpoints')
        ax.set_ylim((0.99*specfit.min(),1.01*specfit.max()))
        plt.legend()
        plt.xlabel('Spectral Pixel')
        plt.ylabel('log(flat counts)')
        plt.title('Spectral Fit for slit={:d}'.format(slit))
        plt.show()
#        # JXP
#        plt.clf()
#        ax = plt.gca()
#        ax.scatter(pix_fit, log_flat_fit)
#        embed(header='265 of flat.py')

    # Evaluate and save
    spec_model = np.ones_like(flat)
    spec_model[thismask], _ = np.exp(spec_set_fine.value(piximg[thismask]))
    norm_spec = np.ones_like(flat)
    norm_spec[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)

    # Flat field pixels for fitting spatial direction
    # Determine maximum counts in median filtered flat spectrum. Only fit pixels > 0.1 of this maximum
    specfit_interp = interp1d(pix_fit, specfit, kind='linear', bounds_error=False, fill_value=-np.inf)
    log_specfit = specfit_interp(pixvec)
    specvec = np.exp(log_specfit)
    spec_sm = utils.fast_running_median(specvec,np.fmax(np.ceil(0.10*nspec).astype(int),10))
    spec_sm_max = spec_sm.max()
    fit_spat = thismask & inmask &  (spec_model > 1.0) & (spec_model > 0.1*spec_sm_max) & \
               (norm_spec > 0.0) & (norm_spec < 1.7)  #& (flat < nonlinear_counts)
    nfit_spat = np.sum(fit_spat)
    spat_frac = nfit_spat/np.sum(thismask)
    msgs.info('Spatial fit to flatfield for {:}'.format(nfit_spec) + ' pixels')
    if spat_frac < 0.5:
        msgs.warn('Spatial flatfield fit is to only {:4.2f}'.format(100*spat_frac) + '% of the pixels on this slit.' +
                  msgs.newline() + '              Something apperas to be wrong here')


    isrt_spat = np.argsort(ximg[fit_spat])
    ximg_fit = ximg[fit_spat][isrt_spat]
    norm_spec_fit = norm_spec[fit_spat][isrt_spat]
    #norm_spec_ivar = np.ones_like(norm_spec_fit)/(spat_illum_thresh**2)
    nfit_spat = np.sum(fit_spat)

    slitwidth = np.median(slit_righ_in - slit_left_in) # How many pixels wide is the slit at each Y?
    ximg_resln = spat_samp/slitwidth

    med_width = (np.ceil(nfit_spat*ximg_resln)).astype(int)
    normimg_raw = utils.fast_running_median(norm_spec_fit,med_width)
    sig_res = np.fmax(med_width/20.0,0.5)
    normimg = scipy.ndimage.filters.gaussian_filter1d(normimg_raw,sig_res, mode='nearest')

    # mask regions where illumination function takes on extreme values
    if np.any(np.invert(np.isfinite(normimg))):
        msgs.error('Inifinities in slit illumination function computation normimg')

    # Determine the breakpoint spacing from the sampling of the ximg
    ximg_samp = np.median(ximg_fit - np.roll(ximg_fit,1))
    ximg_1pix = 1.0/slitwidth
    # Use breakpoints at a spacing of a 1/10th of a pixel, but do not allow a bsp smaller than the typical sampling
    ximg_bsp  = np.fmax(ximg_1pix/10.0, ximg_samp*1.2)
    bsp_set = pydl.bspline(ximg_fit,nord=4, bkspace=ximg_bsp)
    fullbkpt = bsp_set.breakpoints
    spat_set, outmask_spat, spatfit, _, exit_status = \
        utils.bspline_profile(ximg_fit, normimg, np.ones_like(normimg),np.ones_like(normimg),
        nord=4,upper=5.0, lower=5.0,fullbkpt = fullbkpt)

    # Evaluate and save
    illumflat = np.ones_like(flat)
    illumflat[thismask], _ = spat_set.value(ximg[thismask])
    norm_spec_spat = np.ones_like(flat)
    norm_spec_spat[thismask] = flat[thismask]/np.fmax(spec_model[thismask], 1.0)/np.fmax(illumflat[thismask],0.01)

    if tweak_slits:
        slit_left_out, slit_righ_out, tweak_dict = tweak_slit_edges(
            slit_left_in, slit_righ_in, ximg_fit, normimg, tweak_slits_thresh, tweak_slits_maxfrac)
        # Recreate all the quantities we need based on the tweaked slits
        tslits_dict_out = copy.deepcopy(tslits_dict_in)
        tslits_dict_out['slit_left'][:,slit] = slit_left_out
        tslits_dict_out['slit_righ'][:,slit] = slit_righ_out
        slitmask_out = pixels.tslits2mask(tslits_dict_out)
        thismask_out = (slitmask_out == slit)
        ximg_out, edgmask_out = pixels.ximg_and_edgemask(slit_left_out, slit_righ_out, thismask_out, trim_edg=trim_edg)
        # Note that nothing changes with the tilts, since these were already extrapolated across the whole image.
    else:
        # Generate the edgemask using the original slit boundaries and thismask_in
        slit_left_out = np.copy(slit_left_in)
        slit_righ_out = np.copy(slit_righ_in)
        thismask_out = thismask_in
        ximg_out = ximg_in

    # Add an approximate pixel axis at the top
    if debug:
        plt.clf()
        ax = plt.gca()
        ax.plot(ximg_fit, norm_spec_fit, color='k', marker='o', markersize=0.4, mfc='k', fillstyle='full',linestyle='None',
                label = 'all pixels')
        #ax.plot(ximg_fit[~imed], norm_spec_fit[~imed], color='darkred', marker='+',markersize=4.0, mfc='red',
        #        fillstyle='full', linestyle='None', label = 'masked')
        #ax.plot(ximg_fit[imed], normfit[imed], color='orange', label = 'median spatial profile')
        ax.plot(ximg_fit, spatfit, color='cornflowerblue', label = 'final slit illumination function')
        ymin = np.fmax(0.8 * spatfit.min(), 0.5)
        ymax = 1.2*spatfit.max()
        ax.set_ylim((np.fmax(0.8 * spatfit.min(), 0.5), 1.2 * spatfit.max()))
        ax.set_xlim(ximg_fit.min(), ximg_fit.max())
        plt.vlines(0.0, ymin, ymax, color='lightgreen', linestyle=':', linewidth=2.0, label='original left edge',zorder=8)
        plt.vlines(1.0,ymin,ymax, color='red',linestyle=':', linewidth = 2.0, label='original right edge',zorder=9)
        if tweak_slits:
            if tweak_dict['tweak_left']:
                label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) + ' % of max of left illumprofile'

                plt.hlines(tweak_slits_thresh*tweak_dict['norm_max_left'], ximg_fit.min(), 0.5, color='lightgreen',
                           linewidth=3.0,label=label, zorder=10)
                plt.vlines(tweak_dict['xleft'],ymin,ymax, color='lightgreen',linestyle='--', linewidth = 3.0, label='tweaked left edge',zorder=11)
            if tweak_dict['tweak_righ']:
                label = 'threshold = {:5.2f}'.format(tweak_slits_thresh) + ' % of max of right illumprofile'

                plt.hlines(tweak_slits_thresh * tweak_dict['norm_max_righ'], 0.5, ximg_fit.max(), color='red', linewidth=3.0,
                           label=label, zorder=10)
                plt.vlines(tweak_dict['xrigh'],ymin,ymax, color='red',linestyle='--', linewidth = 3.0, label='tweaked right edge',zorder=20)
        plt.legend()
        plt.xlabel('Normalized Slit Position')
        plt.ylabel('Normflat Spatial Profile')
        plt.title('Illumination Function Fit for slit={:d}'.format(slit))
        plt.show()

    msgs.info('Performing illumination + scattembedered light flat field fit')

    # Flat field pixels for fitting spectral direction
    isrt_spec = np.argsort(piximg[thismask_out])
    pix_twod = piximg[thismask_out][isrt_spec]
    ximg_twod = ximg_out[thismask_out][isrt_spec]
    norm_twod = norm_spec_spat[thismask_out][isrt_spec]

    fitmask = inmask[thismask_out][isrt_spec] & (np.abs(norm_twod - 1.0) < 0.30)
    # Here we ignore the formal photon counting errors and simply assume that a typical error per pixel.
    # This guess is somewhat aribtrary. We then set the rejection threshold with sigrej_illum
    var_value = 0.01
    norm_twod_ivar = fitmask.astype(float)/(var_value**2)
    sigrej_illum = 4.0

    poly_basis = pydl.fpoly(2.0*ximg_twod - 1.0, npoly).T

    # Perform the full 2d fit now
    twod_set, outmask_twod, twodfit, _ , exit_status = \
        utils.bspline_profile(pix_twod, norm_twod, norm_twod_ivar,poly_basis,inmask = fitmask, nord = 4,
        upper=sigrej_illum, lower=sigrej_illum,
        kwargs_bspline = {'bkspace':spec_samp_coarse},kwargs_reject={'groupbadpix':True, 'maxrej': 10})

    if debug:
        resid = (norm_twod  - twodfit)
        badpix = np.invert(outmask_twod) & fitmask
        goodpix = outmask_twod & fitmask
        plt.clf()
        ax = plt.gca()
        ax.plot(pix_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',linestyle='None',
                label = 'good points')
        ax.plot(pix_twod[badpix],resid[badpix], color='red', marker='+',markersize=0.5, mfc='red', fillstyle='full', linestyle='None', label='masked')
        plt.hlines(sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
                   label='rejection thresholds',zorder=10,linewidth=2.0)
        plt.hlines(-sigrej_illum*var_value,pix_twod.min(),pix_twod.max(), color='lawngreen',linestyle='--',
                   zorder=10,linewidth=2.0)
        ax.set_ylim((-0.05,0.05))
        ax.set_xlim((pix_twod.min(), pix_twod.max()))
        plt.legend()
        plt.xlabel('Spectral Pixel')
        plt.ylabel('Residuals from pixelflat 2-d fit')
        plt.title('Spectral Residuals for slit={:d}'.format(slit))
        plt.show()

        plt.clf()
        ax = plt.gca()
        ax.plot(ximg_twod[goodpix], resid[goodpix], color='k', marker='o', markersize=0.2, mfc='k', fillstyle='full',
                linestyle='None',
                label='good points')
        ax.plot(ximg_twod[badpix], resid[badpix], color='red', marker='+', markersize=0.5, mfc='red', fillstyle='full',
                linestyle='None', label='masked')
        plt.hlines(sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
                   label='rejection thresholds', zorder=10,linewidth=2.0)
        plt.hlines(-sigrej_illum * var_value, ximg_twod.min(), ximg_twod.max(), color='lawngreen', linestyle='--',
                   zorder=10,linewidth=2.0)
        ax.set_ylim((-0.05, 0.05))
        ax.set_xlim(-0.02, 1.02)
        plt.legend()
        plt.xlabel('Normalized Slit Position')
        plt.ylabel('Residuals from pixelflat 2-d fit')
        plt.title('Spatial Residuals for slit={:d}'.format(slit))
        plt.show()

    # Evaluate and save
    twod_model = np.ones_like(flat)
    twod_this = np.zeros_like(twodfit)
    twod_this[isrt_spec] = twodfit
    twod_model[thismask_out] = twod_this

    # Compute all the final output images output
    pixelflat = np.ones_like(flat)
    flat_model = np.ones_like(flat)
    flat_model[thismask_out] = twod_model[thismask_out]*np.fmax(illumflat[thismask_out],0.05)*np.fmax(spec_model[thismask_out],1.0)
    pixelflat[thismask_out] = flat[thismask_out]/flat_model[thismask_out]

    # ToDo Add some code here to treat the edges and places where fits go bad?
    # Set the pixelflat to 1.0 wherever the flat was nonlinear
    pixelflat[flat >= nonlinear_counts] = 1.0
    # Do not apply pixelflat field corrections that are greater than 100% to avoid creating edge effects, etc.
    # TODO Should we do the same for the illumflat??
    #pixelflat = np.fmax(np.fmin(pixelflat, 2.0), 0.5)
    pixelflat = np.clip(pixelflat, 0.5, 2.0)

    return pixelflat, illumflat, flat_model, tilts, thismask_out, slit_left_out, slit_righ_out



