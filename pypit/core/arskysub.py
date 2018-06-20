""" Module for sky subtraction
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import sys, os
from matplotlib import pyplot as plt

#from pydl.pydlutils.bspline import bspline
from pypit.core import pydl

from pypit import msgs

from pypit import arutils
from pypit import arpixels

from pypit import ardebug as debugger


def bg_subtraction_slit(slit, slitpix, edge_mask, sciframe, varframe, tilts,
                        bpm=None, crmask=None, tracemask=None, bsp=0.6, sigrej=3.):
    """
    Perform sky subtraction on an input slit

    Parameters
    ----------
    slit : int
      Slit number; indexed 1, 2,
    slitpix : ndarray
      Specifies pixels in the slits
    edgemask : ndarray
      Mask edges of the slit
    sciframe : ndarray
      science frame
    varframe : ndarray
      Variance array
    tilts : ndarray
      Tilts of the wavelengths
    bpm : ndarray, optional
      Bad pixel mask
    crmask : ndarray
      Cosmic ray mask
    tracemask : ndarray
      Object mask
    bsp : float
      Break point spacing
    sigrej : float
      rejection

    Returns
    -------
    bgframe : ndarray
      Sky background image

    """

    # Init
    bgframe = np.zeros_like(sciframe)
    ivar = arutils.calc_ivar(varframe)
    ny = sciframe.shape[0]
    piximg = tilts * (ny-1)

    #
    ordpix = slitpix.copy()
    # Masks
    if bpm is not None:
        ordpix *= 1-bpm.astype(np.int)
    if crmask is not None:
        ordpix *= 1-crmask.astype(np.int)
    if tracemask is not None:
        ordpix *= (1-tracemask.astype(np.int))

    # Sky pixels for fitting
    fit_sky = (ordpix == slit) & (ivar > 0.) & (~edge_mask)
    isrt = np.argsort(piximg[fit_sky])
    wsky = piximg[fit_sky][isrt]
    sky = sciframe[fit_sky][isrt]
    sky_ivar = ivar[fit_sky][isrt]

    # All for evaluation
    all_slit = (slitpix == slit) & (~edge_mask)

    # Pre-fit
    pos_sky = (sky > 1.0) & (sky_ivar > 0.)
    if np.sum(pos_sky) > ny:
        lsky = np.log(sky[pos_sky])
        lsky_ivar = lsky * 0. + 0.1

        # Init bspline to get the sky breakpoints (kludgy)
        tmp = pydl.bspline(wsky[pos_sky], nord=4, bkspace=bsp)

        #skybkpt = bspline_bkpts(wsky[pos_sky], nord=4, bkspace=bsp $
        #, / silent)
        if False:
            from matplotlib import pyplot as plt
            plt.clf()
            ax = plt.gca()
            ax.scatter(wsky[pos_sky], lsky)
            #ax.scatter(wsky[~full_out], sky[~full_out], color='red')
            #ax.plot(wsky, yfit, color='green')
            plt.show()
            debugger.set_trace()
        lskyset, outmask, lsky_fit, red_chi = arutils.bspline_profile(
            wsky[pos_sky], lsky, lsky_ivar, np.ones_like(lsky),
            fullbkpt = tmp.breakpoints, upper=sigrej, lower=sigrej,
            kwargs_reject={'groupbadpix':True})
        res = (sky[pos_sky] - np.exp(lsky_fit)) * np.sqrt(sky_ivar[pos_sky])
        lmask = (res < 5.0) & (res > -4.0)
        sky_ivar[pos_sky] = sky_ivar[pos_sky] * lmask

    # Full fit now
    full_bspline = pydl.bspline(wsky, nord=4, bkspace=bsp)
    skyset, full_out, yfit, _ = arutils.bspline_profile(
        wsky, sky, sky_ivar, np.ones_like(sky),
        fullbkpt=full_bspline.breakpoints,
        upper=sigrej, lower=sigrej, kwargs_reject={'groupbadpix':True, 'maxrej': 10})
    # Evaluate and save
    bgframe[all_slit] = skyset.value(piximg[all_slit])[0] #, skyset)

    # Debugging/checking
    if False:
        from matplotlib import pyplot as plt
        plt.clf()
        ax = plt.gca()
        ax.scatter(wsky[full_out], sky[full_out])
        ax.scatter(wsky[~full_out], sky[~full_out], color='red')
        ax.plot(wsky, yfit, color='green')
        plt.show()

    # Return
    return bgframe


# Utility routine used by local_bg_subtraction_slit
def skyoptimal(wave,data,ivar, oprof, sortpix, sigrej = 3.0, npoly = 1, spatial = None, fullbkpt = None):

    nx = data.size
    nc = oprof.shape[0]
    nobj = int(oprof.size / nc)
    if nc != nx:
        raise ValueError('Object profile should have oprof.shape[0] equal to nx')

    msgs.info('Iter     Chi^2     Rejected Pts')
    xmin = 0.0
    xmax = 1.0

    if ((npoly == 1) | (spatial == None)):
        profile_basis = np.column_stack((oprof, np.ones(nx)))
    else:
        xmin = spatial.min()
        xmax = spatial.max()
        x2 = 2.0 * (spatial - xmin) / (xmax - xmin) - 1
        poly_basis = flegendre(x2, npoly).T
        profile_basis = np.column_stack((oprof, poly_basis))

    #ToDo this seems different from IDL
    if nobj == 1:
        relative_mask = (oprof > 1e-10)
    else:
        relative_mask = (np.sum(oprof, axis=1) > 1e-10)

    indx, = np.where(ivar[sortpix] > 0.0)
    ngood = indx.size
    good = sortpix[indx]
    good = good[wave[good].argsort()]
    relative, = np.where(relative_mask[good])

    (sset1, outmask_good1, yfit1, red_chi1) = bspline_longslit(wave[good], data[good], ivar[good], profile_basis[good, :],
                                                               fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                                               relative=relative,
                                                               kwargs_reject={'groupbadpix': True, 'maxrej': 5})

    chi2 = (data[good] - yfit1) ** 2 * ivar[good]
    chi2_srt = np.sort(chi2)
    gauss_prob = 1.0 - 2.0 * ndtr(-1.2 * sigrej)
    sigind = int(np.fmin(np.rint(gauss_prob * float(ngood)), ngood - 1))
    chi2_sigrej = chi2_srt[sigind]
    mask1 = (chi2 < chi2_sigrej)
    msgs.info('2nd round....')
    msgs.info('Iter     Chi^2     Rejected Pts')

    (sset, outmask_good, yfit, red_chi) = bspline_longslit(wave[good], data[good], ivar[good] * mask1,
                                                           profile_basis[good, :],
                                                           fullbkpt=fullbkpt, upper=sigrej, lower=sigrej,
                                                           relative=relative,
                                                           kwargs_reject={'groupbadpix': True, 'maxrej': 1})

    ncoeff = npoly + nobj
    skyset = bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord, npoly=npoly)
    skyset.coeff = sset.coeff[nobj:, :]  # Coefficients for the sky
    skyset.mask = sset.mask
    skyset.xmin = xmin
    skyset.xmax = xmax

    sky_bmodel, _ = skyset.value(wave, x2=spatial)

    obj_bmodel = np.zeros(sky_bmodel.shape)
    objset = bspline(None, fullbkpt=sset.breakpoints, nord=sset.nord)
    objset.mask = sset.mask
    for i in range(nobj):
        objset.coeff = sset.coeff[i, :]
        obj_bmodel1, _ = objset.value(wave)
        obj_bmodel = obj_bmodel + obj_bmodel1 * profile_basis[:, i]

    outmask = np.zeros(wave.shape, dtype=bool)
    outmask[good] = outmask_good

    return (sky_bmodel, obj_bmodel, outmask)



# This code is deprecated and replaced by bg_subtraction_slit
def orig_bg_subtraction_slit(tslits_dict, pixlocn,
                        slit, tilts, sciframe, varframe, bpix, crpix,
                        settings,
                        tracemask=None,
                        rejsigma=3.0, maskval=-999999.9,
                        method='bspline'):
    """ Extract a science target and background flux
    :param slf:
    :param sciframe:
    :param varframe:
    :return:
    """
    # Unpack tslits
    lordloc = tslits_dict['lcen']
    rordloc = tslits_dict['rcen']
    slitpix = tslits_dict['slitpix']
    # Init
    bgframe = np.zeros_like(sciframe)

    # Begin the algorithm
    # Find which pixels are within the order edges
    msgs.info("Identifying pixels within each order")
    ordpix = arpixels.new_order_pixels(pixlocn, lordloc*0.95+rordloc*0.05,
                              lordloc*0.05+rordloc*0.95)

    msgs.info("Applying bad pixel mask")
    ordpix *= (1-bpix.astype(np.int)) * (1-crpix.astype(np.int))
    if tracemask is not None: ordpix *= (1-tracemask.astype(np.int))

    # Construct an array of pixels to be fit with a spline
    whord = np.where(ordpix == slit+1)
    xvpix  = tilts[whord]
    scipix = sciframe[whord]
    xargsrt = np.argsort(xvpix, kind='mergesort')
    sxvpix  = xvpix[xargsrt]
    sscipix = scipix[xargsrt]

    # Reject deviant pixels -- step through every 1.0/sciframe.shape[0] in sxvpix and reject significantly deviant pixels
    edges = np.linspace(min(0.0,np.min(sxvpix)),max(1.0,np.max(sxvpix)),sciframe.shape[0])
    fitcls = np.zeros(sciframe.shape[0])

    # Identify science target
    maskpix = np.zeros(sxvpix.size)
    msgs.info("Identifying pixels containing the science target")
    msgs.work("Speed up this step with multi-processing")
    for i in range(sciframe.shape[0]-1):
        wpix = np.where((sxvpix>=edges[i]) & (sxvpix<=edges[i+1]))
        if (wpix[0].size>5):
            txpix = sxvpix[wpix]
            typix = sscipix[wpix]
            msk, cf = arutils.robust_polyfit(txpix, typix, 0, sigma=rejsigma)
            maskpix[wpix] = msk
            #fitcls[i] = cf[0]
            wgd=np.where(msk == 0)
            szt = np.size(wgd[0])
            if szt > 8:
                fitcls[i] = np.mean(typix[wgd][szt//2-3:szt//2+4]) # Average the 7 middle pixels
                #fitcls[i] = np.mean(np.random.shuffle(typix[wgd])[:5]) # Average the 5 random pixels
            else:
                fitcls[i] = cf[0]
    '''
    else:
        debugger.set_trace()
        msgs.work("Speed up this step in cython")
        for i in range(sciframe.shape[0]-1):
            wpix = np.where((sxvpix >= edges[i]) & (sxvpix <= edges[i+1]))
            typix = sscipix[wpix]
            szt = typix.size
            if szt > 8:
                fitcls[i] = np.mean(typix[szt//2-3:szt//2+4])  # Average the 7 middle pixels
            elif szt != 0:
                fitcls[i] = np.mean(typix)
            else:
                fitcls[i] = 0.0
        # Trace the sky lines to get a better estimate of the tilts
        scicopy = sciframe.copy()
        scicopy[np.where(ordpix == slit)] = maskval
        scitilts, _ = artrace.model_tilt(det, scicopy, guesstilts=tilts.copy(),
                                         censpec=fitcls, maskval=maskval, plotQA=True)
        xvpix  = scitilts[whord]
        scipix = sciframe[whord]
        varpix = varframe[whord]
        mskpix = tracemask[whord]
        xargsrt = np.argsort(xvpix, kind='mergesort')
        sxvpix  = xvpix[xargsrt]
        sscipix = scipix[xargsrt]
        svarpix = varpix[xargsrt]
        maskpix = mskpix[xargsrt]
    '''

    # Check the mask is reasonable
    scimask = sciframe.copy()
    rxargsrt = np.argsort(xargsrt, kind='mergesort')
    scimask[whord] *= (1.0-maskpix)[rxargsrt]

    # Now trace the sky lines to get a better estimate of the spectral tilt during the observations
    scifrcp = scimask.copy()
    scifrcp[whord] += (maskval*maskpix)[rxargsrt]
    scifrcp[np.where(ordpix != slit+1)] = maskval

    # Check tilts? -- Can also be error in flat fielding or slit illumination
    if False:
        idx = 1893
        plt.clf()
        ax = plt.gca()
        ax.scatter(tilts[idx-2,:], scifrcp[idx-2,:], color='green')
        ax.scatter(tilts[idx-1,:], scifrcp[idx-1,:], color='blue')
        ax.scatter(tilts[idx,:], scifrcp[idx,:], color='red')
        ax.scatter(tilts[idx+1,:], scifrcp[idx+1,:], color='orange')
        ax.set_ylim(0., 3000)
        plt.show()
        debugger.set_trace()
    #
    msgs.info("Fitting sky background spectrum")
    if method == 'bspline':
        msgs.info("Using bspline sky subtraction")
        gdp = (scifrcp != maskval) & (ordpix == slit+1) & (varframe > 0.)
        srt = np.argsort(tilts[gdp])
        #bspl = arutils.func_fit(tilts[gdp][srt], scifrcp[gdp][srt], 'bspline', 3,
        #                        **settings.argflag['reduce']['skysub']['bspline'])
        ivar = arutils.calc_ivar(varframe)
        mask, bspl = arutils.robust_polyfit(tilts[gdp][srt], scifrcp[gdp][srt], 3,
                                            function='bspline',
                                            weights=np.sqrt(ivar)[gdp][srt],
                                            sigma=5., maxone=False,
                                            bspline_par=settings['skysub']['bspline'])
        # Just those in the slit
        in_slit = np.where(slitpix == slit+1)
        bgf_flat = arutils.func_val(bspl, tilts[in_slit].flatten(), 'bspline')
        #bgframe = bgf_flat.reshape(tilts.shape)
        bgframe[in_slit] = bgf_flat
        if msgs._debug['sky_sub']:
            plt_bspline_sky(tilts, scifrcp, bgf_flat, in_slit, gdp)
            debugger.set_trace()
    else:
        msgs.error('Not ready for this method for skysub {:s}'.format(method))

    if np.sum(np.isnan(bgframe)) > 0:
        msgs.warn("NAN in bgframe.  Replacing with 0")
        bad = np.isnan(bgframe)
        bgframe[bad] = 0.

    # Plot to make sure that the result is good
    return bgframe



def plt_bspline_sky(tilts, scifrcp, bgf_flat, inslit, gdp):
    # Setup
    srt = np.argsort(tilts[inslit].flatten())
    # Plot
    plt.close()
    plt.clf()
    ax = plt.gca()
    ax.scatter(tilts[gdp]*tilts.shape[0], scifrcp[gdp], marker='o')
    ax.plot(tilts[inslit].flatten()[srt]*tilts.shape[0], bgf_flat[srt], 'r-')
    plt.show()


