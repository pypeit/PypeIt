""" Module for sky subtraction
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from matplotlib import pyplot as plt


from pypit import msgs

from pypit import arutils
from pypit import arpixels

from pypit import ardebug as debugger


def bg_subtraction_slit(tslits_dict, pixlocn,
                        slit, tilts, sciframe, varframe, bpix, crpix,
                        settings_skysub,
                        tracemask=None,
                        rejsigma=3.0, maskval=-999999.9,
                        method='bspline'):
    """ Extract a science target and background flux
    :param slf:
    :param sciframe:
    :param varframe:
    :return:
    """
    # Unpack
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
                                            bspline_par=settings_skysub['bspline'])
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


