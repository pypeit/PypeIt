""" Module for methods related to tracing arc/sky lines across a slit/order
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy

import numpy as np

from scipy import interpolate

import matplotlib.pyplot as plt

from pypit import msgs
from pypit.core import ararc
from pypit import arutils
from pypit import arparse
from pypit import arpca
from pypit import arqa
from pypit import ardebug as debugger

try:
    from pypit import ginga
except ImportError:
    pass


def analyze_lines(msarc, trcdict, slit, pixcen, tilt_settings, maskval=-999999.9):
    # Analyze each spectral line
    aduse = trcdict["aduse"]
    arcdet = trcdict["arcdet"]
    xtfits = trcdict["xtfit"]
    ytfits = trcdict["ytfit"]
    wmasks = trcdict["wmask"]
    badlines = trcdict["badlines"]

    xtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
    ytilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
    mtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
    wtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval

    # For displaying later
    xmodel = []
    ymodel = []

    for j in range(arcdet.size):
        if not aduse[j]:
            continue
        xtfit = xtfits[j]
        ytfit = ytfits[j]
        wmask = wmasks[j]
        xint = int(xtfit[0])
        sz = (xtfit.size-1)//2

        # Trim if we are off the detector
        lastx = min(xint + 2 * sz + 1, msarc.shape[1])
        if (lastx-xint) < xtfit.size: # Cut down
            dx = (lastx-xint)-xtfit.size
            xtfit = xtfit[:dx]
            ytfit = ytfit[:dx]
            wmask = wmask[np.where(wmask < (xtfit.size+dx))]

        # Perform a scanning polynomial fit to the tilts
        wmfit = np.where(ytfit != maskval)
        if wmfit[0].size > tilt_settings['tilts']['order'] + 1:
            cmfit = arutils.func_fit(xtfit[wmfit], ytfit[wmfit],
                                     tilt_settings['tilts']['function'],
                                     tilt_settings['tilts']['order'],
                                     minv=0.0, maxv=msarc.shape[1] - 1.0)
            model = arutils.func_val(cmfit, xtfit, tilt_settings['tilts']['function'],
                                     minv=0.0, maxv=msarc.shape[1] - 1.0)
        else:
            aduse[j] = False
            badlines += 1
            continue

        # Can this actually happen??
        if maskval in model:
            # Model contains masked values
            aduse[j] = False
            badlines += 1
            continue

        # Perform a robust polynomial fit to the traces
        wmsk, mcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask],
                                              tilt_settings['tilts']['order'],
                                              function=tilt_settings['tilts']['function'],
                                              sigma=2.0, minv=0.0, maxv=msarc.shape[1] - 1.0)

        # Save model
        model = arutils.func_val(mcoeff, xtfit, tilt_settings['tilts']['function'],
                                 minv=0.0, maxv=msarc.shape[1] - 1.0)
        xmodel.append(xtfit)
        ymodel.append(model)

        # Save
        xtilt[xint:lastx, j] = xtfit / (msarc.shape[1] - 1.0)
        # These should be un-normalized for now
        pcen = pixcen[arcdet[j], slit]
        ytilt[xint:lastx, j] = model[pcen-int(xtfit[wmask[0]])]
        mtilt[xint:lastx, j] = model

    # Save
    trcdict['xmodel'] = xmodel
    trcdict['ymodel'] = ymodel
    trcdict["aduse"] = aduse

    # Return
    all_tilts = (xtilt, ytilt, mtilt, wtilt)
    return badlines, all_tilts


def new_tilts_image(tilts, lordloc, rordloc, pad, sz_y):
    """
    Using the tilt (assumed to be fit with a first order polynomial)
    generate an image of the tilts for each slit.

    Parameters
    ----------
    tilts : ndarray
      An (m x n) 2D array specifying the tilt (i.e. gradient) at each
      pixel along the spectral direction (m) for each slit (n).
    lordloc : ndarray
      Location of the left slit edges
    rordloc : ndarray
      Location of the right slit edges
    pad : int
      Set the tilts within each slit, and extend a number of pixels
      outside the slit edges (this number is set by pad).
    sz_y : int
      Number of detector pixels in the spatial direction.

    Returns
    -------
    tiltsimg : ndarray
      An image the same size as the science frame, containing values from 0-1.
      0/1 corresponds to the bottom/top of the detector (in the spectral direction),
      and constant wavelength is represented by a single value from 0-1.

    """
    sz_x, sz_o = tilts.shape
    dszx = (sz_x-1.0)

    tiltsimg = np.zeros((sz_x,sz_y), dtype=float)
    for o in range(sz_o):
        for x in range(sz_x):
            ow = (rordloc[x,o]-lordloc[x,o])/2.0
            oc = (rordloc[x,o]+lordloc[x,o])/2.0
            ymin = int(oc-ow) - pad
            ymax = int(oc+ow) + 1 + pad
            # Check we are in bounds
            if ymin < 0:
                ymin = 0
            elif ymax < 0:
                continue
            if ymax > sz_y-1:
                ymax = sz_y-1
            elif ymin > sz_y-1:
                continue
            # Set the tilt value at each pixel in this row
            for y in range(ymin, ymax):
                yv = (y-lordloc[x, o])/ow - 1.0
                tiltsimg[x,y] = (tilts[x,o]*yv + x)/dszx
    return tiltsimg


def trace_tilt(ordcen, rordloc, lordloc, det, msarc, slitnum, settings_spect,
               tilt_settings, censpec=None, maskval=-999999.9,
               trthrsh=1000.0, nsmth=0, method="fweight", wv_calib=None):
    """
    This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
               trthrsh=1000.0, nsmth=0):

    Parameters
    ----------
    slf
    det
    msarc
    slitnum : int
      Slit number, here indexed from 0
    censpec
    maskval
    trthrsh
    nsmth
    method : str (fweight or cc)

    Returns
    -------
    trcdict : dict

    """
    #if settings.argflag['trace']['slits']['tilts']['idsonly']:
    def pad_dict(indict):
        """ If an arc line is considered bad, fill the
        dictionary arrays with null values
        """
        indict["xtfit"].append(None)
        indict["ytfit"].append(None)
        indict["wmask"].append(None)
        return indict

    # from pypit import arcyutils
    dnum = arparse.get_dnum(det)

    msgs.work("Detecting lines for slit {0:d}".format(slitnum+1))
    tampl, tcent, twid, w, _ = ararc.detect_lines(censpec)
    satval = settings_spect[dnum]['saturation']*settings_spect[dnum]['nonlinear']
    # Order of the polynomials to be used when fitting the tilts.
    arcdet = (tcent[w]+0.5).astype(np.int)
    ampl = tampl[w]

    # Determine the best lines to use to trace the tilts
    ncont = 15
    aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    w = np.where(ampl >= trthrsh)
    aduse[w] = 1
    # Remove lines that are within ncont pixels
    nuse = np.sum(aduse)
    detuse = arcdet[aduse]
    idxuse = np.arange(arcdet.size)[aduse]
    olduse = aduse.copy()
    for s in range(nuse):
        w = np.where((np.abs(arcdet-detuse[s]) <= ncont) & (np.abs(arcdet-detuse[s]) >= 1.0))[0]
        for u in range(w.size):
            if ampl[w[u]] > ampl[olduse][s]:
                aduse[idxuse[s]] = False
                break
    # Restricted to ID lines? [introduced to avoid LRIS ghosts]
    if tilt_settings['tilts']['idsonly']:
        ids_pix = np.round(np.array(wv_calib[str(slitnum)]['xfit'])*(msarc.shape[0]-1))
        idxuse = np.arange(arcdet.size)[aduse]
        for s in idxuse:
            if np.min(np.abs(arcdet[s]-ids_pix)) > 2:
                msgs.info("Ignoring line at row={:d}".format(arcdet[s]))
                aduse[s] = False

    # Divide the detector into Nseg segments,
    # and find the brightest lines in each segment.
    # The total number of lines used to trace the tilts will be  = Nseg*Nuse + Nadd
    # Nseg = 4
    # Nuse = 2
    # Nadd = 8
    # segsz = msarc.shape[0]/float(Nseg)
    # aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    # for s in range(Nseg):
    #     w = np.where((arcdet > s*segsz) & (arcdet <= (s+1)*segsz))[0]
    #     segampl = tampl[w]
    #     asrt = np.argsort(segampl)[::-1]
    #     for u in range(Nuse):
    #         aduse[w[asrt[u]]] = True
    # # Now include some additional bright lines
    # asrt = np.argsort(tampl)[::-1]
    # s, u = 0, 0
    # while u < Nadd:
    #     if not aduse[asrt[s]]:
    #         aduse[asrt[s]] = True
    #         u += 1
    #     s += 1

    # Setup the trace dictionary
    trcdict = {"xtfit":[], "ytfit":[], "wmask":[], "arcdet":arcdet, "aduse":aduse, "badlines":0}

    msgs.info("Modelling arc line tilts with {0:d} arc lines".format(np.sum(aduse)))
    if np.sum(aduse) == 0:
        msgs.warn("No arc lines were deemed usable in slit {0:d} for spectral tilt".format(slitnum))
        return None
    # Go along each order and trace the tilts
    # Start by masking every row, then later unmask the rows with usable arc lines
    msgs.work("This next step could be multiprocessed to speed up the reduction")
    nspecfit = 3
    badlines = 0
    for j in range(arcdet.size):
        # For each detection in this order
        #msgs.info("Tracing tilt of arc line {0:d}/{1:d}".format(j+1, arcdet.size))
        # Check if this is a saturated line
        ysat = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1, ordcen[arcdet[j], slitnum]-nsmth:ordcen[arcdet[j], slitnum]+nsmth+1]
        if np.where(ysat > satval)[0].size != 0:
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        # Get the size of the slit
        sz = int(np.floor(np.abs(rordloc[arcdet[j], slitnum]-lordloc[arcdet[j], slitnum])/2.0)) - 2
        xtfit = np.zeros(2*sz+1)
        ytfit = np.ones(2*sz+1)*maskval  # Fitted centroid
        etfit = np.zeros(2*sz+1)  # Fitted centroid error
        mtfit = np.ones(2*sz+1, dtype=np.int)   # Mask of bad fits
        #apfit = np.zeros(2*sz+1)  # Fitted Amplitude
        xfit = np.arange(-nspecfit, nspecfit+1, 1.0)
        wfit = np.ones(xfit.size, dtype=np.float)
        tstcc = True  # A boolean to tell the loop once a good set of pixels has been found to cross-correlate with
        # Fit up
        pcen = arcdet[j]
        if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
            # Too close to the end of the spectrum
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        offchip = False
        centv = None
        for k in range(0, sz+1-nsmth):
            if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                offchip = True
                break
            if ordcen[pcen, slitnum]+k >= msarc.shape[1]:
                offchip = True
                break
            # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]+k]
            yfit = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], slitnum]+k-nsmth:ordcen[arcdet[j], slitnum]+k+nsmth+1]
            if np.size(yfit) == 0:
                offchip = True
                break
            if len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            # wgd = np.where((yfit<satval)&(yfit!=maskval))
            wgd = np.where(yfit == maskval)
            if wgd[0].size != 0:
                continue
            if method == "fweight":
                if centv is None:
                    centv = np.sum(yfit * (pcen+xfit))/np.sum(yfit)
                wfit[0] = 0.5 + (pcen-centv)
                wfit[-1] = 0.5 - (pcen-centv)
                sumxw = yfit * (pcen+xfit) * wfit
                sumw = yfit * wfit
                centv = np.sum(sumxw)/np.sum(sumw)
                fail = False
            elif method == "cc":
                # Get a copy of the array that will be used to cross-correlate
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]+k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],slitnum]+k-nsmth:ordcen[arcdet[j],slitnum]+k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                cc = np.correlate(ccyfit, yfit, mode='same')
                params, fail = arutils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
            xtfit[k+sz] = ordcen[arcdet[j], slitnum] + k
            ytfit[k+sz] = centv
            etfit[k+sz] = 0.02
            #apfit[k+sz] = params[0]
            if fail:
                mtfit[k+sz] = 1
            else:
                pcen = int(0.5 + centv)
                mtfit[k+sz] = 0
        if offchip:
            # Don't use lines that go off the chip (could lead to a bad trace)
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        for k in range(sz+1-nsmth, sz+1):
            xtfit[k+sz] = ordcen[arcdet[j], slitnum]+k
        # Fit down
        pcen = arcdet[j]
        centv = None
        for k in range(1,sz+1-nsmth):
            if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                offchip = True
                break
            if ordcen[pcen, slitnum]-k < 0:
                offchip = True
                break
            # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]-k]
            yfit = msarc[pcen - nspecfit:pcen + nspecfit + 1,
                   ordcen[arcdet[j], slitnum] - k - nsmth:ordcen[arcdet[j], slitnum] - k + nsmth + 1]
            if len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            if np.size(yfit) == 0:
                offchip = True
                break
            wgd = np.where(yfit == maskval)
            if wgd[0].size != 0:
                continue
            if method == "fweight":
                if centv is None:
                    centv = np.sum(yfit * (pcen+xfit))/np.sum(yfit)
                wfit[0] = 0.5 + (pcen-centv)
                wfit[-1] = 0.5 - (pcen-centv)
                sumxw = yfit * (pcen+xfit) * wfit
                sumw = yfit * wfit
                centv = np.sum(sumxw)/np.sum(sumw)
                fail = False
            elif method == "cc":
                # Get a copy of the array that will be used to cross-correlate
                # (testcc is probably already False from the Fit Up part of the loop, but it's best to be sure)
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]-k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],slitnum]-k-nsmth:ordcen[arcdet[j],slitnum]-k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                cc = np.correlate(ccyfit, yfit, mode='same')
                params, fail = arutils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
            xtfit[sz-k] = ordcen[arcdet[j], slitnum] - k
            ytfit[sz-k] = centv
            etfit[sz-k] = 0.02
            #apfit[sz-k] = params[0]
            if fail:
                mtfit[sz-k] = 1
            else:
                #from IPython import embed
                if np.isfinite(centv) is False: debugger.set_trace() #embed()
                pcen = int(0.5 + centv)
                mtfit[sz-k] = 0
        '''
        jxp_fix = False
        if jxp_fix:
            from desispec.bootcalib import trace_crude_init
            from desispec.bootcalib import trace_fweight as dbtf
            from pypit import ginga
            pcen = arcdet[j]
            img = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], slitnum]-sz:ordcen[arcdet[j], slitnum]+sz+1]
            rot_img = np.rot90(img,3) - 980.  # Bias!
            #
            xcen = np.array([6.6]) # 8.5, 173
            ypass = 173
            # Tune-up first
            for ii in range(4):
                xcen, xsig = dbtf(rot_img, xcen, ycen=np.array([ypass]).astype(int), invvar=None, radius=2.)
            xset, xerr = trace_crude_init(rot_img, xcen, ypass, invvar=None, radius=3., maxshift0=0.5, maxshift=0.15, maxerr=0.2)
            #xcen, xsig = dbtf(rot_img, np.array([6.5,6.5]), ycen=np.array([173,173]).astype(int), invvar=None, radius=2.)
            # Convert back
            y0 = pcen+nspecfit
            ycrude = y0 - xset[:,0]
            #debugger.set_trace()
            #cytfit = ytfit.copy()
            #cytfit[np.where(ytfit < 0)] = np.median(cytfit)
            #debugger.xplot(np.arange(img.shape[1]), ycrude, cytfit)
            #
            #trcdict['save_yt'] = ytfit.copy()
            assert ycrude.size == ytfit.size
            ytfit = ycrude
            #debugger.set_trace()
        '''

        if offchip:
            # Don't use lines that go off the chip (could lead to a bad trace)
            aduse[j] = False
            badlines += 1
            trcdict = pad_dict(trcdict)
            continue
        for k in range(sz+1-nsmth, sz+1):
            xtfit[sz-k] = ordcen[arcdet[j], slitnum]-k

        wmask = np.where(mtfit == 0)
        ytfit[np.where(mtfit == 1)] = maskval

        # Append the trace information into the dictionary
        trcdict["xtfit"].append(xtfit.copy())
        trcdict["ytfit"].append(ytfit.copy())
        trcdict["wmask"].append(wmask[0].copy())
    trcdict["aduse"] = aduse
    trcdict["badlines"] = badlines
    msgs.info("Completed spectral tilt tracing".format(np.sum(aduse)))
    return trcdict



def trace_fweight(fimage, xinit, ltrace=None, rtraceinvvar=None, radius=3.):
    """ Python port of trace_fweight.pro from IDLUTILS

    Parameters:
    -----------
    fimage: 2D ndarray
      Image for tracing
    xinit: ndarray
      Initial guesses for x-trace
    invvar: ndarray, optional
      Inverse variance array for the image
    radius: float, optional
      Radius for centroiding; default to 3.0
    """

    # Init
    nx = fimage.shape[1]
    ny = fimage.shape[0]
    ncen = len(xinit)
    xnew = copy.deepcopy(xinit)
    xerr = np.zeros(ncen) + 999.

    ycen = np.arange(ny, dtype=int)
    invvar = 0. * fimage + 1.
    x1 = xinit - radius + 0.5
    x2 = xinit + radius + 0.5
    ix1 = np.floor(x1).astype(int)
    ix2 = np.floor(x2).astype(int)

    fullpix = int(np.maximum(np.min(ix2-ix1)-1, 0))
    sumw = np.zeros(ny)
    sumxw = np.zeros(ny)
    sumwt = np.zeros(ny)
    sumsx1 = np.zeros(ny)
    sumsx2 = np.zeros(ny)
    qbad = np.array([False]*ny) 

    if invvar is None: 
        invvar = np.zeros_like(fimage) + 1. 

    # Compute
    for ii in range(0,fullpix+3):
        spot = ix1 - 1 + ii
        ih = np.clip(spot,0,nx-1)
        xdiff = spot - xinit
        #
        wt = np.clip(radius - np.abs(xdiff) + 0.5,0,1) * ((spot >= 0) & (spot < nx))
        sumw = sumw + fimage[ycen,ih] * wt
        sumwt = sumwt + wt
        sumxw = sumxw + fimage[ycen,ih] * xdiff * wt
        var_term = wt**2 / (invvar[ycen,ih] + (invvar[ycen,ih] == 0))
        sumsx2 = sumsx2 + var_term
        sumsx1 = sumsx1 + xdiff**2 * var_term
        #qbad = qbad or (invvar[ycen,ih] <= 0)
        qbad = np.any([qbad, invvar[ycen,ih] <= 0], axis=0)

    # Fill up
    good = (sumw > 0) & (~qbad)
    if np.sum(good) > 0:
        delta_x = sumxw[good]/sumw[good]
        xnew[good] = delta_x + xinit[good]
        xerr[good] = np.sqrt(sumsx1[good] + sumsx2[good]*delta_x**2)/sumw[good]

    bad = np.any([np.abs(xnew-xinit) > radius + 0.5, xinit < radius - 0.5, xinit > nx - 0.5 - radius], axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Return
    return xnew, xerr


def echelle_tilt(slf, msarc, det, settings_argflag, settings_spect, pcadesc="PCA trace of the spectral tilts", maskval=-999999.9):
    """ Determine the spectral tilts in each echelle order

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    msarc : numpy ndarray
      Wavelength calibration frame that will be used to trace constant wavelength
    det : int
      Index of the detector
    pcadesc : str (optional)
      A description of the tilts
    maskval : float (optional)
      Mask value used in numpy arrays

    Returns
    -------
    tilts : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    satmask : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    arccen, maskslit, satmask = ararc.get_censpec(slf._lordloc[det-1], slf._rordloc[det-1],
                                                      slf._pixlocn[det-1], msarc, det, settings_spect, gen_satmask=True)
    # If the user sets no tilts, return here
    if settings_argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satmask, None
    norders = maskslit.size

    # Now model the tilt for each slit
    tiltang, centval = None, None
    slitcnt = 0
    for o in range(norders):
        if maskslit[o] == 1:
            continue
        # Determine the tilts for this slit
        trcdict = trace_tilt(slf, det, msarc, o, censpec=arccen[:, o])
        slitcnt += 1
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        aduse = trcdict["aduse"]
        arcdet = trcdict["arcdet"]
        if tiltang is None:
            tiltang = maskval * np.ones((aduse.size, norders))
            centval = maskval * np.ones((aduse.size, norders))
            totnum = aduse.size
        else:
            if aduse.size > totnum:
                tiltang = np.append(tiltang, maskval * np.ones((aduse.size-totnum, norders)), axis=0)
                centval = np.append(centval, maskval * np.ones((aduse.size-totnum, norders)), axis=0)
                totnum = aduse.size
        for j in range(aduse.size):
            if not aduse[j]:
                continue
            lrdiff = slf._rordloc[det-1][arcdet[j], o] - slf._lordloc[det-1][arcdet[j], o]
            xtfit = 2.0 * (trcdict["xtfit"][j]-slf._lordloc[det-1][arcdet[j], o])/lrdiff - 1.0
            ytfit = trcdict["ytfit"][j]
            wmask = trcdict["wmask"][j]
            if wmask.size < settings_argflag['trace']['slits']['tilts']['order']+2:
                maskslit[o] = 1
                continue
            null, tcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask],
                                                  settings_argflag['trace']['slits']['tilts']['order'], sigma=2.0)
            # Save the tilt angle
            tiltang[j, o] = tcoeff[1]  # tan(tilt angle)
            centval[j, o] = tcoeff[0]  # centroid of arc line

    msgs.info("Fitting tilt angles")
    tcoeff = np.ones((settings_argflag['trace']['slits']['tilts']['disporder'] + 1, tiltang.shape[1]))
    maskord = np.where(maskslit == 1)[0]
    extrap_ord = np.zeros(norders)
    for o in range(norders):
        if o in maskord:
            extrap_ord[o] = 1
            continue
        w = np.where(tiltang[:, o] != maskval)
        if np.size(w[0]) <= settings_argflag['trace']['slits']['tilts']['disporder'] + 2:
            extrap_ord[o] = 1.0
            maskord = np.append(maskord, o)
        else:
            null, tempc = arutils.robust_polyfit(centval[:, o][w], tiltang[:, o][w],
                                                 settings_argflag['trace']['slits']['tilts']['disporder'],
                                                 function=settings_argflag['trace']['slits']['function'], sigma=2.0,
                                                 minv=0.0, maxv=msarc.shape[0] - 1)
            tcoeff[:, o] = tempc
    # Sort which orders are masked
    maskord.sort()
    xv = np.arange(msarc.shape[0])
    tiltval = arutils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'],
                               minv=0.0, maxv=msarc.shape[0] - 1).T
    ofit = settings_argflag['trace']['slits']['tilts']['params']
    lnpc = len(ofit) - 1
    if np.sum(1.0 - extrap_ord) > ofit[0] + 1:  # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral tilts")
        ordsnd = np.arange(norders) + 1.0
        xcen = xv[:, np.newaxis].repeat(norders, axis=1)
        fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=settings_argflag['trace']['slits']['function'])
        if not msgs._debug['no_qa']:
            #pcadesc = "Spectral Tilt PCA"
#            arqa.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
            arpca.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(norders)
        extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=settings_argflag['trace']['slits']['function'])
        tilts = extrap_tilt
#        arqa.pca_arctilt(slf, tiltang, centval, tilts)
        arpca.pca_arctilt(slf, tiltang, centval, tilts)
    else:
        outpar = None
        msgs.warn("Could not perform a PCA when tracing the order tilts" + msgs.newline() +
                  "Not enough well-traced orders")
        msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
        xtiltfit = np.array([])
        ytiltfit = np.array([])
        for o in range(tiltang.shape[1]):
            w = np.where(tiltang[:, o] != -999999.9)
            if np.size(w[0]) != 0:
                xtiltfit = np.append(xtiltfit, centval[:, o][w])
                ytiltfit = np.append(ytiltfit, tiltang[:, o][w])
        if np.size(xtiltfit) > settings_argflag['trace']['slits']['tilts']['disporder'] + 2:
            tcoeff = arutils.func_fit(xtiltfit, ytiltfit, settings_argflag['trace']['slits']['function'],
                                      settings_argflag['trace']['slits']['tilts']['disporder'],
                                      minv=0.0, maxv=msarc.shape[0] - 1)
            tiltval = arutils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'], minv=0.0,
                                       maxv=msarc.shape[0] - 1)
            tilts = tiltval[:, np.newaxis].repeat(tiltang.shape[1], axis=1)
        else:
            msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
            msgs.info("Assuming there is no tilt")
            tilts = np.zeros_like(slf._lordloc)

    # Generate tilts image
#    print('calling tilts_image')
#    t = time.clock()
#    _tiltsimg = arcytrace.tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1],
#                                     settings.argflag['trace']['slits']['pad'], msarc.shape[1])
#    print('Old tilts_image: {0} seconds'.format(time.clock() - t))
#    t = time.clock()
    tiltsimg = new_tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1],
                                settings_argflag['trace']['slits']['pad'], msarc.shape[1])
#    print('New tilts_image: {0} seconds'.format(time.clock() - t))
#    assert np.sum(_tiltsimg != tiltsimg) == 0, 'Difference between old and new tilts_image'

    return tiltsimg, satmask, outpar


def multislit_tilt(msarc, lordloc, rordloc, pixlocn, pixcen, slitpix, det,
                   maskslits, tilt_settings, settings_spect, setup,
                   maskval=-999999.9, doqa=False, wv_calib=None):
    """ Determine the spectral tilt of each slit in a multislit image

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    msarc : numpy ndarray
      Wavelength calibration frame that will be used to trace constant wavelength
    det : int
      Index of the detector
    maskval : float (optional)
      Mask value used in numpy arrays
    doqa : bool, optional
      Output QA files.  These can be many files and slow for
      lots of slits

    Returns
    -------
    tilts : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    satmask : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    # tilt_settings['tilts']['order'] : settings.argflag['trace']['slits']['tilts']['order']
    # ofit = settings.argflag['trace']['slits']['tilts']['params']
    arccen, arc_maskslit, _ = ararc.get_censpec(lordloc, rordloc, pixlocn, msarc, det,
                                            settings_spect, gen_satmask=False)
    satmask = np.zeros_like(pixcen)

    ordcen = pixcen.copy()

    # maskslit
    if maskslits is not None:
        mask = maskslits & (arc_maskslit==1)
    else:
        mask = arc_maskslit

    # TODO -- NEED TO PASS THIS BACK!?
    #slf._maskslits[det-1] = mask

    gdslits = np.where(mask == 0)[0]

    # Final tilts image
    final_tilts = np.zeros_like(msarc)

    # Now trace the tilt for each slit
    for slit in gdslits:
        # Determine the tilts for this slit
        trcdict = trace_tilt(pixcen, rordloc, lordloc, det, msarc, slit, settings_spect,
                             tilt_settings, censpec=arccen[:, slit], nsmth=3, wv_calib=wv_calib,
                             trthrsh=tilt_settings['tilts']['trthrsh'])
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        if msgs._debug['tilts']:
            debugger.chk_arc_tilts(msarc, trcdict, sedges=(lordloc[:,slit], rordloc[:,slit]))
            debugger.set_trace()

        # Extract information from the trace dictionary
        aduse = trcdict["aduse"]
        arcdet = trcdict["arcdet"]

        # Analyze spec lines
        badlines, maskrows, tcoeff, all_tilts = analyze_spec_lines(msarc, slit, trcdict, ordcen, tilt_settings)
        if badlines != 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")

        # Prepare polytilts
        polytilts, outpar = prepare_polytilts(msarc, slit, maskrows, tcoeff, all_tilts, tilt_settings, setup=setup)

        '''
        if tilt_settings['tilts']['method'].lower() == "interp":
            tilts = tilts_interp(ordcen, slit, all_tilts, polytilts, arcdet, aduse, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spline":
            tilts = tilts_spline(all_tilts, arcdet, aduse, polytilts, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spca":
            tilts = tilts_spca(msarc, polytilts, ordcen, slit, arcdet, aduse, rordloc, lordloc)
        elif tilt_settings['tilts']['method'].lower() == "pca":
            tilts = polytilts.copy()
        '''

        # Save into final_tilts
        word = np.where(slitpix == slit+1)
        #final_tilts[word] = tilts[word]
        final_tilts[word] = polytilts[word]

        # Now do the QA
        if doqa:
            tiltsplot, ztilto, xdat = prep_tilts_qa(msarc, all_tilts, tilts, arcdet, ordcen, slit, setup)
            msgs.info("Plotting arc tilt QA")
            plot_orderfits(setup, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
                   textplt="Arc line", maxp=9, desc="Arc line spectral tilts",
                   maskval=maskval, slit=slit)
    # Finish
    return final_tilts, satmask, outpar


def fit_tilts(msarc, slit, all_tilts, tilt_settings, maskval=-999999.9, setup=None, doqa=True, show_QA=False):
    # Unpack
    xtilt, ytilt, mtilt, wtilt = all_tilts
    #
    fitxy = [tilt_settings['tilts']['order']+1, tilt_settings['tilts']['yorder']]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(tilt_settings['tilts']['func2D']))
    wgd = np.where(xtilt != maskval)
    # Invert
    coeff2 = arutils.polyfit2d_general(xtilt[wgd], mtilt[wgd]/(msarc.shape[0]-1),
                                       mtilt[wgd]-ytilt[wgd], fitxy,
                                              minx=0., maxx=1., miny=0., maxy=1.,
                                              function=tilt_settings['tilts']['func2D'])
    polytilts = arutils.polyval2d_general(coeff2, np.linspace(0.0, 1.0, msarc.shape[1]),
                                          np.linspace(0.0, 1.0, msarc.shape[0]),
                                          minx=0., maxx=1., miny=0., maxy=1.,
                                          function=tilt_settings['tilts']['func2D'])

    # TODO -- Add a rejection iteration (or two)

    # Residuals
    xv2 = arutils.scale_minmax(xtilt[wgd], minx=0., maxx=1)
    yv2 = arutils.scale_minmax(mtilt[wgd]/(msarc.shape[0]-1), minx=0., maxx=1)
    yfit = np.polynomial.legendre.legval2d(xv2, yv2, coeff2)
    res2 = (mtilt[wgd]-ytilt[wgd]) - yfit
    msgs.info("RMS (pixels): {}".format(np.std(res2)))

    # QA
    if doqa:
        plot_tiltres(setup, mtilt[wgd], ytilt[wgd], yfit, slit=slit, show_QA=show_QA)

    # y normalization and subtract
    ynorm = np.outer(np.linspace(0., 1., msarc.shape[0]), np.ones(msarc.shape[1]))
    polytilts = ynorm - polytilts/(msarc.shape[0]-1)

    # Return
    outpar = None
    return polytilts, outpar


def slit_image(slf, det, scitrace, obj, tilts=None):
    """ Generate slit image for a given object
    Ignores changing plate scale (for now)
    The slit is approximated as a straight line in this calculation
    which should be reasonably accurate.  Better, the object profile
    generated from this approximation is applied in the same fashion
    so that the 'error' is compensated for.
    Parameters
    ----------
    slf
    det
    scitrace
    obj
    Returns
    -------
    slit_img : ndarray
    """
    # Setup
    if tilts is None:
        tilts = slf._tilts[det-1]
    ximg = np.outer(np.ones(tilts.shape[0]), np.arange(tilts.shape[1]))
    dypix = 1./tilts.shape[0]
    #  Trace
    xtrc = np.round(scitrace['traces'][:, obj]).astype(int)
    wch = np.where((xtrc >= 0) & (xtrc <= tilts.shape[1]-1))
    msgs.work("Use 2D spline to evaluate tilts")
    trc_tilt = np.zeros(tilts.shape[0], dtype=np.float)
    trc_tilt[wch] = tilts[np.arange(tilts.shape[0])[wch], xtrc[wch]]
    trc_tilt_img = np.outer(trc_tilt, np.ones(tilts.shape[1]))
    # Slit image
    msgs.work("Should worry about changing plate scale")
    dy = (tilts - trc_tilt_img)/dypix  # Pixels
    dx = ximg - np.outer(scitrace['traces'][:, obj], np.ones(tilts.shape[1]))
    slit_img = np.sqrt(dx**2 - dy**2)
    neg = dx < 0.
    slit_img[neg] *= -1
    # Return
    return slit_img


def tilts_interp(ordcen, slit, all_tilts, polytilts, arcdet, aduse, msarc):
    # Unpack
    xtilt, ytilt, ztilt, mtilt, wtilt = all_tilts
    #
    msgs.info("Interpolating and Extrapolating the tilts")
    xspl = np.linspace(0.0, 1.0, msarc.shape[1])
    # yspl = np.append(0.0, np.append(arcdet[np.where(aduse)]/(msarc.shape[0]-1.0), 1.0))
    # yspl = np.append(0.0, np.append(polytilts[arcdet[np.where(aduse)], msarc.shape[1]/2], 1.0))
    ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse), slit]]])
    yspl = np.append(0.0, np.append(ycen, 1.0))
    zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
    zspl[:, 1:-1] = mtilt[:, np.where(aduse)[0]]
    # zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
    zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
    zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
    # Make sure the endpoints are set to 0.0 and 1.0
    zspl[:, 0] -= zspl[ordcen[0, 0], slit]
    zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, slit], -1] + 1.0
    # Prepare the spline variables
    # if False:
    #     pmin = 0
    #     pmax = -1
    # else:
    #     pmin = int(max(0, np.min(slf._lordloc[det-1])))
    #     pmax = int(min(msarc.shape[1], np.max(slf._rordloc[det-1])))
    # xsbs = np.outer(xspl, np.ones(yspl.size))
    # ysbs = np.outer(np.ones(xspl.size), yspl)
    # zsbs = zspl[wgd]
    # Restrict to good portion of the image
    tiltspl = interpolate.RectBivariateSpline(xspl, yspl, zspl, kx=3, ky=3)
    yval = np.linspace(0.0, 1.0, msarc.shape[0])
    tilts = tiltspl(xspl, yval, grid=True).T
    # Return
    return tilts


def tilts_spline(all_tilts, arcdet, aduse, polytilts, msarc, use_mtilt=False, maskval=-999999.9):
    msgs.info("Performing a spline fit to the tilts")
    # Unpack
    xtilt, ytilt, ztilt, mtilt, wtilt = all_tilts
    #
    wgd = np.where((ytilt != maskval) & (ztilt != maskval))
    txsbs = xtilt[wgd]
    tysbs = ytilt[wgd]
    if use_mtilt:
        tzsbs = mtilt[wgd]/(msarc.shape[0]-1)
    else:
        tzsbs = ztilt[wgd]
    twsbs = wtilt[wgd]

    # Append the end points
    wlo = np.where((ytilt == np.min(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
    whi = np.where((ytilt == np.max(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
    xlo = (xtilt[wlo] * (msarc.shape[1] - 1.0)).astype(np.int)
    xhi = (xtilt[whi] * (msarc.shape[1] - 1.0)).astype(np.int)
    xsbs = np.append(xtilt[wlo], np.append(txsbs, xtilt[whi]))
    ysbs = np.append(np.zeros(wlo[0].size), np.append(tysbs, np.ones(whi[0].size)))

    zlo = ztilt[wlo] + polytilts[0, xlo] - polytilts[arcdet[np.where(aduse)[0][0]], xlo]
    zhi = ztilt[whi] + polytilts[-1, xhi] - polytilts[arcdet[np.where(aduse)[0][-1]], xhi]
    zsbs = np.append(zlo, np.append(tzsbs, zhi))
    wsbs = np.append(wtilt[wlo], np.append(twsbs, wtilt[whi]))
    # Generate the spline curve
    tiltspl = interpolate.SmoothBivariateSpline(xsbs, zsbs, ysbs, w=wsbs, kx=3, ky=3,
                                                s=xsbs.size)#, bbox=[0.0, 1.0, min(zsbs.min(),0.0), max(zsbs.max(),1.0)])
    xspl = np.linspace(0.0, 1.0, msarc.shape[1])
    yspl = np.linspace(0.0, 1.0, msarc.shape[0])
    tilts = tiltspl(xspl, yspl, grid=True).T
    #
    tmp = (msarc.shape[0]-1)*tiltspl(np.arange(1490, 1499)/(msarc.shape[1]-1),
                                     np.array([968.5, 969.34])/(msarc.shape[0]-1), grid=True)
    tmp2 = (msarc.shape[0]-1)*ytilt[1497,28]
    #tmp3 = tiltspl(xsbs, zsbs, grid=True)
    print(tmp, tmp2)
    debugger.set_trace()
    '''
    if msgs._debug['tilts']:
        tiltqa = tiltspl(xsbs, zsbs, grid=False)
        plt.clf()
        # plt.imshow((zsbs-tiltqa)/zsbs, origin='lower')
        # plt.imshow((ysbs-tiltqa)/ysbs, origin='lower')
        plt.plot(xsbs, (ysbs - tiltqa) / ysbs, 'bx')
        plt.plot(xsbs, 1.0 / (wsbs * ysbs), 'r-')
        plt.ylim(-5e-2, 5e-2)
        # plt.colorbar()
        plt.show()
        debugger.set_trace()
    '''
    # Return
    return tilts


def tilts_spca(msarc, polytilts, ordcen, slit, arcdet, aduse, rordloc, lordloc):
    msgs.info("Performing a spca analysis of the tilts")

    # Slit position
    xspl = np.linspace(0.0, 1.0, msarc.shape[1])
    # Trace positions down center of the order
    ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)], slit:slit+1]])
    yspl = np.append(0.0, np.append(ycen, 1.0))
    # Trace positions as measured+modeled
    zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
    zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
    #   Pad the ends with a fake measurement to avoid run-away
    zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
    zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
    # Make sure the endpoints are set to be nearly 0.0 and 1.0
    zspl[:, 0] -= zspl[ordcen[0, slit], 0]
    zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, slit], -1] + 1.0
    # Prepare the spline variables
    if False:
        pmin = 0
        pmax = -1
    else:
        pmin = int(max(0, np.min(lordloc)))
        pmax = int(min(msarc.shape[1], np.max(rordloc)))
    xsbs = np.outer(xspl, np.ones(yspl.size))[pmin:pmax, :]
    ysbs = np.outer(np.ones(xspl.size), yspl)[pmin:pmax, :]  * (msarc.shape[0]-1)
    zsbs = zspl[pmin:pmax, :]
    # Spline
    msgs.work('Consider adding weights to SmoothBivariate in spca')
    # Note how y and z are inverted (intentionally)!!
    tiltspl = interpolate.SmoothBivariateSpline(xsbs.flatten(), zsbs.flatten(),
                                                ysbs.flatten(), kx=3, ky=3, s=xsbs.size/2)
    # Finish
    yval = np.linspace(0.0, 1.0, msarc.shape[0])
    tilts = tiltspl(xspl, yval, grid=True).T
    #
    tmp = (msarc.shape[0]-1)*tiltspl(np.arange(1490, 1499)/(msarc.shape[1]-1),
                                     np.array([968.5, 969.34])/(msarc.shape[0]-1), grid=True)
    tmp2 = 940.
    print(tmp, tmp2)
    coeff = arutils.polyfit2d_general(xsbs, zsbs, ysbs, [5,5],
                                      minx=0., maxx=1., miny=0., maxy=1.,
                                      function='legendre')
    #tmpa = (msarc.shape[0]-1)*arutils.polyval2d_general(coeff, np.arange(1490, 1499)/(msarc.shape[1]-1),
    #tmpb = np.polynomial.legendre.legval2d(xsbs.flatten(), zsbs.flatten(), coeff)
    tmpa = arutils.polyval2d_general(coeff, np.arange(1490, 1499)/(msarc.shape[1]-1),
                                     np.array([968.5, 969.34])/(msarc.shape[0]-1),
                                     minx=0., maxx=1., miny=0., maxy=1.,
                                     function='legendre')
    print(tmpa)
    debugger.set_trace()
    # Return
    return tilts



def prep_tilts_qa(msarc, all_tilts, tilts, arcdet, ordcen, slit, maskval=-999999.9):
    # Unpack
    xtilt, ytilt, ztilt, mtilt, wtilt = all_tilts
    #
    msgs.info("Preparing arc tilt QA data")
    tiltsplot = tilts[arcdet, :].T
    tiltsplot *= (msarc.shape[0] - 1.0)
    # Shift the plotted tilts about the centre of the slit
    ztilto = ztilt.copy()
    adj = np.diag(tilts[arcdet, ordcen[arcdet, slit:slit+1]])
    zmsk = np.where(ztilto == maskval)
    # This magic is likely correct  ;)
    ztilto = 2.0 * np.outer(np.ones(ztilto.shape[0]), adj) - ztilto
    ztilto[zmsk] = maskval
    ztilto[np.where(ztilto != maskval)] *= (msarc.shape[0] - 1.0)
    for i in range(arcdet.size):
        w = np.where(ztilto[:, i] != maskval)
        if w[0].size != 0:
            twa = (xtilt[w[0], i] * (msarc.shape[1] - 1.0) + 0.5).astype(np.int)
            # fitcns = np.polyfit(w[0], ztilt[w[0], i] - tiltsplot[twa, i], 0)[0]
            fitcns = np.median(ztilto[w[0], i] - tiltsplot[twa, i])
            # if abs(fitcns) > 1.0:
            #     msgs.warn("The tilt of Arc Line {0:d} might be poorly traced".format(i+1))
            tiltsplot[:, i] += fitcns

    xdat = xtilt.copy()
    xdat[np.where(xdat != maskval)] *= (msarc.shape[1] - 1.0)

    #    arqa.plot_orderfits(slf, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
    #                        textplt="Arc line", maxp=9, desc="Arc line spectral tilts", maskval=maskval)
    return tiltsplot, ztilto, xdat


def plot_orderfits(setup, model, ydata, xdata=None, xmodl=None, textplt="Slit",
                   maxp=4, desc="", maskval=-999999.9, slit=None):
    """ Generate a QA plot for the blaze function fit to each slit
    Or the arc line tilts

    Parameters
    ----------
    slf : class
      Science Exposure class
    model : ndarray
      (m x n) 2D array containing the model blaze function (m) of a flat frame for each slit (n)
    ydata : ndarray
      (m x n) 2D array containing the extracted 1D spectrum (m) of a flat frame for each slit (n)
    xdata : ndarray, optional
      x values of the data points
    xmodl : ndarry, optional
      x values of the model points
    textplt : str, optional
      A string printed above each panel
    maxp : int, (optional)
      Maximum number of panels per page
    desc : str, (optional)
      A description added to the top of each page
    maskval : float, (optional)
      Value used in arrays to indicate a masked value
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Outfil
    method = inspect.stack()[0][3]
    if 'Arc' in desc:
        method += '_Arc'
    elif 'Blaze' in desc:
        method += '_Blaze'
    else:
        msgs.bug("Unknown type of order fits.  Currently prepared for Arc and Blaze")
    outroot = arqa.set_qa_filename(setup, method, slit=slit)
    #
    npix, nord = ydata.shape
    pages, npp = arqa.get_dimen(nord, maxp=maxp)
    if xdata is None: xdata = np.arange(npix).reshape((npix, 1)).repeat(nord, axis=1)
    if xmodl is None: xmodl = np.arange(model.shape[0])
    # Loop through all pages and plot the results
    ndone = 0
    axesIdx = True
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][0] == 1 and pages[i][1] == 1: axesIdx = False
            elif pages[i][1] == 1: ind = (ipx,)
            elif pages[i][0] == 1: ind = (ipy,)
            else: ind = (ipy, ipx)
            if axesIdx:
                axes[ind].plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes[ind].plot(xmodl, model[:,ndone+j], 'r-')
            else:
                axes.plot(xdata[:,ndone+j], ydata[:,ndone+j], 'bx', drawstyle='steps')
                axes.plot(xmodl, model[:,ndone+j], 'r-')
            ytmp = ydata[:,ndone+j]
            gdy = ytmp != maskval
            ytmp = ytmp[gdy]
            if ytmp.size != 0:
                amn = min(np.min(ytmp), np.min(model[gdy,ndone+j]))
            else:
                amn = np.min(model[:,ndone+j])
            if ytmp.size != 0:
                amx = max(np.max(ytmp), np.max(model[gdy,ndone+j]))
            else: amx = np.max(model[:,ndone+j])
            # Restrict to good pixels
            xtmp = xdata[:,ndone+j]
            gdx = xtmp != maskval
            xtmp = xtmp[gdx]
            if xtmp.size == 0:
                xmn = np.min(xmodl)
                xmx = np.max(xmodl)
            else:
                xmn = np.min(xtmp)
                xmx = np.max(xtmp)
                #xmn = min(np.min(xtmp), np.min(xmodl))
                #xmx = max(np.max(xtmp), np.max(xmodl))
            if axesIdx:
                axes[ind].axis([xmn, xmx, amn-1, amx+1])
                axes[ind].set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            else:
                axes.axis([xmn, xmx, amn, amx])
                axes.set_title("{0:s} {1:d}".format(textplt, ndone+j+1))
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Delete the unnecessary axes
        if axesIdx:
            for j in range(npp[i], axes.size):
                if pages[i][1] == 1: ind = (ipx,)
                elif pages[i][0] == 1: ind = (ipy,)
                else: ind = (ipy, ipx)
                f.delaxes(axes[ind])
                if ipx == pages[i][0]:
                    ipx = 0
                    ipy += 1
        ndone += npp[i]
        # Save the figure
        if axesIdx: axsz = axes.size
        else: axsz = 1.0
        if pages[i][1] == 1 or pages[i][0] == 1: ypngsiz = 11.0/axsz
        else: ypngsiz = 11.0*axes.shape[0]/axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        if desc != "":
            pgtxt = ""
            if len(pages) != 1:
                pgtxt = ", page {0:d}/{1:d}".format(i+1, len(pages))
            f.suptitle(desc + pgtxt, y=1.02, size=16)
        f.tight_layout()
        outfile = outroot+'{:03d}.png'.format(i)
        plt.savefig(outfile, dpi=200)
        plt.close()
        f.clf()
    del f

    plt.rcdefaults()

    return



def plot_tiltres(setup, mtilt, ytilt, yfit, slit=None, outfile=None, show_QA=False):
    """ Generate a QA plot of the residuals for the fit to the tilts
    One slit at a time

    Parameters
    ----------
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Outfil
    method = inspect.stack()[0][3]
    if (outfile is None) and (not show_QA):
        outfile = arqa.set_qa_filename(setup, method, slit=slit)

    # Setup
    plt.figure(figsize=(8, 4.0))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (mtilt-ytilt) - yfit
    ax.scatter(mtilt, res)

    rms = np.std(res)
    ax.text(0.90, 0.90, 'Slit {:d}:  RMS (pix) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black')
    # Label
    ax.set_xlabel('Row')
    ax.set_ylabel('Residual (pix)')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if show_QA:
        plt.show()
    else:
        plt.savefig(outfile, dpi=400)
    plt.close()

    plt.rcdefaults()

    return


