""" Module for methods related to tracing arc/sky lines across a slit/order
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import inspect
import copy

import numpy as np

from scipy import interpolate

import matplotlib.pyplot as plt

from pypeit import msgs
from pypeit.core import arc
from pypeit import utils
from pypeit.core import parse
from pypeit.core import pca
from pypeit.core import qa
from pypeit.core import trace_slits
from pypeit.core import extract
from pypeit import debugger

try:
    from pypeit import ginga
except ImportError:
    pass


def analyze_lines(msarc, trcdict, slit, pixcen, order=2, function='legendre', maskval=-999999.9):
    """
    .. todo::
        This needs a docstring!
    """
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
        if wmfit[0].size > order + 1:
            cmfit = utils.func_fit(xtfit[wmfit], ytfit[wmfit], function, order, minx=0.0,
                                     maxx=msarc.shape[1] - 1.0)
            model = utils.func_val(cmfit, xtfit, function, minx=0.0, maxx=msarc.shape[1] - 1.0)
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
        wmsk, mcoeff = utils.robust_polyfit(xtfit[wmask], ytfit[wmask], order, function=function,
                                              sigma=2.0, minx=0.0, maxx=msarc.shape[1] - 1.0)

        # Save model
        model = utils.func_val(mcoeff, xtfit, function, minx=0.0, maxx=msarc.shape[1] - 1.0)
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


def tilts_image(tilts, lordloc, rordloc, pad, sz_y):
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

def trace_tilt(ordcen, rordloc, lordloc, det, msarc, slitnum, satval,
               idsonly=False, censpec=None, maskval=-999999.9, tracethresh=20.0,
               nsmth=0, method="fweight", wv_calib=None, nonlinear_counts = 1e10):

    """
    This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
               tracethresh=1000.0, nsmth=0):
    # TODO Please expand these docs! This is no simple code.
    Parameters
    ----------
    slf
    det
    msarc
    slitnum : int
      Slit number, here indexed from 0
    censpec
    maskval
    tracethresh : float, optional
      This is now a significance whereas before it was an absolute threshold in counts
    nsmth
    method : str (fweight or cc)
    Returns
    -------
    trcdict : dict
    """
    def pad_dict(indict):
        """ If an arc line is considered bad, fill the
        dictionary arrays with null values
        """
        indict["xtfit"].append(None)
        indict["ytfit"].append(None)
        indict["wmask"].append(None)
        return indict

    dnum = parse.get_dnum(det)

    msgs.work("Detecting lines for slit {0:d}".format(slitnum+1))
    tampl, tampl_cont, tcent, twid, _, w, _ , tnsig = arc.detect_lines(censpec, fit_frac_fwhm=1.75, nonlinear_counts=nonlinear_counts)

    # TODO: Validate satval value?
#    satval = settings_det['saturation']*settings_det['nonlinear']
    # Order of the polynomials to be used when fitting the tilts.
    arcdet = (tcent[w]+0.5).astype(np.int)
    nsig = tnsig[w]

    # Determine the best lines to use to trace the tilts
    ncont = 15
    aduse = np.zeros(arcdet.size, dtype=np.bool)  # Which lines should be used to trace the tilts
    w = np.where(nsig >= tracethresh)
    aduse[w] = 1
    # Remove lines that are within ncont pixels
    nuse = np.sum(aduse)
    detuse = arcdet[aduse]
    idxuse = np.arange(arcdet.size)[aduse]
    olduse = aduse.copy()
    for s in range(nuse):
        w = np.where((np.abs(arcdet-detuse[s]) <= ncont) & (np.abs(arcdet-detuse[s]) >= 1.0))[0]
        for u in range(w.size):
            if nsig[w[u]] > nsig[olduse][s]:
                aduse[idxuse[s]] = False
                break
    # TODO Perhaps a more robust version of this code would only use the lines that were used in the wavelength solution. I guess
    # that would filter out these ghosts and it would also filter out blends for which the tracing will be less robust becuase
    # you are trace_fweighting a blended line?

    # Restricted to ID lines? [introduced to avoid LRIS ghosts]
    if idsonly:
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
                params, fail = utils.gauss_lsqfit(xfit, cc, 0.0)
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
            # check whether the yfit is offchip FW
            if np.size(yfit) == 0:
                offchip = True
                break
            elif len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            else:
                pass
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
                #if np.isfinite(centv) == False: # debugging
                #    from IPython import embed
                #    embed()
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
                params, fail = utils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
            xtfit[sz-k] = ordcen[arcdet[j], slitnum] - k
            ytfit[sz-k] = centv
            etfit[sz-k] = 0.02
            #apfit[sz-k] = params[0]
            if fail:
                mtfit[sz-k] = 1
            else:
                #from IPython import embed
                if np.isfinite(centv) == False: debugger.set_trace() #embed()
                pcen = int(0.5 + centv)
                mtfit[sz-k] = 0

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

def echelle_tilt(slf, msarc, det, settings_argflag, settings_spect,
                 pcadesc="PCA trace of the spectral tilts", maskval=-999999.9, doqa=True):
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
    arccen, maskslit, satmask = arc.get_censpec(slf._lordloc[det-1], slf._rordloc[det-1],
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
            null, tcoeff = utils.robust_polyfit(xtfit[wmask], ytfit[wmask],
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
            null, tempc = utils.robust_polyfit(centval[:, o][w], tiltang[:, o][w],
                                                 settings_argflag['trace']['slits']['tilts']['disporder'],
                                                 function=settings_argflag['trace']['slits']['function'], sigma=2.0,
                                                 minx=0.0, maxx=msarc.shape[0] - 1)
            tcoeff[:, o] = tempc
    # Sort which orders are masked
    maskord.sort()
    xv = np.arange(msarc.shape[0])
    tiltval = utils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'],
                               minx=0.0, maxx=msarc.shape[0] - 1).T
    ofit = settings_argflag['trace']['slits']['tilts']['params']
    lnpc = len(ofit) - 1
    if np.sum(1.0 - extrap_ord) > ofit[0] + 1:  # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral tilts")
        ordsnd = np.arange(norders) + 1.0
        xcen = xv[:, np.newaxis].repeat(norders, axis=1)
        fitted, outpar = pca.basis(xcen, tiltval, tcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=settings_argflag['trace']['slits']['function'])
        if doqa:
            #pcadesc = "Spectral Tilt PCA"
#            qa.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
            pca.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(norders)
        extrap_tilt, outpar = pca.extrapolate(outpar, orders, function=settings_argflag['trace']['slits']['function'])
        tilts = extrap_tilt
#        qa.pca_arctilt(slf, tiltang, centval, tilts)
        pca.pca_arctilt(slf, tiltang, centval, tilts)
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
            tcoeff = utils.func_fit(xtiltfit, ytiltfit, settings_argflag['trace']['slits']['function'],
                                      settings_argflag['trace']['slits']['tilts']['disporder'],
                                      minx=0.0, maxx=msarc.shape[0] - 1)
            tiltval = utils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'], minx=0.0,
                                       maxx=msarc.shape[0] - 1)
            tilts = tiltval[:, np.newaxis].repeat(tiltang.shape[1], axis=1)
        else:
            msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
            msgs.info("Assuming there is no tilt")
            tilts = np.zeros_like(slf._lordloc)

    # Generate tilts image
    tiltsimg = tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1],
                           settings_argflag['trace']['slits']['pad'], msarc.shape[1])

    return tiltsimg, satmask, outpar


'''
def multislit_tilt(msarc, lordloc, rordloc, pixlocn, pixcen, slitpix, det,
                   maskslits, tilt_settings, settings_spect, setup,
                   maskval=-999999.9, doqa=False, wv_calib=None):
    """ Determine the spectral tilt of each slit in a multislit image

    Parameters
    ----------
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
    arccen, arc_maskslit, _ = arc.get_censpec(lordloc, rordloc, pixlocn, msarc, det,
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
                             tracethresh=tilt_settings['tilts']['tracethresh'])
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        if doqa:
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

        """
        if tilt_settings['tilts']['method'].lower() == "interp":
            tilts = tilts_interp(ordcen, slit, all_tilts, polytilts, arcdet, aduse, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spline":
            tilts = tilts_spline(all_tilts, arcdet, aduse, polytilts, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spca":
            tilts = tilts_spca(msarc, polytilts, ordcen, slit, arcdet, aduse, rordloc, lordloc)
        elif tilt_settings['tilts']['method'].lower() == "pca":
            tilts = polytilts.copy()
        """

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
'''



'''
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
'''


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


def tilts_spline(all_tilts, arcdet, aduse, polytilts, msarc, use_mtilt=False, maskval=-999999.9, doqa=False):
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
    '''
    if doqa:
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
    coeff = utils.polyfit2d_general(xsbs, zsbs, ysbs, [5,5],
                                      minx=0., maxx=1., miny=0., maxy=1.,
                                      function='legendre')
    #tmpa = (msarc.shape[0]-1)*utils.polyval2d_general(coeff, np.arange(1490, 1499)/(msarc.shape[1]-1),
    #tmpb = np.polynomial.legendre.legval2d(xsbs.flatten(), zsbs.flatten(), coeff)
    tmpa = polyval2d_general(coeff, np.arange(1490, 1499)/(msarc.shape[1]-1),
                                     np.array([968.5, 969.34])/(msarc.shape[0]-1),
                                     minx=0., maxx=1., miny=0., maxy=1.,
                                     function='legendre')
    print(tmpa)
    debugger.set_trace()
    # Return
    return tilts


def polyval2d_general(c, x, y, function="polynomial", minx=None, maxx=None, miny=None, maxy=None):
    """Document me please. I think this evaluates on an image??"""

    if ('2d' in function):
        function = function[:-2]

    if function == "polynomial":
        xx, yy = np.meshgrid(x, y)
        return np.polynomial.polynomial.polyval2d(xx, yy, c)
    elif function in ["legendre", "chebyshev"]:
        # Scale x-direction
        xv = scale_minmax(x, minx=minx, maxx=maxx)
        # Scale y-direction
        yv = scale_minmax(y, minx=miny, maxx=maxy)
        xx, yy = np.meshgrid(xv, yv)
        if function == "legendre":
            return np.polynomial.legendre.legval2d(xx, yy, c)
        elif function == "chebyshev":
            return np.polynomial.chebyshev.chebval2d(xx, yy, c)
    else:
        msgs.error("Function {0:s} has not yet been implemented".format(function))
    return None


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

    #    qa.plot_orderfits(slf, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
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
    outroot = qa.set_qa_filename(setup, method, slit=slit)
    #
    npix, nord = ydata.shape
    pages, npp = qa.get_dimen(nord, maxp=maxp)
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



def plot_tiltres(setup, mtilt, ytilt, yfit, slit=None, outfile=None, show_QA=False, out_dir=None):
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
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

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



# This is useful when you do fits the other way. If you are fitting the tilts, this gives a flipped result
def eval_2d_at_tilts(tilts_spec, tilts_mask, shape, thismask, slit_cen, coeff2, func2d):

    # Compute Tilt model
    nspec = slit_cen.shape[0]
    nspat = tilts_spec.shape[0]
    nuse = tilts_spec.shape[1]

    spec_vec = np.arange(nspec)
    xnspecmin1 = float(nspec-1)
    tilts_img = np.zeros((nspec,nspat))
    tilts_img[thismask] = fit2tilts(shape, thismask, slit_cen, coeff2, func2d)
    piximg = xnspecmin1*tilts_img
    tilts_2dfit = np.zeros_like(tilts_spec)
    for iline in range(nuse):
        indx = np.where(tilts_mask[:, iline])[0]
        for ispat in indx:
            tilts_2dfit[ispat, iline] = np.interp(tilts_spec[0, iline], spec_vec, piximg[:, ispat])

    return tilts_2dfit


def fit_tilts_backup(trc_tilt_dict, slit_cen, spat_order=3, spec_order=4, maxdev = 0.2, sigrej = 3.0, func2d='legendre2d', doqa=True, setup = 'test',
              slit = 0, show_QA=False, out_dir=None, debug=False):
    """

    Parameters
    ----------
    trc_tilt_dict: dict
        Diciontary containing tilt info

    Optional Parameters
    -------------------
        slit:
        all_tilts:
        order:
        yorder:
        func2D:
        maskval:
        setup:
        doqa:
        show_QA:
        out_dir:

    Returns:

    """

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    fwhm = trc_tilt_dict['fwhm']
    maxdev_pix = maxdev*fwhm
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    nspat = trc_tilt_dict['nspat']
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    nuse = np.sum(use_tilt)
    tilts = trc_tilt_dict['tilts'][:,use_tilt]   # legendre polynomial fit
#   JFH Before we were fitting the fits. Now we fit the actual flux weighted centroided tilts.
    tilts_fit = trc_tilt_dict['tilts_fit'][:,use_tilt]   # legendre polynomial fit
    tilts_err = trc_tilt_dict['tilts_err'][:,use_tilt]   # flux weighted centroidding error
    tilts_dspat = trc_tilt_dict['tilts_dspat'][:,use_tilt] # spatial offset from the central trace
    tilts_spec_fit = trc_tilt_dict['tilts_spec_fit'][:,use_tilt] # line spectral pixel position from legendre fit evaluated at slit center
    tilts_mask = trc_tilt_dict['tilts_mask'][:,use_tilt] # Reflects if trace is on the slit
    tilts_mad = trc_tilt_dict['tilts_mad'][:,use_tilt] # Reflects if trace is on the slit

    # Do one last round of rejection here at the pixel level, i.e. we already rejected lines before
    #tot_mask = tilts_mask & (delta_tilt < maxdev) & (tilts_err < 900)
    tot_mask = tilts_mask & (tilts_err < 900)
    fitxy = [spat_order, spec_order]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2d))

    adderr = 0.03
    tilts_invvar = ((tilts_mad < 100.0) & (tilts_mad > 0.0))/(np.abs(tilts_mad) + adderr)**2
    # Previously we were using tilts_spec.flatten/xnspecmin1 in place of tilts_spec_fit
    # Burles was fitting the entire tilt and not the offset.
    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_dspat.flatten()/xnspatmin1, tilts_spec_fit.flatten()/xnspecmin1,
                                               fitxy, x2=tilts.flatten()/xnspecmin1, inmask = tot_mask.flatten(),
                                               invvar = xnspecmin1**2*tilts_invvar.flatten(),
                                               function=func2d, maxiter=100, lower=sigrej, upper=sigrej,
                                               maxdev=maxdev_pix/xnspecmin1,minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0,
                                               use_mad=False, sticky=False)

    fitmask = fitmask.reshape(tilts_dspat.shape)
    spec_fit2d_1 = xnspecmin1*utils.func_val(coeff2, tilts_dspat[tot_mask]/xnspatmin1, func2d, x2=tilts[tot_mask]/xnspecmin1,
                                               minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    spec_fit2d = np.zeros_like(tilts_dspat)
    spec_fit2d[tot_mask] = spec_fit2d_1
    # Residuals in pixels
    res_fit = tilts_spec_fit[fitmask] - spec_fit2d[fitmask]
    rms_fit = np.std(res_fit)
    msgs.info("Residuals: 2D Legendre Fit")
    msgs.info("RMS (pixels): {}".format(rms_fit))
    msgs.info("RMS/FWHM: {}".format(rms_fit/fwhm))

    # These are locations that were fit but were rejected
    rej_mask = tot_mask & np.invert(fitmask)

    # Tilt model
    tilts_img = fit2tilts((nspec, nspat), slit_cen, coeff2, func2d)
    piximg = xnspecmin1*tilts_img
    # Now construct the 2d model of each tilt  (above we modeled the line spectral position at slit center)
    spec_vec = np.arange(nspec)
    nlines = trc_tilt_dict['nlines']
    tilt_2dfit_all = np.zeros((nspat, nlines))
    for iline in range(nlines):
        for ispat in range(nspat): # Could make this faster by only looping over the parts on the actual slit
            tilt_2dfit_all[ispat, iline] = np.interp(trc_tilt_dict['tilts_spec_fit'][0, iline], spec_vec, piximg[:, ispat])

    trc_tilt_dict_out = copy.deepcopy(trc_tilt_dict)
    trc_tilt_dict_out['tilt_2dfit'] = tilt_2dfit_all
    tilts_2dfit = tilt_2dfit_all[:, use_tilt]

    # Actual 2D Model Tilt Residuals
    res = tilts[fitmask] - tilts_2dfit[fitmask]
    rms = np.std(res)
    msgs.info("Residuals: Actual 2D Tilt Residuals")
    msgs.info("RMS (pixels): {}".format(rms))
    msgs.info("RMS/FWHM: {}".format(rms/fwhm))

    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(fitmask == False), coeff2=coeff2, spec_order = spec_order, spat_order = spat_order,
                         minx = -1.0, maxx = 1.0, minx2 = 0.0, maxx2 = 1.0, func=func2d)

    # Now do some QA

    if doqa:
        # We could also be plotting the actual thing we fit below. Right now I'm trying with the 2d tilts themselves
        plot_tilt_2d(tilts_dspat, tilts, tilts_2dfit, tot_mask, rej_mask, spat_order, spec_order, rms, fwhm,
                     slit=slit, setup=setup, show_QA=show_QA, out_dir=out_dir)
        plot_tilt_spat(tilts_dspat, tilts, tilts_2dfit, tilts_spec_fit, tot_mask, rej_mask, spat_order, spec_order, rms, fwhm,
                       slit=slit, setup=setup, show_QA=show_QA, out_dir=out_dir)
        plot_tilt_spec(tilts_spec_fit, tilts, tilts_2dfit, tot_mask, rej_mask, rms, fwhm, slit=slit,
                       setup = setup, show_QA=show_QA, out_dir=out_dir)

    return tilts_img, tilt_fit_dict, trc_tilt_dict_out



def fit_tilts_xavier(trc_tilt_dict, spat_order=3, spec_order=4, maxdev = 0.2, sigrej = 3.0, func2d='legendre2d', doqa=True, setup = 'test',
              slit = 0, show_QA=False, out_dir=None, debug=False):
    """

    Parameters
    ----------
    trc_tilt_dict: dict
        Diciontary containing tilt info

    Optional Parameters
    -------------------
        slit:
        all_tilts:
        order:
        yorder:
        func2D:
        maskval:
        setup:
        doqa:
        show_QA:
        out_dir:

    Returns:

    """

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    fwhm = trc_tilt_dict['fwhm']
    maxdev_pix = maxdev*fwhm
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    nspat = trc_tilt_dict['nspat']
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    nuse = np.sum(use_tilt)
    tilts = trc_tilt_dict['tilts'][:,use_tilt]   # legendre polynomial fit
#   JFH Before we were fitting the fits. Now we fit the actual flux weighted centroided tilts.
#   tilts_fit = trc_tilt_dict['tilts_fit'][:,use_tilt]   # legendre polynomial fit
    tilts_err = trc_tilt_dict['tilts_err'][:,use_tilt]   # flux weighted centroidding error
    tilts_dspat = trc_tilt_dict['tilts_dspat'][:,use_tilt] # spatial offset from the central trace
    tilts_spec = trc_tilt_dict['tilts_spec'][:,use_tilt] # line_spec spectral pixel position
    tilts_mask = trc_tilt_dict['tilts_mask'][:,use_tilt] # Reflects if trace is on the slit

    # Do one last round of rejection here at the pixel level, i.e. we already rejected lines before
    #tot_mask = tilts_mask & (delta_tilt < maxdev) & (tilts_err < 900)
    tot_mask = tilts_mask & (tilts_err < 900)
    fitxy = [spat_order, spec_order]

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with a low order, 2D {:s}".format(func2d))

    # This is a 2-d polynomila, i.e. z(x,y), or following how robust_polyfit works y(x,x2)
    # x = spatial position along image, i.e. np.arange(nspat) for each line
    # y = spectral pixel where arc line was detected, i.e. line_spec replicated everywhere for each line
    # z = (spectral_line_trace - spectral_pos on extracted arc), i.e. tilt - y

    # We fit the function in this way i.e. with independent variable y being the actual tilt, becuase we are trying to
    # determine the mapping which brings the spectral line back to the same spectral_pos on extracted arc for which
    # we determined the wavelength solution

    # Xavier's way of fitting the offset
    # Fits are done in dimensionless coordinates to allow for different binnings between i.e. the arc and the science frame
    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_dspat.flatten()/xnspatmin1,
                                               (tilts.flatten() - tilts_spec.flatten())/xnspecmin1,
                                               fitxy,x2=tilts.flatten()/xnspecmin1, inmask = tot_mask.flatten(),
                                               function=func2d, maxiter=100, lower=sigrej, upper=sigrej,
                                               maxdev=maxdev_pix/xnspecmin1,minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0,
                                               use_mad=True, sticky=False)
    fitmask = fitmask.reshape(tilts_dspat.shape)
    delta_spec_fit1 = xnspecmin1*utils.func_val(coeff2, tilts_dspat[tot_mask]/xnspatmin1, func2d, x2=tilts[tot_mask]/xnspecmin1,
                                               minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    delta_spec_fit = np.zeros_like(tilts_dspat)
    delta_spec_fit[tot_mask] = delta_spec_fit1
    # Residuals in pixels
    res = (tilts[fitmask] - tilts_spec[fitmask]) - delta_spec_fit[fitmask]
    rms = np.std(res)
    msgs.info("RMS (pixels): {}".format(rms))
    msgs.info("RMS/FWHM: {}".format(rms/fwhm))

    # These are locations that were fit but were rejected
    rej_mask = tot_mask & np.invert(fitmask)

    # TODO: Make these plots a part of the standard QA and write to a file
    if debug:

        # QA1
        xmin = 1.1*tilts_dspat[tot_mask].min()
        xmax = 1.1*tilts_dspat[tot_mask].max()

        # Show the fit
        fig, ax = plt.subplots(figsize=(12, 18))
        # dummy mappable shows the spectral pixel
        #dummie_cax = ax.scatter(lines_spec, lines_spec, c=lines_spec, cmap=cmap)
        ax.cla()
        ax.plot(tilts_dspat[tot_mask], tilts[tot_mask], color='black', linestyle=' ', mfc ='None', marker='o',
                markersize = 9.0, markeredgewidth=1.0,zorder=4, label='Good Tilt')
        ax.plot(tilts_dspat[rej_mask], tilts[rej_mask], color ='red',linestyle=' ', mfc = 'None', marker='o',
                markersize=9.0, markeredgewidth=2.0, zorder=5, label='Rejected')
        ax.plot(tilts_dspat[tot_mask], tilts_spec[tot_mask] + delta_spec_fit[tot_mask], color='black', linestyle=' ', marker='o',
                markersize = 2.0,markeredgewidth=1.0,zorder=1)

        ax.set_xlim((xmin,xmax))
        ax.set_xlabel('Spatial Offset from Central Trace (pixels)', fontsize=15)
        ax.set_ylabel('Spectral Pixel',fontsize=15)
        ax.legend()
        ax.set_title('Tilts vs Fit (spat_order, spec_order)=({:d},{:d}) for slit={:d}: RMS = {:5.3f}, '
                     'RMS/FWHM={:5.3f}'.format(spat_order,spec_order,slit,rms, rms/fwhm),fontsize=15)
        plt.show()

        # QA2
        # Show the fit residuals as a function of spatial position
        line_indx = np.outer(np.ones(nspat), np.arange(nuse))
        lines_spec = tilts_spec[0,:]
        cmap = mpl.cm.get_cmap('coolwarm', nuse)

        fig, ax = plt.subplots(figsize=(14, 12))
        # dummy mappable shows the spectral pixel
        dummie_cax = ax.scatter(lines_spec, lines_spec, c=lines_spec, cmap=cmap)
        ax.cla()

        for iline in range(nuse):
            iall = (line_indx == iline) & tot_mask
            irej = (line_indx == iline) & tot_mask & rej_mask
            this_color = cmap(iline)
            # plot the residuals
            ax.plot(tilts_dspat[iall], (tilts[iall] - tilts_spec[iall]) - delta_spec_fit[iall], color=this_color,
                    linestyle='-', linewidth=3.0, marker='None', alpha=0.5)
            ax.plot(tilts_dspat[irej],(tilts[irej] - tilts_spec[irej]) - delta_spec_fit[irej],linestyle=' ',
                    marker='o', color = 'limegreen', mfc='limegreen', markersize=5.0)

        ax.hlines(0.0, xmin, xmax, linestyle='--', linewidth=2.0, color='k', zorder=10)

        legend_elements = [Line2D([0], [0], color='cornflowerblue', linestyle='-', linewidth=3.0, label='residual'),
                           Line2D([0], [0], color='limegreen', linestyle=' ', marker='o', mfc='limegreen', markersize=7.0, label='rejected')]

        ax.set_xlim((xmin,xmax))
        ax.set_xlabel('Spatial Offset from Central Trace (pixels)')
        ax.set_ylabel('Arc Line Tilt Residual (pixels)')
        ax.legend(handles=legend_elements)
        cb = fig.colorbar(dummie_cax, ticks=lines_spec)
        cb.set_label('Spectral Pixel')
        plt.show()

    # TODO Add QA where we overlay the final model of the tilts on the image using ginga!

    # QA
    if doqa:
        plot_tiltres_xavier(setup, tilts[tot_mask], tilts_spec[tot_mask], delta_spec_fit[tot_mask], fitmask[tot_mask], fwhm, slit=slit, show_QA=show_QA, out_dir=out_dir)

    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(fitmask == False), coeff2=coeff2, spec_order = spec_order, spat_order = spat_order,
                         minx = -1.0, maxx = 1.0, minx2 = 0.0, maxx2 = 1.0, func=func2d)

    return tilt_fit_dict



def fit_tilts_joe(trc_tilt_dict, spat_order=3, spec_order=4, maxdev = 0.2, sigrej = 3.0, func2d='legendre2d', doqa=True, setup = 'test',
              slit = 0, show_QA=False, out_dir=None, debug=False):
    """

    Parameters
    ----------
    trc_tilt_dict: dict
        Diciontary containing tilt info

    Optional Parameters
    -------------------
        slit:
        all_tilts:
        order:
        yorder:
        func2D:
        maskval:
        setup:
        doqa:
        show_QA:
        out_dir:

    Returns:

    """

    import matplotlib as mpl
    from matplotlib.lines import Line2D

    nspec = trc_tilt_dict['nspec']
    nspat = trc_tilt_dict['nspat']
    fwhm = trc_tilt_dict['fwhm']
    maxdev_pix = maxdev*fwhm
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    nspat = trc_tilt_dict['nspat']
    use_tilt = trc_tilt_dict['use_tilt']                 # mask for good/bad tilts, based on aggregate fit, frac good pixels
    nuse = np.sum(use_tilt)
    tilts = trc_tilt_dict['tilts'][:,use_tilt]   # legendre polynomial fit
#   JFH Before we were fitting the fits. Now we fit the actual flux weighted centroided tilts.
    tilts_fit = trc_tilt_dict['tilts_fit'][:,use_tilt]   # legendre polynomial fit
    tilts_err = trc_tilt_dict['tilts_err'][:,use_tilt]   # flux weighted centroidding error
    tilts_spat = trc_tilt_dict['tilts_spat'][:,use_tilt] # spatial position of the tilt on the image
    tilts_dspat = trc_tilt_dict['tilts_dspat'][:,use_tilt] # spatial offset from the central trace
    tilts_spec = trc_tilt_dict['tilts_spec'][:,use_tilt] # line_spec spectral pixel position
    tilts_mask = trc_tilt_dict['tilts_mask'][:,use_tilt] # Reflects if trace is on the slit

    # Do one last round of rejection here at the pixel level, i.e. we already rejected lines before
    #tot_mask = tilts_mask & (delta_tilt < maxdev) & (tilts_err < 900)
    tot_mask = tilts_mask & (tilts_err < 900)

    # Fit the inverted model with a 2D polynomial
    msgs.info("Fitting tilts with polynomial warping")
    # Burles used the midpoint of the tilt fit instead of tilts_spec below
    zfit = np.zeros_like(tilts_spec)
    for iline in range(nuse):
        imask = tilts_mask[:,iline]
        zfit[:,iline] = np.full(nspat,np.interp(0.0,tilts_dspat[imask,iline],tilts_fit[imask,iline]))

    tilts_spec=zfit.copy()

    from skimage import transform
    from skimage.measure import ransac

    from IPython import embed
    embed()
    src = np.stack((tilts_spat[tot_mask], tilts_spec[tot_mask]),axis=1)
    dst = np.stack((tilts_spat[tot_mask], tilts[tot_mask]),axis=1)

    poly_order = 3
    #poly_trans = transform.PolynomialTransform()
    #out = poly_trans.estimate(dst,src,order=poly_order)
    #residuals = poly_trans.residuals(dst, src)
    #inliers=np.ones_like(residuals,dtype=bool)


    from skimage.transform import PolynomialTransform

    class PolyTF3(PolynomialTransform):
        def estimate(self, src, dst):
            return super().estimate(src, dst, order=poly_order)

    poly_trans, inliers = ransac((dst, src), PolyTF3, min_samples=20, residual_threshold=0.1,max_trials=1000)
    residuals = poly_trans.residuals(dst[inliers, :], src[inliers, :])
    rms = np.std(residuals)
    msgs.info("RMS (pixels): {}".format(rms))
    msgs.info("RMS/FWHM: {}".format(rms / fwhm))

    # This will generate the tilts image
    #piximg_src = np.outer(np.arange(nspec),np.ones(nspat))
    #piximg_dst = transform.warp(piximg_src, poly_trans, mode='edge',cval=float('nan'))
    #imask = (piximg_dst  > 1077) & (piximg_dst < 1078)
    #ginga.show_image(piximg*imask)


    tilt_fit_dict = dict(nspec = nspec, nspat = nspat, ngood_lines=np.sum(use_tilt), npix_fit = np.sum(tot_mask),
                         npix_rej = np.sum(inliers == False), coeff2=poly_trans.params, poly_order=poly_order,
                         spec_order = spec_order, spat_order = spat_order,
                         minx = -1.0, maxx = 1.0, minx2 = 0.0, maxx2 = 1.0, func=func2d)

    return tilt_fit_dict


def fit2tilts_joe(shape, slit_cen, coeff2, func):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    Returns
    -------
    piximg: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """


    from skimage import transform
    import scipy
    # Compute the tilts image
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    poly_trans = transform.PolynomialTransform(params=coeff2)
    # This will generate the tilts image centered at the midpoit of the image
    piximg_src = np.outer(np.arange(nspec),np.ones(nspat))
    piximg_dst = transform.warp(piximg_src, poly_trans, mode='edge',cval=float('nan'))
    #imask = (piximg  > 1024) & (piximg < 1026)
    #ginga.show_image(piximg*imask)

    #spec_img = np.outer(np.arange(nspec),np.ones(nspat))
    #spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
    #slit_cen_img = np.outer(slit_cen, np.ones(nspat))   # center of the slit replicated spatially

    tilts = np.fmax(np.fmin(piximg_dst/xnspecmin1, 1.2),-0.2)

    # Normalized spatial offset image (from central trace)
    #dspat_img_nrm = (spat_img - slit_cen_img)
    # Normalized spatial offset image offset by nspat//2 to the frame of the piximg
    #dspat_img_nrm_cen = (spat_img - slit_cen_img) + nspat//2

    # For pixels with completely bad profile values, interpolate from trace.
    #piximg_spline=scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, piximg)
    #piximg = np.zeros((nspec,nspat))
    #piximg_test = piximg_spline(spec_img.flatten(), dspat_img_nrm_cen.flatten(), grid=False)
    #piximg_final = piximg_test.reshape(nspec, nspat)
    #tilts = piximg_final/xnspecmin1
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting

    return tilts


def fit2tilts_xavier(shape, slit_cen, coeff2, func):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    Returns
    -------
    piximg: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """

    # Compute the tilts image
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)

    spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
    slit_cen_img = np.outer(slit_cen, np.ones(nspat))   # center of the slit replicated spatially
    # Normalized spatial offset image (from central trace)
    dspat_img_nrm = (spat_img - slit_cen_img)/xnspatmin1
    # normalized spec image
    spec_img_nrm = np.outer(np.arange(nspec), np.ones(nspat))/xnspecmin1

    delta_spec_fit_nrm = utils.func_val(coeff2, dspat_img_nrm, func, x2=spec_img_nrm, minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    tilts = (spec_img_nrm - delta_spec_fit_nrm)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)

    return tilts


def plot_tiltres_xavier(setup, mtilt, ytilt, yfit, fitmask, fwhm, slit=None, outfile=None, show_QA=False, out_dir=None):
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
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

    # Setup
    plt.figure(figsize=(8, 4))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (mtilt - ytilt) - yfit
    xmin = 0.95*np.asarray([ytilt,mtilt]).min()
    xmax = 1.05*np.asarray([ytilt,mtilt]).max()
    ax.hlines(0.0, xmin, xmax,linestyle='--', color='green')
    ax.plot(ytilt[fitmask], (res[fitmask]), 'ko', mfc='k', markersize=5.0, label='Good Points')
    ax.plot(ytilt[~fitmask],(res[~fitmask]), 'ro', mfc='r', markersize=5.0, label='Rejected Points')
    #ax.scatter(mtilt, res)
    rms = np.std(res[fitmask])
    ax.text(0.90, 0.90, 'Slit {:d}:  RMS (pixels) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    ax.text(0.90, 0.80, ' Slit {:d}:  RMS/FWHM = {:0.5f}'.format(slit, rms/fwhm),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    # Label
    ax.set_xlabel('Spectral Pixel')
    ax.set_ylabel('RMS (pixels)')
    ax.set_title('RMS of Each Arc Line Traced')
    ax.set_xlim((xmin,xmax))
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if show_QA:
        plt.show()

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    plt.close()
    plt.rcdefaults()

    return

# Joe version
def plot_tiltres_joe(setup, mtilt, ytilt, yfit, fitmask, fwhm, slit=None, outfile=None, show_QA=False, out_dir=None):
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
        outfile = qa.set_qa_filename(setup, method, slit=slit, out_dir=out_dir)

    # Setup
    plt.figure(figsize=(10, 6))
    plt.clf()
    ax = plt.gca()

    # Scatter plot
    res = (mtilt - yfit)
    xmin = 0.95*np.asarray([ytilt,mtilt]).min()
    xmax = 1.05*np.asarray([ytilt,mtilt]).max()
    ax.hlines(0.0, xmin, xmax,linestyle='--', color='green')
    ax.plot(ytilt[fitmask], (res[fitmask]), 'ko', mfc='k', markersize=5.0, label='Good Points')
    ax.plot(ytilt[~fitmask],(res[~fitmask]), 'ro', mfc='r', markersize=5.0, label='Rejected Points')
    #ax.scatter(mtilt, res)
    rms = np.std(res[fitmask])
    ax.text(0.90, 0.90, 'Slit {:d}:  RMS (pixels) = {:0.5f}'.format(slit, rms),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    ax.text(0.90, 0.80, ' Slit {:d}:  RMS/FWHM = {:0.5f}'.format(slit, rms/fwhm),
            transform=ax.transAxes, size='large', ha='right', color='black',fontsize=16)
    # Label
    ax.set_xlabel('Spectral Pixel')
    ax.set_ylabel('RMS (pixels)')
    ax.set_title('RMS of Each Arc Line Traced')
    ax.set_xlim((xmin,xmax))
    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if show_QA:
        plt.show()

    if outfile is not None:
        plt.savefig(outfile, dpi=400)

    plt.close()
    plt.rcdefaults()

    return


# This did not work
def fit_tilts_rej(tilts_dspat, tilts_spec_fit, tilts, tilts_invvar, tot_mask, slit_cen, spat_order, spec_order,
                  maxdev = 0.2, maxrej=None, sigrej = 3.0, maxiter = 20, func2d='legendre2d'):

    from pypeit.core import pydl

    nspat = tilts_dspat.shape[0]
    nuse = tilts_dspat.shape[1]
    nspec = slit_cen.shape[0]
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    fitxy = [spat_order, spec_order]

    lower=sigrej
    upper=sigrej

    msgs.info('                                  Iteratively Fitting Tilts' + msgs.newline() +
              '                                  -------------------------' + msgs.newline() +
              '******************************  Iter  # rejected   RMS (pixels)  ******************************' + msgs.newline() +
              '                                ----  ----------  ------------- ')
    # Reject bad pixels
    iIter = 0
    qdone = False
    thismask = np.copy(tot_mask.flatten())
    while (not qdone) and (iIter < maxiter):
        # Do a fit with no rejection
        fitmask, coeff2 = utils.robust_polyfit_djs(tilts_dspat.flatten()/xnspatmin1,
                                                   tilts_spec_fit.flatten()/xnspecmin1,
                                                   fitxy, x2=tilts.flatten()/xnspecmin1,
                                                   inmask=(thismask*tot_mask.flatten()),
                                                   invvar=xnspecmin1**2*tilts_invvar.flatten(),
                                                   function=func2d, maxiter=0,
                                                   minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)

        # Compute Tilt model
        tilts_2dfit = eval_2d_at_tilts(tilts_spec_fit, tot_mask, (nspec, nspat), slit_cen, coeff2, func2d)

        # Perform rejection
        thismask, qdone = pydl.djs_reject(tilts.flatten(), tilts_2dfit.flatten(), outmask=thismask,inmask=tot_mask.flatten(),
                                          lower=lower,upper=upper,maxdev=maxdev,maxrej=maxrej,
                                          use_mad=True,sticky=False)
#        invvar = tilts_invvar.flatten(),
        nrej = np.sum(tot_mask.flatten() & np.invert(thismask))
        rms = np.std(tilts.flat[thismask] - tilts_2dfit.flat[thismask])
        msgs.info('                                 {:d}        {:d}  {:6.3f}'.format(iIter, nrej,rms))
        iIter += 1

    if (iIter == maxiter) & (maxiter != 0):
        msgs.warn('Maximum number of iterations maxiter={:}'.format(maxiter) + ' reached in fit_tilts_rej')
    outmask = np.copy(thismask)
    if np.sum(outmask) == 0:
        msgs.warn('All points were rejected!!! The fits will be zero everywhere.')

    # Do the final fit
    fitmask, coeff2 = utils.robust_polyfit_djs(tilts_dspat.flatten()/xnspatmin1,
                                               tilts_spec_fit.flatten()/xnspecmin1,
                                               fitxy, x2=tilts.flatten()/xnspecmin1,
                                               inmask=(outmask*tot_mask.flatten()),
                                               invvar=xnspecmin1 ** 2 * tilts_invvar.flatten(),
                                               function=func2d, maxiter=0,
                                               minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)



    outmask = outmask.reshape((nspat,nuse))
    return outmask, coeff2


# This was too slow
def fit_tilts_func(coeff2, tilts_dspat, tilts_spec_fit, tilts, tilts_invvar, tot_mask, inmask, slit_cen, spat_order, spec_order):

    """

    Parameters
    ----------
    trc_tilt_dict: dict
        Diciontary containing tilt info

    Optional Parameters
    -------------------
        slit:
        all_tilts:
        order:
        yorder:
        func2D:
        maskval:
        setup:
        doqa:
        show_QA:
        out_dir:

    Returns:

    """

    #print('Call:', coeff2)
    coeff2_arr = coeff2.reshape((spat_order + 1, spec_order+1))
    nspat = tilts_dspat.shape[0]
    nuse = tilts_dspat.shape[1]
    nspec = slit_cen.shape[0]
    spec_vec = np.arange(nspec)
    xnspecmin1 = float(nspec-1)
    # Tilt model
    tilts_img = fit2tilts((nspec, nspat), slit_cen, coeff2_arr, 'legendre2d')
    piximg = xnspecmin1*tilts_img
    # Now construct the 2d model of each tilt  (above we modeled the line spectral position at slit center)
    chi2 = 0.0
    for iline in range(nuse):
        indx = np.where(tot_mask[:,iline] & inmask[:,iline])[0]
        for ispat in indx:
            tilt_2dfit = np.interp(tilts_spec_fit[0, iline], spec_vec, piximg[:, ispat])
            chi2 += ((tilts[ispat,iline] - tilt_2dfit)*np.sqrt(tilts_invvar[ispat,iline]))**2

    return chi2




# JFH This is now defunct
def gridtilts(shape, thismask, slit_cen, coeff2, func2d, spec_order, spat_order, pad_spec=30, pad_spat = 5, method='interp'):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    Returns
    -------
    piximg: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """

    # Compute the tilts image
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    # JFH This histogram method is not preferred, since it basically does NGP. It is however super fast, so for big images
    # it is useful to have it
    if 'hist2d' in method:
        oversamp_spec=5
        oversamp_spat=3
        spec_ind, spat_ind = np.where(thismask)
        min_spec = spec_ind.min() - pad_spec
        max_spec = spec_ind.max() + pad_spec
        num_spec = max_spec - min_spec + 1

        min_spat = spat_ind.min() - pad_spat
        max_spat = spat_ind.max() + pad_spat
        num_spat = max_spat - min_spat + 1

        spec_lin = np.linspace(min_spec,max_spec,num = int(np.round(num_spec*oversamp_spec)))
        spat_lin = np.linspace(min_spat,max_spat,num = int(np.round(num_spat*oversamp_spat)))
        spat_img, spec_img = np.meshgrid(spat_lin, spec_lin)
        # Normalized spatial offset image (from central trace)
        slit_cen_lin = (scipy.interpolate.interp1d(np.arange(nspec),slit_cen,bounds_error=False,fill_value='extrapolate'))(spec_lin)
        slit_cen_img = np.outer(slit_cen_lin, np.ones(spat_img.shape[1]))  # center of the slit replicated spatially
        dspat_img_nrm = (spat_img - slit_cen_img)/xnspatmin1
        spec_img_nrm = spec_img/xnspecmin1

        # normalized spec image
        tracepix = spec_img + xnspecmin1*utils.func_val(coeff2, spec_img_nrm, func2d, x2=dspat_img_nrm,
                                                              minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
        norm_img, spec_edges, spat_edges = np.histogram2d(tracepix.flatten(), spat_img.flatten(),
                                                          bins=[np.arange(nspec+1), np.arange(nspat+1)], density=False)
        weigh_img, spec_edges, spat_edges = np.histogram2d(tracepix.flatten(), spat_img.flatten(),
                                                           bins=[np.arange(nspec+1), np.arange(nspat+1)],
                                                           weights = spec_img.flatten(),density=False)
        piximg =(norm_img > 0.0)*weigh_img/(norm_img + (norm_img == 0.0))
        inmask = thismask & (norm_img > 0) & (piximg/xnspecmin1 > -0.2) & (piximg/xnspecmin1 < 1.2)
    # This is the defulat method although scipy.interpolate.griddata is a bit slow
    elif 'interp' in method:
        spec_vec_pad = np.arange(-pad_spec,nspec+pad_spec)
        spat_vec_pad = np.arange(-pad_spat,nspat+pad_spat)
        spat_img, spec_img = np.meshgrid(spat_vec, spec_vec)
        spat_img_pad, spec_img_pad = np.meshgrid(np.arange(-pad_spat,nspat+pad_spat),np.arange(-pad_spec,nspec+pad_spec))
        slit_cen_pad = (scipy.interpolate.interp1d(spec_vec,slit_cen,bounds_error=False,fill_value='extrapolate'))(spec_vec_pad)
        thismask_pad = np.zeros_like(spec_img_pad,dtype=bool)
        ind_spec, ind_spat = np.where(thismask)
        slit_cen_img_pad= np.outer(slit_cen_pad, np.ones(nspat + 2*pad_spat))  # center of the slit replicated spatially
        # Normalized spatial offset image (from central trace)
        dspat_img_nrm = (spat_img_pad - slit_cen_img_pad)/xnspatmin1
        # normalized spec image
        spec_img_nrm = spec_img_pad/xnspecmin1
        # Embed the old thismask in the new larger padded thismask
        thismask_pad[ind_spec + pad_spec,ind_spat + pad_spat] = thismask[ind_spec,ind_spat]
        # Now grow the thismask_pad
        kernel = np.ones((2*pad_spec, 2*pad_spat))/float(4*pad_spec*pad_spat)
        thismask_grow = scipy.ndimage.convolve(thismask_pad.astype(float), kernel, mode='nearest') > 0.0
        # Evaluate the tilts on the padded image grid
        tracepix = spec_img_pad[thismask_grow] + xnspecmin1*utils.func_val(coeff2, spec_img_nrm[thismask_grow], func2d, x2=dspat_img_nrm[thismask_grow],
                                                              minx=0.0, maxx=1.0, minx2=-1.0, maxx2=1.0)
        ## TESTING STARTS
        """
        ikeep = np.isfinite(tracepix)
        sigma = np.full_like(spec_img_pad[thismask_grow], 10.0)/xnspecmin1
        fitxy = [spec_order, spat_order]
        fitmask, coeff2_tilts = utils.robust_polyfit_djs(tracepix/xnspecmin1, spec_img_pad[thismask_grow]/xnspecmin1,
                                                         fitxy, x2=spat_img_pad[thismask_grow]/xnspatmin1,
                                                         sigma=sigma,
                                                         upper=5.0, lower=5.0, maxdev=10.0/xnspecmin1,
                                                         inmask=ikeep, function=func2d, maxiter=20,
                                                         minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0, use_mad=False)
        ## TESTING ENDS
        # values(points) \equiv spec_pos(tilt,spat_pos) which is the piximg that we want to create via griddata interpolation
        """
        ikeep = np.isfinite(tracepix)
        points = np.stack((tracepix[ikeep], spat_img_pad[thismask_grow][ikeep]), axis=1)
        values =spec_img_pad[thismask_grow][ikeep]
        piximg = scipy.interpolate.griddata(points, values, (spec_img, spat_img), method='cubic')
        inmask = thismask & np.isfinite(piximg) & (piximg/xnspecmin1 > -0.2) & (piximg/xnspecmin1 < 1.2)

    # Now simply do a 2d polynomial fit with just rejection of crazy behavior, i.e. 10 pixels
    fitxy = [spec_order, spat_order]
    sigma = np.full_like(spec_img,10.0)/xnspecmin1
    fitmask, coeff2_tilts = utils.robust_polyfit_djs(spec_img.flatten()/xnspecmin1, piximg.flatten()/xnspecmin1,
                                                     fitxy, x2=spat_img.flatten()/xnspatmin1, sigma = sigma.flatten(),
                                                     upper=5.0, lower=5.0, maxdev = 10.0/xnspecmin1,
                                                     inmask=inmask.flatten(), function=func2d, maxiter=20,
                                                     minx=0.0, maxx=1.0, minx2=0.0,maxx2=1.0,use_mad=False)
    irej = np.invert(fitmask) & inmask.flatten()
    msgs.info('Rejected {:d}/{:d} pixels in final tilts image after gridding'.format(np.sum(irej),np.sum(inmask)))
    # normalized tilts image

    tilts = utils.func_val(coeff2_tilts, spec_img/xnspecmin1, func2d, x2=spat_img/xnspatmin1,minx=0.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting

    return coeff2_tilts, tilts



def fit2tilts_backup(shape, slit_cen, coeff2, func):
    """

    Parameters
    ----------
    tilt_fit_dict: dict
        Tilt fit dictioary produced by fit_tilts

    Returns
    -------
    piximg: ndarray, float
       Image indicating how spectral pixel locations move across the image. This output is used in the pipeline.
    """

    # Compute the tilts image
    nspec, nspat = shape
    xnspecmin1 = float(nspec-1)
    xnspatmin1 = float(nspat-1)

    spat_img = np.outer(np.ones(nspec), np.arange(nspat)) # spatial position everywhere along image
    slit_cen_img = np.outer(slit_cen, np.ones(nspat))   # center of the slit replicated spatially
    # Normalized spatial offset image (from central trace)
    dspat_img_nrm = (spat_img - slit_cen_img)/xnspatmin1
    # normalized spec image
    spec_img_nrm = np.outer(np.arange(nspec), np.ones(nspat))/xnspecmin1

    tilts = utils.func_val(coeff2, dspat_img_nrm, func, x2=spec_img_nrm, minx=-1.0, maxx=1.0, minx2=0.0, maxx2=1.0)
    # Added this to ensure that tilts are never crazy values due to extrapolation of fits which can break
    # wavelength solution fitting
    tilts = np.fmax(np.fmin(tilts, 1.2),-0.2)

    return tilts


# This should be deprecated, since there is another version in core.artraceslits
#def trace_fweight(fimage, xinit, ltrace=None, rtraceinvvar=None, radius=3.,maskval=999999.9):
def trace_fweight_deprecated(fimage, xinit, ltrace=None, rtraceinvvar=None, radius=3., maskval=999999.9):
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
    maskval : float, optional
    """

    # Init
    nx = fimage.shape[1]
    ny = fimage.shape[0]
    ncen = len(xinit)
    xnew = copy.deepcopy(xinit)
    xerr = np.zeros(ncen) + maskval

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
        xerr[bad] = maskval

    # Return
    return xnew, xerr


def echelle_tilt(slf, msarc, det, settings_argflag, settings_spect,
                 pcadesc="PCA trace of the spectral tilts", maskval=-999999.9, doqa=True):
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
    arccen, maskslit, satmask = arc.get_censpec(slf._lordloc[det-1], slf._rordloc[det-1],
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
            null, tcoeff = utils.robust_polyfit(xtfit[wmask], ytfit[wmask],
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
            null, tempc = utils.robust_polyfit(centval[:, o][w], tiltang[:, o][w],
                                                 settings_argflag['trace']['slits']['tilts']['disporder'],
                                                 function=settings_argflag['trace']['slits']['function'], sigma=2.0,
                                                 minv=0.0, maxv=msarc.shape[0] - 1)
            tcoeff[:, o] = tempc
    # Sort which orders are masked
    maskord.sort()
    xv = np.arange(msarc.shape[0])
    tiltval = utils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'],
                               minv=0.0, maxv=msarc.shape[0] - 1).T
    ofit = settings_argflag['trace']['slits']['tilts']['params']
    lnpc = len(ofit) - 1
    if np.sum(1.0 - extrap_ord) > ofit[0] + 1:  # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral tilts")
        ordsnd = np.arange(norders) + 1.0
        xcen = xv[:, np.newaxis].repeat(norders, axis=1)
        fitted, outpar = pca.basis(xcen, tiltval, tcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=settings_argflag['trace']['slits']['function'])
        if doqa:
            #pcadesc = "Spectral Tilt PCA"
#            qa.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
            pca.pca_plot(slf, outpar, ofit, 'Arc', pcadesc=pcadesc, addOne=False)
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(norders)
        extrap_tilt, outpar = pca.extrapolate(outpar, orders, function=settings_argflag['trace']['slits']['function'])
        tilts = extrap_tilt
#        qa.pca_arctilt(slf, tiltang, centval, tilts)
        pca.pca_arctilt(slf, tiltang, centval, tilts)
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
            tcoeff = utils.func_fit(xtiltfit, ytiltfit, settings_argflag['trace']['slits']['function'],
                                      settings_argflag['trace']['slits']['tilts']['disporder'],
                                      minv=0.0, maxv=msarc.shape[0] - 1)
            tiltval = utils.func_val(tcoeff, xv, settings_argflag['trace']['slits']['function'], minv=0.0,
                                       maxv=msarc.shape[0] - 1)
            tilts = tiltval[:, np.newaxis].repeat(tiltang.shape[1], axis=1)
        else:
            msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
            msgs.info("Assuming there is no tilt")
            tilts = np.zeros_like(slf._lordloc)

    # Generate tilts image
    tiltsimg = tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1],
                           settings_argflag['trace']['slits']['pad'], msarc.shape[1])

    return tiltsimg, satmask, outpar


def pca_arctilt(slf, tiltang, centval, tilts, maxp=25, maskval=-999999.9):
    """ Generate a QA plot for the blaze function fit to each slit

    Parameters
    ----------
    slf : class
      Science Exposure class
    tiltang : ndarray
      (m x n) 2D array containing the measured tilts (m) for each slit (n)
    centval : ndarray
      (m x n) 2D array containing the pixel location (in the spectral direction) of the measured tilts (m) for each slit (n)
    tilts : ndarray
      (m x n) 2D array containing the model tilts (m) for each slit (n)
    maxp : int, (optional)
      Maximum number of panels per page
    maskval : float, (optional)
      Value used in arrays to indicate a masked value
    """
    plt.rcdefaults()
    plt.rcParams['font.family'] = 'times new roman'

    method = inspect.stack()[0][3]
    outroot = arqa.set_qa_filename(slf.setup, method)
    #
    npc = tiltang.shape[1]
    pages, npp = arqa.get_dimen(npc, maxp=maxp)
    x0 = np.arange(tilts.shape[0])
    # First calculate the min and max values for the plotting axes, to make sure they are all the same
    w = np.where(tiltang != maskval)
    medv = np.median(tiltang[w])
    madv = 1.4826 * np.median(np.abs(medv - tiltang[w]))
    ymin, ymax = medv - 3.0 * madv, medv + 3.0 * madv
    ymin = min(ymin, np.min(tilts))
    ymax = max(ymax, np.max(tilts))
    # Check that ymin and ymax are set, if not, return without plotting
    if ymin is None or ymax is None:
        msgs.warn("Arc tilt fits were not plotted")
        return
    # Generate the plots
    ndone = 0
    for i in range(len(pages)):
        f, axes = plt.subplots(pages[i][1], pages[i][0])
        ipx, ipy = 0, 0
        for j in range(npp[i]):
            if pages[i][1] == 1:
                ind = (ipx,)
            elif pages[i][0] == 1:
                ind = (ipy,)
            else:
                ind = (ipy, ipx)
            w = np.where(tiltang[:, ndone] != maskval)
            if np.size(w[0]) != 0:
                axes[ind].plot(centval[:, ndone][w], tiltang[:, ndone][w], 'bx')
                axes[ind].plot(x0, tilts[:, ndone], 'r-')
            axes[ind].axis([0, tilts.shape[0] - 1, ymin, ymax])
            axes[ind].set_title("Slit {0:d}".format(1 + ndone))
            axes[ind].tick_params(labelsize=8)
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
            ndone += 1
        # Delete the unnecessary axes
        for j in range(npp[i], axes.size):
            if pages[i][1] == 1:
                ind = (ipx,)
            elif pages[i][0] == 1:
                ind = (ipy,)
            else:
                ind = (ipy, ipx)
            f.delaxes(axes[ind])
            ipx += 1
            if ipx == pages[i][0]:
                ipx = 0
                ipy += 1
        # Save the figure
        if pages[i][1] == 1 or pages[i][0] == 1:
            ypngsiz = 11.0 / axes.size
        else:
            ypngsiz = 11.0 * axes.shape[0] / axes.shape[1]
        f.set_size_inches(11.0, ypngsiz)
        f.tight_layout()
        outfile = outroot + '{:02d}.png'.format(i)
        plt.savefig(outfile, dpi=200)
        plt.close()
        # f.clf()
        # del f

    plt.rcdefaults()

    return


'''
def multislit_tilt(msarc, lordloc, rordloc, pixlocn, pixcen, slitpix, det,
                   maskslits, tilt_settings, settings_spect, setup,
                   maskval=-999999.9, doqa=False, wv_calib=None):
    """ Determine the spectral tilt of each slit in a multislit image

    Parameters
    ----------
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
    arccen, arc_maskslit, _ = arc.get_censpec(lordloc, rordloc, pixlocn, msarc, det,
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
                             tracethresh=tilt_settings['tilts']['tracethresh'])
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        if doqa:
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

        """
        if tilt_settings['tilts']['method'].lower() == "interp":
            tilts = tilts_interp(ordcen, slit, all_tilts, polytilts, arcdet, aduse, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spline":
            tilts = tilts_spline(all_tilts, arcdet, aduse, polytilts, msarc)
        elif tilt_settings['tilts']['method'].lower() == "spca":
            tilts = tilts_spca(msarc, polytilts, ordcen, slit, arcdet, aduse, rordloc, lordloc)
        elif tilt_settings['tilts']['method'].lower() == "pca":
            tilts = polytilts.copy()
        """

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
'''


#TODO this should be deleted, but I'm saving some it for now becasue I'm still testing.
'''


        # JFH Code block starts here
        ########
        # def get_tilts(npca = 2, fwhm=4.0, ncoeff=5, maxdev_tracefit=0.1,percentile_reject=0.10, max_badpix_frac=0.20, tcrude_maxerr=1.0,
        # tcrude_maxshift=3.0, tcrude_maxshift0=3.0, tcrude_nave=5,)

        # DEBUGGIN
        slit = 4

        from pypeit.core import pixels
        from pypeit.core import extract


        nspat = self.msarc.shape[1]
        nspec = self.msarc.shape[0]
        arcimg = self.msarc
        arc_spec = self.arccen[:, slit]
        slit_left = self.tslits_dict['lcen'][:,slit].copy()
        slit_righ = self.tslits_dict['rcen'][:,slit].copy()
        inmask = (self.bpm == False)
        # Center of the slit
        slit_cen = (slit_left + slit_righ)/2.0

        slitmask = self.spectrograph.slitmask(self.tslits_dict)
        thismask = slitmask == slit

        # Tilt specific Optional parameters
        tracethresh = 10.0 # significance threshold for an arc line to be traced
        sigdetect = 5.0 # This is the significance threshold for finding neighboring lines. The number of these neighboring lines
        #  determines how many nsig > tracethresh lines that may be removed because they are too close.
        nfwhm_neigh = 3.0
        only_these_lines = None # These are lines from the wavelength solution that are known to be good. If set code only uses these
        # identified lines to avoid problems with ghosts (is this still necessary with new tracing code?)
        maxdev_tracefit = 1.0 # maximum deviation for iterative trace fitting (flux weighted, then gaussian weighted)
        sigrej_trace = 3.0 # From each line we compute a median absolute deviation of the trace from the polynomial fit. We then
        # analyze the distribution of mads for all the lines, and reject sigrej_trace outliers from that distribution.
        npca = 1 # Order of PCA for iteration to improve tracing
        coeff_npoly_pca = 1 # Polynomial order for fit to PCA coefficients to improve tracing
        sigrej_pca = 2.0 # Outlier rejection significance for PCA fit to arc line coefficients
        ncoeff = 5 # Order of legendre polynomial fits to the tilts

        max_badpix_frac = 0.20 # Maximum fraction of total pixels masked by the trace_gweight algorithm (because the residuals are too large)
        # Trace Crude parameters
        tcrude_maxerr = 1.0 #
        tcrude_maxshift = 3.0
        tcrude_maxshift0 = 3.0
        tcrude_nave = 5
        show_tracefits = False # show the trace fits

        # These are parameters for 2D fitting
        spec_order = 4
        maxdev_2d = 1.0 # Maximum absolute deviation used in rejection for 2d legendre polynomial fit to ensemble of all spectral line tilts
        sigrej_2d = 3.0 # Outlier significance threshold for rejection for 2d legendre polynomial fit to ensemble of all spectral line tilts
        show_tilts = True
        debug = True

        trc_tilt_dict0 = tracewave.trace_tilts_work(
            arcimg, lines_spec, lines_spat, thismask, inmask=inmask,
            tilts_guess=None,fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
            sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
            tcrude_maxerr=tcrude_maxerr, tcrude_maxshift=tcrude_maxshift,
            tcrude_maxshift0=tcrude_maxshift0,
            tcrude_nave=tcrude_nave,show_fits=show_tracefits)
        # Now evaluate the model of the tilts for all of our lines
        # Testing!!!!
        # Now perform a fit to the tilts
        tilt_fit_dict0 = tracewave.fit_tilts(trc_tilt_dict0, spat_order=ncoeff, spec_order=spec_order,debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')

        # Do a PCA fit, which rejects some outliers
        iuse = trc_tilt_dict0['use_tilt']
        nuse =np.sum(iuse)
        msgs.info('PCA modeling {:d} good tilts'.format(nuse))
        #pca_explained_var = None
        # TODO Should we truncate this PCA by hand, or just let it explain variance
        # PCA fit good orders, and predict bad orders
        pca_fit, poly_fit_dict, pca_mean, pca_vectors = extract.pca_trace(
            trc_tilt_dict0['tilts_sub_fit'], predict=np.invert(iuse), npca = npca, coeff_npoly=coeff_npoly_pca,
            lower=sigrej_pca, upper=sigrej_pca, order_vec=lines_spec,xinit_mean=lines_spec,
            minv=0.0, maxv=float(trc_tilt_dict0['nsub'] - 1), debug=True)

        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict1 = tracewave.trace_tilts_work(
            arcimg, lines_spec, lines_spat, thismask, inmask=inmask,
            tilts_guess=pca_fit, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
            sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
            show_fits=show_tracefits)


        # Now perform a fit to the tilts
        tilt_fit_dict1 = tracewave.fit_tilts(trc_tilt_dict1, spat_order=ncoeff, spec_order=spec_order, debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg1 = tracewave.fit2piximg(tilt_fit_dict1)

        # Now trace again with the piximg model as the starting crutches


        """
        # Since the y independent variable is the tilts in the way we do the 2d fit, and soem tilts are spurios, it does
        # no good to evaluate the global fit at these spurious tilts to get the new tracing crutches. The following is a
        # hack to determine this from the piximg. Perhaps there is a better way.
        spec_img = np.outer(np.arange(nspec), np.ones(nspat))
        spec_vec = np.arange(nspec)
        nlines=len(lines_spec)
        interp_flag = np.ones(nlines,dtype=bool)
        tilts_crutch = np.zeros((nspat, nlines))
        spat_min = trc_tilt_dict1['spat_min']
        spat_max = trc_tilt_dict1['spat_max']
        piximg1[np.invert(thismask)] = 1e10
        # Is there a faster more clever way to do this?
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            min_spec = int(np.fmax(np.round(lines_spec[iline]) - np.round(0.01*nspec),0))
            max_spec = int(np.fmin(np.round(lines_spec[iline]) + np.round(0.01*nspec),nspec-1))
            piximg_sub = piximg1[min_spec:max_spec,:]
            for ispat in range(min_spat,max_spat):
                if np.any(np.diff(piximg_sub[:,ispat] < 0)):
                    # If we ever encounter an unsorted piximg_sub the logic below makes no sense so just use the
                    # previous polynomial fit as the crutch and exit this loop
                    tilts_crutch[ispat,iline] = trc_tilt_dict1['tilts_fit'][:,iline]
                    msgs.warn('piximg is not monotonic. Your tilts are probably bad')
                    break
                else:
                    tilts_crutch[ispat,iline] = np.interp(lines_spec[iline],piximg_sub[:,ispat],spec_vec[min_spec:max_spec])


        trc_tilt_dict1['tilts_crutch'] = tilts_crutch
        # Now trace again with the PCA predictions as the starting crutches
        trc_tilt_dict2 = tracewave.trace_tilts(arcimg, lines_spec, lines_spat, slit_width, thismask, inmask=inmask,
                                               tilts_guess=tilts_crutch, fwhm=fwhm,ncoeff=ncoeff, maxdev_tracefit=maxdev_tracefit,
                                               sigrej_trace=sigrej_trace, max_badpix_frac=max_badpix_frac,
                                               show_fits=show_tracefits)

        # Now perform a second fit to the tilts
        tilt_fit_dict2 = tracewave.fit_tilts(trc_tilt_dict2, spat_order=ncoeff, spec_order=spec_order, debug=True,
                                             maxdev=maxdev_2d, sigrej=sigrej_2d,
                                             doqa=True,setup='test',slit=0, show_QA=False, out_dir='./')
        # Now evaluate the model of the tilts for all of our lines
        piximg2 = tracewave.fit2piximg(tilt_fit_dict2)
        """

#        from matplotlib import pyplot as plt
#        tilt_mask = trc_tilt_dict1['tilts_mask']
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], tilts_crutch[tilt_mask], 'ko', markersize=2.0)
#        plt.plot(trc_tilt_dict1['tilts_spat'][tilt_mask], trc_tilt_dict1['tilts_fit'][tilt_mask], 'ro', mfc='none', markersize=2.0)


        # Now do a fit


        """
        for iline in range(nlines):
            min_spat = np.fmax(spat_min[iline], 0)
            max_spat = np.fmin(spat_max[iline], nspat - 1)
            sub_img = (piximg*thismask)[:, min_spat:max_spat]
            spec_img_sub = spec_img[  :, min_spat:max_spat]
            ispec_min = np.argmin(np.abs(sub_img - lines_spec[iline]), axis=0)
            tilts_crutch[min_spat:max_spat,iline] = spec_img_sub[ispec_min,np.arange(sub_img.shape[1])]
        """



'''