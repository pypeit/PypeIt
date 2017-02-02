from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
import os
import copy
from pypit import arqa
from pypit import ararc
from pypit import armsgs
from pypit import arutils
from pypit import arpca
from pypit import arplot
from pypit import arparse as settings
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.ndimage as ndimage
from collections import Counter

# Logging
msgs = armsgs.get_logger()

from pypit import ardebug as debugger


def assign_slits(binarr, edgearr, ednum=100000, lor=-1):
    """
    This routine will traces the locations of the slit edges

    Parameters
    ----------
    binarr : numpy ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
    edgearr : numpy ndarray
      An array of negative/positive numbers (left/right edges respectively) and zeros (no edge)
    ednum : int
      A dummy number given to define slit edges
    lor : int (-1 or +1)
      A flag that indicates if the left edge (-1) or right edge (+1) should be assigned

    Returns
    -------
    edgearr : numpy ndarray
      An array of negative/positive numbers (left/right edges respectively) and zeros (no edge)
    """
    from pypit import arcytrace

    if lor == -1:
        lortxt = "left"
    else:
        lortxt = "right"
    outitnm = 1
    oldedgearr = edgearr.copy()
    prevedgearr = edgearr.copy()
    while True:
        msgs.prindent("Outer {0:s} edge loop, Iteration {1:d}".format(lortxt, outitnm))
        labnum = lor*ednum
        itnm = 0
        nslit = 0
        cmnold = None
        firstpass = True
        while True:
            edgehist = np.zeros(binarr.shape[1]*2, dtype=np.int)
            itnm += 1
            # Locate edges relative to the most common edge
            if lor == -1:
                wl = np.where(edgearr <= -2*ednum)
            else:
                wl = np.where(edgearr >= 2*ednum)
            if wl[0].size == 0:
                break
            cl = Counter(edg for edg in edgearr[wl])
            comml = cl.most_common(1)
            ww = np.where(edgearr == comml[0][0])
            if not firstpass:
                if (cmnold[0] == comml[0][0]) and (cmnold[1] == comml[0][1]):
                    # Nothing has changed since the previous iteration, so end the loop
                    break
                if comml[0][1] < binarr.shape[0]/100.0:
                    # Now considering an edge that spans less than 1 per cent of the detector ---> insignificant
                    break
            cmnold = comml[0]
            # Extract just these elements
            tedgearr = edgearr[ww[0], :]
            # Calculate the offset
            offs = binarr.shape[1]
            # Add these into edgehist
            edgehist[offs] = ww[0].size
            # And a fudge to give this edge detection some width (for peak finding, below)
            edgehist[offs-1] = 1 + ww[0].size/2
            edgehist[offs+1] = 1 + ww[0].size/2
            # Find the difference between unknown edges
            if lor == -1:
                www = np.where(tedgearr <= -2*ednum)
            else:
                www = np.where(tedgearr >= 2*ednum)
            if www[0].size == 0:
                break
            shft = www[1] - ww[1][www[0]]  # Calculate the shift between right edges
            shft += offs  # Apply the offset to the edgehist arr
            arcytrace.edge_sum(edgehist, shft)
            # Smooth the histogram with a Gaussian of standard deviation 1 pixel to reduce noise
            smedgehist = ndimage.uniform_filter1d(edgehist, 3)
            # Identify peaks (which indicate the locations of the right slit edges)
            arrlfr = smedgehist[0:-4]
            arrlft = smedgehist[1:-3]
            arrcen = smedgehist[2:-2]
            arrrgt = smedgehist[3:-1]
            arrrfr = smedgehist[4:]
            wpk = np.where((arrcen >= arrlft) & (arrcen > arrrgt) &  # Exactly one of these should be >=
                           ((arrlft > arrlfr) | (arrrgt > arrrfr)))[0]
            if wpk.size == 0:
                # No more peaks
                break
            if wpk.size != 1:
                wpkmsk = arcytrace.prune_peaks(smedgehist, wpk, np.where(wpk+2 == offs)[0][0])
                wpk = wpk[np.where(wpkmsk == 1)]
            if wpk.size == 0:
                # After pruning, there are no more peaks
                break
            pks = wpk+2  # Shifted by 2 because of the peak finding algorithm above
            pedges = arcytrace.find_peak_limits(smedgehist, pks)
            if np.all(pedges[:, 1]-pedges[:, 0] == 0):
                # Remaining peaks have no width
                break
            if msgs._debug['trace'] and False:
                plt.clf()
                plt.plot(arrcen, 'k-', drawstyle='steps')
                plt.plot(wpk, np.zeros(wpk.size), 'ro')
                plt.show()
            # Label all edge ids (in the original edgearr) that are located in each peak with the same number
            for ii in range(pks.size):
                shbad = np.zeros(edgearr.shape)
                wp = np.where((shft >= pedges[ii, 0]) & (shft <= pedges[ii, 1]))
                vals = np.unique(tedgearr[(www[0][wp], www[1][wp])])
                # Fit the edge detections in this edge and calculate the offsets
                strev = "np.where("
                for vv in vals:
                    strev += "(edgearr=={0:d})|".format(vv)
                strev = strev[:-1] + ")"
                widx = eval(strev)
                if widx[0].size < 2*settings.argflag['trace']['slits']['polyorder']:
                    continue
                badmsk, fitcof = arutils.robust_polyfit(widx[0], widx[1],
                                                        settings.argflag['trace']['slits']['polyorder'],
                                                        function=settings.argflag['trace']['slits']['function'],
                                                        minv=0, maxv=binarr.shape[0]-1)
                shbad[widx] = badmsk
                smallhist = np.zeros(101, dtype=np.int)
                meddiff = np.zeros(vals.size)
                for vv in range(vals.size):
                    widx = np.where((edgearr == vals[vv]) & (shbad == 0))
                    if widx[0].size == 0:
                        # These pixels were deemed to be bad
                        continue
                    diff = widx[1] - arutils.func_val(fitcof, widx[0],
                                                      settings.argflag['trace']['slits']['function'],
                                                      minv=0, maxv=binarr.shape[0]-1)
                    diff = 50 + np.round(diff).astype(np.int)
                    arcytrace.edge_sum(smallhist, diff)
                    meddiff[vv] = np.median(diff)
                # Find the peaks of this distribution
                wspk = np.where((smallhist[1:-1] >= smallhist[2:]) & (smallhist[1:-1] > smallhist[:-2]))[0]
                wspk += 1  # Add one here to account for peak finding
                if msgs._debug['trace'] and False:
                    plt.clf()
                    plt.plot(smallhist, 'k-', drawstyle='steps')
                    plt.show()

                for pp in range(wspk.size):  # For all small peaks identified
                    for vv in range(vals.size):
                        if lor == -1 and vals[vv] > -2*ednum:
                            continue
                        elif lor == 1 and vals[vv] < 2*ednum:
                            continue
                        # Make sure this value is within 1 pixel of the peak
                        if meddiff[vv] < wspk[pp]-1:
                            continue
                        if meddiff[vv] > wspk[pp]+1:
                            continue
                        edgearr[np.where(edgearr == vals[vv])] = labnum
                        meddiff[vv] = -1  # Flag this val as done
                    labnum += lor*1
                # Find any vals that weren't applied
                for vv in range(vals.size):
                    if meddiff[vv] == -1:
                        continue
                    edgearr[np.where(edgearr == vals[vv])] = 0
            nslit += pks.size
            msgs.prindent("  Inner loop, Iteration {0:d}, {1:d} {2:s} edges assigned ({3:d} total)".format(itnm, pks.size, lortxt, nslit))
            firstpass = False
        outitnm += 1
        if lor == -1:
            edgearr[np.where(edgearr <= -2*ednum)] = 0
        else:
            edgearr[np.where(edgearr >= 2*ednum)] = 0
        if np.array_equal(edgearr, oldedgearr) or np.array_equal(edgearr, prevedgearr):
            break
        elif outitnm > 10:
            msgs.warn("Edge assignment may not have converged")
            msgs.info("Please check the slit edge traces")
            #debugger.set_trace()
            break
        else:
            oldedgearr = prevedgearr.copy()
            prevedgearr = edgearr.copy()
            if lor == -1:
                edgearr[np.where(edgearr <= -ednum)] -= ednum
            else:
                edgearr[np.where(edgearr >= ednum)] += ednum
    # Ignore any order detections that weren't identified in the loop
    if lor == -1:
        edgearr[np.where(edgearr <= -2*ednum)] = 0
    else:
        edgearr[np.where(edgearr >= 2*ednum)] = 0
    # Sort vals by increasing spatial position on the detector
    # First, determine the model for the most common slit edge
    if lor == -1:
        wcm = np.where(edgearr <= -ednum)
    else:
        wcm = np.where(edgearr >= ednum)
    if wcm[0].size != 0:
        cntr = Counter(edg for edg in edgearr[wcm])
        commn = cntr.most_common(1)
        wedx, wedy = np.where(edgearr == commn[0][0])
        msk, cf = arutils.robust_polyfit(wedx, wedy,
                                         settings.argflag['trace']['slits']['polyorder'],
                                         function=settings.argflag['trace']['slits']['function'],
                                         minv=0, maxv=binarr.shape[0]-1)
        cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                   settings.argflag['trace']['slits']['function'],
                                   minv=0, maxv=binarr.shape[0]-1)
        if lor == -1:
            vals = np.unique(edgearr[np.where(edgearr < 0)])
        else:
            vals = np.unique(edgearr[np.where(edgearr > 0)])
        diffarr = np.zeros(vals.size)
        diffstd = 0.0
        for jj in range(vals.size):
            wedx, wedy = np.where(edgearr == vals[jj])
            diffarr[jj] = np.mean(wedy-cenmodl[wedx])
            diffstd += np.std(wedy-cenmodl[wedx])
        diffstd /= vals.size
        dasrt = np.argsort(diffarr)
        # Relabel the edges from left to right
        edgearr[wcm] += lor*ednum
        labnum = lor*ednum
        diffarrsrt = diffarr[dasrt]
        diffs = diffarrsrt[1:] - diffarrsrt[:-1]
        for jj in range(vals.size):
            wrplc = np.where(edgearr == lor*ednum + vals[dasrt[jj]])
            edgearr[wrplc] = labnum
            if jj != vals.size-1:
                if diffs[jj] > 3.0*diffstd:
                    # The next edge must be a different edge
                    labnum += lor*1
    return


def expand_slits(slf, mstrace, det, ordcen, extord):
    """
    This routine will traces the locations of the slit edges

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    mstrace: numpy ndarray
      Calibration frame that will be used to identify slit edges
    det : int
      Index of the detector
    ordcen : ndarray
      An array providing the physical pixel locations corresponding to the slit centres
    extord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)

    Returns
    -------
    lordloc : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    rordloc : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    """
    from pypit import arcytrace

    # Calculate the pixel locations of th eorder edges
    pixcen = phys_to_pix(ordcen, slf._pixlocn[det - 1], 1)
    msgs.info("Expanding slit traces to slit edges")
    mordwid, pordwid = arcytrace.expand_slits(mstrace, pixcen, extord.astype(np.int))
    # Fit a function for the difference between left edge and the centre trace
    ldiff_coeff, ldiff_fit = arutils.polyfitter2d(mordwid, mask=-1,
                                                  order=settings.argflag['trace']['slits']['diffpolyorder'])
    # Fit a function for the difference between left edge and the centre trace
    rdiff_coeff, rdiff_fit = arutils.polyfitter2d(pordwid, mask=-1,
                                                  order=settings.argflag['trace']['slits']['diffpolyorder'])
    lordloc = ordcen - ldiff_fit.T
    rordloc = ordcen + rdiff_fit.T
    return lordloc, rordloc


def trace_object(slf, det, sciframe, varframe, crmask, trim=2.0,
                 triml=None, trimr=None, sigmin=2.0, bgreg=None,
                 maskval=-999999.9, order=0, doqa=True):
    """ Finds objects, and traces their location on the detector
    Parameters
    ----------
    slf
    sciframe
    varframe
    crmask
    trim
    triml
    trimr
    sigmin
    bgreg
    maskval
    order
    doqa

    Returns
    -------

    """
    from pypit import arcytrace
    from pypit import arcyutils
    smthby = 7
    rejhilo = 1
    bgreg = 20
    traceorder = 2   # Order of polynomial used to trace the objects
    if triml is None: triml = trim
    if trimr is None: trimr = trim
    npix = int(slf._pixwid[det-1][order] - triml - trimr)
    if bgreg is None: bgreg = npix
    # Interpolate the science array onto a new grid (with constant spatial slit length)
    msgs.info("Rectifying science frame")
    xint = np.linspace(0.0, 1.0, sciframe.shape[0])
    yint = np.linspace(0.0, 1.0, sciframe.shape[1])
    scispl = interp.RectBivariateSpline(xint, yint, sciframe, bbox=[0.0, 1.0, 0.0, 1.0], kx=1, ky=1, s=0)
    varspl = interp.RectBivariateSpline(xint, yint, varframe, bbox=[0.0, 1.0, 0.0, 1.0], kx=1, ky=1, s=0)
    crmspl = interp.RectBivariateSpline(xint, yint, crmask, bbox=[0.0, 1.0, 0.0, 1.0], kx=1, ky=1, s=0)
    xx, yy = np.meshgrid(np.linspace(0.0,1.0,sciframe.shape[0]),np.linspace(0.0,1.0,npix), indexing='ij')
    ro = (slf._rordloc[det-1][:,order]-trimr).reshape((-1,1))/(sciframe.shape[1]-1.0)
    lo = (slf._lordloc[det-1][:,order]+triml).reshape((-1,1))/(sciframe.shape[1]-1.0)
    vv = (lo+(ro-lo)*yy).flatten()
    xx = xx.flatten()
    recsh = (sciframe.shape[0],npix)
    rec_sciframe = scispl.ev(xx, vv).reshape(recsh)
    rec_varframe = varspl.ev(xx, vv).reshape(recsh)
    rec_crmask   = crmspl.ev(xx, vv).reshape(recsh)
    # Update the CR mask to ensure it only contains 1's and 0's
    rec_crmask[np.where(rec_crmask>0.2)] = 1.0
    rec_crmask[np.where(rec_crmask<=0.2)] = 0.0
    msgs.info("Estimating object profiles")
    # Smooth the S/N frame
    rec_sigframe_bin = arcyutils.smooth_x(rec_sciframe/np.sqrt(rec_varframe), 1.0-rec_crmask, smthby, rejhilo, maskval)
    #rec_varframe_bin = arcyutils.smooth_x(rec_varframe, 1.0-rec_crmask, smthby, rejhilo, maskval)
    #rec_sigframe_bin = np.sqrt(rec_varframe_bin/(smthby-2.0*rejhilo))
    #sigframe = rec_sciframe_bin*(1.0-rec_crmask)/rec_sigframe_bin
    ww = np.where(rec_crmask==0.0)
    med, mad = arutils.robust_meanstd(rec_sigframe_bin[ww])
    ww = np.where(rec_crmask==1.0)
    rec_sigframe_bin[ww] = maskval
    srtsig = np.sort(rec_sigframe_bin,axis=1)
    ww = np.where(srtsig[:,-2] > med + sigmin*mad)
    mask_sigframe = np.ma.array(rec_sigframe_bin, mask=rec_crmask, fill_value=maskval)
    # Collapse along the spectral direction to get object profile
    trcprof = np.ma.mean(mask_sigframe[ww[0],:], axis=0).filled(0.0)
    trcxrng = np.arange(npix)/(npix-1.0)
    msgs.info("Identifying objects that are significantly detected")
    # Find significantly detected objects
    mskpix, coeff = arutils.robust_polyfit(trcxrng, trcprof, 1+npix//40, function='legendre', sigma=2.0, minv=0.0, maxv=1.0)
    backg = arutils.func_val(coeff, trcxrng, 'legendre', minv=0.0, maxv=1.0)
    trcprof -= backg
    wm = np.where(mskpix==0)
    if wm[0].size==0:
        msgs.warn("No objects found")
        return None
    med, mad = arutils.robust_meanstd(trcprof[wm])
    trcprof -= med
    objl, objr, bckl, bckr = arcytrace.find_objects(trcprof, bgreg, mad)
    #plt.clf()
    #plt.plot(trcxrng,trcprof,'k-')
    #wl = np.where(bckl[:,0]==1)
    #plt.plot(trcxrng[wl],trcprof[wl],'ro')
    #wr = np.where(bckr[:,0]==1)
    #plt.plot(trcxrng[wr],trcprof[wr],'go')
    #plt.plot(trcxrng[objl[0]:objr[0]+1],trcprof[objl[0]:objr[0]+1],'bx')
    #plt.show()
    nobj = objl.size
    if msgs._debug['trace_obj']:
        debugger.set_trace()
        #nobj = 1
    if nobj==1:
        msgs.info("Found {0:d} object".format(objl.size))
        msgs.info("Tracing {0:d} object".format(objl.size))
    else:
        msgs.info("Found {0:d} objects".format(objl.size))
        msgs.info("Tracing {0:d} objects".format(objl.size))
    # Max obj
    if nobj > settings.argflag['science']['extraction']['maxnumber']:
        nobj = settings.argflag['science']['extraction']['maxnumber']
        msgs.warn("Restricting to the brightest {:d} objects found".format(nobj))
    # Trace objects
    cval = np.zeros(nobj)
    allsfit = np.array([])
    allxfit = np.array([])
    for o in range(nobj):
        xfit = np.arange(objl[o],objr[o]).reshape((1,-1))/(npix-1.0)
        cent = np.ma.sum(mask_sigframe[:,objl[o]:objr[o]]*xfit, axis=1)
        wght = np.ma.sum(mask_sigframe[:,objl[o]:objr[o]], axis=1)
        cent /= wght
        centfit = cent.filled(maskval)
        specfit = np.linspace(-1.0, 1.0, sciframe.shape[0])
        w = np.where(centfit != maskval)
        specfit = specfit[w]
        centfit = centfit[w]
        mskbad, coeffs = arutils.robust_polyfit(specfit,centfit,traceorder,function="legendre", minv=-1.0, maxv=1.0)
        cval[o] = arutils.func_val(coeffs, np.array([0.0]), "legendre", minv=-1.0, maxv=1.0)[0]
        w = np.where(mskbad==0.0)
        if w[0].size!=0:
            allxfit = np.append(allxfit, specfit[w])
            allsfit = np.append(allsfit, centfit[w]-cval[o])
    if nobj == 0:
        msgs.warn("No objects detected in slit")
        return dict(nobj=0, traces=None, object=None, background=None)
    # Tracing
    msgs.info("Performing global trace to all objects")
    mskbad, coeffs = arutils.robust_polyfit(allxfit,allsfit,traceorder,function="legendre", minv=-1.0, maxv=1.0)
    trcfunc = arutils.func_val(coeffs, np.linspace(-1.0, 1.0, sciframe.shape[0]), "legendre", minv=-1.0, maxv=1.0)
    msgs.info("Constructing a trace for all objects")
    trcfunc = trcfunc.reshape((-1,1)).repeat(nobj, axis=1)
    trccopy = trcfunc.copy()
    for o in range(nobj): trcfunc[:,o] += cval[o]
    if nobj==1: msgs.info("Converting object trace to detector pixels")
    else: msgs.info("Converting object traces to detector pixels")
    ofst = slf._lordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1) + triml
    diff = (slf._rordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1)
            - slf._lordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1))
    # Convert central trace
    traces = ofst + (diff-triml-trimr)*trcfunc
    # Convert left object trace
    for o in range(nobj): trccopy[:,o] = trcfunc[:,o] - cval[o] + objl[o]/(npix-1.0)
    trobjl = ofst + (diff-triml-trimr)*trccopy
    # Convert right object trace
    for o in range(nobj): trccopy[:,o] = trcfunc[:,o] - cval[o] + objr[o]/(npix-1.0)
    trobjr = ofst + (diff-triml-trimr)*trccopy
    # Make an image of pixel weights for each object
    xint = np.linspace(0.0, 1.0, sciframe.shape[0])
    yint = np.linspace(0.0, 1.0, npix)
    yint = np.append(-yint[1],np.append(yint,2.0-yint[-2]))
    msgs.info("Creating an image weighted by object pixels")
    #plt.plot(trcxrng[objl[0]:objr[0]+1],trcprof[objl[0]:objr[0]+1],'bx')
    rec_obj_img = np.zeros((sciframe.shape[0],sciframe.shape[1],nobj))
    for o in range(nobj):
        obj = np.zeros(npix)
        obj[objl[o]:objr[o]+1]=1
        scitmp = np.append(0.0,np.append(obj,0.0))
        objframe = scitmp.reshape(1,-1).repeat(sciframe.shape[0],axis=0)
        objspl = interp.RectBivariateSpline(xint, yint, objframe, bbox=[0.0, 1.0, yint.min(), yint.max()], kx=1, ky=1, s=0)
        xx, yy = np.meshgrid(np.linspace(0,1.0,sciframe.shape[0]), np.arange(0,sciframe.shape[1]), indexing='ij')
        lo = (slf._lordloc[det-1][:,order]+triml).reshape((-1,1))
        vv = ((yy-lo)/(npix-1.0)).flatten()
        xx = xx.flatten()
        wf = np.where((vv>=yint[0])&(vv<=yint[-1]))
        rec_obj_arr = objspl.ev(xx[wf], vv[wf])
        idxarr = np.zeros(sciframe.shape).flatten()
        idxarr[wf] = 1
        idxarr = idxarr.reshape(sciframe.shape)
        rec_img = np.zeros_like(sciframe)
        rec_img[np.where(idxarr==1)] = rec_obj_arr
        rec_obj_img[:,:,o] = rec_img.copy()
    # Make an image of pixel weights for the background region of each object
    msgs.info("Creating an image weighted by background pixels")
    rec_bg_img = np.zeros((sciframe.shape[0],sciframe.shape[1],nobj))
    for o in range(nobj):
        backtmp = np.append(0.0,np.append(bckl[:,o]+bckr[:,o],0.0))
        bckframe = backtmp.reshape(1,-1).repeat(sciframe.shape[0],axis=0)
        bckspl = interp.RectBivariateSpline(xint, yint, bckframe, bbox=[0.0, 1.0, yint.min(), yint.max()], kx=1, ky=1, s=0)
        xx, yy = np.meshgrid(np.linspace(0,1.0,sciframe.shape[0]), np.arange(0,sciframe.shape[1]), indexing='ij')
        lo = (slf._lordloc[det-1][:,order]+triml).reshape((-1,1))
        vv = ((yy-lo)/(npix-1.0)).flatten()
        xx = xx.flatten()
        wf = np.where((vv>=yint[0])&(vv<=yint[-1]))
        rec_bg_arr = bckspl.ev(xx[wf], vv[wf])
        idxarr = np.zeros(sciframe.shape).flatten()
        idxarr[wf] = 1
        idxarr = idxarr.reshape(sciframe.shape)
        rec_img = np.zeros_like(sciframe)
        rec_img[np.where(idxarr==1)] = rec_bg_arr
        rec_bg_img[:,:,o] = rec_img.copy()
        #arutils.ds9plot(rec_img)
    # Save the quality control
    if doqa and (not msgs._debug['no_qa']):
        arqa.obj_trace_qa(slf, sciframe, trobjl, trobjr, root="object_trace", normalize=False)
    # Trace dict
    tracedict = dict({})
    tracedict['nobj'] = nobj
    tracedict['traces'] = traces
    tracedict['object'] = rec_obj_img
    tracedict['background'] = rec_bg_img
    return tracedict


def trace_slits(slf, mstrace, det, pcadesc="", maskBadRows=False):
    """
    This routine traces the locations of the slit edges

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    mstrace : numpy ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
    det : int
      Index of the detector
    pcadesc : str, optional
      A descriptive string of text to be annotated as a title on the QA PCA plots
    maskBadRows : bool, optional
      Mostly useful for echelle data where the slit edges are bent relative to
      the pixel columns. Do not set this keyword to True if slit edges are
      almost aligned with the pixel columns.

    Returns
    -------
    lcenint : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    rcenint : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    dnum = settings.get_dnum(det)
    ednum = 100000  # A large dummy number used for slit edge assignment. ednum should be larger than the number of edges detected
    from pypit import arcytrace

    msgs.info("Preparing trace frame for slit detection")
    # Generate a binned (or smoothed) version of the trace frame
    binarr = ndimage.uniform_filter(mstrace, size=(3, 1))
    binbpx = slf._bpix[det-1].copy()
    plxbin = slf._pixlocn[det-1][:, :, 0].copy()
    plybin = slf._pixlocn[det-1][:, :, 1].copy()
    if False:
        # Use this for debugging
        binbpx = np.zeros(mstrace.shape, dtype=np.int)
        xs = np.arange(mstrace.shape[0] * 1.0) * settings.spect[dnum]['xgap']
        xt = 0.5 + np.arange(mstrace.shape[0] * 1.0) + xs
        ys = np.arange(mstrace.shape[1]) * settings.spect[dnum]['ygap'] * settings.spect[dnum]['ysize']
        yt = settings.spect[dnum]['ysize'] * (0.5 + np.arange(mstrace.shape[1] * 1.0)) + ys
        xloc, yloc = np.meshgrid(xt, yt)
        plxbin, plybin = xloc.T, yloc.T
    #    binby = 5
#    binarr = arcyutils.bin_x(mstrace, binby, 0)
#    binbpx = arcyutils.bin_x(slf._bpix[det-1], binby, 0)
#    plxbin = arcyutils.bin_x(slf._pixlocn[det-1][:,:,0], binby, 1)
#    plybin = arcyutils.bin_x(slf._pixlocn[det-1][:,:,1], binby, 1)

    # Specify how many times to repeat the median filter
    medrep = 3
    if len(settings.argflag['trace']['slits']['single']) > 0:
        edgearr = np.zeros(binarr.shape, dtype=np.int)
        # Add a user-defined slit?
        # Syntax is a list of values, 2 per detector that define the slit
        # according to column values.  The 2nd value (for the right edge)
        # must be >0 to be applied.  Example for LRISr [-1, -1, 7, 295]
        # which means the code skips user-definition for the first detector
        # but adds one for the 2nd.
        ledge, redge = (det-1)*2, (det-1)*2+1
        if settings.argflag['trace']['slits']['single'][redge] > 0:
            msgs.warn("Using input slit edges on detector {:d}: [{:g},{:g}]".format(
                    det,
                    settings.argflag['trace']['slits']['single'][ledge],
                    settings.argflag['trace']['slits']['single'][redge]))
            msgs.warn("Better know what you are doing!")
            edgearr[:, settings.argflag['trace']['slits']['single'][ledge]] = -1
            edgearr[:, settings.argflag['trace']['slits']['single'][redge]] = +1
    else:
        # Even better would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian near the maximum in each column
        msgs.info("Detecting slit edges")
        sqmstrace = np.sqrt(np.abs(binarr))
        for ii in range(medrep):
            sqmstrace = ndimage.median_filter(sqmstrace, size=(3, 7))
        # Make sure there are no spuriously low pixels
        sqmstrace[(sqmstrace < 1.0) & (sqmstrace >= 0.0)] = 1.0
        sqmstrace[(sqmstrace > -1.0) & (sqmstrace <= 0.0)] = -1.0
        # Apply a Sobel filter
        filt = ndimage.sobel(sqmstrace, axis=1, mode='nearest')
        msgs.info("Applying bad pixel mask")
        filt *= (1.0 - binbpx)  # Apply to the bad pixel mask
        siglev = np.sign(filt)*(filt**2)/sqmstrace
        tedges = np.zeros(binarr.shape, dtype=np.float)
        wl = np.where(siglev > +settings.argflag['trace']['slits']['sigdetect'])  # A positive gradient is a left edge
        wr = np.where(siglev < -settings.argflag['trace']['slits']['sigdetect'])  # A negative gradient is a right edge
        tedges[wl] = -1.0
        tedges[wr] = +1.0
        # import astropy.io.fits as pyfits
        # hdu = pyfits.PrimaryHDU(filt)
        # hdu.writeto("filt_{0:02d}.fits".format(det))
        # hdu = pyfits.PrimaryHDU(sqmstrace)
        # hdu.writeto("sqmstrace_{0:02d}.fits".format(det))
        # hdu = pyfits.PrimaryHDU(binarr)
        # hdu.writeto("binarr_{0:02d}.fits".format(det))
        # hdu = pyfits.PrimaryHDU(siglev)
        # hdu.writeto("siglev_{0:02d}.fits".format(det))
        nedgear = arcytrace.clean_edges(siglev, tedges)
        if maskBadRows:
            msgs.info("Searching for bad pixel rows")
            edgsum = np.sum(nedgear, axis=0)
            sigma = 1.4826*np.median(np.abs(edgsum-np.median(edgsum)))
            w = np.where(np.abs(edgsum) >= 1.5*sigma)[0]
            maskcols = np.unique(np.append(w, np.append(w+1, w-1)))
            msgs.info("Masking {0:d} bad pixel rows".format(maskcols.size))
            for i in range(maskcols.size):
                if maskcols[i] < 0 or maskcols[i] >= nedgear.shape[1]:
                    continue
                nedgear[:, maskcols[i]] = 0
        ######
        msgs.info("Applying bad pixel mask")
        nedgear *= (1-binbpx.astype(np.int))  # Apply to the bad pixel mask
        # eroll = np.roll(binbpx, 1, axis=1)
        # eroll[:,0] = eroll[:,1]
        # nedgear *= (1.0-eroll)  # Apply to the new detection algorithm (with shift)
        # Now roll back
        # nedgear = np.roll(nedgear, -1, axis=1)
        # edgearr = np.zeros_like(nedgear)
        # edgearr[np.where((nedgear == +1) | (tedgear == +1))] = +1
        # edgearr[np.where((nedgear == -1) | (tedgear == -1))] = -1
        sigedg = np.copy(siglev)
        sigedg[np.where(nedgear == 0)] = 0
        sigedg[np.where(np.isinf(sigedg) | np.isnan(sigedg))] = 0
        if settings.argflag['trace']['slits']['number'] > 0:
            # Now identify the number of most significantly detected peaks (specified by the user)
            amnmx = np.argsort(sigedg, axis=1)
            edgearr = np.zeros(binarr.shape, dtype=np.int)
            xsm = np.arange(binarr.shape[0])
            for ii in range(0, settings.argflag['trace']['slits']['number']):
                wset = np.where(sigedg[(xsm, amnmx[:, ii])] != 0)
                edgearr[(wset[0], amnmx[wset[0], ii])] = 1
                wset = np.where(sigedg[(xsm, amnmx[:, amnmx.shape[1] - 1 - ii])] != 0)
                edgearr[(wset[0], amnmx[wset[0], amnmx.shape[1] - 1 - ii])] = -1
        else:
            edgearr = np.copy(nedgear)

    # Assign a number to each of the edges
    msgs.info("Matching slit edges")
    lcnt, rcnt = arcytrace.match_edges(edgearr, ednum)
    if lcnt >= ednum or rcnt >= ednum:
        msgs.error("Found more edges than allowed by ednum. Set ednum to a larger number.")
    if lcnt == 1:
        letxt = "edge"
    else:
        letxt = "edges"
    if rcnt == 1:
        retxt = "edge"
    else:
        retxt = "edges"
    msgs.info("{0:d} left {1:s} and {2:d} right {3:s} were found in the trace".format(lcnt, letxt, rcnt, retxt))
    if (lcnt == 0) and (rcnt == 0):
        if np.median(binarr) > 500:
            msgs.warn("Found flux but no edges.  Assuming they go to the edge of the detector.")
            edgearr[:, -1] = 2*ednum
            rcnt = 1
            edgearr[:, 0] = -2*ednum
            lcnt = 1
        else:
            msgs.error("Unable to trace any edges"+msgs.newline()+"try a different method to trace the order edges")
    elif rcnt == 0:
        msgs.warn("Unable to find a right edge. Adding one in.")
        edgearr[:, -1] = 2*ednum
        rcnt = 1
    elif lcnt == 0:
        msgs.warn("Unable to find a left edge. Adding one in.")
        edgearr[:, 0] = -2*ednum
        lcnt = 1
    msgs.info("Assigning slit edge traces")
    # Find the most common set of edges
    minvf, maxvf = plxbin[0, 0], plxbin[-1, 0]
    edgearrcp = edgearr.copy()
    # If slits are set as "close" by the user, take the absolute value
    # of the detections and ignore the left/right edge detections
    if settings.argflag['trace']['slits']['maxgap'] is not None:
        edgearrcp[np.where(edgearrcp < 0)] += 1 + np.max(edgearrcp) - np.min(edgearrcp)
    # Assign left edges
    msgs.info("Assigning left slit edges")
    if lcnt == 1:
        edgearrcp[np.where(edgearrcp <= -2*ednum)] = -ednum
    else:
        assign_slits(binarr, edgearrcp, lor=-1)
    # Assign right edges
    msgs.info("Assigning right slit edges")
    if rcnt == 1:
        edgearrcp[np.where(edgearrcp >= 2*ednum)] = ednum
    else:
        assign_slits(binarr, edgearrcp, lor=+1)
    if settings.argflag['trace']['slits']['maxgap'] is not None:
        vals = np.sort(np.unique(edgearrcp[np.where(edgearrcp != 0)]))
        hasedge = arcytrace.close_edges(edgearrcp, vals, int(settings.argflag['trace']['slits']['maxgap']))
        # Find all duplicate edges
        edgedup = vals[np.where(hasedge == 1)]
        if edgedup.size > 0:
            for jj in range(edgedup.size):
                # Raise all remaining edges by one
                if jj != edgedup.size-1:
                    edgedup[jj+1:] += 1
                edgearrcp[np.where(edgearrcp > edgedup[jj])] += 1
                # Now investigate the duplicate
                wdup = np.where(edgearrcp == edgedup[jj])
                alldup = edgearr[wdup]
                alldupu = np.unique(alldup)
                cntr = Counter(edg for edg in alldup)
                commn = cntr.most_common(alldupu.size)
                shftsml = np.zeros(len(commn))
                shftarr = np.zeros(wdup[0].size, dtype=np.int)
                wghtarr = np.zeros(len(commn))
                duploc = [None for ii in range(len(commn))]
                for ii in range(len(commn)):
                    wghtarr[ii] = commn[ii][1]
                    duploc[ii] = np.where(edgearr[wdup] == commn[ii][0])
                changesmade = True
                while changesmade:
                    # Keep shifting pixels until the best match is found
                    changesmade = False
                    # First calculate the old model
                    cf = arutils.func_fit(wdup[0], wdup[1]+shftarr,
                                          settings.argflag['trace']['slits']['function'],
                                          settings.argflag['trace']['slits']['polyorder'],
                                          minv=0, maxv=binarr.shape[0]-1)
                    cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                               settings.argflag['trace']['slits']['function'],
                                               minv=0, maxv=binarr.shape[0]-1)
                    chisqold = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr).sum()
                    for ii in range(1, len(commn)):
                        # Shift by +1
                        adj = np.zeros(wdup[0].size)
                        adj[duploc[ii]] += 1
                        cf = arutils.func_fit(wdup[0], wdup[1]+shftarr+adj,
                                              settings.argflag['trace']['slits']['function'],
                                              settings.argflag['trace']['slits']['polyorder'],
                                              minv=0, maxv=binarr.shape[0]-1)
                        cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                                   settings.argflag['trace']['slits']['function'],
                                                   minv=0, maxv=binarr.shape[0]-1)
                        chisqp = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr-adj).sum()
                        # Shift by -1
                        adj = np.zeros(wdup[0].size)
                        adj[duploc[ii]] -= 1
                        cf = arutils.func_fit(wdup[0], wdup[1]+shftarr+adj,
                                              settings.argflag['trace']['slits']['function'],
                                              settings.argflag['trace']['slits']['polyorder'],
                                              minv=0, maxv=binarr.shape[0]-1)
                        cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                                                   settings.argflag['trace']['slits']['function'],
                                                   minv=0, maxv=binarr.shape[0]-1)
                        chisqm = np.abs(cenmodl[wdup[0]]-wdup[1]-shftarr-adj).sum()
                        # Test which solution is best:
                        if chisqold < chisqp and chisqold < chisqm:
                            # No changes are needed
                            continue
                        else:
                            changesmade = True
                            if chisqp < chisqm:
                                shftarr[duploc[ii]] += 1
                                shftsml[ii] += 1
                            else:
                                shftarr[duploc[ii]] -= 1
                                shftsml[ii] -= 1
                # Find the two most common edges
                cntr = Counter(sarr for sarr in shftarr)
                commn = cntr.most_common(2)
                if commn[0][0] > commn[1][0]:  # Make sure that suffix 'a' is assigned the leftmost edge
                    wdda = np.where(shftarr == commn[0][0])
                    wddb = np.where(shftarr == commn[1][0])
                    shadj = commn[0][0] - commn[1][0]
                else:
                    wdda = np.where(shftarr == commn[1][0])
                    wddb = np.where(shftarr == commn[0][0])
                    shadj = commn[1][0] - commn[0][0]
                wvla = np.unique(edgearr[wdup][wdda])
                wvlb = np.unique(edgearr[wdup][wddb])
                # Now generate the dual edge
                arcytrace.dual_edge(edgearr, edgearrcp, wdup[0], wdup[1], wvla, wvlb, shadj,
                                    int(settings.argflag['trace']['slits']['maxgap']), edgedup[jj])
        # Now introduce new edge locations
        vals = np.sort(np.unique(edgearrcp[np.where(edgearrcp != 0)]))
        edgearrcp = arcytrace.close_slits(binarr, edgearrcp, vals, int(settings.argflag['trace']['slits']['maxgap']))
    # Update edgearr
    edgearr = edgearrcp.copy()
    iterate = True
    # Get the final estimates of left and right edges
    eaunq = np.unique(edgearr)
    lcnt = np.where(eaunq < 0)[0].size
    rcnt = np.where(eaunq > 0)[0].size
    if lcnt == 0:
        msgs.warn("Unable to find a left edge. Adding one in.")
        edgearr[:, 0] = -2 * ednum
        lcnt = 1
    if rcnt == 0:
        msgs.warn("Unable to find a right edge. Adding one in.")
        edgearr[:, -1] = 2 * ednum
        rcnt = 1
    if (lcnt == 1) & (rcnt > 1):
        msgs.warn("Only one left edge, and multiple right edges.")
        msgs.info("Restricting right edge detection to the most significantly detected edge.")
        wtst = np.where(eaunq > 0)[0]
        bval, bidx = -np.median(siglev[np.where(edgearr == eaunq[wtst[0]])]), 0
        for r in range(1, rcnt):
            wed = np.where(edgearr == eaunq[wtst[r]])
            tstv = -np.median(siglev[wed])
            if tstv > bval:
                bval = tstv
                bidx = r
        edgearr[np.where((edgearr > 0) & (edgearr != eaunq[wtst[bidx]]))] = 0
        edgearr[np.where(edgearr > 0)] = +ednum  # Reset the edge value
        rcnt = 1
    if (lcnt > 1) & (rcnt == 1):
        msgs.warn("Only one right edge, and multiple left edges.")
        msgs.info("Restricting left edge detection to the most significantly detected edge.")
        wtst = np.where(eaunq < 0)[0]
        bval, bidx = np.median(siglev[np.where(edgearr == eaunq[wtst[0]])]), 0
        for r in range(1, lcnt):
            wed = np.where(edgearr == eaunq[wtst[r]])
            tstv = np.median(siglev[wed])
            if tstv > bval:
                bval = tstv
                bidx = r
        edgearr[np.where((edgearr < 0) & (edgearr != eaunq[wtst[bidx]]))] = 0
        edgearr[np.where(edgearr < 0)] = -ednum  # Reset the edge value
        lcnt = 1
    if lcnt == 1:
        letxt = "edge"
    else:
        letxt = "edges"
    if rcnt == 1:
        retxt = "edge"
    else:
        retxt = "edges"
    msgs.info("{0:d} left {1:s} and {2:d} right {3:s} were found in the trace".format(lcnt, letxt, rcnt, retxt))
    while iterate:
        # Calculate the minimum and maximum left/right edges
        ww = np.where(edgearr < 0)
        lmin, lmax = -np.max(edgearr[ww]), -np.min(edgearr[ww])  # min/max are switched because of the negative signs
        ww = np.where(edgearr > 0)
        rmin, rmax = np.min(edgearr[ww]), np.max(edgearr[ww])  # min/max are switched because of the negative signs
        #msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(settings.argflag['trace']['slits']['fracignore'], int(edgearr.shape[0]*binby)))
        msgs.info("Ignoring any slit that spans < {0:3.2f}x{1:d} pixels on the detector".format(settings.argflag['trace']['slits']['fracignore'], int(edgearr.shape[0])))
        fracpix = int(settings.argflag['trace']['slits']['fracignore']*edgearr.shape[0])
        lnc, lxc, rnc, rxc, ldarr, rdarr = arcytrace.ignore_orders(edgearr, fracpix, lmin, lmax, rmin, rmax)
        lmin += lnc
        rmin += rnc
        lmax -= lxc
        rmax -= rxc
        iterate = False
        if settings.argflag['trace']['slits']['number'] == 1:  # Another check on slits for singleSlit
            if lmax < lmin:
                msgs.warn("Unable to find a left edge2. Adding one in.")
                iterate = True
                edgearr[:,0] = -2*ednum
                lcnt = 1
            if rmax < rmin:
                msgs.warn("Unable to find a right edge2. Adding one in.")
                iterate = True
                edgearr[:,-1] = 2*ednum
                rcnt = 1
    # Trace left slit edges
    # First, determine the model for the most common left slit edge
    wcm = np.where(edgearr < 0)
    cntr = Counter(edg for edg in edgearr[wcm])
    commn = cntr.most_common(1)
    wedx, wedy = np.where(edgearr == commn[0][0])
    msk, cf = arutils.robust_polyfit(wedx, wedy,
                                     settings.argflag['trace']['slits']['polyorder'],
                                     function=settings.argflag['trace']['slits']['function'],
                                     minv=0, maxv=binarr.shape[0]-1)
    cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                               settings.argflag['trace']['slits']['function'],
                               minv=0, maxv=binarr.shape[0]-1)

    msgs.info("Fitting left slit traces")
    lcoeff = np.zeros((1+settings.argflag['trace']['slits']['polyorder'], lmax-lmin+1))
    ldiffarr = np.zeros(lmax-lmin+1)
    lwghtarr = np.zeros(lmax-lmin+1)
    lnmbrarr = np.zeros(lmax-lmin+1)
    offs = cenmodl[int(binarr.shape[0]/2)]
#    lfail = np.array([])
#    minvf, maxvf = slf._pixlocn[det-1][0, 0, 0], slf._pixlocn[det-1][-1, 0, 0]
    for i in range(lmin, lmax+1):
        w = np.where(edgearr == -i)
        if np.size(w[0]) <= settings.argflag['trace']['slits']['polyorder']+2:
            # lfail = np.append(lfail,i-lmin)
            continue
        tlfitx = plxbin[w]
        tlfity = plybin[w]
        ldiffarr[i-lmin] = np.mean(w[1]-cenmodl[w[0]]) + offs
        lwghtarr[i-lmin] = np.size(w[0])/float(binarr.shape[0])
        lnmbrarr[i-lmin] = -i
        #lcoeff[:, i-lmin] = arutils.func_fit(tlfitx, tlfity, settings.argflag['trace']['slits']['function'],
        #                                     settings.argflag['trace']['slits']['polyorder'], minv=minvf, maxv=maxvf)
        msk, lcoeff[:, i-lmin] = arutils.robust_polyfit(tlfitx, tlfity,
                                                        settings.argflag['trace']['slits']['polyorder'],
                                                        function=settings.argflag['trace']['slits']['function'],
                                                        minv=minvf, maxv=maxvf)
#		xv=np.linspace(0,edgearr.shape[0])
#		yv=np.polyval(coeffl[i-lmin,:],xv)
#		plt.plot(w[0],w[1],'ro')
#		plt.plot(xv,yv,'r-')
#		plt.show()
#		plt.clf()
    # Trace right slit edges
    # First, determine the model for the most common right slit edge
    wcm = np.where(edgearr > 0)
    cntr = Counter(edg for edg in edgearr[wcm])
    commn = cntr.most_common(1)
    wedx, wedy = np.where(edgearr == commn[0][0])
    msk, cf = arutils.robust_polyfit(wedx, wedy,
                                     settings.argflag['trace']['slits']['polyorder'],
                                     function=settings.argflag['trace']['slits']['function'],
                                     minv=0, maxv=binarr.shape[0]-1)
    cenmodl = arutils.func_val(cf, np.arange(binarr.shape[0]),
                               settings.argflag['trace']['slits']['function'],
                               minv=0, maxv=binarr.shape[0]-1)

    msgs.info("Fitting right slit traces")
    rcoeff = np.zeros((1+settings.argflag['trace']['slits']['polyorder'], rmax-rmin+1))
    rdiffarr = np.zeros(rmax-rmin+1)
    rwghtarr = np.zeros(rmax-rmin+1)
    rnmbrarr = np.zeros(rmax-rmin+1)
    offs = cenmodl[int(binarr.shape[0]/2)]
#	rfail = np.array([])
    for i in range(rmin, rmax+1):
        w = np.where(edgearr == i)
        if np.size(w[0]) <= settings.argflag['trace']['slits']['polyorder']+2:
#			rfail = np.append(rfail, i-rmin)
            continue
        tlfitx = plxbin[w]
        tlfity = plybin[w]
        #rcoeff[:, i-rmin] = arutils.func_fit(tlfitx, tlfity, settings.argflag['trace']['slits']['function'],
        #                                     settings.argflag['trace']['slits']['polyorder'], minv=minvf, maxv=maxvf)
        rdiffarr[i-rmin] = np.mean(w[1]-cenmodl[w[0]]) + offs
        rwghtarr[i-rmin] = np.size(w[0])/float(binarr.shape[0])
        rnmbrarr[i-rmin] = i
        msk, rcoeff[:, i-rmin] = arutils.robust_polyfit(tlfitx, tlfity,
                                                       settings.argflag['trace']['slits']['polyorder'],
                                                       function=settings.argflag['trace']['slits']['function'],
                                                       minv=minvf, maxv=maxvf)
    # Check if no further work is needed (i.e. there only exists one order)
    if (lmax+1-lmin == 1) and (rmax+1-rmin == 1):
        # Just a single order has been identified (i.e. probably longslit)
        msgs.info("Only one slit was identified. Should be a longslit.")
        xint = slf._pixlocn[det-1][:, 0, 0]
        lcenint = np.zeros((mstrace.shape[0], 1))
        rcenint = np.zeros((mstrace.shape[0], 1))
        lcenint[:, 0] = arutils.func_val(lcoeff[:, 0], xint, settings.argflag['trace']['slits']['function'],
                                         minv=minvf, maxv=maxvf)
        rcenint[:, 0] = arutils.func_val(rcoeff[:, 0], xint, settings.argflag['trace']['slits']['function'],
                                         minv=minvf, maxv=maxvf)
        return lcenint, rcenint, np.zeros(1, dtype=np.bool)
    msgs.info("Synchronizing left and right slit traces")
    # Define the array of pixel values along the dispersion direction
    xv = plxbin[:, 0]
    num = (lmax-lmin)//2  # Should be an int, right?
    lval = lmin + num  # Pick an order, somewhere in between lmin and lmax
    lv = (arutils.func_val(lcoeff[:, lval-lmin], xv, settings.argflag['trace']['slits']['function'], minv=minvf, maxv=maxvf)+0.5).astype(np.int)
    if np.any(lv < 0) or np.any(lv+1 >= binarr.shape[1]):
        msgs.warn("At least one slit is poorly traced")
        msgs.info("Refer to the manual, and adjust the input trace parameters")
        msgs.error("Cannot continue without a successful trace")
    mnvalp = np.median(binarr[:, lv+1])  # Go one row above and one row below an order edge,
    mnvalm = np.median(binarr[:, lv-1])  # then see which mean value is greater.

    """
    lvp = (arutils.func_val(lcoeff[:,lval+1-lmin],xv,settings.argflag['trace']['slits']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
    edgbtwn = arcytrace.find_between(edgearr,lv,lvp,1)
    print lval, edgbtwn
    # edgbtwn is a 3 element array that determines what is between two adjacent left edges
    # edgbtwn[0] is the next right order along, from left order lval
    # edgbtwn[1] is only !=-1 when there's an order overlap.
    # edgebtwn[2] is only used when a left order is found before a right order
    if edgbtwn[0] == -1 and edgbtwn[1] == -1:
        rsub = edgbtwn[2]-(lval) # There's an order overlap
    elif edgbtwn[1] == -1: # No overlap
        rsub = edgbtwn[0]-(lval)
    else: # There's an order overlap
        rsub = edgbtwn[1]-(lval)
    """
    if mnvalp > mnvalm:
        lvp = (arutils.func_val(lcoeff[:, lval+1-lmin], xv, settings.argflag['trace']['slits']['function'],
                                minv=minvf, maxv=maxvf)+0.5).astype(np.int)
        edgbtwn = arcytrace.find_between(edgearr, lv, lvp, 1)
        # edgbtwn is a 3 element array that determines what is between two adjacent left edges
        # edgbtwn[0] is the next right order along, from left order lval
        # edgbtwn[1] is only !=-1 when there's an order overlap.
        # edgebtwn[2] is only used when a left order is found before a right order
        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2]-lval  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0]-lval
        else:  # There's an order overlap
            rsub = edgbtwn[1]-lval
    else:
        lvp = (arutils.func_val(lcoeff[:, lval-1-lmin], xv, settings.argflag['trace']['slits']['function'],
                                minv=minvf, maxv=maxvf)+0.5).astype(np.int)
        edgbtwn = arcytrace.find_between(edgearr, lvp, lv, -1)
        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2]-(lval-1)  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0]-(lval-1)
        else:  # There's an order overlap
            rsub = edgbtwn[1]-(lval-1)

#	rva = arutils.func_val(rcoeff[:,0],midval,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#	for i in range(1,rmax-rmin+1):
#		rvb = arutils.func_val(rcoeff[:,i],midval,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#		if rvb > lv and rva < lv:
#			lvp = arutils.func_val(lcoeff[:,500-lmin],midval*2.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			lvm = arutils.func_val(lcoeff[:,500-lmin],0.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			rvap = arutils.func_val(rcoeff[:,i-1],midval*2.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			rvam = arutils.func_val(rcoeff[:,i-1],0.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			rvbp = arutils.func_val(rcoeff[:,i],midval*2.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			rvbm = arutils.func_val(rcoeff[:,i],0.0,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)
#			mina = np.min([lv-rva,lvp-rvap,lvm-rvam])
#			minb = np.min([rvb-lv,rvbp-lvp,rvbm-lvm])
#			if mina <= 1.0 or minb <= 0.0:
#				msgs.error("Orders are too close or the fitting quality is too poor")
#			yva = np.arange(0.0,mina)
#			yvb = np.arange(0.0,minb)
#			xmga, ymga = np.meshgrid(xv,yva)
#			xmgb, ymgb = np.meshgrid(xv,yvb)
#			ymga += np.array([arutils.func_val(rcoeff[:,i-1],xv,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)])
#			ymgb += np.array([arutils.func_val(rcoeff[:,i],xv,settings.argflag['trace']['slits']['function'],min=xvmn,max=xvmx)])
#			xmga = xmga.flatten().astype(np.int)
#			ymga = ymga.flatten().astype(np.int)
#			xmgb = xmgb.flatten().astype(np.int)
#			ymgb = ymgb.flatten().astype(np.int)
#			meda = np.median(binarr[xmga,ymga])
#			medb = np.median(binarr[xmgb,ymgb])
#			if meda > medb:
#				rsub = rmin+i-1-500
#			else:
#				rsub = rmin+i-500
#			break
#		else:
#			rva = rvb
    msgs.info("Relabelling slit edges")
    rsub = int(round(rsub))
    if lmin < rmin-rsub:
        esub = lmin - (settings.argflag['trace']['slits']['pca']['extrapolate']['neg']+1)
    else:
        esub = (rmin-rsub) - (settings.argflag['trace']['slits']['pca']['extrapolate']['neg']+1)

    wl = np.where(edgearr < 0)
    wr = np.where(edgearr > 0)
    edgearr[wl] += esub
    edgearr[wr] -= (esub+rsub)
    lnmbrarr += esub
    rnmbrarr -= (esub+rsub)

    # Insert new rows into coefficients arrays if rsub != 0 (if orders were not labelled correctly, there will be a mismatch for the lcoeff and rcoeff)
    almin, almax = -np.max(edgearr[wl]), -np.min(edgearr[wl]) # min and max switched because left edges have negative values
    armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])
    nmord = settings.argflag['trace']['slits']['polyorder']+1
    if armin != almin:
        if armin < almin:
            lcoeff = np.append(np.zeros((nmord, almin-armin)), lcoeff, axis=1)
            ldiffarr = np.append(np.zeros(almin-armin), ldiffarr)
            lnmbrarr = np.append(np.zeros(almin-armin), lnmbrarr)
            lwghtarr = np.append(np.zeros(almin-armin), lwghtarr)
        else:
            rcoeff = np.append(np.zeros((nmord, armin-almin)), rcoeff, axis=1)
            rdiffarr = np.append(np.zeros(armin-almin), rdiffarr)
            rnmbrarr = np.append(np.zeros(armin-almin), rnmbrarr)
            rwghtarr = np.append(np.zeros(armin-almin), rwghtarr)
    if armax != almax:
        if armax < almax:
            rcoeff = np.append(rcoeff, np.zeros((nmord, almax-armax)), axis=1)
            rdiffarr = np.append(rdiffarr, np.zeros(almax-armax))
            rnmbrarr = np.append(rnmbrarr, np.zeros(almax-armax))
            rwghtarr = np.append(rwghtarr, np.zeros(almax-armax))
        else:
            lcoeff = np.append(lcoeff, np.zeros((nmord, armax-almax)), axis=1)
            ldiffarr = np.append(ldiffarr, np.zeros(armax-almax))
            lnmbrarr = np.append(lnmbrarr, np.zeros(armax-almax))
            lwghtarr = np.append(lwghtarr, np.zeros(armax-almax))

    # Now consider traces where both the left and right edges are detected
    ordunq = np.unique(edgearr)
    lunqt = ordunq[np.where(ordunq < 0)[0]]
    runqt = ordunq[np.where(ordunq > 0)[0]]
    lunq = np.arange(lunqt.min(), lunqt.max()+1)
    runq = np.arange(runqt.min(), runqt.max()+1)
    # Determine which orders are detected on both the left and right edge
    gord = np.intersect1d(-lunq, runq, assume_unique=True)
    # We need to ignore the orders labelled rfail and lfail.
    lg = np.where(np.in1d(-lunq, gord))[0]
    rg = np.where(np.in1d(runq, gord))[0]
    lgm = np.where(np.in1d(-lunq, gord, invert=True))[0]
    rgm = np.where(np.in1d(runq, gord, invert=True))[0]
    maxord = np.max(np.append(gord, np.append(-lunq[lgm], runq[rgm])))
    lcent = arutils.func_val(lcoeff[:,-lunq[lg][::-1]-1-settings.argflag['trace']['slits']['pca']['extrapolate']['neg']], xv,
                             settings.argflag['trace']['slits']['function'], minv=minvf, maxv=maxvf)
    rcent = arutils.func_val(rcoeff[:,runq[rg]-1-settings.argflag['trace']['slits']['pca']['extrapolate']['neg']], xv,
                             settings.argflag['trace']['slits']['function'], minv=minvf, maxv=maxvf)
    slitcen = 0.5*(lcent+rcent).T
    ##############
    if settings.argflag['trace']['slits']['pca']['type'] == 'order':
        #maskord = np.where((np.all(lcoeff[:,lg],axis=0)==False)|(np.all(rcoeff[:,rg],axis=0)==False))[0]
        maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
        ordsnd = np.arange(min(almin, armin), max(almax, armax)+1)
        totord = ordsnd[-1]+settings.argflag['trace']['slits']['pca']['extrapolate']['pos']
        # Identify the orders to be extrapolated during reconstruction
        extrapord = (1.0-np.in1d(np.linspace(1.0, totord, totord), gord).astype(np.int)).astype(np.bool)
        msgs.info("Performing a PCA on the order edges")
        ofit = settings.argflag['trace']['slits']['pca']['params']
        lnpc = len(ofit)-1
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        coeffs = arutils.func_fit(xv, slitcen, settings.argflag['trace']['slits']['function'],
                                  settings.argflag['trace']['slits']['polyorder'], minv=minvf, maxv=maxvf)
        for i in range(ordsnd.size):
            if i in maskord:
                coeffs = np.insert(coeffs, i, 0.0, axis=1)
                slitcen = np.insert(slitcen, i, 0.0, axis=1)
                lcent = np.insert(lcent, i, 0.0, axis=0)
                rcent = np.insert(rcent, i, 0.0, axis=0)
        xcen = xv[:, np.newaxis].repeat(ordsnd.size, axis=1)
        fitted, outpar = arpca.basis(xcen, slitcen, coeffs, lnpc, ofit, x0in=ordsnd, mask=maskord,
                                     skipx0=False, function=settings.argflag['trace']['slits']['function'])
        if not msgs._debug['no_qa']:
            arpca.pc_plot(slf, outpar, ofit, pcadesc=pcadesc)
        # Extrapolate the remaining orders requested
        orders = 1.0+np.arange(totord)
        extrap_cent, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['slits']['function'])
        # Fit a function for the difference between left and right edges.
        diff_coeff, diff_fit = arutils.polyfitter2d(rcent-lcent, mask=maskord,
                                                    order=settings.argflag['trace']['slits']['diffpolyorder'])
        # Now extrapolate the order difference
        ydet = np.linspace(0.0, 1.0, lcent.shape[0])
        ydetd = ydet[1]-ydet[0]
        lnum = ordsnd[0]-1.0
        ydet = np.append(-ydetd*np.arange(1.0, 1.0+lnum)[::-1], ydet)
        ydet = np.append(ydet, 1.0+ydetd*np.arange(1.0, 1.0+settings.argflag['trace']['slits']['pca']['extrapolate']['pos']))
        xde, yde = np.meshgrid(np.linspace(0.0, 1.0, lcent.shape[1]), ydet)
        extrap_diff = arutils.polyval2d(xde, yde, diff_coeff).T
        msgs.info("Refining the trace for reconstructed and predicted orders")
        # NOTE::  MIGHT NEED TO APPLY THE BAD PIXEL MASK HERE TO BINARR
        msgs.work("Should the bad pixel mask be applied to the frame here?")
        refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff,
                                            [gord[0]-orders[0], orders[-1]-gord[-1]], orders,
                                            ofit[0], slf._pixlocn[det-1],
                                            function=settings.argflag['trace']['slits']['function'])
        # Generate the left and right edges
        lcen = refine_cent - 0.5*extrap_diff
        rcen = refine_cent + 0.5*extrap_diff
        # lcen = extrap_cent - 0.5*extrap_diff
        # rcen = extrap_cent + 0.5*extrap_diff
    elif settings.argflag['trace']['slits']['pca']['type'] == 'pixel':
        maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
        allord = np.arange(ldiffarr.shape[0])
        ww = np.where(np.in1d(allord, maskord) == False)[0]
        # Unmask where an order edge is located
        maskrows = np.ones(binarr.shape[1], dtype=np.int)
        ldiffarr = np.round(ldiffarr[ww]).astype(np.int)
        rdiffarr = np.round(rdiffarr[ww]).astype(np.int)
        maskrows[ldiffarr] = 0
        maskrows[rdiffarr] = 0
        # Extract the slit edge ID numbers associated with the acceptable traces
        lnmbrarr = lnmbrarr[ww]
        rnmbrarr = rnmbrarr[ww]
        # Fill in left/right coefficients
        tcoeff = np.ones((settings.argflag['trace']['slits']['polyorder']+1, binarr.shape[1]))
        tcoeff[:, ldiffarr] = lcoeff[:, ww]
        tcoeff[:, rdiffarr] = rcoeff[:, ww]
        # Weight the PCA fit by the number of detections in each slit edge
        pxwght = np.zeros(binarr.shape[1])
        pxwght[ldiffarr] = lwghtarr[ww]
        pxwght[rdiffarr] = rwghtarr[ww]
        maskrw = np.where(maskrows == 1)[0]
        maskrw.sort()
        extrap_row = maskrows.copy()
        xv = np.arange(binarr.shape[0])
        # trace values
        trcval = arutils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'],
                                  minv=minvf, maxv=maxvf).T
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        ofit = settings.argflag['trace']['slits']['pca']['params']
        lnpc = len(ofit)-1
        # Only do a PCA if there are enough good slits
        if np.sum(1.0-extrap_row) > ofit[0]+1:
            # Perform a PCA on the locations of the slits
            msgs.info("Performing a PCA on the slit traces")
            ordsnd = np.arange(binarr.shape[1])
            xcen = xv[:, np.newaxis].repeat(binarr.shape[1], axis=1)
            fitted, outpar = arpca.basis(xcen, trcval, tcoeff, lnpc, ofit, weights=pxwght,
                                         x0in=ordsnd, mask=maskrw, skipx0=False,
                                         function=settings.argflag['trace']['slits']['function'])
            if not msgs._debug['no_qa']:
                arpca.pc_plot(slf, outpar, ofit, pcadesc=pcadesc, addOne=False)
            # Now extrapolate to the whole detector
            pixpos = np.arange(binarr.shape[1])
            extrap_trc, outpar = arpca.extrapolate(outpar, pixpos,
                                                   function=settings.argflag['trace']['slits']['function'])
            # Extract the resulting edge traces
            lcen = extrap_trc[:, ldiffarr]
            rcen = extrap_trc[:, rdiffarr]
            # Perform a final shift fit to ensure the traces closely follow the edge detections
            for ii in range(lnmbrarr.size):
                wedx, wedy = np.where(edgearr == lnmbrarr[ii])
                shft = np.mean(lcen[wedx, ii]-wedy)
                lcen[:, ii] -= shft
            for ii in range(rnmbrarr.size):
                wedx, wedy = np.where(edgearr == rnmbrarr[ii])
                shft = np.mean(rcen[wedx, ii]-wedy)
                rcen[:, ii] -= shft
        else:
            allord = np.arange(lcent.shape[0])
            maskord = np.where((np.all(lcent, axis=1) == False) | (np.all(rcent, axis=1) == False))[0]
            ww = np.where(np.in1d(allord, maskord) == False)[0]
            lcen = lcent[ww, :].T.copy()
            rcen = rcent[ww, :].T.copy()
        extrapord = np.zeros(lcen.shape[1], dtype=np.bool)
    else:
        allord = np.arange(lcent.shape[0])
        maskord = np.where((np.all(lcent, axis=1) == False) | (np.all(rcent, axis=1) == False))[0]
        ww = np.where(np.in1d(allord, maskord) == False)[0]
        lcen = lcent[ww, :].T.copy()
        rcen = rcent[ww, :].T.copy()
        extrapord = np.zeros(lcen.shape[1], dtype=np.bool)
    # Interpolate the best solutions for all orders with a cubic spline
    #msgs.info("Interpolating the best solutions for all orders with a cubic spline")
    # zmin, zmax = arplot.zscale(binarr)
    # plt.clf()
    # extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
    # implot = plt.imshow(binarr, extent=extnt, origin='lower', interpolation='none', aspect='auto')
    # implot.set_cmap("gray")
    # plt.colorbar()
    # implot.set_clim(zmin, zmax)
    #xint = slf._pixlocn[det-1][:,0,0]
    #lcenint = np.zeros((mstrace.shape[0], lcen.shape[1]))
    #rcenint = np.zeros((mstrace.shape[0], rcen.shape[1]))
    lcenint = lcen.copy()
    rcenint = rcen.copy()
    #for i in range(lcen.shape[1]):
        # if i+1 not in gord: pclr = '--'
        # else: pclr = '-'
        # plt.plot(xv_corr, lcen_corr[:,i], 'g'+pclr, linewidth=0.1)
        # plt.plot(xv_corr, rcen_corr[:,i], 'b'+pclr, linewidth=0.1)
        #lcenint[:,i] = arutils.spline_interp(xint, xv, lcen[:,i])
        #rcenint[:,i] = arutils.spline_interp(xint, xv, rcen[:,i])
#       plt.plot(xint,lcenint[:,i],'k.')
#       plt.plot(xint,rcenint[:,i],'k.')
#       plt.show()
    # if tracedesc != "":
    #     plt.title(tracedesc)
    # slf._qa.savefig(dpi=1200, orientation='portrait', bbox_inches='tight')
    # plt.close()

    # Illustrate where the orders fall on the detector (physical units)
    if msgs._debug['trace']:
        from pypit import ginga
        viewer, ch = ginga.show_image(mstrace)
        ginga.show_slits(viewer, ch, lcenint, rcenint)
        debugger.set_trace()
    if settings.argflag['run']['qa']:
        msgs.work("Not yet setup with ginga")
    return lcenint, rcenint, extrapord


def refine_traces(binarr, outpar, extrap_cent, extrap_diff, extord, orders,
                  fitord, locations, function='polynomial'):
    """
    Parameters
    ----------
    binarr
    outpar
    extrap_cent
    extrap_diff
    extord
    orders
    fitord
    locations
    function

    Returns
    -------

    """
    from pypit import arcytrace
    # Refine the orders in the positive direction
    i = extord[1]
    hiord = phys_to_pix(extrap_cent[:,-i-2], locations, 1)
    nxord = phys_to_pix(extrap_cent[:,-i-1], locations, 1)
    mask = np.ones(orders.size)
    mask[0:extord[0]] = 0.0
    mask[-extord[1]:] = 0.0
    extfit = extrap_cent.copy()
    outparcopy = copy.deepcopy(outpar)
    while i > 0:
        loord = hiord
        hiord = nxord
        nxord = phys_to_pix(extrap_cent[:,-i], locations, 1)
        minarrL = arcytrace.minbetween(binarr, loord, hiord) # Minimum counts between loord and hiord
        minarrR = arcytrace.minbetween(binarr, hiord, nxord)
        minarr = 0.5*(minarrL+minarrR)
        srchz = np.abs(extfit[:,-i]-extfit[:,-i-1])/3.0
        lopos = phys_to_pix(extfit[:,-i]-srchz, locations, 1) # The pixel indices for the bottom of the search window
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,-i])))
        diffarr = np.round(extrap_diff[:,-i]).astype(np.int)
        shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch)
        relshift = np.mean(shift+extrap_diff[:,-i]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[-i])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[-i]), relshift))
        # Renew guess for the next order
        mask[-i] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, -i, relshift, fitord, function=function)
        if fail:
            msgs.warn("Order refinement has large residuals -- check order traces")
            return extrap_cent, outparcopy
        i -= 1
    # Refine the orders in the negative direction
    i = extord[0]
    loord = phys_to_pix(extrap_cent[:,i+1], locations, 1)
    extrap_cent = extfit.copy()
    outparcopy = copy.deepcopy(outpar)
    while i > 0:
        hiord = loord
        loord = phys_to_pix(extfit[:,i], locations, 1)
        minarr = arcytrace.minbetween(binarr,loord, hiord)
        srchz = np.abs(extfit[:,i]-extfit[:,i-1])/3.0
        lopos = phys_to_pix(extfit[:,i-1]-srchz, locations, 1)
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,i-1])))
        diffarr = np.round(extrap_diff[:,i-1]).astype(np.int)
        shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch)
        relshift = np.mean(shift+extrap_diff[:,i-1]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[i-1])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[i-1]),relshift))
        # Renew guess for the next order
        mask[i-1] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, i-1, relshift, fitord, function=function)
        if fail:
            msgs.warn("Order refinement has large residuals -- check order traces")
            return extrap_cent, outparcopy
        i -= 1
    return extfit, outpar


def trace_tilt(slf, det, msarc, slitnum, censpec=None, maskval=-999999.9,
               trthrsh=1000.0, nsmth=0, method = "fweight"):
    """
    This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
               trthrsh=1000.0, nsmth=0):

    Parameters
    ----------
    slf
    det
    msarc
    slitnum
    censpec
    maskval
    trthrsh
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

    # from pypit import arcyutils
    dnum = settings.get_dnum(det)

    msgs.work("Detecting lines for slit {0:d}".format(slitnum+1))
    ordcen = slf._pixcen[det-1].copy()
    tampl, tcent, twid, w, satsnd, _ = ararc.detect_lines(slf, det, msarc, censpec=censpec)
    satval = settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear']
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
    # Restricted to ID lines? [avoid ghosts]
    if settings.argflag['trace']['slits']['tilts']['idsonly']:
        ids_pix = np.round(np.array(slf._wvcalib[det-1]['xfit'])*(msarc.shape[0]-1))
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
        sz = int(np.floor(np.abs(slf._rordloc[det-1][arcdet[j], slitnum]-slf._lordloc[det-1][arcdet[j], slitnum])/2.0)) - 2
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
            if len(yfit.shape) == 2:
                yfit = np.median(yfit, axis=1)
            if np.size(yfit) == 0:
                offchip = True
                break
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


def trace_fweight(fimage, xinit, invvar=None, radius=3.):
    '''Python port of trace_fweight.pro from IDLUTILS

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
    '''
    # Definitions for Cython
    #cdef int nx,ny,ncen

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

    fullpix = int(np.maximum(np.min(ix2-ix1)-1,0))
    sumw = np.zeros(ny)
    sumxw = np.zeros(ny)
    sumwt = np.zeros(ny)
    sumsx1 = np.zeros(ny)
    sumsx2 = np.zeros(ny)
    qbad = np.array([False]*ny) 

    if invvar is None: 
        invvar = np.zeros_like(fimage) + 1. 

    '''
    cdef np.ndarray[ITYPE_t, ndim=1] ycen = np.arange(ny, dtype=ITYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] invvar = 0. * fimage + 1.
    cdef np.ndarray[DTYPE_t, ndim=1] x1 = xinit - radius + 0.5
    cdef np.ndarray[DTYPE_t, ndim=1] x2 = xinit + radius + 0.5
    cdef np.ndarray[ITYPE_t, ndim=1] ix1 = np.fix(x1)
    cdef np.ndarray[ITYPE_t, ndim=1] ix2 = np.fix(x2)
    cdef np.ndarray[DTYPE_t, ndim=1] fullpix = np.maximum(np.min(ix2-ix1)-1),0)

    cdef np.ndarray[DTYPE_t, ndim=1] sumw = np.zeros(ny, dypte=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] sumxw = np.zeros(ny, dypte=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] sumwt = np.zeros(ny, dypte=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] sumsx1 = np.zeros(ny, dypte=DTYPE)
    cdef np.ndarray[DTYPE_t, ndim=1] sumsx2 = np.zeros(ny, dypte=DTYPE)
    cdef np.ndarray[ITYPE_t, ndim=1] qbad = np.zeros(ny, dypte=ITYPE)
    '''

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
    good = (sumw > 0) &  (~qbad)
    if np.sum(good) > 0:
        delta_x = sumxw[good]/sumw[good]
        xnew[good] = delta_x + xinit[good]
        xerr[good] = np.sqrt(sumsx1[good] + sumsx2[good]*delta_x**2)/sumw[good]

    bad = np.any([np.abs(xnew-xinit) > radius + 0.5,xinit < radius - 0.5,xinit > nx - 0.5 - radius],axis=0)
    if np.sum(bad) > 0:
        xnew[bad] = xinit[bad]
        xerr[bad] = 999.0

    # Return
    return xnew, xerr


def echelle_tilt(slf, msarc, det, pcadesc="PCA trace of the spectral tilts", maskval=-999999.9):
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
    from pypit import arcytrace

    arccen, maskslit, satmask = get_censpec(slf, msarc, det, gen_satmask=True)
    # If the user sets no tilts, return here
    if settings.argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satmask, None
    norders = maskslit.size

    # Now model the tilt for each slit
    tiltang, centval = None, None
    for o in range(norders):
        if maskslit[o] == 1:
            continue
        # Determine the tilts for this slit
        trcdict = trace_tilt(slf, det, msarc, o, censpec=arccen[:, o])
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
            if wmask.size < settings.argflag['trace']['slits']['tilts']['order']+2:
                maskslit[o] = 1
                continue
            null, tcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask],
                                                  settings.argflag['trace']['slits']['tilts']['order'], sigma=2.0)
            # Save the tilt angle
            tiltang[j, o] = tcoeff[1]  # tan(tilt angle)
            centval[j, o] = tcoeff[0]  # centroid of arc line

    msgs.info("Fitting tilt angles")
    tcoeff = np.ones((settings.argflag['trace']['slits']['tilts']['disporder'] + 1, tiltang.shape[1]))
    maskord = np.where(maskslit == 1)[0]
    extrap_ord = np.zeros(norders)
    for o in range(norders):
        if o in maskord:
            extrap_ord[o] = 1
            continue
        w = np.where(tiltang[:, o] != maskval)
        if np.size(w[0]) <= settings.argflag['trace']['slits']['tilts']['disporder'] + 2:
            extrap_ord[o] = 1.0
            maskord = np.append(maskord, o)
        else:
            null, tempc = arutils.robust_polyfit(centval[:, o][w], tiltang[:, o][w],
                                                 settings.argflag['trace']['slits']['tilts']['disporder'],
                                                 function=settings.argflag['trace']['slits']['function'], sigma=2.0,
                                                 minv=0.0, maxv=msarc.shape[0] - 1)
            tcoeff[:, o] = tempc
    # Sort which orders are masked
    maskord.sort()
    xv = np.arange(msarc.shape[0])
    tiltval = arutils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'],
                               minv=0.0, maxv=msarc.shape[0] - 1).T
    ofit = settings.argflag['trace']['slits']['tilts']['params']
    lnpc = len(ofit) - 1
    if np.sum(1.0 - extrap_ord) > ofit[0] + 1:  # Only do a PCA if there are enough good orders
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the spectral tilts")
        ordsnd = np.arange(norders) + 1.0
        xcen = xv[:, np.newaxis].repeat(norders, axis=1)
        fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=False,
                                     function=settings.argflag['trace']['slits']['function'])
        if not msgs._debug['no_qa']:
            arpca.pc_plot(slf, outpar, ofit, pcadesc=pcadesc)
        # Extrapolate the remaining orders requested
        orders = 1.0 + np.arange(norders)
        extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['slits']['function'])
        tilts = extrap_tilt
        if not msgs._debug['no_qa']:
            arpca.pc_plot_arctilt(slf, tiltang, centval, tilts)
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
        if np.size(xtiltfit) > settings.argflag['trace']['slits']['tilts']['disporder'] + 2:
            tcoeff = arutils.func_fit(xtiltfit, ytiltfit, settings.argflag['trace']['slits']['function'],
                                      settings.argflag['trace']['slits']['tilts']['disporder'],
                                      minv=0.0, maxv=msarc.shape[0] - 1)
            tiltval = arutils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'], minv=0.0,
                                       maxv=msarc.shape[0] - 1)
            tilts = tiltval[:, np.newaxis].repeat(tiltang.shape[1], axis=1)
        else:
            msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
            msgs.info("Assuming there is no tilt")
            tilts = np.zeros_like(slf._lordloc)

    # Generate tilts image
    tiltsimg = arcytrace.tilts_image(tilts, slf._lordloc[det-1], slf._rordloc[det-1], msarc.shape[1])
    return tiltsimg, satmask, outpar


def multislit_tilt(slf, msarc, det, maskval=-999999.9):
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

    Returns
    -------
    tilts : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    satmask : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    arccen, maskslit, satmask = get_censpec(slf, msarc, det, gen_satmask=True)
    # If the user sets no tilts, return here
    if settings.argflag['trace']['slits']['tilts']['method'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satmask, None

    ordcen = slf._pixcen[det - 1].copy()
    fitxy = [settings.argflag['trace']['slits']['tilts']['order'], 1]

    # Now trace the tilt for each slit
    for o in range(arccen.shape[1]):
        # Determine the tilts for this slit
        trcdict = trace_tilt(slf, det, msarc, o, censpec=arccen[:, o], nsmth=3)
        if trcdict is None:
            # No arc lines were available to determine the spectral tilt
            continue
        if msgs._debug['tilts']:
            from pypit import ginga
            ginga.chk_arc_tilts(msarc, trcdict, sedges=(slf._lordloc[det-1][:,o], slf._rordloc[det-1][:,o]))
            #debugger.set_trace()
        # Extract information from the trace dictionary
        aduse = trcdict["aduse"]
        arcdet = trcdict["arcdet"]
        xtfits = trcdict["xtfit"]
        ytfits = trcdict["ytfit"]
        wmasks = trcdict["wmask"]
        badlines = trcdict["badlines"]
        # Initialize some arrays
        maskrows = np.ones(msarc.shape[0], dtype=np.int)
        tcoeff = np.ones((settings.argflag['trace']['slits']['tilts']['order'] + 1, msarc.shape[0]))
        xtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        ytilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        ztilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        mtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        wtilt = np.ones((msarc.shape[1], arcdet.size)) * maskval
        # Analyze each spectral line
        for j in range(arcdet.size):
            if not aduse[j]:
                continue
            xtfit = xtfits[j]
            ytfit = ytfits[j]
            wmask = wmasks[j]
            xint = int(xtfit[0])
            sz = (xtfit.size-1)/2

            # Perform a scanning polynomial fit to the tilts
            # model = arcyutils.polyfit_scan_intext(xtfit, ytfit, np.ones(ytfit.size, dtype=np.float), mtfit,
            #                                       2, sz/6, 3, maskval)
            wmfit = np.where(ytfit != maskval)
            if wmfit[0].size > settings.argflag['trace']['slits']['tilts']['order'] + 1:
                cmfit = arutils.func_fit(xtfit[wmfit], ytfit[wmfit], settings.argflag['trace']['slits']['function'],
                                         settings.argflag['trace']['slits']['tilts']['order'],
                                         minv=0.0, maxv=msarc.shape[1] - 1.0)
                model = arutils.func_val(cmfit, xtfit, settings.argflag['trace']['slits']['function'],
                                         minv=0.0, maxv=msarc.shape[1] - 1.0)
            else:
                aduse[j] = False
                badlines += 1
                continue

            if maskval in model:
                # Model contains masked values
                aduse[j] = False
                badlines += 1
                continue

            # Perform a robust polynomial fit to the traces
            if settings.argflag['trace']['slits']['tilts']['method'].lower() == "spca":
                yfit = ytfit[wmask] / (msarc.shape[0] - 1.0)
            else:
                yfit = (2.0 * model[sz] - ytfit[wmask]) / (msarc.shape[0] - 1.0)
            wmsk, mcoeff = arutils.robust_polyfit(xtfit[wmask], yfit,
                                                  settings.argflag['trace']['slits']['tilts']['order'],
                                                  function=settings.argflag['trace']['slits']['function'],
                                                  sigma=2.0, minv=0.0, maxv=msarc.shape[1] - 1.0)
            # Update the mask
            wmask = wmask[np.where(wmsk == 0)]

            # Save the tilt angle, and unmask the row
            factr = (msarc.shape[0] - 1.0) * arutils.func_val(mcoeff, ordcen[arcdet[j], 0],
                                                              settings.argflag['trace']['slits']['function'],
                                                              minv=0.0, maxv=msarc.shape[1] - 1.0)
            idx = int(factr + 0.5)
            if (idx > 0) and (idx < msarc.shape[0]):
                maskrows[idx] = 0
                tcoeff[:, idx] = mcoeff.copy()
            # Restrict to good IDs?
            if settings.argflag['trace']['slits']['tilts']['idsonly']:
                if not aduse[j]:
                    maskrows[idx] = 1

            xtilt[xint:xint + 2 * sz + 1, j] = xtfit / (msarc.shape[1] - 1.0)
            ytilt[xint:xint + 2 * sz + 1, j] = arcdet[j] / (msarc.shape[0] - 1.0)
            ztilt[xint:xint + 2 * sz + 1, j] = ytfit / (msarc.shape[0] - 1.0)
            if settings.argflag['trace']['slits']['tilts']['method'].lower() == "spline":
                mtilt[xint:xint + 2 * sz + 1, j] = model / (msarc.shape[0] - 1.0)
            elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "interp":
                mtilt[xint:xint + 2 * sz + 1, j] = (2.0 * model[sz] - model) / (msarc.shape[0] - 1.0)
            else:
                mtilt[xint:xint + 2 * sz + 1, j] = (2.0 * model[sz] - model) / (msarc.shape[0] - 1.0)
            wbad = np.where(ytfit == maskval)[0]
            ztilt[xint + wbad, j] = maskval
            if wmask.size != 0:
                sigg = max(1.4826 * np.median(np.abs(ytfit - model)[wmask]) / np.sqrt(2.0), 1.0)
                wtilt[xint:xint + 2 * sz + 1, j] = 1.0 / sigg
            # Extrapolate off the slit to the edges of the chip
            nfit = 6  # Number of pixels to fit a linear function to at the end of each trace
            xlof, xhif = np.arange(xint, xint + nfit), np.arange(xint + 2 * sz + 1 - nfit, xint + 2 * sz + 1)
            xlo, xhi = np.arange(xint), np.arange(xint + 2 * sz + 1, msarc.shape[1])
            glon = np.mean(xlof * mtilt[xint:xint + nfit, j]) - np.mean(xlof) * np.mean(mtilt[xint:xint + nfit, j])
            glod = np.mean(xlof ** 2) - np.mean(xlof) ** 2
            clo = np.mean(mtilt[xint:xint + nfit, j]) - (glon / glod) * np.mean(xlof)
            yhi = mtilt[xint + 2 * sz + 1 - nfit:xint + 2 * sz + 1, j]
            ghin = np.mean(xhif * yhi) - np.mean(xhif) * np.mean(yhi)
            ghid = np.mean(xhif ** 2) - np.mean(xhif) ** 2
            chi = np.mean(yhi) - (ghin / ghid) * np.mean(xhif)
            mtilt[0:xint, j] = (glon / glod) * xlo + clo
            mtilt[xint + 2 * sz + 1:, j] = (ghin / ghid) * xhi + chi
        if badlines != 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")

        # Masking
        maskrw = np.where(maskrows == 1)[0]
        maskrw.sort()
        extrap_row = maskrows.copy()
        xv = np.arange(msarc.shape[1])
        # Tilt values
        tiltval = arutils.func_val(tcoeff, xv, settings.argflag['trace']['slits']['function'],
                                   minv=0.0, maxv=msarc.shape[1] - 1.0).T
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        ofit = settings.argflag['trace']['slits']['tilts']['params']
        lnpc = len(ofit) - 1
        # Only do a PCA if there are enough good orders
        if np.sum(1.0 - extrap_row) > ofit[0] + 1:
            # Perform a PCA on the tilts
            msgs.info("Performing a PCA on the tilts")
            ordsnd = np.linspace(0.0, 1.0, msarc.shape[0])
            xcen = xv[:, np.newaxis].repeat(msarc.shape[0], axis=1)
            fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, weights=None,
                                         x0in=ordsnd, mask=maskrw, skipx0=False,
                                         function=settings.argflag['trace']['slits']['function'])
            if not msgs._debug['no_qa']:
                arpca.pc_plot(slf, outpar, ofit, pcadesc="Spectral Tilts PCA", addOne=False)
            # Extrapolate the remaining orders requested
            orders = np.linspace(0.0, 1.0, msarc.shape[0])
            extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=settings.argflag['trace']['slits']['function'])
            polytilts = extrap_tilt.T
        else:
            # Fit the model with a 2D polynomial
            msgs.warn("Could not perform a PCA when tracing the spectral tilt" + msgs.newline() +
                      "Not enough well-traced arc lines")
            msgs.info("Fitting tilts with a low order, 2D polynomial")
            wgd = np.where(xtilt != maskval)
            coeff = arutils.polyfit2d_general(xtilt[wgd], ytilt[wgd], mtilt[wgd], fitxy)
            polytilts = arutils.polyval2d_general(coeff, np.linspace(0.0, 1.0, msarc.shape[1]),
                                                  np.linspace(0.0, 1.0, msarc.shape[0]))

        if settings.argflag['trace']['slits']['tilts']['method'].lower() == "interp":
            msgs.info("Interpolating and Extrapolating the tilts")
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # yspl = np.append(0.0, np.append(arcdet[np.where(aduse)]/(msarc.shape[0]-1.0), 1.0))
            # yspl = np.append(0.0, np.append(polytilts[arcdet[np.where(aduse)], msarc.shape[1]/2], 1.0))
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)]]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
            zspl[:, 1:-1] = mtilt[:, np.where(aduse)[0]]
            # zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
            zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
            zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
            # Make sure the endpoints are set to 0.0 and 1.0
            zspl[:, 0] -= zspl[ordcen[0, 0], 0]
            zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, 0], -1] + 1.0
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
            tiltspl = interp.RectBivariateSpline(xspl, yspl, zspl, kx=3, ky=3)
            yval = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yval, grid=True).T
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "spline":
            msgs.info("Performing a spline fit to the tilts")
            wgd = np.where((ytilt != maskval) & (ztilt != maskval))
            txsbs = xtilt[wgd]
            tysbs = ytilt[wgd]
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
            tiltspl = interp.SmoothBivariateSpline(xsbs, zsbs, ysbs, w=wsbs, kx=3, ky=3, s=xsbs.size,
                                                   bbox=[0.0, 1.0, min(zsbs.min(), 0.0), max(zsbs.max(), 1.0)])
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            yspl = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yspl, grid=True).T
            # QA
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
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "spca":
            # Slit position
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # Trace positions down center of the order
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)]]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            # Trace positions as measured+modeled
            zspl = np.zeros((msarc.shape[1], np.sum(aduse) + 2))
            zspl[:, 1:-1] = polytilts[arcdet[np.where(aduse)[0]], :].T
            zspl[:, 0] = zspl[:, 1] + polytilts[0, :] - polytilts[arcdet[np.where(aduse)[0][0]], :]
            zspl[:, -1] = zspl[:, -2] + polytilts[-1, :] - polytilts[arcdet[np.where(aduse)[0][-1]], :]
            # Make sure the endpoints are set to 0.0 and 1.0
            zspl[:, 0] -= zspl[ordcen[0, 0], 0]
            zspl[:, -1] = zspl[:, -1] - zspl[ordcen[-1, 0], -1] + 1.0
            # Prepare the spline variables
            if False:
                pmin = 0
                pmax = -1
            else:
                pmin = int(max(0, np.min(slf._lordloc[det - 1])))
                pmax = int(min(msarc.shape[1], np.max(slf._rordloc[det - 1])))
            xsbs = np.outer(xspl, np.ones(yspl.size))[pmin:pmax, :]
            ysbs = np.outer(np.ones(xspl.size), yspl)[pmin:pmax, :]
            zsbs = zspl[pmin:pmax, :]
            # Spline
            msgs.work('Consider adding weights to SmoothBivariate in spca')
            tiltspl = interp.SmoothBivariateSpline(xsbs.flatten(),
                                                   zsbs.flatten(),
                                                   ysbs.flatten(), kx=3, ky=3, s=xsbs.size)
            # Finish
            yval = np.linspace(0.0, 1.0, msarc.shape[0])
            tilts = tiltspl(xspl, yval, grid=True).T
            if False:
                tiltqa = tiltspl(xsbs.flatten(), zsbs.flatten(), grid=False).reshape(xsbs.shape)
                plt.clf()
                # plt.imshow((zsbs-tiltqa)/zsbs, origin='lower')
                plt.imshow((ysbs - tiltqa) / ysbs, origin='lower')
                plt.colorbar()
                plt.show()
                debugger.set_trace()
        elif settings.argflag['trace']['slits']['tilts']['method'].lower() == "pca":
            tilts = polytilts.copy()


    # Now do the QA
    msgs.info("Preparing arc tilt QA data")
    tiltsplot = tilts[arcdet, :].T
    tiltsplot *= (msarc.shape[0] - 1.0)
    # Shift the plotted tilts about the centre of the slit
    ztilto = ztilt.copy()
    adj = np.diag(tilts[arcdet, ordcen[arcdet]])
    zmsk = np.where(ztilto == maskval)
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

    msgs.info("Plotting arc tilt QA")
    if not msgs._debug['no_qa']:
        arqa.plot_orderfits(slf, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
                        textplt="Arc line", maxp=9, desc="Arc line spectral tilts", maskval=maskval)
    return tilts, satmask, outpar


def get_censpec(slf, frame, det, gen_satmask=False):
    """
    The value of "tilts" returned by this function is of the form:
    tilts = tan(tilt angle), where "tilt angle" is the angle between
    (1) the line representing constant wavelength and
    (2) the column of pixels that is most closely parallel with the spatial direction of the slit.

    The angle is determined relative to the axis defined by ...

    In other words, tilts = y/x according to the docs/get_locations_orderlength.JPG file.

    """
    from pypit import arcyarc
    dnum = settings.get_dnum(det)

    ordcen = 0.5*(slf._lordloc[det-1]+slf._rordloc[det-1])
    ordwid = 0.5*np.abs(slf._lordloc[det-1]-slf._rordloc[det-1])
    if gen_satmask:
        msgs.info("Generating a mask of arc line saturation streaks")
        satmask = arcyarc.saturation_mask(frame, settings.spect[dnum]['saturation']*settings.spect[dnum]['nonlinear'])
        satsnd = arcyarc.order_saturation(satmask, (ordcen+0.5).astype(np.int), (ordwid+0.5).astype(np.int))
    # Extract a rough spectrum of the arc in each slit
    msgs.info("Extracting an approximate arc spectrum at the centre of each slit")
    tordcen = None
    maskslit = np.zeros(ordcen.shape[1], dtype=np.int)
    for i in range(ordcen.shape[1]):
        wl = np.size(np.where(ordcen[:, i] <= slf._pixlocn[det-1][0,0,1])[0])
        wh = np.size(np.where(ordcen[:, i] >= slf._pixlocn[det-1][0,-1,1])[0])
        if wl == 0 and wh == 0:
            # The center of the slit is always on the chip
            if tordcen is None:
                tordcen = np.zeros((ordcen.shape[0], 1), dtype=np.float)
                tordcen[:, 0] = ordcen[:, i]
            else:
                tordcen = np.append(tordcen, ordcen[:, i].reshape((ordcen.shape[0], 1)), axis=1)
        else:
            # A slit isn't always on the chip
            if tordcen is None:
                tordcen = np.zeros((ordcen.shape[0], 1), dtype=np.float)
            else:
                tordcen = np.append(tordcen, ordcen[:, i].reshape((ordcen.shape[0], 1)), axis=1)
            maskslit[i] = 1
    w = np.where(maskslit == 0)[0]
    if tordcen is None:
        msgs.warn("Could not determine which slits are fully on the detector")
        msgs.info("Assuming all slits are fully on the detector")
        ordcen = phys_to_pix(ordcen, slf._pixlocn[det-1], 1)
    else:
        ordcen = phys_to_pix(tordcen[:,w], slf._pixlocn[det-1], 1)

    pixcen = np.arange(frame.shape[0])
    temparr = pixcen.reshape(frame.shape[0], 1).repeat(ordcen.shape[1], axis=1)
    # Average over three pixels to remove some random fluctuations, and increase S/N
    op1 = ordcen+1
    op2 = ordcen+2
    om1 = ordcen-1
    om2 = ordcen-2
    w = np.where(om1 < 0)
    om1[w] += 1
    w = np.where(om2 == -1)
    om2[w] += 1
    w = np.where(om2 == -2)
    om2[w] += 2
    w = np.where(op1 >= frame.shape[1])
    op1[w] -= 1
    w = np.where(op2 == frame.shape[1])
    op2[w] -= 1
    w = np.where(op2 == frame.shape[1]+1)
    op2[w] -= 2
    arccen = (1.0/5.0) * (frame[temparr, ordcen] +
              frame[temparr, op1] + frame[temparr, op2] +
              frame[temparr, om1] + frame[temparr, om2])
    del temparr

    if gen_satmask:
        return arccen, maskslit, satsnd
    else:
        return arccen, maskslit


def gen_pixloc(frame, det, gen=True):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    frame : ndarray
      uniformly illuminated and normalized flat field frame
    det : int
      Index of the detector

    Returns
    -------
    locations : ndarray
      A 3D array containing the x center, y center, x width and y width of each pixel.
      The returned array has a shape:   frame.shape + (4,)
    """
    dnum = settings.get_dnum(det)
    msgs.info("Deriving physical pixel locations on the detector")
    locations = np.zeros((frame.shape[0],frame.shape[1],4))
    if gen:
        msgs.info("Pixel gap in the dispersion direction = {0:4.3f}".format(settings.spect[dnum]['xgap']))
        msgs.info("Pixel size in the dispersion direction = {0:4.3f}".format(1.0))
        xs = np.arange(frame.shape[0]*1.0)*settings.spect[dnum]['xgap']
        xt = 0.5 + np.arange(frame.shape[0]*1.0) + xs
        msgs.info("Pixel gap in the spatial direction = {0:4.3f}".format(settings.spect[dnum]['ygap']))
        msgs.info("Pixel size in the spatial direction = {0:4.3f}".format(settings.spect[dnum]['ysize']))
        ys = np.arange(frame.shape[1])*settings.spect[dnum]['ygap']*settings.spect[dnum]['ysize']
        yt = settings.spect[dnum]['ysize']*(0.5 + np.arange(frame.shape[1]*1.0)) + ys
        xloc, yloc = np.meshgrid(xt, yt)
#		xwid, ywid = np.meshgrid(xs,ys)
        msgs.info("Saving pixel locations")
        locations[:,:,0] = xloc.T
        locations[:,:,1] = yloc.T
        locations[:,:,2] = 1.0
        locations[:,:,3] = settings.spect[dnum]['ysize']
    else:
        msgs.error("Have not yet included an algorithm to automatically generate pixel locations")
    return locations


def phys_to_pix(array, pixlocn, axis):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    array : ndarray
      An array of physical pixel locations
    pixlocn : ndarray
      A 3D array containing the x center, y center, x width and y width of each pixel.
    axis : int
      The axis that the input array probes

    Returns
    -------
    pixarr : ndarray
      The pixel locations of the input array (as seen on a computer screen)
    """
    from pypit import arcytrace

    if axis == 0:
        diff = pixlocn[:,0,0]
    else:
        diff = pixlocn[0,:,1]
    if len(np.shape(array)) == 1:
        pixarr = arcytrace.phys_to_pix(np.array([array]).T, diff).flatten()
    else:
        pixarr = arcytrace.phys_to_pix(array, diff)
    return pixarr


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
    xtrc = np.round(scitrace['traces'][:,obj]).astype(int)
    msgs.work("Use 2D spline to evaluate tilts")
    trc_tilt = tilts[np.arange(tilts.shape[0]), xtrc]
    trc_tilt_img = np.outer(trc_tilt, np.ones(tilts.shape[1]))
    # Slit image
    msgs.work("Should worry about changing plate scale")
    dy = (tilts - trc_tilt_img)/dypix  # Pixels
    dx = ximg - np.outer(scitrace['traces'][:,obj],np.ones(tilts.shape[1]))
    slit_img = np.sqrt(dx**2 - dy**2)
    neg = dx < 0.
    slit_img[neg] *= -1
    # Return
    return slit_img
