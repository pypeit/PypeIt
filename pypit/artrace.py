import numpy as np
import os
import copy
import arqa
import ararc
import arcyarc
import arcytrace
import arcyutils
import armsgs
import arutils
import arpca
import arplot
import matplotlib.pyplot as plt
import scipy.interpolate as interp
import scipy.ndimage as ndimage

# Logging
msgs = armsgs.get_logger()

try:
    from xastropy.xutils import xdebug as debugger
except:
    import pdb as debugger


def dispdir(msframe, dispwin=None, mode=0):
    """
    Estimate which axis is predominantly the dispersion direction of the data

    Parameters
    ----------
    msframe : ndarray
      Master calibration frame used to estimate the dispersion direction
    dispwin : list, optional
      A user-specified window to determine the dispersion (formatted as a section)
    mode : int
      mode = 0 for longslit data where msframe=msarc
      mode = 1 for echelle data where msframe=msflat

    Returns
    -------
    dispaxis : int
      The predominant dispersion axis of the data
    """
    msgs.info("Determining the dispersion direction")
    ds1, ds2 = msframe.shape
    if dispwin is None:
        #min1, max1 = ds1/2-10, ds1/2+10
        #min2, max2 = ds2/2-10, ds2/2+10
        min1, max1 = ds1/4, 3*(ds1/4)
        min2, max2 = ds2/4, 3*(ds2/4)
    elif type(dispwin) is list:  # User has specified the location of the window (x1:x2,y1:y2)
        min1, max1 = dispwin[0]
        min2, max2 = dispwin[1]
    else: # User has specified the size of the window
        min1, max1 = ds1/2-dispwin, ds1/2+dispwin
        min2, max2 = ds2/2-dispwin, ds2/2+dispwin
    # Generate the two test statistics
    #test1 = np.median(msframe[min1:max1,:],axis=0)
    #test2 = np.median(msframe[:,min2:max2],axis=1)
    test1 = np.mean(msframe[min1:max1,:], axis=0)   # Using mean for LRIS
    test2 = np.mean(msframe[:,min2:max2], axis=1)
    # Calculate the step difference
    htst1 = test1[1:]-test1[:-1]
    htst2 = test2[1:]-test2[:-1]
    # Get the standard deviation of the step difference
    std1, std2 = np.std(htst1), np.std(htst2)
    # Return the dispersion axis
    if std1 > std2:
        if mode == 0:
            msgs.info("Dispersion axis is predominantly along a column")
            return 1
        else:
            msgs.info("Dispersion axis is predominantly along a row")
            return 0
    else:
        if mode == 0:
            msgs.info("Dispersion axis is predominantly along a row")
            return 0
        else:
            msgs.info("Dispersion axis is predominantly along a column")
            return 1


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
    mskpix, coeff = arutils.robust_polyfit(trcxrng, trcprof, 1+npix/40, function='legendre', sigma=2.0, minv=0.0, maxv=1.0)
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
        nobj = 1
    if nobj==1:
        msgs.info("Found {0:d} object".format(objl.size))
        msgs.info("Tracing {0:d} object".format(objl.size))
    else:
        msgs.info("Found {0:d} objects".format(objl.size))
        msgs.info("Tracing {0:d} objects".format(objl.size))
    # Trace objects
    cval = np.zeros(nobj)
    allsfit = np.array([])
    allxfit = np.array([])
    for o in xrange(nobj):
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
    for o in xrange(nobj): trcfunc[:,o] += cval[o]
    if nobj==1: msgs.info("Converting object trace to detector pixels")
    else: msgs.info("Converting object traces to detector pixels")
    ofst = slf._lordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1) + triml
    diff = (slf._rordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1)
            - slf._lordloc[det-1][:,order].reshape((-1,1)).repeat(nobj,axis=1))
    # Convert central trace
    traces = ofst + (diff-triml-trimr)*trcfunc
    # Convert left object trace
    for o in xrange(nobj): trccopy[:,o] = trcfunc[:,o] - cval[o] + objl[o]/(npix-1.0)
    trobjl = ofst + (diff-triml-trimr)*trccopy
    # Convert right object trace
    for o in xrange(nobj): trccopy[:,o] = trcfunc[:,o] - cval[o] + objr[o]/(npix-1.0)
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
    if doqa:
        arqa.obj_trace_qa(slf, sciframe, trobjl, trobjr, root="object_trace", normalize=False)
    # Trace dict
    tracedict = dict({})
    tracedict['nobj'] = nobj
    tracedict['traces'] = traces
    tracedict['object'] = rec_obj_img
    tracedict['background'] = rec_bg_img
    return tracedict


def trace_orders(slf, mstrace, det, pcadesc="", maskBadRows=False, singleSlit=False):
    """
    This routine will traces the locations of the slit edges

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
    singleSlit : bool, optional
      If True, only the most significant slit edge identified will be returned
      Set singleSlit=True for longslit data.

    Returns
    -------
    lcenint : ndarray
      Locations of the left slit edges (in physical pixel coordinates)
    rcenint : ndarray
      Locations of the right slit edges (in physical pixel coordinates)
    extrapord : ndarray
      A boolean mask indicating if an order was extrapolated (True = extrapolated)
    """
    msgs.info("Preparing trace frame for order edge detection")
    # Generate a binned version of the trace frame
    msgs.work("binby=1 makes this slow and ineffective -- increase this to 10, and add as a parameter of choice by the user")
    binby = 5
    binarr = arcyutils.bin_x(mstrace, binby, 0)
    binbpx = arcyutils.bin_x(slf._bpix[det-1], binby, 0)
    plxbin = arcyutils.bin_x(slf._pixlocn[det-1][:,:,0], binby, 1)
    plybin = arcyutils.bin_x(slf._pixlocn[det-1][:,:,1], binby, 1)
    msgs.info("Detecting order edges")
    #debugger.set_trace()
    if singleSlit:
        edgearr = np.zeros(binarr.shape, dtype=np.int)
        detect = True
        # Add a user-defined slit?
        # Syntax is a list of values, 2 per detector that define the slit
        # according to column values.  The 2nd value (for the right edge)
        # must be >0 to be applied.  Example for LRISr [-1, -1, 7, 295]
        # which means the code skips user-definition for the first detector
        # but adds one for the 2nd.
        if len(slf._argflag['trace']['orders']['sng_slit']) > 0:
            ledge, redge = (det-1)*2, (det-1)*2+1
            if slf._argflag['trace']['orders']['sng_slit'][redge] > 0:
                msgs.warn("Using input slit edges on detector {:d}: [{:g},{:g}]".format(
                        det,
                        slf._argflag['trace']['orders']['sng_slit'][ledge],
                        slf._argflag['trace']['orders']['sng_slit'][redge]))
                msgs.warn("Better know what you are doing!")
                edgearr[:, slf._argflag['trace']['orders']['sng_slit'][ledge]] = -1
                edgearr[:, slf._argflag['trace']['orders']['sng_slit'][redge]] = +1
                detect = False
        if detect:
            msgs.info("Detecting slit edges")
            filt = ndimage.sobel(binarr, axis=1, mode='constant')
            msgs.info("Applying bad pixel mask")
            filt *= (1.0-binbpx)  # Apply to the old detection algorithm
            amin = np.argmin(filt, axis=1)
            amax = np.argmax(filt, axis=1)
            edgearr[np.arange(edgearr.shape[0]), amin] = +1
            edgearr[np.arange(edgearr.shape[0]), amax] = -1
        # Even better would be to fit the filt/sqrt(abs(binarr)) array with a Gaussian near the maximum in each column
    else:
        ######
        # Old detection algorithm
        tedgear = arcytrace.detect_edges(binarr, slf._dispaxis)
        #arutils.ds9plot(tedgear)
        ######
        # New detection algorithm
        # First do left edges
        troll = np.roll(binarr, 1, axis=1-slf._dispaxis)
        if slf._dispaxis == 0:
            troll[:,0] = troll[:,1]
        else:
            troll[0,:] = troll[1,:]
        # First test for an edge
        diff = np.zeros_like(binarr)
        w = np.where(troll != 0.0)
        diff[w] = (troll[w]-binarr[w])/binarr[w]
        siglev = 1.4826*np.median(np.abs(diff))
        ttedges = np.zeros_like(binarr)
        wr = np.where(diff > +6.0*siglev)
        wl = np.where(diff < -6.0*siglev)
        ttedges[wr] = +1.0
        ttedges[wl] = -1.0
        # Second test for an edge
        diff = (troll-binarr)
        siglev = 1.4826*np.median(np.abs(diff))
        tedges = np.zeros_like(binarr)
        wr = np.where((diff > +6.0*siglev) & (ttedges == +1))
        wl = np.where((diff < -6.0*siglev) & (ttedges == -1))
        tedges[wr] = +1.0
        tedges[wl] = -1.0
        nedgear = arcytrace.clean_edges(diff, tedges, slf._dispaxis)
        #arutils.ds9plot(nedgear)
        if maskBadRows:
            msgs.info("Searching for bad pixel rows")
            edgsum = np.sum(nedgear, axis=0)
            sigma = 1.4826*np.median(np.abs(edgsum-np.median(edgsum)))
            w = np.where(np.abs(edgsum) >= 1.5*sigma)[0]
        #   maskcols = np.unique(np.append(w,np.append(np.append(w+2,w+1),np.append(w-2,w-1))))
            maskcols = np.unique(np.append(w,np.append(w+1, w-1)))
            msgs.info("Masking {0:d} bad pixel rows".format(maskcols.size))
            for i in xrange(maskcols.size):
                if maskcols[i] < 0 or maskcols[i] >= nedgear.shape[1]: continue
                nedgear[:, maskcols[i]] = 0
        ######
        msgs.info("Applying bad pixel mask")
        tedgear *= (1.0-binbpx)  # Apply to the old detection algorithm
        nedgear *= (1.0-binbpx)  # Apply to the new detection algorithm
        eroll = np.roll(binbpx, 1, axis=1)
        eroll[:,0] = eroll[:,1]
        nedgear *= (1.0-eroll)  # Apply to the new detection algorithm (with shift)
        # Now roll back
        nedgear = np.roll(nedgear, -1, axis=1)
        edgearr = np.zeros_like(nedgear)
        edgearr[np.where((nedgear == +1) | (tedgear == +1))] = +1
        edgearr[np.where((nedgear == -1) | (tedgear == -1))] = -1
    #arutils.ds9plot(edgearr)
    # Assign a number to each of the edges
    msgs.info("Matching order edges")
    lcnt, rcnt = arcytrace.match_edges(edgearr, 0)
    if lcnt == 1: letxt = "edge"
    else: letxt = "edges"
    if rcnt == 1: retxt = "edge"
    else: retxt = "edges"
    msgs.info("{0:d} left {1:s} and {2:d} right {3:s} were found in the trace".format(lcnt, letxt, rcnt, retxt))
    if (lcnt == 0) & (rcnt == 0):
        if np.median(binarr) > 500:
            msgs.warn("Found flux but no edges.  Assuming they go to the edge of the detector.")
            edgearr[:,-1] = 1000
            rcnt = 1
            edgearr[:,0] = -1000
            lcnt = 1
        else:
            msgs.error("Unable to trace any edges"+msgs.newline()+"try a different method to trace the order edges")
    elif (rcnt == 0) & (lcnt == 1):
        msgs.warn("Unable to find a right edge. Adding one in.")
        edgearr[:,-1] = 1000
        rcnt = 1
    elif (lcnt == 0) & (rcnt == 1):
        msgs.warn("Unable to find a left edge. Adding one in.")
        edgearr[:,0] = -1000
        lcnt = 1
    msgs.info("Assigning orders")
    #arutils.ds9plot(edgearr)
    iterate = True
    while iterate:
        lmin, lmax, rmin, rmax = arcytrace.assign_orders(edgearr, slf._dispaxis, lcnt, rcnt)
        msgs.info("Ignoring orders that span < {0:3.2f}x{1:d} pixels on the detector".format(slf._argflag['trace']['orders']['fracignore'],int(edgearr.shape[slf._dispaxis]*binby)))
        fracpix = int(slf._argflag['trace']['orders']['fracignore']*edgearr.shape[slf._dispaxis])
        lnc, lxc, rnc, rxc, ldarr, rdarr = arcytrace.ignore_orders(edgearr, slf._dispaxis, fracpix, lmin, lmax, rmin, rmax)
        lmin += lnc
        rmin += rnc
        lmax -= lxc
        rmax -= rxc
        iterate = False
        if singleSlit: # Another check on slits for singleSlit
            if lmax < lmin:
                msgs.warn("Unable to find a left edge2. Adding one in.")
                iterate = True
                edgearr[:,0] = -1000
                lcnt = 1
            if rmax < rmin:
                msgs.warn("Unable to find a right edge2. Adding one in.")
                iterate = True
                edgearr[:,-1] = 1000
                rcnt = 1
    # Left order traces
    msgs.info("Fitting left order traces")
    lcoeff = np.zeros((1+slf._argflag['trace']['orders']['polyorder'],lmax-lmin+1))
#	lfail = np.array([])
    minvf, maxvf = slf._pixlocn[det-1][0,0,0], slf._pixlocn[det-1][-1,0,0]
    for i in xrange(lmin, lmax+1):
        w = np.where(edgearr == -i)
        if np.size(w[0]) <= slf._argflag['trace']['orders']['polyorder']+2:
#			lfail = np.append(lfail,i-lmin)
            continue
        tlfitx = plxbin[w]
        tlfity = plybin[w]
        lcoeff[:, i-lmin] = arutils.func_fit(tlfitx, tlfity, slf._argflag['trace']['orders']['function'],
                                             slf._argflag['trace']['orders']['polyorder'], minv=minvf, maxv=maxvf)
#		xv=np.linspace(0,edgearr.shape[slf._dispaxis-0])
#		yv=np.polyval(coeffl[i-lmin,:],xv)
#		plt.plot(w[slf._dispaxis-0],w[1-slf._dispaxis],'ro')
#		plt.plot(xv,yv,'r-')
#		plt.show()
#		plt.clf()
    msgs.info("Fitting right order traces")
    rcoeff = np.zeros((1+slf._argflag['trace']['orders']['polyorder'], rmax-rmin+1))
#	rfail = np.array([])
    for i in xrange(rmin, rmax+1):
        w = np.where(edgearr == i)
        if np.size(w[0]) <= slf._argflag['trace']['orders']['polyorder']+2:
#			rfail = np.append(rfail, i-rmin)
            continue
        tlfitx = plxbin[w]
        tlfity = plybin[w]
        rcoeff[:, i-rmin] = arutils.func_fit(tlfitx, tlfity, slf._argflag['trace']['orders']['function'],
                                             slf._argflag['trace']['orders']['polyorder'], minv=minvf, maxv=maxvf)
    # Check if no further work is needed (i.e. there only exists one order)
    if (lmax+1-lmin == 1) and (rmax+1-rmin == 1):
        # Just a single order has been identified (i.e. probably longslit)
        msgs.info("Only one order was identified.  Should be a longslit.")
        xint = slf._pixlocn[det-1][:,0,0]
        lcenint = np.zeros((mstrace.shape[0], 1))
        rcenint = np.zeros((mstrace.shape[0], 1))
        lcenint[:,0] = arutils.func_val(lcoeff[:,0], xint, slf._argflag['trace']['orders']['function'],
                                        minv=minvf, maxv=maxvf)
        rcenint[:,0] = arutils.func_val(rcoeff[:,0], xint, slf._argflag['trace']['orders']['function'],
                                        minv=minvf, maxv=maxvf)
        return lcenint, rcenint, np.zeros(1, dtype=np.bool)
    msgs.info("Synchronizing left and right order traces")
    # Define the array of pixel values along the dispersion direction
    xv = plxbin[:,0]
    #midval = np.mean(xv)
    num = (lmax-lmin)/2
    lval = lmin + num  # Pick an order, somewhere in between lmin and lmax
    lv = (arutils.func_val(lcoeff[:,lval-lmin], xv, slf._argflag['trace']['orders']['function'],
                           minv=minvf, maxv=maxvf)+0.5).astype(np.int)
    mnvalp = np.median(binarr[:, lv+1])  # Go one row above and one row below an order edge,
    mnvalm = np.median(binarr[:, lv-1])  # then see which mean value is greater.

    """
    lvp = (arutils.func_val(lcoeff[:,lval+1-lmin],xv,slf._argflag['trace']['orders']['function'],min=minvf,max=maxvf)+0.5).astype(np.int)
    edgbtwn = arcytrace.find_between(edgearr,lv,lvp,slf._dispaxis,1)
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
        lvp = (arutils.func_val(lcoeff[:,lval+1-lmin], xv, slf._argflag['trace']['orders']['function'],
                                minv=minvf, maxv=maxvf)+0.5).astype(np.int)
        edgbtwn = arcytrace.find_between(edgearr, lv, lvp, slf._dispaxis, 1)
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
        lvp = (arutils.func_val(lcoeff[:,lval-1-lmin], xv, slf._argflag['trace']['orders']['function'],
                                minv=minvf, maxv=maxvf)+0.5).astype(np.int)
        edgbtwn = arcytrace.find_between(edgearr, lvp, lv, slf._dispaxis, -1)
        if edgbtwn[0] == -1 and edgbtwn[1] == -1:
            rsub = edgbtwn[2]-(lval-1)  # There's an order overlap
        elif edgbtwn[1] == -1:  # No overlap
            rsub = edgbtwn[0]-(lval-1)
        else:  # There's an order overlap
            rsub = edgbtwn[1]-(lval-1)

#	rva = arutils.func_val(rcoeff[:,0],midval,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#	for i in xrange(1,rmax-rmin+1):
#		rvb = arutils.func_val(rcoeff[:,i],midval,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#		if rvb > lv and rva < lv:
#			lvp = arutils.func_val(lcoeff[:,500-lmin],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			lvm = arutils.func_val(lcoeff[:,500-lmin],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvap = arutils.func_val(rcoeff[:,i-1],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvam = arutils.func_val(rcoeff[:,i-1],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvbp = arutils.func_val(rcoeff[:,i],midval*2.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			rvbm = arutils.func_val(rcoeff[:,i],0.0,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)
#			mina = np.min([lv-rva,lvp-rvap,lvm-rvam])
#			minb = np.min([rvb-lv,rvbp-lvp,rvbm-lvm])
#			if mina <= 1.0 or minb <= 0.0:
#				msgs.error("Orders are too close or the fitting quality is too poor")
#			yva = np.arange(0.0,mina)
#			yvb = np.arange(0.0,minb)
#			xmga, ymga = np.meshgrid(xv,yva)
#			xmgb, ymgb = np.meshgrid(xv,yvb)
#			ymga += np.array([arutils.func_val(rcoeff[:,i-1],xv,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)])
#			ymgb += np.array([arutils.func_val(rcoeff[:,i],xv,slf._argflag['trace']['orders']['function'],min=xvmn,max=xvmx)])
#			xmga = xmga.flatten().astype(np.int)
#			ymga = ymga.flatten().astype(np.int)
#			xmgb = xmgb.flatten().astype(np.int)
#			ymgb = ymgb.flatten().astype(np.int)
#			if slf._dispaxis == 0:
#				meda = np.median(binarr[xmga,ymga])
#				medb = np.median(binarr[xmgb,ymgb])
#			else:
#				meda = np.median(binarr[ymga,xmga])
#				medb = np.median(binarr[ymgb,xmgb])
#			if meda > medb:
#				rsub = rmin+i-1-500
#			else:
#				rsub = rmin+i-500
#			break
#		else:
#			rva = rvb
    msgs.info("Relabelling order edges")
    if lmin < rmin-rsub:
        esub = (lmin) - (slf._argflag['trace']['orders']['pcxneg']+1)
    else:
        esub = (rmin-rsub) - (slf._argflag['trace']['orders']['pcxneg']+1)

    wl = np.where(edgearr < 0)
    wr = np.where(edgearr > 0)
    edgearr[wl] += esub
    edgearr[wr] -= (esub+rsub)

    # Insert new rows into coefficients arrays if rsub != 0 (if orders were not labelled correctly, there will be a mismatch for the lcoeff and rcoeff)
    almin, almax = -np.max(edgearr[wl]), -np.min(edgearr[wl]) # min and max switched because left edges have negative values
    armin, armax = np.min(edgearr[wr]), np.max(edgearr[wr])
    nmord = slf._argflag['trace']['orders']['polyorder']+1
    if armin != almin:
        if armin < almin:
            lcoeff = np.append(np.zeros((nmord,almin-armin)),lcoeff,axis=1)
        else:
            rcoeff = np.append(np.zeros((nmord,armin-almin)),rcoeff,axis=1)
    if armax != almax:
        if armax < almax:
            rcoeff = np.append(rcoeff,np.zeros((nmord,almax-armax)),axis=1)
        else:
            lcoeff = np.append(lcoeff,np.zeros((nmord,armax-almax)),axis=1)

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
    maxord = np.max(np.append(gord,np.append(-lunq[lgm],runq[rgm])))
    #addnbad = maxord-np.max(gord)
    lcent = arutils.func_val(lcoeff[:,-lunq[lg][::-1]-1-slf._argflag['trace']['orders']['pcxneg']], xv,
                             slf._argflag['trace']['orders']['function'], minv=minvf, maxv=maxvf)
    rcent = arutils.func_val(rcoeff[:,runq[rg]-1-slf._argflag['trace']['orders']['pcxneg']], xv,
                             slf._argflag['trace']['orders']['function'], minv=minvf, maxv=maxvf)
    slitcen = 0.5*(lcent+rcent).T
    ##############
#	zmin, zmax = arplot.zscale(binarr)
#	if slf._dispaxis == 0:
#		extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
#	else:
#		extnt = (slf._pixlocn[0,0,0], slf._pixlocn[0,-1,0], slf._pixlocn[0,0,1], slf._pixlocn[-1,0,1])
#	implot = plt.imshow(binarr, extent=extnt, origin='lower', interpolation='none', aspect='auto')
#	implot.set_cmap("gray")
#	plt.colorbar()
#	implot.set_clim(zmin,zmax)
#	# Interpolate the best solutions for all orders with a cubic spline
#	if slf._dispaxis == 0:
#		xint = slf._pixlocn[:,0,0]
#	else:
#		xint = slf._pixlocn[0,:,0]
#	for i in xrange(rcent.shape[0]):
#		pclr = '-'
#		plt.plot(np.arange(lcent.shape[1]),lcent[i,:],'g'+pclr,linewidth=3.1)
#		plt.plot(np.arange(rcent.shape[1]),rcent[i,:],'b'+pclr,linewidth=3.1)
#	plt.show()
#	null = raw_input("wait...")
    ##############
    #maskord = np.where((np.all(lcoeff[:,lg],axis=0)==False)|(np.all(rcoeff[:,rg],axis=0)==False))[0]
    maskord = np.where((np.all(lcoeff, axis=0) == False) | (np.all(rcoeff, axis=0) == False))[0]
# 	print almin, armin, almax, armax
    ordsnd = np.arange(min(almin, armin), max(almax, armax)+1)
    totord = ordsnd[-1]+slf._argflag['trace']['orders']['pcxpos']
    # Identify the orders to be extrapolated during reconstruction
    extrapord = (1.0-np.in1d(np.linspace(1.0, totord, totord), gord).astype(np.int)).astype(np.bool)
    msgs.info("Performing a PCA on the order edges")
    ofit = slf._argflag['trace']['orders']['pca']
    lnpc = len(ofit)-1
    msgs.work("May need to do a check here to make sure ofit is reasonable")
    coeffs = arutils.func_fit(xv, slitcen, slf._argflag['trace']['orders']['function'],
                              slf._argflag['trace']['orders']['polyorder'], minv=minvf, maxv=maxvf)
    for i in xrange(ordsnd.size):
        if i in maskord:
            coeffs = np.insert(coeffs, i, 0.0, axis=1)
            slitcen = np.insert(slitcen, i, 0.0, axis=1)
            lcent = np.insert(lcent, i, 0.0, axis=0)
            rcent = np.insert(rcent, i, 0.0, axis=0)
    xcen = xv[:,np.newaxis].repeat(ordsnd.size, axis=1)
    fitted, outpar = arpca.basis(xcen, slitcen, coeffs, lnpc, ofit, x0in=ordsnd, mask=maskord, skipx0=True,
                                 function=slf._argflag['trace']['orders']['function'])
    # If the PCA worked OK, do the following
    msgs.work("Should something be done here inbetween the two basis calls?")
    fitted, outpar = arpca.basis(xcen, slitcen, coeffs, lnpc, ofit, x0in=ordsnd, mask=maskord,
                                 skipx0=False, function=slf._argflag['trace']['orders']['function'])
    arpca.pc_plot(slf, outpar, ofit, pcadesc=pcadesc)
    # Extrapolate the remaining orders requested
    orders = 1.0+np.arange(totord)
    extrap_cent, outpar = arpca.extrapolate(outpar, orders, function=slf._argflag['trace']['orders']['function'])
    # Fit a function for the difference between left and right edges.
    diff_coeff, diff_fit = arutils.polyfitter2d(rcent-lcent, mask=maskord,
                                                order=slf._argflag['trace']['orders']['diffpolyorder'])
    # Now extrapolate the order difference
    ydet = np.linspace(0.0, 1.0, lcent.shape[0])
    ydetd = ydet[1]-ydet[0]
    lnum = ordsnd[0]-1.0
    ydet = np.append(-ydetd*np.arange(1.0, 1.0+lnum)[::-1], ydet)
    ydet = np.append(ydet, 1.0+ydetd*np.arange(1.0, 1.0+slf._argflag['trace']['orders']['pcxpos']))
    xde, yde = np.meshgrid(np.linspace(0.0, 1.0, lcent.shape[1]), ydet)
    extrap_diff = arutils.polyval2d(xde, yde, diff_coeff).T
    msgs.info("Refining the trace for reconstructed and predicted orders")
    #refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff, [slf._argflag['trace']['orders']['pcxneg'],slf._argflag['trace']['orders']['pcxpos']], orders, slf._dispaxis, ofit[0], slf._pixlocn, function=slf._argflag['trace']['orders']['function'])
    ######
    ## NOTE::  MIGHT NEED TO APPLY THE BAD PIXEL MASK HERE TO BINARR
    msgs.work("Should the bad pixel mask be applied to the frame here?")
    refine_cent, outpar = refine_traces(binarr, outpar, extrap_cent, extrap_diff,
                                        [gord[0]-orders[0], orders[-1]-gord[-1]], orders, slf._dispaxis,
                                        ofit[0], slf._pixlocn[det-1], msgs,
                                        function=slf._argflag['trace']['orders']['function'])
    # Generate the left and right edges
    lcen = refine_cent - 0.5*extrap_diff
    rcen = refine_cent + 0.5*extrap_diff
#	lcen = extrap_cent - 0.5*extrap_diff
#	rcen = extrap_cent + 0.5*extrap_diff
    # Interpolate the best solutions for all orders with a cubic spline
    msgs.info("Interpolating the best solutions for all orders with a cubic spline")
    # zmin, zmax = arplot.zscale(binarr)
    # plt.clf()
    # extnt = (slf._pixlocn[0,0,1], slf._pixlocn[0,-1,1], slf._pixlocn[0,0,0], slf._pixlocn[-1,0,0])
    # implot = plt.imshow(binarr, extent=extnt, origin='lower', interpolation='none', aspect='auto')
    # implot.set_cmap("gray")
    # plt.colorbar()
    # implot.set_clim(zmin, zmax)
    xint = slf._pixlocn[det-1][:,0,0]
    lcenint = np.zeros((mstrace.shape[0], lcen.shape[1]))
    rcenint = np.zeros((mstrace.shape[0], rcen.shape[1]))
    for i in xrange(lcen.shape[1]):
        # if i+1 not in gord: pclr = '--'
        # else: pclr = '-'
        # plt.plot(xv_corr, lcen_corr[:,i], 'g'+pclr, linewidth=0.1)
        # plt.plot(xv_corr, rcen_corr[:,i], 'b'+pclr, linewidth=0.1)
        lcenint[:,i] = arutils.spline_interp(xint, xv, lcen[:,i])
        rcenint[:,i] = arutils.spline_interp(xint, xv, rcen[:,i])
#			plt.plot(xint,lcenint[:,i],'k.')
#			plt.plot(xint,rcenint[:,i],'k.')
#		plt.show()
    # if tracedesc != "":
    #     plt.title(tracedesc)
    # slf._qa.savefig(dpi=1200, orientation='portrait', bbox_inches='tight')
    # plt.close()

    # Illustrate where the orders fall on the detector (physical units)
    if slf._argflag['run']['qcontrol']:
        msgs.work("Not yet setup with ginga")
        # # Set up a ds9 instance
        # d = ds9.ds9()
        # # Load the image
        # d.set_np2arr(mstrace)
        # # Zoom to fit
        # d.set('zoom to fit')
        # # Change the colormap and scaling
        # d.set('cmap gray')
        # d.set('scale log')
        # # Plot the regions
        # if tracedesc != "":
        #     d.set('regions load ' + '"' + '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'],tracedesc) + '"')
        # else:
        #     d.set('regions load ' + '"' + '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir']) + '"')
        # # Save the image
        # # Check if the user wants to peruse the output as it becomes available
        # if slf._argflag['run']['stopcheck']:
        #     null=raw_input(msgs.input()+"Press enter to continue...")
        # else:
        #     msgs.info("DS9 window was updated")
    return lcenint, rcenint, extrapord


def refine_traces(binarr, outpar, extrap_cent, extrap_diff, extord, orders, dispaxis,
                  fitord, locations, function='polynomial'):
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
        minarrL = arcytrace.minbetween(binarr, loord, hiord, dispaxis) # Minimum counts between loord and hiord
        minarrR = arcytrace.minbetween(binarr, hiord, nxord, dispaxis)
        minarr = 0.5*(minarrL+minarrR)
        srchz = np.abs(extfit[:,-i]-extfit[:,-i-1])/3.0
        lopos = phys_to_pix(extfit[:,-i]-srchz, locations, 1) # The pixel indices for the bottom of the search window
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,-i])))
        diffarr = np.round(extrap_diff[:,-i]).astype(np.int)
        shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch, dispaxis)
        relshift = np.mean(shift+extrap_diff[:,-i]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[-i])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[-i]), relshift))
        # Renew guess for the next order
        mask[-i] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, -i, relshift, fitord, msgs, function=function)
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
        minarr = arcytrace.minbetween(binarr,loord, hiord, dispaxis)
        srchz = np.abs(extfit[:,i]-extfit[:,i-1])/3.0
        lopos = phys_to_pix(extfit[:,i-1]-srchz, locations, 1)
        numsrch = np.int(np.max(np.round(2.0*srchz-extrap_diff[:,i-1])))
        diffarr = np.round(extrap_diff[:,i-1]).astype(np.int)
        shift = arcytrace.find_shift(binarr, minarr, lopos, diffarr, numsrch, dispaxis)
        relshift = np.mean(shift+extrap_diff[:,i-1]/2-srchz)
        if shift == -1:
            msgs.info("  Refining order {0:d}: NO relative shift applied".format(int(orders[i-1])))
            relshift = 0.0
        else:
            msgs.info("  Refining order {0:d}: relative shift = {1:+f}".format(int(orders[i-1]),relshift))
        # Renew guess for the next order
        mask[i-1] = 1.0
        extfit, outpar, fail = arpca.refine_iter(outpar, orders, mask, i-1, relshift, fitord, msgs, function=function)
        if fail:
            msgs.warn("Order refinement has large residuals -- check order traces")
            return extrap_cent, outparcopy
        i -= 1
    return extfit, outpar


def model_tilt(slf, det, msarc, censpec=None, maskval=-999999.9,
               trthrsh=1000.0):
    """
    This function performs a PCA analysis on the arc tilts for a single spectrum (or order)
    """
    msgs.work("Detecting lines")
    tampl, tcent, twid, w, satsnd, _ = ararc.detect_lines(slf, det, msarc, censpec=censpec)
    satval = slf._spect['det'][det-1]['saturation']*slf._spect['det'][det-1]['nonlinear']
    if slf._argflag['trace']['orders']['tilts'].lower() == "zero":
        # Assuming there is no spectral tilt
        tilts = np.outer(np.linspace(0.0, 1.0, msarc.shape[0]), np.ones(msarc.shape[1]))
        return tilts, satsnd
    # Order of the polynomials to be used when fitting the tilts.
    fitxy = [slf._argflag['trace']['orders']['tiltorder'], 1]
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
    for s in xrange(nuse):
        w = np.where((np.abs(arcdet-detuse[s]) <= ncont) & (np.abs(arcdet-detuse[s]) >= 1.0))[0]
        for u in xrange(w.size):
            if ampl[w[u]] > ampl[olduse][s]:
                aduse[idxuse[s]] = False
                break
    # Restricted to ID lines? [avoid ghosts]
    if slf._argflag['trace']['orders']['use_ids_only']:
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
    # for s in xrange(Nseg):
    #     w = np.where((arcdet > s*segsz) & (arcdet <= (s+1)*segsz))[0]
    #     segampl = tampl[w]
    #     asrt = np.argsort(segampl)[::-1]
    #     for u in xrange(Nuse):
    #         aduse[w[asrt[u]]] = True
    # # Now include some additional bright lines
    # asrt = np.argsort(tampl)[::-1]
    # s, u = 0, 0
    # while u < Nadd:
    #     if not aduse[asrt[s]]:
    #         aduse[asrt[s]] = True
    #         u += 1
    #     s += 1

    if slf._argflag['trace']['orders']['tilts'].lower() == "pca":
        nitertilts = 2
    else:
        nitertilts = 1
    ordcen = slf._pixcen[det-1].copy()
    guesstilts = None
    outpar = {}
    msgs.info("Modelling arc line tilts with {0:d} arc lines".format(np.sum(aduse)))
    # Go along each order and trace the tilts
    for tt in xrange(nitertilts):
        # Start by masking every row, then later unmask the rows with usable arc lines
        maskrows = np.ones(msarc.shape[0], dtype=np.int)
        if nitertilts > 1:
            msgs.info("Iterating on spectral tilts -- Iteration {0:d}/{1:d}".format(tt+1, nitertilts))
        tcoeff = np.ones((slf._argflag['trace']['orders']['tiltorder']+1, msarc.shape[0]))
        weights = np.zeros(msarc.shape[0])
        msgs.work("This next step could be multiprocessed to speed up the reduction")
        nspecfit = 7
        nsmth = 3  # Number of pixels +/- in the spatial direction to include in the fit (0=no smoothing, 1=3 pixels, 2=5 pixels...)
        xtilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        ytilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        ztilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        stilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        mtilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        wtilt = np.ones((msarc.shape[1], arcdet.size))*maskval
        badlines = 0
        for j in xrange(arcdet.size):
            # For each detection in this order
            msgs.info("Tracing tilt of arc line {0:d}/{1:d}".format(j+1, arcdet.size))
            # Check if this is a saturated line
            ysat = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1, ordcen[arcdet[j], 0]-nsmth:ordcen[arcdet[j], 0]+nsmth+1]
            if np.where(ysat > satval)[0].size != 0:
                aduse[j] = False
                badlines += 1
                continue
            # Get the size of the slit
            sz = int(np.floor(np.abs(slf._rordloc[det-1][arcdet[j], 0]-slf._lordloc[det-1][arcdet[j], 0])/2.0)) - 2
            xtfit = np.zeros(2*sz+1)
            ytfit = np.ones(2*sz+1)*maskval  # Fitted centroid
            etfit = np.zeros(2*sz+1)  # Fitted centroid error
            mtfit = np.ones(2*sz+1, dtype=np.int)   # Mask of bad fits
            apfit = np.zeros(2*sz+1)  # Fitted Amplitude
            shfit = np.zeros(2*sz+1)  # Shift applied
            xfit = np.arange(-nspecfit, nspecfit+1, 1.0)
            tstcc = True  # A boolean to tell the loop once a good set of pixels has been found to cross-correlate with
            # Fit up
            pcen = arcdet[j]
            if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                # Too close to the end of the spectrum
                aduse[j] = False
                badlines += 1
                continue
            offchip = False
            for k in xrange(0, sz+1-nsmth):
                if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                    offchip = True
                    break
                if ordcen[pcen, 0]+k >= msarc.shape[1]:
                    offchip = True
                    break
                # Get a copy of the array that will be used to cross-correlate
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]+k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]+k-nsmth:ordcen[arcdet[j],0]+k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]+k]
                yfit = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], 0]+k-nsmth:ordcen[arcdet[j], 0]+k+nsmth+1]
                if len(yfit.shape) == 2:
                    yfit = np.median(yfit, axis=1)
                if np.size(yfit) == 0:
                    offchip = True
                    break
                # wgd = np.where((yfit<satval)&(yfit!=maskval))
                wgd = np.where(yfit == maskval)
                if wgd[0].size != 0:
                    continue
                cc = np.correlate(ccyfit, yfit, mode='same')
                params, fail = arutils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
                if guesstilts is not None:
                    shfit[k+sz] = guesstilts[arcdet[j], ordcen[arcdet[j], 0]+k]*float(msarc.shape[0]-1.0)
                xtfit[k+sz] = ordcen[arcdet[j], 0] + k
                ytfit[k+sz] = centv - shfit[k+sz]
                etfit[k+sz] = 0.02
                apfit[k+sz] = params[0]
                if fail:
                    mtfit[k+sz] = 1
                else:
                    pcen = int(0.5 + centv)
                    mtfit[k+sz] = 0
            if offchip:
                # Don't use lines that go off the chip (could lead to a bad trace)
                aduse[j] = False
                badlines += 1
                continue
            for k in range(sz+1-nsmth, sz+1):
                xtfit[k+sz] = ordcen[arcdet[j], 0]+k
            # Fit down
            pcen = arcdet[j]
            for k in xrange(1,sz+1-nsmth):
                if (pcen < nspecfit) or (pcen > msarc.shape[0]-(nspecfit+1)):
                    offchip = True
                    break
                if ordcen[pcen, 0]-k < 0:
                    offchip = True
                    break
                # Get a copy of the array that will be used to cross-correlate
                # (testcc is probably already False from the Fit Up part of the loop, but it's best to be sure)
                if tstcc:
                    # ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]-k]
                    ccyfit = msarc[arcdet[j]-nspecfit:arcdet[j]+nspecfit+1,ordcen[arcdet[j],0]-k-nsmth:ordcen[arcdet[j],0]-k+nsmth+1]
                    if len(ccyfit.shape) == 2:
                        ccyfit = np.median(ccyfit, axis=1)
                    wgd = np.where(ccyfit == maskval)
                    if wgd[0].size != 0:
                        continue
                    ccval = arcdet[j] + np.sum(xfit*ccyfit)/np.sum(ccyfit)
                    tstcc = False  # Once we have an array, there's no need to keep looking
                # yfit = msarc[pcen-nspecfit:pcen+nspecfit+1,ordcen[arcdet[j],0]-k]
                yfit = msarc[pcen-nspecfit:pcen+nspecfit+1, ordcen[arcdet[j], 0]-k-nsmth:ordcen[arcdet[j], 0]-k+nsmth+1]
                if len(yfit.shape) == 2:
                    yfit = np.median(yfit, axis=1)
                if np.size(yfit) == 0:
                    offchip = True
                    break
                wgd = np.where(yfit == maskval)
                if wgd[0].size != 0:
                    continue
                cc = np.correlate(ccyfit, yfit, mode='same')
                params, fail = arutils.gauss_lsqfit(xfit, cc, 0.0)
                centv = ccval + pcen - arcdet[j] - params[1]
                if guesstilts is not None:
                    shfit[sz-k] = guesstilts[arcdet[j], ordcen[arcdet[j], 0]-k]*float(msarc.shape[0]-1.0)
                xtfit[sz-k] = ordcen[arcdet[j], 0] - k
                ytfit[sz-k] = centv - shfit[sz-k]
                etfit[sz-k] = 0.02
                apfit[sz-k] = params[0]
                if fail:
                    mtfit[sz-k] = 1
                else:
                    pcen = int(0.5 + centv)
                    mtfit[sz-k] = 0
            if offchip:
                # Don't use lines that go off the chip (could lead to a bad trace)
                aduse[j] = False
                badlines += 1
                continue
            for k in range(sz+1-nsmth, sz+1):
                xtfit[sz-k] = ordcen[arcdet[j], 0]-k

            wmask = np.where(mtfit == 0)
            if guesstilts is not None:
                ytfit[wmask] -= np.median(ytfit[wmask])

            # Perform a scanning polynomial fit to the tilts
            model = arcyutils.polyfit_scan_intext(xtfit, ytfit, np.ones(ytfit.size, dtype=np.float), mtfit,
                                                  2, sz/6, 3, maskval)
            # if np.sum(np.isnan(model)) > 0:
            #     debugger.set_trace()
            # if len(np.where(np.abs(model) > 1e10)[0]):
            #     debugger.set_trace()

            if maskval in model:
                # Model contains masked values
                aduse[j] = False
                badlines += 1
                continue

            ytfit[np.where(mtfit == 1)] = maskval

            # Perform a robust polynomial fit to the traces
            if slf._argflag['trace']['orders']['tilts'].lower() == "spca":
                yfit = ytfit[wmask]/(msarc.shape[0]-1.0)
            else:
                yfit = (2.0*model[sz]-ytfit[wmask])/(msarc.shape[0]-1.0)
            wmsk, mcoeff = arutils.robust_polyfit(xtfit[wmask], yfit,
                                                  slf._argflag['trace']['orders']['tiltorder'],
                                                  function=slf._argflag['trace']['orders']['function'],
                                                  sigma=2.0, minv=0.0, maxv=msarc.shape[1]-1.0)
            # Update the mask
            wmask = wmask[0][np.where(wmsk == 0)]

            # Save the tilt angle, and unmask the row
            factr = (msarc.shape[0]-1.0)*arutils.func_val(mcoeff, ordcen[arcdet[j],0],
                                                          slf._argflag['trace']['orders']['function'],
                                                          minv=0.0, maxv=msarc.shape[1]-1.0)
            if guesstilts is not None:
                factr += guesstilts[arcdet[j], ordcen[arcdet[j], 0]]*float(msarc.shape[0]-1.0)
            idx = int(factr+0.5)
            if (idx > 0) and (idx < msarc.shape[0]):
                maskrows[idx] = 0
                tcoeff[:, idx] = mcoeff.copy()
                weights[idx] = (np.abs(np.median(apfit[wmask])))**0.25
            # Restrict to good IDs?
            if slf._argflag['trace']['orders']['use_ids_only']:
                if not aduse[j]:
                    maskrows[idx] = 1

            # if True:
            #     #debugger.set_trace()
            #     plt.clf()
            #     plt.plot(xtfit[wmask], ytfit[wmask], 'bx')
            #     plt.plot(xtfit, model, 'r-')
            #     plt.show()
            #     continue

            xint = int(xtfit[0])
            xtilt[xint:xint+2*sz+1, j] = xtfit/(msarc.shape[1]-1.0)
            ytilt[xint:xint+2*sz+1, j] = arcdet[j]/(msarc.shape[0]-1.0)
            ztilt[xint:xint+2*sz+1, j] = ytfit/(msarc.shape[0]-1.0)
            stilt[xint:xint+2*sz+1, j] = shfit
            if slf._argflag['trace']['orders']['tilts'].lower() == "spline":
                mtilt[xint:xint+2*sz+1, j] = model/(msarc.shape[0]-1.0)
            elif slf._argflag['trace']['orders']['tilts'].lower() == "interp":
                mtilt[xint:xint+2*sz+1, j] = (2.0*model[sz]-model)/(msarc.shape[0]-1.0)
            else:
                mtilt[xint:xint+2*sz+1, j] = (2.0*model[sz]-model)/(msarc.shape[0]-1.0)
            wbad = np.where(ytfit == maskval)[0]
            ztilt[xint+wbad, j] = maskval
            if wmask.size != 0:
                sigg = max(1.4826 * np.median(np.abs(ytfit-model)[wmask])/np.sqrt(2.0), 1.0)
                wtilt[xint:xint+2*sz+1, j] = 1.0 / sigg
            # Extrapolate off the slit to the edges of the chip
            nfit = 6  # Number of pixels to fit a linear function to at the end of each trace
            xlof, xhif = np.arange(xint, xint+nfit), np.arange(xint+2*sz+1-nfit, xint+2*sz+1)
            xlo, xhi = np.arange(xint), np.arange(xint+2*sz+1, msarc.shape[1])
            glon = np.mean(xlof*mtilt[xint:xint+nfit, j]) - np.mean(xlof)*np.mean(mtilt[xint:xint+nfit, j])
            glod = np.mean(xlof**2) - np.mean(xlof)**2
            clo = np.mean(mtilt[xint:xint+nfit, j]) - (glon/glod)*np.mean(xlof)
            yhi = mtilt[xint+2*sz+1-nfit:xint+2*sz+1, j]
            ghin = np.mean(xhif*yhi) - np.mean(xhif)*np.mean(yhi)
            ghid = np.mean(xhif**2) - np.mean(xhif)**2
            chi = np.mean(yhi) - (ghin/ghid)*np.mean(xhif)
            mtilt[0:xint, j] = (glon/glod)*xlo + clo
            mtilt[xint+2*sz+1:, j] = (ghin/ghid)*xhi + chi
        if badlines != 0:
            msgs.warn("There were {0:d} additional arc lines that should have been traced".format(badlines) +
                      msgs.newline() + "(perhaps lines were saturated?). Check the spectral tilt solution")

        # Masking
        weights /= np.max(weights)
        maskrw = np.where(maskrows == 1)[0]
        maskrw.sort()
        extrap_row = maskrows.copy()
        xv = np.arange(msarc.shape[1])
        # Tilt values
        tiltval = arutils.func_val(tcoeff, xv, slf._argflag['trace']['orders']['function'],
                                   minv=0.0, maxv=msarc.shape[1]-1.0).T
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        ofit = slf._argflag['trace']['orders']['pcatilt']
        lnpc = len(ofit)-1
        # Only do a PCA if there are enough good orders
        if np.sum(1.0-extrap_row) > ofit[0]+1:
            # Perform a PCA on the tilts
            msgs.info("Performing a PCA on the tilts")
            ordsnd = np.linspace(0.0, 1.0, msarc.shape[0])
            xcen = xv[:, np.newaxis].repeat(msarc.shape[0], axis=1)
            fitted, outpar = arpca.basis(xcen, tiltval, tcoeff, lnpc, ofit, weights=None,
                                         x0in=ordsnd, mask=maskrw, skipx0=False,
                                         function=slf._argflag['trace']['orders']['function'])
            arpca.pc_plot(slf, outpar, ofit, pcadesc="Spectral Tilts PCA", addOne=False)
            # Extrapolate the remaining orders requested
            orders = np.linspace(0.0, 1.0, msarc.shape[0])
            extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=slf._argflag['trace']['orders']['function'])
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

        # Now apply the starting guess for the tilts
        if guesstilts is not None:
            polytilts += guesstilts

        if slf._argflag['trace']['orders']['tilts'].lower() == "interp":
            msgs.info("Interpolating and Extrapolating the tilts")
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # yspl = np.append(0.0, np.append(arcdet[np.where(aduse)]/(msarc.shape[0]-1.0), 1.0))
            # yspl = np.append(0.0, np.append(polytilts[arcdet[np.where(aduse)], msarc.shape[1]/2], 1.0))
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)]]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            zspl = np.zeros((msarc.shape[1], np.sum(aduse)+2))
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
        elif slf._argflag['trace']['orders']['tilts'].lower() == "spline":
            msgs.info("Performing a spline fit to the tilts")
            wgd = np.where((ytilt != maskval) & (ztilt != maskval))
            txsbs = xtilt[wgd]
            tysbs = ytilt[wgd]
            tzsbs = ztilt[wgd]
            twsbs = wtilt[wgd]
            # Append the end points
            wlo = np.where((ytilt == np.min(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
            whi = np.where((ytilt == np.max(tysbs)) & (ytilt != maskval) & (ztilt != maskval))
            xlo = (xtilt[wlo]*(msarc.shape[1]-1.0)).astype(np.int)
            xhi = (xtilt[whi]*(msarc.shape[1]-1.0)).astype(np.int)
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
                plt.plot(xsbs, (ysbs-tiltqa)/ysbs, 'bx')
                plt.plot(xsbs, 1.0/(wsbs*ysbs), 'r-')
                plt.ylim(-5e-2, 5e-2)
                # plt.colorbar()
                plt.show()
                debugger.set_trace()
        elif slf._argflag['trace']['orders']['tilts'].lower() == "spca":
            # Slit position
            xspl = np.linspace(0.0, 1.0, msarc.shape[1])
            # Trace positions down center of the order
            ycen = np.diag(polytilts[arcdet[np.where(aduse)], ordcen[arcdet[np.where(aduse)]]])
            yspl = np.append(0.0, np.append(ycen, 1.0))
            # Trace positions as measured+modeled
            zspl = np.zeros((msarc.shape[1], np.sum(aduse)+2))
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
                pmin = int(max(0, np.min(slf._lordloc[det-1])))
                pmax = int(min(msarc.shape[1], np.max(slf._rordloc[det-1])))
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
                tiltqa = tiltspl(xsbs.flatten(), zsbs.flatten(), grid=False).reshape( xsbs.shape)
                plt.clf()
                #plt.imshow((zsbs-tiltqa)/zsbs, origin='lower')
                plt.imshow((ysbs-tiltqa)/ysbs, origin='lower')
                plt.colorbar()
                plt.show()
                debugger.set_trace()
        elif slf._argflag['trace']['orders']['tilts'].lower() == "pca":
            tilts = polytilts.copy()
        if tt == 0:
            ztilto = ztilt.copy()
        elif tt < nitertilts-1:
            # Provide a starting guess of the tilts for the next iteration
            guesstilts = tilts.copy()

    # Now do the QA
    msgs.info("Preparing arc tilt QA data")
    tiltsplot = tilts[arcdet, :].T
    tiltsplot *= (msarc.shape[0]-1.0)
    # Shift the plotted tilts about the centre of the slit
    adj = np.diag(tilts[arcdet, ordcen[arcdet]])
    zmsk = np.where(ztilto == maskval)
    ztilto = 2.0*np.outer(np.ones(ztilto.shape[0]), adj) - ztilto
    ztilto[zmsk] = maskval
    ztilto[np.where(ztilto != maskval)] *= (msarc.shape[0]-1.0)
    for i in xrange(arcdet.size):
        w = np.where(ztilto[:, i] != maskval)
        if w[0].size != 0:
            twa = (xtilt[w[0], i]*(msarc.shape[1]-1.0)+0.5).astype(np.int)
            # fitcns = np.polyfit(w[0], ztilt[w[0], i] - tiltsplot[twa, i], 0)[0]
            fitcns = np.median(ztilto[w[0], i] - tiltsplot[twa, i])
            # if abs(fitcns) > 1.0:
            #     msgs.warn("The tilt of Arc Line {0:d} might be poorly traced".format(i+1))
            tiltsplot[:, i] += fitcns

    xdat = xtilt.copy()
    xdat[np.where(xdat != maskval)] *= (msarc.shape[1]-1.0)

    msgs.info("Plotting arc tilt QA")
    arplot.plot_orderfits(slf, tiltsplot, ztilto, xdata=xdat, xmodl=np.arange(msarc.shape[1]),
                          textplt="Arc line", maxp=9, desc="Arc line spectral tilts", maskval=maskval)
    return tilts, satsnd, outpar


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


def trace_tilt(slf, msarc, prefix="", tltprefix="", trcprefix=""):
    """
    The value of "tilts" returned by this function is of the form:
    tilts = tan(tilt angle), where "tilt angle" is the angle between
    (1) the line representing constant wavelength and
    (2) the column of pixels that is most closely parallel with the spatial direction of the slit.

    The angle is determined relative to the axis defined by ...

    In other words, tilts = y/x according to the docs/get_locations_orderlength.JPG file.

    """

    msgs.work("Haven't used physical pixel locations in this routine")
    if slf._argflag['trace']['orders']['tilts'] == 'zero':
        # No calculation is required, simply return the appropriate array of tilts.
        tilts = np.zeros_like(slf._lordloc[det-1])
    # Calculate the perpendicular tilts for each order (this can be used as a first guess if the user asks for the tilts to be traced)
    msgs.info("Calculating perpendicular tilts for each order")
    ocen = 0.5*(slf._lordloc[det-1]+slf._rordloc[det-1])
    dervt = ocen[1:,:]-ocen[:-1,:]
    derv = np.append(dervt[0,:].reshape((1,dervt.shape[1])), 0.5*(dervt[1:,:]+dervt[:-1,:]),axis=0)
    derv = np.append(derv, dervt[-1,:].reshape((1,dervt.shape[1])),axis=0)
#	tilts = np.arctan(-1.0/derv)*180.0/np.pi
    derv = -derv
    if slf._argflag['trace']['orders']['tilts'] == 'perp':
        tilts = derv
#	plt.plot(np.arange(msarc.shape[slf._dispaxis]),slf._lordloc[:,10],'g-')
#	plt.plot(np.arange(msarc.shape[slf._dispaxis]),slf._rordloc[:,10],'b-')
#	showtilts=np.array([10,50,800,1000,2000,3000,3300,3600])
#	for j in xrange(showtilts.size):
#		xplt = np.arange(-5,5)+showtilts[j]
#		yplt = ocen[showtilts[j],10] + np.tan(tilts[showtilts[j],10]*np.pi/180.0)*(xplt-showtilts[j])
#		yplt = np.arange(-5,5)+ocen[showtilts[j],10]
#		xplt = showtilts[j] + (yplt-ocen[showtilts[j],10])/np.tan(tilts[showtilts[j],10]*np.pi/180.0)
#		plt.plot(xplt,yplt,'r-')
#	plt.show()
    # Extract a rough spectrum of the arc in each order
    msgs.info("Extracting an approximate arc spectrum at the centre of each order")
    tordcen = None
    maskorder = np.zeros(ocen.shape[1],dtype=np.int)
    for i in xrange(ocen.shape[1]):
        if slf._dispaxis == 0:
            wl = np.size(np.where(ocen[:,i] <= slf._pixlocn[det-1][0,0,1])[0])
            wh = np.size(np.where(ocen[:,i] >= slf._pixlocn[det-1][0,-1,1])[0])
        else:
            wl = np.size(np.where(ocen[:,i] <= slf._pixlocn[det-1][0,0,1])[0])
            wh = np.size(np.where(ocen[:,i] >= slf._pixlocn[det-1][-1,0,1])[0])
        if wl==0 and wh==0: # The center of the order is always on the chip
            if tordcen is None:
                tordcen = np.zeros((ocen.shape[0],1),dtype=np.int)
                tordcen[:,0] = ocen[:,i]
            else:
                tordcen = np.append(tordcen,ocen[:,i].reshape((ocen.shape[0],1)),axis=1)
        else: # An order isn't on the chip
            if tordcen is None:
                tordcen = np.zeros((ocen.shape[0],1),dtype=np.int)
            else:
                tordcen = np.append(tordcen,ocen[:,i].reshape((ocen.shape[0],1)),axis=1)
            maskorder[i] = 1
    w = np.where(maskorder==0)[0]
    if tordcen is None:
        msgs.warn("Could not determine which full orders are on the detector")
        msgs.info("Assuming all orders are fully on the detector")
        ordcen = phys_to_pix(ocen, slf._pixlocn[det-1], 1)
    else:
        ordcen = phys_to_pix(tordcen[:,w], slf._pixlocn[det-1], 1)

    pixcen = np.arange(msarc.shape[0])
    temparr = pixcen.reshape(msarc.shape[0],1).repeat(ordcen.shape[1],axis=1)
    # Average over three pixels to remove some random fluctuations, and increase S/N
    op1 = ordcen+1
    op2 = ordcen+2
    om1 = ordcen-1
    om2 = ordcen-2
    w = np.where(om1<0)
    om1[w] += 1
    w = np.where(om2==-1)
    om2[w] += 1
    w = np.where(om2==-2)
    om2[w] += 2
    if slf._dispaxis == 0:
        w = np.where(op1>=msarc.shape[1])
        op1[w] -= 1
        w = np.where(op2==msarc.shape[1])
        op2[w] -= 1
        w = np.where(op2==msarc.shape[1]+1)
        op2[w] -= 2
        arccen = (msarc[temparr,ordcen]+msarc[temparr,op1]+msarc[temparr,op2]+msarc[temparr,om1]+msarc[temparr,om2])/5.0
#		arccel = (msarc[temparr,ordcen-2])#+msarc[temparr,ordcen+1]+msarc[temparr,ordcen-1])/3.0
#		arccer = (msarc[temparr,ordcen+2])#+msarc[temparr,ordcen+1]+msarc[temparr,ordcen-1])/3.0
    else:
        w = np.where(op1>=msarc.shape[0])
        op1[w] -= 1
        w = np.where(op2==msarc.shape[0])
        op2[w] -= 1
        w = np.where(op2==msarc.shape[0]+1)
        op2[w] -= 2
        arccen = (msarc[ordcen,temparr]+msarc[op1,temparr]+msarc[op2,temparr]+msarc[om1,temparr]+msarc[om2,temparr])/5.0
    del temparr
#		arccel = (msarc[ordcen-2,temparr])#+msarc[ordcen+1,temparr]+msarc[ordcen-1,temparr])/3.0
#		arccer = (msarc[ordcen+2,temparr])#+msarc[ordcen+1,temparr]+msarc[ordcen-1,temparr])/3.0
#	for i in xrange(arccen.shape[1]):
#		plt.clf()
#		plt.plot(pixcen,arccen[:,i],'g-')
#		plt.plot(pixcen,arccer[:,i],'r-')
#		plt.plot(pixcen,arccel[:,i],'b-')
#		plt.show()
    msgs.info("Generating a mask of arc line saturation streaks")
    satmask = arcyarc.saturation_mask(msarc, slf._spect['det']['saturation']*slf._spect['det']['nonlinear'])
    ordwid = 0.5*np.abs(slf._lordloc-slf._rordloc)
    satsnd = arcyarc.order_saturation(satmask,ordcen,(ordwid+0.5).astype(np.int),slf._dispaxis)
#	arutils.ds9plot((1.0-satmask)*msarc)
#	plt.plot(pixcen,arccen[:,10],'k-',drawstyle='steps')
#	plt.show()
    # Detect the location of the arc lines
    msgs.info("Detecting the strongest, nonsaturated arc lines in each order")
    #####
    # Old algorithm for arc line detection
#	arcdet = arcyarc.detections_allorders(arccen, satsnd)
    #####
    # New algorithm for arc line detection
    pixels=[]
    totnum = 0
    siglev = 2.0*slf._argflag['arc']['calibrate']['detection']
    bpfit = 5 # order of the polynomial used to fit the background 'continuum'
    fitp = slf._argflag['arc']['calibrate']['nfitpix']
    for o in xrange(arccen.shape[1]):
        pixels.append([])
        detns = arccen[:,o]
        xrng = np.arange(float(detns.size))
        mask = np.zeros(detns.size,dtype=np.int)
        mskcnt=0
        while True:
            w = np.where(mask==0)
            xfit = xrng[w]
            yfit = detns[w]
            ct = np.polyfit(xfit,yfit,bpfit)
            yrng = np.polyval(ct,xrng)
            sigmed = 1.4826*np.median(np.abs(detns[w]-yrng[w]))
            w = np.where(detns>yrng+1.5*sigmed)
            mask[w] = 1
            if mskcnt == np.sum(mask): break # No new values have been included in the mask
            mskcnt = np.sum(mask)
# 		plt.plot(xrng,detns,'k-',drawstyle='steps')
# 		plt.plot(xrng,yrng,'r-')
# 		plt.show()
# 		plt.clf()
        w = np.where(mask==0)
        xfit = xrng[w]
        yprep = detns - yrng
        sfit = 1.4826*np.abs(detns[w]-yrng[w])
        ct = np.polyfit(xfit,sfit,bpfit)
        yerr = np.polyval(ct,xrng)
        myerr = np.median(np.sort(yerr)[:yerr.size/2])
        yerr[np.where(yerr < myerr)] = myerr
        # Find all significant detections
        tpixt, num = arcyarc.detections_sigma(yprep,yerr,np.zeros(satsnd.shape[0],dtype=np.int),siglev/2.0,siglev) # The last argument is the overall minimum significance level of an arc line detection and the second last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
        pixt = arcyarc.remove_similar(tpixt, num)
        pixt = pixt[np.where(pixt!=-1)].astype(np.int)
        tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng,yprep,pixt,fitp)
        w = np.where((np.isnan(twid)==False) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent>0.0) & (tcent<xrng[-1]))
        pixels[o] = (tcent[w]+0.5).astype(np.int)
        if np.size(w[0])>totnum:
            totnum = np.size(w[0])
    # Convert this into an arcdet array
    arcdet = -1.0*np.ones((totnum,arccen.shape[1]))
    for o in xrange(arccen.shape[1]):
        arcdet[:pixels[o].size,o] = pixels[o]
#	y = arccen[:,10]
#	pixt = arcdet[np.where(arcdet[:,10]!=-1)[0],10].astype(np.float)
#	plt.plot(pixcen, y, 'k-', drawstyle='steps')
#	ym=(np.max(y)-np.min(y))/100.0
#	for i in xrange(len(pixt)):
#		yp = y[np.argmin(np.abs(pixt[i]-pixcen))]
#		plt.plot([pixt[i]-0.5,pixt[i]-0.5],[yp+ym/2.0,yp + ym],'b-')
#	plt.show()
    if slf._argflag['trace']['orders']['tilts'] in ['perp','zero']:
        centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        nm=0
        for i in xrange(maskorder.size):
            if maskorder[i] == 1: continue
            w = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
            if np.size(w) != 0: centval[:np.size(w),i] = w
            nm += 1
    elif slf._argflag['trace']['orders']['tilts'] == 'fit1D':
        # Go along each order and fit the tilts in 1D
        tiltang=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        msgs.work("This next step could be multiprocessed to speed up the reduction")
        nm = 0
        for i in xrange(maskorder.size): # For each order
            msgs.info("Tracing tilts -- order {0:d}/{1:d}".format(i+1,maskorder.size))
            if maskorder[i] == 1: continue
            pixt = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
            for j in xrange(pixt.size): # For each detection in this order
                sz = int(np.floor(np.abs(slf._rordloc[pixt[j],i]-slf._lordloc[pixt[j],i])/2.0))-1
                if slf._dispaxis == 0:
                    xtfit = np.arange(-sz,sz+1,1.0) # pixel along the arc line
                    ytfit = np.zeros(2*sz+1) # Fitted centroid
                    etfit = np.zeros(2*sz+1) # Fitted centroid error
                    mtfit = np.zeros(2*sz+1) # Mask of bad fits
                    # Fit up
                    pcen = pixt[j]
                    if (pcen < 3) or (pcen > msarc.shape[0]-4): continue
                    offchip = False
                    for k in xrange(0,sz+1):
                        if (pcen < 3) or (pcen > msarc.shape[0]-4):
                            offchip == True
                            break
                        if ordcen[pcen,nm]+k >= msarc.shape[1]:
                            mtfit[k+sz] = 1.0
                            offchip = True
                            break
                        xfit = np.arange(pcen-3, pcen+3+1, 1.0)  # 2x4 + 1 = 9 pixels total along the spectral dimension
                        yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,nm]+k]
                        if np.size(yfit) == 0:
                            mtfit[k+sz] = 1.0
                            offchip = True
                            break
                        #params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
                        params, fail = arutils.gauss_fit(xfit,yfit,pcen)
                        ytfit[k+sz] = params[1]
                        etfit[k+sz] = 0.02
                        if fail: mtfit[k+sz] = 1.0
                        else: pcen = int(0.5+params[1])
                    if offchip: continue
                    # Fit down
                    pcen = int(0.5+ytfit[sz]) # Start with the best-fitting centroid at arccen
                    for k in xrange(1,sz+1):
                        if (pcen < 3) or (pcen > msarc.shape[0]-4):
                            offchip == True
                            break
                        if ordcen[pcen,nm]-k < 0:
                            mtfit[sz-k] = 1.0
                            offchip = True
                            break
                        xfit = np.arange(pcen-3, pcen+3+1, 1.0)
                        yfit = msarc[pcen-3:pcen+3+1,ordcen[pcen,nm]-k]
                        if np.size(yfit) == 0:
                            mtfit[sz-k] = 1.0
                            offchip = True
                            break
                        #params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
                        params, fail = arutils.gauss_fit(xfit,yfit,pcen)
                        ytfit[sz-k] = params[1]
                        etfit[sz-k] = 0.02
                        if fail: mtfit[sz-k] = 1.0
                        else: pcen = int(0.5+params[1])
                    if offchip: continue
                else:
                    xtfit = np.arange(-sz,sz+1,1.0) # pixel along the arc line
                    ytfit = np.zeros(2*sz+1) # Fitted centroid
                    etfit = np.zeros(2*sz+1) # Fitted centroid error
                    mtfit = np.zeros(2*sz+1) # Mask of bad fits
                    # Fit up
                    pcen = pixt[j]
                    if (pcen < 3) or (pcen > msarc.shape[1]-4): continue
                    offchip = False
                    for k in xrange(0,sz+1):
                        if (pcen < 3) or (pcen > msarc.shape[0]-4):
                            offchip == True
                            break
                        if ordcen[pcen,nm]+k >= msarc.shape[0]:
                            mtfit[k+sz] = 1.0
                            offchip = True
                            break
                        xfit = np.arange(pcen-3, pcen+3+1, 1.0)  # 2x3 + 1 = 7 pixels total along the spectral dimension
                        yfit = msarc[ordcen[pcen,nm]+k,pcen-3:pcen+3+1]
                        if np.size(yfit) == 0:
                            mtfit[k+sz] = 1.0
                            offchip = True
                            break
                        #params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
                        params, fail = arutils.gauss_fit(xfit,yfit,pcen)
                        ytfit[k+sz] = params[1]
                        etfit[k+sz] = 0.02
                        if fail: mtfit[k+sz] = 1.0
                        else: pcen = int(0.5+params[1])
                    if offchip: continue
                    # Fit down
                    pcen = int(0.5+ytfit[sz]) # Start with the best-fitting centroid at arccen
                    for k in xrange(1,sz+1):
                        if (pcen < 3) or (pcen > msarc.shape[0]-4):
                            offchip == True
                            break
                        if ordcen[pcen,nm]-k < 0:
                            mtfit[sz-k] = 1.0
                            offchip = True
                            break
                        xfit = np.arange(pcen-3, pcen+3+1, 1.0)
                        yfit = msarc[ordcen[pcen,nm]-k,pcen-3:pcen+3+1]
                        if np.size(yfit) == 0:
                            mtfit[sz-k] = 1.0
                            offchip = True
                            break
                        #params, perror, fail = arfitbase.fit_gauss(xfit,yfit)
                        params, fail = arutils.gauss_fit(xfit,yfit,pcen)
                        ytfit[sz-k] = params[1]
                        etfit[sz-k] = 0.02
                        if fail: mtfit[sz-k] = 1.0
                        else: pcen = int(0.5+params[1])
                    if offchip: continue
                wmask = np.where(mtfit==0.0)
#				try:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1,w=1.0/mt)
#				except:
#					tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1)
                null, tcoeff = arutils.robust_polyfit(xtfit[wmask], ytfit[wmask], slf._argflag['trace']['orders']['tiltdisporder'], sigma=2.0)
                #tcoeff = np.polynomial.polynomial.polyfit(xtfit[wmask],ytfit[wmask],1)
                # Save the tilt angle
                tiltang[j,i] = tcoeff[1] # tan(tilt angle)
                centval[j,i] = tcoeff[0] # centroid of arc line
            nm += 1
        msgs.info("Fitting tilt angles")
        tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
        maskord = np.where(maskorder==1)[0]
        extrap_ord = np.zeros(maskorder.size)
        for o in xrange(maskorder.size):
            if o in maskord:
                extrap_ord[o] = 1
                continue
            w = np.where(tiltang[:,o]!=-999999.9)
            if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
                extrap_ord[o] = 1.0
                maskord = np.append(maskord,o)
            else:
                null, tempc = arutils.robust_polyfit(centval[:,o][w],tiltang[:,o][w], slf._argflag['trace']['orders']['tiltdisporder'], function=slf._argflag['trace']['orders']['function'],sigma=2.0,
                                                     minv=0.0, maxv=msarc.shape[slf._dispaxis]-1)
                tcoeff[:,o] = tempc
#				tcoeff[:,o] = arutils.func_fit(centval[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
#				plt.clf()
#				plt.plot(centval[:,o][w],tiltang[:,o][w],'bx')
#				xmod = np.linspace(np.min(centval[:,o][w]),np.max(centval[:,o][w]),1000)
#				ymod = arutils.func_val(tcoeff[:,o],xmod,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
#				plt.plot(xmod,ymod,'r-')
#		plt.show()
        maskord.sort()
        xv = np.arange(msarc.shape[slf._dispaxis])
        tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'], minv=0.0,
                                   maxv=msarc.shape[slf._dispaxis]-1).T
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        ofit = slf._argflag['trace']['orders']['pcatilt']
        lnpc = len(ofit)-1
        if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
            # Perform a PCA on the tilts
            msgs.info("Performing a PCA on the order edges")
            ordsnd = np.arange(tiltang.shape[1])+1.0
            xcen = xv[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
            fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=slf._argflag['trace']['orders']['function'])
            # If the PCA worked OK, do the following
            msgs.work("Should something be done here inbetween the two basis calls?")
            fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=slf._argflag['trace']['orders']['function'])
            arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
            # Extrapolate the remaining orders requested
            orders = 1.0+np.arange(maskorder.size)
            extrap_tilt, outpar = arpca.extrapolate(outpar, orders, msgs, function=slf._argflag['trace']['orders']['function'])
            tilts = extrap_tilt
            arpca.pc_plot_arctilt(tiltang, centval, tilts, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
        else:
            msgs.warn("Could not perform a PCA when tracing the order tilts"+msgs.newline()+"Not enough well-traced orders")
            msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
            xtiltfit = np.array([])
            ytiltfit = np.array([])
            for o in xrange(tiltang.shape[1]):
                w = np.where(tiltang[:,o]!=-999999.9)
                if np.size(w[0]) != 0:
                    xtiltfit = np.append(xtiltfit,centval[:,o][w])
                    ytiltfit = np.append(ytiltfit,tiltang[:,o][w])
            if np.size(xtiltfit) > slf._argflag['trace']['orders']['tiltdisporder']+2:
                tcoeff = arutils.func_fit(xtiltfit,ytiltfit,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],
                                          minv=0.0, maxv=msarc.shape[slf._dispaxis]-1)
                tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'], minv=0.0,
                                           maxv=msarc.shape[slf._dispaxis]-1)
                tilts = tiltval[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
            else:
                msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
                msgs.info("Assuming there is no tilt")
                tilts = np.zeros_like(slf._lordloc)
    elif slf._argflag['trace']['orders']['tilts'] == 'fit2D':
        # Go along each order and fit the tilts in 2D
        tiltang=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        msgs.work("This next step could be multiprocessed to speed up the reduction")

        for i in xrange(arccen.shape[1]):
            pixt = arcdet[np.where(arcdet[:,i]!=-1)[0],i]
            for j in xrange(pixt.size):
#				plt.plot(pixcen, arccen[:,i], 'k-', drawstyle='steps')
#				plt.plot([pixt[j],pixt[j]],[0.0,70000.0],'r-')
#				plt.show()
#				plt.clf()
                # Extract a small region
                sz = int(np.floor(np.abs(slf._rordloc[pixt[j],i]-slf._lordloc[pixt[j],i])/2.0))-1
                if slf._dispaxis == 0:
                    minx, maxx = pixt[j]-4, pixt[j]+4+1
                    miny, maxy = ordcen[pixt[j],i]-sz, ordcen[pixt[j],i]+sz+1
                    # Only use arc lines that are entirely on the chip
                    if (minx < 0) or (maxx > msarc.shape[0]) or (miny < 0) or (maxy > msarc.shape[1]): continue
# 					if minx < 0: minx = 0
# 					elif maxx > msarc.shape[0]: maxx = msarc.shape[0]
# 					if miny < 0: miny = 0
# 					elif maxy > msarc.shape[1]: maxy = msarc.shape[1]
                    sqrdata = msarc[minx:maxx,miny:maxy]
                    # Fit the 2D square region
                    pos = [miny-ordcen[pixt[j],i],minx-pixt[j],maxy-ordcen[pixt[j],i],maxx-pixt[j]]
#					pos = [minx-pixt[j],miny-ordcen[pixt[j],i],maxx-pixt[j],maxy-ordcen[pixt[j],i]]
#					angle, error, fail = arfitbase.fit_tilt(np.rot90(sqrdata), pos, -derv[pixt[j],i])
                    angle, error, fail = arfitbase.fit_tilt(np.rot90(sqrdata), pos, 0.0)
                    if not fail: angle *= -1.0
                else:
                    minx, maxx = ordcen[pixt[j],i]-sz, ordcen[pixt[j],i]+sz+1
                    miny, maxy = pixt[j]-4, pixt[j]+4+1
                    if (minx < 0) or (maxx > msarc.shape[0]) or (miny < 0) or (maxy > msarc.shape[1]): continue
# 					if minx < 0: minx = 0
# 					elif maxx > msarc.shape[0]: maxx = msarc.shape[0]
# 					if miny < 0: miny = 0
# 					elif maxy > msarc.shape[1]: maxy = msarc.shape[1]
                    sqrdata = msarc[minx:maxx,miny:maxy]
                    # Fit the 2D square region
                    pos = [minx-ordcen[pixt[j],i],miny-pixt[j],maxx-ordcen[pixt[j],i],maxy-pixt[j]]
#					angle, error, fail = arfitbase.fit_tilt(sqrdata, pos, derv[pixt[j],i])
                    angle, error, fail = arfitbase.fit_tilt(sqrdata, pos, 0.0)
                # Save the tilt angle
                if not fail:
                    tiltang[j,i] = angle
                    centval[j,i] = float(pixt[j])
            w = np.where(tiltang[:,i]!=-999999.9)
            plt.plot(arcdet[:,i][w],np.arctan(tiltang[:,i][w])*180.0/np.pi,'bx')
            plt.show()
            plt.clf()
        # Perform a PCA on the tilts
        tiltang[w] = derv[w]
        msgs.info("Fitting tilt angles")
        tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
        ogd = np.ones(tiltang.shape[1])
        for o in xrange(tiltang.shape[1]):
            w = np.where(tiltang[:,o]!=-999999.9)
            if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
                ogd[o] = 0.0
            tcoeff[:,o] = arutils.func_fit(arcdet[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],
                                           minv=0.0, maxv=msarc.shape[slf._dispaxis]-1)
        gd = np.where(ogd==1.0)[0]
        xv = np.arange(msarc.shape[slf._dispaxis])
        tiltval = arutils.func_val(tcoeff[:,gd],xv,slf._argflag['trace']['orders']['function'], minv=0.0,
                                   maxv=msarc.shape[slf._dispaxis]-1).T
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the tilt angles")
        ofit = slf._argflag['trace']['orders']['pcatilt']
        lnpc = len(ofit)-1
        msgs.work("May need to do a check here to make sure tilt ofit is reasonable")
        coeffs = arutils.func_fit(xv,tiltval,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],
                                  minv=0.0, maxv=msarc.shape[slf._dispaxis])
        xcen = xv[:,np.newaxis].repeat(gd.size,axis=1)
        fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=True,function=slf._argflag['trace']['orders']['function'])
        # If the PCA worked OK, do the following
        msgs.work("Should something be done here inbetween the two basis calls?")
        fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=False,function=slf._argflag['trace']['orders']['function'])
        arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
        # Extrapolate the remaining orders requested
        orders = 1.0+np.arange(arcdet.shape[1])
        extrap_tilt, outpar = arpca.extrapolate(outpar, orders, function=slf._argflag['trace']['orders']['function'])
        tilts = extrap_tilt

    elif slf._argflag['trace']['orders']['tilts'] == 'trace':
        osiz = 0.5*(slf._rordloc-slf._lordloc)
        ordsiz = np.round(np.abs(0.5*(slf._rordloc-slf._lordloc))).astype(np.int)
        # Try to find a trace at every arc line
        msgs.info("Tracing tilts")
        w = np.where(maskorder==0)[0]
        ttiltang, tcentval = arcytrace.trace_tilts(msarc,ordcen,ordsiz[:,w],arcdet,slf._dispaxis,1+2*np.max(ordsiz))
        wB = np.where(ttiltang!=-999999.9)
        centval=-999999.9*np.ones((tcentval.shape[0],maskorder.size))
        tiltang=-999999.9*np.ones((ttiltang.shape[0],maskorder.size))
        nm=0
        for i in xrange(maskorder.size):
            if maskorder[i] == 1: continue
            tx, tya, tyb = np.arange(tcentval.shape[0]), np.array([nm]), np.array([i])
            xya = np.ix_(tx,tya)
            xyb = np.ix_(tx,tyb)
            centval[xyb] = tcentval[xya]
            tiltang[xyb] = ttiltang[xya]
            nm += 1
        # Convert this to tilt angles and fit the tilts along the orders
        #msgs.info("Plotting tilts")
        #zmin, zmax = arplot.zscale(msarc)
        #implot = plt.imshow(msarc, extent=(0, msarc.shape[1], 0, msarc.shape[0]), origin='lower', interpolation='none', aspect='auto')
        #implot.set_cmap("gray")
        #implot.set_clim(zmin,zmax)
        #for o in xrange(derv.shape[1]):
        #	for l in xrange(derv.shape[0]):
        #		if derv[l,o] == -999999.9: break
        #		yplt = np.arange(-5,6)+ocen[arcdet[l,o],o]
        #		xplt = arcdet[l,o] + derv[l,o]*(yplt-ocen[arcdet[l,o],o])
        #		plt.plot(xplt,yplt,'r-')
        #plt.show()
        #plt.clf()
        msgs.info("Calculating tilt angles")
        #tiltang = -999999.9*np.ones_like(derv)
        #w = np.where(derv!=-999999.9)
#		tiltang[w] = np.arctan(derv[w])*180.0/np.pi
        #tiltang[w] = derv[w]
        msgs.info("Fitting tilt angles")
        tcoeff = np.zeros((slf._argflag['trace']['orders']['tiltdisporder']+1,maskorder.size))
        maskord = np.array([],dtype=np.int)
        extrap_ord = np.zeros(maskorder.size)
#		o = 0
        for ow in xrange(maskorder.size):
            if maskorder[ow] == 1:
                maskord = np.append(maskord,ow)
                extrap_ord[ow] = 1.0
                continue
            w = np.where(tiltang[:,ow]!=-999999.9)
            if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
                extrap_ord[ow] = 1.0
                maskord = np.append(maskord,ow)
            else:
                tcoeff[:,ow] = arutils.func_fit(centval[:,ow][w],tiltang[:,ow][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],
                                                minv=0.0, maxv=msarc.shape[slf._dispaxis]-1)
#				o += 1
        xv = np.arange(msarc.shape[slf._dispaxis])
        tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'], minv=0.0,
                                   maxv=msarc.shape[slf._dispaxis]-1).T
        plt.clf()
        for ow in xrange(maskorder.size):
            if maskorder[ow] == 1:
                continue
            w = np.where(tiltang[:,ow]!=-999999.9)
            if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
                continue
            plt.plot(centval[:,ow][w],tiltang[:,ow][w],'bx')
            plt.plot(xv,tiltval[:,ow],'r-')
        plt.show()
        msgs.work("May need to do a check here to make sure ofit is reasonable")
        maskord.sort()
        ofit = slf._argflag['trace']['orders']['pcatilt']
        lnpc = len(ofit)-1
        if np.sum(1.0-extrap_ord) > ofit[0]+1: # Only do a PCA if there are enough good orders
            # Perform a PCA on the tilts
            msgs.info("Performing a PCA on the order tilts")
            ordsnd = np.arange(maskorder.size)+1.0
            xcen = xv[:,np.newaxis].repeat(maskorder.size,axis=1)
            fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=True,function=slf._argflag['trace']['orders']['function'])
            # If the PCA worked OK, do the following
            msgs.work("Should something be done here inbetween the two basis calls?")
            fitted, outpar = arpca.basis(xcen,tiltval,tcoeff,lnpc,ofit,x0in=ordsnd,mask=maskord,skipx0=False,function=slf._argflag['trace']['orders']['function'])
            arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts", prefix=prefix)
            # Extrapolate the remaining orders requested
            orders = 1.0+np.arange(maskorder.size)
            extrap_tilt, outpar = arpca.extrapolate(outpar, orders, msgs, function=slf._argflag['trace']['orders']['function'])
            tilts = extrap_tilt
        else:
            msgs.warn("Could not perform a PCA when tracing the order tilts"+msgs.newline()+"Not enough well-traced orders")
            msgs.info("Attempting to fit tilts by assuming the tilt is order-independent")
            xtiltfit = np.array([])
            ytiltfit = np.array([])
            for o in xrange(tiltang.shape[1]):
                w = np.where(tiltang[:,o]!=-999999.9)
                if np.size(w[0]) != 0:
                    xtiltfit = np.append(xtiltfit,centval[:,o][w])
                    ytiltfit = np.append(ytiltfit,tiltang[:,o][w])
            if np.size(xtiltfit) > slf._argflag['trace']['orders']['tiltdisporder']+2:
                tcoeff = arutils.func_fit(xtiltfit,ytiltfit,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],
                                          minv=0.0, maxv=msarc.shape[slf._dispaxis]-1)
                tiltval = arutils.func_val(tcoeff,xv,slf._argflag['trace']['orders']['function'], minv=0.0,
                                           maxv=msarc.shape[slf._dispaxis]-1)
                tilts = tiltval[:,np.newaxis].repeat(tiltang.shape[1],axis=1)
            else:
                msgs.warn("There are still not enough detections to obtain a reliable tilt trace")
                msgs.info("Assuming there is no tilt")
                tilts = np.zeros_like(slf._lordloc)
        """
        msgs.info("Fitting tilt angles")
        tcoeff = np.ones((slf._argflag['trace']['orders']['tiltdisporder']+1,tiltang.shape[1]))
        ogd = np.ones(tiltang.shape[1])
        for o in xrange(tiltang.shape[1]):
            w = np.where(tiltang[:,o]!=-999999.9)
            if np.size(w[0]) <= slf._argflag['trace']['orders']['tiltdisporder']+2:
                ogd[o] = 0.0
                continue
            tcoeff[:,o] = arutils.func_fit(arcdet[:,o][w],tiltang[:,o][w],slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis]-1)
        gd = np.where(ogd==1.0)[0]
        xv = np.arange(msarc.shape[slf._dispaxis])
        tiltval = arutils.func_val(tcoeff[:,gd],xv,slf._argflag['trace']['orders']['function'],min=0.0,max=msarc.shape[slf._dispaxis]-1).T
        # Perform a PCA on the tilts
        msgs.info("Performing a PCA on the tilt angles")
        ofit = slf._argflag['trace']['orders']['pcatilt']
        lnpc = len(ofit)-1
        msgs.bug("May need to do a check here to make sure tilt ofit is reasonable")
        coeffs = arutils.func_fit(xv,tiltval,slf._argflag['trace']['orders']['function'],slf._argflag['trace']['orders']['tiltdisporder'],min=0.0,max=msarc.shape[slf._dispaxis])
        xcen = xv[:,np.newaxis].repeat(gd.size,axis=1)
        fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=True,function=slf._argflag['trace']['orders']['function'])
        # If the PCA worked OK, do the following
        msgs.bug("Should something be done here inbetween the two basis calls?")
        fitted, outpar = arpca.basis(xcen,tiltval,coeffs,lnpc,ofit,x0in=gd+1.0,skipx0=False,function=slf._argflag['trace']['orders']['function'])
        arpca.pc_plot(outpar, ofit, plotsdir=slf._argflag['run']['plotsdir'], pcatype="tilts")
        # Extrapolate the remaining orders requested
        orders = 1.0+np.arange(arcdet.shape[1])
        extrap_tilt, outpar = arpca.extrapolate(outpar,orders,function=slf._argflag['trace']['orders']['function'])
        tilts = extrap_tilt
        """
    else:
        msgs.warn("I don't know how to deal with '{0:s}' tilts".format(slf._argflag['trace']['orders']['tilts']))
        msgs.info("Assuming there is no tilt")
        centval=-999999.9*np.ones((arcdet.shape[0],maskorder.size))
        nm=0
        for i in xrange(maskorder.size):
            if maskorder[i] == 1: continue
            w = arcdet[np.where(arcdet[:,nm]!=-1)[0],nm]
            if np.size(w) != 0: centval[:np.size(w),i] = w
            nm += 1
        tilts = np.zeros_like(slf._lordloc)
    # Write out ds9 regions file for slit tilts.
    msgs.info("Writing QC files")
    if tltprefix != "":
        tracereg = open("{0:s}/{1:s}_trace_tilts.reg".format(slf._argflag['run']['plotsdir'],tltprefix),'w')
    else:
        tracereg = open("{0:s}/trace_tilts.reg".format(slf._argflag['run']['plotsdir']),'w')
    tracereg.write("# Region file format: DS9 version 4.1\n")
    tracereg.write('global color=red dashlist=8 3 width=3 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    tracereg.write("image\n")
    # Correct the real pixel locations to image pixel locations
    cv_corr = (centval+0.5).astype(np.int)
    nm = 0
    for i in xrange(tilts.shape[1]):
        if maskorder[i] == 1: continue
        for j in xrange(centval.shape[0]):
            if centval[j,i] == -999999.9: break
            incpt = cv_corr[j,i] - tilts[cv_corr[j,i],i]*ordcen[cv_corr[j,i],nm]
            xmin = ordcen[cv_corr[j,i],nm] - ordwid[cv_corr[j,i],i]
            xmax = ordcen[cv_corr[j,i],nm] + ordwid[cv_corr[j,i],i]
            ymin = incpt + xmin*tilts[cv_corr[j,i],i]
            ymax = incpt + xmax*tilts[cv_corr[j,i],i]
            if slf._dispaxis == 0:
                tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(xmin+1,ymin+1,xmax+1,ymax+1))
            else:
                tracereg.write('line({0:.2f},{1:.2f},{2:.2f},{3:.2f}) # line=0 0\n'.format(ymin+1,xmin+1,ymax+1,xmax+1))
        nm += 1
    tracereg.close()
    # Plot the tilts in real time if the user requests
    if slf._argflag['run']['qcontrol']:
        import ds9
        # Set up a ds9 instance
        d = ds9.ds9()
        # Load the image
        d.set_np2arr(msarc)
        # Zoom to fit
        d.set('zoom to fit')
        # Change the colormap and scaling
        d.set('cmap gray')
        d.set('scale log')
        # Plot the regions
        if tltprefix != "":
            d.set('regions load ' + '"' + '{0:s}/{1:s}_trace_tilts.reg'.format(slf._argflag['run']['plotsdir'], tltprefix) + '"')
            if trcprefix != "":
                tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
                if os.path.exists(tfil):
                    d.set('regions load ' + '"' + tfil + '"')
                else:
                    msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
            else:
                tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
                if os.path.exists(tfil):
                    d.set('regions load ' + '"' + tfil + '"')
                else:
                    msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
        else:
            d.set('regions load ' + '"' + '{0:s}/trace_tilts.reg'.format(slf._argflag['run']['plotsdir']) + '"')
            if trcprefix != "":
                tfil = '{0:s}/{1:s}_trace_orders.reg'.format(slf._argflag['run']['plotsdir'], trcprefix)
                if os.path.exists(tfil):
                    d.set('regions load ' + '"' + tfil + '"')
                else:
                    msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
            else:
                tfil = '{0:s}/trace_orders.reg'.format(slf._argflag['run']['plotsdir'])
                if os.path.exists(tfil):
                    d.set('regions load ' + '"' + tfil + '"')
                else:
                    msgs.warn("Couldn't open order locations DS9 region file:"+msgs.newline()+tfil)
        # Save the image
        if slf._argflag['run']['stopcheck']:
            null=raw_input(msgs.input()+"Press enter to continue...")
        else:
            msgs.info("DS9 window was updated")
    return tilts, satsnd


def gen_pixloc(slf, frame, det, gen=True):
    """
    Generate an array of physical pixel coordinates

    Parameters
    ----------
    slf : class
      Science Exposure class
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
    msgs.info("Deriving physical pixel locations on the detector")
    locations = np.zeros((frame.shape[0],frame.shape[1],4))
    if gen:
        msgs.info("Pixel gap in the dispersion direction = {0:4.3f}".format(slf._spect['det'][det-1]['xgap']))
        msgs.info("Pixel size in the dispersion direction = {0:4.3f}".format(1.0))
        xs = np.arange(frame.shape[slf._dispaxis]*1.0)*slf._spect['det'][det-1]['xgap']
        xt = 0.5 + np.arange(frame.shape[slf._dispaxis]*1.0) + xs
        msgs.info("Pixel gap in the spatial direction = {0:4.3f}".format(slf._spect['det'][det-1]['ygap']))
        msgs.info("Pixel size in the spatial direction = {0:4.3f}".format(slf._spect['det'][det-1]['ysize']))
        ys = np.arange(frame.shape[1-slf._dispaxis])*slf._spect['det'][det-1]['ygap']*slf._spect['det'][det-1]['ysize']
        yt = slf._spect['det'][det-1]['ysize']*(0.5 + np.arange(frame.shape[1-slf._dispaxis]*1.0)) + ys
        xloc, yloc = np.meshgrid(xt, yt)
#		xwid, ywid = np.meshgrid(xs,ys)
        msgs.info("Saving pixel locations")
        if slf._dispaxis == 0:
            locations[:,:,0] = xloc.T
            locations[:,:,1] = yloc.T
            locations[:,:,2] = 1.0
            locations[:,:,3] = slf._spect['det'][det-1]['ysize']
        else:
            locations[:,:,0] = xloc
            locations[:,:,1] = yloc
            locations[:,:,2] = 1.0
            locations[:,:,3] = slf._spect['det'][det-1]['ysize']
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

    if axis == 0:
        diff = pixlocn[:,0,0]
    else:
        diff = pixlocn[0,:,1]
    # if axis == dispaxis:
    #     if axis == 0:
    #         diff = pixlocn[:,0,0]
    #     else:
    #         diff = pixlocn[0,:,0]
    # else:
    #     if axis == 0:
    #         diff = pixlocn[:,0,1]
    #     else:
    #         diff = pixlocn[0,:,1]
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
