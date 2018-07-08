from __future__ import (print_function, absolute_import, division, unicode_literals)

import os
import time
import inspect

import numpy as np
from matplotlib import pyplot as plt
from matplotlib import gridspec

#from pypit import arparse as settings
from pypit import arparse
from pypit import msgs
from pypit import arutils
from pypit import ararclines
from pypit import arqa
from pypit import arpixels
from pypit import ardebug as debugger


def setup_param(spectro_class, msarc_shape, fitstbl, arc_idx,
                calibrate_lamps=None):
    """ Setup for arc analysis

    Parameters
    ----------
    spectrograph : str
    msarc_shape : tuple
    fitstbl : Table
      Contains relevant information from fits header files
    arc_idx : int
      Index of the arc frame in the fitstbl
    calibrate_lamps : str, optional
       List of lamps used

    """
    # Defaults
    arcparam = dict(llist='',
                    disp=0.,             # Ang/unbinned pixel
                    b1=0.,               # Pixel fit term (binning independent)
                    b2=0.,               # Pixel fit term
                    lamps=[],            # Line lamps on
                    wv_cen=0.,           # Estimate of central wavelength
                    wvmnx=[2900.,12000.],# Guess at wavelength range
                    disp_toler=0.1,      # 10% tolerance
                    match_toler=3.,      # Matching tolerance (pixels)
                    min_ampl=300.,       # Minimum amplitude
                    func='legendre',     # Function for fitting
                    n_first=1,           # Order of polynomial for first fit
                    n_final=4,           # Order of polynomial for final fit
                    nsig_rej=2.,         # Number of sigma for rejection
                    nsig_rej_final=3.0,  # Number of sigma for rejection (final fit)
                    Nstrong=13)          # Number of lines for auto-analysis

    # Instrument/disperser specific
    disperser = fitstbl["dispname"][arc_idx]
    binspatial, binspectral = arparse.parse_binning(fitstbl['binning'][arc_idx])
    modify_dict = spectro_class.setup_arcparam(arcparam, disperser=disperser, fitstbl=fitstbl,
                                               arc_idx=arc_idx, binspatial=binspatial,
                                               binspectral=binspectral, msarc_shape=msarc_shape)
    # Load linelist
    #if settings.argflag['arc']['calibrate']['lamps'] is not None:
    if calibrate_lamps is not None:
        arcparam['lamps'] = calibrate_lamps
    slmps = arcparam['lamps'][0]
    for lamp in arcparam['lamps'][1:]:
        slmps=slmps+','+lamp
    msgs.info('Loading line list using {:s} lamps'.format(slmps))

    arcparam['llist'] = ararclines.load_arcline_list(arcparam['lamps'], disperser,
                                                     spectro_class.spectrograph,
                                                     wvmnx=arcparam['wvmnx'],
                                                     modify_parse_dict=modify_dict)
    # Binning
    arcparam['disp'] *= binspectral

    # Return
    return arcparam

def get_censpec(lordloc, rordloc, pixlocn, frame, det, nonlinear_counts=None, gen_satmask=False):
    """ Extract a simple spectrum down the center of each slit
    Parameters
    ----------
    frame : ndarray
      Image
    det : int
    gen_satmask : bool, optional
      Generate a saturation mask?

    Returns
    -------
    arccen : ndarray
      Extracted arcs.  This *need* not be one per slit/order,
      although I wish it were (with `rejected` ones padded with zeros)
    maskslit : bool array
      1 = Bad slit/order for extraction (incomplete)
      0 = Ok
    satmask : ndarray
      Saturation mask
      None if gen_satmask=False
    """
    dnum = arparse.get_dnum(det)

    ordcen = 0.5*(lordloc+rordloc)
    ordwid = 0.5*np.abs(lordloc-rordloc)
    satsnd = None
    if gen_satmask and nonlinear_counts is not None:
        msgs.info("Generating a mask of arc line saturation streaks")
        satmask = saturation_mask(frame, nonlinear_counts)
        satsnd = order_saturation(satmask, (ordcen+0.5).astype(int), (ordwid+0.5).astype(int))

    # Extract a rough spectrum of the arc in each slit
    msgs.info("Extracting an approximate arc spectrum at the centre of each slit")
    tordcen = None
    maskslit = np.zeros(ordcen.shape[1], dtype=np.int)
    for i in range(ordcen.shape[1]):
        wl = np.size(np.where(ordcen[:, i] <= pixlocn[0,0,1])[0])
        wh = np.size(np.where(ordcen[:, i] >= pixlocn[0,-1,1])[0])
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
        ordcen = arpixels.phys_to_pix(ordcen, pixlocn, 1)
    else:
        ordcen = arpixels.phys_to_pix(tordcen[:,w], pixlocn, 1)

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
    # Construct the good ones
    gd_arccen = (1.0/5.0) * (frame[temparr, ordcen] +
                             frame[temparr, op1] + frame[temparr, op2] +
                             frame[temparr, om1] + frame[temparr, om2])
    # Pad masked ones with zeros
    if np.sum(maskslit) > 0:
        arccen = np.zeros((gd_arccen.shape[0], maskslit.size))
        gd = maskslit == 0
        arccen[:,gd] = gd_arccen
    else:
        arccen = gd_arccen
    del temparr

    return arccen, maskslit, satsnd

def detect_lines(censpec, nfitpix=5, nonlinear=None):
    """
    Extract an arc down the center of the chip and identify
    statistically significant lines for analysis.

    Parameters
    ----------
    det : int
      Index of the detector
    msarc : ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
    censpec : ndarray, optional
      A 1D spectrum to be searched for significant detections

    Returns
    -------
    tampl : ndarray
      The amplitudes of the line detections
    tcent : ndarray
      The centroids of the line detections
    twid : ndarray
      The 1sigma Gaussian widths of the line detections
    w : ndarray
      An index array indicating which detections are the most reliable.
    detns : ndarray
      The spectrum used to find detections. This spectrum has
      had any "continuum" emission subtracted off
    """
    # TODO I don't see that there is any continuum subtraction being done here contrary to what the docs say.

    # Extract a rough spectrum of the arc in each order
    msgs.info("Detecting lines")
    '''
    if MK_SATMASK:
        ordwid = 0.5*np.abs(slf._lordloc[det-1] - slf._rordloc[det-1])
        msgs.info("Generating a mask of arc line saturation streaks")
        satmask = saturation_mask(msarc, nonlinear)
        satsnd = order_saturation(satmask, ordcen, (ordwid+0.5).astype(np.int))
    else:
        satsnd = np.zeros_like(ordcen)
    '''
    # Detect the location of the arc lines
    msgs.info("Detecting the strongest, nonsaturated lines")

    fitp = nfitpix # settings.argflag['arc']['calibrate']['nfitpix']
    if len(censpec.shape) == 3:
        detns = censpec[:, 0].flatten()
    else:
        detns = censpec.copy()
    detns = detns.astype(np.float)
    xrng = np.arange(detns.size, dtype=np.float)

    # Find all significant detections
    # TODO -- Need to add nonlinear back in here
    pixt = np.where((detns > 0.0) &  # (detns < slf._nonlinear[det-1]) &
                    (detns > np.roll(detns, 1)) & (detns >= np.roll(detns, -1)) &
                    (np.roll(detns, 1) > np.roll(detns, 2)) & (np.roll(detns, -1) > np.roll(detns, -2)) &#)[0]
                    (np.roll(detns, 2) > np.roll(detns, 3)) & (np.roll(detns, -2) > np.roll(detns, -3)))[0]
#                    (np.roll(detns, 3) > np.roll(detns, 4)) & (np.roll(detns, -3) > np.roll(detns, -4)) & # )[0]
#                    (np.roll(detns, 4) > np.roll(detns, 5)) & (np.roll(detns, -4) > np.roll(detns, -5)))[0]
    tampl, tcent, twid = fit_arcspec(xrng, detns, pixt, fitp)
    w = np.where((~np.isnan(twid)) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent > 0.0) & (tcent < xrng[-1]))
    # Check the results
    #plt.clf()
    #plt.plot(xrng,detns,'k-')
    #plt.plot(tcent,tampl,'ro')
    #plt.show()
    # Return
    return tampl, tcent, twid, w, detns


def fit_arcspec(xarray, yarray, pixt, fitp):

    # Setup the arrays with fit parameters
    sz_p = pixt.size
    sz_a = yarray.size
    ampl, cent, widt = -1.0*np.ones(sz_p, dtype=np.float),\
                       -1.0*np.ones(sz_p, dtype=np.float),\
                       -1.0*np.ones(sz_p, dtype=np.float)

    for p in range(sz_p):
        pmin = pixt[p]-(fitp-1)//2
        pmax = pixt[p]-(fitp-1)//2 + fitp
        if pmin < 0:
            pmin = 0
        if pmax > sz_a:
            pmax = sz_a
        if pmin == pmax:
            continue
        if pixt[p]-pmin <= 1 or pmax-pixt[p] <= 1:
            continue  # Probably won't be a good solution
        # Fit the gaussian
        try:
            popt = arutils.func_fit(xarray[pmin:pmax], yarray[pmin:pmax], "gaussian", 3)
            ampl[p] = popt[0]
            cent[p] = popt[1]
            widt[p] = popt[2]
        except RuntimeError:
            pass
    return ampl, cent, widt


def simple_calib(msarc, aparm, censpec, nfitpix=5, get_poly=False,
                 IDpixels=None, IDwaves=None):
    """Simple calibration algorithm for longslit wavelengths

    Uses slf._arcparam to guide the analysis

    Parameters
    ----------
    msarc : ndarray
    aparm : dict
    censpec : ndarray
    get_poly : bool, optional
      Pause to record the polynomial pix = b0 + b1*lambda + b2*lambda**2

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """

    # Extract the arc
    msgs.work("Detecting lines..")
    tampl, tcent, twid, w, yprep = detect_lines(censpec, nfitpix=nfitpix)

    # Cut down to the good ones
    tcent = tcent[w]
    tampl = tampl[w]
    msgs.info('Detected {:d} lines in the arc spectrum.'.format(len(w[0])))

    # Read Arc linelist
    llist = aparm['llist']

    # IDs were input by hand
    if IDpixels is not None:
        # Check that there are at least 5 values
        pixels = np.array(IDpixels) # settings.argflag['arc']['calibrate']['IDpixels'])
        if np.sum(pixels > 0.) < 5:
            msgs.error("Need to give at least 5 pixel values!")
        #
        msgs.info("Using input lines to seed the wavelength solution")
        # Calculate median offset
        mdiff = [np.min(np.abs(tcent-pix)) for pix in pixels]
                 #settings.argflag['arc']['calibrate']['IDpixels']]
        med_poff = np.median(np.array(mdiff))
        msgs.info("Will apply a median offset of {:g} pixels".format(med_poff))

        # Match input lines to observed spectrum
        nid = pixels.size # len(settings.argflag['arc']['calibrate']['IDpixels'])
        idx_str = np.ones(nid).astype(int)
        ids = np.zeros(nid)
        idsion = np.array(['     ']*nid)
        gd_str = np.arange(nid).astype(int)
        for jj,pix in enumerate(pixels): #settings.argflag['arc']['calibrate']['IDpixels']):
            diff = np.abs(tcent-pix-med_poff)
            if np.min(diff) > 2.:
                debugger.set_trace()
                msgs.error("No match with input pixel {:g}!".format(pix))
            else:
                imn = np.argmin(diff)
            # Set
            idx_str[jj] = imn
            # Take wavelength from linelist instead of input value
            wdiff = np.abs(llist['wave']-IDwaves[jj]) # settings.argflag['arc']['calibrate']['IDwaves'][jj])
            imnw = np.argmin(wdiff)
            if wdiff[imnw] > 0.015:  # Arbitrary tolerance
                msgs.error("Input IDwaves={:g} is not in the linelist.  Fix".format(
                    IDwaves[jj]))
                        #settings.argflag['arc']['calibrate']['IDwaves'][jj]))
            else:
                ids[jj] = llist['wave'][imnw]
                idsion[jj] = llist['Ion'][imnw]
                msgs.info("Identifying arc line: {:s} {:g}".format(idsion[jj],ids[jj]))
    else:
        # Generate dpix pairs
        msgs.info("Using pair algorithm for wavelength solution")
        nlist = len(llist)
        dpix_list = np.zeros((nlist,nlist))
        for kk,row in enumerate(llist):
            #dpix_list[kk,:] = (np.array(row['wave'] - llist['wave']))/disp
            dpix_list[kk,:] = msarc.shape[0]*(aparm['b1']*(np.array(row['wave'] - llist['wave'])) + aparm['b2']*np.array(row['wave']**2 - llist['wave']**2) )

        # Lambda pairs for the strongest N lines
        srt = np.argsort(tampl)
        idx_str = srt[-aparm['Nstrong']:]
        idx_str.sort()
        dpix_obs = np.zeros((aparm['Nstrong'], aparm['Nstrong']))
        for kk,idx in enumerate(idx_str):
            dpix_obs[kk,:] = np.array(tcent[idx] - tcent[idx_str])

        # Match up (ugly loops)
        ids = np.zeros(aparm['Nstrong'])
        idsion = np.array(['     ']*aparm['Nstrong'])
        for kk in range(aparm['Nstrong']):
            med_off = np.zeros(nlist)
            for ss in range(nlist):
                dpix = dpix_list[ss]
                min_off = []
                for jj in range(aparm['Nstrong']):
                    min_off.append(np.min(np.abs(dpix_obs[kk,jj]-dpix)))
                med_off[ss] = np.median(min_off)
            # Set by minimum
            idm = np.argmin(med_off)
            ids[kk] = llist['wave'][idm]
            idsion[kk] = llist['Ion'][idm]

        # Calculate disp of the strong lines
        disp_str = np.zeros(aparm['Nstrong'])
        for kk in range(aparm['Nstrong']):
            disp_val = (ids[kk]-ids)/(tcent[idx_str[kk]]-tcent[idx_str])
            isf = np.isfinite(disp_val)
            disp_str[kk] = np.median(disp_val[isf])
        # Consider calculating the RMS with clipping
        gd_str = np.where( np.abs(disp_str-aparm['disp'])/aparm['disp'] < aparm['disp_toler'])[0]
        msgs.info('Found {:d} lines within the dispersion threshold'.format(len(gd_str)))
        if len(gd_str) < 5:
            if msgs._debug['arc']:
                msgs.warn('You should probably try your best to ID lines now.')
                debugger.set_trace()
                debugger.plot1d(yprep)
            else:
                msgs.error('Insufficient lines to auto-fit.')

    # Debug
    #debug=True
    if msgs._debug['arc']:
        #tmp = list(gd_str)
        #tmp.pop(1)
        #gd_str = np.array(tmp)
        #xdb.xpcol(tcent[idx_str[gd_str]],ids[gd_str])
        #xdb.xplot(tcent[idx_str[gd_str]],ids[gd_str],scatter=True)
        # debugger.xplot(yprep)
        debugger.set_trace()

    msgs.work('Cross correlate here?')

    # Setup for fitting
    ifit = idx_str[gd_str]
    sv_ifit = list(ifit) # Keep the originals
    all_ids = -999.*np.ones(len(tcent))
    all_idsion = np.array(['12345']*len(tcent))
    all_ids[ifit] = ids[gd_str]
    all_idsion[ifit] = idsion[gd_str]
    # Fit
    n_order = aparm['n_first']
    flg_quit = False
    fmin, fmax = -1., 1.
    msgs.info('Iterative wavelength fitting..')
    while (n_order <= aparm['n_final']) and (flg_quit is False):
        #msgs.info('n_order={:d}'.format(n_order))
        # Fit with rejection
        xfit, yfit = tcent[ifit], all_ids[ifit]
        mask, fit = arutils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej'], minv=fmin, maxv=fmax)
        # Reject but keep originals (until final fit)
        ifit = list(ifit[mask == 0]) + sv_ifit
        # Find new points (should we allow removal of the originals?)
        twave = arutils.func_val(fit, tcent, aparm['func'], minv=fmin, maxv=fmax)
        for ss,iwave in enumerate(twave):
            mn = np.min(np.abs(iwave-llist['wave']))
            if mn/aparm['disp'] < aparm['match_toler']:
                imn = np.argmin(np.abs(iwave-llist['wave']))
                #if msgs._debug['arc']:
                #    print('Adding {:g} at {:g}'.format(llist['wave'][imn],tcent[ss]))
                # Update and append
                all_ids[ss] = llist['wave'][imn]
                all_idsion[ss] = llist['Ion'][imn]
                ifit.append(ss)
        # Keep unique ones
        ifit = np.unique(np.array(ifit,dtype=int))
        #if msgs._debug['arc']:
        #    debugger.set_trace()
        # Increment order
        if n_order < aparm['n_final']:
            n_order += 1
        else:
            # This does 2 iterations at the final order
            flg_quit = True

    # Final fit (originals can now be rejected)
    fmin, fmax = 0., 1.
    xfit, yfit = tcent[ifit]/(msarc.shape[0]-1), all_ids[ifit]
    mask, fit = arutils.robust_polyfit(xfit, yfit, n_order, function=aparm['func'], sigma=aparm['nsig_rej_final'], minv=fmin, maxv=fmax)#, debug=True)
    irej = np.where(mask==1)[0]
    if len(irej) > 0:
        xrej = xfit[irej]
        yrej = yfit[irej]
        for imask in irej:
            msgs.info('Rejecting arc line {:g}'.format(yfit[imask]))
    else:
        xrej = []
        yrej = []
    xfit = xfit[mask==0]
    yfit = yfit[mask==0]
    ions = all_idsion[ifit][mask==0]
    #
    if msgs._debug['arc']:
        wave = arutils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
            'legendre', minv=fmin, maxv=fmax)
        debugger.set_trace()

        #debugger.xplot(xfit, np.ones(len(xfit)), scatter=True,
        #    xtwo=np.arange(msarc.shape[0]),ytwo=yprep)
        #debugger.xplot(xfit,yfit, scatter=True, xtwo=np.arange(msarc.shape[0]),
        #    ytwo=wave)
        #debugger.set_trace()
        #wave = arutils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
        #    'legendre', min=fmin, max=fmax)

    # 2nd order Poly fit for archival
    #get_poly=True
    if get_poly:
        poly_fit = arutils.func_fit(yfit,xfit, 'polynomial',2, minv=fmin, maxv=fmax)
        print(' Most likely you with to record these values:')
        print(poly_fit)
        debugger.set_trace()
    # Pack up fit
    final_fit = dict(fitc=fit, function=aparm['func'], xfit=xfit, yfit=yfit,
        ions=ions, fmin=fmin, fmax=fmax, xnorm=float(msarc.shape[0]),
        xrej=xrej, yrej=yrej, mask=mask, spec=yprep, nrej=aparm['nsig_rej_final'],
        shift=0., tcent=tcent)
    # QA
    #arc_fit_qa(slf, final_fit, slit)
    # RMS
    rms_ang = arutils.calc_fit_rms(xfit, yfit, fit, aparm['func'], minv=fmin, maxv=fmax)
    wave = arutils.func_val(fit, np.arange(msarc.shape[0])/float(msarc.shape[0]),
                            aparm['func'], minv=fmin, maxv=fmax)
    rms_pix = rms_ang/np.median(np.abs(wave-np.roll(wave,1)))
    msgs.info("Fit RMS = {} pix".format(rms_pix))
    # Return
    return final_fit


def calib_with_arclines(aparm, spec, use_method="general"):
    """Holy grail algorithms for wavelength calibration

    Uses arcparam to guide the analysis

    Parameters
    ----------
    aparm
    spec
    use_method : str, optional

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """
    import arclines.holy.grail
    # Extract the arc
    #msgs.work("Detecting lines")
    #tampl, tcent, twid, w, satsnd, spec = detect_lines( slf, det, msarc, censpec=censpec)

    if use_method == "semi-brute":
        best_dict, final_fit = arclines.holy.grail.semi_brute(spec, aparm['lamps'], aparm['wv_cen'], aparm['disp'], fit_parm=aparm, min_ampl=aparm['min_ampl'])
    elif use_method == "basic":
        stuff = arclines.holy.grail.basic(spec, aparm['lamps'], aparm['wv_cen'], aparm['disp'])
        status, ngd_match, match_idx, scores, final_fit = stuff
    else:
        # Now preferred
        best_dict, final_fit = arclines.holy.grail.general(spec, aparm['lamps'], fit_parm=aparm, min_ampl=aparm['min_ampl'])
    #
    return final_fit


def order_saturation(satmask, ordcen, ordwid):
    """
    .. todo::
        Document this!
    """
    sz_y, sz_x = satmask.shape
    sz_o = ordcen.shape[1]

    xmin = ordcen - ordwid
    xmax = ordcen + ordwid + 1
    xmin[xmin < 0] = 0
    xmax[xmax >= sz_x] = sz_x

    ordsat = np.zeros((sz_y, sz_o), dtype=int)
    for o in range(sz_o):
        for y in range(sz_y):
            ordsat[y,o] = (xmax[y,o] > xmin[y,o]) & np.any(satmask[y,xmin[y,o]:xmax[y,o]] == 1)

    return ordsat


def search_for_saturation_edge(a, x, y, sy, dx, satdown, satlevel, mask):
    sx = dx
    localx = a[x+sx,y+sy]
    while True:
        mask[x+sx,y+sy] = True
        sx += dx
        if x+sx > a.shape[0]-1 or x+sx < 0:
            break
        if a[x+sx,y+sy] >= localx/satdown and a[x+sx,y+sy]<satlevel:
            break
        localx = a[x+sx,y+sy]
    return mask


def determine_saturation_region(a, x, y, sy, dy, satdown, satlevel, mask):
    localy = a[x,y+sy]
    while True:
        mask[x,y+sy] = True
        mask = search_for_saturation_edge(a, x, y, sy, 1, satdown, satlevel, mask)
        mask = search_for_saturation_edge(a, x, y, sy, -1, satdown, satlevel, mask)

        sy += dy
        if y+sy > a.shape[1]-1 or y+sy < 0:
            return mask
        if a[x,y+sy] >= localy/satdown and a[x,y+sy] < satlevel:
            return mask
        localy = a[x,y+sy]


def saturation_mask(a, satlevel):
    """
    ... todo::
        Document this!
    """
    mask = np.zeros(a.shape, dtype=bool)
    a_is_saturated = a >= satlevel
    if not np.any(a_is_saturated):
        return mask.astype(int)

    satdown = 1.001
    sz_x, sz_y = a.shape

    for y in range (0,sz_y):
        for x in range(0,sz_x):
            if a_is_saturated[x,y] and not mask[x,y]:
                mask[x,y] = True
                mask = determine_saturation_region(a, x, y, 0, 1, satdown, satlevel, mask)
                mask = determine_saturation_region(a, x, y, -1, -1, satdown, satlevel, mask)

    return mask.astype(int)


def arc_fit_qa(setup, fit, slit, outfile=None, ids_only=False, title=None):
    """
    QA for Arc spectrum

    Parameters
    ----------
    fit : Wavelength fit
    arc_spec : ndarray
      Arc spectrum
    outfile : str, optional
      Name of output file
      or 'show' to show on screen
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    # Grab the named of the method
    method = inspect.stack()[0][3]
    # Outfil
    if outfile is None:
        outfile = arqa.set_qa_filename(setup, method, slit=slit)
    #
    arc_spec = fit['spec']

    # Begin
    if not ids_only:
        plt.figure(figsize=(8, 4.0))
        plt.clf()
        gs = gridspec.GridSpec(2, 2)
        idfont = 'xx-small'
    else:
        plt.figure(figsize=(11, 8.5))
        plt.clf()
        gs = gridspec.GridSpec(1, 1)
        idfont = 'small'

    # Simple spectrum plot
    ax_spec = plt.subplot(gs[:,0])
    ax_spec.plot(np.arange(len(arc_spec)), arc_spec)
    ymin, ymax = 0., np.max(arc_spec)
    ysep = ymax*0.03
    for kk, x in enumerate(fit['xfit']*fit['xnorm']):
        yline = np.max(arc_spec[int(x)-2:int(x)+2])
        # Tick mark
        ax_spec.plot([x,x], [yline+ysep*0.25, yline+ysep], 'g-')
        # label
        ax_spec.text(x, yline+ysep*1.3,
            '{:s} {:g}'.format(fit['ions'][kk], fit['yfit'][kk]), ha='center', va='bottom',
            size=idfont, rotation=90., color='green')
    ax_spec.set_xlim(0., len(arc_spec))
    ax_spec.set_ylim(ymin, ymax*1.2)
    ax_spec.set_xlabel('Pixel')
    ax_spec.set_ylabel('Flux')
    if title is not None:
        ax_spec.text(0.04, 0.93, title, transform=ax_spec.transAxes,
                     size='x-large', ha='left')#, bbox={'facecolor':'white'})
    if ids_only:
        plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
        plt.savefig(outfile, dpi=800)
        plt.close()
        return

    # Arc Fit
    ax_fit = plt.subplot(gs[0, 1])
    # Points
    ax_fit.scatter(fit['xfit']*fit['xnorm'], fit['yfit'], marker='x')
    if len(fit['xrej']) > 0:
        ax_fit.scatter(fit['xrej']*fit['xnorm'], fit['yrej'], marker='o',
            edgecolor='gray', facecolor='none')
    # Solution
    xval = np.arange(len(arc_spec))
    wave = arutils.func_val(fit['fitc'], xval/fit['xnorm'], 'legendre',
        minv=fit['fmin'], maxv=fit['fmax'])
    ax_fit.plot(xval, wave, 'r-')
    xmin, xmax = 0., len(arc_spec)
    ax_fit.set_xlim(xmin, xmax)
    ymin,ymax = np.min(wave)*.95,  np.max(wave)*1.05
    ax_fit.set_ylim(np.min(wave)*.95,  np.max(wave)*1.05)
    ax_fit.set_ylabel('Wavelength')
    ax_fit.get_xaxis().set_ticks([]) # Suppress labeling
    # Stats
    wave_fit = arutils.func_val(fit['fitc'], fit['xfit'], 'legendre',
        minv=fit['fmin'], maxv=fit['fmax'])
    rms = np.sqrt(np.sum((fit['yfit']-wave_fit)**2)/len(fit['xfit'])) # Ang
    dwv_pix = np.median(np.abs(wave-np.roll(wave,1)))
    ax_fit.text(0.1*len(arc_spec), 0.90*ymin+(ymax-ymin),
        r'$\Delta\lambda$={:.3f}$\AA$ (per pix)'.format(dwv_pix), size='small')
    ax_fit.text(0.1*len(arc_spec), 0.80*ymin+(ymax-ymin),
        'RMS={:.3f} (pixels)'.format(rms/dwv_pix), size='small')
    # Arc Residuals
    ax_res = plt.subplot(gs[1,1])
    res = fit['yfit']-wave_fit
    ax_res.scatter(fit['xfit']*fit['xnorm'], res/dwv_pix, marker='x')
    ax_res.plot([xmin,xmax], [0.,0], 'k--')
    ax_res.set_xlim(xmin, xmax)
    ax_res.set_xlabel('Pixel')
    ax_res.set_ylabel('Residuals (Pix)')

    # Finish
    plt.tight_layout(pad=0.2, h_pad=0.0, w_pad=0.0)
    if outfile == 'show':
        plt.show()
    else:
        plt.savefig(outfile, dpi=800)
    plt.close()

    plt.rcdefaults()


    return



