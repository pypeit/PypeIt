from __future__ import (print_function, absolute_import, division, unicode_literals)

import numpy as np
from pypit import arpca
from pypit import arparse as settings
from pypit import armsgs
from pypit import arsave
from pypit import arutils
from pypit import ararclines
from pypit import arqa
from matplotlib import pyplot as plt
import os

from pypit import ardebug as debugger

# Logging
msgs = armsgs.get_logger()


def detect_lines(slf, det, msarc, censpec=None, MK_SATMASK=False):
    """
    Extract an arc down the center of the chip and identify
    statistically significant lines for analysis.

    Parameters
    ----------
    slf : Class instance
      An instance of the Science Exposure class
    det : int
      Index of the detector
    msarc : ndarray
      Calibration frame that will be used to identify slit traces (in most cases, the slit edge)
    censpec : ndarray, optional
      A 1D spectrum to be searched for significant detections
    MK_SATMASK : bool, optional
      Generate a mask of arc line saturation streaks? Mostly used for echelle data
      when saturation in one order can cause bleeding into a neighbouring order.

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
    satsnd : ndarray
      A mask indicating where which pixels contain saturation streaks
    yprep : ndarray
      The spectrum used to find detections. This spectrum has
      had any "continuum" emission subtracted off
    """
    from pypit import arcyarc
    # Extract a rough spectrum of the arc in each order
    msgs.info("Detecting lines")
    msgs.info("Extracting an approximate arc spectrum at the centre of the chip")
    if msgs._debug['flexure'] or (msarc is None):
        ordcen = slf._pixcen
    else:
        ordcen = slf.GetFrame(slf._pixcen, det)
    if censpec is None:
        #pixcen = np.arange(msarc.shape[0], dtype=np.int)
        #ordcen = (msarc.shape[1]/2)*np.ones(msarc.shape[0],dtype=np.int)
        #if len(ordcen.shape) != 1: msgs.error("The function artrace.model_tilt should only be used for"+msgs.newline()+"a single spectrum (or order)")
        #ordcen = ordcen.reshape((ordcen.shape[0],1))
        msgs.work("No orders being masked at the moment")
        # Average over several pixels to remove some random fluctuations, and increase S/N
        op1 = ordcen+1
        op2 = ordcen+2
        om1 = ordcen-1
        om2 = ordcen-2
        censpec = (msarc[:,ordcen]+msarc[:,op1]+msarc[:,op2]+msarc[:,om1]+msarc[:,om2])/5.0
    # Generate a saturation mask
    if MK_SATMASK:
        ordwid = 0.5*np.abs(slf._lordloc[det-1] - slf._rordloc[det-1])
        msgs.info("Generating a mask of arc line saturation streaks")
        satmask = arcyarc.saturation_mask(msarc, slf._nonlinear[det-1])
        satsnd = arcyarc.order_saturation(satmask, ordcen, (ordwid+0.5).astype(np.int))
    else:
        satsnd = np.zeros_like(ordcen)
    # Detect the location of the arc lines
    msgs.info("Detecting the strongest, nonsaturated lines")
    #####
    # Old algorithm for arc line detection
#   arcdet = arcyarc.detections_allorders(censpec, satsnd)
    #####
    # New algorithm for arc line detection
    #pixels=[]
    siglev = settings.argflag['arc']['calibrate']['detection']
    bpfit = 5  # order of the polynomial used to fit the background 'continuum'
    fitp = settings.argflag['arc']['calibrate']['nfitpix']
    if len(censpec.shape) == 3: detns = censpec[:, 0].flatten()
    else: detns = censpec.copy()
    xrng = np.arange(float(detns.size))
    yrng = np.zeros(detns.size)
    mask = np.zeros(detns.size, dtype=np.int)
    mask[np.where(detns < 0.0)] = 1.0
    mskcnt = np.sum(mask)
    while True:
        w = np.where(mask == 0)
        xfit = xrng[w]
        yfit = detns[w]
        ct = np.polyfit(xfit, yfit, bpfit)
        yrng = np.polyval(ct, xrng)
        sigmed = 1.4826*np.median(np.abs(detns[w]-yrng[w]))
        w = np.where(detns > yrng + siglev*sigmed)
        mask[w] = 1
        if mskcnt == np.sum(mask):
            break  # No new values have been included in the mask
        mskcnt = np.sum(mask)
    w = np.where(mask == 0)
    xfit = xrng[w]
    yprep = detns - yrng
    sfit = 1.4826*np.abs(detns[w]-yrng[w])
    ct = np.polyfit(xfit, sfit, bpfit)
    yerr = np.polyval(ct, xrng)
    myerr = np.median(np.sort(yerr)[:yerr.size//2])
    yerr[np.where(yerr < myerr)] = myerr
    # Find all significant detections
    # The last argument is the overall minimum significance level of an arc line detection and the second
    # last argument is the level required by an individual pixel before the neighbourhood of this pixel is searched.
    pixt = np.where((yprep/yerr > 0.0) &  # (yprep < slf._nonlinear[det-1]) &
                    (yprep > np.roll(yprep, 1)) & (yprep >= np.roll(yprep, -1)) &
                    (np.roll(yprep, 1) > np.roll(yprep, 2)) & (np.roll(yprep, -1) > np.roll(yprep, -2)) &#)[0]
                    (np.roll(yprep, 2) > np.roll(yprep, 3)) & (np.roll(yprep, -2) > np.roll(yprep, -3)))[0]
#                    (np.roll(yprep, 3) > np.roll(yprep, 4)) & (np.roll(yprep, -3) > np.roll(yprep, -4)) & # )[0]
#                    (np.roll(yprep, 4) > np.roll(yprep, 5)) & (np.roll(yprep, -4) > np.roll(yprep, -5)))[0]
    tampl, tcent, twid, ngood = arcyarc.fit_arcorder(xrng, yprep, pixt, fitp)
    w = np.where((~np.isnan(twid)) & (twid > 0.0) & (twid < 10.0/2.35) & (tcent > 0.0) & (tcent < xrng[-1]))
    # Check the results
    #plt.clf()
    #plt.plot(xrng,yprep,'k-')
    #plt.plot(tcent,tampl,'ro')
    #plt.show()
    # Return
    return tampl, tcent, twid, w, satsnd, yprep


def setup_param(slf, sc, det, fitsdict):
    """ Setup for arc analysis

    Parameters
    ----------
    det : int
      detctor index
    fitsdict : dict
      Contains relevant information from fits header files

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

    modify_dict = None
    # Instrument/disperser specific
    sname = settings.argflag['run']['spectrograph']
    idx = settings.spect['arc']['index'][sc]
    disperser = fitsdict["dispname"][idx[0]]
    binspatial, binspectral = settings.parse_binning(fitsdict['binning'][idx[0]])
    if sname == 'shane_kast_blue':
        # Could have the following depend on lamps that were turned on
        lamps = ['CdI','HgI','HeI']
        #arcparam['llist'] = settings.argflag['run']['pypitdir'] + 'data/arc_lines/kast_blue.lst'
        if disperser == '600/4310':
            arcparam['disp']=1.02
            arcparam['b1']=6.88935788e-04
            arcparam['b2']=-2.38634231e-08
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4250.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname=='shane_kast_red':
        lamps = ['NeI']
        #arcparam['llist'] = settings.argflag['run']['pypitdir'] + 'data/arc_lines/kast_red.lst'
        if disperser == '600/7500':
            arcparam['disp']=1.30
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=0.63
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
            arcparam['wv_cen'] = 6600.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname=='shane_kast_red_ret':
        lamps = ['NeI']
        #arcparam['llist'] = settings.argflag['run']['pypitdir'] + 'data/arc_lines/kast_red.lst'
        if disperser == '600/7500':
            arcparam['disp']=2.35
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        elif disperser == '1200/5000':
            arcparam['disp']=1.17
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][0] = 5000.
            arcparam['n_first']=2 # Should be able to lock on
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname=='keck_lris_blue':
        lamps = ['NeI', 'ArI', 'CdI', 'KrI', 'XeI', 'ZnI','CdI','HgI']
        if disperser == '600/4000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.63 # Ang per pixel (unbinned)
            arcparam['b1']= 4.54698031e-04
            arcparam['b2']= -6.86414978e-09
            arcparam['wvmnx'][1] = 6000.
            arcparam['wv_cen'] = 4000.
        elif disperser == '400/3400':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.02
            #arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0]
            arcparam['b1']= 2.72694493e-04
            arcparam['b2']= -5.30717321e-09
            arcparam['wvmnx'][1] = 6000.
        elif disperser == '300/5000':
            arcparam['n_first'] = 2
            arcparam['wv_cen'] = 4500.
            arcparam['disp'] = 1.43
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname=='keck_lris_red':
        arcparam['wv_cen'] = fitsdict['headers'][idx[0]][0]['WAVELEN']
        lamps = ['ArI','NeI','HgI','KrI','XeI']  # Should set according to the lamps that were on
        if disperser == '600/7500':
            arcparam['n_first']=3 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
        elif disperser == '600/10000':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.80 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][1] = 12000.
        elif disperser == '400/8500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=1.19 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][1] = 11000.
            arcparam['min_ampl'] = 3000.  # Lines tend to be very strong
            arcparam['nsig_rej_final'] = 5.
        elif disperser == '900/5500':
            arcparam['n_first']=2 # Too much curvature for 1st order
            arcparam['disp']=0.53 # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0] / binspectral
            arcparam['wvmnx'][1] = 7000.
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname=='wht_isis_blue':
        modify_dict = dict(NeI={'min_wave': 3000.,'min_intensity': 299,
                                'min_Aki': 0.},ArI={'min_intensity': 399.})
        lamps=['CuI','NeI','ArI']
        if fitsdict["dichroic"][idx[0]].strip() == '5300':
            arcparam['wvmnx'][1] = 6000.
        else:
            msgs.error('Not ready for this dichroic {:s}!'.format(disperser))
        if disperser == 'R300B':
            arcparam['n_first']=1  #
            arcparam['disp']=0.80  # Ang per pixel (unbinned)
            arcparam['b1']= 1./arcparam['disp']/slf._msarc[det-1].shape[0]
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    elif sname == 'tng_dolores':
        lamps = ['NeI', 'HgI']
        if disperser == 'LR-R':
            arcparam['n_first'] = 2  # Too much curvature for 1st order
            arcparam['disp'] = 2.61  # Ang per pixel (unbinned)
            arcparam['disp_toler'] = 0.1  # Ang per pixel (unbinned)
            arcparam['wvmnx'][0] = 4470.0
            arcparam['wvmnx'][1] = 10073.0
            arcparam['wv_cen'] = 7400.
            arcparam['b1'] = 1. / arcparam['disp'] / slf._msarc[det - 1].shape[0] / binspectral
        else:
            msgs.error('Not ready for this disperser {:s}!'.format(disperser))
    else:
        msgs.error('ararc.setup_param: Not ready for this instrument {:s}!'.format(sname))
    # Load linelist
    if settings.argflag['arc']['calibrate']['lamps'] is not None:
        arcparam['lamps'] = settings.argflag['arc']['calibrate']['lamps']
    else:
        arcparam['lamps'] = lamps
    slmps = lamps[0]
    for lamp in lamps[1:]:
        slmps=slmps+','+lamp
    msgs.info('Loading line list using {:s} lamps'.format(slmps))
    arcparam['llist'] = ararclines.load_arcline_list(slf, idx, lamps, disperser,
        wvmnx=arcparam['wvmnx'], modify_parse_dict=modify_dict)
    # Binning
    arcparam['disp'] *= binspectral

    # Return
    return arcparam


def simple_calib(slf, det, get_poly=False):
    """Simple calibration algorithm for longslit wavelengths

    Uses slf._arcparam to guide the analysis

    Parameters
    ----------
    get_poly : bool, optional
      Pause to record the polynomial pix = b0 + b1*lambda + b2*lambda**2

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """

    # Extract the arc
    msgs.work("Detecting lines..")
    tampl, tcent, twid, w, satsnd, yprep = detect_lines(slf, det, slf._msarc[det-1])

    # Cut down to the good ones
    tcent = tcent[w]
    tampl = tampl[w]
    msgs.info('Detected {:d} lines in the arc spectrum.'.format(len(w[0])))

    # Parameters (just for convenience)
    aparm = slf._arcparam[det-1]

    # Read Arc linelist
    llist = aparm['llist']

    # IDs were input by hand
    if len(settings.argflag['arc']['calibrate']['IDpixels']) > 0:
        # Check that there are at least 5 values
        pixels = np.array(settings.argflag['arc']['calibrate']['IDpixels'])
        if np.sum(pixels > 0.) < 5:
            msgs.error("Need to give at least 5 pixel values!")
        #
        msgs.info("Using input lines to seed the wavelength solution")
        # Calculate median offset
        mdiff = [np.min(np.abs(tcent-pix)) for pix in
                 settings.argflag['arc']['calibrate']['IDpixels']]
        med_poff = np.median(np.array(mdiff))
        msgs.info("Will apply a median offset of {:g} pixels".format(med_poff))

        # Match input lines to observed spectrum
        nid = len(settings.argflag['arc']['calibrate']['IDpixels'])
        idx_str = np.ones(nid).astype(int)
        ids = np.zeros(nid)
        idsion = np.array(['     ']*nid)
        gd_str = np.arange(nid).astype(int)
        for jj,pix in enumerate(settings.argflag['arc']['calibrate']['IDpixels']):
            diff = np.abs(tcent-pix-med_poff)
            if np.min(diff) > 2.:
                debugger.set_trace()
                msgs.error("No match with input pixel {:g}!".format(pix))
            else:
                imn = np.argmin(diff)
            # Set
            idx_str[jj] = imn
            # Take wavelength from linelist instead of input value
            wdiff = np.abs(llist['wave']-settings.argflag['arc']['calibrate']['IDwaves'][jj])
            imnw = np.argmin(wdiff)
            if wdiff[imnw] > 0.015:  # Arbitrary tolerance
                msgs.error("Input IDwaves={:g} is not in the linelist.  Fix".format(
                        settings.argflag['arc']['calibrate']['IDwaves'][jj]))
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
            dpix_list[kk,:] = slf._msarc[det-1].shape[0]*(aparm['b1']*(np.array(row['wave'] - llist['wave'])) + aparm['b2']*np.array(row['wave']**2 - llist['wave']**2) )

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
    xfit, yfit = tcent[ifit]/(slf._msarc[det-1].shape[0]-1), all_ids[ifit]
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
        msarc = slf._msarc[det-1]
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
        ions=ions, fmin=fmin, fmax=fmax, xnorm=float(slf._msarc[det-1].shape[0]),
        xrej=xrej, yrej=yrej, mask=mask, spec=yprep, nrej=aparm['nsig_rej_final'],
        shift=0., tcent=tcent)
    # QA
    arqa.arc_fit_qa(slf, final_fit)
    # RMS
    rms_ang = arutils.calc_fit_rms(xfit, yfit, fit, aparm['func'], minv=fmin, maxv=fmax)
    wave = arutils.func_val(fit, np.arange(slf._msarc[det-1].shape[0])/float(slf._msarc[det-1].shape[0]),
                            aparm['func'], minv=fmin, maxv=fmax)
    rms_pix = rms_ang/np.median(np.abs(wave-np.roll(wave,1)))
    msgs.info("Fit RMS = {} pix".format(rms_pix))
    # Return
    return final_fit


def calib_with_arclines(slf, det, get_poly=False, use_method="general"):
    """Simple calibration algorithm for longslit wavelengths

    Uses slf._arcparam to guide the analysis

    Parameters
    ----------
    get_poly : bool, optional
      Pause to record the polynomial pix = b0 + b1*lambda + b2*lambda**2

    Returns
    -------
    final_fit : dict
      Dict of fit info
    """
    from arclines.holy.grail import basic, semi_brute, general
    # Parameters (just for convenience)
    aparm = slf._arcparam[det-1]
    # Extract the arc
    msgs.work("Detecting lines")
    tampl, tcent, twid, w, satsnd, spec = detect_lines(slf, det, slf._msarc[det-1])

    if use_method == "semi-brute":
        best_dict, final_fit = semi_brute(spec, aparm['lamps'], aparm['wv_cen'], aparm['disp'], fit_parm=aparm, min_ampl=aparm['min_ampl'])
    elif use_method == "basic":
        stuff = basic(spec, aparm['lamps'], aparm['wv_cen'], aparm['disp'])
        status, ngd_match, match_idx, scores, final_fit = stuff
    else:
        # Now preferred
        best_dict, final_fit = general(spec, aparm['lamps'], fit_parm=aparm, min_ampl=aparm['min_ampl'])
    arqa.arc_fit_qa(slf, final_fit)
    #
    return final_fit
