""" Module for PYPIT extraction code
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import time
import copy
import inspect

import numpy as np

from matplotlib import pyplot as plt
from matplotlib import gridspec, font_manager

from astropy import units
from astropy.stats import sigma_clip

from pypit import msgs
from pypit import arqa
from pypit import artrace
from pypit import arutils
from pypit import ardebug as debugger

# MASK VALUES FROM EXTRACTION
# 0 
# 2**0 = Flagged as bad detector pixel
# 2**1 = Flagged as affected by Cosmic Ray 
# 2**5 = Flagged as NAN (from something gone wrong)
# 2**6 = Entire region masked

mask_flags = dict(bad_pix=2**0, CR=2**1, NAN=2**5, bad_row=2**6)


# This routine is deprecated, replaced by extract_boxcar
def boxcar(specobjs, sciframe, varframe, bpix, skyframe, crmask, scitrace, mswave,
           maskslits, slitpix):
    """ Perform boxcar extraction on the traced objects.
    Also perform a local sky subtraction

    Parameters
    ----------
    specobjs : list of dict
      list of SpecObj objects
    sciframe : ndarray
      science frame, should be sky subtracted
    varframe : ndarray
      variance image
    bgframe : ndarray
      sky background frame
    crmask : int ndarray
        mask of cosmic ray hits
    scitrace : list
      traces, object and background trace images

    Returns
    -------
    bgcorr : ndarray
      Correction to the sky background in the object window
    """
    # Subtract the global sky
    skysub = sciframe-skyframe

    bgfitord = 1  # Polynomial order used to fit the background
    slits = range(len(scitrace))
    nslit = len(slits)
    gdslits = np.where(~maskslits)[0]
    cr_mask = 1.0-crmask
    bgfit = np.linspace(0.0, 1.0, sciframe.shape[1])
    bgcorr = np.zeros_like(cr_mask)
    # Loop on Slits
    for sl in slits:
        if sl not in gdslits:
            continue
        word = np.where((slitpix == sl + 1) & (varframe > 0.))
        if word[0].size == 0:
            continue
        mask_slit = np.zeros(sciframe.shape, dtype=np.float)
        mask_slit[word] = 1.0
        # Loop on Objects
        nobj = scitrace[sl]['nobj']
        for o in range(nobj):
            msgs.info("Performing boxcar extraction of object {0:d}/{1:d} in slit {2:d}/{3:d}".format(o+1, nobj, sl+1, nslit))
            if scitrace[sl]['object'] is None:
                # The object for all slits is provided in the first extension
                objreg = np.copy(scitrace[0]['object'][:, :, o])
                wzro = np.where(slitpix != sl + 1)
                objreg[wzro] = 0.0
            else:
                objreg = scitrace[sl]['object'][:, :, o]
            # Fit the background
            msgs.info("   Fitting the background")
            if scitrace[sl]['background'] is None:
                # The background for all slits is provided in the first extension
                bckreg = np.copy(scitrace[0]['background'][:, :, o])
                wzro = np.where(slitpix != sl + 1)
                bckreg[wzro] = 0.0
            else:
                bckreg = scitrace[sl]['background'][:, :, o]
            # Trim CRs further
            bg_mask = np.zeros_like(sciframe)
            bg_mask[np.where((bckreg*cr_mask <= 0.))] = 1.
            bg_mask[np.where((slitpix != sl + 1))] = 1.
            mask_sci = np.ma.array(skysub, mask=bg_mask, fill_value=0.)
            clip_image = sigma_clip(mask_sci, axis=1, sigma=3.)  # For the mask only
            # Fit
#            print('calling func2d_fit_val')
#            t = time.clock()
#            _bgframe = arcyutils.func2d_fit_val(bgfit, sciframe,
#                                                (~clip_image.mask)*bckreg*cr_mask, bgfitord)
#            print('Old func2d_fit_val: {0} seconds'.format(time.clock() - t))
#            t = time.clock()
            bgframe = new_func2d_fit_val(skysub, bgfitord, x=bgfit,
                                         w=(~clip_image.mask)*bckreg*cr_mask)
#            print('New func2d_fit_val: {0} seconds'.format(time.clock() - t))
            # Some fits are really wonky ... in both methods
#            if np.sum(bgframe != _bgframe) != 0:
#                plt.imshow(np.ma.log10(sciframe), origin='lower', interpolation='nearest',
#                           aspect='auto')
#                plt.colorbar()
#                plt.show()
#                w=(~clip_image.mask)*bckreg*cr_mask
#                plt.imshow(np.ma.log10(sciframe*w), origin='lower', interpolation='nearest',
#                           aspect='auto')
#                plt.colorbar()
#                plt.show()
#                plt.imshow(np.ma.log10(_bgframe), origin='lower',
#                           interpolation='nearest', aspect='auto')
#                plt.colorbar()
#                plt.show()
#                plt.imshow(np.ma.log10(bgframe), origin='lower',
#                           interpolation='nearest', aspect='auto')
#                plt.colorbar()
#                plt.show()
#                plt.imshow(np.ma.log10(np.absolute(bgframe-_bgframe)), origin='lower',
#                           interpolation='nearest', aspect='auto')
#                plt.colorbar()
#                plt.show()
#                plt.imshow(np.ma.log10(np.ma.divide(bgframe,_bgframe)), origin='lower',
#                           interpolation='nearest', aspect='auto')
#                plt.colorbar()
#                plt.show()
#
#                d = np.amax(np.absolute(bgframe-_bgframe), axis=1)
#                i = np.argmax(d)
#                plt.plot(bgfit, sciframe[i,:])
#                plt.plot(bgfit, sciframe[i,:]*w[i,:])
#                plt.plot(bgfit, bgframe[i,:])
#                plt.plot(bgfit, _bgframe[i,:])
#                plt.show()
#
#            assert np.sum(bgframe != _bgframe) == 0, 'Difference between old and new func2d_fit_val'

            # Weights
            weight = objreg*mask_slit
            sumweight = np.sum(weight, axis=1)
            # Deal with fully masked regions
            fully_masked = sumweight == 0.
            if np.any(fully_masked):
                weight[fully_masked,:] = 1.
                sumweight[fully_masked] = weight.shape[1]
            # Generate wavelength array (average over the pixels)
            wvsum = np.sum(mswave*weight, axis=1)
            wvsum /= sumweight
            # Generate sky spectrum (flux per pixel)
            skysum = np.sum(skyframe*weight, axis=1)
            skysum /= sumweight
            # Total the object flux
            msgs.info("   Summing object counts")
            scisum = np.sum((skysub-bgframe)*weight, axis=1)
            # Total the variance array
            msgs.info("   Summing variance array")
            varsum = np.sum(varframe*weight, axis=1)
            # Update background correction image
            tmp = bckreg + objreg
            gdp = np.where((tmp > 0) & (slitpix == sl + 1))
            bgcorr[gdp] = bgframe[gdp]
            # Mask
            boxmask = np.zeros(wvsum.shape, dtype=np.int)
            # Bad detector pixels
            BPs = np.sum(weight*bpix, axis=1)
            bp = BPs > 0.
            boxmask[bp] += mask_flags['bad_pix']
            # CR
            CRs = np.sum(weight*cr_mask, axis=1)
            cr = CRs > 0.
            boxmask[cr] += mask_flags['CR']
            # Fully masked
            if np.any(fully_masked):
                boxmask[fully_masked] += mask_flags['bad_row']
                scisum[fully_masked] = 0.
                varsum[fully_masked] = 0.
                skysum[fully_masked] = 0.
            # NAN
            NANs = np.isnan(scisum)
            if np.any(NANs):
                msgs.warn("   NANs in the spectrum somehow...")
                boxmask[NANs] += mask_flags['NANs']
                scisum[NANs] = 0.
                varsum[NANs] = 0.
                skysum[NANs] = 0.
            # Check on specobjs
            if not specobjs[sl][o].check_trace(scitrace[sl]['traces'][:, o]):
                debugger.set_trace()
                msgs.error("Bad match to specobj in boxcar!")
            # Fill
            specobjs[sl][o].boxcar['wave'] = wvsum.copy()*units.AA  # Yes, units enter here
            specobjs[sl][o].boxcar['counts'] = scisum.copy()
            specobjs[sl][o].boxcar['var'] = varsum.copy()
            if np.sum(specobjs[sl][o].boxcar['var']) == 0.:
                debugger.set_trace()
            specobjs[sl][o].boxcar['sky'] = skysum.copy()  # per pixel
            specobjs[sl][o].boxcar['mask'] = boxmask.copy()
            # Find boxcar size
            slit_sz = []
            inslit = np.where(weight == 1.)
            for ii in range(weight.shape[0]):
                inrow = inslit[0] == ii
                if np.sum(inrow) > 0:
                    slit_sz.append(np.max(inslit[1][inrow])-np.min(inslit[1][inrow]))
            slit_pix = np.median(slit_sz)  # Pixels
            specobjs[sl][o].boxcar['size'] = slit_pix
    # Return
    return bgcorr

# This routine is deprecated
def obj_profiles(det, specobjs, sciframe, varframe, crmask,
                 scitrace, tilts, maskslits, slitpix,
                 extraction_profile='gaussian',
                 COUNT_LIM=25., doqa=True, pickle_file=None):
    """ Derive spatial profiles for each object
    Parameters
    ----------
    det : int
    specobjs : list
    sciframe : ndarray
      Sky subtracted science frame
    varframe : ndarray
    crmask : ndarray
    scitrace : list
    tilts : ndarray
    maskslits : ndarray

    Returns
    -------
    All internal on specobj and scitrace objects
    """
    # Init QA
    #
    sigframe = np.sqrt(varframe)
    slits = range(len(specobjs))
    gdslits = np.where(~maskslits)[0]
    # Loop on slits
    for sl in slits:
        if sl not in gdslits:
            continue
        # Loop on objects
        nobj = scitrace[sl]['nobj']
        if nobj == 0:
            continue
        scitrace[sl]['opt_profile'] = []
        msgs.work("Should probably loop on S/N")
        for o in range(nobj):
            msgs.info("Deriving spatial profile of object {0:d}/{1:d} in slit {2:d}/{3:d}".format(o+1, nobj, sl+1, len(specobjs)))
            # Get object pixels
            if scitrace[sl]['background'] is None:
                # The object for all slits is provided in the first extension
                objreg = np.copy(scitrace[0]['object'][:, :, o])
                wzro = np.where(slitpix != sl + 1)
                objreg[wzro] = 0.0
            else:
                objreg = scitrace[sl]['object'][:, :, o]
            # Calculate slit image
            slit_img = artrace.slit_image(scitrace[sl], o, tilts)
            # Object pixels
            weight = objreg.copy()
            # Identify good rows
            gdrow = np.where(specobjs[sl][o].boxcar['counts'] > COUNT_LIM)[0]
            # Normalized image
            norm_img = sciframe / np.outer(specobjs[sl][o].boxcar['counts'], np.ones(sciframe.shape[1]))
            # Eliminate rows with CRs (wipes out boxcar)
            crspec = np.sum(crmask*weight, axis=1)
            cr_rows = np.where(crspec > 0)[0]
            weight[cr_rows, :] = 0.
            #
            if len(gdrow) > 100:  # Good S/N regime
                msgs.info("Good S/N for profile")
                # Eliminate low count regions
                badrow = np.where(specobjs[sl][o].boxcar['counts'] < COUNT_LIM)[0]
                weight[badrow, :] = 0.
                # Extract profile
                gdprof = (weight > 0) & (sigframe > 0.) & (~np.isnan(slit_img))  # slit_img=nan if the slit is partially on the chip
                slit_val = slit_img[gdprof]
                flux_val = norm_img[gdprof]
                #weight_val = sciframe[gdprof]/sigframe[gdprof]  # S/N
                weight_val = 1./sigframe[gdprof]  # 1/N
                msgs.work("Weight by S/N in boxcar extraction? [avoid CRs; smooth?]")
                # Fit
                fdict = dict(func=extraction_profile, deg=3, extrap=False)
                if fdict['func'] == 'gaussian':
                    fdict['deg'] = 2
                elif fdict['func'] == 'moffat':
                    fdict['deg'] = 3
                else:
                    msgs.error("Not ready for this type of object profile")
                msgs.work("Might give our own guess here instead of using default")
                guess = None
                # Check if there are enough pixels in the slit to perform fit
                if slit_val.size <= fdict['deg'] + 1:
                    msgs.warn("Not enough pixels to determine profile of object={0:s} in slit {1:d}".format(specobjs[sl][o].idx, sl+1) + msgs.newline() +
                              "Skipping Optimal")
                    fdict['extrap'] = True
                    scitrace[sl]['opt_profile'].append(copy.deepcopy(fdict))
                    continue
                # Fit the profile
                try:
                    mask, gfit = arutils.robust_polyfit(slit_val, flux_val, fdict['deg'], function=fdict['func'], weights=weight_val, maxone=False, guesses=guess)
                except RuntimeError:
                    msgs.warn("Bad Profile fit for object={:s}." + msgs.newline() +
                              "Skipping Optimal".format(specobjs[sl][o].idx))
                    fdict['extrap'] = True
                    scitrace[sl]['opt_profile'].append(copy.deepcopy(fdict))
                    continue
                except ValueError:
                    debugger.set_trace()  # NaNs in the values?  Check
                msgs.work("Consider flagging/removing CRs here")
                # Record
                fdict['param'] = gfit.copy()
                fdict['mask'] = mask
                fdict['slit_val'] = slit_val
                fdict['flux_val'] = flux_val
                scitrace[sl]['opt_profile'].append(copy.deepcopy(fdict))
                specobjs[sl][o].optimal['fwhm'] = fdict['param'][1]  # Pixels
                if msgs._debug['obj_profile']:
                    gdp = mask == 0
                    mn = np.min(slit_val[gdp])
                    mx = np.max(slit_val[gdp])
                    xval = np.linspace(mn, mx, 1000)
                    model = arutils.func_val(gfit, xval, fdict['func'])
                    plt.clf()
                    ax = plt.gca()
                    ax.scatter(slit_val[gdp], flux_val[gdp], marker='.', s=0.7, edgecolor='none', facecolor='black')
                    ax.plot(xval, model, 'b')
                    # Gaussian too?
                    if False:
                        fdictg = dict(func='gaussian', deg=2)
                        maskg, gfitg = arutils.robust_polyfit(slit_val, flux_val, fdict['deg'], function=fdictg['func'], weights=weight_val, maxone=False)
                        modelg = arutils.func_val(gfitg, xval, fdictg['func'])
                        ax.plot(xval, modelg, 'r')
                    plt.show()
                    debugger.set_trace()
            elif len(gdrow) > 10:  #
                msgs.warn("Low extracted flux for obj={:s} in slit {:d}.  Not ready for Optimal".format(specobjs[sl][o].idx,sl+1))
                scitrace[sl]['opt_profile'].append({})
                continue
            elif len(gdrow) >= 0:  # limit is ">= 0" to avoid crash for gdrow=0
                msgs.warn("Low extracted flux for obj={:s} in slit {:d}.  Not ready for Optimal".format(specobjs[sl][o].idx, sl + 1))
                scitrace[sl]['opt_profile'].append({})
                continue
    # QA
    if doqa: #not msgs._debug['no_qa'] and doqa:
        msgs.info("Preparing QA for spatial object profiles")
#        arqa.obj_profile_qa(slf, specobjs, scitrace, det)
        debugger.set_trace()  # Need to avoid slf
        obj_profile_qa(specobjs, scitrace, det)
    return


# This routine is deprecated, replaced by fit_profile_qa
def obj_profile_qa(slf, specobjs, scitrace, det):
    """ Generate a QA plot for the object spatial profile
    Parameters
    ----------
    """

    plt.rcdefaults()
    plt.rcParams['font.family']= 'times new roman'

    method = inspect.stack()[0][3]
    gdslits = np.where(~slf._maskslits[det-1])[0]
    for sl in range(len(specobjs)):
        if sl not in gdslits:
            continue
        # Setup
        nobj = scitrace[sl]['nobj']
        if nobj == 0:
            continue
        ncol = min(3, nobj)
        nrow = nobj // ncol + ((nobj % ncol) > 0)
        # Outfile
        outfile = arqa.set_qa_filename(slf._basename, method, det=det, slit=specobjs[sl][0].slitid)
        # Plot
        plt.figure(figsize=(8, 5.0))
        plt.clf()
        gs = gridspec.GridSpec(nrow, ncol)

        # Plot
        for o in range(nobj):
            fdict = scitrace[sl]['opt_profile'][o]
            if 'param' not in fdict.keys():  # Not optimally extracted
                continue
            ax = plt.subplot(gs[o//ncol, o % ncol])

            # Data
            gdp = fdict['mask'] == 0
            ax.scatter(fdict['slit_val'][gdp], fdict['flux_val'][gdp], marker='.',
                       s=0.5, edgecolor='none')

            # Fit
            mn = np.min(fdict['slit_val'][gdp])
            mx = np.max(fdict['slit_val'][gdp])
            xval = np.linspace(mn, mx, 1000)
            fit = arutils.func_val(fdict['param'], xval, fdict['func'])
            ax.plot(xval, fit, 'r')
            # Axes
            ax.set_xlim(mn,mx)
            # Label
            ax.text(0.02, 0.90, 'Obj={:s}'.format(specobjs[sl][o].idx),
                    transform=ax.transAxes, ha='left', size='small')

        plt.savefig(outfile, dpi=500)
        plt.close()

    plt.rcdefaults()

# This routine is deprecated, replaced by extract_optimal
def optimal_extract(specobjs, sciframe, varframe,
                    crmask, scitrace, tilts, mswave,
                    maskslits, slitpix, calib_wavelength='vacuum',
                    pickle_file=None, profiles=None):
    """ Preform optimal extraction
    Standard Horne approach

    Parameters
    ----------
    slf
    det
    specobjs
    sciframe
    varframe
    crmask
    scitrace
    COUNT_LIM
    pickle_file

    Returns
    -------
    newvar : ndarray
      Updated variance array that includes object model
    """
    # Inverse variance
    model_ivar = np.zeros_like(varframe)
    cr_mask = 1.0-crmask
    gdvar = (varframe > 0.) & (cr_mask == 1.)
    model_ivar[gdvar] = arutils.calc_ivar(varframe[gdvar])
    # Object model image
    obj_model = np.zeros_like(varframe)
    gdslits = np.where(~maskslits)[0]
    # Loop on slits
    for sl in range(len(specobjs)):
        if sl not in gdslits:
            continue
        # Loop on objects
        nobj = scitrace[sl]['nobj']
        if nobj == 0:
            continue
        for o in range(nobj):
            msgs.info("Performing optimal extraction of object {0:d}/{1:d} in slit {2:d}/{3:d}".format(o+1, nobj, sl+1, len(specobjs)))
            # Get object pixels
            if scitrace[sl]['background'] is None:
                # The object for all slits is provided in the first extension
                objreg = np.copy(scitrace[0]['object'][:, :, o])
                wzro = np.where(slitpix != sl + 1)
                objreg[wzro] = 0.0
            else:
                objreg = scitrace[sl]['object'][:, :, o]
            # Fit dict
            fit_dict = scitrace[sl]['opt_profile'][o]
            if 'param' not in fit_dict.keys():
                continue
            # Slit image
            slit_img = artrace.slit_image(scitrace[sl], o, tilts)
            #msgs.warn("Turn off tilts")
            # Object pixels
            weight = objreg.copy()
            gdo = (weight > 0) & (model_ivar > 0)
            # Profile image
            prof_img = np.zeros_like(weight)
            prof_img[gdo] = arutils.func_val(fit_dict['param'], slit_img[gdo],
                                             fit_dict['func'])
            # Normalize
            norm_prof = np.sum(prof_img, axis=1)
            prof_img /= np.outer(norm_prof + (norm_prof == 0.), np.ones(prof_img.shape[1]))
            # Mask (1=good)
            mask = np.zeros_like(prof_img)
            mask[gdo] = 1.
            mask *= cr_mask

            # Optimal flux
            opt_num = np.sum(mask * sciframe * model_ivar * prof_img, axis=1)
            opt_den = np.sum(mask * model_ivar * prof_img**2, axis=1)
            opt_flux = opt_num / (opt_den + (opt_den == 0.))
            # Optimal wave
            opt_num = np.sum(mswave * model_ivar * prof_img**2, axis=1)
            opt_den = np.sum(model_ivar * prof_img**2, axis=1)
            opt_wave = opt_num / (opt_den + (opt_den == 0.))
            # Replace fully masked rows with mean wavelength (they will have zero flux and ivar)
            full_mask = opt_num == 0.
            if np.any(full_mask):
                msgs.warn("Replacing fully masked regions with mean wavelengths")
                mnwv = np.mean(mswave, axis=1)
                opt_wave[full_mask] = mnwv[full_mask]
            if (np.sum(opt_wave < 1.) > 0) and calib_wavelength != "pixel":
                debugger.set_trace()
                msgs.error("Zero value in wavelength array. Uh-oh")
            # Optimal ivar
            opt_num = np.sum(mask * model_ivar * prof_img**2, axis=1)
            ivar_den = np.sum(mask * prof_img, axis=1)
            opt_ivar = opt_num * arutils.calc_ivar(ivar_den)

            # Save
            specobjs[sl][o].optimal['wave'] = opt_wave.copy()*units.AA  # Yes, units enter here
            specobjs[sl][o].optimal['counts'] = opt_flux.copy()
            gdiv = (opt_ivar > 0.) & (ivar_den > 0.)
            opt_var = np.zeros_like(opt_ivar)
            opt_var[gdiv] = arutils.calc_ivar(opt_ivar[gdiv])
            specobjs[sl][o].optimal['var'] = opt_var.copy()
            #specobjs[o].boxcar['sky'] = skysum  # per pixel

            # Update object model
            counts_image = np.outer(opt_flux, np.ones(prof_img.shape[1]))
            obj_model += prof_img * counts_image
            '''
            if 'OPTIMAL' in msgs._debug:
                debugger.set_trace()
                debugger.xplot(opt_wave, opt_flux, np.sqrt(opt_var))
            '''

    # KBW: Using variance_frame here produces a circular import.  I've
    # changed this function to return the object model, then this last
    # step is done in arproc.

#    # Generate new variance image
#    newvar = arproc.variance_frame(slf, det, sciframe, -1,
#                                   skyframe=skyframe, objframe=obj_model)
#    # Return
#    return newvar

    return obj_model


def boxcar_cen(slf, det, img):
    """ Simple boxcar down center of the slit

    Parameters
    ----------
    slf
    det
    img

    Returns
    -------
    spec : ndarray

    """
    # Extract a slit down the center (as in ararc, or should be!)
    ordcen = slf.GetFrame(slf._pixcen, det)
    op1 = ordcen+1
    op2 = ordcen+2
    om1 = ordcen-1
    om2 = ordcen-2
    # Extract
    censpec = (img[:,ordcen]+img[:,op1]+img[:,op2]+img[:,om1]+img[:,om2])/5.0
    if len(censpec.shape) == 3:
        censpec = censpec[:, 0].flatten()
    # Return
    return censpec


def new_func2d_fit_val(y, order, x=None, w=None):
    """
    Fit a polynomial to each column in y.

    if y is 2D, always fit along columns

    test if y is a MaskedArray
    """
    # Check input
    if y.ndim > 2:
        msgs.error('y cannot have more than 2 dimensions.')

    _y = np.atleast_2d(y)
    ny, npix = _y.shape

    # Set the x coordinates
    if x is None:
        _x = np.linspace(-1,1,npix)
    else:
        if x.ndim != 1:
            msgs.error('x must be a vector')
        if x.size != npix:
            msgs.error('Input x must match y vector or column length.')
        _x = x.copy()

    # Generate the Vandermonde matrix
    vand = np.polynomial.polynomial.polyvander(_x, order)

    # Fit with appropriate weighting
    if w is None:
        # Fit without weights
        c = np.linalg.lstsq(vand, _y.T)[0]
        ym = np.sum(c[:,:,None] * vand.T[:,None,:], axis=0)
    elif w.ndim == 1:
        # Fit with the same weight for each vector
        _vand = w[:,None] * vand
        __y = w[None,:] * _y
        c = np.linalg.lstsq(_vand, __y.T)[0]
        ym = np.sum(c[:,:,None] * vand.T[:,None,:], axis=0)
    else:
        # Fit with different weights for each vector
        if w.shape != y.shape:
            msgs.error('Input w must match y axis length or y shape.')
        # Prep the output model
        ym = np.empty(_y.shape, dtype=float)
        # Weight the data
        __y = w * _y
        for i in range(ny):
            # Weight the Vandermonde matrix for this y vector
            _vand = w[i,:,None] * vand
            c = np.linalg.lstsq(_vand, __y[i,:])[0]
            ym[i,:] = np.sum(c[:,None] * vand.T[:,:], axis=0)

    # Return the model with the appropriate shape
    return ym if y.ndim == 2 else ym[0,:]



