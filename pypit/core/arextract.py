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






def extract_asymbox2(image,left_in,right_in,ycen = None,weight_image = None):
    """ Extract the total flux within a variable window at many positions. This routine will accept an asymmetric/variable window
    specified by the left_in and right_in traces.  The ycen position is optional. If it is not provied, it is assumed to be integers
    in the spectral direction (as is typical for traces). Traces are expected to run vertically to be consistent with other
    extract_  routines. Based on idlspec2d/spec2d/extract_asymbox2.pro

    Parameters
    ----------
    image :  float ndarray
        Image to extract from. It is a 2-d array with shape (nspec, nspat)
    left  :  float ndarray
        Left boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

    right  :  float ndarray
        Right boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.


    Optional Parameters
    -------------------
    ycen :  float ndarray
        Y positions corresponding to "Left"  and "Right" (expected as integers). Will be cast to an integer if floats
        are provided. This needs to have the same shape as left and right broundarys provided above. In other words,
        either a  2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

    weight_image: float ndarray
        Weight map to be applied to image before boxcar. It is a 2-d array with shape (nspec, nspat)

    Returns
    -------
    fextract:   ndarray
       Extracted flux at positions specified by (left<-->right, ycen). The output will have the same shape as
       Left and Right, i.e.  an 2-d  array with shape (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for
       the case of a single trace.


    Revision History
    ----------------
    24-Mar-1999  Written by David Schlegel, Princeton.
    17-Feb-2003  Written with slow IDL routine, S. Burles, MIT
    22-Apr-2018  Ported to python by Joe Hennawi
    """

    # ToDO it would be nice to avoid this transposing, but I got confused during the IDL port
    left = left_in.T
    right = right_in.T

    dim = left.shape
    ndim = left.ndim
    if (ndim == 1):
        nTrace = 1
        npix = dim[0]
    else:
        nTrace = dim[0]
        npix = dim[1]

    if ycen is None:
        if ndim == 1:
            ycen_out = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen_out = np.outer(np.ones(nTrace, dtype='int'), np.arange(npix, dtype='int'), )
        else:
            raise ValueError('trace is not 1 or 2 dimensional')
    else:
        ycen_out = ycen.T
        ycen_out = np.rint(ycen_out).astype(int)

    if ((np.size(left) != np.size(ycen_out)) | (np.shape(left) != np.shape(ycen_out))):
        raise ValueError('Number of elements and left of trace and ycen must be equal')

    idims = image.shape
    nspat = idims[1]
    nspec = idims[0]

    maxwindow = np.max(right - left)
    tempx = np.int(maxwindow + 3.0)

    bigleft = np.outer(left[:], np.ones(tempx))
    bigright = np.outer(right[:], np.ones(tempx))
    spot = np.outer(np.ones(npix * nTrace), np.arange(tempx)) + bigleft - 1
    bigy = np.outer(ycen_out[:], np.ones(tempx, dtype='int'))

    fullspot = np.array(np.fmin(np.fmax(np.round(spot + 1) - 1, 0), nspat - 1), int)
    fracleft = np.fmax(np.fmin(fullspot - bigleft, 0.5), -0.5)
    fracright = np.fmax(np.fmin(bigright - fullspot, 0.5), -0.5)
    del bigleft
    del bigright
    bool_mask1 = (spot >= -0.5) & (spot < (nspat - 0.5))
    bool_mask2 = (bigy >= 0) & (bigy <= (nspec - 1))
    weight = (np.fmin(np.fmax(fracleft + fracright, 0), 1)) * bool_mask1 * bool_mask2
    del spot
    del fracleft
    del fracright
    bigy = np.fmin(np.fmax(bigy, 0), nspec - 1)

    if weight_image != None:
        temp = np.array([weight_image[x1, y1] * image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)
        temp_wi = np.array([weight_image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2_wi = np.reshape(weight.flatten() * temp_wi, (nTrace, npix, tempx))
        f_ivar = np.sum(temp2_wi, axis=2)
        fextract = fextract / (f_ivar + (f_ivar == 0)) * (f_ivar > 0)
    else:
        # Might be more pythonic way to code this. I needed to switch the flattening order in order to get
        # this to work
        temp = np.array([image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)

    # IDL version model functionality not implemented yet
    # At the moment I'm not reutnring the f_ivar for the weight_image mode. I'm not sure that this functionality is even
    # ever used

    if(nTrace ==1):
        fextract = fextract.reshape(npix)
    return fextract.T


def extract_boxcar(image,trace_in, radius_in, ycen = None):
    """ Extract the total flux within a boxcar window at many positions. The ycen position is optional. If it is not provied, it is assumed to be integers
     in the spectral direction (as is typical for traces). Traces are expected to run vertically to be consistent with other
     extract_  routines. Based on idlspec2d/spec2d/extract_boxcar.pro

     Parameters
     ----------
     image :  float ndarray
         Image to extract from. It is a 2-d array with shape (nspec, nspat)

     trace_in :  float ndarray
         Trace for the region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
         (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

     radius :  float or ndarray
         boxcar radius in floating point pixels. This can be either be in put as a scalar or as an array to perform
         boxcar extraction a varaible radius. If an array is input it must have the same size and shape as trace_in, i.e.
         a 2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.


     Optional Parameters
     -------------------
     ycen :  float ndarray
         Y positions corresponding to trace_in (expected as integers). Will be rounded to the nearest integer if floats
         are provided. This needs to have the same shape as trace_in  provided above. In other words,
         either a  2-d  array with shape (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.


     Returns
     -------
     fextract:   ndarray
         Extracted flux at positions specified by (left<-->right, ycen). The output will have the same shape as
         Left and Right, i.e.  an 2-d  array with shape (nspec, nTrace) array if multiple traces were input, or a 1-d array with shape (nspec) for
         the case of a single trace.

     Revision History
     ----------------
     24-Mar-1999  Written by David Schlegel, Princeton.
     22-Apr-2018  Ported to python by Joe Hennawi
     """


    # Checks on radius
    if (isinstance(radius_in,int) or isinstance(radius_in,float)):
        radius = radius_in
    elif ((np.size(radius)==np.size(trace_in)) & (np.shape(radius) == np.shape(trace_in))):
        radius = radius_in.T
    else:
        raise ValueError('Boxcar radius must a be either an integer, a floating point number, or an ndarray '
                         'with the same shape and size as trace_in')

    trace = trace_in.T

    dim = trace.shape
    ndim = len(dim)
    if (ndim == 1):
        nTrace = 1
        npix = dim[0]
    else:
        nTrace = dim[0]
        npix = dim[1]

    if ycen is None:
        if ndim == 1:
            ycen_out = np.arange(npix, dtype='int')
        elif ndim == 2:
            ycen_out = np.outer(np.ones(nTrace, dtype='int'), np.arange(npix, dtype='int'), )
        else:
            raise ValueError('trace is not 1 or 2 dimensional')
    else:
        ycen_out = ycen.T
        ycen_out = np.rint(ycen_out).astype(int)

    if ((np.size(trace) != np.size(ycen_out)) | (np.shape(trace) != np.shape(ycen_out))):
        raise ValueError('Number of elements and shape of trace and ycen must be equal')



    left = trace - radius
    right = trace + radius
    fextract = extract_asymbox2(image, left, right, ycen_out)

    return fextract


def extract_optimal(waveimg, imgminsky, ivar, mask, oprof, skyimg, rn2_img, box_radius, specobj):

    nspat = imgminsky.shape[1]
    nspec = imgminsky.shape[0]

    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    var_no = np.abs(skyimg - np.sqrt(2.0) * np.sqrt(rn2_img)) + rn2_img

    ispec, ispat = np.where(oprof > 0.0)
    mincol = np.min(ispat)
    maxcol = np.max(ispat) + 1
    nsub = maxcol - mincol

    mask_sub = mask[:,mincol:maxcol]
    wave_sub = waveimg[:,mincol:maxcol]
    ivar_sub = np.fmax(ivar[:,mincol:maxcol],0.0) # enforce positivity since these are used as weights
    vno_sub = np.fmax(var_no[:,mincol:maxcol],0.0)

    rn2_sub = rn2_img[:,mincol:maxcol]
    img_sub = imgminsky[:,mincol:maxcol]
    sky_sub = skyimg[:,mincol:maxcol]
    oprof_sub = oprof[:,mincol:maxcol]
    # enforce normalization and positivity of object profiles
    norm = np.nansum(oprof_sub,axis = 1)
    norm_oprof = np.outer(norm, np.ones(nsub))
    oprof_sub = np.fmax(oprof_sub/norm_oprof, 0.0)

    ivar_denom = np.nansum(mask_sub*oprof_sub, axis=1)
    mivar_num = np.nansum(mask_sub*ivar_sub*oprof_sub**2, axis=1)
    mivar_opt = mivar_num/(ivar_denom + (ivar_denom == 0.0))
    flux_opt = np.nansum(mask_sub*ivar_sub*img_sub*oprof_sub, axis=1)/(mivar_num + (mivar_num == 0.0))
    # Optimally extracted noise variance (sky + read noise) only. Since
    # this variance is not the same as that used for the weights, we
    # don't get the usual cancellation. Additional denom factor is the
    # analog of the numerator in Horne's variance formula. Note that we
    # are only weighting by the profile (ivar_sub=1) because
    # otherwise the result depends on the signal (bad).
    nivar_num =np.nansum(mask_sub*oprof_sub**2, axis=1) # Uses unit weights
    nvar_opt = ivar_denom*((mask_sub*vno_sub*oprof_sub**2).sum(axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    nivar_opt = 1.0/(nvar_opt + (nvar_opt == 0.0))
    # Optimally extract sky and (read noise)**2 in a similar way
    sky_opt = ivar_denom*(np.nansum(mask_sub*sky_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    rn2_opt = ivar_denom*(np.nansum(mask_sub*rn2_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    rn_opt = np.sqrt(rn2_opt)
    rn_opt[np.isnan(rn_opt)]=0.0

    tot_weight = np.nansum(mask_sub*ivar_sub*oprof_sub, axis=1)
    mask_opt = (tot_weight > 0.0) & (mivar_num > 0.0) & (ivar_denom > 0.0)
    frac_use = np.nansum((mask_sub*ivar_sub > 0.0)*oprof_sub, axis=1)
    # Use the same weights = oprof^2*mivar for the wavelenghts as the flux.
    # Note that for the flux, one of the oprof factors cancels which does
    # not for the wavelengths.
    wave_opt = np.nansum(mask_sub*ivar_sub*wave_sub*oprof_sub**2, axis=1)/(mivar_num + (mivar_num == 0.0))
    # Interpolate wavelengths over masked pixels
    badwvs = (mivar_num <= 0) | (np.isfinite(wave_opt) == False) | (wave_opt <= 0.0)
    if badwvs.any():
        oprof_smash = np.nansum(oprof_sub**2, axis=1)
        # Can we use the profile average wavelengths instead?
        oprof_good = badwvs & (oprof_smash > 0.0)
        if oprof_good.any():
            wave_opt[oprof_good] = np.nansum(wave_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)/np.nansum(oprof_sub[oprof_good,:]**2, axis=1)
        oprof_bad = badwvs & ((oprof_smash <= 0.0) | (np.isfinite(oprof_smash) == False) | (wave_opt <= 0.0) | (np.isfinite(wave_opt) == False))
        if oprof_bad.any():
            # For pixels with completely bad profile values, interpolate from trace.
            f_wave = RectBivariateSpline(spec_vec,spat_vec, waveimg)
            wave_opt[oprof_bad] = f_wave(specobj.trace_spec[oprof_bad], specobj.trace_spat[oprof_bad],grid=False)

    flux_model = np.outer(flux_opt,np.ones(nsub))*oprof_sub
    chi2_num = np.nansum((img_sub - flux_model)**2*ivar_sub*mask_sub,axis=1)
    chi2_denom = np.fmax(np.nansum(ivar_sub*mask_sub > 0.0, axis=1) - 1.0, 1.0)
    chi2 = chi2_num/chi2_denom

    # Fill in the optimally extraction tags
    specobj.optimal['WAVE_OPT'] = wave_opt    # Optimally extracted wavelengths
    specobj.optimal['FLUX_OPT'] = flux_opt    # Optimally extracted flux
    specobj.optimal['IVAR_OPT'] = mivar_opt   # Inverse variance of optimally extracted flux using modelivar image
    specobj.optimal['NIVAR_OPT'] = nivar_opt  # Optimally extracted noise variance (sky + read noise) only
    specobj.optimal['MASK_OPT'] = mask_opt    # Mask for optimally extracted flux
    specobj.optimal['SKY_OPT'] = sky_opt      # Optimally extracted sky
    specobj.optimal['RN_OPT'] = rn_opt        # Square root of optimally extracted read noise squared
    specobj.optimal['FRAC_USE'] = frac_use    # Fraction of pixels in the object profile subimage used for this extraction
    specobj.optimal['CHI2'] = chi2            # Reduced chi2 of the model fit for this spectral pixel

    # Fill in the boxcar extraction tags
    flux_box  = extract_boxcar(imgminsky*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    # Denom is computed in case the trace goes off the edge of the image
    box_denom = extract_boxcar(waveimg*mask > 0.0, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    wave_box  = extract_boxcar(waveimg*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)/(box_denom + (box_denom == 0.0))
    varimg = 1.0/(ivar + (ivar == 0.0))
    var_box  = extract_boxcar(varimg*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    nvar_box  = extract_boxcar(var_no*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    sky_box  = extract_boxcar(skyimg*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    rn2_box  = extract_boxcar(rn2_img*mask, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    rn_posind = (rn2_box > 0.0)
    rn_box = np.zeros(rn2_box.shape,dtype=float)
    rn_box[rn_posind] = np.sqrt(rn2_box[rn_posind])
    pixtot  = extract_boxcar(ivar*0 + 1.0, specobj.trace_spat,box_radius, ycen = specobj.trace_spec)
    # If every pixel is masked then mask the boxcar extraction
    mask_box = (extract_boxcar(ivar*mask == 0.0, specobj.trace_spat,box_radius, ycen = specobj.trace_spec) != pixtot)

    bad_box = (wave_box <= 0.0) | (np.isfinite(wave_box) == False) | (box_denom == 0.0)
    # interpolate bad wavelengths over masked pixels
    if bad_box.any():
        f_wave = RectBivariateSpline(spec_vec, spat_vec, waveimg)
        wave_box[bad_box] = f_wave(specobj.trace_spec[bad_box], specobj.trace_spat[bad_box],grid=False)

    ivar_box = 1.0/(var_box + (var_box == 0.0))
    nivar_box = 1.0/(nvar_box + (nvar_box == 0.0))

    specobj.boxcar['WAVE_BOX'] = wave_box
    specobj.boxcar['FLUX_BOX'] = flux_box*mask_box
    specobj.boxcar['IVAR_BOX'] = ivar_box*mask_box
    specobj.boxcar['NIVAR_BOX'] = nivar_box*mask_box
    specobj.boxcar['MASK_BOX'] = mask_box
    specobj.boxcar['SKY_BOX'] = sky_box
    specobj.boxcar['RN_BOX'] = rn_box

    return




def findfwhm(model, sig_x):
    """ Calculate the spatial FWHM from an object profile. Utitlit routine for fit_profile

    Parameters
    ----------
    model :   numpy float 2-d array [nspec, nspat]
    x :

    Returns
    -------
    peak :  Peak value of the profile model
    peak_x:  sig_x location where the peak value is obtained
    lwhm:   Value of sig_x at the left width at half maximum
    rwhm:   Value of sig_x at the right width at half maximum

    Revision History
    ----------------
    11-Mar-2005  Written by J. Hennawi and S. Burles David Schlegel, Princeton.
    28-May-2018  Ported to python by J. Hennawi
    """


    peak = (model*(np.abs(sig_x) < 1.)).max()
    peak_x = sig_x[(model*(np.abs(sig_x) < 1.)).argmax()]

    lrev = ((sig_x < peak_x) & (model < 0.5*peak))[::-1]
    lind, = np.where(lrev)
    if(lind.size > 0):
        lh = lind.min()
        lwhm = (sig_x[::-1])[lh]
    else:
        lwhm = -0.5*2.3548

    rind, = np.where((sig_x > peak_x) & (model < 0.5*peak))
    if(rind.size > 0):
        rh = rind.min()
        rwhm = sig_x[rh]
    else:
        rwhm = 0.5 * 2.3548

    return (peak, peak_x, lwhm, rwhm)



def fit_profile_qa(x_tot,y_tot, model_tot, l_limit = None, r_limit = None, ind = None,
                   title =' ', xtrunc = 1e6, xlim = None, ylim = None, qafile = None):


    # Plotting pre-amble
    plt.close("all")
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    width = 10.0 # Golden ratio 1.618
    fig, ax = plt.subplots(1, figsize=(width, width/1.618))

    if ind is None:
        ind = np.slice(x_tot.size)

    x = x_tot.flat[ind]
    y = y_tot.flat[ind]
    model = model_tot.flat[ind]

    ax.plot(x,y,color='k',marker='o',markersize= 0.3, mfc='k',fillstyle='full',linestyle='None')


    max_model = np.fmin(model.max(), 2.0)
    if xlim is None:
        goodpix = model > 0.001*max_model
        if goodpix.any():
            minx = np.fmax(1.1*(x[goodpix]).min(), -xtrunc)
            maxx = np.fmin(1.1*(x[goodpix]).max(),  xtrunc)
        else:
            minx = -5.0
            maxx = 5.0
    else:
        minx = -xlim
        maxx = xlim
    xlimit = (minx, maxx)

    nsamp = 150
    half_bin = (maxx - minx)/nsamp/2.0
    if ylim is None:
        ymax = np.fmax(1.5*model.max(), 0.3)
        ymin = np.fmin(-0.1*ymax, -0.05)
        ylim = (ymin, ymax)
    plot_mid = (np.arange(nsamp) + 0.5)/nsamp*(maxx - minx) + minx

    y20 = np.zeros(nsamp)
    y80 = np.zeros(nsamp)
    y50 = np.zeros(nsamp)
    model_samp = np.zeros(nsamp)
    nbin = np.zeros(nsamp)

    for i in range(nsamp):
        dist = np.abs(x - plot_mid[i])
        close = dist < half_bin
        yclose = y[close]
        nclose = close.sum()
        nbin[i] = nclose
        if close.any():
            closest = (dist[close]).argmin()
            model_samp[i] = (model[close])[closest]
        if nclose > 3:
            s = yclose.argsort()
            y50[i] = yclose[s[int(np.rint((nclose - 1)*0.5))]]
            y80[i] = yclose[s[int(np.rint((nclose - 1)*0.8))]]
            y20[i] = yclose[s[int(np.rint((nclose - 1)*0.2))]]
            ax.plot([plot_mid[i],plot_mid[i]], [y20[i],y80[i]], linewidth=1.0, color='cornflowerblue')

    icl = nbin > 3
    if icl.any():
        ax.plot(plot_mid[icl],y50[icl],marker = 'o', color='lime', markersize=2, fillstyle='full', linestyle='None')
    else:
        ax.plot(plot_mid, y50, marker='o', color='lime', markersize=2, fillstyle = 'full', linestyle='None')

    isort = x.argsort()
    ax.plot(x[isort], model[isort], color='red', linewidth=1.0)



    if l_limit is not None:
        ax.axvline(x =l_limit, color='orange',linewidth=1.0)
    if r_limit is not None:
        ax.axvline(x=r_limit, color='orange',linewidth=1.0)

    ax.set_xlim(xlimit)
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.set_xlabel(r'$x/\sigma$')
    ax.set_ylabel('Normalized Profile')

    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.show()
    return


def fit_profile(image, ivar, waveimg, trace_in, wave, flux, fluxivar,
                thisfwhm=4.0, MAX_TRACE_CORR = 2.0, SN_GAUSS = 3.0, wvmnx = [2900.0,30000.0],
                hwidth = None, PROF_NSIGMA = None, NO_DERIV = False, GAUSS = False):

    """Fit a non-parametric object profile to an object spectrum, unless the S/N ratio is low (> SN_GAUSS) in which
    fit a simple Gaussian. Port of IDL LOWREDUX long_gprofile.pro

     Parameters
     ----------
     image : numpy float 2-d array [nspec, nspat]
         sky-subtracted image
     ivar : numpy float 2-d array [nspec, nspat]
         inverse variance of sky-subtracted image
     waveimg numpy float 2-d array [nspec, nspat]
         2-d wavelength map
     trace_in : numpy 1-d array [nspec]
         object trace
     wave : numpy 1-d array [nspec]
         extracted wavelength of spectrum
     flux : numpy 1-d array [nspec]
         extracted flux of spectrum
     fluxivar : numpy 1-d array [nspec]
         inverse variance of extracted flux spectrum


    Optional Parameters
    ----------
    thisfwhm : float
         fwhm of the object trace
    MAX_TRACE_CORR : float [default = 2.0]
         maximum trace correction to apply
    SN_GAUSS : float [default = 3.0]
         S/N ratio below which code just uses a Gaussian
    wvmnx : float [default = [2900.0,30000.0]
         wavelength range of usable part of spectrum
    hwidth : float [default = None]
         object maskwidth determined from object finding algorithm. If = None,
         code defaults to use 3.0*(np.max(thisfwhm) + 1.0)
    PROF_NSIGMA : float [default = None]
         Number of sigma to include in the profile fitting. This option is only needed for bright objects that are not
         point sources, which allows the profile fitting to fit the high S/N wings (rather than the default behavior
         which truncates exponentially). This allows for extracting all the flux and results in better sky-subtraction
         for bright extended objects.
    NO_DERIV : boolean [default = False]
         disables determination of derivatives and exponential tr

     Returns
     -------
     :func:`tuple`
         A tuple containing the (sset, outmask, yfit, reduced_chi), where

            sset: object
               bspline object
            outmask: : :class:`numpy.ndarray`
               output mask which the same size as xdata
            yfit  : :class:`numpy.ndarray`
               result of the bspline fit (same size as xdata)
            reduced_chi: float
               value of the reduced chi^2
     """

    from scipy.interpolate import interp1d
    from scipy.ndimage.filters import median_filter
    from pypit.core.pydl import bspline
    from pypit.core.pydl import iterfit as  bspline_iterfit
    from pypit.core.pydl import djs_maskinterp
    from pypit.arutils import bspline_profile
    from astropy.stats import sigma_clipped_stats
    from scipy.special import erfcinv

    if hwidth is None: 3.0*(np.max(thisfwhm) + 1.0)
    if PROF_NSIGMA is not None:
        NO_DERIV = True

    thisfwhm = np.fmax(thisfwhm,1.0) # require the FWHM to be greater than 1 pixel

    xnew = trace_in

    nspat = image.shape[1]
    nspec = image.shape[0]

    # create some images we will need
    profile_model = np.zeros((nspec,nspat))
    sub_obj = image
    sub_ivar = ivar
    sub_wave = waveimg
    sub_trace = trace_in
    sub_x = np.arange(nspat)
    sn2_sub = np.zeros((nspec,nspat))
    spline_sub = np.zeros((nspec,nspat))


    flux_sm = median_filter(flux, size=5, mode = 'reflect')
    fluxivar_sm =  median_filter(fluxivar, size = 5, mode = 'reflect')
#    flux_sm = djs_median(flux, width = 5, boundary = 'reflect')
    #    fluxivar_sm =  djs_median(fluxivar, width = 5, boundary = 'reflect')
    fluxivar_sm = fluxivar_sm*(fluxivar > 0.0)

    indsp = (wave > wvmnx[0]) & (wave < wvmnx[1]) & \
             np.isfinite(flux_sm) & (flux_sm < 5.0e5) &  \
             (flux_sm > -1000.0) & (fluxivar_sm > 0.0)


    b_answer, bmask   = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp],kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    b_answer, bmask2  = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask, kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    c_answer, cmask   = bspline_iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask2,kwargs_bspline={'everyn': 30}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    spline_flux, _ = b_answer.value(wave[indsp])
    cont_flux, _ = c_answer.value(wave[indsp])

    sn2 = (np.fmax(spline_flux*(np.sqrt(np.fmax(fluxivar_sm[indsp], 0))*bmask2),0))**2
    ind_nonzero = (sn2 > 0)
    nonzero = np.sum(ind_nonzero)
    if(nonzero >0):
        (mean, med_sn2, stddev) = sigma_clipped_stats(sn2,sigma_lower=3.0,sigma_upper=5.0)
    else: med_sn2 = 0.0
    sn2_med = median_filter(sn2, size=9, mode='reflect')
    #sn2_med = djs_median(sn2, width = 9, boundary = 'reflect')
    igood = (ivar > 0.0)
    ngd = np.sum(igood)
    if(ngd > 0):
        isrt = np.argsort(wave[indsp])
        sn2_interp = interp1d((wave[indsp])[isrt],sn2_med[isrt],assume_sorted=False, bounds_error=False,fill_value = 'extrapolate')
        sn2_sub[igood] = sn2_interp(sub_wave[igood])
    msgs.info('sqrt(med(S/N)^2) = ' + "{:5.2f}".format(np.sqrt(med_sn2)))

    min_wave = np.min(wave[indsp])
    max_wave = np.max(wave[indsp])
    spline_flux1 = np.zeros(nspec)
    cont_flux1 = np.zeros(nspec)
    sn2_1 = np.zeros(nspec)
    ispline = (wave >= min_wave) & (wave <= max_wave)
    spline_tmp, _ = b_answer.value(wave[ispline])
    spline_flux1[ispline] = spline_tmp
    cont_tmp, _ = c_answer.value(wave[ispline])
    cont_flux1[ispline] = cont_tmp
    s2_1_interp = interp1d((wave[indsp])[isrt], sn2[isrt],assume_sorted=False, bounds_error=False,fill_value = 'extrapolate')
    sn2_1[ispline] = s2_1_interp(wave[ispline])
    bmask = np.zeros(nspec,dtype='bool')
    bmask[indsp] = bmask2
    spline_flux1 = djs_maskinterp(spline_flux1,(bmask == False))
    cmask2 = np.zeros(nspec,dtype='bool')
    cmask2[indsp] = cmask
    cont_flux1 = djs_maskinterp(cont_flux1,(cmask2 == False))

    (_, _, sigma1) = sigma_clipped_stats(flux[indsp],sigma_lower=3.0,sigma_upper=5.0)

    if(med_sn2 <= 2.0):
        spline_sub[igood]= np.fmax(sigma1,0)
    else:
        if((med_sn2 <=5.0) and (med_sn2 > 2.0)):
            spline_flux1 = cont_flux1
        # Interp over points <= 0 in boxcar flux or masked points using cont model
        badpix = (spline_flux1 <= 0.5) | (bmask == False)
        goodval = (cont_flux1 > 0.0) & (cont_flux1 < 5e5)
        indbad1 = badpix & goodval
        nbad1 = np.sum(indbad1)
        if(nbad1 > 0):
            spline_flux1[indbad1] = cont_flux1[indbad1]
        indbad2 = badpix & ~goodval
        nbad2 = np.sum(indbad2)
        ngood0 = np.sum(~badpix)
        if((nbad2 > 0) or (ngood0 > 0)):
            spline_flux1[indbad2] = np.median(spline_flux1[~badpix])
        # take a 5-pixel median to filter out some hot pixels
        spline_flux1 = median_filter(spline_flux1,size=5,mode ='reflect')
        # Create the normalized object image
        if(ngd > 0):
            isrt = np.argsort(wave)
            spline_sub_interp = interp1d(wave[isrt],spline_flux1[isrt],assume_sorted=False, bounds_error=False,fill_value = 'extrapolate')
            spline_sub[igood] = spline_sub_interp(sub_wave[igood])
        else:
            spline_sub[igood] = np.fmax(sigma1, 0)

    norm_obj = (spline_sub != 0.0)*sub_obj/(spline_sub + (spline_sub == 0.0))
    norm_ivar = sub_ivar*spline_sub**2

    # Cap very large inverse variances
    ivar_mask = (norm_obj > -0.2) & (norm_obj < 0.7) & (sub_ivar > 0.0) & np.isfinite(norm_obj) & np.isfinite(norm_ivar)
    norm_ivar = norm_ivar*ivar_mask
    good = (norm_ivar.flatten() > 0.0)
    ngood = np.sum(good)


    xtemp = (np.cumsum(np.outer(4.0 + np.sqrt(np.fmax(sn2_1, 0.0)),np.ones(nspat)))).reshape((nspec,nspat))
    xtemp = xtemp/xtemp.max()


    # norm_x is the x position along the image centered on the object  trace
    norm_x = np.outer(np.ones(nspec), sub_x) - np.outer(sub_trace,np.ones(nspat))

    sigma = np.full(nspec, thisfwhm/2.3548)
    fwhmfit = sigma*2.3548
    trace_corr = np.zeros(nspec)

    # If we have too few pixels to fit a profile or S/N is too low, just use a Gaussian profile
    # TODO Clean up the logic below. It is formally correct but
    # redundant since no trace correction has been created or applied yet.  I'm only postponing doing it
    # to preserve this syntax for later when it is needed
    if((ngood < 10) or (med_sn2 < SN_GAUSS**2) or (GAUSS is True)):
        msgs.info("Too few good pixels or S/N <" + "{:5.1f}".format(SN_GAUSS) + " or GAUSS flag set")
        msgs.info("Returning Gaussian profile")
        sigma_x = norm_x/(np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr,np.ones(nspat)))
        profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0*np.pi)*(sigma_x**2 < 25.)
        msgs.info("FWHM="  + "{:6.2f}".format(thisfwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)))
        nxinf = np.sum(np.isfinite(xnew) == False)
        if(nxinf != 0):
            msgs.warn("Nan pixel values in trace correction")
            msgs.warn("Returning original trace....")
            xnew = trace_in
        inf = np.isfinite(profile_model) == False
        ninf = np.sum(inf)
        if (ninf != 0):
            msgs.warn("Nan pixel values in object profile... setting them to zero")
            profile_model[inf] = 0.0
        # Normalize profile
        norm = np.outer(np.sum(profile_model,1),np.ones(nspat))
        if(np.sum(norm) > 0.0):
            profile_model = profile_model/norm
        if ngood > 0:
            title_string = ' '
            indx = good
        else:
            title_string = 'No good pixels, showing all'
            indx = None

        fit_profile_qa(sigma_x, norm_obj, profile_model, title = title_string, ind=indx, xtrunc= 7.0)
        return (profile_model, xnew, fwhmfit, med_sn2)


    msgs.info("Gaussian vs b-spline of width " + "{:6.2f}".format(thisfwhm) + " pixels")
    area = 1.0
    sigma_x = norm_x / (np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr, np.ones(nspat)))

    mask = np.full(nspec*nspat, False, dtype=bool)

    # The following lines set the limits for the b-spline fit
    limit = erfcinv(0.1/np.sqrt(med_sn2))*np.sqrt(2.0)
    if(PROF_NSIGMA is None):
        sinh_space = 0.25*np.log10(np.fmax((1000./np.sqrt(med_sn2)),10.))
        abs_sigma = np.fmin((np.abs(sigma_x.flat[good])).max(),2.0*limit)
        min_sigma = np.fmax(sigma_x.flat[good].min(), (-abs_sigma))
        max_sigma = np.fmin(sigma_x.flat[good].max(), (abs_sigma))
        nb = (np.arcsinh(abs_sigma)/sinh_space).astype(int) + 1
    else:
        msgs.info("Using PROF_NSIGMA= " + "{:6.2f}".format(PROF_NSIGMA) + " for extended/bright objects")
        nb = np.round(PROF_NSIGMA > 10)
        max_sigma = PROF_NSIGMA
        min_sigma = -1*PROF_NSIGMA
        sinh_space = np.arcsinh(PROF_NSIGMA)/nb

    rb = np.sinh((np.arange(nb) + 0.5) * sinh_space)
    bkpt = np.concatenate([(-rb)[::-1], rb])
    keep = ((bkpt >= min_sigma) & (bkpt <= max_sigma))
    bkpt = bkpt[keep]

    # Attempt B-spline first
    GOOD_PIX = (sn2_sub > SN_GAUSS**2) & (norm_ivar > 0)
    IN_PIX   = (sigma_x >= min_sigma) & (sigma_x <= max_sigma) & (norm_ivar > 0)
    ngoodpix = np.sum(GOOD_PIX)
    ninpix     = np.sum(IN_PIX)

    if (ngoodpix >= 0.2*ninpix):
        inside,  = np.where((GOOD_PIX & IN_PIX).flatten())
    else:
        inside, = np.where(IN_PIX.flatten())


    si = inside[np.argsort(sigma_x.flat[inside])]
    sr = si[::-1]

    bset, bmask = bspline_iterfit(sigma_x.flat[si],norm_obj.flat[si], invvar = norm_ivar.flat[si]
                           , nord = 4, bkpt = bkpt, maxiter = 15, upper = 1, lower = 1)
    mode_fit, _ = bset.value(sigma_x.flat[si])
    median_fit = np.median(norm_obj[norm_ivar > 0.0])

    # TODO I don't follow the logic behind this statement but I'm leaving it for now. If the median is large it is used, otherwise we  user zero???
    if (np.abs(median_fit) > 0.01):
        msgs.info("Median flux level in profile is not zero: median = " + "{:7.4f}".format(median_fit))
    else:
        median_fit = 0.0

    # Find the peak and FWHM based this profile fit
    (peak, peak_x, lwhm, rwhm) = findfwhm(mode_fit - median_fit, sigma_x.flat[si])
    trace_corr = np.full(nspec, peak_x)
    min_level = peak*np.exp(-0.5*limit**2)

    bspline_fwhm = (rwhm - lwhm)*thisfwhm/2.3548
    msgs.info("Bspline FWHM: " + "{:7.4f}".format(bspline_fwhm) + ", compared to initial object finding FWHM: " + "{:7.4f}".format(thisfwhm) )
    sigma = sigma * (rwhm-lwhm)/2.3548

    limit = limit * (rwhm-lwhm)/2.3548

    rev_fit = mode_fit[::-1]
    lind, = np.where(((rev_fit < (min_level+median_fit)) & (sigma_x.flat[sr] < peak_x)) | (sigma_x.flat[sr] < (peak_x-limit)))
    if (lind.size > 0):
        lp = lind.min()
        l_limit = sigma_x.flat[sr[lp]]
    else:
        l_limit = min_sigma

    rind, = np.where(((mode_fit < (min_level+median_fit)) & (sigma_x.flat[si] > peak_x)) | (sigma_x.flat[si] > (peak_x+limit)))
    if (rind.size > 0):
        rp = rind.min()
        r_limit = sigma_x.flat[si[rp]]
    else:
        r_limit = max_sigma

    msgs.info("Trace limits: limit = " + "{:7.4f}".format(limit) + ", min_level = " + "{:7.4f}".format(min_level) +
              ", l_limit = " + "{:7.4f}".format(l_limit) + ", r_limit = " + "{:7.4f}".format(r_limit))

    # Just grab the data points within the limits
    mask[si]=((norm_ivar.flat[si] > 0) & (np.abs(norm_obj.flat[si] - mode_fit) < 0.1))
    inside, = np.where((sigma_x.flat[si] > l_limit) & (sigma_x.flat[si] < r_limit) & mask[si])
    ninside = inside.size


    # If we have too few pixels after this step, then again just use a Gaussian profile and return. Note that
    # we are following the original IDL code here and not using the improved sigma and trace correction for the
    # profile at this stage since ninside is so small.
    if(ninside < 10):
        msgs.info("Too few pixels inside l_limit and r_limit")
        msgs.info("Returning Gaussian profile")
        profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0*np.pi)*(sigma_x**2 < 25.)
        msgs.info("FWHM="  + "{:6.2f}".format(thisfwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)))
        nxinf = np.sum(np.isfinite(xnew) == False)
        if(nxinf != 0):
            msgs.warn("Nan pixel values in trace correction")
            msgs.warn("Returning original trace....")
            xnew = trace_in
        inf = np.isfinite(profile_model) == False
        ninf = np.sum(inf)
        if (ninf != 0):
            msgs.warn("Nan pixel values in object profile... setting them to zero")
            profile_model[inf] = 0.0
        # Normalize profile
        norm = np.outer(np.sum(profile_model,1),np.ones(nspat))
        if(np.sum(norm) > 0.0):
            profile_model = profile_model/norm

        fit_profile_qa(sigma_x, norm_obj, profile_model, l_limit = l_limit, r_limit = r_limit, ind =good, xlim = 7.0)
        return (profile_model, xnew, fwhmfit, med_sn2)

    sigma_iter = 3
    isort =  (xtemp.flat[si[inside]]).argsort()
    inside = si[inside[isort]]
    pb =np.ones(inside.size)

    # ADD the loop here later
    for iiter in range(1,sigma_iter + 1):
        mode_zero, _ = bset.value(sigma_x.flat[inside])
        mode_zero = mode_zero*pb

        mode_min05, _ = bset.value(sigma_x.flat[inside]-0.5)
        mode_plu05, _ = bset.value(sigma_x.flat[inside]+0.5)
        mode_shift = (mode_min05  - mode_plu05)*pb*((sigma_x.flat[inside] > (l_limit + 0.5)) &
                                                (sigma_x.flat[inside] < (r_limit - 0.5)))

        mode_by13, _ = bset.value(sigma_x.flat[inside]/1.3)
        mode_stretch = mode_by13*pb/1.3 - mode_zero

        nbkpts = (np.log10(np.fmax(med_sn2, 11.0))).astype(int)

        xx = np.sum(xtemp, 1)/nspat
        profile_basis = np.column_stack((mode_zero,mode_shift))

        mode_shift_out = bspline_profile(xtemp.flat[inside], norm_obj.flat[inside], norm_ivar.flat[inside], profile_basis
                                      ,maxiter=1,kwargs_bspline= {'nbkpts':nbkpts})
        mode_shift_set = mode_shift_out[0]
        temp_set = bspline(None, fullbkpt = mode_shift_set.breakpoints,nord=mode_shift_set.nord)
        temp_set.coeff = mode_shift_set.coeff[0, :]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_shift_set.coeff[1, :]
        h1, _ = temp_set.value(xx)
        ratio_10 = (h1/(h0 + (h0 == 0.0)))
        delta_trace_corr = ratio_10/(1.0 + np.abs(ratio_10)/0.1)
        trace_corr = trace_corr + delta_trace_corr

        profile_basis = np.column_stack((mode_zero,mode_stretch))
        mode_stretch_out = bspline_profile(xtemp.flat[inside], norm_obj.flat[inside], norm_ivar.flat[inside], profile_basis,
                                            maxiter=1,fullbkpt = mode_shift_set.breakpoints)
        mode_stretch_set = mode_stretch_out[0]
        temp_set = bspline(None, fullbkpt = mode_stretch_set.breakpoints,nord=mode_stretch_set.nord)
        temp_set.coeff = mode_stretch_set.coeff[0, :]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_stretch_set.coeff[1, :]
        h2, _ = temp_set.value(xx)
        h0 = np.fmax(h0 + h2*mode_stretch.sum()/mode_zero.sum(),0.1)
        ratio_20 = (h2 / (h0 + (h0 == 0.0)))
        sigma_factor = 0.3 * ratio_20 / (1.0 + np.abs(ratio_20))

        msgs.info("Iteration# " + "{:3d}".format(iiter))
        msgs.info("Median abs value of trace correction = " + "{:8.3f}".format(np.median(np.abs(delta_trace_corr))))
        msgs.info("Median abs value of width correction = " + "{:8.3f}".format(np.median(np.abs(sigma_factor))))

        sigma = sigma*(1.0 + sigma_factor)
        area = area * h0/(1.0 + sigma_factor)

        sigma_x = norm_x / (np.outer(sigma, np.ones(nspat)) - np.outer(trace_corr, np.ones(nspat)))

        # Update the profile B-spline fit for the next iteration
        if iiter < sigma_iter-1:
            ss = sigma_x.flat[inside].argsort()
            pb = (np.outer(area, np.ones(nspat,dtype=float))).flat[inside]
            keep = (bkpt >= sigma_x.flat[inside].min()) & (bkpt <= sigma_x.flat[inside].max())
            if keep.sum() == 0:
                keep = np.ones(bkpt.size, type=bool)
            bset_out = bspline_profile(sigma_x.flat[inside[ss]],norm_obj.flat[inside[ss]],norm_ivar.flat[inside[ss]],pb[ss],
                                    nord = 4, bkpt=bkpt[keep],maxiter=2)
            bset = bset_out[0] # This updated bset used for the next set of trace corrections

    # Apply trace corrections only if they are small (added by JFH)
    if np.median(np.abs(trace_corr*sigma)) < MAX_TRACE_CORR:
        xnew = trace_corr * sigma + trace_in
    else:
        xnew = trace_in

    fwhmfit = sigma*2.3548
    ss=sigma_x.flatten().argsort()
    inside, = np.where((sigma_x.flat[ss] >= min_sigma) &
                       (sigma_x.flat[ss] <= max_sigma) &
                       mask[ss] &
                       np.isfinite(norm_obj.flat[ss]) &
                       np.isfinite(norm_ivar.flat[ss]))
    pb = (np.outer(area, np.ones(nspat,dtype=float)))
    bset_out = bspline_profile(sigma_x.flat[ss[inside]],norm_obj.flat[ss[inside]], norm_ivar.flat[ss[inside]], pb.flat[ss[inside]],
                            nord=4, bkpt = bkpt, upper = 10, lower=10)
    bset = bset_out[0]
    outmask = bset_out[1]

    # inmask = False for pixels within (min_sigma, max_sigma), True outside
    inmask = (sigma_x.flatten() > min_sigma) & (sigma_x.flatten() < max_sigma)
    full_bsp = np.zeros(nspec*nspat, dtype=float)
    sigma_x_inmask = sigma_x.flat[inmask]
    yfit_out, _  = bset.value(sigma_x_inmask)
    full_bsp[inmask] = yfit_out
    isrt = sigma_x_inmask.argsort()
    (peak, peak_x, lwhm, rwhm) = findfwhm(yfit_out[isrt] - median_fit, sigma_x_inmask[isrt])


    left_bool = (((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] < peak_x)) | (sigma_x.flat[ss] < (peak_x-limit)))[::-1]
    ind_left, = np.where(left_bool)
    lp = np.fmax(ind_left.min(), 0)
    righ_bool = ((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] > peak_x))  | (sigma_x.flat[ss] > (peak_x+limit))
    ind_righ, = np.where(righ_bool)
    rp = np.fmax(ind_righ.min(), 0)
    l_limit = ((sigma_x.flat[ss])[::-1])[lp] - 0.1
    r_limit = sigma_x.flat[ss[rp]] + 0.1

    while True:
        l_limit += 0.1
        l_fit, _ = bset.value(np.asarray([l_limit]))
        l2, _ = bset.value(np.asarray([l_limit])* 0.9)
        l_deriv = (np.log(l2[0]) - np.log(l_fit[0]))/(0.1*l_limit)
        if (l_deriv < -1.0) | (l_limit >= -1.0):
            break

    while True:
        r_limit -= 0.1
        r_fit, _ = bset.value(np.asarray([r_limit]))
        r2, _ = bset.value(np.asarray([r_limit])* 0.9)
        r_deriv = (np.log(r2[0]) - np.log(r_fit[0]))/(0.1*r_limit)
        if (r_deriv > 1.0) | (r_limit <= 1.0):
            break


    # JXP kludge
    if PROF_NSIGMA is not None:
       #By setting them to zero we ensure QA won't plot them in the profile QA.
       l_limit = 0.0
       r_limit = 0.0


    # Hack to fix degenerate profiles which have a positive derivative
    if (l_deriv < 0) and (r_deriv > 0) and NO_DERIV is False:
        left = sigma_x.flatten() < l_limit
        full_bsp[left] =  np.exp(-(sigma_x.flat[left]-l_limit)*l_deriv) * l_fit
        right = sigma_x.flatten() > r_limit
        full_bsp[right] = np.exp(-(sigma_x.flat[right] - r_limit) * r_deriv) * r_fit

    # Final object profile
    full_bsp = full_bsp.reshape(nspec,nspat)
    profile_model = full_bsp*pb
    res_mode = (norm_obj.flat[ss[inside]] - profile_model.flat[ss[inside]])*np.sqrt(norm_ivar.flat[ss[inside]])
    chi_good = (outmask == True) & (norm_ivar.flat[ss[inside]] > 0)
    chi_med = np.median(res_mode[chi_good]**2)
    chi_zero = np.median(norm_obj.flat[ss[inside]]**2*norm_ivar.flat[ss[inside]])

    msgs.info("----------  Results of Profile Fit ----------")
    msgs.info(" min(fwhmfit)={:5.2f}".format(fwhmfit.min()) +
              " max(fwhmfit)={:5.2f}".format(fwhmfit.max()) + " median(chi)={:5.2f}".format(chi_med) +
              " nbkpts={:2d}".format(bkpt.size))

    nxinf = np.sum(np.isfinite(xnew) == False)
    if (nxinf != 0):
        msgs.warn("Nan pixel values in trace correction")
        msgs.warn("Returning original trace....")
        xnew = trace_in
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if (ninf != 0):
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    # Normalize profile
    norm = np.outer(np.sum(profile_model, 1), np.ones(nspat))
    if (np.sum(norm) > 0.0):
        profile_model = profile_model / norm
    msgs.info("FWHM="  + "{:6.2f}".format(thisfwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)))
    fit_profile_qa(sigma_x, norm_obj/(pb + (pb == 0.0)), full_bsp, l_limit = l_limit, r_limit = r_limit, ind = ss[inside], xlim = PROF_NSIGMA)

    return (profile_model, xnew, fwhmfit, med_sn2)




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



