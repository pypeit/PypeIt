""" Module for PypeIt extraction code
"""
import time
import copy
import inspect

import numpy as np
import scipy


#from matplotlib import gridspec, font_manager

from astropy import stats

from pypeit import msgs
from pypeit.core import pydl
from pypeit import utils
from pypeit.core import pixels
from pypeit import ginga
from matplotlib import pyplot as plt
from pypeit.core import trace_slits
from pypeit.core import arc
from scipy import interpolate

from sklearn.decomposition import PCA
from pypeit import newspecobj
from pypeit import newspecobjs
#from pypeit import tracepca
from pypeit.core.pydl import spheregroup

from IPython import embed

# MASK VALUES FROM EXTRACTION
# 0 
# 2**0 = Flagged as bad detector pixel
# 2**1 = Flagged as affected by Cosmic Ray 
# 2**5 = Flagged as NAN (from something gone wrong)
# 2**6 = Entire region masked

mask_flags = dict(bad_pix=2**0, CR=2**1, NAN=2**5, bad_row=2**6)


def extract_asymbox2(image,left_in,right_in, ycen=None, weight_image=None):
    """ Extract the total flux within a variable window at many positions. This routine will accept an asymmetric/variable window
    specified by the left_in and right_in traces.  The ycen position is optional. If it is not provied, it is assumed to be integers
    in the spectral direction (as is typical for traces). Traces are expected to run vertically to be consistent with other
    extract_  routines. Based on idlspec2d/spec2d/extract_asymbox2.pro

    Args:
    image :  float ndarray
        Image to extract from. It is a 2-d array with shape (nspec, nspat)
    left  :  float ndarray
        Left boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.

    right  :  float ndarray
        Right boundary of region to be extracted (given as floating pt pixels). This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) forthe case of a single trace.


    Returns:
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
            ycen_out = np.arange(npix, dtype=int)
        elif ndim == 2:
            ycen_out = np.outer(np.ones(nTrace, dtype=int), np.arange(npix, dtype=int))
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

    if weight_image is not None:
        temp = np.array([weight_image[x1, y1] * image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2 = np.reshape(weight.flatten() * temp, (nTrace, npix, tempx))
        fextract = np.sum(temp2, axis=2)
        temp_wi = np.array([weight_image[x1, y1] for (x1, y1) in zip(bigy.flatten(), fullspot.flatten())])
        temp2_wi = np.reshape(weight.flatten() * temp_wi, (nTrace, npix, tempx))
        f_ivar = np.sum(temp2_wi, axis=2)
        fextract = fextract / (f_ivar + (f_ivar == 0)) * (f_ivar > 0)
    else:
        # Might be a more pythonic way to code this. I needed to switch the flattening order in order to get
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
    elif ((np.size(radius_in)==np.size(trace_in)) & (np.shape(radius_in) == np.shape(trace_in))):
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
            ycen_out = np.outer(np.ones(nTrace, dtype=int), np.arange(npix, dtype=int))
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




def extract_optimal(sciimg,ivar, mask, waveimg, skyimg, rn2_img, thismask, oprof, box_radius, specobj, min_frac_use = 0.05):

    """ Calculate the spatial FWHM from an object profile. Utility routine for fit_profile

    Parameters
    ----------
    sciimg : float ndarray shape (nspec, nspat)
       Science frame
    ivar: float ndarray shape (nspec, nspat)
       inverse variance of science frame. Can be a model or deduced from the image itself.
    mask: boolean ndarray
       mask indicating which pixels are good. Good pixels = True, Bad Pixels = False
    waveimg :  float ndarray
        Wavelength image. float 2-d array with shape (nspec, nspat)
    skyimg: float ndarray shape (nspec, nspat)
        Image containing our model of the sky
    rn2_img: float ndarray shape (nspec, nspat)
        Image containing the read noise squared (including digitization noise due to gain, i.e. this is an effective read noise)
    thismask: bool ndarray shape (nspec, nspat)
        Image indicating which pixels are on the slit/order in question. True=Good.
    oprof: float ndarray shape (nspec, nspat)
        Image containing the profile of the object that we are extracting
    box_radius: float
        Size of boxcar window in floating point pixels in the spatial direction.
    specobj: SpecObj object (from the SpecObj class in specobj.py).
         This is the container that holds object, trace,
         and extraction information for the object in question. This routine operates one object at a time.
    min_frac_use: float, optional, default = 0.05. If the sum of object profile arcoss the spatial direction
           are less than this value, the optimal extraction of this spectral pixel is masked because the majority of the
           object profile has been masked

    Returns
    -------
    Return value is None. The specobj object is changed in place with the boxcar and optimal dictionaries being filled
    with the extraction parameters.

    Revision History
    ----------------
    11-Mar-2005  Written by J. Hennawi and S. Burles.
    28-May-2018  Ported to python by J. Hennawi
    """
    imgminsky = sciimg - skyimg
    nspat = imgminsky.shape[1]
    nspec = imgminsky.shape[0]

    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    # TODO This makes no sense for difference imaging? Not sure we need NIVAR anyway
    var_no = np.abs(skyimg - np.sqrt(2.0) * np.sqrt(rn2_img)) + rn2_img

    ispec, ispat = np.where(oprof > 0.0)

    # Exit gracefully if we have no positive object profiles, since that means something was wrong with object fitting
    if not np.any(oprof > 0.0):
        msgs.warn('Object profile is zero everywhere. This aperture is junk.')
        junk = np.zeros(nspec)
        # Fill in the optimally extraction tags
        specobj.optimal['WAVE'] = junk
        specobj.optimal['COUNTS'] = junk
        specobj.optimal['COUNTS_IVAR'] = junk
        specobj.optimal['COUNTS_SIG'] = junk
        specobj.optimal['COUNTS_NIVAR'] = junk
        specobj.optimal['MASK'] = junk
        specobj.optimal['COUNTS_SKY'] = junk
        specobj.optimal['COUNTS_RN'] = junk
        specobj.optimal['FRAC_USE'] = junk
        specobj.optimal['CHI2'] = junk
        # Fill in the boxcar tags
        specobj.boxcar['WAVE'] = junk
        specobj.boxcar['COUNTS'] = junk
        specobj.boxcar['COUNTS_SIG'] = junk
        specobj.boxcar['COUNTS_IVAR'] = junk
        specobj.boxcar['COUNTS_NIVAR'] = junk
        specobj.boxcar['MASK'] = junk
        specobj.boxcar['COUNTS_SKY'] = junk
        specobj.boxcar['COUNTS_RN'] = junk
        specobj.boxcar['BOX_RADIUS'] = 0.0

        return None

    mincol = np.min(ispat)
    maxcol = np.max(ispat) + 1
    nsub = maxcol - mincol

    mask_sub = mask[:,mincol:maxcol]
    thismask_sub = thismask[:, mincol:maxcol]
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
    nivar_num = np.nansum(mask_sub*oprof_sub**2, axis=1) # Uses unit weights
    nvar_opt = ivar_denom*((mask_sub*vno_sub*oprof_sub**2).sum(axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    nivar_opt = 1.0/(nvar_opt + (nvar_opt == 0.0))
    # Optimally extract sky and (read noise)**2 in a similar way
    sky_opt = ivar_denom*(np.nansum(mask_sub*sky_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    rn2_opt = ivar_denom*(np.nansum(mask_sub*rn2_sub*oprof_sub**2, axis=1))/(nivar_num**2 + (nivar_num**2 == 0.0))
    rn_opt = np.sqrt(rn2_opt)
    rn_opt[np.isnan(rn_opt)]=0.0

    tot_weight = np.nansum(mask_sub*ivar_sub*oprof_sub, axis=1)
    prof_norm = np.nansum(oprof_sub, axis=1)
    frac_use = (prof_norm > 0.0)*np.nansum((mask_sub*ivar_sub > 0.0)*oprof_sub, axis=1)/(prof_norm + (prof_norm == 0.0))

    # Use the same weights = oprof^2*mivar for the wavelenghts as the flux.
    # Note that for the flux, one of the oprof factors cancels which does
    # not for the wavelengths.
    wave_opt = np.nansum(mask_sub*ivar_sub*wave_sub*oprof_sub**2, axis=1)/(mivar_num + (mivar_num == 0.0))
    mask_opt = (tot_weight > 0.0) & (frac_use > min_frac_use) & (mivar_num > 0.0) & (ivar_denom > 0.0) & \
               np.isfinite(wave_opt) & (wave_opt > 0.0)

    # Interpolate wavelengths over masked pixels
    badwvs = (mivar_num <= 0) | np.invert(np.isfinite(wave_opt)) | (wave_opt <= 0.0)
    if badwvs.any():
        oprof_smash = np.nansum(thismask_sub*oprof_sub**2, axis=1)
        # Can we use the profile average wavelengths instead?
        oprof_good = badwvs & (oprof_smash > 0.0)
        if oprof_good.any():
            wave_opt[oprof_good] = np.nansum(
                wave_sub[oprof_good,:]*thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)/\
                                   np.nansum(thismask_sub[oprof_good,:]*oprof_sub[oprof_good,:]**2, axis=1)
        oprof_bad = badwvs & ((oprof_smash <= 0.0) | (np.isfinite(oprof_smash) == False) | (wave_opt <= 0.0) | (np.isfinite(wave_opt) == False))
        if oprof_bad.any():
            # For pixels with completely bad profile values, interpolate from trace.
            f_wave = scipy.interpolate.RectBivariateSpline(spec_vec,spat_vec, waveimg*thismask)
            wave_opt[oprof_bad] = f_wave(specobj.trace_spec[oprof_bad], specobj.TRACE_SPAT[oprof_bad],grid=False)

    flux_model = np.outer(flux_opt,np.ones(nsub))*oprof_sub
    chi2_num = np.nansum((img_sub - flux_model)**2*ivar_sub*mask_sub,axis=1)
    chi2_denom = np.fmax(np.nansum(ivar_sub*mask_sub > 0.0, axis=1) - 1.0, 1.0)
    chi2 = chi2_num/chi2_denom

    # Fill in the optimally extraction tags
    '''
    specobj.optimal['WAVE'] = wave_opt    # Optimally extracted wavelengths
    specobj.optimal['COUNTS'] = flux_opt    # Optimally extracted flux
    specobj.optimal['COUNTS_IVAR'] = mivar_opt   # Inverse variance of optimally extracted flux using modelivar image
    specobj.optimal['COUNTS_SIG'] = np.sqrt(utils.inverse(mivar_opt))
    specobj.optimal['COUNTS_NIVAR'] = nivar_opt  # Optimally extracted noise variance (sky + read noise) only
    specobj.optimal['MASK'] = mask_opt    # Mask for optimally extracted flux
    specobj.optimal['COUNTS_SKY'] = sky_opt      # Optimally extracted sky
    specobj.optimal['COUNTS_RN'] = rn_opt        # Square root of optimally extracted read noise squared
    specobj.optimal['FRAC_USE'] = frac_use    # Fraction of pixels in the object profile subimage used for this extraction
    specobj.optimal['CHI2'] = chi2            # Reduced chi2 of the model fit for this spectral pixel
    '''
    specobj.OPT_WAVE = wave_opt    # Optimally extracted wavelengths
    specobj.OPT_COUNTS = flux_opt    # Optimally extracted flux
    specobj.OPT_COUNTS_IVAR = mivar_opt   # Inverse variance of optimally extracted flux using modelivar image
    specobj.OPT_COUNTS_SIG = np.sqrt(utils.inverse(mivar_opt))
    specobj.OPT_COUNTS_NIVAR = nivar_opt  # Optimally extracted noise variance (sky + read noise) only
    specobj.OPT_MASK = mask_opt    # Mask for optimally extracted flux
    specobj.OPT_COUNTS_SKY = sky_opt      # Optimally extracted sky
    specobj.OPT_COUNTS_RN = rn_opt        # Square root of optimally extracted read noise squared
    specobj.OPT_FRAC_USE = frac_use    # Fraction of pixels in the object profile subimage used for this extraction
    specobj.OPT_CHI2 = chi2            # Reduced chi2 of the model fit for this spectral pixel

    # Fill in the boxcar extraction tags
    flux_box  = extract_boxcar(imgminsky*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    # Denom is computed in case the trace goes off the edge of the image
    box_denom = extract_boxcar(waveimg*mask > 0.0, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    wave_box  = extract_boxcar(waveimg*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)/(box_denom + (box_denom == 0.0))
    varimg = 1.0/(ivar + (ivar == 0.0))
    var_box  = extract_boxcar(varimg*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    nvar_box  = extract_boxcar(var_no*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    sky_box  = extract_boxcar(skyimg*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    rn2_box  = extract_boxcar(rn2_img*mask, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    rn_posind = (rn2_box > 0.0)
    rn_box = np.zeros(rn2_box.shape,dtype=float)
    rn_box[rn_posind] = np.sqrt(rn2_box[rn_posind])
    pixtot  = extract_boxcar(ivar*0 + 1.0, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec)
    # If every pixel is masked then mask the boxcar extraction
    mask_box = (extract_boxcar(ivar*mask == 0.0, specobj.TRACE_SPAT,box_radius, ycen = specobj.trace_spec) != pixtot) & \
               np.isfinite(wave_box) & (wave_box > 0.0)

    bad_box = (wave_box <= 0.0) | np.invert(np.isfinite(wave_box)) | (box_denom == 0.0)
    # interpolate bad wavelengths over masked pixels
    if bad_box.any():
        f_wave = scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, waveimg)
        wave_box[bad_box] = f_wave(specobj.trace_spec[bad_box], specobj.TRACE_SPAT[bad_box],grid=False)

    ivar_box = 1.0/(var_box + (var_box == 0.0))
    nivar_box = 1.0/(nvar_box + (nvar_box == 0.0))

    specobj.BOX_WAVE = wave_box
    specobj.BOX_COUNTS = flux_box*mask_box
    specobj.BOX_COUNTS_IVAR = ivar_box*mask_box
    specobj.BOX_COUNTS_SIG = np.sqrt(utils.inverse(ivar_box*mask_box))
    specobj.BOX_COUNTS_NIVAR = nivar_box*mask_box
    specobj.BOX_MASK = mask_box
    specobj.BOX_COUNTS_SKY = sky_box
    specobj.BOX_COUNTS_RN = rn_box
    specobj.BOX_RADIUS = box_radius

    return None




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



def qa_fit_profile(x_tot,y_tot, model_tot, l_limit = None, r_limit = None, ind = None,
                   title =' ', xtrunc = 1e6, xlim = None, ylim = None, qafile = None):

    # Plotting pre-amble
    plt.close("all")
    #plt.clf()
#    plt.rc('text', usetex=True)
#    plt.rc('font', family='serif')
    width = 10.0 # Golden ratio 1.618
    fig, ax = plt.subplots(1, figsize=(width, width/1.618))

    if ind is None:
        indx = np.slice(x_tot.size)
    else:
        if len(ind) == 0:
            indx = np.slice(x_tot.size)
            title = title + ': no good pixels, showing all'
        else:
            indx = ind

    x = x_tot.flat[indx]
    y = y_tot.flat[indx]
    model = model_tot.flat[indx]

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
        ymax = np.fmax(1.5*model.max(), 0.1)
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
            ax.plot([plot_mid[i],plot_mid[i]], [y20[i],y80[i]], linewidth=1.2, color='orange')

    icl = nbin > 3
    if icl.any():
        ax.plot(plot_mid[icl],y50[icl],marker = 'o', color='lime', markersize=2, fillstyle='full', linestyle='None')
    else:
        ax.plot(plot_mid, y50, marker='o', color='lime', markersize=2, fillstyle = 'full', linestyle='None')

    isort = x.argsort()
    ax.plot(x[isort], model[isort], color='red', linewidth=1.0)



    if l_limit is not None:
        ax.axvline(x =l_limit, color='cornflowerblue',linewidth=2.0)
    if r_limit is not None:
        ax.axvline(x=r_limit, color='cornflowerblue',linewidth=2.0)

    ax.set_xlim(xlimit)
    ax.set_ylim(ylim)
    ax.set_title(title)
    ax.set_xlabel(r'$x/\sigma$')
    ax.set_ylabel('Normalized Profile')

    fig.subplots_adjust(left=0.15, right=0.9, top=0.9, bottom=0.15)
    plt.show()
    return


def return_gaussian(sigma_x, norm_obj, fwhm, med_sn2, obj_string, show_profile,
                    ind = None, l_limit = None, r_limit=None, xlim = None, xtrunc = 1e6):

    """
    Utility function to return Gaussian object profile in the case of low S/N ratio or too many rejected pixels.
    Args:
        sigma_x: ndarray (nspec, nspat)
            Spatial of gaussian
        norm_obj: ndarray (nspec, nspat)
            Normalized 2-d spectrum.
        fwhm: float
            FWHM parameter for Gaussian profile
        med_sn2: float
            Median (S/N)^2 used only for QA
        obj_string: str
            String identifying object. Used only for QA.
        show_profile: bool
            Is set, qa plot will be shown to screen
        ind: array, int
            Good indices in object profile. Used by QA routine.
        l_limit: float
            Left limit of profile fit where derivative is evaluated for Gaussian apodization. Used by QA routine.
        r_limit: float
            Right limit of profile fit where derivative is evaluated for Gaussian apodization. Used by QA routine.
        xlim: float
            Spatial location to trim object profile for plotting in QA routine.
        xtrunc: float
            Spatial nsigma to truncate object profile for plotting in QA routine.
    Returns:

    """

    profile_model = np.exp(-0.5*sigma_x**2)/np.sqrt(2.0 * np.pi)*(sigma_x ** 2 < 25.)
    info_string = "FWHM=" + "{:6.2f}".format(fwhm) + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2))
    title_string = obj_string + ', ' + info_string
    msgs.info(title_string)
    inf = np.isfinite(profile_model) == False
    ninf = np.sum(inf)
    if (ninf != 0):
        msgs.warn("Nan pixel values in object profile... setting them to zero")
        profile_model[inf] = 0.0
    # JFH Not sure we need this normalization.
    #nspat = profile_model.shape[1]
    #norm_spec = np.sum(profile_model, 1)
    #norm = np.outer(norm_spec, np.ones(nspat))
    #if np.sum(norm) > 0.0:
    #    profile_model = profile_model*(norm > 0.0)/(norm + (norm <= 0.0))
    if show_profile:
        qa_fit_profile(sigma_x, norm_obj, profile_model, title = title_string, l_limit = l_limit, r_limit = r_limit,
                       ind=ind, xlim = xlim, xtrunc=xtrunc)
    return profile_model


def fit_profile(image, ivar, waveimg, thismask, spat_img, trace_in, wave, flux, fluxivar,
                inmask=None, thisfwhm=4.0, max_trace_corr=2.0, sn_gauss=4.0, #, wvmnx = (2900.0,30000.0),
                maskwidth=None, prof_nsigma=None, no_deriv=False, gauss=False, obj_string='',
                show_profile=False):

    """Fit a non-parametric object profile to an object spectrum, unless the S/N ratio is low (> sn_gauss) in which
    fit a simple Gaussian. Port of IDL LOWREDUX long_gprofile.pro

     Parameters
     ----------
     image : numpy float 2-d array (nspec, nspat)
         sky-subtracted image
     ivar : numpy float 2-d array (nspec, nspat)
         inverse variance of sky-subtracted image
     waveimg numpy float 2-d array (nspec, nspat)
         2-d wavelength map
     spat_img: float ndarray, shape (nspec, nspat)
         Image containing the spatial location of pixels. If not input,
         it will be computed via spat_img = np.outer(np.ones(nspec), np.arange(nspat))
     trace_in : numpy 1-d array (nspec,)
         object trace
     wave : numpy 1-d array (nspec,)
         extracted wavelength of spectrum
     flux : numpy 1-d array (nspec,)
         extracted flux of spectrum
     fluxivar : numpy 1-d array (nspec,)
         inverse variance of extracted flux spectrum


    Optional Parameters
    ----------
    thisfwhm : float
         fwhm of the object trace
    max_trace_corr : float [default = 2.0]
         maximum trace correction to apply
    sn_gauss : float [default = 3.0]
         S/N ratio below which code just uses a Gaussian
    wvmnx : float [default = [2900.0,30000.0]
         wavelength range of usable part of spectrum
    maskwidth : float [default = None]
         object maskwidth determined from object finding algorithm. If = None,
         code defaults to use 3.0*(np.max(thisfwhm) + 1.0)
    prof_nsigma : float [default = None]
         Number of sigma to include in the profile fitting. This option is only needed for bright objects that are not
         point sources, which allows the profile fitting to fit the high S/N wings (rather than the default behavior
         which truncates exponentially). This allows for extracting all the flux and results in better sky-subtraction
         for bright extended objects.
    no_deriv : boolean [default = False]
         disables determination of derivatives and exponential apodization

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

    if inmask is None:
        inmask = (ivar > 0.0) & thismask

    totmask = inmask & (ivar > 0.0) & thismask

    if maskwidth is None: 3.0*(np.max(thisfwhm) + 1.0)
    if prof_nsigma is not None:
        no_deriv = True


    thisfwhm = np.fmax(thisfwhm,1.0) # require the FWHM to be greater than 1 pixel

    nspat = image.shape[1]
    nspec = image.shape[0]

    # dspat is the spatial position along the image centered on the object trace
    dspat = spat_img - np.outer(trace_in, np.ones(nspat))
    # create some images we will need
    sn2_img = np.zeros((nspec,nspat))
    spline_img = np.zeros((nspec,nspat))

    flux_sm = scipy.ndimage.filters.median_filter(flux, size=5, mode = 'reflect')
    fluxivar_sm0 =  scipy.ndimage.filters.median_filter(fluxivar, size = 5, mode = 'reflect')
    fluxivar_sm0 = fluxivar_sm0*(fluxivar > 0.0)
    wave_min = waveimg[thismask].min()
    wave_max = waveimg[thismask].max()

    # This adds an error floor to the fluxivar_sm, preventing too much rejection at high-S/N (i.e. standard stars)
    # TODO implement the ivar_cap function here in utils.
    sn_cap = 100.0
    fluxivar_sm = utils.clip_ivar(flux_sm,fluxivar_sm0, sn_cap)
    indsp = (wave >= wave_min) & (wave <= wave_max) & \
             np.isfinite(flux_sm) & \
             (flux_sm > -1000.0) & (fluxivar_sm > 0.0)
    b_answer, bmask   = pydl.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp],
                                     kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    b_answer, bmask2  = pydl.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask,
                                     kwargs_bspline={'everyn': 1.5}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    c_answer, cmask   = pydl.iterfit(wave[indsp], flux_sm[indsp], invvar = fluxivar_sm[indsp]*bmask2,
                                     kwargs_bspline={'everyn': 30}, kwargs_reject={'groupbadpix':True,'maxrej':1})
    spline_flux, _ = b_answer.value(wave[indsp])
    cont_flux, _ = c_answer.value(wave[indsp])

    sn2 = (np.fmax(spline_flux*(np.sqrt(np.fmax(fluxivar_sm[indsp], 0))*bmask2),0))**2
    ind_nonzero = (sn2 > 0)
    nonzero = np.sum(ind_nonzero)
    if(nonzero >0):
        (mean, med_sn2, stddev) = stats.sigma_clipped_stats(sn2,sigma_lower=3.0,sigma_upper=5.0)
    else:
        med_sn2 = 0.0

    #min_wave = np.min(wave[indsp])
    #max_wave = np.max(wave[indsp])
    # TODO -- JFH document this
    spline_flux1 = np.zeros(nspec)
    cont_flux1 = np.zeros(nspec)
    sn2_1 = np.zeros(nspec)
    ispline = (wave >= wave_min) & (wave <= wave_max)
    spline_tmp, _ = b_answer.value(wave[ispline])
    spline_flux1[ispline] = spline_tmp
    cont_tmp, _ = c_answer.value(wave[ispline])
    cont_flux1[ispline] = cont_tmp
    isrt = np.argsort(wave[indsp])
    s2_1_interp = scipy.interpolate.interp1d(wave[indsp][isrt], sn2[isrt],assume_sorted=False, bounds_error=False,fill_value = 0.0)
    sn2_1[ispline] = s2_1_interp(wave[ispline])
    bmask = np.zeros(nspec,dtype='bool')
    bmask[indsp] = bmask2
    spline_flux1 = pydl.djs_maskinterp(spline_flux1,(bmask == False))
    cmask2 = np.zeros(nspec,dtype='bool')
    cmask2[indsp] = cmask
    cont_flux1 = pydl.djs_maskinterp(cont_flux1,(cmask2 == False))
    (_, _, sigma1) = stats.sigma_clipped_stats(flux[indsp],sigma_lower=3.0,sigma_upper=5.0)

    sn2_med_filt = scipy.ndimage.filters.median_filter(sn2, size=9, mode='reflect')
    if np.any(totmask):
        sn2_interp = scipy.interpolate.interp1d(wave[indsp][isrt],sn2_med_filt[isrt],assume_sorted=False,
                                                bounds_error=False,fill_value = 'extrapolate')
        sn2_img[totmask] = sn2_interp(waveimg[totmask])
    else:
        msgs.warn('All pixels are masked')

    msgs.info('sqrt(med(S/N)^2) = ' + "{:5.2f}".format(np.sqrt(med_sn2)))

    # TODO -- JFH document this
    if(med_sn2 <= 2.0):
        spline_img[totmask]= np.fmax(sigma1,0)
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
        indbad2 = badpix & np.invert(goodval)
        nbad2 = np.sum(indbad2)
        ngood0 = np.sum(np.invert(badpix))
        if((nbad2 > 0) or (ngood0 > 0)):
            spline_flux1[indbad2] = np.median(spline_flux1[np.invert(badpix)])
        # take a 5-pixel median to filter out some hot pixels
        spline_flux1 = scipy.ndimage.filters.median_filter(spline_flux1,size=5,mode ='reflect')

        # Create the normalized object image
        if np.any(totmask):
            igd = (wave >= wave_min) & (wave <= wave_max)
            isrt1 = np.argsort(wave[igd])
            #plt.plot(wave[igd][isrt1], spline_flux1[igd][isrt1])
            #plt.show()
            spline_img_interp = scipy.interpolate.interp1d(wave[igd][isrt1],spline_flux1[igd][isrt1],assume_sorted=False,
                                                           bounds_error=False,fill_value = 'extrapolate')
            spline_img[totmask] = spline_img_interp(waveimg[totmask])
        else:
            spline_img[totmask] = np.fmax(sigma1, 0)

    # TODO -- JFH document this
    norm_obj = (spline_img != 0.0)*image/(spline_img + (spline_img == 0.0))
    norm_ivar = ivar*spline_img**2

    # Cap very large inverse variances
    ivar_mask = (norm_obj > -0.2) & (norm_obj < 0.7) & totmask & np.isfinite(norm_obj) & np.isfinite(norm_ivar)
    norm_ivar = norm_ivar*ivar_mask
    good = (norm_ivar.flatten() > 0.0)
    ngood = np.sum(good)

    xtemp = (np.cumsum(np.outer(4.0 + np.sqrt(np.fmax(sn2_1, 0.0)),np.ones(nspat)))).reshape((nspec,nspat))
    xtemp = xtemp/xtemp.max()

    sigma = np.full(nspec, thisfwhm/2.3548)
    fwhmfit = sigma*2.3548
    trace_corr = np.zeros(nspec)
    msgs.info("Gaussian vs b-spline of width " + "{:6.2f}".format(thisfwhm) + " pixels")
    area = 1.0
    # sigma_x represents the profile argument, i.e. (x-x0)/sigma
    sigma_x = (dspat/(np.outer(sigma, np.ones(nspat))) - np.outer(trace_corr, np.ones(nspat)))

    # If we have too few pixels to fit a profile or S/N is too low, just use a Gaussian profile
    if((ngood < 10) or (med_sn2 < sn_gauss**2) or (gauss is True)):
        msgs.info("Too few good pixels or S/N <" + "{:5.1f}".format(sn_gauss) + " or gauss flag set")
        msgs.info("Returning Gaussian profile")
        profile_model = return_gaussian(sigma_x, norm_obj, thisfwhm, med_sn2, obj_string,show_profile,ind=good,xtrunc=7.0)
        return (profile_model, trace_in, fwhmfit, med_sn2)

    mask = np.full(nspec*nspat, False, dtype=bool)

    # The following lines set the limits for the b-spline fit
    limit = scipy.special.erfcinv(0.1/np.sqrt(med_sn2))*np.sqrt(2.0)
    if(prof_nsigma is None):
        sinh_space = 0.25*np.log10(np.fmax((1000./np.sqrt(med_sn2)),10.))
        abs_sigma = np.fmin((np.abs(sigma_x.flat[good])).max(),2.0*limit)
        min_sigma = np.fmax(sigma_x.flat[good].min(), (-abs_sigma))
        max_sigma = np.fmin(sigma_x.flat[good].max(), (abs_sigma))
        nb = (np.arcsinh(abs_sigma)/sinh_space).astype(int) + 1
    else:
        msgs.info("Using prof_nsigma= " + "{:6.2f}".format(prof_nsigma) + " for extended/bright objects")
        nb = np.round(prof_nsigma > 10)
        max_sigma = prof_nsigma
        min_sigma = -1*prof_nsigma
        sinh_space = np.arcsinh(prof_nsigma)/nb

    rb = np.sinh((np.arange(nb) + 0.5) * sinh_space)
    bkpt = np.concatenate([(-rb)[::-1], rb])
    keep = ((bkpt >= min_sigma) & (bkpt <= max_sigma))
    bkpt = bkpt[keep]

    # Attempt B-spline first
    GOOD_PIX = (sn2_img > sn_gauss**2) & (norm_ivar > 0)
    IN_PIX   = (sigma_x >= min_sigma) & (sigma_x <= max_sigma) & (norm_ivar > 0)
    ngoodpix = np.sum(GOOD_PIX)
    ninpix     = np.sum(IN_PIX)

    if (ngoodpix >= 0.2*ninpix):
        inside,  = np.where((GOOD_PIX & IN_PIX).flatten())
    else:
        inside, = np.where(IN_PIX.flatten())


    si = inside[np.argsort(sigma_x.flat[inside])]
    sr = si[::-1]

    bset, bmask = pydl.iterfit(sigma_x.flat[si],norm_obj.flat[si], invvar = norm_ivar.flat[si],
                                   nord = 4, bkpt = bkpt, maxiter = 15, upper = 1, lower = 1)
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

    # If we have too few pixels after this step, then again just use a Gaussian profile and return.
    if(ninside < 10):
        msgs.info("Too few pixels inside l_limit and r_limit")
        msgs.info("Returning Gaussian profile")
        profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,show_profile,
                                                          ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
        return (profile_model, trace_in, fwhmfit, med_sn2)

    sigma_iter = 3
    isort = (xtemp.flat[si[inside]]).argsort()
    inside = si[inside[isort]]
    pb = np.ones(inside.size)

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

        mode_shift_out = utils.bspline_profile(xtemp.flat[inside], norm_obj.flat[inside], norm_ivar.flat[inside], profile_basis
                                      ,maxiter=1,kwargs_bspline= {'nbkpts':nbkpts})
        # Check to see if the mode fit failed, if so punt and return a Gaussian
        if not np.any(mode_shift_out[1]):
            msgs.info('B-spline fit to trace correction failed for fit to ninside = {:}'.format(ninside) + ' pixels')
            msgs.info("Returning Gaussian profile")
            profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                            show_profile, ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
            return (profile_model, trace_in, fwhmfit, med_sn2)


        mode_shift_set = mode_shift_out[0]
        temp_set = pydl.bspline(None, fullbkpt = mode_shift_set.breakpoints,nord=mode_shift_set.nord)
        temp_set.coeff = mode_shift_set.coeff[0, :]
        h0, _ = temp_set.value(xx)
        temp_set.coeff = mode_shift_set.coeff[1, :]
        h1, _ = temp_set.value(xx)
        ratio_10 = (h1/(h0 + (h0 == 0.0)))
        delta_trace_corr = ratio_10/(1.0 + np.abs(ratio_10)/0.1)
        trace_corr = trace_corr + delta_trace_corr

        profile_basis = np.column_stack((mode_zero,mode_stretch))
        mode_stretch_out = utils.bspline_profile(xtemp.flat[inside], norm_obj.flat[inside], norm_ivar.flat[inside], profile_basis,
                                            maxiter=1,fullbkpt = mode_shift_set.breakpoints)
        if not np.any(mode_stretch_out[1]):
            msgs.info('B-spline fit to width correction failed for fit to ninside = {:}'.format(ninside) + ' pixels')
            msgs.info("Returning Gaussian profile")
            profile_model  = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                                              show_profile,ind=good, l_limit=l_limit, r_limit=r_limit, xlim=7.0)
            return (profile_model, trace_in, fwhmfit, med_sn2)

        mode_stretch_set = mode_stretch_out[0]
        temp_set = pydl.bspline(None, fullbkpt = mode_stretch_set.breakpoints,nord=mode_stretch_set.nord)
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

        sigma_x = dspat/(np.outer(sigma, np.ones(nspat))) - np.outer(trace_corr, np.ones(nspat))

        # Update the profile B-spline fit for the next iteration
        if iiter < sigma_iter-1:
            ss = sigma_x.flat[inside].argsort()
            pb = (np.outer(area, np.ones(nspat,dtype=float))).flat[inside]
            keep = (bkpt >= sigma_x.flat[inside].min()) & (bkpt <= sigma_x.flat[inside].max())
            if keep.sum() == 0:
                keep = np.ones(bkpt.size, type=bool)
            bset_out = utils.bspline_profile(sigma_x.flat[inside[ss]],norm_obj.flat[inside[ss]],norm_ivar.flat[inside[ss]],pb[ss],
                                    nord = 4, bkpt=bkpt[keep],maxiter=2)
            if not np.any(bset_out[1]):
                msgs.info('B-spline to profile in trace and width correction loop failed for fit to ninside = {:}'.format(ninside) + ' pixels')
                msgs.info("Returning Gaussian profile")
                profile_model = return_gaussian(sigma_x, norm_obj, bspline_fwhm, med_sn2, obj_string,
                                                                  show_profile, ind=good, l_limit=l_limit,r_limit=r_limit, xlim=7.0)
                return (profile_model, trace_in, fwhmfit, med_sn2)

            bset = bset_out[0] # This updated bset used for the next set of trace corrections

    # Apply trace corrections only if they are small (added by JFH)
    if np.median(np.abs(trace_corr*sigma)) < max_trace_corr:
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
    bset_out = utils.bspline_profile(sigma_x.flat[ss[inside]],norm_obj.flat[ss[inside]], norm_ivar.flat[ss[inside]], pb.flat[ss[inside]],
                            nord=4, bkpt = bkpt, upper = 10, lower=10)
    bset = bset_out[0]
    outmask = bset_out[1]

    # igood = False for pixels within (min_sigma, max_sigma), True outside
    igood = (sigma_x.flatten() > min_sigma) & (sigma_x.flatten() < max_sigma)
    full_bsp = np.zeros(nspec*nspat, dtype=float)
    sigma_x_igood = sigma_x.flat[igood]
    yfit_out, _  = bset.value(sigma_x_igood)
    full_bsp[igood] = yfit_out
    isrt2 = sigma_x_igood.argsort()
    (peak, peak_x, lwhm, rwhm) = findfwhm(yfit_out[isrt2] - median_fit, sigma_x_igood[isrt2])


    left_bool = (((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] < peak_x)) | (sigma_x.flat[ss] < (peak_x-limit)))[::-1]
    ind_left, = np.where(left_bool)
    if ind_left.size == 0:
        lp = 0
    else:
        lp = np.fmax(ind_left.min(), 0)
    righ_bool = ((full_bsp[ss] < (min_level+median_fit)) & (sigma_x.flat[ss] > peak_x))  | (sigma_x.flat[ss] > (peak_x+limit))
    ind_righ, = np.where(righ_bool)
    if ind_righ.size == 0:
        rp = 0
    else:
        rp = np.fmax(ind_righ.min(), 0)

    l_limit = ((sigma_x.flat[ss])[::-1])[lp] - 0.1
    r_limit = sigma_x.flat[ss[rp]] + 0.1

    # Determine the left and right locations (l_limit and r_limit) where the profile logarithmic derivative crosses 1.0 for apodization
    # of the object profiles. If the profile slope is never actually one, then just find the maximum value in the interval between
    # (l_limit, -1.0) or (1.0, r_limit)
    l_lim_vec = np.arange(l_limit+0.1,-1.0, 0.1)
    l_lim_vec = np.asarray([-1.1]) if len(l_lim_vec) == 0 else l_lim_vec
    l_fit1, _ = bset.value(l_lim_vec)
    l_fit2, _ = bset.value(l_lim_vec*0.9)
    l_deriv_vec = (np.log(l_fit2) - np.log(l_fit1))/(0.1*l_lim_vec)
    l_deriv_max = np.fmax(l_deriv_vec.min(), -1.0)

    r_lim_vec = np.arange(r_limit-0.1,1.0, -0.1)
    r_lim_vec = np.asarray([1.1]) if len(r_lim_vec) == 0 else r_lim_vec

    r_fit1, _ = bset.value(r_lim_vec)
    r_fit2, _ = bset.value(r_lim_vec*0.9)
    r_deriv_vec = (np.log(r_fit2) - np.log(r_fit1))/(0.1*r_lim_vec)
    r_deriv_max = np.fmin(r_deriv_vec.max(), 1.0)


    while True:
        l_limit += 0.1
        l_fit, _ = bset.value(np.asarray([l_limit]))
        l2, _ = bset.value(np.asarray([l_limit])* 0.9)
        l_deriv = (np.log(l2[0]) - np.log(l_fit[0]))/(0.1*l_limit)
        if (l_deriv <= l_deriv_max) | (l_limit >= -1.0):
            break

    while True:
        r_limit -= 0.1
        r_fit, _ = bset.value(np.asarray([r_limit]))
        r2, _ = bset.value(np.asarray([r_limit])* 0.9)
        r_deriv = (np.log(r2[0]) - np.log(r_fit[0]))/(0.1*r_limit)
        if (r_deriv >= r_deriv_max) | (r_limit <= 1.0):
            break

    # JXP kludge
    if prof_nsigma is not None:
       #By setting them to zero we ensure QA won't plot them in the profile QA.
       l_limit = 0.0
       r_limit = 0.0

    # Apodization of object profiles with exponential, ensuring continuity of first derivative
    if (l_deriv < 0) and (r_deriv > 0) and no_deriv is False:
        left = sigma_x.flatten() < l_limit
        full_bsp[left] =  np.exp(-(sigma_x.flat[left]-l_limit)*l_deriv) * l_fit[0]
        right = sigma_x.flatten() > r_limit
        full_bsp[right] = np.exp(-(sigma_x.flat[right] - r_limit) * r_deriv) * r_fit[0]

    # Final object profile
    full_bsp = full_bsp.reshape(nspec,nspat)
    profile_model = full_bsp*pb
    res_mode = (norm_obj.flat[ss[inside]] - profile_model.flat[ss[inside]])*np.sqrt(norm_ivar.flat[ss[inside]])
    chi_good = (outmask == True) & (norm_ivar.flat[ss[inside]] > 0)
    chi_med = np.median(res_mode[chi_good]**2)
    chi_zero = np.median(norm_obj.flat[ss[inside]]**2*norm_ivar.flat[ss[inside]])

    msgs.info("--------------------  Results of Profile Fit --------------------")
    msgs.info(" min(fwhmfit)={:5.2f}".format(fwhmfit.min()) +
              " max(fwhmfit)={:5.2f}".format(fwhmfit.max()) + " median(chi^2)={:5.2f}".format(chi_med) +
              " nbkpts={:2d}".format(bkpt.size))
    msgs.info("-----------------------------------------------------------------")

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

    info_string = "FWHM range:" + "{:5.2f}".format(fwhmfit.min()) + " - {:5.2f}".format(fwhmfit.max()) \
                  + ", S/N=" + "{:8.3f}".format(np.sqrt(med_sn2)) + ", median(chi^2)={:8.3f}".format(chi_med)
    title_string = obj_string + ' ' + info_string
    if(show_profile):
        qa_fit_profile(sigma_x, norm_obj/(pb + (pb == 0.0)), full_bsp,
                       l_limit = l_limit, r_limit = r_limit, ind = ss[inside], xlim = prof_nsigma, title = title_string)

    return (profile_model, xnew, fwhmfit, med_sn2)


def parse_hand_dict(hand_extract_dict):
    """ Utility routine for objfind to parase the hand_extract_dict dictionary for hand selected apertures

    Parameters
    ----------
    hand_extract_dict:   dictionary

    Returns
    -------
    hand_spec:  spectral pixel location, numpy float 1-d array with size equal to number of hand aperatures requested

    hand_spat:  spatial pixel location, numpy float 1-d array with size equal to number of hand aperatures requested

    hand_det:   Detector for hand apertures. This should either be ba numpy float 1-d array with size equal to number of hand
                apertures requested, or a single number which applies to all the hand apertures provied by hand_spec, hand_spat

    hand_fwhm:  hand aperture fwhm for extraction. This should either be ba numpy float 1-d array with size equal to number of hand
                apertures requested, or a single number which applies to all the hand apertures provied by hand_spec, hand_spat

    Revision History
    ----------------
    23-June-2018  Written by J. Hennawi
    """


    if ('hand_extract_spec' not in hand_extract_dict.keys() | 'hand_extract_spat' not in hand_extract_dict.keys()):
        raise ValueError('hand_extract_spec and hand_extract_spat must be set in the hand_extract_dict')

    hand_extract_spec=np.asarray(hand_extract_dict['hand_extract_spec'])
    hand_extract_spat=np.asarray(hand_extract_dict['hand_extract_spat'])
    hand_extract_det = np.asarray(hand_extract_dict['hand_extract_det'])
    if(hand_extract_spec.size == hand_extract_spat.size == hand_extract_det.size) == False:
        raise ValueError('hand_extract_spec, hand_extract_spat, and hand_extract_det must have the same size in the hand_extract_dict')
    nhand = hand_extract_spec.size

    hand_extract_fwhm = np.asarray(hand_extract_dict['hand_extract_fwhm'])
    # hand_extract_fwhm = hand_extract_dict.get('hand_extract_fwhm')
    # if hand_extract_fwhm is not None:
    #     hand_extract_fwhm = np.asarray(hand_extract_fwhm)
    #     if(hand_extract_fwhm.size==hand_extract_spec.size):
    #         pass
    #     elif (hand_extract_fwhm.size == 1):
    #         hand_extract_fwhm = np.full(nhand, hand_extract_fwhm)
    #     else:
    #         raise ValueError('hand_extract_fwhm must either be a number of have the same size as hand_extract_spec and hand_extract_spat')
    # else:
    #     hand_extract_fwhm = np.full(nhand, None)
    return hand_extract_spec, hand_extract_spat, hand_extract_det, hand_extract_fwhm



def iter_tracefit(image, xinit_in, ncoeff, inmask = None, trc_inmask = None, fwhm = 3.0, maxdev = 2.0, maxiter = 25,
                  niter=9, gweight=False, show_fits=False, idx = None, verbose=False, xmin= None, xmax = None):
    """ Utility routine for object find to iteratively trace and fit. Used by both objfind and ech_objfind

    Parameters
    ----------
    image: ndaarray, float
        Image of objects to be traced
    xinit_in: ndarray, float
        Initial guesses for spatial direction trace. This can either be an 2-d  array with shape
        (nspec, nTrace) array, or a 1-d array with shape (nspec) for the case of a single trace.
    ncoeff: int
        Order of polynomial fits to trace

    Optional Parameter
    ------------------
    inmask: ndarray, bool
        Input mask for the image
    trc_inmask: ndarray, bool
        Input mask for the trace, i.e. places where you know the trace is going to be bad that you always want to mask in the
        fits. Same size as xinit_in (nspec, nTrace)
    fwhm: float
        fwhm width parameter for the flux or gaussian weighted tracing. For flux weighted trace the code does a third
        of the iterations with window 1.3*fwhm, a third with 1.1*fwhm, and a third with fwhm. For Gaussian weighted tracing
        it uses the fwhm to determine the sigma of the Gausisan which is used for all iterations.
    gweight: bool, default = False
        If gweight is True the routine does Gaussian weighted tracing, if it is False it will do flux weighted tracing.
        Normally the approach is to do a round of flux weighted tracing first, and then refine the traces with Gaussian
        weighted tracing.
    show_fits: bool, default = False
        Will plot the data and the fits.
    idx: ndarray of strings, default = None
        Array of idx IDs for each object. Used only if show_fits is true for the plotting.
    xmin: float, default = None
        Lower reference for robust_polyfit polynomial fitting. Default is to use zero
    xmax: float, defualt = None
        Upper refrence for robust_polyfit polynomial fitting.  Default is to use the image size in nspec direction
    Returns
    -------
    xpos: ndarray, float
       The output has the same size as xinit_in and contains the fit to the spatial direction of trace for each
       object.



    Revision History
    ----------------
    23-June-2018  Written by J. Hennawi
    """


    if inmask is None:
        inmask = np.ones_like(image,dtype=bool)

    # Allow for single vectors as input as well:
    nspec = xinit_in.shape[0]

    if xmin is None:
        xmin = 0.0
    if xmax is None:
        xmax = float(nspec-1)

    # Deal with the possibility of vectors as inputs instead of 2d arrays
    if xinit_in.ndim == 1:
        nobj = 1
        xinit = xinit_in.reshape(nspec,1)
        if trc_inmask is not None:
            trc_inmask_out = trc_inmask.reshape(nspec,1)
        else:
            trc_inmask_out = np.ones_like(xinit,dtype=bool)
    else:
        nobj = xinit_in.shape[1]
        xinit = xinit_in
        if trc_inmask is not None:
            trc_inmask_out = trc_inmask
        else:
            trc_inmask_out = np.ones_like(xinit,dtype=bool)

    spec_vec = np.arange(nspec)

    if verbose:
        msgs.info('Fitting the object traces')

    # Iterate flux weighted centroiding
    fwhm_vec = np.zeros(niter)
    fwhm_vec[0:niter//3] = 1.3*fwhm
    fwhm_vec[niter//3:2*niter//3] = 1.1*fwhm
    fwhm_vec[2*niter//3:] = fwhm

    if gweight:
        title_text = 'Gaussian Weighted'
    else:
        title_text = 'Flux Weighted'

    xfit1 = np.copy(xinit)

    for iiter in range(niter):
        if gweight:
            xpos1, xerr1 = trace_slits.trace_gweight(image*inmask,xfit1, invvar=inmask.astype(float),sigma=fwhm/2.3548)
        else:
            xpos1, xerr1 = trace_slits.trace_fweight(image*inmask,xfit1, invvar = inmask.astype(float), radius = fwhm_vec[iiter])

        # If a trc_inmask was input, always set the masked values to the initial input crutch. The point is that the crutch
        # initially comes from either the standard or the slit boundaries, and if we continually replace it for all iterations
        # we will naturally always extraplate the trace to match the shape of a high S/N ratio fit (i.e. either the standard)
        # or the flat which was used to determine the slit edges.
        xpos1[np.invert(trc_inmask_out)] = xinit[np.invert(trc_inmask_out)]

        # Do not do any kind of masking based on xerr1. Trace fitting is much more robust when masked pixels are simply
        # replaced by the tracing crutch. We thus do not do weighted fits, i.e. uniform weights, but we set the relative
        # weight of the trc_inmask pixels to be lower. This way they still get a say but do not overly influence the fit.
        xinvvar = np.ones_like(xpos1.T)
        xinvvar[np.invert(trc_inmask_out.T)] = 0.1
        pos_set1 = pydl.xy2traceset(np.outer(np.ones(nobj),spec_vec), xpos1.T,
                                    #inmask = trc_inmask_out.T,
                                    ncoeff=ncoeff, maxdev=maxdev,
                                    maxiter=maxiter, invvar=xinvvar, xmin=xmin, xmax =xmax)
        xfit1 = pos_set1.yfit.T
        # bad pixels have errors set to 999 and are returned to lie on the input trace. Use this only for plotting below
        #errmask = (xerr1 > 990.0)  # bad pixels have errors set to 999 and are returned to lie on the input trace
        outmask = pos_set1.outmask.T
        # Plot all the points that were not masked initially
        if(show_fits) & (iiter == niter - 1):
            for iobj in range(nobj):
                # The sum of all these masks adds up to the number of pixels.
                inmask_trc = np.invert(trc_inmask_out[:,iobj]) # masked on the way in
                errmask = xerr1[:,iobj] > 990.0 # masked by fweight or gweight, was set to input trace and still fit
                rejmask = np.invert(outmask[:, iobj]) & np.invert(inmask_trc) # was good on the way in, masked by the poly fit
                nomask = outmask[:, iobj] & np.invert(errmask) # was actually fit and not replaced to input trace
                plt.plot(spec_vec[nomask],xpos1[nomask,iobj],marker='o', c='k', markersize=3.0,linestyle='None',label=title_text + ' Centroid')
                plt.plot(spec_vec,xinit[:,iobj],c='g', zorder = 25, linewidth=2.0,linestyle='--', label='initial guess')
                plt.plot(spec_vec,xfit1[:,iobj],c='red',zorder=30,linewidth = 2.0, label ='fit to trace')
                if np.any(errmask):
                    plt.plot(spec_vec[errmask],xfit1[errmask,iobj], c='blue',marker='+',
                             markersize=5.0,linestyle='None',zorder= 20, label='masked by tracing, set to init guess')
                if np.any(rejmask):
                    plt.plot(spec_vec[rejmask],xpos1[rejmask,iobj], c='cyan',marker='v',
                             markersize=5.0,linestyle='None',zorder= 20, label='masked by polynomial fit')
                if np.any(inmask_trc):
                    plt.plot(spec_vec[inmask_trc],xpos1[inmask_trc,iobj],
                             c='orange',marker='s',markersize=3.0,linestyle='None',zorder= 20, label='input masked points, not fit')
                try:
                    plt.title(title_text + ' Centroid to object {:s}.'.format(idx[iobj]))
                except TypeError:
                    plt.title(title_text + ' Centroid to object {:d}.'.format(iobj))
                plt.ylim((0.995*xfit1[:, iobj].min(), 1.005*xfit1[:, iobj].max()))
                plt.xlabel('Spectral Pixel')
                plt.ylabel('Spatial Pixel')
                plt.legend()
                plt.show()

    # Returns the fit, the actual weighted traces, and the pos_set1 object
    return xfit1, xpos1, xerr1, pos_set1


def create_skymask_fwhm(sobjs, thismask):
    """
    Creates a skymask from a SpecObjs object using the fwhm of each object
    Args:
        sobjs: SpecObjs object
            Objects for which you would like to create the mask
        thismask: ndarray, bool, shape (nspec, nspat)
            Boolean image indicating pixels which are on the slit

    Returns:
        skymask:  ndarray, bool, shape (nspec, nspat)
           Boolean image with the same size as thismask indicating which pixels are usable for global sky subtraction.
           True = usable for sky subtraction, False = should be masked when sky subtracting.

    """
    nobj = len(sobjs)
    skymask = np.copy(thismask)
    if nobj == 0:
        return skymask
    else:
        nspec, nspat = thismask.shape
        # spatial position everywhere along image
        spat_img = np.outer(np.ones(nspec, dtype=int),np.arange(nspat, dtype=int))
        all_fwhm = sobjs.FWHM
        med_fwhm = np.median(all_fwhm)
        for iobj in range(nobj):
            # Create a mask for the pixels that will contribue to the arc
            slit_img = np.outer(sobjs[iobj].TRACE_SPAT, np.ones(nspat))  # central trace replicated spatially
            objmask_now = thismask & (spat_img > (slit_img - med_fwhm)) & (spat_img < (slit_img + med_fwhm))
            skymask = skymask & np.invert(objmask_now)

        # Check that we have not performed too much masking
        if (np.sum(skymask)/np.sum(thismask) < 0.10):
            msgs.warn('More than 90% of  usable area on this slit would be masked and not used by global sky subtraction. '
                      'Something is probably wrong with object finding for this slit. Not masking object for global sky subtraction.')
            return np.copy(thismask)
        else:
            return skymask

def objfind(image, thismask, slit_left, slit_righ, inmask=None, fwhm=3.0, maxdev=2.0, ir_redux=False, spec_min_max=None,
            hand_extract_dict=None, std_trace=None, extrap_npoly=3, ncoeff=5, nperslit=None, bg_smth=5.0,
            extract_maskwidth=4.0, sig_thresh=10.0, peak_thresh=0.0, abs_thresh=0.0, trim_edg=(5,5),
            skymask_nthresh=1.0, specobj_dict=None, cont_fit=True, npoly_cont=1,
            show_peaks=False, show_fits=False, show_trace=False, show_cont=False, debug_all=False, qa_title='objfind'):

    """ Find the location of objects in a slitmask slit or a echelle order.


    Args:
    image :  float ndarray
        Image to search for objects from. This image has shape (nspec, nspat) image.shape where the first dimension (nspec)
        is spectral, and second dimension (nspat) is spatial. Note this image can either have the sky background in it, or have already been sky subtracted.
        Object finding works best on sky-subtracted images, but often one runs on the frame with sky first to identify the brightest
        objects which are then masked (see skymask below) in sky subtraction.
    thismask:  boolean ndarray
        Boolean mask image specifying the pixels which lie on the slit/order to search for objects on.
        The convention is: True = on the slit/order, False  = off the slit/order

    slit_left:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)
    slit_righ:  float ndarray
        Left boundary of slit/order to be extracted (given as floating pt pixels). This a 1-d array with shape (nspec, 1)
        or (nspec)
    inmask: float ndarray, default = None
        Input mask image.
    spec_min_max: tuple of float or int, (2), default=None.
        This is tuple of two elements which defines the minimum and maximum of the slit in the
        spectral direction on the detector. If not passed in it will be determined automatically from the thismask
    fwhm: float, default = 3.0
        Estimated fwhm of the objects in pixels
    maxdev (float): default=2.0
        Maximum deviation of pixels from polynomial fit to trace used to reject bad pixels in trace fitting.
    hand_extract_dict: dict, default = None
        Dictionary containing information about apertures requested by user that should be place by hand in the object list.
        This option is useful for cases like an emission line obect that the code fails to find with its significance threshold
    std_trace: float ndarray, shape = (nspec,), default = None
        This is a one dimensioal array containing the standard star trace which is used as a crutch for tracing. If the
        no standard star is provided the code uses the the slit boundaries as the crutch.
    ncoeff: int, default = 5
        Order of legendre polynomial fits to the trace
    nperslit: int, default = 10
        Maximum number of objects allowed per slit. The code will take the nperslit most significant detections.
    bg_smth: float, default = 5.0
        Size of the smoothing kernel in units of fwhm used to determine the background level from the smash of the image
        along the curved traces. This background subtracted smashed image is used for peak finding to identify objects
    extract_maskwidth: float, default = 3.0
        This parameter determines the initial size of the region in units of fwhm that will be used for local sky subtraction in the routine
        skysub.local_skysub_extract.
    sig_thresh: float, default = 5.0
        Significance threshold for object detection. The code uses the maximum of the thresholds defined by sig_thresh,
        peak_thresh, and abs_thresh.
         For the default behavior peak_thresh and abs_thresh are zero, so sig_thresh defines the threshold.
    peak_thresh: float, default = 0.0
        Peak threshold for object detection. This is a number between 0 and 1 and represents the fraction of the brightest
        object on the slit that will be kept as an object, i.e. if ymax is the brightest object of the spectrum smashed
        out in the spectral direction, all objects with ypeak > peak_thresh*ymak are kept. The code uses the maximum of the
        thresholds defined by sig_thresh, peak_thers, and abs_thresh.
    abs_thresh: float, defalt = 0.0.
        Absolute threshold for object detection.  Objects are found by smashing out the spectral direction along the curved
        slit/order traces, and abs_thresh is in the units of this smashed profile.  The code uses the maximum of the
        thresholds defined by sig_thresh, peak_thers, and abs_thresh.
    trim_edg: tuple of integers or float, default = (3,3)
        Ignore objects within this many pixels of the left and right slit boundaries, where the first element refers to the left
        and second refers to the right.
    skymask_nthresh: float, default = 2.0
        The multiple of the final object finding threshold (see above) which is used to create the skymask using the value
        of the peak flux in the slit profile (image with the spectral direction smashed out).
    cont_fit (bool): default=True:
        Fit a continuum to the illumination pattern across the slit when peak finding
    npoly_cont (int): default=1
        Order of polynomial fit to the illumination pattern across the slit when peak finding
    specobj_dict: dict, default = None
         Dictionary containing meta-data for the objects that will be propgated into the SpecObj objects, i.e. setup,
         slitid, detector, object type, and pipeline. The default is None, in which case the following dictionary will be used.
         specobj_dict = {'setup': None, 'slitid': 999, 'det': 1, 'objtype': 'unknown', 'pypeline': 'unknown'}

    Returns:
    sobjs:   SpecoObjs object
         Object containing the information about the objects found on the slit/order
    np.ndarray: Skymask image



    Revision History
    ----------------
    10-Mar-2005  First version written by D. Schlegel, LBL
    2005-2018    Improved by J. F. Hennawi and J. X. Prochaska
    23-June-2018 Ported to python by J. F. Hennawi and significantly improved
    """

    #debug_all=True
    if debug_all:
        show_peaks=True
        #show_fits = True
        show_trace = True
        show_cont = True

    if specobj_dict is None:
        specobj_dict = dict(setup=None, slitid=999, det=1, objtype='unknown', pypeline='unknown', orderindx=999)

    # Check that peak_thresh values make sense
    if ((peak_thresh >=0.0) & (peak_thresh <=1.0)) == False:
        msgs.error('Invalid value of peak_thresh. It must be between 0.0 and 1.0')

    frameshape = image.shape
    nspec = frameshape[0]
    nspat = frameshape[1]
    specmid = nspec//2

    # Some information about this slit we need for later when we instantiate specobj objects
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)
    slit_spec_pos = nspec/2.0

    slit_spat_pos = (np.interp(slit_spec_pos, spec_vec, slit_left), np.interp(slit_spec_pos, spec_vec, slit_righ))

    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)

    # If a mask was not passed in, create it
    if inmask is None:
        inmask = thismask
    # If spec_min_max was not passed in, determine it from the thismask
    if spec_min_max is None:
        ispec, ispat = np.where(thismask)
        spec_min_max = (ispec.min(), ispec.max())


    totmask = thismask & inmask & np.invert(edgmask)
    thisimg =image*totmask
    #  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    #nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())
    # Mask skypixels with 2 fwhm of edge
    left_asym = np.outer(slit_left,np.ones(int(nsamp))) + np.outer(xsize/nsamp, np.arange(nsamp))
    righ_asym = left_asym + np.outer(xsize/nsamp, np.ones(int(nsamp)))
    # This extract_asymbox2 call smashes the image in the spectral direction along the curved object traces
    # TODO Should we be passing the mask here with extract_asymbox or not?
    flux_spec = extract_asymbox2(thisimg, left_asym, righ_asym, weight_image=totmask.astype(float))
    mask_spec = extract_asymbox2(totmask, left_asym, righ_asym, weight_image=totmask.astype(float)) < 0.3
    flux_mean, flux_median, flux_sig = stats.sigma_clipped_stats(flux_spec, mask = mask_spec, axis=0, sigma = 3.0,
                                                           cenfunc='median', stdfunc=utils.nan_mad_std)
##   New CODE
    # 1st iteration
    smash_mask = np.isfinite(flux_mean)
    flux_mean_med0 = np.median(flux_mean[smash_mask])
    flux_mean[np.invert(smash_mask)] = flux_mean_med0
    fluxsub0 = flux_mean - flux_mean_med0
    fluxconv0 = scipy.ndimage.filters.gaussian_filter1d(fluxsub0, fwhm/2.3548, mode='nearest')

    cont_samp = np.fmin(int(np.ceil(nsamp/(fwhm/2.3548))), 30)
    cont, cont_mask0 = arc.iter_continuum(
        fluxconv0, inmask=smash_mask, fwhm=fwhm,cont_frac_fwhm=2.0, sigthresh=2.0, sigrej=2.0, cont_samp=cont_samp,
        npoly=(0 if (nsamp/fwhm < 20.0) else npoly_cont), cont_mask_neg=ir_redux, debug=show_cont,
        qa_title='Smash Image Background, 1st iteration: Slit# {:d}'.format(specobj_dict['slitid']))

    # Second iteration
    flux_mean_med = np.median(flux_mean[cont_mask0])
    fluxsub = flux_mean - flux_mean_med
    fluxconv = scipy.ndimage.filters.gaussian_filter1d(fluxsub, fwhm/2.3548, mode='nearest')

    cont, cont_mask = arc.iter_continuum(
        fluxconv, inmask=smash_mask, fwhm=fwhm, cont_frac_fwhm=2.0, sigthresh=2.0, sigrej=2.0, cont_samp=cont_samp,
        npoly=(0 if (nsamp/fwhm < 20.0) else npoly_cont), cont_mask_neg=ir_redux, debug=show_cont,
        qa_title='Smash Image Background: 2nd iteration: Slit# {:d}'.format(specobj_dict['slitid']))
    fluxconv_cont = (fluxconv - cont) if cont_fit else fluxconv
    # JFH TODO Do we need a running median as was done in the OLD code? Maybe needed for long slits. We could use
    #  use the cont_mask to isolate continuum pixels, and then interpolate the unmasked pixels.

## New CODE

##   OLD CODE
#    smash_mask = np.isfinite(flux_mean)
#    flux_mean_med = np.median(flux_mean[smash_mask])
#    flux_mean[np.invert(smash_mask)] = 0.0
#    if (nsamp < 3.0*bg_smth*fwhm):
#        # This may lead to many negative fluxsub values..
#        # TODO: Calculate flux_mean_med by avoiding the peak
#        fluxsub = flux_mean - flux_mean_med
#    else:
#        kernel_size= int(np.ceil(bg_smth*fwhm) // 2 * 2 + 1) # This ensure kernel_size is odd
#        # TODO should we be using  scipy.ndimage.filters.median_filter to better control the boundaries?
#        fluxsub = flux_mean - scipy.signal.medfilt(flux_mean, kernel_size=kernel_size)
#        # This little bit below deals with degenerate cases for which the slit gets brighter toward the edge, i.e. when
#        # alignment stars saturate and bleed over into other slits. In this case the median smoothed profile is the nearly
#        # everywhere the same as the profile itself, and fluxsub is full of zeros (bad!). If 90% or more of fluxsub is zero,
#        # default to use the unfiltered case
#        isub_bad = (fluxsub == 0.0)
#        frac_bad = np.sum(isub_bad)/nsamp
#        if frac_bad > 0.9:
#            fluxsub = flux_mean - flux_mean_med
#
#    fluxconv = scipy.ndimage.filters.gaussian_filter1d(fluxsub, fwhm/2.3548, mode='nearest')
#
#    cont_samp = np.fmin(int(np.ceil(nsamp/(fwhm/2.3548))), 30)
#    cont, cont_mask = arc.iter_continuum(fluxconv, inmask=smash_mask, fwhm=fwhm,
#                                         cont_frac_fwhm=2.0, sigthresh=2.0,
#                                         sigrej=2.0, cont_samp=cont_samp,
#                                         npoly=(0 if (nsamp/fwhm < 20.0) else npoly_cont),
#                                         cont_mask_neg=ir_redux, debug=debug_all)
#    fluxconv_cont = (fluxconv - cont) if cont_fit else fluxconv
## OLD CODE

    if not np.any(cont_mask):
        cont_mask = np.ones(int(nsamp),dtype=bool) # if all pixels are masked for some reason, don't mask

    (mean, med, skythresh) = stats.sigma_clipped_stats(fluxconv_cont[cont_mask], sigma=1.5)
    (mean, med, sigma)     = stats.sigma_clipped_stats(fluxconv_cont[cont_mask], sigma=2.5)

    if(skythresh == 0.0) & (sigma != 0.0):
        skythresh = sigma
    elif(skythresh == 0.0) & (sigma == 0.0):  # if both SKYTHRESH and sigma are zero mask out the zero pixels and reavaluate
        good = fluxconv_cont > 0.0
        if np.any(good) == True:
            (mean, med_sn2, skythresh) = stats.sigma_clipped_stats(fluxconv_cont[good], sigma=1.5)
            (mean, med_sn2, sigma) = stats.sigma_clipped_stats(fluxconv_cont[good], sigma=2.5)
        else:
            msgs.error('Object finding failed. All the elements of the fluxconv_cont spatial profile array are zero')
    else:
        pass

    # Now find all the peaks without setting any threshold
    ypeak, _, xcen, sigma_pk, _, _, _, _ = arc.detect_lines(fluxconv_cont, cont_subtract = False, fwhm = fwhm,
                                                            input_thresh = 'None', debug=False)

    # Get rid of peaks within trim_edg of slit edge which are almost always spurious, this should have been handled
    # with the edgemask, but we do it here anyway
    not_near_edge = (xcen > trim_edg[0]) & (xcen < (nsamp - trim_edg[1]))
    if np.any(np.invert(not_near_edge)):
        msgs.warn('Discarding {:d}'.format(np.sum(np.invert(not_near_edge))) +
                  ' at spatial pixels spat = {:}'.format(xcen[np.invert(not_near_edge)]) +
                  ' which land within trim_edg = (left, right) = {:}'.format(trim_edg) +
                  ' pixels from the slit boundary for this nsamp = {:5.2f}'.format(nsamp) + ' wide slit')
        msgs.warn('You must decrease from the current value of trim_edg in order to keep them')
        msgs.warn('Such edge objects are often spurious')

    xcen = xcen[not_near_edge]
    ypeak = ypeak[not_near_edge]

    # If the user requested the nperslit most significant peaks have been requested, then grab and return only these lines
    if nperslit is not None:
        #if nperslit > len(ypeak):
        #    msgs.warn('Requested nperslit = {:}'.format(nperslit) + ' most significant objects but only npeak = {:}'.format(len(xcen)) +
        #              ' were found. Returning all the objects found.')
        #else:
        ikeep = (ypeak.argsort()[::-1])[0:nperslit]
        xcen = xcen[ikeep]
        ypeak = ypeak[ikeep]

    npeak = len(xcen)
    # Instantiate a null specobj
    sobjs = newspecobjs.SpecObjs()
    # Choose which ones to keep and discard based on threshold params. Create SpecObj objects

    # Possible thresholds    [significance,  fraction of brightest, absolute]
    thresh_peak = peak_thresh * ypeak.max() if len(ypeak) > 0 else 0.0
    threshvec = np.array([sig_thresh * sigma, thresh_peak, abs_thresh])
    threshold = threshvec.max()
    if npeak > 0:
        if threshvec.argmax() == 0:
            msgs.info('Used SIGNIFICANCE threshold: sig_thresh = {:3.1f}'.format(sig_thresh) +
                      ' * sigma = {:5.2f}'.format(sigma))
        elif threshvec.argmax() == 1:
            msgs.info('Used FRACTION of BRIGHTEST threshold: peak_thresh = {:3.1f}'.format(peak_thresh) +
                      ' * ypeak_max = {:5.2f}'.format(ypeak.max()))
        elif threshvec.argmax() == 2:
            msgs.info('Used ABSOLUTE threshold of abs_thresh = {:5.2f}'.format(abs_thresh))
        msgs.info('Object finding threshold of: {:5.2f}'.format(threshold))
        # Trim to only objects above this threshold
        ikeep = (ypeak >= threshold)
        xcen = xcen[ikeep]
        ypeak = ypeak[ikeep]
        nobj_reg = len(xcen)
        # Now create SpecObj objects for all of these
        for iobj in range(nobj_reg):
            # ToDo Label with objid and objind here?
            if specobj_dict['pypeline'] == 'MultiSlit':
                thisobj = newspecobj.SpecObj('MultiSlit', specobj_dict['det'],
                                             (slit_left[specmid], slit_righ[specmid]),
                                             slitid=specobj_dict['slitid'],
                                             objtype=specobj_dict['objtype'])
            else:
                embed('1711 of extract')
            thisobj.spat_fracpos = xcen[iobj]/nsamp
            thisobj.smash_peakflux = ypeak[iobj]
            thisobj.smash_nsig = ypeak[iobj]/sigma
            sobjs.add_sobj(thisobj)
    else:
        nobj_reg = 0


    # ToDo Also plot the edge trimming boundaries on the QA here.
    if show_peaks:
        spat_approx_vec = slit_left[specmid] + xsize[specmid]*np.arange(nsamp)/nsamp
        spat_approx = slit_left[specmid] + xsize[specmid]*xcen/nsamp
        # Define the plotting function
        #plt.plot(spat_approx_vec, fluxsub/sigma, color ='cornflowerblue',linestyle=':', label='Collapsed Flux')
        plt.plot(spat_approx_vec, fluxconv_cont/sigma, color='black', label = 'Collapsed flux (FWHM convol)')
        plt.plot(spat_approx_vec[cont_mask], fluxconv_cont[cont_mask]/sigma, color='red', markersize=3.0,
                 mfc='red', linestyle='None', fillstyle='full',
                 zorder=9, marker='o', label = 'Used for threshold')
        plt.hlines(threshold/sigma,spat_approx_vec.min(),spat_approx_vec.max(), color='red',linestyle='--', label='Threshold')
        plt.hlines(1.0,spat_approx_vec.min(),spat_approx_vec.max(), color='green',linestyle=':', label='+- 1 sigma')
        plt.hlines(-1.0,spat_approx_vec.min(),spat_approx_vec.max(), color='green',linestyle=':')

        plt.plot(spat_approx, ypeak/sigma, color='red', marker='o', markersize=10.0, mfc='lawngreen', fillstyle='full',
        linestyle='None', zorder = 10,label='Object Found')
        plt.legend()
        plt.xlabel('Approximate Spatial Position (pixels)')
        plt.ylabel('F/sigma (significance)')
        plt.title(qa_title + ': Slit# {:d}'.format(specobj_dict['slitid']))
        plt.show()

    # Now loop over all the regular apertures and assign preliminary traces to them.
    for iobj in range(nobj_reg):
        # Was a standard trace provided? If so, use that as a crutch.
        if std_trace is not None:
            if iobj == 0:
                msgs.info('Using input STANDARD star trace as crutch for object tracing')
            x_trace = np.interp(specmid, spec_vec, std_trace)
            shift = np.interp(specmid, spec_vec,
                              slit_left + xsize * sobjs[iobj].spat_fracpos) - x_trace
            sobjs[iobj].TRACE_SPAT = std_trace + shift
        else:    # If no standard is provided shift left slit boundary over to be initial trace
            # ToDO make this the average left and right boundary instead. That would be more robust.
            sobjs[iobj].TRACE_SPAT = slit_left  + xsize*sobjs[iobj].spat_fracpos
        sobjs[iobj].trace_spec = spec_vec
        sobjs[iobj].spat_pixpos = sobjs[iobj].TRACE_SPAT[specmid]
        # Set the idx for any prelminary outputs we print out. These will be updated shortly
        sobjs[iobj].set_name()

        # Determine the fwhm max
        yhalf = 0.5*sobjs[iobj].smash_peakflux
        xpk = sobjs[iobj].spat_fracpos*nsamp
        x0 = int(np.rint(xpk))
        # TODO It seems we have two codes that do similar things, i.e. findfwhm in arextract.py. Could imagine having one
        # Find right location where smash profile croses yhalf
        if x0 < (int(nsamp)-1):
            ind_righ, = np.where(fluxconv_cont[x0:] < yhalf)
            if len(ind_righ) > 0:
                i2 = ind_righ[0]
                if i2 == 0:
                    xrigh = None
                else:
                    xrigh_int = scipy.interpolate.interp1d(fluxconv_cont[x0 + i2-1:x0 + i2 + 1], x0 + np.array([i2-1,i2],dtype=float),assume_sorted=False)
                    xrigh = xrigh_int([yhalf])[0]
            else:
                xrigh = None
        else:
            xrigh = None
        # Find left location where smash profile crosses yhalf
        if x0 > 0:
            ind_left, = np.where(fluxconv_cont[0:np.fmin(x0+1,int(nsamp)-1)] < yhalf)
            if len(ind_left) > 0:
                i1 = (ind_left[::-1])[0]
                if i1 == (int(nsamp)-1):
                    xleft = None
                else:
                    xleft_int = scipy.interpolate.interp1d(fluxconv_cont[i1:i1+2],np.array([i1,i1+1],dtype=float), assume_sorted= False)
                    xleft = xleft_int([yhalf])[0]
            else:
                xleft = None
        else:
            xleft = None

        if (xleft is None) & (xrigh is None):
            fwhm_measure = None
        elif xrigh is None:
            fwhm_measure = 2.0*(xpk- xleft)
        elif xleft is None:
            fwhm_measure = 2.0*(xrigh - xpk)
        else:
            fwhm_measure = (xrigh - xleft)

        if fwhm_measure is not None:
            sobjs[iobj].FWHM = np.sqrt(np.fmax(fwhm_measure**2 - fwhm**2, (fwhm/2.0)**2)) # Set a floor of fwhm/2 on fwhm
        else:
            sobjs[iobj].FWHM = fwhm


    if (len(sobjs) == 0) & (hand_extract_dict is None):
        msgs.info('No objects found')
        skymask = create_skymask_fwhm(sobjs,thismask)
        return newspecobjs.SpecObjs(), skymask[thismask]


    msgs.info('Fitting the object traces')

    if len(sobjs) > 0:
        # Note the transpose is here to pass in the TRACE_SPAT correctly.
        xinit_fweight = np.copy(sobjs.TRACE_SPAT.T)
        spec_mask = (spec_vec >= spec_min_max[0]) & (spec_vec <= spec_min_max[1])
        trc_inmask = np.outer(spec_mask, np.ones(len(sobjs), dtype=bool))
        xfit_fweight, _, _, _= iter_tracefit(image, xinit_fweight, ncoeff, inmask=inmask, trc_inmask=trc_inmask, fwhm=fwhm,
                                             maxdev=maxdev, idx=sobjs.name,
                                             show_fits=show_fits)
        xinit_gweight = np.copy(xfit_fweight)
        xfit_gweight, _ , _, _= iter_tracefit(image, xinit_gweight, ncoeff,inmask=inmask, trc_inmask=trc_inmask, fwhm=fwhm,
                                              maxdev=maxdev, gweight = True,
                                              idx=sobjs.name, show_fits=show_fits)
        # Linearly extrapolate the traces where they are bad
        #xfit_final = trace_slits.extrapolate_trace(xfit_gweight, spec_min_max, npoly=extrap_npoly)

        # assign the final trace
        for iobj in range(nobj_reg):
            sobjs[iobj].TRACE_SPAT = xfit_gweight[:, iobj]
            sobjs[iobj].spat_pixpos = sobjs[iobj].TRACE_SPAT[specmid]
            sobjs[iobj].set_name()


    # Now deal with the hand apertures if a hand_extract_dict was passed in. Add these to the SpecObj objects
    if hand_extract_dict is not None:
        # First Parse the hand_dict
        hand_extract_spec, hand_extract_spat, hand_extract_det, hand_extract_fwhm = parse_hand_dict(hand_extract_dict)
        # Determine if these hand apertures land on the slit in question
        hand_on_slit = np.where(np.array(thismask[np.rint(hand_extract_spec).astype(int), np.rint(hand_extract_spat).astype(int)]))[0]
        hand_extract_spec = hand_extract_spec[hand_on_slit]
        hand_extract_spat = hand_extract_spat[hand_on_slit]
        hand_extract_det  = hand_extract_det[hand_on_slit]
        hand_extract_fwhm = hand_extract_fwhm[hand_on_slit]
        nobj_hand = len(hand_extract_spec)

        # Decide how to assign a trace to the hand objects
        if nobj_reg > 0:  # Use brightest object on slit?
            smash_peakflux = sobjs.smash_peakflux
            ibri = smash_peakflux.argmax()
            trace_model = sobjs[ibri].TRACE_SPAT
            med_fwhm_reg = np.median(sobjs.FWHM)
        elif std_trace is not None:   # If no objects found, use the standard?
            trace_model = std_trace
        else:  # If no objects or standard use the slit boundary
            trace_model = slit_left
        # Loop over hand_extract apertures and create and assign specobj
        for iobj in range(nobj_hand):
            thisobj = newspecobj.SpecObj(frameshape,
                                       det=specobj_dict['det'],
                                       setup=specobj_dict['setup'], slitid=specobj_dict['slitid'],
                                       orderindx = specobj_dict['orderindx'],
                                       objtype=specobj_dict['objtype'])
            thisobj.hand_extract_spec = hand_extract_spec[iobj]
            thisobj.hand_extract_spat = hand_extract_spat[iobj]
            thisobj.hand_extract_det = hand_extract_det[iobj]
            thisobj.hand_extract_fwhm = hand_extract_fwhm[iobj]
            thisobj.hand_extract_flag = True
            f_ximg = scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, ximg)
            thisobj.spat_fracpos = f_ximg(thisobj.hand_extract_spec, thisobj.hand_extract_spat, grid=False) # interpolate from ximg
            thisobj.smash_peakflux = np.interp(thisobj.spat_fracpos*nsamp,np.arange(nsamp),fluxconv_cont) # interpolate from fluxconv
            # assign the trace
            spat_0 = np.interp(thisobj.hand_extract_spec, spec_vec, trace_model)
            shift = thisobj.hand_extract_spat - spat_0
            thisobj.TRACE_SPAT = trace_model + shift
            thisobj.trace_spec = spec_vec
            thisobj.spat_pixpos = thisobj.TRACE_SPAT[specmid]
            thisobj.set_idx()
            if hand_extract_fwhm[iobj] is not None: # If a hand_extract_fwhm was input use that for the fwhm
                thisobj.FWHM = hand_extract_fwhm[iobj]
            elif nobj_reg > 0: # Otherwise is None was input, then use the median of objects on this slit if they are present
                thisobj.FWHM = med_fwhm_reg
            else:  # Otherwise just use the FWHM parameter input to the code (or the default value)
                thisobj.FWHM = fwhm
            # Finish
            sobjs.add_sobj(thisobj)

    nobj = len(sobjs)
    # If there are no regular aps and no hand aps, just return
    #if nobj == 0:
    #    return (None, skymask, objmask)

    ## Okay now loop over all the regular aps and exclude any which within the fwhm of the hand_extract_APERTURES
    if nobj_reg > 0 and hand_extract_dict is not None:
        spat_pixpos = sobjs.spat_pixpos
        hand_flag = sobjs.hand_extract_flag
        spec_fwhm = sobjs.FWHM
        #spat_pixpos = np.array([spec.spat_pixpos for spec in specobjs])
        #hand_flag = np.array([spec.hand_extract_flag for spec in specobjs])
        #spec_fwhm = np.array([spec.FWHM for spec in specobjs])
        reg_ind, = np.where(~hand_flag)
        hand_ind, = np.where(hand_flag)
        #med_fwhm = np.median(spec_fwhm[~hand_flag])
        #spat_pixpos_hand = spat_pixpos[hand_ind]
        keep = np.ones(nobj,dtype=bool)
        for ihand in hand_ind:
            close = np.abs(sobjs[reg_ind].spat_pixpos - spat_pixpos[ihand]) <= 0.6*spec_fwhm[ihand]
            if np.any(close):
                # Print out a warning
                msgs.warn('Deleting object(s) {}'.format(sobjs[reg_ind[close]].name) +
                          ' because it collides with a user specified hand_extract aperture')
                #for ihand in range(len(close)):
                #    if close[ihand] is True:
                #        msgs.warn('Hand aperture at (hand_extract_spec, hand_extract_spat) = ({:6.2f}'.format(sobjs[hand_ind[ihand]].hand_extract_spec) +
                #                  ',{:6.2f})'.format(sobjs[hand_ind[ihand]].hand_extract_spat) +
                #                  ' lands within 0.6*med_fwhm = {:4.2f}'.format(0.6*med_fwhm) + ' pixels of this object')
                keep[reg_ind[close]] = False

        sobjs = sobjs[keep]

    # Sort objects according to their spatial location
    nobj = len(sobjs)
    spat_pixpos = sobjs.spat_pixpos
    sobjs = sobjs[spat_pixpos.argsort()]
    # Assign integer objids
    sobjs[:].objid = np.arange(nobj)

    # Assign the maskwidth and compute some inputs for the object mask
    xtmp = (np.arange(nsamp) + 0.5)/nsamp
    qobj = np.zeros_like(xtmp)
    for iobj in range(nobj):
        if skythresh > 0.0:
            sobjs[iobj].maskwidth = extract_maskwidth*sobjs[iobj].FWHM*(1.0 + 0.5*np.log10(np.fmax(sobjs[iobj].smash_peakflux/skythresh,1.0)))
        else:
            sobjs[iobj].maskwidth = extract_maskwidth*sobjs[iobj].FWHM
        sep = np.abs(xtmp-sobjs[iobj].spat_fracpos)
        sep_inc = sobjs[iobj].maskwidth/nsamp
        close = sep <= sep_inc
        qobj[close] += sobjs[iobj].smash_peakflux*np.exp(np.fmax(-2.77*(sep[close]*nsamp)**2/sobjs[iobj].FWHM**2,-9.0))

    # Create an objmask. This is created here in case we decide to use it later, but it is not currently used
    skymask_objflux = np.copy(thismask)
    skymask_objflux[thismask] = np.interp(ximg[thismask],xtmp,qobj) < (skymask_nthresh*threshold)
    # Still have to make the skymask
    skymask_fwhm = create_skymask_fwhm(sobjs,thismask)
    skymask = skymask_objflux | skymask_fwhm
    # If requested display the resulting traces on top of the image
    if show_trace:
        viewer, ch = ginga.show_image(image*(thismask*inmask))
        ginga.show_slits(viewer, ch, slit_left.T, slit_righ.T, slit_ids = sobjs[0].slitid)
        for iobj in range(nobj):
            if sobjs[iobj].hand_extract_flag == False:
                color = 'orange'
            else:
                color = 'blue'
            ginga.show_trace(viewer, ch,sobjs[iobj].TRACE_SPAT, trc_name = sobjs[iobj].name, color=color)


    return sobjs, skymask[thismask]

def remap_orders(xinit, spec_min_max, inverse=False):

    """
    This code remaps echelle orders to all extend over the same numer of pixels. It is meant to deal with cases
    where the echelle orders do not completely span the detector. Initial tests with PCA using this remapping
    indicate (for vlt_xshooter_nir) that the remapping does not work as well as simply linear extroplating the traces
    in the PCA. The traces appear to compress better onto a PCA when they are not remapped. So this functionality
    is experimemtnal and not currently used.


    Args:
        xinit: ndarray, (nspec, norders)
           Array of input traces that one wants to remap.
        spec_min_max: ndarray, (2, norders)
           Array containing the spec_min and spec_max defined for each order
        inverse: bool, default = False,
           If True, perform the inverse re-mapping rather than the re-mapping.

    Returns:

    """

    nspec, norders = xinit.shape
    spec_vec = np.arange(nspec)
    spec_vec_norm = spec_vec/float(nspec-1)
    xinit_remap = np.zeros_like(xinit)
    for iord in range(norders):
        igood = (spec_vec >= spec_min_max[0,iord]) & (spec_vec <= spec_min_max[1,iord])
        spec_norm_iord = (spec_vec - spec_vec[igood].min())/(spec_vec[igood].max()  - spec_vec[igood].min())
        if inverse:
            xinit_remap[:, iord] = scipy.interpolate.interp1d(spec_vec_norm, xinit[:, iord], kind='linear',
                                                              bounds_error=False, fill_value='extrapolate')(spec_norm_iord)
        else:
            xinit_remap[:,iord] = scipy.interpolate.interp1d(spec_norm_iord[igood], xinit[igood,iord], kind='linear',
                                                             bounds_error=False, fill_value='extrapolate')(spec_vec_norm)
    return xinit_remap

# TODO: JFH It would be really ideal if we could replace this pca with a weighted PCA!!
def pca_trace(xinit_in, spec_min_max=None, predict = None, npca = None, pca_explained_var=99.0,
              coeff_npoly = None, coeff_weights=None, debug=True, order_vec = None, lower = 3.0,
              upper = 3.0, minv = None,maxv = None, maxrej=1,
              xinit_mean = None):

    """
    Use a PCA model to determine the best object (or slit edge) traces for echelle spectrographs.

    Args:
      xinit:  ndarray, (nspec, norders)
         Array of input traces that one wants to PCA model. For object finding this will be the traces for orders where
         an object was detected. If an object was not detected on some orders (see ech_objfind), the standard star
         (or order boundaries)  will be  assigned to these orders at the correct fractional slit position, and a joint PCA
         fit will be performed to the detected traces and the standard/slit traces.

    spec_min_max: float or int ndarray, (2, norders), default=None.
         This is a 2-d array which defines the minimum and maximum of each order in the
         spectral direction on the detector. This should only be used for echelle spectrographs for which the orders do not
         entirely cover the detector, and each order passed in for xinit_in is a succession of orders on the detector.
         The code will re-map the traces such that they all have the same length, compute the PCA, and then re-map the orders
         back. This improves performanc for echelle spectrographs by removing the nonlinear shrinking of the orders so that
         the linear pca operation can better predict the traces. THIS IS AN EXPERIMENTAL FEATURE. INITIAL TESTS WITH
         XSHOOTER-NIR INDICATED THAT IT DID NOT IMPROVE PERFORMANCE AND SIMPLY LINEAR EXTRAPOLATION OF THE ORDERS INTO THE
         REGIONS THAT ARE NOT ILLUMINATED PERFORMED SIGNIFICANTLY BETTER. DO NOT USE UNTIL FURTHER TESTING IS PERFORMED. IT
         COULD HELP WITH OTHER MORE NONLINEAR SPECTROGRAPHS.
    predict: ndarray, bool (norders,), default = None
         Orders which have True are those that will be predicted by extrapolating the fit of the PCA coefficents for those
         orders which have False set in this array. The default is None, which means that the coefficients of all orders
         will be fit simultaneously and no extrapolation will be performed. For object finding, we use the standard star
         (or slit boundaries) as the input for orders for which a trace is not identified and fit the coefficients of all
         simultaneously. Thus no extrapolation is performed. For tracing slit boundaries it may be useful to perform
          extrapolations.
    npca: int, default = None
         number of PCA components to be kept. The maximum number of possible PCA components would be = norders, which is to say
         that no PCA compression woulud be performed. For the default of None, npca will be automatically determinedy by
         calculating the minimum number of components required to explain 99% (pca_explained_var) of the variance in the different orders.
    pca_explained_var: float, default = 99
         Amount of explained variance cut used to determine where to truncate the PCA, i.e. to determine npca.
    coeff_npoly: int, default = None
         Order of polynomial fits used for PCA coefficients fitting. The defualt is None, which means that coeff_noly
         will be automatically determined by taking the number of orders into account. PCA components that explain
         less variance (and are thus much noiser) are fit with lower order.
    coeff_weights (np.ndarray): shape = (norders,), default=None
         If input these weights will be used for the polynomial fit to the PCA coefficients. Even if you are predicting
         orders and hence only fitting a subset of the orders != norders, the shape of coeff_weights must be norders.
         Just give the orders you don't plan to fit a weight of zero. This option is useful for fitting object
         traces since the weights can be set to (S/N)^2 of each order.
         TODO: Perhaps we should get rid of the predict option and simply allow the user to set the weights of the orders
         they want predicted to be zero. That would be more straightforward, but would require a rework of the code.

    debug: bool, default = False
         Show plots useful for debugging.

    Returns:
    --------
    pca_fit:  ndarray, float (nspec, norders)
        Array with the same size as xinit, which contains the pca fitted orders.
    """

    nspec, norders = xinit_in.shape

    if order_vec is None:
        order_vec = np.arange(norders,dtype=float)

    if predict is None:
        predict = np.zeros(norders,dtype=bool)

    # use_order = True orders used to predict the predict = True bad orders
    use_order = np.invert(predict)
    ngood = np.sum(use_order)

    if ngood < 2:
        msgs.warn('There are no good traces to PCA fit. There is probably a bug somewhere. Exiting and returning input traces.')
        return xinit_in, {}, None, None

    if spec_min_max is not None:
        xinit = remap_orders(xinit_in, spec_min_max)
    else:
        xinit = xinit_in

    # Take out the mean position of each input trace
    if xinit_mean is None:
        xinit_mean = np.mean(xinit, axis=0)

    xpca = xinit - xinit_mean
    xpca_use = xpca[:, use_order].T
    pca_full = PCA()
    pca_full.fit(xpca_use)
    var = np.cumsum(np.round(pca_full.explained_variance_ratio_, decimals=6) * 100)
    npca_full = var.size
    if npca is None:
        if var[0]>=pca_explained_var:
            npca = 1
            msgs.info('The first PCA component contains more than {:5.3f} of the information'.format(pca_explained_var))
        else:
            npca = int(np.ceil(np.interp(pca_explained_var, var,np.arange(npca_full)+1)))
            msgs.info('Truncated PCA to contain {:5.3f}'.format(pca_explained_var) + '% of the total variance. ' +
                      'Number of components to keep is npca = {:d}'.format(npca))
    else:
        npca = int(npca)
        var_trunc = np.interp(float(npca),np.arange(npca_full)+1.0, var)
        msgs.info('Truncated PCA with npca={:d} components contains {:5.3f}'.format(npca, var_trunc) + '% of the total variance.')

    if npca_full < npca:
        msgs.warn('Not enough good traces for a PCA fit of the requested dimensionality. The full (non-compressing) PCA has size: '
                  'npca_full = {:d}'.format(npca_full) + ' is < npca = {:d}'.format(npca))
        msgs.warn('Using the input trace for now. But you should lower npca <= npca_full')
        return xinit_in, {}, None, None

    if coeff_npoly is None:
        coeff_npoly = int(np.fmin(np.fmax(np.floor(3.3*ngood/norders),1.0),3.0))


    # Polynomial coefficient for PCA coefficients
    npoly_vec =np.zeros(npca, dtype=int)
    # Fit first pca dimension (with largest variance) with a higher order npoly depending on number of good orders.
    # Fit all higher dimensions (with lower variance) with a line
    # Cascade down and use lower order polynomial for PCA directions that contain less variance
    for ipoly in range(npca):
        npoly_vec[ipoly] = np.fmax(coeff_npoly - ipoly,1)

        pca = PCA(n_components=npca)
        pca_coeffs_use = pca.fit_transform(xpca_use)
        pca_vectors = pca.components_

    pca_coeffs_new = np.zeros((norders, npca))
    fit_dict = {}
    # Now loop over the dimensionality of the compression and perform a polynomial fit to
    for idim in range(npca):
        # Only fit the use_order orders, then use this to predict the others
        xfit = order_vec[use_order]
        yfit = pca_coeffs_use[:,idim]
        ncoeff = npoly_vec[idim]
        # Apply a 10% relative error to each coefficient. This performs better than use_mad, since larger coefficients
        # will always be considered inliers, if the coefficients vary rapidly with order as they sometimes do.
        sigma = np.fmax(0.1*np.abs(yfit), 0.1)
        invvar = utils.inverse(sigma**2)
        use_weights =  coeff_weights[use_order] if coeff_weights is not None else None
        # TODO Note that we are doing a weighted fit using the coeff_weights, but the rejection is still done
        # usnig the ad-hoc invvar created in the line above. I cannot think of a better way.
        msk_new, poly_out = utils.robust_polyfit_djs(xfit, yfit, ncoeff, invvar = invvar, weights=use_weights,
                                                     function='polynomial', maxiter=25,
                                                     lower=lower, upper=upper,
                                                     maxrej=maxrej,
                                                     sticky=False, use_mad=False, minx = minv, maxx = maxv)
        # ToDO robust_poly_fit needs to return minv and maxv as outputs for the fits to be usable downstream
        pca_coeffs_new[:,idim] = utils.func_val(poly_out, order_vec, 'polynomial')
        fit_dict[str(idim)] = {}
        fit_dict[str(idim)]['coeffs'] = poly_out
        fit_dict[str(idim)]['minv'] = minv
        fit_dict[str(idim)]['maxv'] = maxv
        if debug:
            # Evaluate the fit
            xvec = np.linspace(order_vec.min(),order_vec.max(),num=100)
            robust_mask_new = msk_new == 1
            plt.plot(xfit, yfit, 'ko', mfc='None', markersize=8.0, label='pca coeff')
            plt.plot(xfit[~robust_mask_new], yfit[~robust_mask_new], 'r+', markersize=20.0,label='robust_polyfit_djs rejected')
            plt.plot(xvec, utils.func_val(poly_out, xvec, 'polynomial'),ls='-.', color='steelblue',
                     label='Polynomial fit of order={:d}'.format(ncoeff))
            plt.xlabel('Order Number', fontsize=14)
            plt.ylabel('PCA Coefficient', fontsize=14)
            plt.title('PCA Fit for Dimension #{:d}/{:d}'.format(idim + 1,npca))
            plt.legend()
            plt.show()

    pca_model = np.outer(pca.mean_,np.ones(norders)) + (np.dot(pca_coeffs_new, pca_vectors)).T
#   pca_model_mean = np.mean(pca_model,0)
#   pca_fit = np.outer(np.ones(nspec), xinit_mean) + (pca_model - pca_model_mean)
#   JFH which is correct?
    pca_fit = np.outer(np.ones(nspec), xinit_mean) + (pca_model)

    if spec_min_max is not None:
        pca_out = remap_orders(pca_fit, spec_min_max, inverse=True)
    else:
        pca_out = pca_fit

    return pca_out, fit_dict, pca.mean_, pca_vectors



def ech_objfind(image, ivar, slitmask, slit_left, slit_righ, inmask=None, spec_min_max=None,
                fof_link=1.5, order_vec=None, plate_scale=0.2, ir_redux=False,
                std_trace=None, extrap_npoly=3, ncoeff=5, npca=None, coeff_npoly=None, max_snr=2.0, min_snr=1.0, nabove_min_snr=2,
                pca_explained_var=99.0, box_radius=2.0, fwhm=3.0, maxdev=2.0, hand_extract_dict=None, nperslit=5, bg_smth=5.0,
                extract_maskwidth=3.0, sig_thresh = 10.0, peak_thresh=0.0, abs_thresh=0.0, specobj_dict=None,
                trim_edg=(5,5), cont_fit=True, npoly_cont=1, show_peaks=False, show_fits=False, show_single_fits=False,
                show_trace=False, show_single_trace=False, debug=False, show_pca=False, debug_all=False):
    """
    Object finding routine for Echelle spectrographs. This routine:
       1) runs object finding on each order individually
       2) Links the objects found together using a friends-of-friends algorithm on fractional order position.
       3) For objects which were only found on some orders, the standard (or the slit boundaries) are placed at the appropriate
          fractional position along the order.
       4) A PCA fit to the traces is performed using the routine above pca_fit

    Args:
    image:  float ndarray, shape (nspec, nspat)
        Image to search for objects from. This image has shape (nspec, nspat) where the first dimension (nspec) is spectral,
        and second dimension (nspat) is spatial. Note this image can either have the sky background in it, or have already been sky subtracted.
        Object finding works best on sky-subtracted images. Ideally objfind would be run in another routine, global sky-subtraction performed, and
        then this code should be run. However, it is also possible to run this code on non sky subtracted images.
    ivar: float ndarray, shape (nspec, nspat)
       Inverse variance image for the input image.
    slitmask: int ndarray, shape (nspec, nspat)
       Integer image indicating the pixels that belong to each order. Pixels that are not on an order have value -1, and those
       that are on an order have a value equal to the slit number (i.e. 0 to nslits-1 from left to right on the image)
    slit_left:  float ndarray
        Left boundary of orders to be extracted (given as floating pt pixels). This a 2-d array with shape (nspec, norders)
    slit_righ:  float ndarray
        Left boundary of orders to be extracted (given as floating pt pixels). This a 2-d array with shape (nspec, norders)
    inmask: ndarray, bool, shape (nspec, nspat), default = None
        Input mask for the input image.
    fwhm: float, default = 3.0
        Estimated fwhm of the objects in pixels
    maxdev (float): default=2.0
        Maximum deviation of pixels from polynomial fit to trace used to reject bad pixels in trace fitting.
    spec_min_max: float or int ndarray, (2, norders), default=None. This is a 2-d array which defines the minimum and maximum of each order in the
       spectral direction on the detector. This should only be used for echelle spectrographs for which the orders do not
       entirely cover the detector. The pca_trace code will re-map the traces such that they all have the same length,
       compute the PCA, and then re-map the orders back. This improves performanc for echelle spectrographs by removing
       the nonlinear shrinking of the orders so that the linear pca operation can better predict the traces. If not
       passed in it will be determined automitically from the slitmask
    fof_link: float, default = 1.5"
        Friends-of-friends linking length in arcseconds used to link together traces across orders. The routine links
        together at the same fractional slit position and links them together with a friends-of-friends algorithm using
        this linking length.
    plate_scale: float or ndarray, if an ndarray shape is (norders,) default = 0.2
       plate scale of your detector, in unit of arcsec/pix. This can either be a single float for every order, or an array
       with size norders indicating the plate scale of each order.
    ncoeff: int, default = 5
       Order of polynomial fit to traces
    npca: int, default = None
       Nmber of PCA components you want to keep. default is None and it will be assigned automatically by calculating
       the number of components contains approximately 99% of the variance
    coeff_npoly: int, default = None,
       order of polynomial used for PCA coefficients fitting. Default is None and this will be determined automatically.
    min_snr: float, default = 0.2
       Minimum SNR for keeping an object. For an object to be kept it must have a median S/N ratio above min_snr for
       at least nabove_min_snr orders.
    nabove_min_snr: int, default = 1
       The required number of orders that an object must have with median SNR>min_snr in order to be kept.
    box_radius: float,
      box_car extraction radius in arcseconds for SNR calculation and trimming
    sig_thresh: threshord for finding objects
    show_peaks: whether plotting the QA of peak finding of your object in each order
    show_fits: Plot trace fitting for final fits using PCA as crutch
    show_single_fits: Plot trace fitting for single order fits
    show_trace: whether display the resulting traces on top of the image
    debug:

    Returns
    -------
    (sobjs, skymask): tuple

    sobjs: object
      Specobjs object containing the objects detected
    skymask: float ndarray, same shape as image
      Skymask indicating which pixels can be used for global sky subtraction
    """

    #debug_all=True
    if debug_all:
        show_peaks = True
        #show_fits = True
        #show_single_fits = True
        show_trace = True
        show_pca = True
        #show_single_trace = True
        debug = True


    if specobj_dict is None:
        specobj_dict = {'setup': 'unknown', 'slitid': 999, 'orderindx': 999,
                        'det': 1, 'objtype': 'unknown', 'pypeline': 'Echelle'}


    # TODO Update FOF algorithm here with the one from scikit-learn.

    allmask = slitmask > -1
    if inmask is None:
        inmask = allmask

    frameshape = image.shape
    nspec = frameshape[0]
    nspat = frameshape[1]
    norders = slit_left.shape[1]

    if spec_min_max is None:
        spec_min_max = np.zeros(2,norders)
        for iord in range(norders):
            ispec, ispat = np.where(slitmask == iord)
            spec_min_max[:,iord] = ispec.min(), ispec.max()

    if order_vec is None:
        order_vec = np.arange(norders)

    if isinstance(plate_scale,(float, int)):
        plate_scale_ord = np.full(norders, plate_scale)
    elif isinstance(plate_scale,(np.ndarray, list, tuple)):
        if len(plate_scale) == norders:
            plate_scale_ord = plate_scale
        elif len(plate_scale) == 1:
            plate_scale_ord = np.full(norders, plate_scale[0])
        else:
            msgs.error('Invalid size for plate_scale. It must either have one element or norders elements')
    else:
        msgs.error('Invalid type for plate scale')

    specmid = nspec // 2
    spec_vec = np.arange(nspec)
    slit_width = slit_righ - slit_left
    slit_spec_pos = nspec/2.0
    slit_spat_pos = np.zeros((norders, 2))
    for iord in range(norders):
        slit_spat_pos[iord, :] = (np.interp(slit_spec_pos, spec_vec, slit_left[:,iord]), np.interp(slit_spec_pos, spec_vec, slit_righ[:,iord]))

    # create the ouptut images skymask and objmask
    skymask_objfind = np.copy(allmask)
    # Loop over orders and find objects
    sobjs = newspecobjs.SpecObjs()
    # ToDo replace orderindx with the true order number here? Maybe not. Clean up slitid and orderindx!
    for iord in range(norders):
        msgs.info('Finding objects on order # {:d}'.format(order_vec[iord]))
        thismask = slitmask == iord
        inmask_iord = inmask & thismask
        specobj_dict['slitid'] = iord
        specobj_dict['orderindx'] = iord
        try:
            std_in = std_trace[:,iord]
        except TypeError:
            std_in = None
        sobjs_slit, skymask_objfind[thismask] = \
            objfind(image, thismask, slit_left[:,iord], slit_righ[:,iord], spec_min_max=spec_min_max[:,iord],
                    inmask=inmask_iord,std_trace=std_in, extrap_npoly=extrap_npoly, ncoeff=ncoeff, fwhm=fwhm, maxdev=maxdev,
                    hand_extract_dict=hand_extract_dict, ir_redux=ir_redux,
                    nperslit=nperslit, bg_smth=bg_smth, extract_maskwidth=extract_maskwidth, sig_thresh=sig_thresh,
                    peak_thresh=peak_thresh, abs_thresh=abs_thresh, trim_edg=trim_edg, cont_fit=cont_fit,
                    npoly_cont=npoly_cont, show_peaks=show_peaks,
                    show_fits=show_single_fits, show_trace=show_single_trace, specobj_dict=specobj_dict)
        # ToDO make the specobjs _set_item_ work with expressions like this spec[:].orderindx = iord
        for spec in sobjs_slit:
            spec.ech_orderindx = iord
            spec.ech_order = order_vec[iord]
            _ = spec.set_idx()
        sobjs.add_sobj(sobjs_slit)

    nfound = len(sobjs)

    if nfound == 0:
        return sobjs, skymask_objfind[allmask]

    FOF_frac = fof_link/(np.median(np.median(slit_width,axis=0)*plate_scale_ord))
    # Run the FOF. We use fake coordinaes
    fracpos = sobjs.spat_fracpos
    ra_fake = fracpos/1000.0 # Divide all angles by 1000 to make geometry euclidian
    dec_fake = 0.0*fracpos
    if nfound>1:
        (inobj_id, multobj_id, firstobj_id, nextobj_id) = spheregroup(ra_fake, dec_fake, FOF_frac/1000.0)
        obj_id_init = inobj_id.copy()
    elif nfound==1:
        obj_id_init = np.zeros(1,dtype='int')

    uni_obj_id_init, uni_ind_init = np.unique(obj_id_init, return_index=True)

    # Now loop over the unique objects and check that there is only one object per order. If FOF
    # grouped > 1 objects on the same order, then this will be popped out as its own unique object
    obj_id = obj_id_init.copy()
    nobj_init = len(uni_obj_id_init)
    for iobj in range(nobj_init):
        for iord in range(norders):
            on_order = (obj_id_init == uni_obj_id_init[iobj]) & (sobjs.ech_orderindx == iord)
            if (np.sum(on_order) > 1):
                msgs.warn('Found multiple objects in a FOF group on order iord={:d}'.format(iord) + msgs.newline() +
                          'Spawning new objects to maintain a single object per order.')
                off_order = (obj_id_init == uni_obj_id_init[iobj]) & (sobjs.ech_orderindx != iord)
                ind = np.where(on_order)[0]
                if np.any(off_order):
                    # Keep the closest object to the location of the rest of the group (on other orders)
                    # as corresponding to this obj_id, and spawn new obj_ids for the others.
                    frac_mean = np.mean(fracpos[off_order])
                    min_dist_ind = np.argmin(np.abs(fracpos[ind] - frac_mean))
                else:
                    # If there are no other objects with this obj_id to compare to, then we simply have multiple
                    # objects grouped together on the same order, so just spawn new object IDs for them to maintain
                    # one obj_id per order
                    min_dist_ind = 0
                ind_rest = np.setdiff1d(ind,ind[min_dist_ind])
                obj_id[ind_rest] = (np.arange(len(ind_rest)) + 1) + obj_id_init.max()

    uni_obj_id, uni_ind = np.unique(obj_id, return_index=True)
    nobj = len(uni_obj_id)
    msgs.info('FOF matching found {:d}'.format(nobj) + ' unique objects')

    gfrac = np.zeros(nfound)
    for jj in range(nobj):
        this_obj_id = obj_id == uni_obj_id[jj]
        gfrac[this_obj_id] = np.median(fracpos[this_obj_id])

    uni_frac = gfrac[uni_ind]

    # Sort with respect to fractional slit location to guarantee that we have a similarly sorted list of objects later
    isort_frac = uni_frac.argsort()
    uni_obj_id = uni_obj_id[isort_frac]
    uni_frac = uni_frac[isort_frac]

    sobjs_align = sobjs.copy()
    # Loop over the orders and assign each specobj a fractional position and a obj_id number
    for iobj in range(nobj):
        for iord in range(norders):
            on_order = (obj_id == uni_obj_id[iobj]) & (sobjs_align.ech_orderindx == iord)
            # ToDO fix specobjs set_item to get rid of these crappy loops
            for spec in sobjs_align[on_order]:
                spec.ech_fracpos = uni_frac[iobj]
                spec.ech_objid = uni_obj_id[iobj]
                spec.objid = uni_obj_id[iobj]
                spec.ech_frac_was_fit = False

    # Now loop over objects and fill in the missing objects and their traces. We will fit the fraction slit position of
    # the good orders where an object was found and use that fit to predict the fractional slit position on the bad orders
    # where no object was found
    ## TODO
    for iobj in range(nobj):
        # Grab all the members of this obj_id from the object list
        indx_obj_id = sobjs_align.ech_objid == uni_obj_id[iobj]
        nthisobj_id = np.sum(indx_obj_id)
        # Perform the fit if this objects shows up on more than three orders
        if (nthisobj_id > 3) and (nthisobj_id<norders):
            thisorderindx = sobjs_align[indx_obj_id].ech_orderindx
            goodorder = np.zeros(norders, dtype=bool)
            goodorder[thisorderindx] = True
            badorder = np.invert(goodorder)
            xcen_good = (sobjs_align[indx_obj_id].TRACE_SPAT).T
            slit_frac_good = (xcen_good-slit_left[:,goodorder])/slit_width[:,goodorder]
            # Fractional slit position averaged across the spectral direction for each order
            frac_mean_good = np.mean(slit_frac_good, 0)
            # Perform a  linear fit to fractional slit position
            #TODO Do this as a S/N weighted fit similar to what is now in the pca_trace algorithm?
            msk_frac, poly_coeff_frac = utils.robust_polyfit_djs(order_vec[goodorder], frac_mean_good, 1,
                                                                 function='polynomial', maxiter=20, lower=2, upper=2,
                                                                 use_mad= True, sticky=False,
                                                                 minx = order_vec.min(), maxx=order_vec.max())
            frac_mean_new = np.zeros(norders)
            frac_mean_new[badorder] = utils.func_val(poly_coeff_frac, order_vec[badorder], 'polynomial',
                                                     minx = order_vec.min(),maxx=order_vec.max())
            frac_mean_new[goodorder] = frac_mean_good
            # TODO This QA needs some work
            if show_pca:
                frac_mean_fit = utils.func_val(poly_coeff_frac, order_vec, 'polynomial')
                plt.plot(order_vec[goodorder][msk_frac], frac_mean_new[goodorder][msk_frac], 'ko', mfc='k', markersize=8.0, label='Good Orders Kept')
                plt.plot(order_vec[goodorder][np.invert(msk_frac)], frac_mean_new[goodorder][np.invert(msk_frac)], 'ro', mfc='k', markersize=8.0, label='Good Orders Rejected')
                plt.plot(order_vec[badorder], frac_mean_new[badorder], 'ko', mfc='None', markersize=8.0, label='Predicted Bad Orders')
                plt.plot(order_vec,frac_mean_new,'+',color='cyan',markersize=12.0,label='Final Order Fraction')
                plt.plot(order_vec, frac_mean_fit, 'r-', label='Fractional Order Position Fit')
                plt.xlabel('Order Index', fontsize=14)
                plt.ylabel('Fractional Slit Position', fontsize=14)
                plt.title('Fractional Slit Position Fit')
                plt.legend()
                plt.show()
        else:
            frac_mean_new = np.full(norders, uni_frac[iobj])

        # Now loop over the orders and add objects on the ordrers for which the current object was not found
        for iord in range(norders):
            # Is the current object detected on this order?
            on_order = (sobjs_align.ech_objid == uni_obj_id[iobj]) & (sobjs_align.ech_orderindx == iord)
            if not np.any(on_order):
                # Add this to the sobjs_align, and assign required tags
                thisobj = newspecobj.SpecObj(frameshape, slit_spat_pos[iord,:], det = sobjs_align[0].det,
                                           setup = sobjs_align[0].setup, slitid = iord,
                                           objtype=sobjs_align[0].objtype, pypeline=sobjs_align[0].pypeline)
                thisobj.ech_orderindx = iord
                thisobj.ech_order = order_vec[iord]
                thisobj.spat_fracpos = uni_frac[iobj]
                # Assign traces using the fractional position fit above
                if std_trace is not None:
                    x_trace = np.interp(slit_spec_pos, spec_vec, std_trace[:,iord])
                    shift = np.interp(slit_spec_pos, spec_vec,slit_left[:,iord] + slit_width[:,iord]*frac_mean_new[iord]) - x_trace
                    thisobj.TRACE_SPAT = std_trace[:,iord] + shift
                else:
                    thisobj.TRACE_SPAT = slit_left[:, iord] + slit_width[:, iord] * frac_mean_new[iord]  # new trace
                thisobj.trace_spec = spec_vec
                thisobj.spat_pixpos = thisobj.TRACE_SPAT[specmid]
                # Use the real detections of this objects for the FWHM
                this_obj_id = obj_id == uni_obj_id[iobj]
                # Assign to the fwhm of the nearest detected order
                imin = np.argmin(np.abs(sobjs_align[this_obj_id].ech_orderindx - iord))
                thisobj.FWHM = sobjs_align[imin].FWHM
                thisobj.maskwidth = sobjs_align[imin].maskwidth
                thisobj.ech_fracpos = uni_frac[iobj]
                thisobj.ech_objid = uni_obj_id[iobj]
                #thisobj.objid = uni_obj_id[iobj]
                thisobj.ech_frac_was_fit = True
                thisobj.set_idx()
                sobjs_align.add_sobj(thisobj)
                obj_id = np.append(obj_id, uni_obj_id[iobj])
                gfrac = np.append(gfrac, uni_frac[iobj])

    # Loop over the objects and perform a quick and dirty extraction to assess S/N.
    varimg = utils.calc_ivar(ivar)
    flux_box = np.zeros((nspec, norders, nobj))
    ivar_box = np.zeros((nspec, norders, nobj))
    mask_box = np.zeros((nspec, norders, nobj))
    SNR_arr = np.zeros((norders, nobj))
    for iobj in range(nobj):
        for iord in range(norders):
            indx = (sobjs_align.ech_objid == uni_obj_id[iobj]) & (sobjs_align.ech_orderindx == iord)
            spec = sobjs_align[indx][0]
            thismask = slitmask == iord
            inmask_iord = inmask & thismask
            # TODO make the snippet below its own function quick_extraction()
            box_rad_pix = box_radius/plate_scale_ord[iord]
            flux_tmp  = extract_boxcar(image*inmask_iord, spec.TRACE_SPAT,box_rad_pix, ycen = spec.trace_spec)
            var_tmp  = extract_boxcar(varimg*inmask_iord, spec.TRACE_SPAT,box_rad_pix, ycen = spec.trace_spec)
            ivar_tmp = utils.calc_ivar(var_tmp)
            pixtot  = extract_boxcar(ivar*0 + 1.0, spec.TRACE_SPAT,box_rad_pix, ycen = spec.trace_spec)
            mask_tmp = (extract_boxcar(ivar*inmask_iord == 0.0, spec.TRACE_SPAT,box_rad_pix, ycen = spec.trace_spec) != pixtot)
            flux_box[:,iord,iobj] = flux_tmp*mask_tmp
            ivar_box[:,iord,iobj] = np.fmax(ivar_tmp*mask_tmp,0.0)
            mask_box[:,iord,iobj] = mask_tmp
            (mean, med_sn, stddev) = stats.sigma_clipped_stats(flux_box[mask_tmp,iord,iobj]*np.sqrt(ivar_box[mask_tmp,iord,iobj]),
            sigma_lower=5.0,sigma_upper=5.0)
            # ToDO assign this to sobjs_align for use in the extraction
            SNR_arr[iord,iobj] = med_sn
            spec.ech_snr = med_sn

    # Purge objects with low SNR that don't show up in enough orders, sort the list of objects with respect to obj_id
    # and orderindx
    keep_obj = np.zeros(nobj,dtype=bool)
    sobjs_trim = newspecobjs.SpecObjs()
    # objids are 1 based so that we can easily asign the negative to negative objects
    iobj_keep = 1
    for iobj in range(nobj):
        if (SNR_arr[:,iobj].max() > max_snr) or (np.sum(SNR_arr[:,iobj] > min_snr) >= nabove_min_snr):
            keep_obj[iobj] = True
            ikeep = sobjs_align.ech_objid == uni_obj_id[iobj]
            sobjs_keep = sobjs_align[ikeep].copy()
            for spec in sobjs_keep:
                spec.ech_objid = iobj_keep
                #spec.objid = iobj_keep
            sobjs_trim.add_sobj(sobjs_keep[np.argsort(sobjs_keep.ech_orderindx)])
            iobj_keep += 1
        else:
            msgs.info('Purging object #{:d}'.format(iobj) + ' which does not satisfy max_snr > {:5.2f} OR min_snr > {:5.2f}'.format(max_snr, min_snr) +
                      ' on at least nabove_min_snr >= {:d}'.format(nabove_min_snr) + ' orders')

    nobj_trim = np.sum(keep_obj)

    if nobj_trim == 0:
        sobjs_final = newspecobjs.SpecObjs()
        skymask = create_skymask_fwhm(sobjs_final, allmask)
        return sobjs_final, skymask[allmask]

    SNR_arr_trim = SNR_arr[:,keep_obj]


    sobjs_final = sobjs_trim.copy()
    # Loop over the objects one by one and adjust/predict the traces
    pca_fits = np.zeros((nspec, norders, nobj_trim))

    # Create the trc_inmask for iterative fitting below
    trc_inmask = np.zeros((nspec, norders), dtype=bool)
    for iord in range(norders):
        trc_inmask[:,iord] = (spec_vec >= spec_min_max[0,iord]) & (spec_vec <= spec_min_max[1,iord])


    for iobj in range(nobj_trim):
        indx_obj_id = sobjs_final.ech_objid == (iobj + 1)
        # PCA predict all the orders now (where we have used the standard or slit boundary for the bad orders above)
        msgs.info('Fitting echelle object finding PCA for object {:d}\{:d} with median SNR = {:5.3f}'.format(
                iobj + 1,nobj_trim,np.median(sobjs_final[indx_obj_id].ech_snr)))
        pca_fits[:, :, iobj], _, _, _ = pca_trace((sobjs_final[indx_obj_id].TRACE_SPAT).T,
                                                 npca=npca, pca_explained_var=pca_explained_var,
                                                 coeff_npoly=coeff_npoly,
                                                 coeff_weights = np.fmax(sobjs_final[indx_obj_id].ech_snr, 1.0),
                                                 debug=show_pca)
        # Trial and error shows weighting by S/N instead of S/N^2 performs better

        # TODO: Check transpose
#        pca_fits[:, :, iobj] = tracepca.pca_trace_object(sobjs_final[indx_obj_id].TRACE_SPAT.T,
#                                                         order=coeff_npoly, npca=npca,
#                                                         pca_explained_var=pca_explained_var,
#                                                         debug=debug)

        # Perform iterative flux weighted centroiding using new PCA predictions
        xinit_fweight = pca_fits[:,:,iobj].copy()
        inmask_now = inmask & allmask
        xfit_fweight, _, _, _= iter_tracefit(image, xinit_fweight, ncoeff, inmask = inmask_now, trc_inmask=trc_inmask,
                                             fwhm=fwhm, maxdev=maxdev, show_fits=show_fits)
        # Perform iterative Gaussian weighted centroiding
        xinit_gweight = xfit_fweight.copy()
        xfit_gweight, _ , _, _= iter_tracefit(image, xinit_gweight, ncoeff, inmask = inmask_now,  trc_inmask=trc_inmask,
                                              gweight=True, fwhm=fwhm, maxdev=maxdev, show_fits=show_fits)
        #xfit_final = trace_slits.extrapolate_trace(xfit_gweight, spec_min_max, npoly=extrap_npoly)
        #TODO  Assign the new traces. Only assign the orders that were not orginally detected and traced. If this works
        # well, we will avoid doing all of the iter_tracefits above to make the code faster.
        for iord, spec in enumerate(sobjs_final[indx_obj_id]):
            # JFH added the condition on ech_frac_was_fit with S/N cut on 7-7-19.
            # TODO is this robust against half the order being masked?
            if spec.ech_frac_was_fit & (spec.ech_snr > 1.0):
                    spec.TRACE_SPAT = xfit_gweight[:,iord]
                    spec.spat_pixpos = spec.TRACE_SPAT[specmid]




    #TODO Put in some criterion here that does not let the fractional position change too much during the iterative
    # tracefitting. The problem is spurious apertures identified on one slit can be pulled over to the center of flux
    # resulting in a bunch of objects landing on top of each other.

    # Set the IDs
    for spec in sobjs_final:
        spec.ech_order = order_vec[spec.ech_orderindx]
    sobjs_final.set_idx()

    skymask_fwhm = create_skymask_fwhm(sobjs_final,allmask)
    skymask = skymask_objfind | skymask_fwhm

    if show_trace:
        viewer, ch = ginga.show_image(image*allmask)

        for spec in sobjs_trim:
            color = 'red' if spec.ech_frac_was_fit else 'magenta'
            ## Showing the final flux weighted centroiding from PCA predictions
            ginga.show_trace(viewer, ch, spec.TRACE_SPAT, spec.name, color=color)


        for iobj in range(nobj_trim):
            for iord in range(norders):
                ## Showing PCA predicted locations before recomputing flux/gaussian weighted centroiding
                ginga.show_trace(viewer, ch, pca_fits[:,iord, iobj], str(uni_frac[iobj]), color='yellow')
                ## Showing the final traces from this routine
                ginga.show_trace(viewer, ch, sobjs_final.TRACE_SPAT[iord].T, sobjs_final.name, color='cyan')


        # Labels for the points
        text_final = [dict(type='text', args=(nspat / 2 -40, nspec / 2, 'final trace'),
                           kwargs=dict(color='cyan', fontsize=20))]

        text_pca = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 30, 'PCA fit'),kwargs=dict(color='yellow', fontsize=20))]

        text_fit = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 60, 'predicted'),kwargs=dict(color='red', fontsize=20))]

        text_notfit = [dict(type='text', args=(nspat / 2 -40, nspec / 2 - 90, 'originally found'),kwargs=dict(color='magenta', fontsize=20))]

        canvas = viewer.canvas(ch._chname)
        canvas_list = text_final + text_pca + text_fit + text_notfit
        canvas.add('constructedcanvas', canvas_list)
    # TODO two things need to be debugged. 1) For objects which were found and traced, i don't think we should be updating the tracing with
    # the PCA. This just adds a failutre mode. 2) The PCA fit is going wild for X-shooter. Debug that.
    return sobjs_final, skymask[allmask]
