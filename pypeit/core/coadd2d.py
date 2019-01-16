""" Module for image processing core methods
"""
from __future__ import (print_function, absolute_import, division, unicode_literals)

import astropy.stats
import numpy as np
import glob
import os
from pypeit import msgs
from pypeit import utils
from pypeit import masterframe

def get_coadd2d_files(redux_path, objprefix):
    """
    Grabs the files necessary to perform a 2d coadd.

    Args:
        redux_path: str
           redux path, e.g. '/Users/joe/python/PypeIt-development-suite/REDUX_OUT/Gemini_GNIRS/'
        objprefix:
           object prefix, i.e. if the files are spec1d_pisco_*, then this would be 'pisco'

    Returns:
        (spec2d_files, spec1d_files, tracefiles, tiltfiles, waveimgfiles)

    """
    spec2d_files = glob.glob(redux_path + 'Science/spec2d_' + objprefix + '*')
    spec1d_files = [files.replace('spec2d', 'spec1d') for files in spec2d_files]
    # Get the master dir
    head0 = fits.getheader(spec2d_files[0])
    mdir = os.path.basename(head0['PYPMFDIR'])+'/'
    master_path = redux_path + mdir
    tiltfiles = []
    waveimgfiles = []
    tracefiles = []
    for file in spec2d_files:
        head = fits.getheader(file)
        trace_key = '{:s}'.format(head['TRACMKEY'])
        wave_key = '{:s}'.format(head['ARCMKEY'])
        tracefiles.append(masterframe.master_name('trace', trace_key, master_path))
        waveimgfiles.append(masterframe.master_name('wave', wave_key, master_path))
        tiltfiles.append(masterframe.master_name('tilts', wave_key, master_path))

    return spec2d_files, spec1d_files, tracefiles, tiltfiles, waveimgfiles


def get_wave_ind(wave_grid, wave_min, wave_max):
    """
    Utility routine used by coadd2d to determine the starting and ending indices of a wavelength grid.

    Args:
        wave_grid: float ndarray
          Wavelength grid.
        wave_min: float
          Minimum wavelength covered by the data in question.
        wave_max: float
          Maximum wavelength covered by the data in question.

    Returns:
        (ind_lower, ind_upper): tuple, int
          Integer lower and upper indices into the array wave_grid that cover the interval (wave_min, wave_max)
    """

    diff = wave_grid - wave_min
    diff[diff > 0] = np.inf
    ind_lower = np.argmin(np.abs(diff))
    diff = wave_max - wave_grid
    diff[diff > 0] = np.inf
    ind_upper = np.argmin(np.abs(diff))

    return ind_lower, ind_upper


def coadd2d(sciimg_stack, sciivar_stack, skymodel_stack, inmask_stack, tilts_stack, waveimg_stack, trace_stack,
            thismask_stack, weights=None, loglam_grid=None, wave_grid=None):
    """
    This routine will perform a 2d co-add of a stack of PypeIt spec2d reduction outputs. The slit will be
    'rectified' onto a spatial and spectral grid, which encompasses the spectral and spatial coverage of the image stacks.
    The rectification uses nearest grid point interpolation to avoid covariant errors.
    Dithering is supported as all images are centered relative to a set of reference traces in trace_stack.

    Args:
    -----
        sciimg_stack: float ndarray, shape = (nimgs, nspec, nspat)
           Stack of science images
        sciivar_stack: float ndarray, shape = (nimgs, nspec, nspat)
           Stack of inverse variance images
        skymodel_stack: float ndarray, shape = (nimgs, nspec, nspat)
           Stack of the model sky
        inmask_stack: bool ndarray, shape = (nimgs, nspec, nspat)
           Stack of input masks. True = Good, False=Bad
        tilts_stack: float ndarray, shape = (nimgs, nspec, nspat)
           Stack of the wavelength tilts
        waveimg_stack: float ndarray, shape = (nimgs, nspec, nspat)
           Stack of the wavelength images
        trace_stack: float ndarray, shape (nimgs, nspec)
           Stack of reference traces about which the images are rectified and coadded.
           If the images were not dithered then this reference trace can simply be the center
           of the slit slitcen = (slit_left + slit_righ)/2. If the images were dithered, then this could trace_stack
           could either be the slitcen appropriately shifted with the dither pattern, or it could the trace of the object
           of interest in each exposure determined by running PypeIt on the individual images.
        thismask_stack: bool ndarray, shape = (nimgs, nspec, nspat)
           Stack of masks indicating which pixels are on the slit in question. True = On slit, False=Off slit.

    Optional Args:
    --------------
        weights: float ndarray, shape = (nimgs,), default = None
            Set of weights used for the weighted combined of the rectified images using weighted_combine. If the
            weights are not provided then a uniform weighting will be adopted.
        loglam_grid: float ndarray, shape = any, default = None
            Wavelength grid in log10(wave) onto which the image stacks will be rectified. The code will automatically
            choose the subset of this grid encompassing the wavelength coverage of the image stacks proviced
            (using waveimg_stack). Either loglam grid or wave_grid need to be provided
        wave_grid:  float ndarray, shape = any, default = None
            Same as loglam_grid but allowing a grid in wave to be provided instead of the log10(wave).

    Returns:
        (sciimg, sciivar, imgminsky, outmask, nused, tilts, waveimg, dspat, thismask, tslits_dict)

        sciimg: float ndarray shape = (nspec_coadd, nspat_coadd)
            Rectified and coadded science image
        sciivar: float ndarray shape = (nspec_coadd, nspat_coadd)
            Rectified and coadded inverse variance image with correct error propagation
        imgminsky: float ndarray shape = (nspec_coadd, nspat_coadd)
            Rectified and coadded sky subtracted image
        outmask: bool ndarray shape = (nspec_coadd, nspat_coadd)
            Output mask for rectified and coadded images. True = Good, False=Bad.
        nused: int ndarray shape = (nspec_coadd, nspat_coadd)
            Image of integers indicating the number of images from the image stack that contributed to each pixel
        tilts: float ndarray shape = (nspec_coadd, nspat_coadd)
            The averaged tilts image corresponding to the rectified and coadded data.
        waveimg: float ndarray shape = (nspec_coadd, nspat_coadd)
            The averaged wavelength image corresponding to the rectified and coadded data.
        dspat: float ndarray shape = (nspec_coadd, nspat_coadd)
            The average spatial offsets in pixels from the reference trace trace_stack corresponding to the rectified
            and coadded data.
        thismask: bool ndarray shape = (nspec_coadd, nspat_coadd)
            Output mask for rectified and coadded images. True = Good, False=Bad. This image is trivial, and
            is simply an image of True values the same shape as the rectified and coadded data.
        tslits_dict: dict
            tslits_dict dictionary containing the information about the slits boundaries. The slit boundaries
            are trivial and are simply vertical traces at 0 and nspat_coadd-1.
    """
    nimgs, nspec, nspat = sciimg_stack.shape

    # Determine the wavelength grid that we will use for the current slit/order
    wave_min = waveimg_stack[thismask_stack].min()
    wave_max = waveimg_stack[thismask_stack].max()
    if loglam_grid is not None:
        ind_lower, ind_upper = get_wave_ind(loglam_grid, np.log10(wave_min), np.log10(wave_max))
        loglam_bins = loglam_grid[ind_lower:ind_upper + 1]
        wave_bins = np.power(10.0, loglam_bins)
    elif wave_grid is not None:
        ind_lower, ind_upper = get_wave_ind(wave_grid, wave_min, wave_max)
        wave_bins = wave_grid[ind_lower:ind_upper + 1]
        loglam_bins = np.log10(wave_bins)
    else:
        msgs.error('You must input either a uniformly space loglam grid or wave grid')

    if weights is None:
        msgs.info('No weights were provided. Using uniform weights.')
        weights = np.ones(nimgs) / float(nimgs)

    # Create the slit_cen_stack and determine the minimum and maximum spatial offsets that we need to cover to determine
    # the spatial bins
    spat_img = np.outer(np.ones(nspec), np.arange(nspat))
    dspat_stack = np.zeros_like(sciimg_stack)
    spat_min = np.inf
    spat_max = -np.inf
    for img in range(nimgs):
        # center of the slit replicated spatially
        slit_cen_img = np.outer(trace_stack[img, :], np.ones(nspat))
        dspat_iexp = (spat_img - slit_cen_img)
        dspat_stack[img, :, :] = dspat_iexp
        thismask_now = thismask_stack[img, :, :]
        spat_min = np.fmin(spat_min, dspat_iexp[thismask_now].min())
        spat_max = np.fmax(spat_max, dspat_iexp[thismask_now].max())

    spat_min_int = int(np.floor(spat_min))
    spat_max_int = int(np.ceil(spat_max))
    dspat_bins = np.arange(spat_min_int, spat_max_int + 1, 1)

    sci_list = [sciimg_stack, sciimg_stack - skymodel_stack, tilts_stack, waveimg_stack, dspat_stack]
    var_list = [utils.calc_ivar(sciivar_stack)]

    sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = rebin2d(
        wave_bins, dspat_bins, waveimg_stack, dspat_stack, thismask_stack, inmask_stack, sci_list, var_list)

    # Now compute the final stack with sigma clipping
    sigrej = 3.0
    maxiters = 10
    sci_list_out, var_list_out, outmask, nused = weighted_combine(
        weights, sci_list_rebin, var_list_rebin, (norm_rebin_stack != 0),
        sigma_clip=True, sigma_clip_stack=sci_list_rebin[1], sigrej=sigrej, maxiters=maxiters)
    sciimg, imgminsky, tilts, waveimg, dspat = sci_list_out
    sciivar = utils.calc_ivar(var_list_out[0])

    nspec_coadd, nspat_coadd = imgminsky.shape

    # Now let's try to do an extraction
    thismask = np.ones_like(imgminsky, dtype=bool)
    slit_left = np.full(nspec_coadd, 0.0)
    slit_righ = np.full(nspec_coadd, nspat_coadd - 1)
    tslits_dict = dict(lcen=slit_left, rcen=slit_righ,
                             nspec=nspec_coadd, nspat=nspat_coadd,
                             pad=0, binning=None, spectrograph=None)

    return sciimg, sciivar, imgminsky, outmask, nused, tilts, waveimg, dspat, thismask, tslits_dict


def img_list_error_check(sci_list, var_list):
    """
    Utility routine for dealing dealing with lists of image stacks for rebin2d and weigthed_combine routines below. This
    routine checks that the images sizes are correct and routines the shape of the image stacks.
    Args:
        sci_list: list
            List of  float ndarray images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with the  weights, inmask_stack, and possibly sigma clipping
        var_list: list
            List of  float ndarray variance images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with proper erorr propagation, i.e.
            using the  weights**2, inmask_stack, and possibly sigma clipping

    Returns:
        shape: tuple
            The shapes of the image stacks, (nimgs, nspec, nspat)

    """
    shape_sci_list = []
    for img in sci_list:
        shape_sci_list.append(img.shape)
        if img.ndim != 3:
            msgs.error('Dimensionality of an image in sci_list is not 3')

    shape_var_list = []
    for img in var_list:
        shape_var_list.append(img.shape)
        if img.ndim != 3:
            msgs.error('Dimensionality of an image in var_list is not 3')

    for isci in shape_sci_list:
        if isci != shape_sci_list[0]:
            msgs.error('An image in sci_list have different dimensions')
        for ivar in shape_var_list:
            if ivar != shape_var_list[0]:
                msgs.error('An image in var_list have different dimensions')
            if isci != ivar:
                msgs.error('An image in sci_list had different dimensions than an image in var_list')

    shape = shape_sci_list[0]

    return shape

def rebin2d(spec_bins, spat_bins, waveimg_stack, spatimg_stack, thismask_stack, inmask_stack, sci_list, var_list):
    """
    Rebin a set of images and propagate variance onto a new spectral and spatial grid. This routine effectively
    performs "recitifies" images using np.histogram2d which is extremely fast and effectiveluy performs
    nearest grid point interpolation.

    Args:
        spec_bins: float ndarray, shape = (nspec_rebin)
           Spectral bins to rebin to.
        spat_bins: float ndarray, shape = (nspat_rebin)
           Spatial bins to rebin to.
        waveimg_stack: float ndarray, shape = (nimgs, nspec, nspat)
            Stack of nimgs wavelength images with shape = (nspec, nspat) each
        spatimg_stack: float ndarray, shape = (nimgs, nspec, nspat)
            Stack of nimgs spatial position images with shape = (nspec, nspat) each
        thismask_stack: bool ndarray, shape = (nimgs, nspec, nspat)
            Stack of nimgs images with shape = (nspec, nspat) indicating the locatons on the pixels on an image that
            are on the slit in question.
        inmask_stack: bool ndarray, shape = (nimgs, nspec, nspat)
            Stack of nimgs images with shape = (nspec, nspat) indicating which pixels on an image are masked.
            True = Good, False = Bad
        sci_list: list
            List of  float ndarray images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be rebinned onto the new spec_bins, spat_bins
        var_list: list
            List of  float ndarray variance images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be rebbinned with proper erorr propagation

    Returns:
        sci_list_out: list
           The list of ndarray rebinned images with new shape (nimgs, nspec_rebin, nspat_rebin)
        var_list_out: list
           The list of ndarray rebinned variance images with correct error propagation with shape
           (nimgs, nspec_rebin, nspat_rebin)
        norm_rebin_stack: int ndarray, shape (nimgs, nspec_rebin, nspat_rebin)
           An image stack indicating the integer occupation number of a given pixel. In other words, this number would be zero
           for empty bins, one for bins that were populated by a single pixel, etc. This image takes the input
           inmask_stack into account. The output mask for each image can be formed via
           outmask_rebin_satck = (norm_rebin_stack > 0)
        nsmp_rebin_stack: int ndarray, shape (nimgs, nspec_rebin, nspat_rebin)
           An image stack indicating the integer occupation number of a given pixel taking only the thismask_stack into
           account, but taking the inmask_stack into account. This image is mainly constructed for bookeeping purposes,
           as it represents the number of times each pixel in the rebin image was populated taking only the "geometry"
           of the rebinning into account (i.e. the thismask_stack), but not the masking (inmask_stack).

    """

    shape = img_list_error_check(sci_list, var_list)
    nimgs = shape[0]
    # allocate the output mages
    nspec_rebin = spec_bins.size - 1
    nspat_rebin = spat_bins.size - 1
    shape_out = (nimgs, nspec_rebin, nspat_rebin)
    nsmp_rebin_stack = np.zeros(shape_out)
    norm_rebin_stack = np.zeros(shape_out)
    sci_list_out = []
    for ii in range(len(sci_list)):
        sci_list_out.append(np.zeros(shape_out))
    var_list_out = []
    for jj in range(len(var_list)):
        var_list_out.append(np.zeros(shape_out))

    for img in range(nimgs):
        # This fist image is purely for bookeeping purposes to determine the number of times each pixel
        # could have been sampled
        thismask = thismask_stack[img, :, :]
        spec_rebin_this = waveimg_stack[img, :, :][thismask]
        spat_rebin_this = spatimg_stack[img, :, :][thismask]

        nsmp_rebin_stack[img, :, :], spec_edges, spat_edges = np.histogram2d(spec_rebin_this, spat_rebin_this,
                                                               bins=[spec_bins, spat_bins], density=False)

        finmask = thismask & inmask_stack[img,:,:]
        spec_rebin = waveimg_stack[img, :, :][finmask]
        spat_rebin = spatimg_stack[img, :, :][finmask]
        norm_img, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                          bins=[spec_bins, spat_bins], density=False)
        norm_rebin_stack[img, :, :] = norm_img

        # Rebin the science images
        for indx, sci in enumerate(sci_list):
            weigh_sci, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                               bins=[spec_bins, spat_bins], density=False,
                                                               weights=sci[img,:,:][finmask])
            sci_list_out[indx][img, :, :] = (norm_img > 0.0) * weigh_sci/(norm_img + (norm_img == 0.0))

        # Rebin the variance images, note the norm_img**2 factor for correct error propagation
        for indx, var in enumerate(var_list):
            weigh_var, spec_edges, spat_edges = np.histogram2d(spec_rebin, spat_rebin,
                                                               bins=[spec_bins, spat_bins], density=False,
                                                               weights=var[img, :, :][finmask])
            var_list_out[indx][img, :, :] = (norm_img > 0.0)*weigh_var/(norm_img + (norm_img == 0.0))**2



    return sci_list_out, var_list_out, norm_rebin_stack.astype(int), nsmp_rebin_stack.astype(int)

# TODO make weights optional and do uniform weighting without.
def weighted_combine(weights, sci_list, var_list, inmask_stack,
                     sigma_clip=False, sigma_clip_stack = None, sigrej=None, maxiters=5):
    """

    Args:
        weights: ndarray, float shape (nimgs)
        sci_list: list
            List of  float ndarray images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with the  weights, inmask_stack, and possibly sigma clipping
        var_list: list
            List of  float ndarray variance images (each being an image stack with shape (nimgs, nspec, nspat))
            which are to be combined with proper erorr propagation, i.e.
            using the  weights**2, inmask_stack, and possibly sigma clipping
        inmask_stack: ndarray, boolean, shape (nimgs, nspec, nspat)
            Array of input masks for the images. True = Good, False=Bad
        sigma_clip: bool, default = False
            Combine with a mask by sigma clipping the image stack. Only valid if nimgs > 2
        sigma_clip_stack: ndarray, float, shape (nimgs, nspec, nspat), default = None
            The image stack to be used for the sigma clipping. For example if
            if the list of images to be combined with the weights is [sciimg_stack, waveimg_stack, tilts_stack] you
            would be sigma clipping with sciimg_stack, and would set sigma_clip_stack = sciimg_stack
        sigrej: int or float, default = None
            Rejection threshold for sigma clipping. Code defaults to determining this automatically based
            on the numberr of images provided.
        maxiters:
            Maximum number of iterations for sigma clipping using astropy.stats.SigmaClip

    Returns:
        sci_list_out: list
           The list of ndarray float combined images with shape (nspec, nspat)
        var_list_out: list
           The list of ndarray propagated variance images with shape (nspec, nspat)
        outmask: bool ndarray, shape (nspec, nspat)
           Mask for combined image. True=Good, False=Bad
        nused: int ndarray, shape (nspec, nspat)
           Image of integers indicating the number of images that contributed to each pixel
    """

    shape = img_list_error_check(sci_list, var_list)
    nimgs = shape[0]
    nspec = shape[1]
    nspat = shape[2]

    if nimgs == 1:
        # If only one image is passed in, simply return the input lists of images, but reshaped
        # to be (nspec, nspat)
        msgs.warn('Cannot combine a single image. Returning input images')
        sci_list_out = []
        for sci_stack in sci_list:
            sci_list_out.append(sci_stack.reshape((nspec, nspat)))
        var_list_out = []
        for var_stack in var_list:
            var_list_out.append(var_stack.reshape((nspec, nspat)))
        outmask = inmask_stack.reshape((nspec, nspat))
        nused = outmask.astype(int)
        return sci_list_out, var_list_out, outmask, nused

    if sigma_clip and nimgs >= 3:
        if sigma_clip_stack is None:
            msgs.error('You must specify sigma_clip_stack, i.e. which quantity to use for sigma clipping')
        if sigrej is None:
            if nimgs <= 2:
                sigrej = 100.0  # Irrelevant for only 1 or 2 files, we don't sigma clip below
            elif nimgs == 3:
                sigrej = 1.1
            elif nimgs == 4:
                sigrej = 1.3
            elif nimgs == 5:
                sigrej = 1.6
            elif nimgs == 6:
                sigrej = 1.9
            else:
                sigrej = 2.0
        # sigma clip if we have enough images
        # mask_stack > 0 is a masked value. numpy masked arrays are True for masked (bad) values
        data = np.ma.MaskedArray(sigma_clip_stack, np.invert(inmask_stack))
        sigclip = astropy.stats.SigmaClip(sigma=sigrej, maxiters=maxiters, cenfunc='median')
        data_clipped = sigclip(data, axis=0, masked=True)
        mask_stack = np.invert(data_clipped.mask)  # mask_stack = True are good values
    else:
        if sigma_clip and nimgs < 3:
            msgs.warn('Sigma clipping requested, but you cannot sigma clip with less than 3 images. '
                      'Proceeding without sigma clipping')
        mask_stack = inmask_stack  # mask_stack = True are good values

    nused = np.sum(mask_stack, axis=0)
    weights_stack = np.einsum('i,ijk->ijk', weights, mask_stack)
    weights_sum = np.sum(weights_stack, axis=0)
    sci_list_out = []
    for sci_stack in sci_list:
        sci_list_out.append(np.sum(sci_stack*weights_stack, axis=0)/(weights_sum + (weights_sum == 0.0)))
    var_list_out = []
    for var_stack in var_list:
        var_list_out.append(np.sum(var_stack * weights_stack**2, axis=0) / (weights_sum + (weights_sum == 0.0))**2)
    # Was it masked everywhere?
    outmask = np.any(mask_stack, axis=0)

    return sci_list_out, var_list_out, outmask, nused
