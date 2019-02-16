""" Module for image processing core methods
"""

import astropy.stats
import numpy as np
import os
import scipy
from astropy.io import fits
from pypeit import msgs
from pypeit import utils
from pypeit import masterframe
from pypeit import waveimage
from pypeit import wavetilts
from pypeit import reduce
from pypeit import processimages

from pypeit.core import load, coadd, pixels
from pypeit.core import parse
from pypeit import traceslits
from pypeit.spectrographs import util


def get_brightest_obj(specobjs_list, echelle=True):
    """
    Utility routine to find the brightest object in each exposure given a specobjs_list. This currently only works
    for echelle.

    Parameters:
        specobjs_list: list
           List of SpecObjs objects.
    Optional Parameters:
        echelle: bool, default=True

    Returns:
    (objid, snr_bar), tuple

        objid: ndarray, int, shape (len(list),)
            Array of object ids representing the brightest object in each exposure
        snr_bar: ndarray, float, shape (len(list),)
            Average S/N over all the orders for this object

    """
    nexp = len(specobjs_list)
    objid = np.zeros(nexp, dtype=int)
    snr_bar = np.zeros(nexp)
    if echelle:
        norders = specobjs_list[0].ech_orderindx.max() + 1
        for iexp, sobjs in enumerate(specobjs_list):
            uni_objid = np.unique(sobjs.ech_objid)
            nobjs = len(uni_objid)
            order_snr = np.zeros((norders, nobjs))
            for iord in range(norders):
                for iobj in range(nobjs):
                    ind = (sobjs.ech_orderindx == iord) & (sobjs.ech_objid == uni_objid[iobj])
                    flux = sobjs[ind][0].optimal['COUNTS']
                    sig = sobjs[ind][0].optimal['COUNTS_SIG']
                    wave = sobjs[ind][0].optimal['WAVE']
                    mask = sobjs[ind][0].optimal['MASK']
                    rms_sn, weights = coadd.sn_weights(flux, sig, mask, wave, const_weights=True)
                    order_snr[iord, iobj] = rms_sn

            # Compute the average SNR and find the brightest object
            snr_bar_vec = np.mean(order_snr, axis=0)
            objid[iexp] = uni_objid[snr_bar_vec.argmax()]
            snr_bar[iexp] = snr_bar_vec[snr_bar_vec.argmax()]
    else:
        msgs.error('Automated objid selection is not yet implemented for multislit')
        # -- for multislit this routine will:
        #  Determine which object on the slitmask is brightest and attempt to match them up somehow.  The problem is
        #  how to associate the objects together on the slits if they have moved.

    return objid, snr_bar

def optimal_weights(specobjs_list, slitid, objid):
    """
    Determine optimal weights for a 2d coadds. This script grabs the information from SpecObjs list for the
    object with specified slitid and objid and passes to coadd.sn_weights to determine the optimal weights for
    each exposure. This routine will also pass back the trace and the wavelengths (optimally extracted) for each
    exposure.

    Args:
        specobjs_list: list
           list of SpecObjs objects contaning the objects that were extracted from each frame that will contribute
           to the coadd.
        slitid: int
           The slitid that has the brightest object whose S/N will be used to determine the weight for each frame.
        objid: int
           The objid index of the brightest object whose S/N will be used to determine the weight for each frame.

    Returns:
        (rms_sn, weights, trace_stack, wave_stack)
        rms_sn : ndarray, shape = (len(specobjs_list),)
            Root mean square S/N value for each input spectra
        weights : ndarray, shape (len(specobjs_list),)
            Weights to be applied to the spectra. These are signal-to-noise squared weights.
        trace_stack: ndarray, shape = (len(specobs_list), nspec)
            Traces for each exposure
        wave_stack: ndarray, shape = (len(specobs_list), nspec)
            Wavelengths (optimally extracted) for each exposure.

    """

    nexp = len(specobjs_list)
    nspec = specobjs_list[0][0].trace_spat.shape[0]
    # Grab the traces, flux, wavelength and noise for this slit and objid.
    trace_stack = np.zeros((nexp, nspec), dtype=float)
    flux_stack = np.zeros((nexp, nspec), dtype=float)
    sig_stack = np.zeros((nexp, nspec), dtype=float)
    wave_stack = np.zeros((nexp, nspec), dtype=float)
    mask_stack = np.zeros((nexp, nspec), dtype=bool)

    for iexp, sobjs in enumerate(specobjs_list):
        ithis = (sobjs.slitid == slitid) & (sobjs.objid == objid[iexp])
        trace_stack[iexp, :] = sobjs[ithis].trace_spat
        flux_stack[iexp,:] = sobjs[ithis][0].optimal['COUNTS']
        sig_stack[iexp,:] = sobjs[ithis][0].optimal['COUNTS_SIG']
        wave_stack[iexp,:] = sobjs[ithis][0].optimal['WAVE']
        mask_stack[iexp,:] = sobjs[ithis][0].optimal['MASK']

    # TODO For now just use the zero as the reference for the wavelengths? Perhaps we should be rebinning the data though?
    rms_sn, weights = coadd.sn_weights(flux_stack, sig_stack, mask_stack, wave_stack)

    return rms_sn, weights, trace_stack, wave_stack

def det_error_msg(exten, sdet):
    # Print out error message if extension is not found
    msgs.error("Extension {:s} for requested detector {:s} was not found.\n".format(exten)  +
               " Maybe you chose the wrong detector to coadd? "
               "Set with --det= or check file contents with pypeit_show_2dspec Science/spec2d_XXX --list".format(sdet))

def load_coadd2d_stacks(spec2d_files, det):
    """

    Args:
        spec2d_files: list
           List of spec2d filenames
        det: int
           detector in question

    Returns:
        stack_dict: dict
           Dictionary containing all the images and keys required for perfomring 2d coadds.

    """

    # Get the detector string
    sdet = parse.get_dnum(det, prefix=False)

    # Grab the files
    spec1d_files = [files.replace('spec2d', 'spec1d') for files in spec2d_files]
    # Get the master dir
    head0 = fits.getheader(spec2d_files[0])
    master_dir = os.path.basename(head0['PYPMFDIR'])+'/'
    redux_path =  os.cwd()
    master_path = redux_path + master_dir
    tiltfiles = []
    waveimgfiles = []
    tracefiles = []

    head2d_list=[]
    for file in spec2d_files:
        head = fits.getheader(file)
        head2d_list.append(head)
        trace_key = '{:s}'.format(head['TRACMKEY']) + '_{:02d}'.format(det)
        wave_key = '{:s}'.format(head['ARCMKEY']) + '_{:02d}'.format(det)
        tracefiles.append(masterframe.master_name('trace', trace_key, master_path))
        waveimgfiles.append(masterframe.master_name('wave', wave_key, master_path))
        tiltfiles.append(masterframe.master_name('tilts', wave_key, master_path))

    nfiles = len(spec2d_files)

    specobjs_list = []
    head1d_list=[]
    # TODO Sort this out with the correct detector extensions etc.
    # Read in the image stacks
    for ifile in range(nfiles):
        waveimg = waveimage.load_waveimage(waveimgfiles[ifile])
        tilts = wavetilts.load_tilts(tiltfiles[ifile])
        hdu = fits.open(spec2d_files[ifile])
        # One detector, sky sub for now
        names = [hdu[i].name for i in range(len(hdu))]
        # science image
        try:
            exten = names.index('DET{:s}-PROCESSED'.format(sdet))
        except:  # Backwards compatability
            det_error_msg(exten, sdet)
        sciimg = hdu[exten].data
        # skymodel
        try:
            exten = names.index('DET{:s}-SKY'.format(sdet))
        except:  # Backwards compatability
            det_error_msg(exten, sdet)
        skymodel = hdu[exten].data
        # Inverse variance model
        try:
            exten = names.index('DET{:s}-IVARMODEL'.format(sdet))
        except ValueError:  # Backwards compatability
            det_error_msg(exten, sdet)
        sciivar = hdu[exten].data
        # Mask
        try:
            exten = names.index('DET{:s}-MASK'.format(sdet))
        except ValueError:  # Backwards compatability
            det_error_msg(exten, sdet)
        mask = hdu[exten].data
        if ifile == 0:
            # the two shapes accomodate the possibility that waveimg and tilts are binned differently
            shape_wave = (nfiles,waveimg.shape[0],waveimg.shape[1])
            shape_sci = (nfiles,sciimg.shape[0],sciimg.shape[1])
            waveimg_stack = np.zeros(shape_wave,dtype=float)
            tilts_stack = np.zeros(shape_wave,dtype=float)
            sciimg_stack = np.zeros(shape_sci,dtype=float)
            skymodel_stack = np.zeros(shape_sci,dtype=float)
            sciivar_stack = np.zeros(shape_sci,dtype=float)
            mask_stack = np.zeros(shape_sci,dtype=float)

        waveimg_stack[ifile,:,:] = waveimg
        tilts_stack[ifile,:,:] = tilts
        sciimg_stack[ifile,:,:] = sciimg
        sciivar_stack[ifile,:,:] = sciivar
        mask_stack[ifile,:,:] = mask
        skymodel_stack[ifile,:,:] = skymodel

        sobjs, head = load.load_specobjs(spec1d_files[ifile])
        head1d_list.append(head)
        specobjs_list.append(sobjs)

    # Right now we assume there is a single tslits_dict for all images and read in the first one
    # TODO this needs to become a tslits_dict for each file to accomodate slits defined by flats taken on different
    # nights
    tslits_dict = traceslits.load_tslits_dict(tracefiles[0])
    spectrograph = util.load_spectrograph(tslits_dict['spectrograph'])
    slitmask = pixels.tslits2mask(tslits_dict)
    slitmask_stack = np.einsum('i,jk->ijk', np.ones(nfiles), slitmask)

    # Fill the master key dict
    head2d = head2d_list[0]
    master_key_dict = {}
    master_key_dict['frame'] = head2d['FRAMMKEY']  + '_{:02d}'.format(det)
    master_key_dict['bpm']   = head2d['BPMMKEY']   + '_{:02d}'.format(det)
    master_key_dict['bias']  = head2d['BIASMKEY']  + '_{:02d}'.format(det)
    master_key_dict['arc']   = head2d['ARCMKEY']   + '_{:02d}'.format(det)
    master_key_dict['trace'] = head2d['TRACMKEY']  + '_{:02d}'.format(det)
    master_key_dict['flat']  = head2d['FLATMKEY']  + '_{:02d}'.format(det)
    stack_dict = dict(specobjs_list=specobjs_list, tslits_dict=tslits_dict,
                      slitmask_stack=slitmask_stack,
                      sciimg_stack=sciimg_stack, sciivar_stack=sciivar_stack,
                      skymodel_stack=skymodel_stack, mask_stack=mask_stack,
                      tilts_stack=tilts_stack, waveimg_stack=waveimg_stack,
                      head1d_list = head1d_list, head2d_list=head2d_list,
                      redux_path=redux_path, master_path=master_path, master_dir=master_dir,
                      master_key_dict=master_key_dict,
                      spectrograph = tslits_dict['spectrograph'])

    return stack_dict




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

def broadcast_weights(weights, shape):
    """
    Utility routine to broadcast weights to be the size of image stacks specified by shape
    Args:
        weights: float ndarray of weights.
            Options for the shape of weights are:
                (nimgs,)              -- a single weight per image
                (nimgs, nspec)        -- wavelength dependent weights per image
                (nimgs, nspec, nspat) -- weights already have the shape of the image stack and are simply returned
        shape: tuple of integers
            Shape of the image stacks for weighted coadding (nimgs, nspec, nspat)

    Returns:

    """
    # Create the weights stack images from the wavelength dependent weights, i.e. propagate these
    # weights to the spatial direction
    if weights.ndim == 1:
        # One float per image
        weights_stack = np.einsum('i,ijk->ijk', weights, np.ones(shape))
    elif weights.ndim == 2:
        # Wavelength dependent weights per image
        weights_stack = np.einsum('ij,k->ijk', weights, np.ones(shape[2]))
    elif weights.ndim == 3:
        # Full image stack of weights
        if weights.shape != shape:
            msgs.error('The shape of weights does not match the shape of the image stack')
        weights_stack = weights
    else:
        msgs.error('Unrecognized dimensionality for weights')

    return weights_stack

def coadd2d(trace_stack, sciimg_stack, sciivar_stack, skymodel_stack, inmask_stack, tilts_stack, waveimg_stack,
            thismask_stack, weights=None, loglam_grid=None, wave_grid=None):
    """
    This routine will perform a 2d co-add of a stack of PypeIt spec2d reduction outputs. The slit will be
    'rectified' onto a spatial and spectral grid, which encompasses the spectral and spatial coverage of the image stacks.
    The rectification uses nearest grid point interpolation to avoid covariant errors.
    Dithering is supported as all images are centered relative to a set of reference traces in trace_stack.

    Args:
    -----
        trace_stack: float ndarray, shape (nimgs, nspec)
           Stack of reference traces about which the images are rectified and coadded.
           If the images were not dithered then this reference trace can simply be the center
           of the slit slitcen = (slit_left + slit_righ)/2. If the images were dithered, then this could trace_stack
           could either be the slitcen appropriately shifted with the dither pattern, or it could the trace of the object
           of interest in each exposure determined by running PypeIt on the individual images.
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
        thismask_stack: bool ndarray, shape = (nimgs, nspec, nspat)
           Stack of masks indicating which pixels are on the slit in question. True = On slit, False=Off slit.

    Optional Args:
    --------------
        weights: float ndarray, shape = (nimgs,), (nimgs, nspec), or (nimgs, nspec, nspat), default = None
            Set of weights used for the weighted combined of the rectified images using weighted_combine.
            If the weights are not provided then a uniform weighting will be adopted. Weights will be broadast
            to the correct size of the image stacks using the broadcast_weights utility function.
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
    wavemask = thismask_stack & (waveimg_stack > 1.0) # TODO This cut on waveimg_stack should not be necessary
    wave_lower = waveimg_stack[wavemask].min()
    wave_upper = waveimg_stack[wavemask].max()
    if loglam_grid is not None:
        ind_lower, ind_upper = get_wave_ind(loglam_grid, np.log10(wave_lower), np.log10(wave_upper))
        loglam_bins = loglam_grid[ind_lower:ind_upper + 1]
        wave_bins = np.power(10.0, loglam_bins)
    elif wave_grid is not None:
        ind_lower, ind_upper = get_wave_ind(wave_grid, wave_lower, wave_upper)
        wave_bins = wave_grid[ind_lower:ind_upper + 1]
        loglam_bins = np.log10(wave_bins)
    else:
        msgs.error('You must input either a uniformly space loglam grid or wave grid')

    if weights is None:
        msgs.info('No weights were provided. Using uniform weights.')
        weights = np.ones(nimgs)/float(nimgs)

    weights_stack = broadcast_weights(weights, sciimg_stack.shape)

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

    sci_list = [weights_stack, sciimg_stack, sciimg_stack - skymodel_stack, tilts_stack, waveimg_stack, dspat_stack]
    var_list = [utils.calc_ivar(sciivar_stack)]

    sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = rebin2d(
        wave_bins, dspat_bins, waveimg_stack, dspat_stack, thismask_stack, inmask_stack, sci_list, var_list)

    # Now compute the final stack with sigma clipping
    sigrej = 3.0
    maxiters = 10
    # sci_list_rebin[0] = rebinned weights image stack
    # sci_list_rebin[1:] = stacks of images that we want to weighted combine
    # sci_list_rebin[2] = rebinned sciimg-sky_model images that we used for the sigma clipping
    sci_list_out, var_list_out, outmask, nused = weighted_combine(
        sci_list_rebin[0], sci_list_rebin[1:], var_list_rebin, (norm_rebin_stack != 0),
        sigma_clip=True, sigma_clip_stack=sci_list_rebin[2], sigrej=sigrej, maxiters=maxiters)
    sciimg, imgminsky, tilts, waveimg, dspat = sci_list_out
    sciivar = utils.calc_ivar(var_list_out[0])

    # Compute the midpoints vectors, and lower/upper bins of the rectified image
    wave_mid = ((wave_bins + np.roll(wave_bins,1))/2.0)[1:]
    wave_min = wave_bins[:-1]
    wave_max = wave_bins[1:]
    dspat_mid = ((dspat_bins + np.roll(dspat_bins,1))/2.0)[1:]

    # Interpolate the dspat images wherever the coadds are masked because a given pixel was not sampled. This is done
    # because the dspat image is not allowed to have holes if it is going to work with local_skysub_extract
    nspec_coadd, nspat_coadd = imgminsky.shape
    spat_img_coadd, spec_img_coadd = np.meshgrid(np.arange(nspat_coadd), np.arange(nspec_coadd))
    if np.any(np.invert(outmask)):
        points_good = np.stack((spec_img_coadd[outmask], spat_img_coadd[outmask]), axis=1)
        points_bad = np.stack((spec_img_coadd[np.invert(outmask)], spat_img_coadd[np.invert(outmask)]), axis=1)
        values_dspat = dspat[outmask]
        dspat_bad = scipy.interpolate.griddata(points_good, values_dspat, points_bad, method='cubic')
        dspat[np.invert(outmask)] = dspat_bad
        # Points outside the convex hull of the data are set to nan. We identify those and simply assume them values from
        # the dspat_img_fake, which is what dspat would be on a regular perfectly rectified image grid.
        nanpix = np.isnan(dspat)
        if np.any(nanpix):
            dspat_img_fake = spat_img_coadd + dspat_mid[0]
            dspat[nanpix] = dspat_img_fake[nanpix]

    coadd_dict = dict(wave_bins=wave_bins, dspat_bins=dspat_bins,
                      wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max, dspat_mid=dspat_mid,
                      sciimg=sciimg, sciivar=sciivar, imgminsky=imgminsky, outmask=outmask,
                      nused=nused, tilts=tilts, waveimg=waveimg, dspat=dspat,
                      nspec=imgminsky.shape[0], nspat=imgminsky.shape[1])
    return coadd_dict


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

# TODO Break up into separate methods?
def extract_coadd2d(stack_dict, master_dir, ir_redux=False, par=None, show=False, show_peaks=False):
    """
    Main routine to run the extraction for 2d coadds. This performs 2d coadd specific tasks, and then also performs
    some of the tasks analogous to the pypeit.extract_one method. Docs coming soon....
    Args:
        stack_dict:
        master_dir:
        ir_redux:
        par:
        show:
        show_peaks:

    Returns:

    """

    # Find the objid of the brighest object, and the average snr across all orders
    nslits = stack_dict['tslits_dict']['slit_left'].shape[1]
    objid, snr_bar = get_brightest_obj(stack_dict['specobjs_list'], echelle=True)
    # TODO Print out a report here on the image stack, i.e. S/N of each image

    spectrograph = util.load_spectrograph(stack_dict['spectrograph'])
    par = spectrograph.default_pypeit_par() if par is None else par

    # Grab the wavelength grid that we will rectify onto
    wave_grid = spectrograph.wavegrid()

    coadd_list = []
    nspec_vec = np.zeros(nslits,dtype=int)
    nspat_vec = np.zeros(nslits,dtype=int)
    # ToDO Generalize this to be a loop over detectors, such tha the coadd_list is an ordered dict (perhaps) with
    # all the slits on all detectors
    for islit in range(nslits):
        msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(islit,nslits-1))
        # Determine the wavelength dependent optimal weights and grab the reference trace
        rms_sn, weights, trace_stack, wave_stack = optimal_weights(stack_dict['specobjs_list'], islit, objid)
        thismask_stack = stack_dict['slitmask_stack'] == islit
        # Perform the 2d coadd
        coadd_dict = coadd2d(trace_stack, stack_dict['sciimg_stack'], stack_dict['sciivar_stack'],
                             stack_dict['skymodel_stack'], (stack_dict['mask_stack'] == 0),
                             stack_dict['tilts_stack'], stack_dict['waveimg_stack'], thismask_stack,
                             weights = weights, wave_grid=wave_grid)
        coadd_list.append(coadd_dict)
        nspec_vec[islit]=coadd_dict['nspec']
        nspat_vec[islit]=coadd_dict['nspat']

    # Determine the size of the psuedo image
    nspat_pad = 10
    nspec_psuedo = nspec_vec.max()
    nspat_psuedo = np.sum(nspat_vec) + (nslits + 1)*nspat_pad
    spec_vec_psuedo = np.arange(nspec_psuedo)
    shape_psuedo = (nspec_psuedo, nspat_psuedo)
    imgminsky_psuedo = np.zeros(shape_psuedo)
    sciivar_psuedo = np.zeros(shape_psuedo)
    waveimg_psuedo = np.zeros(shape_psuedo)
    tilts_psuedo = np.zeros(shape_psuedo)
    spat_psuedo = np.zeros(shape_psuedo)
    nused_psuedo = np.zeros(shape_psuedo, dtype=int)
    inmask_psuedo = np.zeros(shape_psuedo, dtype=bool)
    wave_mid = np.zeros((nspec_psuedo, nslits))
    wave_min = np.zeros((nspec_psuedo, nslits))
    wave_max = np.zeros((nspec_psuedo, nslits))
    dspat_mid = np.zeros((nspat_psuedo, nslits))

    spat_left = nspat_pad
    slit_left = np.zeros((nspec_psuedo, nslits))
    slit_righ = np.zeros((nspec_psuedo, nslits))
    spec_min1 = np.zeros(nslits)
    spec_max1 = np.zeros(nslits)
    for islit, coadd_dict in enumerate(coadd_list):
        spat_righ = spat_left + nspat_vec[islit]
        ispec = slice(0,nspec_vec[islit])
        ispat = slice(spat_left,spat_righ)
        imgminsky_psuedo[ispec, ispat] = coadd_dict['imgminsky']
        sciivar_psuedo[ispec, ispat] = coadd_dict['sciivar']
        waveimg_psuedo[ispec, ispat] = coadd_dict['waveimg']
        tilts_psuedo[ispec, ispat] = coadd_dict['tilts']
        # spat_psuedo is the sub-pixel image position on the rebinned psuedo image
        inmask_psuedo[ispec, ispat] = coadd_dict['outmask']
        image_temp = (coadd_dict['dspat'] -  coadd_dict['dspat_mid'][0] + spat_left)*coadd_dict['outmask']
        spat_psuedo[ispec, ispat] = image_temp
        nused_psuedo[ispec, ispat] = coadd_dict['nused']
        wave_mid[ispec, islit] = coadd_dict['wave_mid']
        wave_min[ispec, islit] = coadd_dict['wave_min']
        wave_max[ispec, islit] = coadd_dict['wave_max']

        dspat_mid[ispat, islit] = coadd_dict['dspat_mid']
        slit_left[:,islit] = np.full(nspec_psuedo, spat_left)
        slit_righ[:,islit] = np.full(nspec_psuedo, spat_righ)
        spec_max1[islit] = nspec_vec[islit]-1
        spat_left = spat_righ + nspat_pad

    tslits_dict_psuedo = dict(slit_left=slit_left, slit_righ=slit_righ, nspec=nspec_psuedo, nspat=nspat_psuedo, pad=0,
                              nslits = nslits, binspectral=1, binspatial=1, spectrograph=spectrograph.spectrograph,
                              spec_min=spec_min1, spec_max=spec_max1)

    slitmask_psuedo = pixels.tslits2mask(tslits_dict_psuedo)
    # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
    # which wreaks havoc on the local sky-subtraction
    min_slit_frac = 0.70
    spec_min = np.zeros(nslits)
    spec_max = np.zeros(nslits)
    for islit in range(nslits):
        slit_width = np.sum(inmask_psuedo*(slitmask_psuedo == islit),axis=1)
        slit_width_img = np.outer(slit_width, np.ones(nspat_psuedo))
        med_slit_width = np.median(slit_width_img[slitmask_psuedo==islit])
        nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
        nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
        slit_width_sm = scipy.ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
        igood = (slit_width_sm > min_slit_frac*med_slit_width)
        spec_min[islit] = spec_vec_psuedo[igood].min()
        spec_max[islit] = spec_vec_psuedo[igood].max()
        bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_psuedo == islit)
        inmask_psuedo[bad_pix] = False

    # Update with tslits_dict_psuedo
    tslits_dict_psuedo['spec_min'] = spec_min
    tslits_dict_psuedo['spec_max'] = spec_max
    slitmask_psuedo = pixels.tslits2mask(tslits_dict_psuedo)

    # Make a fake bitmask from the outmask. We are kludging the crmask to be the outmask_psuedo here, and setting the bpm to
    # be good everywhere
    mask = processimages.ProcessImages.build_mask(imgminsky_psuedo, sciivar_psuedo, np.invert(inmask_psuedo),
                                                  np.zeros_like(inmask_psuedo), slitmask=slitmask_psuedo)

    redux = reduce.instantiate_me(spectrograph, tslits_dict_psuedo, mask, ir_redux=ir_redux, par=par, objtype = 'science')

    if show:
        redux.show('image', image=imgminsky_psuedo*(mask == 0), chname = 'imgminsky', slits=True, clear=True)
    # Object finding
    sobjs_obj, nobj, skymask_init = redux.find_objects(imgminsky_psuedo, sciivar_psuedo, ir_redux=ir_redux,
                                                       show_peaks=show_peaks, show=show)
    # Local sky-subtraction
    global_sky_psuedo = np.zeros_like(imgminsky_psuedo) # No global sky for co-adds since we go straight to local
    rn2img_psuedo = global_sky_psuedo # No rn2img for co-adds since we go do not model noise
    skymodel_psuedo, objmodel_psuedo, ivarmodel_psuedo, outmask_psuedo, sobjs = \
        redux.local_skysub_extract(imgminsky_psuedo, sciivar_psuedo, tilts_psuedo, waveimg_psuedo, global_sky_psuedo,
                                   rn2img_psuedo, sobjs_obj, spat_pix=spat_psuedo,
                                   model_noise=False, show_profile=show, show=show)

    if ir_redux:
        sobjs.purge_neg()

    # Add the information about the fixed wavelength grid to the sobjs
    for spec in sobjs:
        spec.boxcar['WAVE_GRID'] = wave_mid[:,spec.slitid]
        spec.boxcar['WAVE_GRID_MIN'] = wave_min[:,spec.slitid]
        spec.boxcar['WAVE_GRID_MAX'] = wave_max[:,spec.slitid]

        spec.optimal['WAVE_GRID'] = wave_mid[:,spec.slitid]
        spec.optimal['WAVE_GRID_MIN'] = wave_min[:,spec.slitid]
        spec.optimal['WAVE_GRID_MAX'] = wave_max[:,spec.slitid]


    # TODO Implement flexure and heliocentric corrections on the single exposure 1d reductions and apply them to the
    # waveimage. Change the data model to accomodate a wavelength model for each image.
    # Using the same implementation as in core/pypeit

    # Write out the psuedo master files to disk
    master_key_dict = stack_dict['master_key_dict']
    waveImage = waveimage.WaveImage(None, None, None, None, None,
                                      master_key=master_key_dict['arc'], master_dir=master_dir)
    waveImage.save_master(waveimg_psuedo)
    traceSlits = traceslits.TraceSlits(None, None, None, master_key=master_key_dict['trace'], master_dir=master_dir)
    traceSlits.save_master(tslits_dict_psuedo)

    return imgminsky_psuedo, sciivar_psuedo, skymodel_psuedo, objmodel_psuedo, ivarmodel_psuedo, outmask_psuedo, sobjs


# TODO make weights optional and do uniform weighting without.
def weighted_combine(weights, sci_list, var_list, inmask_stack,
                     sigma_clip=False, sigma_clip_stack = None, sigrej=None, maxiters=5):
    """

    Args:
        weights: float ndarray of weights.
            Options for the shape of weights are:
                (nimgs,)              -- a single weight per image in the stack
                (nimgs, nspec)        -- wavelength dependent weights per image in the stack
                (nimgs, nspec, nspat) -- weights input with the shape of the image stack

             Note that the weights are distinct from the mask which is dealt with via inmask_stack argument so there
             should not be any weights that are set to zero (although in principle this would still work).

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
    weights_stack = broadcast_weights(weights, shape)
    weights_mask_stack = weights_stack*mask_stack

    weights_sum = np.sum(weights_mask_stack, axis=0)
    sci_list_out = []
    for sci_stack in sci_list:
        sci_list_out.append(np.sum(sci_stack*weights_mask_stack, axis=0)/(weights_sum + (weights_sum == 0.0)))
    var_list_out = []
    for var_stack in var_list:
        var_list_out.append(np.sum(var_stack * weights_mask_stack**2, axis=0) / (weights_sum + (weights_sum == 0.0))**2)
    # Was it masked everywhere?
    outmask = np.any(mask_stack, axis=0)

    return sci_list_out, var_list_out, outmask, nused
