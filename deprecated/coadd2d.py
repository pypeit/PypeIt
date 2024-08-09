"""
Module for performing two-dimensional coaddition of spectra.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os
import copy

from IPython import embed

import numpy as np
import scipy
from matplotlib import pyplot as plt

from astropy.io import fits

from pypeit import msgs
from pypeit import utils
from pypeit import ginga
from pypeit import specobjs
from pypeit.masterframe import MasterFrame
from pypeit.waveimage import WaveImage
from pypeit.wavetilts import WaveTilts
from pypeit import specobjs
from pypeit import slittrace
from pypeit import reduce
from pypeit.core import extract
from pypeit.core import load, coadd, pixels
from pypeit.core import parse
from pypeit.core import combine
from pypeit.images import scienceimage
from pypeit.spectrographs import util
from pypeit import calibrations
from pypeit.par import PypeItPar


# TODO: Is this commented code still needed?
#def reference_trace_stack(slitid, stack_dict, offsets=None, objid=None):
#    """
#    Utility function for determining the reference trace about which 2d coadds are performed.
#    There are two modes of operation to determine the reference trace for the 2d coadd of a given slit/order:
#
#     1) offsets: we stack about the center of the slit for the slit in question with the input offsets added
#     2) ojbid: we stack about the trace ofa reference object for this slit given for each exposure by the input objid
#
#    Either offsets or objid must be provided, but the code will raise an exception if both are provided.
#
#    Args:
#        slitid (int):
#           The slit or order that we are currently considering
#        stack_dict (dict):
#           Dictionary containing all the images and keys required for perfomring 2d coadds.
#        offsets (list or np.ndarray):
#           An array of offsets with the same dimensionality as the nexp, the numer of images being coadded.
#        objid: (list or np.ndarray):
#           An array of objids with the same dimensionality as the nexp, the number of images being coadded.
#
#    Returns:
#        ref_trace_stack
#
#        ref_trace_stack (np.ndarray):
#            An array with shape (nspec, nexp) containing the reference trace for each of the nexp exposures.
#
#    """
#
#    if offsets is not None and objid is not None:
#        msgs.errror('You can only input offsets or an objid, but not both')
#    nexp = len(offsets) if offsets is not None else len(objid)
#    if offsets is not None:
#        tslits_dict_list = stack_dict['tslits_dict_list']
#        nspec, nslits = tslits_dict_list[0]['slit_left'].shape
#        ref_trace_stack = np.zeros((nspec, nexp))
#        for iexp, tslits_dict in enumerate(tslits_dict_list):
#            ref_trace_stack[:, iexp] = (tslits_dict[:, slitid]['slit_left'] + tslits_dict[:, slitid]['slit_righ'])/2.0 + offsets[iexp]
#    elif objid is not None:
#        specobjs_list = stack_dict['specobjs_list']
#        nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]
#        # Grab the traces, flux, wavelength and noise for this slit and objid.
#        ref_trace_stack = np.zeros((nspec, nexp), dtype=float)
#        for iexp, sobjs in enumerate(specobjs_list):
#            # TODO Should be as in optimal_weights
#            ithis = (sobjs.SLITID == slitid) & (sobjs.OBJID == objid[iexp])
#            ref_trace_stack[:, iexp] = sobjs[ithis].TRACE_SPAT
#    else:
#        msgs.error('You must input either offsets or an objid to determine the stack of reference traces')
#
#    return ref_trace_stack

#def optimal_weights(specobjs_list, slitid, objid, sn_smooth_npix, const_weights=False):
#    """
#    Determine optimal weights for 2d coadds. This script grabs the information from SpecObjs list for the
#    object with specified slitid and objid and passes to coadd.sn_weights to determine the optimal weights for
#    each exposure.
#
#    Args:
#        specobjs_list (list):
#           list of SpecObjs objects contaning the objects that were extracted from each frame that will contribute
#           to the coadd.
#        slitid (int):
#           The slitid that has the brightest object whose S/N will be used to determine the weight for each frame.
#        objid (int):
#           The objid index of the brightest object whose S/N will be used to determine the weight for each frame.
#        sn_smooth_npix (float):
#           Number of pixels used for determining smoothly varying S/N ratio weights.
#        const_weights (bool):
#           Use constant weights for coadding the exposures. Default=False
#
#    Returns:
#        rms_sn, weights
#
#        rms_sn : ndarray, shape = (len(specobjs_list),)
#            Root mean square S/N value for each input spectra
#        weights : ndarray, shape (len(specobjs_list),)
#            Weights to be applied to the spectra. These are signal-to-noise squared weights.
#    """
#
#    nexp = len(specobjs_list)
#    nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]
#    # Grab the traces, flux, wavelength and noise for this slit and objid.
#    flux_stack = np.zeros((nspec, nexp), dtype=float)
#    ivar_stack = np.zeros((nspec, nexp), dtype=float)
#    wave_stack = np.zeros((nspec, nexp), dtype=float)
#    mask_stack = np.zeros((nspec, nexp), dtype=bool)
#
#    for iexp, sobjs in enumerate(specobjs_list):
#        embed()
#        try:
#            ithis = (sobjs.SLITID == slitid) & (sobjs.OBJID == objid[iexp])
#        except AttributeError:
#            ithis = (sobjs.ECH_ORDERINDX == slitid) & (sobjs.ECH_OBJID == objid[iexp])
#        try:
#            flux_stack[:,iexp] = sobjs[ithis][0].OPT_COUNTS
#        except:
#            embed(header='104')
#        ivar_stack[:,iexp] = sobjs[ithis][0].OPT_COUNTS_IVAR
#        wave_stack[:,iexp] = sobjs[ithis][0].OPT_WAVE
#        mask_stack[:,iexp] = sobjs[ithis][0].OPT_MASK
#
#    # TODO For now just use the zero as the reference for the wavelengths? Perhaps we should be rebinning the data though?
#    rms_sn, weights = coadd1d.sn_weights(wave_stack, flux_stack, ivar_stack, mask_stack, sn_smooth_npix,
#                                         const_weights=const_weights)
#    return rms_sn, weights.T

def det_error_msg(exten, sdet):
    # Print out error message if extension is not found
    msgs.error("Extension {:s} for requested detector {:s} was not found.\n".format(exten)  +
               " Maybe you chose the wrong detector to coadd? "
               "Set with --det= or check file contents with pypeit_show_2dspec Science/spec2d_XXX --list".format(sdet))


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
        tuple: Returns (ind_lower, ind_upper), Integer lower and upper
        indices into the array wave_grid that cover the interval
        (wave_min, wave_max)
    """

    diff = wave_grid - wave_min
    diff[diff > 0] = np.inf
    if not np.any(diff < 0):
        ind_lower = 0
        msgs.warn('Your wave grid does not extend blue enough. Taking bluest point')
    else:
        ind_lower = np.argmin(np.abs(diff))
    diff = wave_max - wave_grid
    diff[diff > 0] = np.inf
    if not np.any(diff < 0):
        ind_upper = wave_grid.size-1
        msgs.warn('Your wave grid does not extend red enough. Taking reddest point')
    else:
        ind_upper = np.argmin(np.abs(diff))

    return ind_lower, ind_upper



def get_wave_bins(thismask_stack, waveimg_stack, wave_grid):

    # Determine the wavelength grid that we will use for the current slit/order
    # TODO This cut on waveimg_stack should not be necessary
    wavemask = thismask_stack & (waveimg_stack > 1.0)
    wave_lower = waveimg_stack[wavemask].min()
    wave_upper = waveimg_stack[wavemask].max()
    ind_lower, ind_upper = get_wave_ind(wave_grid, wave_lower, wave_upper)
    wave_bins = wave_grid[ind_lower:ind_upper + 1]

    return wave_bins


def get_spat_bins(thismask_stack, trace_stack):

    nimgs, nspec, nspat = thismask_stack.shape
    # Create the slit_cen_stack and determine the minimum and maximum
    # spatial offsets that we need to cover to determine the spatial
    # bins
    spat_img = np.outer(np.ones(nspec), np.arange(nspat))
    dspat_stack = np.zeros_like(thismask_stack,dtype=float)
    spat_min = np.inf
    spat_max = -np.inf
    for img in range(nimgs):
        # center of the slit replicated spatially
        slit_cen_img = np.outer(trace_stack[:, img], np.ones(nspat))
        dspat_iexp = (spat_img - slit_cen_img)
        dspat_stack[img, :, :] = dspat_iexp
        thismask_now = thismask_stack[img, :, :]
        spat_min = np.fmin(spat_min, dspat_iexp[thismask_now].min())
        spat_max = np.fmax(spat_max, dspat_iexp[thismask_now].max())

    spat_min_int = int(np.floor(spat_min))
    spat_max_int = int(np.ceil(spat_max))
    dspat_bins = np.arange(spat_min_int, spat_max_int + 1, 1,dtype=float)

    return dspat_bins, dspat_stack


def compute_coadd2d(ref_trace_stack, sciimg_stack, sciivar_stack, skymodel_stack, inmask_stack, tilts_stack,
                    thismask_stack, waveimg_stack, wave_grid, weights='uniform'):
    """
    Construct a 2d co-add of a stack of PypeIt spec2d reduction outputs.

    Slits are 'rectified' onto a spatial and spectral grid, which
    encompasses the spectral and spatial coverage of the image stacks.
    The rectification uses nearest grid point interpolation to avoid
    covariant errors.  Dithering is supported as all images are centered
    relative to a set of reference traces in trace_stack.

    Args:
        trace_stack (`numpy.ndarray`_):
            Stack of reference traces about which the images are
            rectified and coadded.  If the images were not dithered then
            this reference trace can simply be the center of the slit::

                slitcen = (slit_left + slit_righ)/2

            If the images were dithered, then this object can either be
            the slitcen appropriately shifted with the dither pattern,
            or it could be the trace of the object of interest in each
            exposure determined by running PypeIt on the individual
            images.  Shape is (nimgs, nspec).
        sciimg_stack (`numpy.ndarray`_):
            Stack of science images.  Shape is (nimgs, nspec, nspat).
        sciivar_stack (`numpy.ndarray`_):
            Stack of inverse variance images.  Shape is (nimgs, nspec,
            nspat).
        skymodel_stack (`numpy.ndarray`_):
            Stack of the model sky.  Shape is (nimgs, nspec, nspat).
        inmask_stack (`numpy.ndarray`_):
            Boolean array with the input masks for each image; `True`
            values are *good*, `False` values are *bad*.  Shape is
            (nimgs, nspec, nspat).
        tilts_stack (`numpy.ndarray`_):
           Stack of the wavelength tilts traces.  Shape is (nimgs,
           nspec, nspat).
        waveimg_stack (`numpy.ndarray`_):
           Stack of the wavelength images.  Shape is (nimgs, nspec,
           nspat).
        thismask_stack (`numpy.ndarray`_):
            Boolean array with the masks indicating which pixels are on
            the slit in question.  `True` values are on the slit;
            `False` values are off the slit.  Shape is (nimgs, nspec,
            nspat).
        weights (`numpy.ndarray`_, optional):
            The weights used when combining the rectified images (see
            :func:`weighted_combine`).  If no weights are provided,
            uniform weighting is used.  Weights are broadast to the
            correct size of the image stacks (see
            :func:`broadcast_weights`), as necessary.  Shape must be
            (nimgs,), (nimgs, nspec), or (nimgs, nspec, nspat).
        loglam_grid (`numpy.ndarray`_, optional):
            Wavelength grid in log10(wave) onto which the image stacks
            will be rectified.  The code will automatically choose the
            subset of this grid encompassing the wavelength coverage of
            the image stacks provided (see :func:`waveimg_stack`).
            Either `loglam_grid` or `wave_grid` must be provided.
        wave_grid (`numpy.ndarray`_, optional):
            Same as `loglam_grid` but in angstroms instead of
            log(angstroms). (TODO: Check units...)

    Returns:
        tuple: Returns the following (TODO: This needs to be updated):
            - sciimg: float ndarray shape = (nspec_coadd, nspat_coadd):
              Rectified and coadded science image
            - sciivar: float ndarray shape = (nspec_coadd, nspat_coadd):
              Rectified and coadded inverse variance image with correct
              error propagation
            - imgminsky: float ndarray shape = (nspec_coadd,
              nspat_coadd): Rectified and coadded sky subtracted image
            - outmask: bool ndarray shape = (nspec_coadd, nspat_coadd):
              Output mask for rectified and coadded images. True = Good,
              False=Bad.
            - nused: int ndarray shape = (nspec_coadd, nspat_coadd):
              Image of integers indicating the number of images from the
              image stack that contributed to each pixel
            - tilts: float ndarray shape = (nspec_coadd, nspat_coadd):
              The averaged tilts image corresponding to the rectified
              and coadded data.
            - waveimg: float ndarray shape = (nspec_coadd, nspat_coadd):
              The averaged wavelength image corresponding to the
              rectified and coadded data.
            - dspat: float ndarray shape = (nspec_coadd, nspat_coadd):
              The average spatial offsets in pixels from the reference
              trace trace_stack corresponding to the rectified and
              coadded data.
            - thismask: bool ndarray shape = (nspec_coadd, nspat_coadd):
              Output mask for rectified and coadded images. True = Good,
              False=Bad. This image is trivial, and is simply an image
              of True values the same shape as the rectified and coadded
              data.
            - tslits_dict: dict: tslits_dict dictionary containing the
              information about the slits boundaries. The slit
              boundaries are trivial and are simply vertical traces at 0
              and nspat_coadd-1.
    """
    nimgs, nspec, nspat = sciimg_stack.shape

    if 'uniform' in weights:
        msgs.info('No weights were provided. Using uniform weights.')
        weights = np.ones(nimgs)/float(nimgs)

    weights_stack = combine.broadcast_weights(weights, sciimg_stack.shape)

    # Determine the wavelength grid that we will use for the current slit/order
    wave_bins = get_wave_bins(thismask_stack, waveimg_stack, wave_grid)
    dspat_bins, dspat_stack = get_spat_bins(thismask_stack, ref_trace_stack)

    sci_list = [weights_stack, sciimg_stack, sciimg_stack - skymodel_stack, tilts_stack,
                waveimg_stack, dspat_stack]
    var_list = [utils.calc_ivar(sciivar_stack)]

    sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack \
            = rebin2d(wave_bins, dspat_bins, waveimg_stack, dspat_stack, thismask_stack,
                      inmask_stack, sci_list, var_list)
    # Now compute the final stack with sigma clipping
    sigrej = 3.0
    maxiters = 10
    # sci_list_rebin[0] = rebinned weights image stack
    # sci_list_rebin[1:] = stacks of images that we want to weighted combine
    # sci_list_rebin[2] = rebinned sciimg-sky_model images that we used for the sigma clipping
    sci_list_out, var_list_out, outmask, nused \
            = combine.weighted_combine(sci_list_rebin[0], sci_list_rebin[1:], var_list_rebin,
                               norm_rebin_stack != 0, sigma_clip=True,
                               sigma_clip_stack=sci_list_rebin[2], sigrej=sigrej,
                               maxiters=maxiters)
    sciimg, imgminsky, tilts, waveimg, dspat = sci_list_out
    sciivar = utils.calc_ivar(var_list_out[0])

    # Compute the midpoints vectors, and lower/upper bins of the rectified image
    wave_mid = ((wave_bins + np.roll(wave_bins,1))/2.0)[1:]
    wave_min = wave_bins[:-1]
    wave_max = wave_bins[1:]
    dspat_mid = ((dspat_bins + np.roll(dspat_bins,1))/2.0)[1:]

    # Interpolate the dspat images wherever the coadds are masked
    # because a given pixel was not sampled. This is done because the
    # dspat image is not allowed to have holes if it is going to work
    # with local_skysub_extract
    nspec_coadd, nspat_coadd = imgminsky.shape
    spat_img_coadd, spec_img_coadd = np.meshgrid(np.arange(nspat_coadd), np.arange(nspec_coadd))

    if np.any(np.invert(outmask)):
        points_good = np.stack((spec_img_coadd[outmask], spat_img_coadd[outmask]), axis=1)
        points_bad = np.stack((spec_img_coadd[np.invert(outmask)],
                                spat_img_coadd[np.invert(outmask)]), axis=1)
        values_dspat = dspat[outmask]
        dspat_bad = scipy.interpolate.griddata(points_good, values_dspat, points_bad,
                                               method='cubic')
        dspat[np.invert(outmask)] = dspat_bad
        # Points outside the convex hull of the data are set to nan. We
        # identify those and simply assume them values from the
        # dspat_img_fake, which is what dspat would be on a regular
        # perfectly rectified image grid.
        nanpix = np.isnan(dspat)
        if np.any(nanpix):
            dspat_img_fake = spat_img_coadd + dspat_mid[0]
            dspat[nanpix] = dspat_img_fake[nanpix]

    return dict(wave_bins=wave_bins, dspat_bins=dspat_bins, wave_mid=wave_mid, wave_min=wave_min,
                wave_max=wave_max, dspat_mid=dspat_mid, sciimg=sciimg, sciivar=sciivar,
                imgminsky=imgminsky, outmask=outmask, nused=nused, tilts=tilts, waveimg=waveimg,
                dspat=dspat, nspec=imgminsky.shape[0], nspat=imgminsky.shape[1])



def rebin2d(spec_bins, spat_bins, waveimg_stack, spatimg_stack, thismask_stack, inmask_stack, sci_list, var_list):
    """
    Rebin a set of images and propagate variance onto a new spectral and spatial grid. This routine effectively
    "recitifies" images using np.histogram2d which is extremely fast and effectiveluy performs
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
        tuple: Returns the following:
            - sci_list_out: list: The list of ndarray rebinned images
              with new shape (nimgs, nspec_rebin, nspat_rebin)
            - var_list_out: list: The list of ndarray rebinned variance
              images with correct error propagation with shape (nimgs,
              nspec_rebin, nspat_rebin)
            - norm_rebin_stack: int ndarray, shape (nimgs, nspec_rebin,
              nspat_rebin): An image stack indicating the integer
              occupation number of a given pixel. In other words, this
              number would be zero for empty bins, one for bins that
              were populated by a single pixel, etc. This image takes
              the input inmask_stack into account. The output mask for
              each image can be formed via outmask_rebin_satck =
              (norm_rebin_stack > 0)
            - nsmp_rebin_stack: int ndarray, shape (nimgs, nspec_rebin,
              nspat_rebin): An image stack indicating the integer
              occupation number of a given pixel taking only the
              thismask_stack into account, but taking the inmask_stack
              into account. This image is mainly constructed for
              bookeeping purposes, as it represents the number of times
              each pixel in the rebin image was populated taking only
              the "geometry" of the rebinning into account (i.e. the
              thismask_stack), but not the masking (inmask_stack).
    """

    shape = combine.img_list_error_check(sci_list, var_list)
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
# TODO: Move this out of core

class CoAdd2d(object):

    """
    Main routine to run the extraction for 2d coadds.

    Algorithm steps are as follows:
        - Fill this in.

    This performs 2d coadd specific tasks, and then also performs some
    of the tasks analogous to the pypeit.extract_one method. Docs coming
    soon....

    Args:
        stack_dict:
        master_dir:
        det (int):
        samp_fact: float
           sampling factor to make the wavelength grid finer or coarser.  samp_fact > 1.0 oversamples (finer),
           samp_fact < 1.0 undersamples (coarser)
        ir_redux:
        par:
        show:
        show_peaks:

    """


    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                     ir_redux=False, master_dir=None, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """
        Instantiate the CoAdd2d subclass appropriate for the provided spectrograph.

        The class must be subclassed this class CoAdd2d.

        Args:
            spec2dfiles (list):
                List of spec2d files
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                The instrument used to collect the data to be reduced.
            par (:class:`pypeit.par.parset.ParSet`):
                Parset object

            master_dir (:obj:`str`, optional):
                Directory for the coadding master frames. If None,
                set by ``os.getcwd()``.

        Returns:
            :class:`CoAdd2d`: One of the subclasses with :class:`CoAdd2d` as its
            base.
        """

        return next(c for c in cls.__subclasses__() if c.__name__ == spectrograph.pypeline)(
            spec2dfiles, spectrograph, par, det=det, offsets=offsets, weights=weights, sn_smooth_npix=sn_smooth_npix,
            ir_redux=ir_redux, master_dir=master_dir, show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug, **kwargs_wave)

    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, master_dir=None, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """

        Args:
            spec2d_files:
            det:
            offsets (ndarray): default=None
                Spatial offsets to be applied to each image before coadding. For the default mode of None, images
                are registered automatically using the trace of the brightest object.
            weights (str, list or ndarray):
                Mode for the weights used to coadd images. Options are 'auto' (default), 'uniform', or list/array of
                weights with shape = (nexp,) can be input and will be applied to the image. Note 'auto' is not allowed
                if offsets are input, and if set this will cause an exception.
            sn_smooth_npix:
            ir_redux:
            master_dir (:obj:`str`, optional):
                Directory for the coadding master frames. If None,
                set by ``os.getcwd()``.
            par:
            std:
            show:
            show_peaks:
            debug:
            **kwargs_wave:
        """

        ## Use Cases:
        #  1) offsets is None -- auto compute offsets from brightest object, so then default to auto_weights=True
        #  2) offsets not None, weights = None (uniform weighting) or weights is not None (input weights)
        #  3) offsets not None, auto_weights=True (Do not support)
        if offsets is not None and 'auto' in weights:
            msgs.error("Automatic weights cannot be computed for input offsets. "
                       "Set weights='uniform' or input an array of weights with shape (nexp,)")
        self.spec2d_files = spec2d_files
        self.spectrograph = spectrograph
        self.par = par
        self.det = det
        self.offsets = offsets
        self.weights = weights
        self.ir_redux = ir_redux
        self.master_dir = os.getcwd() if master_dir is None else master_dir
        self.show = show
        self.show_peaks = show_peaks
        self.debug_offsets = debug_offsets
        self.debug = debug
        self.stack_dict = None
        self.pseudo_dict = None

        self.objid_bri = None
        self.slitid_bri  = None
        self.snr_bar_bri = None


        # Load the stack_dict
        self.stack_dict = self.load_coadd2d_stacks(self.spec2d_files)
        self.pypeline = self.spectrograph.pypeline

        # Check that there are the same number of slits on every exposure
        nslits_list = [slits.nslits for slits in self.stack_dict['slits_list']]
        if not len(set(nslits_list))==1:
            msgs.error('Not all of your exposures have the same number of slits. Check your inputs')
        # TODO: Do the same check above but for the shape and binning
        # of the input images?
        self.nslits = nslits_list[0]
        self.nexp = len(self.stack_dict['specobjs_list'])
        self.nspec = self.stack_dict['slits_list'][0].nspec
        self.binning = np.array([self.stack_dict['slits_list'][0].binspec,
                                 self.stack_dict['slits_list'][0].binspat])

        # If smoothing is not input, smooth by 10% of the spectral dimension
        self.sn_smooth_npix = sn_smooth_npix if sn_smooth_npix is not None else 0.1*self.nspec

    def optimal_weights(self, slitorderid, objid, const_weights=False):
        """
        Determine optimal weights for 2d coadds. This script grabs the information from SpecObjs list for the
        object with specified slitid and objid and passes to coadd.sn_weights to determine the optimal weights for
        each exposure.

        Args:
            slitorderid (int):
               The slit or order id that has the brightest object whose S/N will be used to determine the weight for each frame.
            objid (np.ndarray):
               Array of object indices with  shape = (nexp,) of the brightest object whose S/N will be used to determine the weight for each frame.
            const_weights (bool):
               Use constant weights for coadding the exposures. Default=False

        Returns:
            rms_sn : ndarray, shape = (len(specobjs_list),)
                Root mean square S/N value for each input spectra
            weights : ndarray, shape (len(specobjs_list),)
                Weights to be applied to the spectra. These are signal-to-noise squared weights.
        """

        nexp = len(self.stack_dict['specobjs_list'])
        nspec = self.stack_dict['specobjs_list'][0][0].TRACE_SPAT.shape[0]
        # Grab the traces, flux, wavelength and noise for this slit and objid.
        flux_stack = np.zeros((nspec, nexp), dtype=float)
        ivar_stack = np.zeros((nspec, nexp), dtype=float)
        wave_stack = np.zeros((nspec, nexp), dtype=float)
        mask_stack = np.zeros((nspec, nexp), dtype=bool)

        for iexp, sobjs in enumerate(self.stack_dict['specobjs_list']):
            ithis = sobjs.slitorder_objid_indices(slitorderid, objid[iexp])
            flux_stack[:, iexp] = sobjs[ithis].OPT_COUNTS
            ivar_stack[:, iexp] = sobjs[ithis].OPT_COUNTS_IVAR
            wave_stack[:, iexp] = sobjs[ithis].OPT_WAVE
            mask_stack[:, iexp] = sobjs[ithis].OPT_MASK

        # TODO For now just use the zero as the reference for the wavelengths? Perhaps we should be rebinning the data though?
        rms_sn, weights = coadd.sn_weights(wave_stack, flux_stack, ivar_stack, mask_stack, self.sn_smooth_npix,
                                           const_weights=const_weights)
        return rms_sn, weights.T


    def coadd(self, only_slits=None):

        only_slits = [only_slits] if (only_slits is not None and
                                      isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        good_slits = np.arange(self.nslits) if only_slits is None else only_slits

        coadd_list = []
        for islit in good_slits:
            msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(islit, self.nslits - 1))
            ref_trace_stack = self.reference_trace_stack(islit, offsets=self.offsets, objid=self.objid_bri)
            thismask_stack = self.stack_dict['slitmask_stack'] == islit
            # TODO Can we get rid of this one line simply making the weights returned by parse_weights an
            # (nslit, nexp) array?
            # This one line deals with the different weighting strategies between MultiSlit echelle. Otherwise, we
            # would need to copy this method twice in the subclasses
            if 'auto_echelle' in self.use_weights:
                rms_sn, weights = self.optimal_weights(islit, self.objid_bri)
            else:
                weights = self.use_weights
            # Perform the 2d coadd
            coadd_dict = compute_coadd2d(ref_trace_stack, self.stack_dict['sciimg_stack'],
                                           self.stack_dict['sciivar_stack'],
                                           self.stack_dict['skymodel_stack'], self.stack_dict['mask_stack'] == 0,
                                           self.stack_dict['tilts_stack'], thismask_stack,
                                           self.stack_dict['waveimg_stack'],
                                           self.wave_grid, weights=weights)
            coadd_list.append(coadd_dict)

        return coadd_list


    def create_pseudo_image(self, coadd_list):
        """ THIS UNDOCUMENTED CODE PROBABLY SHOULD GENERATE AND RETURN
        STANDARD PYPEIT OBJCTS INSTEAD OF SOME UNDEFINED DICT"""



        nspec_vec = np.zeros(self.nslits,dtype=int)
        nspat_vec = np.zeros(self.nslits,dtype=int)
        for islit, cdict in enumerate(coadd_list):
            nspec_vec[islit]=cdict['nspec']
            nspat_vec[islit]=cdict['nspat']

        # Determine the size of the pseudo image
        nspat_pad = 10
        nspec_pseudo = nspec_vec.max()
        nspat_pseudo = np.sum(nspat_vec) + (self.nslits + 1)*nspat_pad
        spec_vec_pseudo = np.arange(nspec_pseudo)
        shape_pseudo = (nspec_pseudo, nspat_pseudo)
        imgminsky_pseudo = np.zeros(shape_pseudo)
        sciivar_pseudo = np.zeros(shape_pseudo)
        waveimg_pseudo = np.zeros(shape_pseudo)
        tilts_pseudo = np.zeros(shape_pseudo)
        spat_img_pseudo = np.zeros(shape_pseudo)
        nused_pseudo = np.zeros(shape_pseudo, dtype=int)
        inmask_pseudo = np.zeros(shape_pseudo, dtype=bool)
        wave_mid = np.zeros((nspec_pseudo, self.nslits))
        wave_mask = np.zeros((nspec_pseudo, self.nslits),dtype=bool)
        wave_min = np.zeros((nspec_pseudo, self.nslits))
        wave_max = np.zeros((nspec_pseudo, self.nslits))
        dspat_mid = np.zeros((nspat_pseudo, self.nslits))

        spat_left = nspat_pad
        slit_left = np.zeros((nspec_pseudo, self.nslits))
        slit_righ = np.zeros((nspec_pseudo, self.nslits))
        spec_min1 = np.zeros(self.nslits)
        spec_max1 = np.zeros(self.nslits)

        nspec_grid = self.wave_grid_mid.size
        for islit, coadd_dict in enumerate(coadd_list):
            spat_righ = spat_left + nspat_vec[islit]
            ispec = slice(0,nspec_vec[islit])
            ispat = slice(spat_left,spat_righ)
            imgminsky_pseudo[ispec, ispat] = coadd_dict['imgminsky']
            sciivar_pseudo[ispec, ispat] = coadd_dict['sciivar']
            waveimg_pseudo[ispec, ispat] = coadd_dict['waveimg']
            tilts_pseudo[ispec, ispat] = coadd_dict['tilts']
            # spat_img_pseudo is the sub-pixel image position on the rebinned pseudo image
            inmask_pseudo[ispec, ispat] = coadd_dict['outmask']
            image_temp = (coadd_dict['dspat'] -  coadd_dict['dspat_mid'][0] + spat_left)*coadd_dict['outmask']
            spat_img_pseudo[ispec, ispat] = image_temp
            nused_pseudo[ispec, ispat] = coadd_dict['nused']
            wave_min[ispec, islit] = coadd_dict['wave_min']
            wave_max[ispec, islit] = coadd_dict['wave_max']
            wave_mid[ispec, islit] = coadd_dict['wave_mid']
            wave_mask[ispec, islit] = True
            # Fill in the rest of the wave_mid with the corresponding points in the wave_grid
            #wave_this = wave_mid[wave_mask[:,islit], islit]
            #ind_upper = np.argmin(np.abs(self.wave_grid_mid - wave_this.max())) + 1
            #if nspec_vec[islit] != nspec_pseudo:
            #    wave_mid[nspec_vec[islit]:, islit] = self.wave_grid_mid[ind_upper:ind_upper + (nspec_pseudo-nspec_vec[islit])]


            dspat_mid[ispat, islit] = coadd_dict['dspat_mid']
            slit_left[:,islit] = np.full(nspec_pseudo, spat_left)
            slit_righ[:,islit] = np.full(nspec_pseudo, spat_righ)
            spec_max1[islit] = nspec_vec[islit]-1
            spat_left = spat_righ + nspat_pad

        slits_pseudo \
                = slittrace.SlitTraceSet(slit_left, slit_righ, nspat=nspat_pseudo,
                                         spectrograph=self.spectrograph.spectrograph,
                                         specmin=spec_min1, specmax=spec_max1,
                                         master_key=self.stack_dict['master_key_dict']['trace'],
                                         master_dir=self.master_dir)
        slitmask_pseudo = slits_pseudo.slit_img()
        # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
        # which wreaks havoc on the local sky-subtraction
        min_slit_frac = 0.70
        spec_min = np.zeros(self.nslits)
        spec_max = np.zeros(self.nslits)
        for islit in range(self.nslits):
            slit_width = np.sum(inmask_pseudo*(slitmask_pseudo == islit),axis=1)
            slit_width_img = np.outer(slit_width, np.ones(nspat_pseudo))
            med_slit_width = np.median(slit_width_img[slitmask_pseudo == islit])
            nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
            nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
            slit_width_sm = scipy.ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
            igood = (slit_width_sm > min_slit_frac*med_slit_width)
            spec_min[islit] = spec_vec_pseudo[igood].min()
            spec_max[islit] = spec_vec_pseudo[igood].max()
            bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_pseudo == islit)
            inmask_pseudo[bad_pix] = False

        # Update with tslits_dict_pseudo
        slits_pseudo.specmin = spec_min
        slits_pseudo.specmax = spec_max

        return dict(nspec=nspec_pseudo, nspat=nspat_pseudo, imgminsky=imgminsky_pseudo, sciivar=sciivar_pseudo,
                           inmask=inmask_pseudo, tilts=tilts_pseudo,
                           waveimg=waveimg_pseudo, spat_img = spat_img_pseudo,
                           slits=slits_pseudo,
                           wave_mask=wave_mask, wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max)

    def reduce(self, pseudo_dict, show=None, show_peaks=None):

        show = self.show if show is None else show
        show_peaks = self.show_peaks if show_peaks is None else show_peaks

        # Generate a ScienceImage
        sciImage = scienceimage.ScienceImage(self.spectrograph, self.det,
                                                      self.par['scienceframe']['process'],
                                                      pseudo_dict['imgminsky'],
                                                      pseudo_dict['sciivar'],
                                                      np.zeros_like(pseudo_dict['inmask']),  # Dummy bpm
                                                      rn2img=np.zeros_like(pseudo_dict['inmask']),  # Dummy rn2img
                                                      crmask=np.invert(pseudo_dict['inmask']))
        slitmask_pseudo = pseudo_dict['slits'].slit_img()
        sciImage.build_mask(slitmask=slitmask_pseudo)

        # Make changes to parset specific to 2d coadds
        parcopy = copy.deepcopy(self.par)
        parcopy['scienceimage']['findobj']['trace_npoly'] = 3        # Low order traces since we are rectified
        #parcopy['scienceimage']['find_extrap_npoly'] = 1  # Use low order for trace extrapolation
        # Instantiate Calibrations class
        caliBrate = calibrations.MultiSlitCalibrations(None, parcopy['calibrations'], self.spectrograph)
        caliBrate.slits = pseudo_dict['slits']
        caliBrate.tilts_dict = dict(tilts=pseudo_dict['tilts'])
        caliBrate.mswave = pseudo_dict['waveimg']
        #
        redux=reduce.instantiate_me(sciImage, self.spectrograph, parcopy, caliBrate,
                                    ir_redux=self.ir_redux, objtype='science_coadd2d',
                                    det=self.det, binning=self.binning, show=show)

        if show:
            redux.show('image', image=pseudo_dict['imgminsky']*(sciImage.mask == 0), chname = 'imgminsky', slits=True, clear=True)
        # Object finding
        sobjs_obj, nobj, skymask_init = redux.find_objects(sciImage.image, show_peaks=show_peaks)
        # Local sky-subtraction
        global_sky_pseudo = np.zeros_like(pseudo_dict['imgminsky']) # No global sky for co-adds since we go straight to local
        skymodel_pseudo, objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs = redux.local_skysub_extract(
            caliBrate.mswave, global_sky_pseudo, sobjs_obj, spat_pix=pseudo_dict['spat_img'], model_noise=False,
            show_profile=show, show=show)

        if self.ir_redux:
            sobjs.purge_neg()

        # TODO: Removed this, but I'm not sure that's what you want...
#        # Add the information about the fixed wavelength grid to the sobjs
#        for spec in sobjs:
#            idx = spec.slit_orderindx
#            # Fill
#            spec.BOX_WAVE_GRID_MASK, spec.OPT_WAVE_GRID_MASK = [pseudo_dict['wave_mask'][:,idx]]*2
#            spec.BOX_WAVE_GRID, spec.OPT_WAVE_GRID = [pseudo_dict['wave_mid'][:,idx]]*2
#            spec.BOX_WAVE_GRID_MIN, spec.OPT_WAVE_GRID_MIN = [pseudo_dict['wave_min'][:,idx]]*2
#            spec.BOX_WAVE_GRID_MAX, spec.OPT_WAVE_GRID_MAX = [pseudo_dict['wave_max'][:,idx]]*2

        # Add the rest to the pseudo_dict
        pseudo_dict['skymodel'] = skymodel_pseudo
        pseudo_dict['objmodel'] = objmodel_pseudo
        pseudo_dict['ivarmodel'] = ivarmodel_pseudo
        pseudo_dict['outmask'] = outmask_pseudo
        pseudo_dict['sobjs'] = sobjs
        self.pseudo_dict=pseudo_dict

        return pseudo_dict['imgminsky'], pseudo_dict['sciivar'], skymodel_pseudo, objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs


    def save_masters(self):

        # Write out the pseudo master files to disk
        master_key_dict = self.stack_dict['master_key_dict']

        # TODO: These saving operations are a temporary kludge
        waveImage = WaveImage(None, None, None, self.spectrograph,  # spectrograph is needed for header
                              None, None, master_key=master_key_dict['arc'],
                              master_dir=self.master_dir)
        waveImage.save(image=self.pseudo_dict['waveimg'])

        # TODO: Assumes overwrite=True
        self.pseudo_dict['slits'].to_master()

    def snr_report(self, snr_bar, slitid=None):

        # Print out a report on the SNR
        msg_string = msgs.newline() + '-------------------------------------'
        msg_string += msgs.newline() + '  Summary for highest S/N object'
        if slitid is not None:
            msg_string += msgs.newline() + '      found on slitid = {:d}            '.format(slitid)
        msg_string += msgs.newline() + '-------------------------------------'
        msg_string += msgs.newline() + '           exp#        S/N'
        for iexp, snr in enumerate(snr_bar):
            msg_string += msgs.newline() + '            {:d}         {:5.2f}'.format(iexp, snr)

        msg_string += msgs.newline() + '-------------------------------------'
        msgs.info(msg_string)

    def get_good_slits(self, only_slits):

        only_slits = [only_slits] if (only_slits is not None and
                                        isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        good_slits = np.arange(self.nslits) if only_slits is None else only_slits
        return good_slits

    def offset_slit_cen(self, slitid, offsets):
        # TODO: Check that slitid is available for all slit objects
        # TODO: Check that all slits have the same nspec
        ref_trace_stack = np.zeros((self.stack_dict['slits_list'][0].nspec, len(offsets)),
                                   dtype=float)
        for iexp, slits in enumerate(self.stack_dict['slits_list']):
            ref_trace_stack[:, iexp] = slits.center[:,slitid] - offsets[iexp]
        return ref_trace_stack

    def get_wave_grid(self, **kwargs_wave):
        """
        Routine to create a wavelength grid for 2d coadds using all of the wavelengths of the extracted objects. Calls
        coadd1d.get_wave_grid.

        Args:
            **kwargs_wave (dict):
                Optional argumments for coadd1d.get_wve_grid function

        Returns:
            tuple: Returns the following:
                - wave_grid (np.ndarray): New wavelength grid, not
                  masked
                - wave_grid_mid (np.ndarray): New wavelength grid
                  evaluated at the centers of the wavelength bins, that
                  is this grid is simply offset from wave_grid by
                  dsamp/2.0, in either linear space or log10 depending
                  on whether linear or (log10 or velocity) was
                  requested.  For iref or concatenate the linear
                  wavelength sampling will be calculated.
                - dsamp (float): The pixel sampling for wavelength grid
                  created.
        """

        nobjs_tot = np.array([len(spec) for spec in self.stack_dict['specobjs_list']]).sum()
        waves = np.zeros((self.nspec, nobjs_tot))
        masks = np.zeros_like(waves, dtype=bool)
        indx = 0
        for spec_this in self.stack_dict['specobjs_list']:
            for spec in spec_this:
                waves[:, indx] = spec.OPT_WAVE
                masks[:, indx] = spec.OPT_MASK
                indx += 1

        wave_grid, wave_grid_mid, dsamp = coadd.get_wave_grid(waves, masks=masks, **kwargs_wave)

        return wave_grid, wave_grid_mid, dsamp

    def load_coadd2d_stacks(self, spec2d_files):
        """
        Routine to read in required images for 2d coadds given a list of spec2d files.

        Args:
            spec2d_files: list
               List of spec2d filenames
            det: int
               detector in question

        Returns:
            dict: Dictionary containing all the images and keys required
            for perfomring 2d coadds.
        """

        # Get the detector string
        sdet = parse.get_dnum(self.det, prefix=False)

        # Get the master dir

        redux_path = os.getcwd()

        # Grab the files
        head2d_list = []
        tracefiles = []
        waveimgfiles = []
        tiltfiles = []
        spec1d_files = []
        for f in spec2d_files:
            head = fits.getheader(f)
            if os.path.exists(head['PYPMFDIR']):
                master_path = head['PYPMFDIR']
            else:
                master_dir = os.path.basename(head['PYPMFDIR'])
                master_path = os.path.join(os.path.split(os.path.split(f)[0])[0], master_dir)

            trace_key = '{0}_{1:02d}'.format(head['TRACMKEY'], self.det)
            wave_key = '{0}_{1:02d}'.format(head['ARCMKEY'], self.det)

            head2d_list.append(head)
            spec1d_files.append(f.replace('spec2d', 'spec1d'))
            tracefiles.append(os.path.join(master_path,
                            '{0}.gz'.format(MasterFrame.construct_file_name('Slits', trace_key))))
#                                           MasterFrame.construct_file_name('Trace', trace_key)))
            waveimgfiles.append(os.path.join(master_path,
                                             MasterFrame.construct_file_name('Wave', wave_key)))
            tiltfiles.append(os.path.join(master_path,
                                          MasterFrame.construct_file_name('Tilts', wave_key)))

        nfiles = len(spec2d_files)

        specobjs_list = []
        head1d_list = []
        slits_list = []
        # TODO Sort this out with the correct detector extensions etc.
        # Read in the image stacks
        waveimgfile, tiltfile, tracefile = None, None, None
        for ifile in range(nfiles):
            # Load up the calibs, if needed
            if waveimgfiles[ifile] != waveimgfile:
                waveimg = WaveImage.from_master_file(waveimgfiles[ifile]).image
            if tiltfile != tiltfiles[ifile]:
                tilts = WaveTilts.from_master_file(tiltfiles[ifile]).tilts_dict
            # Save
            waveimgfile = waveimgfiles[ifile]
            tiltfile = tiltfiles[ifile]
            #
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
                shape_wave = (nfiles, waveimg.shape[0], waveimg.shape[1])
                shape_sci = (nfiles, sciimg.shape[0], sciimg.shape[1])
                waveimg_stack = np.zeros(shape_wave, dtype=float)
                tilts_stack = np.zeros(shape_wave, dtype=float)
                sciimg_stack = np.zeros(shape_sci, dtype=float)
                skymodel_stack = np.zeros(shape_sci, dtype=float)
                sciivar_stack = np.zeros(shape_sci, dtype=float)
                mask_stack = np.zeros(shape_sci, dtype=float)
                slitmask_stack = np.zeros(shape_sci, dtype=float)

            # Slit Traces and slitmask
            # TODO: Don't understand this if statement
            if tracefile != tracefiles[ifile]:
                slits = slittrace.SlitTraceSet.from_file(tracefiles[ifile])
                # Check the spectrograph names
                # TODO: Should this be done here?
                if slits.spectrograph != self.spectrograph.spectrograph:
                    msgs.error('Spectrograph read from {0} is not correct.  Expected {1}.'.format(
                                tracefiles[ifile], self.spectrograph.spectrograph))

            tracefile = tracefiles[ifile]
            #
            slits_list.append(slits)
            slitmask_stack[ifile, :, :] = slits.slit_img()
            waveimg_stack[ifile, :, :] = waveimg
            tilts_stack[ifile, :, :] = tilts['tilts']
            sciimg_stack[ifile, :, :] = sciimg
            sciivar_stack[ifile, :, :] = sciivar
            mask_stack[ifile, :, :] = mask
            skymodel_stack[ifile, :, :] = skymodel

            # Specobjs
            if os.path.isfile(spec1d_files[ifile]):
                sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_files[ifile])
                head1d_list.append(sobjs.header)
                this_det = sobjs.DET == self.det
                specobjs_list.append(sobjs[this_det])

        # slitmask_stack = np.einsum('i,jk->ijk', np.ones(nfiles), slitmask)

        # Fill the master key dict
        head2d = head2d_list[0]
        master_key_dict = {}
        master_key_dict['frame'] = head2d['FRAMMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['bpm'] = head2d['BPMMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['bias'] = head2d['BIASMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['arc'] = head2d['ARCMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['trace'] = head2d['TRACMKEY'] + '_{:02d}'.format(self.det)
        master_key_dict['flat'] = head2d['FLATMKEY'] + '_{:02d}'.format(self.det)

        # TODO In the future get this stuff from the headers once data model finalized
        # TODO: spectrograph already exists in self and is a required
        # argument of the init. So I use it here, and force it to be
        # the same as what's read by the SlitTraceSet file above.
        #spectrograph = util.load_spectrograph(tslits_dict['spectrograph'])

        return dict(specobjs_list=specobjs_list, slits_list=slits_list,
                    slitmask_stack=slitmask_stack,
                    sciimg_stack=sciimg_stack, sciivar_stack=sciivar_stack,
                    skymodel_stack=skymodel_stack, mask_stack=mask_stack,
                    tilts_stack=tilts_stack, waveimg_stack=waveimg_stack,
                    head1d_list=head1d_list, head2d_list=head2d_list,
                    redux_path=redux_path,
                    master_key_dict=master_key_dict,
                    spectrograph=self.spectrograph.spectrograph,
                    pypeline=self.spectrograph.pypeline)

# Multislit can coadd with:
# 1) input offsets or if offsets is None, it will find the brightest trace and compute them
# 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object

# Echelle can either stack with:
# 1) input offsets or if offsets is None, it will find the objid of brightest trace and stack all orders relative to the trace of this object.
# 2) specified weights, or if weights is None and auto_weights=True,
#    it will use wavelength dependent weights determined from the spectrum of the brightest objects objid on each order

class MultiSlit(CoAdd2d):
    """
    Child of Coadd2d for Multislit and Longslit reductions

        # Multislit can coadd with:
        # 1) input offsets or if offsets is None, it will find the brightest trace and compute them
        # 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object


    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, master_dir=None, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(MultiSlit, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                        sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, master_dir=master_dir, par=par,
                                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                                        debug=debug, **kwargs_wave)


        ## Use Cases:
        #  1) offsets is None -- auto compute offsets from brightest object, so then default to auto_weights=True
        #  2) offsets not None, weights = None (uniform weighting) or weights is not None (input weights)
        #  3) offsets not None, auto_weights=True (Do not support)

        # Default wave_method for Multislit is linear
        kwargs_wave['wave_method'] = 'linear' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        if offsets is None:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri, self.offsets = self.compute_offsets()

        self.use_weights = self.parse_weights(weights)

    def parse_weights(self, weights):

        if 'auto' in weights:
            rms_sn, use_weights = self.optimal_weights(self.slitid_bri, self.objid_bri, const_weights=True)
            return use_weights
        elif 'uniform' in weights:
            return 'uniform'
        elif isinstance(weights, (list, np.ndarray)):
            if len(weights) != self.nexp:
                msgs.error('If weights are input it must be a list/array with same number of elements as exposures')
            return weights
        else:
            msgs.error('Unrecognized format for weights')

    # TODO When we run multislit, we actually compute the rebinned images twice. Once here to compute the offsets
    # and another time to weighted_combine the images in compute2d. This could be sped up
    def compute_offsets(self):

        objid_bri, slitid_bri, snr_bar_bri = self.get_brightest_obj(self.stack_dict['specobjs_list'], self.nslits)
        msgs.info('Determining offsets using brightest object on slit: {:d} with avg SNR={:5.2f}'.format(slitid_bri,np.mean(snr_bar_bri)))
        thismask_stack = self.stack_dict['slitmask_stack'] == slitid_bri
        trace_stack_bri = np.zeros((self.nspec, self.nexp))
        # TODO Need to think abbout whether we have multiple tslits_dict for each exposure or a single one
        for iexp in range(self.nexp):
            trace_stack_bri[:,iexp] = self.stack_dict['slits_list'][iexp].center[:,slitid_bri]
#            trace_stack_bri[:,iexp] = (self.stack_dict['tslits_dict_list'][iexp]['slit_left'][:,slitid_bri] +
#                                       self.stack_dict['tslits_dict_list'][iexp]['slit_righ'][:,slitid_bri])/2.0
        # Determine the wavelength grid that we will use for the current slit/order
        wave_bins = get_wave_bins(thismask_stack, self.stack_dict['waveimg_stack'], self.wave_grid)
        dspat_bins, dspat_stack = get_spat_bins(thismask_stack, trace_stack_bri)

        sci_list = [self.stack_dict['sciimg_stack'] - self.stack_dict['skymodel_stack']]
        var_list = []

        msgs.info('Rebinning Images')
        sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = rebin2d(
            wave_bins, dspat_bins, self.stack_dict['waveimg_stack'], dspat_stack, thismask_stack,
            (self.stack_dict['mask_stack'] == 0), sci_list, var_list)
        thismask = np.ones_like(sci_list_rebin[0][0,:,:],dtype=bool)
        nspec_pseudo, nspat_pseudo = thismask.shape
        slit_left = np.full(nspec_pseudo, 0.0)
        slit_righ = np.full(nspec_pseudo, nspat_pseudo)
        inmask = norm_rebin_stack > 0
        traces_rect = np.zeros((nspec_pseudo, self.nexp))
        sobjs = specobjs.SpecObjs()
        #specobj_dict = {'setup': 'unknown', 'slitid': 999, 'orderindx': 999, 'det': self.det, 'objtype': 'unknown',
        #                'pypeline': 'MultiSLit' + '_coadd_2d'}
        for iexp in range(self.nexp):
            sobjs_exp, _ = extract.objfind(sci_list_rebin[0][iexp,:,:], thismask, slit_left, slit_righ,
                                           inmask=inmask[iexp,:,:], ir_redux=self.ir_redux,
                                           fwhm=self.par['scienceimage']['findobj']['find_fwhm'],
                                           trim_edg=self.par['scienceimage']['findobj']['find_trim_edge'],
                                           npoly_cont=self.par['scienceimage']['findobj']['find_npoly_cont'],
                                           maxdev=self.par['scienceimage']['findobj']['find_maxdev'],
                                           ncoeff=3, sig_thresh=self.par['scienceimage']['findobj']['sig_thresh'], nperslit=1,
                                           show_trace=self.debug_offsets, show_peaks=self.debug_offsets)
            sobjs.add_sobj(sobjs_exp)
            traces_rect[:, iexp] = sobjs_exp.TRACE_SPAT
        # Now deterimine the offsets. Arbitrarily set the zeroth trace to the reference
        med_traces_rect = np.median(traces_rect,axis=0)
        offsets = med_traces_rect[0] - med_traces_rect
        # Print out a report on the offsets
        msg_string = msgs.newline()  + '---------------------------------------------'
        msg_string += msgs.newline() + ' Summary of offsets for highest S/N object   '
        msg_string += msgs.newline() + '         found on slitid = {:d}              '.format(slitid_bri)
        msg_string += msgs.newline() + '---------------------------------------------'
        msg_string += msgs.newline() + '           exp#      offset                  '
        for iexp, off in enumerate(offsets):
            msg_string += msgs.newline() + '            {:d}        {:5.2f}'.format(iexp, off)

        msg_string += msgs.newline() + '-----------------------------------------------'
        msgs.info(msg_string)
        if self.debug_offsets:
            for iexp in range(self.nexp):
                plt.plot(traces_rect[:, iexp], linestyle='--', label='original trace')
                plt.plot(traces_rect[:, iexp] + offsets[iexp], label='shifted traces')
                plt.legend()
            plt.show()

        return objid_bri, slitid_bri, snr_bar_bri, offsets

    def get_brightest_obj(self, specobjs_list, nslits):

        """
        Utility routine to find the brightest object in each exposure given a specobjs_list for MultiSlit reductions.

        Args:
            specobjs_list: list
               List of SpecObjs objects.
            echelle: bool, default=True, optional

        Returns:
            tuple: Returns the following:
                - objid: ndarray, int, shape (len(specobjs_list),):
                  Array of object ids representing the brightest object
                  in each exposure
                - slitid (int): Slit that highest S/N ratio object is on
                  (only for pypeline=MultiSlit)
                - snr_bar: ndarray, float, shape (len(list),): Average
                  S/N over all the orders for this object
        """
        nexp = len(specobjs_list)
        nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]

        slit_snr_max = np.full((nslits, nexp), -np.inf)
        objid_max = np.zeros((nslits, nexp), dtype=int)
        # Loop over each exposure, slit, find the brighest object on that slit for every exposure
        for iexp, sobjs in enumerate(specobjs_list):
            for islit in range(nslits):
                ithis = sobjs.SLITID == islit
                nobj_slit = np.sum(ithis)
                if np.any(ithis):
                    objid_this = sobjs[ithis].OBJID
                    flux = np.zeros((nspec, nobj_slit))
                    ivar = np.zeros((nspec, nobj_slit))
                    wave = np.zeros((nspec, nobj_slit))
                    mask = np.zeros((nspec, nobj_slit), dtype=bool)
                    for iobj, spec in enumerate(sobjs[ithis]):
                        flux[:, iobj] = spec.OPT_COUNTS
                        ivar[:, iobj] = spec.OPT_COUNTS_IVAR
                        wave[:, iobj] = spec.OPT_WAVE
                        mask[:, iobj] = spec.OPT_MASK
                    rms_sn, weights = coadd.sn_weights(wave, flux, ivar, mask, None, const_weights=True)
                    imax = np.argmax(rms_sn)
                    slit_snr_max[islit, iexp] = rms_sn[imax]
                    objid_max[islit, iexp] = objid_this[imax]
        # Find the highest snr object among all the slits
        slit_snr = np.mean(slit_snr_max, axis=1)
        slitid = slit_snr.argmax()
        snr_bar_mean = slit_snr[slitid]
        snr_bar = slit_snr_max[slitid, :]
        objid = objid_max[slitid, :]
        if (snr_bar_mean == -np.inf):
            msgs.error('You do not appear to have a unique reference object that was traced as the highest S/N '
                       'ratio on the same slit of every exposure')

        self.snr_report(snr_bar, slitid=slitid)

        return objid, slitid, snr_bar

    # TODO add an option here to actually use the reference trace for cases where they are on the same slit and it is
    # single slit???
    def reference_trace_stack(self, slitid, offsets=None, objid=None):

        return self.offset_slit_cen(slitid, offsets)


class Echelle(CoAdd2d):
    """
    Child of Coadd2d for Multislit and Longslit reductions

        # Echelle can either stack with:
        # 1) input offsets or if offsets is None, it will find the objid of brightest trace and stack all orders relative to the trace of this object.
        # 2) specified weights, or if weights is None and auto_weights=True,
        #    it will use wavelength dependent weights determined from the spectrum of the brightest objects objid on each order


    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, master_dir=None, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(Echelle, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                      sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, master_dir=master_dir, par=par,
                                      show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug,
                                      **kwargs_wave)

        # Default wave_method for Echelle is log10
        kwargs_wave['wave_method'] = 'log10' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        self.objid_bri = None
        self.slitid_bri  = None
        self.snr_bar_bri = None
        if offsets is None:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri = self.get_brightest_obj(self.stack_dict['specobjs_list'], self.nslits)

        self.use_weights = self.parse_weights(weights)

    def parse_weights(self, weights):

        if 'auto' in weights:
            return 'auto_echelle'
        elif 'uniform' in weights:
            return 'uniform'
        elif isinstance(weights, (list, np.ndarray)):
            if len(weights) != self.nexp:
                msgs.error('If weights are input it must be a list/array with same number of elements as exposures')
            return weights
        else:
            msgs.error('Unrecognized format for weights')

    def get_brightest_obj(self, specobjs_list, nslits):
        """
        Utility routine to find the brightest object in each exposure given a specobjs_list for Echelle reductions.

        Args:
            specobjs_list: list
               List of SpecObjs objects.
            echelle: bool, default=True, optional

        Returns:
            tuple: Returns the following:
                - objid: ndarray, int, shape (len(specobjs_list),):
                  Array of object ids representing the brightest object
                  in each exposure
                - snr_bar: ndarray, float, shape (len(list),): Average
                  S/N over all the orders for this object
        """
        nexp = len(specobjs_list)

        objid = np.zeros(nexp, dtype=int)
        snr_bar = np.zeros(nexp)
        # norders = specobjs_list[0].ech_orderindx.max() + 1
        for iexp, sobjs in enumerate(specobjs_list):
            uni_objid = np.unique(sobjs.ECH_OBJID)
            nobjs = len(uni_objid)
            order_snr = np.zeros((nslits, nobjs))
            for iord in range(nslits):
                for iobj in range(nobjs):
                    ind = (sobjs.ECH_ORDERINDX == iord) & (sobjs.ECH_OBJID == uni_objid[iobj])
                    flux = sobjs[ind][0].OPT_COUNTS
                    ivar = sobjs[ind][0].OPT_COUNTS_IVAR
                    wave = sobjs[ind][0].OPT_WAVE
                    mask = sobjs[ind][0].OPT_MASK
                    rms_sn, weights = coadd.sn_weights(wave, flux, ivar, mask, self.sn_smooth_npix, const_weights=True)
                    order_snr[iord, iobj] = rms_sn

            # Compute the average SNR and find the brightest object
            snr_bar_vec = np.mean(order_snr, axis=0)
            objid[iexp] = uni_objid[snr_bar_vec.argmax()]
            snr_bar[iexp] = snr_bar_vec[snr_bar_vec.argmax()]

        self.snr_report(snr_bar)

        return objid, None, snr_bar

    def reference_trace_stack(self, slitid, offsets=None, objid=None):
        """
        Utility function for determining the reference trace about which 2d coadds are performed.
        There are two modes of operation to determine the reference trace for the 2d coadd of a given slit/order:

         1) offsets: we stack about the center of the slit for the slit in question with the input offsets added
         2) ojbid: we stack about the trace ofa reference object for this slit given for each exposure by the input objid

        Either offsets or objid must be provided, but the code will raise an exception if both are provided.

        Args:
            slitid (int):
               The slit or order that we are currently considering
            stack_dict (dict):
               Dictionary containing all the images and keys required for perfomring 2d coadds.
            offsets (list or np.ndarray):
               An array of offsets with the same dimensionality as the nexp, the numer of images being coadded.
            objid: (list or np.ndarray):
               An array of objids with the same dimensionality as the nexp, the number of images being coadded.

        Returns:
            ref_trace_stack

            ref_trace_stack (np.ndarray):
                An array with shape (nspec, nexp) containing the reference trace for each of the nexp exposures.

        """

        if offsets is not None and objid is not None:
            msgs.errror('You can only input offsets or an objid, but not both')
        nexp = len(offsets) if offsets is not None else len(objid)
        if offsets is not None:
            return self.offset_slit_cen(slitid, offsets)
        elif objid is not None:
            specobjs_list = self.stack_dict['specobjs_list']
            nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]
            # Grab the traces, flux, wavelength and noise for this slit and objid.
            ref_trace_stack = np.zeros((nspec, nexp), dtype=float)
            for iexp, sobjs in enumerate(specobjs_list):
                ithis = (sobjs.ECH_ORDERINDX == slitid) & (sobjs.ECH_OBJID == objid[iexp])
                ref_trace_stack[:, iexp] = sobjs[ithis].TRACE_SPAT
            return ref_trace_stack
        else:
            msgs.error('You must input either offsets or an objid to determine the stack of reference traces')
            return None

