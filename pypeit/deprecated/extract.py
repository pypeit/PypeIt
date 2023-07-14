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
from pypeit import specobjs
#from pypeit import tracepca
from pypeit.core.pydl import spheregroup

from IPython import embed


def objs_in_slit(image, thismask, slit_left, slit_righ, inmask=None, fwhm=3.0, use_user_fwhm=False, boxcar_rad=7.,
                 maxdev=2.0, has_negative=False, spec_min_max=None, hand_extract_dict=None, std_trace=None,
                 ncoeff=5, nperslit=None, sig_thresh=10.0, peak_thresh=0.0, abs_thresh=0.0, trim_edg=(5, 5),
                 cont_sig_thresh=2.0, extract_maskwidth=4.0, specobj_dict=None, cont_fit=True, npoly_cont=1,
                 find_min_max=None, show_peaks=False, show_fits=False, show_trace=False, show_cont=False,
                 debug_all=False, qa_title='objfind', objfindQA_filename=None):
    """
    Find the location of objects in a slitmask slit or a echelle order.
    Args:
        image (`numpy.ndarray`_):
            Image to search for objects from. This image has shape
            (nspec, nspat) image.shape where the first dimension (nspec)
            is spectral, and second dimension (nspat) is spatial. Note
            this image can either have the sky background in it, or have
            already been sky subtracted.  Object finding works best on
            sky-subtracted images, but often one runs on the frame with
            sky first to identify the brightest objects which are then
            masked (see skymask below) in sky subtraction.
        thismask (`numpy.ndarray`_):
            Boolean mask image specifying the pixels which lie on the
            slit/order to search for objects on.  The convention is:
            True = on the slit/order, False = off the slit/order
        slit_left (`numpy.ndarray`_):
            Left boundary of slit/order to be extracted (given as
            floating pt pixels). This a 1-d array with shape (nspec, 1)
            or (nspec)
        slit_righ (`numpy.ndarray`_):
            Right boundary of slit/order to be extracted (given as
            floating pt pixels). This a 1-d array with shape (nspec, 1)
            or (nspec)
        det (:obj:`int`):
            Dectector number of slit to be extracted.
        inmask (`numpy.ndarray`_):
            Floating-point Input mask image.
        spec_min_max (:obj:`tuple`):
            This is tuple (float or int) of two elements which defines the minimum and
            maximum of the SLIT in the spectral direction on the
            detector. If not passed in it will be determined
            automatically from the thismask
        find_min_max (:obj:`tuple`):
            Tuple of integers that defines the minimum and maximum of your OBJECT
            in the spectral direction on the detector. It is only used for object finding.
            This parameter is helpful if your object only has emission lines or at high redshift
            and the trace only shows in part of the detector.
        fwhm (:obj:`float`):
            Estimated fwhm of the objects in pixels
        use_user_fwhm (:obj:`bool`):
            If True PypeIt will use the spatial profile fwm input by the user (i.e. the fwhm parameter above)
            rather than determine the spatial fwhm from the smashed spatial profile via the automated algorithm.
            Default = False.
        boxcar_rad (:obj:`float`, :obj:`int`):
            Boxcar radius in *pixels* to assign to each detected object and to be used later for boxcar extraction.
        maxdev (:obj:`float`):
            Maximum deviation of pixels from polynomial fit to trace
            used to reject bad pixels in trace fitting.
        hand_extract_dict(:obj:`dict`):
            Dictionary containing information about apertures requested
            by user that should be place by hand in the object list.
            This option is useful for cases like an emission line obect
            that the code fails to find with its significance threshold
        std_trace (`numpy.ndarray`_):
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        ncoeff (:obj:`int`):
            Order of legendre polynomial fits to the trace
        nperslit (:obj:`int`):
            Maximum number of objects allowed per slit. The code will
            take the nperslit most significant detections.
        sig_thresh (:obj:`float`):
            Significance threshold for object detection. The code uses
            the maximum of the thresholds defined by sig_thresh,
            peak_thresh, and abs_thresh.  For the default behavior
            peak_thresh and abs_thresh are zero, so sig_thresh defines
            the threshold.
        peak_thresh (:obj:`float`):
            Peak threshold for object detection. This is a number
            between 0 and 1 and represents the fraction of the brightest
            object on the slit that will be kept as an object, i.e. if
            ymax is the brightest object of the spectrum smashed out in
            the spectral direction, all objects with ypeak >
            peak_thresh*ymak are kept. The code uses the maximum of the
            thresholds defined by sig_thresh, peak_thers, and
            abs_thresh.
        abs_thresh (:obj:`float`):
            Absolute threshold for object detection.  Objects are found
            by smashing out the spectral direction along the curved
            slit/order traces, and abs_thresh is in the units of this
            smashed profile.  The code uses the maximum of the
            thresholds defined by sig_thresh, peak_thers, and
            abs_thresh.
        extract_maskwidth (:obj:`float`,optional):
            This parameter determines the initial size of the region in
            units of fwhm that will be used for local sky subtraction in
            the routine skysub.local_skysub_extract.
        cont_sig_thresh (:obj:`float`, optional):
            Significance threshold for peak detection for determinining which pixels to use for the iteratively
            fit continuum of the spectral direction smashed image. This is passed as the sigthresh parameter
            to core.arc.iter_continum. For extremely narrow slits that are almost filled by the object trace set
            this to a smaller number like 1.0 or disable continuum fitting altogether with cont_fit=False below.
            Default = 2.0.
        trim_edg (:obj:`tuple`):
            Ignore objects within this many pixels of the left and right
            slit boundaries, where the first element refers to the left
            and second refers to the right. This is tuple of 2 integers of floats
        has_negative (:obj:`bool`, optional):
            Image has negative object traces, i.e. for IR difference imaging. This impacts how the
            iterative conntinuum is fit to the spectral direction smashed image for object finding. Default=False
        cont_fit (:obj:`bool`):
            Fit a continuum to the illumination pattern across the slit when peak finding
        npoly_cont (:obj:`int`):
            Order of polynomial fit to the illumination pattern across the slit when peak finding
        specobj_dict (:obj:`dict`):
            Dictionary containing meta-data for the objects that will be
            propgated into the SpecObj objects, i.e. SLITID,
            detector, object type, and pipeline. The default is None, in
            which case the following dictionary will be used::

                specobj_dict = {'SLITID': 999, 'DET': 'DET01',
                                'OBJTYPE': 'unknown', 'PYPELINE': 'unknown'}
        show_peaks (:obj:`bool`):
            Whether plotting the QA of peak finding of your object in each order
        show_fits (:obj:`bool`):
            Plot trace fitting for final fits using PCA as crutch
        show_trace (:obj:`bool`):
            Whether display the resulting traces on top of the image
        show_cont (:obj:`bool`):
            Show debugging plot of the routine used to determine the spectrum continuum
        debug_all (:obj:`bool`):
            Show all the debugging plots?
        qa_title (:obj:`str`, optional):
            Title to be printed in the QA plots
        objfindQA_filename: (:obj:`str`, optional):
            Directory + filename of the object profile QA
    Returns:
        :class:`pypeit.specobjs.SpecObjs`: class containing the
              information about the objects found on the slit/order
    Note:
        Revision History:
            - 10-Mar-2005 -- First version written by D. Schlegel, LBL
            - 2005-2018 -- Improved by J. F. Hennawi and J. X. Prochaska
            - 23-June-2018 -- Ported to python by J. F. Hennawi and
              significantly improved
            - 01-Feb-2022 -- Skymask stripped out by JXP
    """

    # debug_all=True
    if debug_all:
        show_peaks = True
        show_fits = True
        show_trace = True
        show_cont = True

    if specobj_dict is None:
        specobj_dict = dict(SLITID=999, DET='DET01', OBJTYPE='unknown', PYPELINE='MultiSlit')

    # Check that peak_thresh values make sense
    if peak_thresh < 0 or peak_thresh > 1:
        msgs.error('Invalid value of peak_thresh. It must be between 0.0 and 1.0')

    nspec, nspat = image.shape
    specmid = nspec // 2

    # Some information about this slit we need for later when we instantiate specobj objects
    spec_vec = np.arange(nspec)
    spat_vec = np.arange(nspat)

    ximg, edgmask = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)

    # If a mask was not passed in, create it
    if inmask is None:
        inmask = thismask

    # If spec_min_max was not passed in, determine it from the thismask
    if spec_min_max is None or np.any([s is None for s in spec_min_max]):
        if spec_min_max is None:
            spec_min_max = [None, None]
        ispec, ispat = np.where(thismask)
        if spec_min_max[0] is None:
            spec_min_max[0] = ispec.min()
        if spec_min_max[1] is None:
            spec_min_max[1] = ispec.max()

    totmask = thismask & inmask & np.invert(edgmask)
    thisimg = image * totmask
    #  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    # nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())
    # Mask skypixels with 2 fwhm of edge
    left_asym = slit_left[:, None] + np.outer(xsize / nsamp, np.arange(nsamp))
    righ_asym = left_asym + np.outer(xsize / nsamp, np.ones(int(nsamp)))
    # This extract_asymbox2 call smashes the image in the spectral direction along the curved object traces
    # TODO Should we be passing the mask here with extract_asymbox or not?
    flux_spec = moment1d(thisimg, (left_asym + righ_asym) / 2, (righ_asym - left_asym),
                         fwgt=totmask.astype(float))[0]
    mask_spec = moment1d(totmask, (left_asym + righ_asym) / 2, (righ_asym - left_asym),
                         fwgt=totmask.astype(float))[0] < 0.3
    if find_min_max is not None:
        find_spec_min, find_spec_max = int(find_min_max[0]), int(find_min_max[1])
        flux_spec = flux_spec[find_spec_min:find_spec_max, :]
        mask_spec = mask_spec[find_spec_min:find_spec_max, :]

    flux_mean, flux_median, flux_sig \
        = stats.sigma_clipped_stats(flux_spec, mask=mask_spec, axis=0, sigma=3.0,
                                    cenfunc='median', stdfunc=utils.nan_mad_std)
    # In some cases flux_spec can be totally masked and the result of sigma_clipped_stats is "masked"
    # and that would crush in the following lines
    # TODO investigate and fix this bug
    if flux_mean is np.ma.core.MaskedConstant():
        msgs.info('No objects found')
        # Instantiate a null specobj
        return specobjs.SpecObjs()

    ##   New CODE
    # 1st iteration
    smash_mask = np.isfinite(flux_mean)
    flux_mean_med0 = np.median(flux_mean[smash_mask])
    flux_mean[np.invert(smash_mask)] = flux_mean_med0
    fluxsub0 = flux_mean - flux_mean_med0
    fluxconv0 = scipy.ndimage.filters.gaussian_filter1d(fluxsub0, fwhm / 2.3548, mode='nearest')

    cont_samp = np.fmin(int(np.ceil(nsamp / (fwhm / 2.3548))), 30)
    cont, cont_mask0 = arc.iter_continuum(
        fluxconv0, inmask=smash_mask, fwhm=fwhm, cont_frac_fwhm=2.0, sigthresh=cont_sig_thresh, sigrej=2.0,
        cont_samp=cont_samp,
        npoly=(0 if (nsamp / fwhm < 20.0) else npoly_cont), cont_mask_neg=has_negative, debug=show_cont,
        qa_title='Smash Image Background, 1st iteration: Slit# {:d}'.format(specobj_dict['SLITID']))

    # Second iteration
    flux_mean_med = np.median(flux_mean[cont_mask0])
    fluxsub = flux_mean - flux_mean_med
    fluxconv = scipy.ndimage.filters.gaussian_filter1d(fluxsub, fwhm / 2.3548, mode='nearest')

    cont, cont_mask = arc.iter_continuum(
        fluxconv, inmask=smash_mask, fwhm=fwhm, cont_frac_fwhm=2.0, sigthresh=cont_sig_thresh, sigrej=2.0,
        cont_samp=cont_samp,
        npoly=(0 if (nsamp / fwhm < 20.0) else npoly_cont), cont_mask_neg=has_negative, debug=show_cont,
        qa_title='Smash Image Background: 2nd iteration: Slit# {:d}'.format(specobj_dict['SLITID']))
    fluxconv_cont = (fluxconv - cont) if cont_fit else fluxconv
    # JFH TODO Do we need a running median as was done in the OLD code? Maybe needed for long slits. We could use
    #  use the cont_mask to isolate continuum pixels, and then interpolate the unmasked pixels.
    ##   New CODE

    # TODO: Leave this in!
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
    #                                         cont_mask_neg=has_negative, debug=debug_all)
    #    fluxconv_cont = (fluxconv - cont) if cont_fit else fluxconv
    ## OLD CODE

    if not np.any(cont_mask):
        cont_mask = np.ones(int(nsamp), dtype=bool)  # if all pixels are masked for some reason, don't mask

    mean_sky, med_sky, skythresh = stats.sigma_clipped_stats(fluxconv_cont[cont_mask], sigma=1.5)
    mean, med, sigma = stats.sigma_clipped_stats(fluxconv_cont[cont_mask], sigma=2.5)

    if skythresh == 0.0 and sigma != 0.0:
        skythresh = sigma
    elif skythresh == 0.0 and sigma == 0.0:  # if both SKYTHRESH and sigma are zero mask out the zero pixels and reavaluate
        good = fluxconv_cont > 0.0
        if np.any(good):
            mean_sky, med_sn2_sky, skythresh = stats.sigma_clipped_stats(fluxconv_cont[good], sigma=1.5)
            mean, med_sn2, sigma = stats.sigma_clipped_stats(fluxconv_cont[good], sigma=2.5)
        else:
            msgs.error('Object finding failed. All the elements of the fluxconv_cont spatial profile array are zero')

    # Now find all the peaks without setting any threshold
    ypeak, _, xcen, sigma_pk, _, good_indx, _, _ = arc.detect_lines(fluxconv_cont, cont_subtract=False, fwhm=fwhm,
                                                                    max_frac_fwhm=5.0, input_thresh='None', debug=False)
    ypeak = ypeak[good_indx]
    xcen = xcen[good_indx]
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
        ikeep = (ypeak.argsort()[::-1])[0:nperslit]
        xcen = xcen[ikeep]
        ypeak = ypeak[ikeep]

    npeak = len(xcen)
    # Instantiate a null specobj
    sobjs = specobjs.SpecObjs()
    # Choose which ones to keep and discard based on threshold params. Create SpecObj objects

    # Possible thresholds    [significance,  fraction of brightest, absolute]
    thresh_peak = peak_thresh * ypeak.max() if len(ypeak) > 0 else 0.0
    threshvec = np.array([mean + sig_thresh * sigma, thresh_peak, abs_thresh])
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
            thisobj = specobj.SpecObj(**specobj_dict)
            #
            thisobj.SPAT_FRACPOS = xcen[iobj] / nsamp
            thisobj.smash_peakflux = ypeak[iobj]
            thisobj.smash_nsig = ypeak[iobj] / sigma
            sobjs.add_sobj(thisobj)
    else:
        nobj_reg = 0

    # ToDo Also plot the edge trimming boundaries on the QA here.
    if show_peaks or objfindQA_filename is not None:
        spat_approx_vec = slit_left[specmid] + xsize[specmid] * np.arange(nsamp) / nsamp
        spat_approx = slit_left[specmid] + xsize[specmid] * xcen / nsamp
        # Define the plotting function
        # plt.plot(spat_approx_vec, fluxsub/sigma, color ='cornflowerblue',linestyle=':', label='Collapsed Flux')
        plt.plot(spat_approx_vec, fluxconv_cont / sigma, color='black', label='Collapsed flux (FWHM convol)')
        plt.plot(spat_approx_vec[cont_mask], fluxconv_cont[cont_mask] / sigma, color='red', markersize=3.0,
                 mfc='red', linestyle='None', fillstyle='full',
                 zorder=9, marker='o', label='Used for threshold')
        plt.hlines(threshold / sigma, spat_approx_vec.min(), spat_approx_vec.max(), color='red', linestyle='--',
                   label='Threshold')
        plt.hlines(1.0, spat_approx_vec.min(), spat_approx_vec.max(), color='green', linestyle=':', label='+- 1 sigma')
        plt.hlines(-1.0, spat_approx_vec.min(), spat_approx_vec.max(), color='green', linestyle=':')

        plt.plot(spat_approx, ypeak / sigma, color='red', marker='o', markersize=10.0, mfc='lawngreen',
                 fillstyle='full',
                 linestyle='None', zorder=10, label='Object Found')
        plt.legend()
        plt.xlabel('Approximate Spatial Position (pixels)')
        plt.ylabel('F/sigma (significance)')
        # plt.title(qa_title + ': Slit# {:d}'.format(objfindQA_dict['SLITORD_ID']))
        plt.title(qa_title)
        if objfindQA_filename is not None:
            plt.savefig(objfindQA_filename, dpi=400)
        if show_peaks:
            viewer, ch = display.show_image(image * (thismask * inmask))
            plt.show()
        plt.close('all')

    # Now loop over all the regular apertures and assign preliminary traces to them.
    for iobj in range(nobj_reg):
        # Was a standard trace provided? If so, use that as a crutch.
        if std_trace is not None:
            if iobj == 0:
                msgs.info('Using input STANDARD star trace as crutch for object tracing')
            x_trace = np.interp(specmid, spec_vec, std_trace)
            shift = np.interp(specmid, spec_vec,
                              slit_left + xsize * sobjs[iobj].SPAT_FRACPOS) - x_trace
            sobjs[iobj].TRACE_SPAT = std_trace + shift
        else:  # If no standard is provided shift left slit boundary over to be initial trace
            # ToDO make this the average left and right boundary instead. That would be more robust.
            sobjs[iobj].TRACE_SPAT = slit_left + xsize * sobjs[iobj].SPAT_FRACPOS
        sobjs[iobj].trace_spec = spec_vec
        sobjs[iobj].SPAT_PIXPOS = sobjs[iobj].TRACE_SPAT[specmid]
        # Set the idx for any prelminary outputs we print out. These will be updated shortly
        sobjs[iobj].set_name()

        # assign FWHM
        if use_user_fwhm:
            sobjs[iobj].FWHM = fwhm

        else:
            # Determine the fwhm max
            yhalf = 0.5 * sobjs[iobj].smash_peakflux
            xpk = sobjs[iobj].SPAT_FRACPOS * nsamp
            x0 = int(np.rint(xpk))
            # TODO It seems we have two codes that do similar things, i.e. findfwhm in arextract.py. Could imagine having one
            # Find right location where smash profile croses yhalf
            if x0 < (int(nsamp) - 1):
                ind_righ, = np.where(fluxconv_cont[x0:] < yhalf)
                if len(ind_righ) > 0:
                    i2 = ind_righ[0]
                    if i2 == 0:
                        xrigh = None
                    else:
                        xrigh_int = scipy.interpolate.interp1d(fluxconv_cont[x0 + i2 - 1:x0 + i2 + 1],
                                                               x0 + np.array([i2 - 1, i2], dtype=float),
                                                               assume_sorted=False)
                        xrigh = xrigh_int([yhalf])[0]
                else:
                    xrigh = None
            else:
                xrigh = None
            # Find left location where smash profile crosses yhalf
            if x0 > 0:
                ind_left, = np.where(fluxconv_cont[0:np.fmin(x0 + 1, int(nsamp) - 1)] < yhalf)
                if len(ind_left) > 0:
                    i1 = (ind_left[::-1])[0]
                    if i1 == (int(nsamp) - 1):
                        xleft = None
                    else:
                        xleft_int = scipy.interpolate.interp1d(fluxconv_cont[i1:i1 + 2],
                                                               np.array([i1, i1 + 1], dtype=float), assume_sorted=False)
                        xleft = xleft_int([yhalf])[0]
                else:
                    xleft = None
            else:
                xleft = None

            # Set FWHM for the object
            if (xleft is None) & (xrigh is None):
                fwhm_measure = None
            elif xrigh is None:
                fwhm_measure = 2.0 * (xpk - xleft)
            elif xleft is None:
                fwhm_measure = 2.0 * (xrigh - xpk)
            else:
                fwhm_measure = (xrigh - xleft)

            if fwhm_measure is not None:
                sobjs[iobj].FWHM = np.sqrt(
                    np.fmax(fwhm_measure ** 2 - fwhm ** 2, (fwhm / 2.0) ** 2))  # Set a floor of fwhm/2 on fwhm
            else:
                sobjs[iobj].FWHM = fwhm

        # assign BOX_RADIUS
        sobjs[iobj].BOX_RADIUS = boxcar_rad

    if (len(sobjs) == 0) & (hand_extract_dict is None):
        msgs.info('No objects found')
        return specobjs.SpecObjs()
    else:
        msgs.info("Automatic finding routine found {0:d} objects".format(len(sobjs)))

    msgs.info('Fitting the object traces')

    if len(sobjs) > 0:
        # Note the transpose is here to pass in the TRACE_SPAT correctly.
        xinit_fweight = np.copy(sobjs.TRACE_SPAT.T)
        spec_mask = (spec_vec >= spec_min_max[0]) & (spec_vec <= spec_min_max[1])
        trc_inmask = np.outer(spec_mask, np.ones(len(sobjs), dtype=bool))
        xfit_fweight = fit_trace(image, xinit_fweight, ncoeff, bpm=np.invert(inmask),
                                 trace_bpm=np.invert(trc_inmask), fwhm=fwhm, maxdev=maxdev,
                                 idx=sobjs.NAME, debug=show_fits)[0]
        xinit_gweight = np.copy(xfit_fweight)
        xfit_gweight = fit_trace(image, xinit_gweight, ncoeff, bpm=np.invert(inmask),
                                 trace_bpm=np.invert(trc_inmask), fwhm=fwhm, maxdev=maxdev,
                                 weighting='gaussian', idx=sobjs.NAME, debug=show_fits)[0]

        # assign the final trace
        for iobj in range(nobj_reg):
            sobjs[iobj].TRACE_SPAT = xfit_gweight[:, iobj]
            sobjs[iobj].SPAT_PIXPOS = sobjs[iobj].TRACE_SPAT[specmid]
            sobjs[iobj].set_name()

    # Now deal with the hand apertures if a hand_extract_dict was passed in. Add these to the SpecObj objects
    if hand_extract_dict is not None:
        # First Parse the hand_dict
        hand_extract_spec, hand_extract_spat, hand_extract_det, hand_extract_fwhm = [
            hand_extract_dict[key] for key in ['spec', 'spat', 'det', 'fwhm']]

        # Determine if these hand apertures land on the slit in question
        hand_on_slit = np.where(np.array(thismask[np.rint(hand_extract_spec).astype(int),
                                                  np.rint(hand_extract_spat).astype(int)]))
        hand_extract_spec = hand_extract_spec[hand_on_slit]
        hand_extract_spat = hand_extract_spat[hand_on_slit]
        hand_extract_det = hand_extract_det[hand_on_slit]
        hand_extract_fwhm = hand_extract_fwhm[hand_on_slit]
        nobj_hand = len(hand_extract_spec)
        msgs.info("Implementing hand apertures for {} sources on the slit".format(nobj_hand))

        # Decide how to assign a trace to the hand objects
        if nobj_reg > 0:  # Use brightest object on slit?
            smash_peakflux = sobjs.smash_peakflux
            ibri = smash_peakflux.argmax()
            trace_model = sobjs[ibri].TRACE_SPAT
            med_fwhm_reg = np.median(sobjs.FWHM)
        elif std_trace is not None:  # If no objects found, use the standard?
            trace_model = std_trace
        else:  # If no objects or standard use the slit boundary
            msgs.warn("No source to use as a trace.  Using the slit boundary")
            trace_model = slit_left

        # Loop over hand_extract apertures and create and assign specobj
        for iobj in range(nobj_hand):
            # Proceed
            thisobj = specobj.SpecObj(**specobj_dict)
            thisobj.hand_extract_spec = hand_extract_spec[iobj]
            thisobj.hand_extract_spat = hand_extract_spat[iobj]
            thisobj.hand_extract_det = hand_extract_det[iobj]
            thisobj.hand_extract_fwhm = hand_extract_fwhm[iobj]
            thisobj.hand_extract_flag = True
            # SPAT_FRACPOS
            f_ximg = scipy.interpolate.RectBivariateSpline(spec_vec, spat_vec, ximg)
            thisobj.SPAT_FRACPOS = float(
                f_ximg(thisobj.hand_extract_spec, thisobj.hand_extract_spat, grid=False))  # interpolate from ximg
            thisobj.smash_peakflux = np.interp(thisobj.SPAT_FRACPOS * nsamp, np.arange(nsamp),
                                               fluxconv_cont)  # interpolate from fluxconv
            # assign the trace
            spat_0 = np.interp(thisobj.hand_extract_spec, spec_vec, trace_model)
            shift = thisobj.hand_extract_spat - spat_0
            thisobj.TRACE_SPAT = trace_model + shift
            thisobj.trace_spec = spec_vec
            thisobj.SPAT_PIXPOS = thisobj.TRACE_SPAT[specmid]
            thisobj.set_name()
            # assign FWHM
            # TODO -- I think FWHM *has* to be input
            if hand_extract_fwhm[iobj] is not None:  # If a hand_extract_fwhm was input use that for the fwhm
                thisobj.FWHM = hand_extract_fwhm[iobj]
            elif nobj_reg > 0:  # Otherwise is None was input, then use the median of objects on this slit if they are present
                thisobj.FWHM = med_fwhm_reg
            else:  # Otherwise just use the FWHM parameter input to the code (or the default value)
                thisobj.FWHM = fwhm
            # assign BOX_RADIUS
            thisobj.BOX_RADIUS = boxcar_rad
            # Finish
            sobjs.add_sobj(thisobj)

    nobj = len(sobjs)

    ## Okay now loop over all the regular aps and exclude any which within the fwhm of the hand_extract_APERTURES
    if nobj_reg > 0 and hand_extract_dict is not None:
        spat_pixpos = sobjs.SPAT_PIXPOS
        hand_flag = sobjs.hand_extract_flag
        spec_fwhm = sobjs.FWHM
        # spat_pixpos = np.array([spec.SPAT_PIXPOS for spec in specobjs])
        # hand_flag = np.array([spec.hand_extract_flag for spec in specobjs])
        # spec_fwhm = np.array([spec.FWHM for spec in specobjs])
        reg_ind, = np.where(np.invert(hand_flag))
        hand_ind, = np.where(hand_flag)
        # med_fwhm = np.median(spec_fwhm[~hand_flag])
        # spat_pixpos_hand = spat_pixpos[hand_ind]
        keep = np.ones(nobj, dtype=bool)
        for ihand in hand_ind:
            close = np.abs(sobjs[reg_ind].SPAT_PIXPOS - spat_pixpos[ihand]) <= 0.6 * spec_fwhm[ihand]
            if np.any(close):
                # Print out a warning
                msgs.warn('Deleting object(s) {}'.format(sobjs[reg_ind[close]].NAME) +
                          ' because it collides with a user specified hand_extract aperture')
                keep[reg_ind[close]] = False

        sobjs = sobjs[keep]

    if len(sobjs) == 0:
        msgs.info('No hand or normal objects found on this slit. Returning')
        return specobjs.SpecObjs()

    # Sort objects according to their spatial location
    nobj = len(sobjs)
    spat_pixpos = sobjs.SPAT_PIXPOS
    sobjs = sobjs[spat_pixpos.argsort()]
    # Assign integer objids
    sobjs.OBJID = np.arange(nobj) + 1

    # Assign the maskwidth and compute some inputs for the object mask
    for iobj in range(nobj):
        # TODO -- This parameter may not be used anywhere
        if skythresh > 0.0:
            sobjs[iobj].maskwidth = extract_maskwidth * sobjs[iobj].FWHM * (
                        1.0 + 0.5 * np.log10(np.fmax(sobjs[iobj].smash_peakflux / skythresh, 1.0)))
        else:
            sobjs[iobj].maskwidth = extract_maskwidth * sobjs[iobj].FWHM

    # If requested display the resulting traces on top of the image
    if show_trace:
        viewer, ch = display.show_image(image * (thismask * inmask))
        display.show_slits(viewer, ch, slit_left.T, slit_righ.T, slit_ids=sobjs[0].SLITID)
        for iobj in range(nobj):
            if sobjs[iobj].hand_extract_flag == False:
                color = 'orange'
            else:
                color = 'blue'
            display.show_trace(viewer, ch, sobjs[iobj].TRACE_SPAT, trc_name=sobjs[iobj].NAME, color=color)

    msgs.info("Successfully traced a total of {0:d} objects".format(len(sobjs)))

    # Finish
    for sobj in sobjs:
        # Add in more info
        sobj.THRESHOLD = threshold
        # Vet
        if not sobj.ready_for_extraction():
            # embed(header=utils.embed_header())
            msgs.error("Bad SpecObj.  Can't proceed")

    # Return
    return sobjs


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



