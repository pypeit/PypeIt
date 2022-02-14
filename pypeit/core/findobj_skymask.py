""" Methods to find objects

.. include:: ../include/links.rst

"""
import copy

import numpy as np
import scipy
from matplotlib import pyplot as plt

from astropy import stats

from pypeit import msgs
from pypeit import utils
from pypeit import specobj
from pypeit import specobjs
from pypeit.core import pydl
from pypeit.core import fitting
from pypeit.core.moment import moment1d
from pypeit import tracepca
from pypeit.core.trace import fit_trace
from pypeit.core import arc
from pypeit.display import display
from pypeit.core import pixels

from IPython import embed


def create_skymask(sobjs, thismask, slit_left, slit_righ, box_rad_pix=None, trim_edg=(5,5), skymask_nthresh=1.0):
    r"""
    Creates a skymask from a SpecObjs object using the fwhm of each object
    and or the boxcar radius

    Args:
        sobjs (:class:`pypeit.specobjs.SpecObjs`):
            Objects for which you would like to create the mask
        thismask (`numpy.ndarray`_):
            Boolean image indicating pixels which are on the slit.
            Shape is :math:`(N_{\rm spec}, N_{\rm spat})`.
        slit_left (`numpy.ndarray`_):
            Left boundary of slit/order to be extracted (given as
            floating pt pixels). This a 1-d array with shape :math:`(N_{\rm spec}, 1)`
            or :math:`(N_{\rm spec},)`
        slit_righ (`numpy.ndarray`_):
            Right boundary of slit/order to be extracted (given as
            floating pt pixels). This a 1-d array with shape :math:`(N_{\rm spec}, 1)`
            or :math:`(N_{\rm spec},)`
        box_rad_pix (:obj:`float`, optional):
            If set, the skymask will be as wide as this radius in pixels.
        skymask_nthresh (:obj:`float`, optional): default = 2.0
            The multiple of the final object finding threshold (see
            above) which is used to create the skymask using the value
            of the peak flux in the slit profile (image with the
            spectral direction smashed out).
        trim_edg (:obj:`tuple`, optional): of integers or float, default = (5,5)
            Ignore objects within this many pixels of the left and right
            slit boundaries, where the first element refers to the left
            and second refers to the right.

    Returns:
        `numpy.ndarray`_: Boolean image with shape :math:`(N_{\rm spec}, N_{\rm spat})`
            (same as thismask) indicating which pixels are usable for
            global sky subtraction.  True = usable for sky subtraction,
            False = should be masked when sky subtracting.
    """
    nobj = len(sobjs)
    ximg, _ = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)
    # How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    #nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())

    # Objmask
    skymask_objflux = np.copy(thismask)
    if nobj == 0:
        msgs.info('No objects were detected. The entire slit will be used to determine the sky subtraction.')
    else:
        # Compute some inputs for the object mask
        xtmp = (np.arange(nsamp) + 0.5)/nsamp
        # threshold for object finding
        for iobj in range(nobj):
            # this will skip also sobjs with THRESHOLD=0, because are the same that have smash_peakflux=0.
            if (sobjs[iobj].smash_peakflux != 0.) and (sobjs[iobj].smash_peakflux != None):
                qobj = np.zeros_like(xtmp)
                sep = np.abs(xtmp-sobjs[iobj].SPAT_FRACPOS)
                sep_inc = sobjs[iobj].maskwidth/nsamp
                close = sep <= sep_inc
                qobj[close] = sobjs[iobj].smash_peakflux * \
                               np.exp(np.fmax(-2.77*(sep[close]*nsamp)**2/sobjs[iobj].FWHM**2, -9.0))
                if sobjs[iobj].THRESHOLD > 0.:
                    skymask_objflux[thismask] &= \
                        np.interp(ximg[thismask], xtmp, qobj) < (skymask_nthresh * sobjs[iobj].THRESHOLD)
    # FWHM
    skymask_fwhm = np.copy(thismask)
    if nobj > 0:
        nspec, nspat = thismask.shape
        # spatial position everywhere along image
        spat_img = np.outer(np.ones(nspec, dtype=int),np.arange(nspat, dtype=int))
        # Boxcar radius?
        if box_rad_pix is not None:
            msgs.info("Using boxcar radius for masking")
        # Loop me
        for iobj in range(nobj):
            # Create a mask for the pixels that will contribute to the object
            skymask_radius = box_rad_pix if box_rad_pix is not None else sobjs[iobj].FWHM
            msgs.info(f"Masking around object {iobj+1} within a radius = {skymask_radius} pixels")
            slit_img = np.outer(sobjs[iobj].TRACE_SPAT, np.ones(nspat))  # central trace replicated spatially
            objmask_now = thismask & (spat_img > (slit_img - skymask_radius)) & (spat_img < (slit_img + skymask_radius))
            skymask_fwhm &= np.invert(objmask_now)

        # Check that we have not performed too much masking
        if (np.sum(skymask_fwhm)/np.sum(thismask) < 0.10):
            msgs.warn('More than 90% of  usable area on this slit would be masked and not used by global sky subtraction. '
                      'Something is probably wrong with object finding for this slit. Not masking object for global sky subtraction.')
            skymask_fwhm = np.copy(thismask)

    # Still have to make the skymask
    # # TODO -- Make sure this is right
    # if box_rad_pix is None:
    #     skymask = skymask_objflux | skymask_fwhm
    # else:  # Enforces boxcar radius masking
    #     skymask = skymask_objflux & skymask_fwhm
    # DP: I think skymask should always be skymask_objflux & skymask_fwhm (i.e., not only when box_rad_pix is not None).
    # In the case of skymask_objflux | skymask_fwhm, if skymask_objflux cannot be computed, the entire slit
    # is used for sky calculation (i.e., skymask_fwhm will not have effect).

    # DP's change which I don't think we should adopt at this time.
    #skymask = skymask_objflux & skymask_fwhm

    # JFH restored old behavior after seeing spurious results for X-shooter. I think the issue here is that the fwhm
    # computation from objs_in_slit is not necessarily that reliable and when large amounts of masking are performed
    # on narrow slits/orders, we have problems. We should revisit this after object finding is refactored since
    # maybe then the fwhm estimates will be more robust.
    if box_rad_pix is None and np.all([sobj.smash_peakflux is not None for sobj in sobjs]) \
            and np.all([sobj.smash_peakflux != 0. for sobj in sobjs]):
        # TODO This is a kludge until we refactor this routine. Basically mask design objects that are not auto-ID
        # always have smash_peakflux undefined. If there is a hybrid situation of auto-ID and maskdesign, the logic
        # here does not really make sense. Soution would be to compute thershold and smash_peakflux for all objects.
        skymask = skymask_objflux | skymask_fwhm
    else:  # Enforces boxcar radius masking
        skymask = skymask_objflux & skymask_fwhm

    # Return
    return skymask[thismask]


def ech_objfind(image, ivar, slitmask, slit_left, slit_righ, order_vec, maskslits, det=1,
                inmask=None, spec_min_max=None, fof_link=1.5, plate_scale=0.2, has_negative=False,
                std_trace=None, ncoeff=5, npca=None, coeff_npoly=None, max_snr=2.0, min_snr=1.0,
                nabove_min_snr=2, pca_explained_var=99.0, box_radius=2.0, fwhm=3.0,
                use_user_fwhm=False, maxdev=2.0, hand_extract_dict=None, nperorder=2,
                extract_maskwidth=3.0, sig_thresh=10.0, peak_thresh=0.0, abs_thresh=0.0,
                cont_sig_thresh=2.0, specobj_dict=None, trim_edg=(5,5), cont_fit=True,
                npoly_cont=1, show_peaks=False, show_fits=False, show_single_fits=False,
                show_trace=False, show_single_trace=False, show_pca=False,
                debug_all=False, objfindQA_filename=None):
    """
    Object finding routine for Echelle spectrographs. This routine:
       1) runs object finding on each order individually
       2) Links the objects found together using a friends-of-friends algorithm on fractional order position.
       3) For objects which were only found on some orders, the standard (or the slit boundaries) are placed at the appropriate
          fractional position along the order.
       4) A PCA fit to the traces is performed using the routine above pca_fit

    Args:
        image (`numpy.ndarray`_):
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial. Note this
            image can either have the sky background in it, or have
            already been sky subtracted.  Object finding works best on
            sky-subtracted images. Ideally objfind would be run in
            another routine, global sky-subtraction performed, and then
            this code should be run. However, it is also possible to run
            this code on non sky subtracted images.
        ivar (`numpy.ndarray`_): float ndarray, shape (nspec, nspat)
            Floating-point inverse variance image for the input image.
            Shape is (nspec, nspat).
        slitmask (`numpy.ndarray`_):
            Integer image indicating the pixels that belong to each
            order. Pixels that are not on an order have value -1, and
            those that are on an order have a value equal to the slit
            number (i.e. 0 to nslits-1 from left to right on the image).
            Shape is (nspec, nspat).
        slit_left (`numpy.ndarray`_):
            Left boundary of orders to be extracted (given as floating
            pt pixels). This a 2-d float array with shape (nspec, norders)
        slit_righ (`numpy.ndarray`_):
            Left boundary of orders to be extracted (given as floating
            pt pixels). This a 2-d float array with shape (nspec, norders)
        order_vec (`numpy.ndarray`_):
            Echelle orders.  This is written to the SpecObj objects.
            It is ok, but not recommended to provide np.arange(norders)
        maskslits (`numpy.ndarray`_):
        det (:obj:`int`, optional):
            Need for hand object
        inmask (`numpy.ndarray`_):
            Boolean input mask for the input image. Shape is (nspec, nspat).
        fwhm (:obj:`float`):
            Estimated fwhm of the objects in pixels
        use_user_fwhm (:obj:`bool`):
            If True PypeIt will use the spatial profile fwm input by the user (i.e. the fwhm parameter above)
            rather than determine the spatial fwhm from the smashed spatial profile via the automated algorithm.
            Default = False.
        hand_extract_dict (:obj:`dict`, optional):
            Dictionary with info on manual extraction.
        maxdev (:obj:`float`):
            Maximum deviation of pixels from polynomial fit to trace
            used to reject bad pixels in trace fitting.
        spec_min_max (`numpy.ndarray`_, optional):
            This is a 2-d array of shape (2, norders) which defines the minimum and maximum of
            each order in the spectral direction on the detector. This
            should only be used for echelle spectrographs for which the
            orders do not entirely cover the detector. The pca_trace
            code will re-map the traces such that they all have the same
            length, compute the PCA, and then re-map the orders back.
            This improves performanc for echelle spectrographs by
            removing the nonlinear shrinking of the orders so that the
            linear pca operation can better predict the traces. If not
            passed in it will be determined automitically from the
            slitmask.
        fof_link (:obj:`float`):
            Friends-of-friends linking length in arcseconds used to link
            together traces across orders. The routine links together at
            the same fractional slit position and links them together
            with a friends-of-friends algorithm using this linking
            length.
        plate_scale (:obj:`float`, `numpy.ndarray`_):
            Plate scale of your detector, in unit of arcsec/pix. This
            can either be a single float for every order, or an array
            with shape (norders,) indicating the plate scale of each order.
        ncoeff (:obj:`int`):
            Order of polynomial fit to traces.
        npca (:obj:`int`):
            Nmber of PCA components you want to keep. default is None
            and it will be assigned automatically by calculating the
            number of components contains approximately 99% of the
            variance. Default = None
        coeff_npoly (:obj:`int`):
            order of polynomial used for PCA coefficients fitting.
            Default is None and this will be determined automatically.
        min_snr (:obj:`float`):
            Minimum SNR for keeping an object. For an object to be kept
            it must have a median S/N ratio above min_snr for at least
            nabove_min_snr orders.
        max_snr (:obj:`float`):
            Required SNR for keeping an object. For an object to be kept
            it must have a max S/N ratio above max_snr.
        nabove_min_snr (:obj:`int`):
            The required number of orders that an object must have with
            median SNR>min_snr in order to be kept.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance
            in the data accounted for by the PCA used to truncate the
            number of PCA coefficients to keep (see `npca`). Ignored
            if `npca` is provided directly. See :func:`pypeit.core.pca.pca_decomposition`.
        cont_sig_thresh (:obj:`float`, optional):
            Significance threshold for peak detection for determinining which pixels to use for the iteratively
            fit continuum of the spectral direction smashed image. This is passed to objfind which is then passed as
            the sigthresh parameter to core.arc.iter_continum. For extremely narrow slits that are almost filled by
            the object trace set this to a smaller number like 1.0 or disable continuum fitting altogether with
            cont_fit=False below. Default = 1.0
        npoly_cont (:obj:`int`):
            Order of polynomial fit to the illumination pattern across the slit when peak finding
        cont_fit (:obj:`bool`):
            Fit a continuum to the illumination pattern across the slit when peak finding
        trim_edg (:obj:`tuple`):
            Ignore objects within this many pixels of the left and right
            slit boundaries, where the first element refers to the left
            and second refers to the right. This is tuple of 2 integers of floats
        specobj_dict (:obj:`dict`):
            Dictionary containing meta-data for the objects that will be
            propgated into the SpecObj objects, i.e. SLITID,
            detector, object type, and pipeline. The default is None, in
            which case the following dictionary will be used::

                specobj_dict = {'SLITID': 999, 'DET': 'DET01',
                                'OBJTYPE': 'unknown', 'PYPELINE': 'unknown'}
        abs_thresh (:obj:`float`):
            Absolute threshold for object detection.  Objects are found
            by smashing out the spectral direction along the curved
            slit/order traces, and abs_thresh is in the units of this
            smashed profile.  The code uses the maximum of the
            thresholds defined by sig_thresh, peak_thers, and
            abs_thresh.
        peak_thresh (:obj:`float`):
            Peak threshold for object detection. This is a number
            between 0 and 1 and represents the fraction of the brightest
            object on the slit that will be kept as an object, i.e. if
            ymax is the brightest object of the spectrum smashed out in
            the spectral direction, all objects with ypeak >
            peak_thresh*ymak are kept. The code uses the maximum of the
            thresholds defined by sig_thresh, peak_thers, and
            abs_thresh.
        extract_maskwidth (:obj:`float`,optional):
            This parameter determines the initial size of the region in
            units of fwhm that will be used for local sky subtraction in
            the routine skysub.local_skysub_extract.
        nperorder (:obj:`int`):
            Maximum number of objects allowed per order. The code will
            take the nperorder most significant detections. However hand apertures will always be returned
            and do not count against this budget.
        std_trace (`numpy.ndarray`_):
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        has_negative (:obj:`bool`, optional):
            Image has negative object traces, i.e. for IR difference imaging. This impacts how the
            iterative conntinuum is fit to the spectral direction smashed image for object finding. Default=False
        box_radius (:obj:`float`):
            Box_car extraction radius in arcseconds to assign to each detected object and to be
            used later for boxcar extraction. In this method box_radius is converted into pixels
            by using the plate scale for the particular order.
            box_radius is also used for SNR calculation and trimming.
        sig_thresh (:obj:`float`):
            Threshold for finding objects
        show_peaks (:obj:`bool`):
            Whether plotting the QA of peak finding of your object in each order
        show_fits (:obj:`bool`):
            Plot trace fitting for final fits using PCA as crutch
        show_single_fits (:obj:`bool`):
            Plot trace fitting for single order fits
        show_trace (:obj:`bool`):
            Whether display the resulting traces on top of the image
        show_single_trace (:obj:`bool`):
            Whether display the resulting traces on top of the single order
        show_pca (:obj:`bool`):
            Whether display debugging plots for pca
        debug_all (:obj:`bool`):
            Show all the debugging plots?
        objfindQA_filename (:obj:`str`, optional):
            Directory + filename of the object profile QA. Default = None.

    Returns:
        :class:`pypeit.specobjs.SpecObjs`: class containing the objects detected
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
        specobj_dict = {'SLITID': 999, 'ECH_ORDERINDX': 999,
                        'DET': det, 'OBJTYPE': 'unknown', 'PYPELINE': 'Echelle'}


    # TODO Update FOF algorithm here with the one from scikit-learn.

    allmask = slitmask > -1
    if inmask is None:
        inmask = allmask

    nspec, nspat = image.shape
    norders = len(order_vec)

    # Find the spat IDs
    gdslit_spat = np.unique(slitmask[slitmask >= 0]).astype(int)  # Unique sorts
    if gdslit_spat.size != norders:
        msgs.error('There is a mismatch between the number of valid orders found by PypeIt and '
                   'the number expected for this spectrograph.  Unable to continue.  Please '
                   'submit an issue on Github: https://github.com/pypeit/PypeIt/issues .')

    if spec_min_max is None:
        spec_min_max = np.zeros((2,norders))
        for iord in range(norders):
            ispec, ispat = np.where(slitmask == gdslit_spat[iord])
            spec_min_max[:,iord] = ispec.min(), ispec.max()

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

    # TODO JFH This hand apertures in echelle needs to be completely refactored.
    # Hand prep
    #   Determine the location of the source on *all* of the orders
    if hand_extract_dict is not None:
        f_spats = []
        for ss, spat, spec in zip(range(len(hand_extract_dict['spec'])),
                                  hand_extract_dict['spat'],
                                  hand_extract_dict['spec']):
            # Find the input slit
            ispec = int(np.clip(np.round(spec),0,nspec-1))
            ispat = int(np.clip(np.round(spat),0,nspec-1))
            slit = slitmask[ispec, ispat]
            # Fractions
            iord_hand = gdslit_spat.tolist().index(slit)
            f_spat = (spat - slit_left[ispec, iord_hand]) / (
                slit_righ[ispec, iord_hand] - slit_left[ispec, iord_hand])
            f_spats.append(f_spat)

    # Loop over orders and find objects
    sobjs = specobjs.SpecObjs()
    # ToDo replace orderindx with the true order number here? Maybe not. Clean up SLITID and orderindx!
    gdorders = np.arange(norders)[np.invert(maskslits)]
    for iord in gdorders: #range(norders):
        qa_title = 'Finding objects on order # {:d}'.format(order_vec[iord])
        msgs.info(qa_title)
        thisslit_gpm = slitmask == gdslit_spat[iord]
        inmask_iord = inmask & thisslit_gpm
        specobj_dict['SLITID'] = gdslit_spat[iord]
        specobj_dict['ECH_ORDERINDX'] = iord
        specobj_dict['ECH_ORDER'] = order_vec[iord]
        std_in = None if std_trace is None else std_trace[:, iord]

        # TODO JFH: Fix this. The way this code works, you should only need to create a single hand object,		
        # not one at every location on the order            
        if hand_extract_dict is not None:
            new_hand_extract_dict = copy.deepcopy(hand_extract_dict)
            for ss, spat, spec, f_spat in zip(range(len(hand_extract_dict['spec'])),
                                              hand_extract_dict['spat'],
                                              hand_extract_dict['spec'], f_spats):
                ispec = int(spec)
                new_hand_extract_dict['spec'][ss] = ispec
                new_hand_extract_dict['spat'][ss] = slit_left[ispec,iord] + f_spat*(
                    slit_righ[ispec,iord]-slit_left[ispec,iord])
        else:
            new_hand_extract_dict = None

        # Get SLTIORD_ID for the objfind QA
        ech_objfindQA_filename = objfindQA_filename.replace('S0999', 'S{:04d}'.format(order_vec[iord])) \
            if objfindQA_filename is not None else None
        # Run
        sobjs_slit = \
            objs_in_slit(image, thisslit_gpm, slit_left[:,iord], slit_righ[:,iord], spec_min_max=spec_min_max[:,iord],
                    inmask=inmask_iord,std_trace=std_in, ncoeff=ncoeff, fwhm=fwhm, use_user_fwhm=use_user_fwhm, maxdev=maxdev,
                    hand_extract_dict=new_hand_extract_dict, has_negative=has_negative,
                    nperslit=nperorder, extract_maskwidth=extract_maskwidth, sig_thresh=sig_thresh,
                    peak_thresh=peak_thresh, abs_thresh=abs_thresh, cont_sig_thresh=cont_sig_thresh,
                    trim_edg=trim_edg, boxcar_rad=box_radius/plate_scale_ord[iord], cont_fit=cont_fit,
                    npoly_cont=npoly_cont, show_peaks=show_peaks, show_fits=show_single_fits,
                    show_trace=show_single_trace, qa_title=qa_title, specobj_dict=specobj_dict,
                    objfindQA_filename=ech_objfindQA_filename)
        sobjs.add_sobj(sobjs_slit)

    nfound = len(sobjs)

    if nfound == 0:
        msgs.warn('No objects found')
        return sobjs

    FOF_frac = fof_link/(np.median(np.median(slit_width,axis=0)*plate_scale_ord))
    # Run the FOF. We use fake coordinates
    fracpos = sobjs.SPAT_FRACPOS
    ra_fake = fracpos/1000.0  # Divide all angles by 1000 to make geometry euclidian
    dec_fake = np.zeros_like(fracpos)
    if nfound>1:
        inobj_id, multobj_id, firstobj_id, nextobj_id \
                = pydl.spheregroup(ra_fake, dec_fake, FOF_frac/1000.0)
        # TODO spheregroup returns zero based indices but we use one based. We should probably add 1 to inobj_id here,
        # i.e. obj_id_init = inobj_id + 1
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
            on_order = (obj_id_init == uni_obj_id_init[iobj]) & (sobjs.ECH_ORDERINDX == iord)
            if (np.sum(on_order) > 1):
                msgs.warn('Found multiple objects in a FOF group on order iord={:d}'.format(iord) + msgs.newline() +
                          'Spawning new objects to maintain a single object per order.')
                off_order = (obj_id_init == uni_obj_id_init[iobj]) & (sobjs.ECH_ORDERINDX != iord)
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
                # JFH OLD LINE with bug
                #obj_id[ind_rest] = (np.arange(len(ind_rest)) + 1) + obj_id_init.max()
                obj_id[ind_rest] = (np.arange(len(ind_rest)) + 1) + obj_id.max()

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
            on_order = (obj_id == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDERINDX == iord)
            sobjs_align[on_order].ECH_FRACPOS = uni_frac[iobj]
            sobjs_align[on_order].ECH_OBJID = uni_obj_id[iobj]
            sobjs_align[on_order].OBJID = uni_obj_id[iobj]
            sobjs_align[on_order].ech_frac_was_fit = False

    # Reset names (just in case)
    sobjs_align.set_names()
    # Now loop over objects and fill in the missing objects and their traces. We will fit the fraction slit position of
    # the good orders where an object was found and use that fit to predict the fractional slit position on the bad orders
    # where no object was found
    for iobj in range(nobj):
        # Grab all the members of this obj_id from the object list
        indx_obj_id = sobjs_align.ECH_OBJID == uni_obj_id[iobj]
        nthisobj_id = np.sum(indx_obj_id)
        # Perform the fit if this objects shows up on more than three orders
        if (nthisobj_id > 3) and (nthisobj_id<norders):
            thisorderindx = sobjs_align[indx_obj_id].ECH_ORDERINDX
            goodorder = np.zeros(norders, dtype=bool)
            goodorder[thisorderindx] = True
            badorder = np.invert(goodorder)
            xcen_good = (sobjs_align[indx_obj_id].TRACE_SPAT).T
            slit_frac_good = (xcen_good-slit_left[:,goodorder])/slit_width[:,goodorder]
            # Fractional slit position averaged across the spectral direction for each order
            frac_mean_good = np.mean(slit_frac_good, 0)
            # Perform a  linear fit to fractional slit position
            #TODO Do this as a S/N weighted fit similar to what is now in the pca_trace algorithm?
            #msk_frac, poly_coeff_frac = fitting.robust_fit(order_vec[goodorder], frac_mean_good, 1,
            pypeitFit = fitting.robust_fit(order_vec[goodorder], frac_mean_good, 1,
                                           function='polynomial', maxiter=20, lower=2, upper=2,
                                           use_mad= True, sticky=False,
                                           minx = order_vec.min(), maxx=order_vec.max())
            frac_mean_new = np.zeros(norders)
            frac_mean_new[badorder] = pypeitFit.eval(order_vec[badorder])#, minx = order_vec.min(),maxx=order_vec.max())
            frac_mean_new[goodorder] = frac_mean_good
            # TODO This QA needs some work
            if show_pca:
                frac_mean_fit = pypeitFit.eval(order_vec)
                plt.plot(order_vec[goodorder][pypeitFit.bool_gpm], frac_mean_new[goodorder][pypeitFit.bool_gpm], 'ko', mfc='k', markersize=8.0, label='Good Orders Kept')
                plt.plot(order_vec[goodorder][np.invert(pypeitFit.bool_gpm)], frac_mean_new[goodorder][np.invert(pypeitFit.bool_gpm)], 'ro', mfc='k', markersize=8.0, label='Good Orders Rejected')
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
            on_order = (sobjs_align.ECH_OBJID == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDERINDX == iord)
            num_on_order = np.sum(on_order)
            if num_on_order == 0:
                # If it is not, create a new sobjs and add to sobjs_align and assign required tags
                thisobj = specobj.SpecObj('Echelle', sobjs_align[0].DET,
                                             OBJTYPE=sobjs_align[0].OBJTYPE,
                                             ECH_ORDERINDX=iord,
                                             ECH_ORDER=order_vec[iord])
                #thisobj.ECH_ORDERINDX = iord
                #thisobj.ech_order = order_vec[iord]
                thisobj.SPAT_FRACPOS = uni_frac[iobj]
                # Assign traces using the fractional position fit above
                if std_trace is not None:
                    x_trace = np.interp(slit_spec_pos, spec_vec, std_trace[:,iord])
                    shift = np.interp(slit_spec_pos, spec_vec,slit_left[:,iord] + slit_width[:,iord]*frac_mean_new[iord]) - x_trace
                    thisobj.TRACE_SPAT = std_trace[:,iord] + shift
                else:
                    thisobj.TRACE_SPAT = slit_left[:, iord] + slit_width[:, iord] * frac_mean_new[iord]  # new trace
                thisobj.trace_spec = spec_vec
                thisobj.SPAT_PIXPOS = thisobj.TRACE_SPAT[specmid]
                # Use the real detections of this objects for the FWHM
                this_obj_id = obj_id == uni_obj_id[iobj]
                # Assign to the fwhm of the nearest detected order
                imin = np.argmin(np.abs(sobjs_align[this_obj_id].ECH_ORDERINDX - iord))
                thisobj.FWHM = sobjs_align[imin].FWHM
                thisobj.maskwidth = sobjs_align[imin].maskwidth
                thisobj.smash_peakflux = sobjs_align[imin].smash_peakflux
                thisobj.THRESHOLD = sobjs_align[imin].THRESHOLD
                thisobj.BOX_RADIUS = sobjs_align[imin].BOX_RADIUS
                thisobj.ECH_FRACPOS = uni_frac[iobj]
                thisobj.ECH_OBJID = uni_obj_id[iobj]
                thisobj.OBJID = uni_obj_id[iobj]
                thisobj.SLITID = gdslit_spat[iord]
                thisobj.ech_frac_was_fit = True
                thisobj.set_name()
                sobjs_align.add_sobj(thisobj)
                obj_id = np.append(obj_id, uni_obj_id[iobj])
                gfrac = np.append(gfrac, uni_frac[iobj])
            elif num_on_order == 1:
                # Object is already on this order so no need to do anything
                pass
            elif num_on_order > 1:
                msgs.error('Problem in echelle object finding. The same objid={:d} appears {:d} times on echelle orderindx ={:d}'
                           ' even after duplicate obj_ids the orders were removed. '
                           'Report this bug to PypeIt developers'.format(uni_obj_id[iobj],num_on_order, iord))



    # Loop over the objects and perform a quick and dirty extraction to assess S/N.
    varimg = utils.calc_ivar(ivar)
    flux_box = np.zeros((nspec, norders, nobj))
    ivar_box = np.zeros((nspec, norders, nobj))
    mask_box = np.zeros((nspec, norders, nobj))
    SNR_arr = np.zeros((norders, nobj))
    slitfracpos_arr = np.zeros((norders, nobj))
    for iobj in range(nobj):
        for iord in range(norders):
            iorder_vec = order_vec[iord]
            indx = sobjs_align.slitorder_objid_indices(iorder_vec, uni_obj_id[iobj])
            #indx = (sobjs_align.ECH_OBJID == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDERINDX == iord)
            #spec = sobjs_align[indx][0]
            inmask_iord = inmask & (slitmask == gdslit_spat[iord])
            # TODO make the snippet below its own function quick_extraction()
            box_rad_pix = box_radius/plate_scale_ord[iord]

            # TODO -- We probably shouldn't be operating on a SpecObjs but instead a SpecObj
            flux_tmp  = moment1d(image*inmask_iord, sobjs_align[indx][0].TRACE_SPAT, 2*box_rad_pix,
                                 row=sobjs_align[indx][0].trace_spec)[0]
            var_tmp  = moment1d(varimg*inmask_iord, sobjs_align[indx][0].TRACE_SPAT, 2*box_rad_pix,
                                row=sobjs_align[indx][0].trace_spec)[0]
            ivar_tmp = utils.calc_ivar(var_tmp)
            pixtot  = moment1d(ivar*0 + 1.0, sobjs_align[indx][0].TRACE_SPAT, 2*box_rad_pix,
                               row=sobjs_align[indx][0].trace_spec)[0]
            mask_tmp = moment1d(ivar*inmask_iord == 0.0, sobjs_align[indx][0].TRACE_SPAT, 2*box_rad_pix,
                                row=sobjs_align[indx][0].trace_spec)[0] != pixtot

            flux_box[:,iord,iobj] = flux_tmp*mask_tmp
            ivar_box[:,iord,iobj] = np.fmax(ivar_tmp*mask_tmp,0.0)
            mask_box[:,iord,iobj] = mask_tmp
            mean, med_sn, stddev = stats.sigma_clipped_stats(flux_box[mask_tmp,iord,iobj]*np.sqrt(ivar_box[mask_tmp,iord,iobj]),
                                                             sigma_lower=5.0,sigma_upper=5.0)
            # ToDO assign this to sobjs_align for use in the extraction
            SNR_arr[iord,iobj] = med_sn
            sobjs_align[indx][0].ech_snr = med_sn
            # For hand extractions
            slitfracpos_arr[iord,iobj] = sobjs_align[indx][0].SPAT_FRACPOS

    # Purge objects with low SNR that don't show up in enough orders, sort the list of objects with respect to obj_id
    # and orderindx
    keep_obj = np.zeros(nobj,dtype=bool)
    sobjs_trim = specobjs.SpecObjs()
    # objids are 1 based so that we can easily asign the negative to negative objects
    iobj_keep = 1
    iobj_keep_not_hand = 1

    # TODO JFH: Fix this ugly and dangerous hack that was added to accomodate hand apertures
    hand_frac = [-1000] if hand_extract_dict is None else [int(np.round(ispat*1000)) for ispat in f_spats]

    ## Loop over objects from highest SNR to lowest SNR. Apply the S/N constraints. Once we hit the maximum number
    # objects requested exit, except keep the hand apertures that were requested.
    isort_SNR_max = np.argsort(np.median(SNR_arr,axis=0))[::-1]
    for iobj in isort_SNR_max:
        hand_ap_flag = int(np.round(slitfracpos_arr[0, iobj]*1000)) in hand_frac
        SNR_constraint = (SNR_arr[:,iobj].max() > max_snr) or (np.sum(SNR_arr[:,iobj] > min_snr) >= nabove_min_snr)
        nperorder_constraint = (iobj_keep-1) < nperorder
        if (SNR_constraint and nperorder_constraint) or hand_ap_flag:
            keep_obj[iobj] = True
            ikeep = sobjs_align.ECH_OBJID == uni_obj_id[iobj]
            sobjs_keep = sobjs_align[ikeep].copy()
            sobjs_keep.ECH_OBJID = iobj_keep
            sobjs_keep.OBJID = iobj_keep
#            for spec in sobjs_keep:
#                spec.ECH_OBJID = iobj_keep
#                #spec.OBJID = iobj_keep
            sobjs_trim.add_sobj(sobjs_keep[np.argsort(sobjs_keep.ECH_ORDERINDX)])
            iobj_keep += 1
            if not hand_ap_flag:
                iobj_keep_not_hand += 1
        else:
            if not nperorder_constraint:
                msgs.info('Purging object #{:d}'.format(iobj) +
                          ' since there are already {:d} objects automatically identified '
                          'and you set nperorder={:d}'.format(iobj_keep_not_hand-1, nperorder))
            else:
                msgs.info('Purging object #{:d}'.format(iobj) + ' which does not satisfy max_snr > {:5.2f} OR min_snr > {:5.2f}'.format(max_snr, min_snr) +
                ' on at least nabove_min_snr >= {:d}'.format(nabove_min_snr) + ' orders')


    nobj_trim = np.sum(keep_obj)

    if nobj_trim == 0:
        msgs.warn('No objects found')
        sobjs_final = specobjs.SpecObjs()
        return sobjs_final

    # TODO JFH: We need to think about how to implement returning a maximum number of objects, where the objects
    # returned are the highest S/N ones. It is a bit complicated with regards to the individual object finding and then
    # the linking that is performed above, and also making sure the hand apertures don't get removed.
    SNR_arr_trim = SNR_arr[:,keep_obj]


    sobjs_final = sobjs_trim.copy()
    # Loop over the objects one by one and adjust/predict the traces
    pca_fits = np.zeros((nspec, norders, nobj_trim))

    # Create the trc_inmask for iterative fitting below
    trc_inmask = np.zeros((nspec, norders), dtype=bool)
    for iord in range(norders):
        trc_inmask[:,iord] = (spec_vec >= spec_min_max[0,iord]) & (spec_vec <= spec_min_max[1,iord])

    for iobj in range(nobj_trim):
        indx_obj_id = sobjs_final.ECH_OBJID == (iobj + 1)
        # PCA predict all the orders now (where we have used the standard or slit boundary for the bad orders above)
        msgs.info('Fitting echelle object finding PCA for object {:d}\{:d} with median SNR = {:5.3f}'.format(
            iobj + 1,nobj_trim,np.median(sobjs_final[indx_obj_id].ech_snr)))
        pca_fits[:,:,iobj] \
                = tracepca.pca_trace_object(sobjs_final[indx_obj_id].TRACE_SPAT.T,
                                            order=coeff_npoly, npca=npca,
                                            pca_explained_var=pca_explained_var,
                                            trace_wgt=np.fmax(sobjs_final[indx_obj_id].ech_snr, 1.0)**2,
                                            debug=show_pca)

        # Trial and error shows weighting by S/N instead of S/N^2 performs better
        # JXP -- Updated to now be S/N**2, i.e. inverse variance, with fitting fit

        # Perform iterative flux weighted centroiding using new PCA predictions
        xinit_fweight = pca_fits[:,:,iobj].copy()
        inmask_now = inmask & allmask
        xfit_fweight = fit_trace(image, xinit_fweight, ncoeff, bpm=np.invert(inmask_now),
                                 trace_bpm=np.invert(trc_inmask), fwhm=fwhm, maxdev=maxdev,
                                 debug=show_fits)[0]

        # Perform iterative Gaussian weighted centroiding
        xinit_gweight = xfit_fweight.copy()
        xfit_gweight = fit_trace(image, xinit_gweight, ncoeff, bpm=np.invert(inmask_now),
                                 trace_bpm=np.invert(trc_inmask), weighting='gaussian', fwhm=fwhm,
                                 maxdev=maxdev, debug=show_fits)[0]

        #TODO  Assign the new traces. Only assign the orders that were not orginally detected and traced. If this works
        # well, we will avoid doing all of the iter_tracefits above to make the code faster.
        for iord, spec in enumerate(sobjs_final[indx_obj_id]):
            # JFH added the condition on ech_frac_was_fit with S/N cut on 7-7-19.
            # TODO is this robust against half the order being masked?
            if spec.ech_frac_was_fit & (spec.ech_snr > 1.0):
                    spec.TRACE_SPAT = xfit_gweight[:,iord]
                    spec.SPAT_PIXPOS = spec.TRACE_SPAT[specmid]

    #TODO Put in some criterion here that does not let the fractional position change too much during the iterative
    # tracefitting. The problem is spurious apertures identified on one slit can be pulled over to the center of flux
    # resulting in a bunch of objects landing on top of each other.

    # Set the IDs
    sobjs_final[:].ECH_ORDER = order_vec[sobjs_final[:].ECH_ORDERINDX]
    #for spec in sobjs_final:
    #    spec.ech_order = order_vec[spec.ECH_ORDERINDX]
    sobjs_final.set_names()

    if show_trace:
        viewer, ch = display.show_image(image*allmask)

        for spec in sobjs_trim:
            color = 'red' if spec.ech_frac_was_fit else 'magenta'
            ## Showing the final flux weighted centroiding from PCA predictions
            display.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        for iobj in range(nobj_trim):
            for iord in range(norders):
                ## Showing PCA predicted locations before recomputing flux/gaussian weighted centroiding
                display.show_trace(viewer, ch, pca_fits[:,iord, iobj], str(uni_frac[iobj]), color='yellow')
                ## Showing the final traces from this routine
                display.show_trace(viewer, ch, sobjs_final.TRACE_SPAT[iord].T, sobjs_final.NAME, color='cyan')

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
    # Vette
    for sobj in sobjs_final:
        if not sobj.ready_for_extraction():
            msgs.error("Bad SpecObj.  Can't proceed")

    return sobjs_final


def objs_in_slit(image, thismask, slit_left, slit_righ, inmask=None, fwhm=3.0, use_user_fwhm=False, boxcar_rad=7.,
                 maxdev=2.0, has_negative=False, spec_min_max=None, hand_extract_dict=None, std_trace=None,
                 ncoeff=5, nperslit=None, sig_thresh=10.0, peak_thresh=0.0, abs_thresh=0.0, trim_edg=(5,5),
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

    #debug_all=True
    if debug_all:
        show_peaks=True
        show_fits = True
        show_trace = True
        show_cont = True

    if specobj_dict is None:
        specobj_dict = dict(SLITID=999, DET='DET01', OBJTYPE='unknown', PYPELINE='MultiSlit')

    # Check that peak_thresh values make sense
    if peak_thresh < 0 or peak_thresh > 1:
        msgs.error('Invalid value of peak_thresh. It must be between 0.0 and 1.0')

    nspec, nspat = image.shape
    specmid = nspec//2

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
    thisimg = image*totmask
    #  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    #nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())
    # Mask skypixels with 2 fwhm of edge
    left_asym = slit_left[:,None] + np.outer(xsize/nsamp, np.arange(nsamp))
    righ_asym = left_asym + np.outer(xsize/nsamp, np.ones(int(nsamp)))
    # This extract_asymbox2 call smashes the image in the spectral direction along the curved object traces
    # TODO Should we be passing the mask here with extract_asymbox or not?
    flux_spec = moment1d(thisimg, (left_asym+righ_asym)/2, (righ_asym-left_asym),
                         fwgt=totmask.astype(float))[0]
    mask_spec = moment1d(totmask, (left_asym+righ_asym)/2, (righ_asym-left_asym),
                         fwgt=totmask.astype(float))[0] < 0.3
    if find_min_max is not None:
        find_spec_min,find_spec_max = int(find_min_max[0]), int(find_min_max[1])
        flux_spec = flux_spec[find_spec_min:find_spec_max,:]
        mask_spec = mask_spec[find_spec_min:find_spec_max,:]

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

    gauss_smth_sigma = (fwhm/2.3548) # JFH Reduced by two
    smash_mask = np.isfinite(flux_mean)
    flux_mean_med0 = np.median(flux_mean[smash_mask])
    flux_mean[np.invert(smash_mask)] = flux_mean_med0
    fluxsub0 = flux_mean - flux_mean_med0
    fluxconv0 = scipy.ndimage.filters.gaussian_filter1d(fluxsub0, gauss_smth_sigma, mode='nearest')

    #show_cont=True
    cont_samp = np.fmin(int(np.ceil(nsamp/(fwhm/2.3548))), 30)
    cont, cont_mask0 = arc.iter_continuum(
        fluxconv0, inmask=smash_mask, fwhm=fwhm, cont_frac_fwhm=2.0, sigthresh=cont_sig_thresh, sigrej=2.0, cont_samp=cont_samp,
        npoly=(0 if (nsamp/fwhm < 20.0) else npoly_cont), cont_mask_neg=has_negative, debug=show_cont, debug_peak_find=False,
        qa_title='Smash Image Background, 1st iteration: Slit# {:d}'.format(specobj_dict['SLITID']))

    # Second iteration
    flux_mean_med = np.median(flux_mean[cont_mask0])
    fluxsub = flux_mean - flux_mean_med
    fluxconv = scipy.ndimage.filters.gaussian_filter1d(fluxsub, gauss_smth_sigma, mode='nearest')

    cont, cont_mask = arc.iter_continuum(
        fluxconv, inmask=smash_mask, fwhm=fwhm, cont_frac_fwhm=2.0, sigthresh=cont_sig_thresh, sigrej=2.0, cont_samp=cont_samp,
        npoly=(0 if (nsamp/fwhm < 20.0) else npoly_cont), cont_mask_neg=has_negative, debug=show_cont, debug_peak_find=False,
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
        cont_mask = np.ones(int(nsamp),dtype=bool) # if all pixels are masked for some reason, don't mask

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
    ypeak, _, xcen, sigma_pk, _, good_indx, _, _ = arc.detect_lines(fluxconv_cont, cont_subtract = False, fwhm = fwhm,
                                                                    max_frac_fwhm = 5.0, input_thresh = 'None', debug=False)
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
            thisobj.SPAT_FRACPOS = xcen[iobj]/nsamp
            thisobj.smash_peakflux = ypeak[iobj]
            thisobj.smash_nsig = ypeak[iobj]/sigma
            sobjs.add_sobj(thisobj)
    else:
        nobj_reg = 0

    # ToDo Also plot the edge trimming boundaries on the QA here.
    if show_peaks or objfindQA_filename is not None:
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
        # plt.title(qa_title + ': Slit# {:d}'.format(objfindQA_dict['SLITORD_ID']))
        plt.title(qa_title)
        if objfindQA_filename is not None:
            plt.savefig(objfindQA_filename, dpi=400)
        if show_peaks:
            viewer, ch = display.show_image(image*(thismask*inmask))
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
        else:    # If no standard is provided shift left slit boundary over to be initial trace
            # ToDO make this the average left and right boundary instead. That would be more robust.
            sobjs[iobj].TRACE_SPAT = slit_left + xsize*sobjs[iobj].SPAT_FRACPOS
        sobjs[iobj].trace_spec = spec_vec
        sobjs[iobj].SPAT_PIXPOS = sobjs[iobj].TRACE_SPAT[specmid]
        # Set the idx for any prelminary outputs we print out. These will be updated shortly
        sobjs[iobj].set_name()

        # assign FWHM
        if use_user_fwhm:
            sobjs[iobj].FWHM = fwhm

        else:
            # Determine the fwhm max
            yhalf = 0.5*sobjs[iobj].smash_peakflux
            xpk = sobjs[iobj].SPAT_FRACPOS*nsamp
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

            # Set FWHM for the object
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
        hand_extract_det  = hand_extract_det[hand_on_slit]
        hand_extract_fwhm = hand_extract_fwhm[hand_on_slit]
        nobj_hand = len(hand_extract_spec)
        msgs.info("Implementing hand apertures for {} sources on the slit".format(nobj_hand))

        # Decide how to assign a trace to the hand objects
        if nobj_reg > 0:  # Use brightest object on slit?
            smash_peakflux = sobjs.smash_peakflux
            ibri = smash_peakflux.argmax()
            trace_model = sobjs[ibri].TRACE_SPAT
            med_fwhm_reg = np.median(sobjs.FWHM)
        elif std_trace is not None:   # If no objects found, use the standard?
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
            thisobj.SPAT_FRACPOS = float(f_ximg(thisobj.hand_extract_spec, thisobj.hand_extract_spat, grid=False)) # interpolate from ximg
            thisobj.smash_peakflux = np.interp(thisobj.SPAT_FRACPOS*nsamp,np.arange(nsamp),fluxconv_cont) # interpolate from fluxconv
            # assign the trace
            spat_0 = np.interp(thisobj.hand_extract_spec, spec_vec, trace_model)
            shift = thisobj.hand_extract_spat - spat_0
            thisobj.TRACE_SPAT = trace_model + shift
            thisobj.trace_spec = spec_vec
            thisobj.SPAT_PIXPOS = thisobj.TRACE_SPAT[specmid]
            thisobj.set_name()
            # assign FWHM
            # TODO -- I think FWHM *has* to be input
            if hand_extract_fwhm[iobj] is not None: # If a hand_extract_fwhm was input use that for the fwhm
                thisobj.FWHM = hand_extract_fwhm[iobj]
            elif nobj_reg > 0: # Otherwise is None was input, then use the median of objects on this slit if they are present
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
        #spat_pixpos = np.array([spec.SPAT_PIXPOS for spec in specobjs])
        #hand_flag = np.array([spec.hand_extract_flag for spec in specobjs])
        #spec_fwhm = np.array([spec.FWHM for spec in specobjs])
        reg_ind, = np.where(np.invert(hand_flag))
        hand_ind, = np.where(hand_flag)
        #med_fwhm = np.median(spec_fwhm[~hand_flag])
        #spat_pixpos_hand = spat_pixpos[hand_ind]
        keep = np.ones(nobj, dtype=bool)
        for ihand in hand_ind:
            close = np.abs(sobjs[reg_ind].SPAT_PIXPOS - spat_pixpos[ihand]) <= 0.6*spec_fwhm[ihand]
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
            sobjs[iobj].maskwidth = extract_maskwidth*sobjs[iobj].FWHM*(1.0 + 0.5*np.log10(np.fmax(sobjs[iobj].smash_peakflux/skythresh,1.0)))
        else:
            sobjs[iobj].maskwidth = extract_maskwidth*sobjs[iobj].FWHM

    # If requested display the resulting traces on top of the image
    if show_trace:
        viewer, ch = display.show_image(image*(thismask*inmask))
        display.show_slits(viewer, ch, slit_left.T, slit_righ.T, slit_ids = sobjs[0].SLITID)
        for iobj in range(nobj):
            if sobjs[iobj].hand_extract_flag == False:
                color = 'orange'
            else:
                color = 'blue'
            display.show_trace(viewer, ch,sobjs[iobj].TRACE_SPAT, trc_name = sobjs[iobj].NAME, color=color)

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

