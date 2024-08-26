
def orig_ech_objfind(image, ivar, slitmask, slit_left, slit_righ, order_vec, maskslits, det='DET01',
                     inmask=None, spec_min_max=None, fof_link=1.5, plate_scale=0.2,
                     std_trace=None, ncoeff=5, npca=None, coeff_npoly=None, max_snr=2.0, min_snr=1.0,
                     nabove_min_snr=2, pca_explained_var=99.0, box_radius=2.0, fwhm=3.0,
                     use_user_fwhm=False, maxdev=2.0, hand_extract_dict=None, nperorder=2,
                     extract_maskwidth=3.0, snr_thresh=10.0,
                     specobj_dict=None, trim_edg=(5, 5),
                     show_peaks=False, show_fits=False, show_single_fits=False,
                     show_trace=False, show_single_trace=False, show_pca=False,
                     debug_all=False, objfindQA_filename=None):
    """
    Object finding routine for Echelle spectrographs.

    This routine:

        #. Runs object finding on each order individually

        #. Links the objects found together using a friends-of-friends algorithm
           on fractional order position.

        #. For objects which were only found on some orders, the standard (or
           the slit boundaries) are placed at the appropriate fractional
           position along the order.

        #. A PCA fit to the traces is performed using the routine above pca_fit

    Args:
        image (`numpy.ndarray`_):
            (Floating-point) Image to use for object search with shape (nspec,
            nspat).  The first dimension (nspec) is spectral, and second
            dimension (nspat) is spatial. Note this image can either have the
            sky background in it, or have already been sky subtracted.  Object
            finding works best on sky-subtracted images. Ideally, object finding
            is run in another routine, global sky-subtraction performed, and
            then this code should be run. However, it is also possible to run
            this code on non-sky-subtracted images.
        ivar (`numpy.ndarray`_):
            Floating-point inverse variance image for the input image.  Shape
            must match ``image``, (nspec, nspat).
        slitmask (`numpy.ndarray`_):
            Integer image indicating the pixels that belong to each order.
            Pixels that are not on an order have value -1, and those that are on
            an order have a value equal to the slit number (i.e. 0 to nslits-1
            from left to right on the image).  Shape must match ``image``,
            (nspec, nspat).
        slit_left (`numpy.ndarray`_):
            Left boundary of orders to be extracted (given as floating point
            pixels).  Shape is (nspec, norders), where norders is the total
            number of traced echelle orders.
        slit_righ (`numpy.ndarray`_):
            Right boundary of orders to be extracted (given as floating point
            pixels).  Shape is (nspec, norders), where norders is the total
            number of traced echelle orders.
        order_vec (`numpy.ndarray`_):
            Vector identifying the Echelle orders for each pair of order edges
            found.  This is saved to the output :class:`~pypeit.specobj.SpecObj`
            objects.  If the orders are not known, this can be
            ``np.arange(norders)`` (but this is *not* recommended).
        maskslits (`numpy.ndarray`_):
            Boolean array selecting orders that should be ignored (i.e., good
            orders are False, bad orders are True).  Shape must be (norders,).
        det (:obj:`str`, optional):
            The name of the detector containing the object.  Only used if
            ``specobj_dict`` is None.
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Must have the same shape as
            ``image``.  If None, all pixels in ``slitmask`` with non-negative
            values are considered good.
        spec_min_max (`numpy.ndarray`_, optional):
            2D array defining the minimum and maximum pixel in the spectral
            direction with useable data for each order.  Shape must be (2,
            norders).  This should only be used for echelle spectrographs for
            which the orders do not entirely cover the detector. PCA tracing
            will re-map the traces such that they all have the same length,
            compute the PCA, and then re-map the orders back.  This improves
            performance for echelle spectrographs by removing the nonlinear
            shrinking of the orders so that the linear pca operation can better
            predict the traces. If None, the minimum and maximum values will be
            determined automatically from ``slitmask``.
        fof_link (:obj:`float`, optional):
            Friends-of-friends linking length in arcseconds used to link
            together traces across orders. The routine links together at
            the same fractional slit position and links them together
            with a friends-of-friends algorithm using this linking
            length.
        plate_scale (:obj:`float`, `numpy.ndarray`_, optional):
            Plate scale in arcsec/pix. This can either be a single float for
            every order, or an array with shape (norders,) providing the plate
            scale of each order.
        std_trace (`numpy.ndarray`_, optional):
            Vector with the standard star trace, which is used as a crutch for
            tracing.  Shape must be (nspec,).  If None, the slit boundaries are
            used as the crutch.
        ncoeff (:obj:`int`, optional):
            Order of polynomial fit to traces.
        npca (:obj:`int`, optional):
            Number of PCA components to keep during PCA decomposition of the
            object traces.  If None, the number of components set by requiring
            the PCA accounts for approximately 99% of the variance.
        coeff_npoly (:obj:`int`, optional):
            Order of polynomial used for PCA coefficients fitting.  If None,
            value set automatically, see
            :func:`~pypeit.tracepca.pca_trace_object`.
        max_snr (:obj:`float`, optional):
            For an object to be included in the output object, it must have a
            max S/N ratio above this value.
        min_snr (:obj:`float`, optional):
            For an object to be included in the output object, it must have a
            a median S/N ratio above this value for at least
            ``nabove_min_snr`` orders (see below).
        nabove_min_snr (:obj:`int`, optional):
            The required number of orders that an object must have with median
            SNR greater than ``min_snr`` in order to be included in the output
            object.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance in the data
            accounted for by the PCA used to truncate the number of PCA
            coefficients to keep (see ``npca``). Ignored if ``npca`` is provided
            directly; see :func:`~pypeit.tracepca.pca_trace_object`.
        box_radius (:obj:`float`, optional):
            Box_car extraction radius in arcseconds to assign to each detected
            object and to be used later for boxcar extraction. In this method
            ``box_radius`` is converted into pixels using ``plate_scale``.
            ``box_radius`` is also used for SNR calculation and trimming.
        fwhm (:obj:`float`, optional):
            Estimated fwhm of the objects in pixels
        use_user_fwhm (:obj:`bool`, optional):
            If True, ``PypeIt`` will use the spatial profile FWHM input by the
            user (see ``fwhm``) rather than determine the spatial FWHM from the
            smashed spatial profile via the automated algorithm.
        maxdev (:obj:`float`, optional):
            Maximum deviation of pixels from polynomial fit to trace
            used to reject bad pixels in trace fitting.
        hand_extract_dict (:obj:`dict`, optional):
            Dictionary with info on manual extraction; see
            :class:`~pypeit.manual_extract.ManualExtractionObj`.
        nperorder (:obj:`int`, optional):
            Maximum number of objects allowed per order.  If there are more
            detections than this number, the code will select the ``nperorder``
            most significant detections. However, hand apertures will always be
            returned and do not count against this budget.
        extract_maskwidth (:obj:`float`, optional):
            Determines the initial size of the region in units of FWHM that will
            be used for local sky subtraction; See :func:`objs_in_slit` and
            :func:`~pypeit.core.skysub.local_skysub_extract`.
        snr_thresh (:obj:`float`, optional):
            SNR threshold for finding objects
        specobj_dict (:obj:`dict`, optional):
            Dictionary containing meta-data for the objects that will be
            propagated into the :class:`~pypeit.specobj.SpecObj` objects.  The
            expected components are:

                - SLITID: The slit ID number
                - DET: The detector identifier
                - OBJTYPE: The object type
                - PYPELINE: The class of pipeline algorithms applied

            If None, the dictionary is filled with the following placeholders::

                specobj_dict = {'SLITID': 999, 'DET': 'DET01',
                                'OBJTYPE': 'unknown', 'PYPELINE': 'unknown'}

        trim_edg (:obj:`tuple`, optional):
            A two-tuple of integers or floats used to ignore objects within this
            many pixels of the left and right slit boundaries, respectively.
        show_peaks (:obj:`bool`, optional):
            Plot the QA of the object peak finding in each order.
        show_fits (:obj:`bool`, optional):
            Plot trace fitting for final fits using PCA as crutch.
        show_single_fits (:obj:`bool`, optional):
            Plot trace fitting for single order fits.
        show_trace (:obj:`bool`, optional):
            Display the object traces on top of the image.
        show_single_trace (:obj:`bool`, optional):
            Display the object traces on top of the single order.
        show_pca (:obj:`bool`, optional):
            Display debugging plots for the PCA decomposition.
        debug_all (:obj:`bool`, optional):
            Show all the debugging plots.  If True, this also overrides any
            provided values for ``show_peaks``, ``show_trace``, and
            ``show_pca``, setting them to True.
        objfindQA_filename (:obj:`str`, optional):
            Full path (directory and filename) for the object profile QA plot.
            If None, not plot is produced and saved.

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`: Object containing the objects
        detected.
    """
    raise DeprecationWarning
    msgs.error("This ginormous method as been Deprecated")

    # debug_all=True
    if debug_all:
        show_peaks = True
        # show_fits = True
        # show_single_fits = True
        show_trace = True
        show_pca = True
        # show_single_trace = True
        # TODO: This isn't used, right?
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
    if gdslit_spat.size != np.sum(np.invert(maskslits)):
        msgs.error('Masking of slitmask not in sync with that of maskslits.  This is a bug')
        # msgs.error('There is a mismatch between the number of valid orders found by PypeIt and '
        #           'the number expected for this spectrograph.  Unable to continue.  Please '
        #           'submit an issue on Github: https://github.com/pypeit/PypeIt/issues .')

    if spec_min_max is None:
        spec_min_max = np.zeros((2, norders), dtype=int)
        for iord in range(norders):
            ispec, ispat = np.where(slitmask == gdslit_spat[iord])
            spec_min_max[:, iord] = ispec.min(), ispec.max()

    # Setup the plate scale
    if isinstance(plate_scale, (float, int)):
        plate_scale_ord = np.full(norders, plate_scale)
    elif isinstance(plate_scale, (np.ndarray, list, tuple)):
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
    slit_spec_pos = nspec / 2.0

    # TODO JFH This hand apertures in echelle needs to be completely refactored.
    # Hand prep
    #   Determine the location of the source on *all* of the orders
    if hand_extract_dict is not None:
        f_spats = []
        for ss, spat, spec in zip(range(len(hand_extract_dict['spec'])),
                                  hand_extract_dict['spat'],
                                  hand_extract_dict['spec']):
            # Find the input slit
            ispec = int(np.clip(np.round(spec), 0, nspec - 1))
            ispat = int(np.clip(np.round(spat), 0, nspat - 1))
            slit = slitmask[ispec, ispat]
            if slit == -1:
                msgs.error('You are requesting a manual extraction at a position ' +
                           f'(spat, spec)={spat, spec} that is not on one of the echelle orders. Check your pypeit file.')
            # Fractions
            iord_hand = gdslit_spat.tolist().index(slit)
            f_spat = (spat - slit_left[ispec, iord_hand]) / (
                    slit_righ[ispec, iord_hand] - slit_left[ispec, iord_hand])
            f_spats.append(f_spat)

    # Loop over orders and find objects
    sobjs = specobjs.SpecObjs()
    # TODO: replace orderindx with the true order number here? Maybe not. Clean
    # up SLITID and orderindx!
    gdorders = np.arange(norders)[np.invert(maskslits)]
    for iord in gdorders:  # range(norders):
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
                new_hand_extract_dict['spat'][ss] = slit_left[ispec, iord] + f_spat * (
                        slit_righ[ispec, iord] - slit_left[ispec, iord])
        else:
            new_hand_extract_dict = None

        # Get SLTIORD_ID for the objfind QA
        ech_objfindQA_filename = objfindQA_filename.replace('S0999', 'S{:04d}'.format(order_vec[iord])) \
            if objfindQA_filename is not None else None
        # Run
        sobjs_slit = \
            objs_in_slit(image, ivar, thisslit_gpm, slit_left[:, iord], slit_righ[:, iord],
                         spec_min_max=spec_min_max[:, iord],
                         inmask=inmask_iord, std_trace=std_in, ncoeff=ncoeff, fwhm=fwhm, use_user_fwhm=use_user_fwhm,
                         maxdev=maxdev,
                         hand_extract_dict=new_hand_extract_dict, nperslit=nperorder,
                         extract_maskwidth=extract_maskwidth,
                         snr_thresh=snr_thresh, trim_edg=trim_edg, boxcar_rad=box_radius / plate_scale_ord[iord],
                         show_peaks=show_peaks, show_fits=show_single_fits,
                         show_trace=show_single_trace, qa_title=qa_title, specobj_dict=specobj_dict,
                         objfindQA_filename=ech_objfindQA_filename)
        sobjs.add_sobj(sobjs_slit)

    nfound = len(sobjs)

    if nfound == 0:
        msgs.warn('No objects found')
        return sobjs

    FOF_frac = fof_link / (np.median(np.median(slit_width, axis=0) * plate_scale_ord))
    # Run the FOF. We use fake coordinates
    fracpos = sobjs.SPAT_FRACPOS
    ra_fake = fracpos / 1000.0  # Divide all angles by 1000 to make geometry euclidian
    dec_fake = np.zeros_like(fracpos)
    if nfound > 1:
        inobj_id, multobj_id, firstobj_id, nextobj_id \
            = pydl.spheregroup(ra_fake, dec_fake, FOF_frac / 1000.0)
        # TODO spheregroup returns zero based indices but we use one based. We should probably add 1 to inobj_id here,
        # i.e. obj_id_init = inobj_id + 1
        obj_id_init = inobj_id.copy()
    elif nfound == 1:
        obj_id_init = np.zeros(1, dtype='int')

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
                ind_rest = np.setdiff1d(ind, ind[min_dist_ind])
                # JFH OLD LINE with bug
                # obj_id[ind_rest] = (np.arange(len(ind_rest)) + 1) + obj_id_init.max()
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
    isort_frac = uni_frac.argsort(kind='stable')
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
        if (nthisobj_id > 3) and (nthisobj_id < norders):
            thisorderindx = sobjs_align[indx_obj_id].ECH_ORDERINDX
            goodorder = np.zeros(norders, dtype=bool)
            goodorder[thisorderindx] = True
            badorder = np.invert(goodorder)
            xcen_good = (sobjs_align[indx_obj_id].TRACE_SPAT).T
            slit_frac_good = (xcen_good - slit_left[:, goodorder]) / slit_width[:, goodorder]
            # Fractional slit position averaged across the spectral direction for each order
            frac_mean_good = np.mean(slit_frac_good, 0)
            # Perform a  linear fit to fractional slit position
            # TODO Do this as a S/N weighted fit similar to what is now in the pca_trace algorithm?
            # msk_frac, poly_coeff_frac = fitting.robust_fit(order_vec[goodorder], frac_mean_good, 1,
            pypeitFit = fitting.robust_fit(order_vec[goodorder], frac_mean_good, 1,
                                           function='polynomial', maxiter=20, lower=2, upper=2,
                                           use_mad=True, sticky=False,
                                           minx=order_vec.min(), maxx=order_vec.max())
            frac_mean_new = np.zeros(norders)
            frac_mean_new[badorder] = pypeitFit.eval(
                order_vec[badorder])  # , minx = order_vec.min(),maxx=order_vec.max())
            frac_mean_new[goodorder] = frac_mean_good
            # TODO This QA needs some work
            if show_pca:
                frac_mean_fit = pypeitFit.eval(order_vec)
                plt.plot(order_vec[goodorder][pypeitFit.bool_gpm], frac_mean_new[goodorder][pypeitFit.bool_gpm], 'ko',
                         mfc='k', markersize=8.0, label='Good Orders Kept')
                plt.plot(order_vec[goodorder][np.invert(pypeitFit.bool_gpm)],
                         frac_mean_new[goodorder][np.invert(pypeitFit.bool_gpm)], 'ro', mfc='k', markersize=8.0,
                         label='Good Orders Rejected')
                plt.plot(order_vec[badorder], frac_mean_new[badorder], 'ko', mfc='None', markersize=8.0,
                         label='Predicted Bad Orders')
                plt.plot(order_vec, frac_mean_new, '+', color='cyan', markersize=12.0, label='Final Order Fraction')
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
                # thisobj.ECH_ORDERINDX = iord
                # thisobj.ech_order = order_vec[iord]
                thisobj.SPAT_FRACPOS = uni_frac[iobj]
                # Assign traces using the fractional position fit above
                if std_trace is not None:
                    x_trace = np.interp(slit_spec_pos, spec_vec, std_trace[:, iord])
                    shift = np.interp(slit_spec_pos, spec_vec,
                                      slit_left[:, iord] + slit_width[:, iord] * frac_mean_new[iord]) - x_trace
                    thisobj.TRACE_SPAT = std_trace[:, iord] + shift
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
                thisobj.smash_snr = sobjs_align[imin].smash_snr
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
                msgs.error(
                    'Problem in echelle object finding. The same objid={:d} appears {:d} times on echelle orderindx ={:d}'
                    ' even after duplicate obj_ids the orders were removed. '
                    'Report this bug to PypeIt developers'.format(uni_obj_id[iobj], num_on_order, iord))

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
            # indx = (sobjs_align.ECH_OBJID == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDERINDX == iord)
            # spec = sobjs_align[indx][0]
            inmask_iord = inmask & (slitmask == gdslit_spat[iord])
            # TODO make the snippet below its own function quick_extraction()
            box_rad_pix = box_radius / plate_scale_ord[iord]

            # TODO -- We probably shouldn't be operating on a SpecObjs but instead a SpecObj
            flux_tmp = moment1d(image * inmask_iord, sobjs_align[indx][0].TRACE_SPAT, 2 * box_rad_pix,
                                row=sobjs_align[indx][0].trace_spec)[0]
            var_tmp = moment1d(varimg * inmask_iord, sobjs_align[indx][0].TRACE_SPAT, 2 * box_rad_pix,
                               row=sobjs_align[indx][0].trace_spec)[0]
            ivar_tmp = utils.calc_ivar(var_tmp)
            pixtot = moment1d(ivar * 0 + 1.0, sobjs_align[indx][0].TRACE_SPAT, 2 * box_rad_pix,
                              row=sobjs_align[indx][0].trace_spec)[0]
            mask_tmp = moment1d(ivar * inmask_iord == 0.0, sobjs_align[indx][0].TRACE_SPAT, 2 * box_rad_pix,
                                row=sobjs_align[indx][0].trace_spec)[0] != pixtot

            flux_box[:, iord, iobj] = flux_tmp * mask_tmp
            ivar_box[:, iord, iobj] = np.fmax(ivar_tmp * mask_tmp, 0.0)
            mask_box[:, iord, iobj] = mask_tmp
            mean, med_sn, stddev = astropy.stats.sigma_clipped_stats(
                flux_box[mask_tmp, iord, iobj] * np.sqrt(ivar_box[mask_tmp, iord, iobj]),
                sigma_lower=5.0, sigma_upper=5.0
            )
            # ToDO assign this to sobjs_align for use in the extraction
            SNR_arr[iord, iobj] = med_sn
            sobjs_align[indx][0].ech_snr = med_sn
            # For hand extractions
            slitfracpos_arr[iord, iobj] = sobjs_align[indx][0].SPAT_FRACPOS

    # Purge objects with low SNR that don't show up in enough orders, sort the list of objects with respect to obj_id
    # and orderindx
    keep_obj = np.zeros(nobj, dtype=bool)
    sobjs_trim = specobjs.SpecObjs()
    # objids are 1 based so that we can easily asign the negative to negative objects
    iobj_keep = 1
    iobj_keep_not_hand = 1

    # TODO JFH: Fix this ugly and dangerous hack that was added to accomodate hand apertures
    hand_frac = [-1000] if hand_extract_dict is None else [int(np.round(ispat * 1000)) for ispat in f_spats]

    ## Loop over objects from highest SNR to lowest SNR. Apply the S/N constraints. Once we hit the maximum number
    # objects requested exit, except keep the hand apertures that were requested.
    isort_SNR_max = np.argsort(np.median(SNR_arr, axis=0), kind='stable')[::-1]
    for iobj in isort_SNR_max:
        hand_ap_flag = int(np.round(slitfracpos_arr[0, iobj] * 1000)) in hand_frac
        SNR_constraint = (SNR_arr[:, iobj].max() > max_snr) or (np.sum(SNR_arr[:, iobj] > min_snr) >= nabove_min_snr)
        nperorder_constraint = (iobj_keep - 1) < nperorder
        if (SNR_constraint and nperorder_constraint) or hand_ap_flag:
            keep_obj[iobj] = True
            ikeep = sobjs_align.ECH_OBJID == uni_obj_id[iobj]
            sobjs_keep = sobjs_align[ikeep].copy()
            sobjs_keep.ECH_OBJID = iobj_keep
            sobjs_keep.OBJID = iobj_keep
            #            for spec in sobjs_keep:
            #                spec.ECH_OBJID = iobj_keep
            #                #spec.OBJID = iobj_keep
            sobjs_trim.add_sobj(sobjs_keep[np.argsort(sobjs_keep.ECH_ORDERINDX, kind='stable')])
            iobj_keep += 1
            if not hand_ap_flag:
                iobj_keep_not_hand += 1
        else:
            if not nperorder_constraint:
                msgs.info('Purging object #{:d}'.format(iobj) +
                          ' since there are already {:d} objects automatically identified '
                          'and you set nperorder={:d}'.format(iobj_keep_not_hand - 1, nperorder))
            else:
                msgs.info('Purging object #{:d}'.format(
                    iobj) + ' which does not satisfy max_snr > {:5.2f} OR min_snr > {:5.2f}'.format(max_snr, min_snr) +
                          ' on at least nabove_min_snr >= {:d}'.format(nabove_min_snr) + ' orders')

    nobj_trim = np.sum(keep_obj)

    if nobj_trim == 0:
        msgs.warn('No objects found')
        sobjs_final = specobjs.SpecObjs()
        return sobjs_final

    # TODO JFH: We need to think about how to implement returning a maximum number of objects, where the objects
    # returned are the highest S/N ones. It is a bit complicated with regards to the individual object finding and then
    # the linking that is performed above, and also making sure the hand apertures don't get removed.
    SNR_arr_trim = SNR_arr[:, keep_obj]

    sobjs_final = sobjs_trim.copy()
    # Loop over the objects one by one and adjust/predict the traces
    pca_fits = np.zeros((nspec, norders, nobj_trim))

    # Create the trc_inmask for iterative fitting below
    trc_inmask = np.zeros((nspec, norders), dtype=bool)
    for iord in range(norders):
        trc_inmask[:, iord] = (spec_vec >= spec_min_max[0, iord]) & (spec_vec <= spec_min_max[1, iord])

    for iobj in range(nobj_trim):
        indx_obj_id = sobjs_final.ECH_OBJID == (iobj + 1)
        # PCA predict all the orders now (where we have used the standard or slit boundary for the bad orders above)
        msgs.info('Fitting echelle object finding PCA for object {:d}/{:d} with median SNR = {:5.3f}'.format(
            iobj + 1, nobj_trim, np.median(sobjs_final[indx_obj_id].ech_snr)))
        pca_fits[:, :, iobj] \
            = tracepca.pca_trace_object(sobjs_final[indx_obj_id].TRACE_SPAT.T,
                                        order=coeff_npoly, npca=npca,
                                        pca_explained_var=pca_explained_var,
                                        trace_wgt=np.fmax(sobjs_final[indx_obj_id].ech_snr, 1.0) ** 2,
                                        debug=show_pca)

        # Trial and error shows weighting by S/N instead of S/N^2 performs better
        # JXP -- Updated to now be S/N**2, i.e. inverse variance, with fitting fit

        # Perform iterative flux weighted centroiding using new PCA predictions
        xinit_fweight = pca_fits[:, :, iobj].copy()
        inmask_now = inmask & allmask
        xfit_fweight = fit_trace(image, xinit_fweight, ncoeff, bpm=np.invert(inmask_now),
                                 trace_bpm=np.invert(trc_inmask), fwhm=fwhm, maxdev=maxdev,
                                 debug=show_fits)[0]

        # Perform iterative Gaussian weighted centroiding
        xinit_gweight = xfit_fweight.copy()
        xfit_gweight = fit_trace(image, xinit_gweight, ncoeff, bpm=np.invert(inmask_now),
                                 trace_bpm=np.invert(trc_inmask), weighting='gaussian', fwhm=fwhm,
                                 maxdev=maxdev, debug=show_fits)[0]

        # TODO  Assign the new traces. Only assign the orders that were not orginally detected and traced. If this works
        # well, we will avoid doing all of the iter_tracefits above to make the code faster.
        for iord, spec in enumerate(sobjs_final[indx_obj_id]):
            # JFH added the condition on ech_frac_was_fit with S/N cut on 7-7-19.
            # TODO is this robust against half the order being masked?
            if spec.ech_frac_was_fit & (spec.ech_snr > 1.0):
                spec.TRACE_SPAT = xfit_gweight[:, iord]
                spec.SPAT_PIXPOS = spec.TRACE_SPAT[specmid]

    # TODO Put in some criterion here that does not let the fractional position change too much during the iterative
    # tracefitting. The problem is spurious apertures identified on one slit can be pulled over to the center of flux
    # resulting in a bunch of objects landing on top of each other.

    # Set the IDs
    sobjs_final[:].ECH_ORDER = order_vec[sobjs_final[:].ECH_ORDERINDX]
    # for spec in sobjs_final:
    #    spec.ech_order = order_vec[spec.ECH_ORDERINDX]
    sobjs_final.set_names()

    if show_trace:
        viewer, ch = display.show_image(image * allmask)

        for spec in sobjs_trim:
            color = 'red' if spec.ech_frac_was_fit else 'magenta'
            ## Showing the final flux weighted centroiding from PCA predictions
            display.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        for iobj in range(nobj_trim):
            for iord in range(norders):
                ## Showing PCA predicted locations before recomputing flux/gaussian weighted centroiding
                display.show_trace(viewer, ch, pca_fits[:, iord, iobj], str(uni_frac[iobj]), color='yellow')
                ## Showing the final traces from this routine
                display.show_trace(viewer, ch, sobjs_final.TRACE_SPAT[iord].T, sobjs_final.NAME, color='cyan')

        # Labels for the points
        text_final = [dict(type='text', args=(nspat / 2 - 40, nspec / 2, 'final trace'),
                           kwargs=dict(color='cyan', fontsize=20))]

        text_pca = [dict(type='text', args=(nspat / 2 - 40, nspec / 2 - 30, 'PCA fit'),
                         kwargs=dict(color='yellow', fontsize=20))]

        text_fit = [dict(type='text', args=(nspat / 2 - 40, nspec / 2 - 60, 'predicted'),
                         kwargs=dict(color='red', fontsize=20))]

        text_notfit = [dict(type='text', args=(nspat / 2 - 40, nspec / 2 - 90, 'originally found'),
                            kwargs=dict(color='magenta', fontsize=20))]

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