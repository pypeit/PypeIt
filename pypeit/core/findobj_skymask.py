""" Methods to find objects

.. include:: ../include/links.rst

"""
import copy

import numpy as np
import scipy.interpolate
import scipy.ndimage
import matplotlib.pyplot as plt

import astropy.stats

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
from pypeit.core import extract
from pypeit.utils import fast_running_median
from IPython import embed


def create_skymask(sobjs, thismask, slit_left, slit_righ, box_rad_pix=None, trim_edg=(5,5),
                   skymask_snr_thresh=1.0):
    r"""
    Creates a sky mask from a :class:`~pypeit.specobjs.SpecObjs` using the fwhm
    of each object and/or the boxcar radius.

    Args:
        sobjs (:class:`~pypeit.specobjs.SpecObjs`):
            Objects for which you would like to create the mask.
        thismask (`numpy.ndarray`_):
            Boolean image selecting pixels that are on the slit.  Shape is
            :math:`(N_{\rm spec}, N_{\rm spat})`.
        slit_left (`numpy.ndarray`_):
            Left boundary of slit/order to be extracted (given as floating point
            pixels). This a 1-d array with shape :math:`(N_{\rm spec}, 1)` or
            :math:`(N_{\rm spec},)`
        slit_righ (`numpy.ndarray`_):
            Right boundary of slit/order to be extracted (given as floating
            point pixels). This a 1-d array with shape :math:`(N_{\rm spec}, 1)`
            or :math:`(N_{\rm spec},)`
        box_rad_pix (:obj:`float`, optional):
            If set, the sky mask will be as wide as this radius in pixels.
        trim_edg (:obj:`tuple`, optional):
            A two-tuple of integers or floats used to ignore objects within this
            many pixels of the left and right slit boundaries, respectively.
        skymask_snr_thresh (:obj:`float`, optional):
            The multiple of the final object finding SNR threshold used to
            create the sky mask via a Gaussian model of the SNR profile.  The
            (spatial) Gaussian model is determined from the image after
            collapsing along the spectral dimension.

    Returns:
        `numpy.ndarray`_: Boolean image with shape :math:`(N_{\rm spec}, N_{\rm
        spat})` (same as ``thismask``) indicating which pixels are usable for
        global sky subtraction (True means the pixel is usable for sky
        subtraction, False means it should be masked when subtracting sky).
    """
    nobj = len(sobjs)
    ximg, _ = pixels.ximg_and_edgemask(slit_left, slit_righ, thismask, trim_edg=trim_edg)
    # How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    #nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())

    # Objmask
    skymask_objsnr = np.copy(thismask)
    if nobj == 0:
        msgs.info('No objects were detected. The entire slit will be used to determine the sky subtraction.')
    else:
        # Compute some inputs for the object mask
        xtmp = (np.arange(nsamp) + 0.5)/nsamp
        # threshold for object finding
        for iobj in range(nobj):
            # this will skip also sobjs with THRESHOLD=0, because are the same that have smash_snr=0.
            if (sobjs[iobj].smash_snr != 0.) and (sobjs[iobj].smash_snr != None):
                qobj = np.zeros_like(xtmp)
                sep = np.abs(xtmp-sobjs[iobj].SPAT_FRACPOS)
                sep_inc = sobjs[iobj].maskwidth/nsamp
                close = sep <= sep_inc
                # This is an analytical SNR profile with a Gaussian shape.
                # JFH modified to use SNR here instead of smash peakflux. I believe that the 2.77 is supposed to be
                # 2.355**2/2, i.e. the argument of a gaussian with sigma = FWHM/2.35
                qobj[close] = sobjs[iobj].smash_snr * \
                               np.exp(np.fmax(-2.77*(sep[close]*nsamp)**2/sobjs[iobj].FWHM**2, -9.0))
                skymask_objsnr[thismask] &= np.interp(ximg[thismask], xtmp, qobj) < skymask_snr_thresh
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
    #     skymask = skymask_objsnr | skymask_fwhm
    # else:  # Enforces boxcar radius masking
    #     skymask = skymask_objsnr & skymask_fwhm
    # DP: I think skymask should always be skymask_objsnr & skymask_fwhm (i.e., not only when box_rad_pix is not None).
    # In the case of skymask_objsnr | skymask_fwhm, if skymask_objsnr cannot be computed, the entire slit
    # is used for sky calculation (i.e., skymask_fwhm will not have effect).

    # DP's change which I don't think we should adopt at this time.
    #skymask = skymask_objsnr & skymask_fwhm

    # JFH restored old behavior after seeing spurious results for X-shooter. I think the issue here is that the fwhm
    # computation from objs_in_slit is not necessarily that reliable and when large amounts of masking are performed
    # on narrow slits/orders, we have problems. We should revisit this after object finding is refactored since
    # maybe then the fwhm estimates will be more robust.
    if box_rad_pix is None and np.all([sobj.smash_snr is not None for sobj in sobjs]) \
            and np.all([sobj.smash_snr != 0. for sobj in sobjs]) and not np.all(skymask_objsnr == thismask):
        # TODO This is a kludge until we refactor this routine. Basically mask design objects that are not auto-ID
        # always have smash_snr undefined. If there is a hybrid situation of auto-ID and maskdesign, the logic
        # here does not really make sense. Soution would be to compute thershold and smash_snr for all objects.
        skymask = skymask_objsnr | skymask_fwhm
    else:  # Enforces boxcar radius masking
        skymask = skymask_objsnr & skymask_fwhm

    # Return
    return skymask[thismask]


def ech_findobj_ineach_order(
    image, ivar, slitmask, slit_left, slit_righ, slit_spats,
    order_vec, orders_gpm, spec_min_max, plate_scale_ord,
    det='DET01', inmask=None, std_trace=None, ncoeff=5, 
    hand_extract_dict=None,
    box_radius=2.0, fwhm=3.0,
    use_user_fwhm=False, maxdev=2.0, nperorder=2,
    extract_maskwidth=3.0, snr_thresh=10.0,
    specobj_dict=None, trim_edg=(5,5),
    show_peaks=False, show_single_fits=False,
    show_single_trace=False, objfindQA_filename=None):
    """
    Find objects in each echelle order, individually.

    This routine:

        - Loops over the good orders

        - Calls the :func:`objs_in_slit` method to find objects in the order.
          See that method for further details.

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
        slit_spats (`numpy.ndarray`_):
            slit_spat values (spatial position 1/2 way up the detector)
            for the orders 
        order_vec (`numpy.ndarray`_):
            Vector identifying the Echelle orders for each pair of order edges
            found.  This is saved to the output :class:`~pypeit.specobj.SpecObj`
            objects.  If the orders are not known, this can be 
            ``np.arange(norders)`` (but this is *not* recommended).
        order_gpm (`numpy.ndarray`_):
            Boolean array indicating which orders are good (True),
            i.e. have good calibrations (wavelengths, etc.).  Shape is
        spec_min_max (`numpy.ndarray`_):
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
        plate_scale_ord (`numpy.ndarray`_):
            An array with shape (norders,) providing the plate 
            scale of each order in arcsec/pix, 
        det (:obj:`str`, optional):
            The name of the detector containing the object.  Only used if
            ``specobj_dict`` is None.
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Must have the same shape as
            ``image``.  If None, all pixels in ``slitmask`` with non-negative
            values are considered good.
        std_trace (`numpy.ndarray`_, optional):
            Vector with the standard star trace, which is used as a crutch for
            tracing.  Shape must be (nspec,).  If None, the slit boundaries are
            used as the crutch.
        ncoeff (:obj:`int`, optional):
            Order of polynomial fit to traces.
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
        nperorder (:obj:`int`, optional):
            Maximum number of objects allowed per order.  If there are more
            detections than this number, the code will select the ``nperorder``
            most significant detections. However, hand apertures will always be
            returned and do not count against this budget.
        specobj_dict (:obj:`dict`, optional):
            Dictionary containing meta-data for the objects that will be
            propagated into the :class:`~pypeit.specobj.SpecObj` objects.  The
            expected components are:
            
                - SLITID: The slit ID number
                - DET: The detector identifier
                - OBJTYPE: The object type
                - PYPELINE: The class of pipeline algorithms applied

            If None, the dictionary is filled with the following placeholders::

                specobj_dict = {'SLITID': 999, 'DET': det, 'ECH_ORDERINDX': 999,
                                'OBJTYPE': 'unknown', 'PYPELINE': 'Echelle'}

        trim_edg (:obj:`tuple`, optional):
            A two-tuple of integers or floats used to ignore objects within this
            many pixels of the left and right slit boundaries, respectively.
        show_peaks (:obj:`bool`, optional):
            Plot the QA of the object peak finding in each order.
        show_single_fits (:obj:`bool`, optional):
            Plot trace fitting for single order fits.
        show_single_trace (:obj:`bool`, optional):
            Display the object traces on top of the single order.
        objfindQA_filename (:obj:`str`, optional):
            Full path (directory and filename) for the object profile QA plot.
            If None, not plot is produced and saved.
     
    Returns:
        :class:`~pypeit.specobjs.SpecObjs`: Object containing the objects
        detected.
    """
    if specobj_dict is None:
        specobj_dict = {'SLITID': 999, 'ECH_ORDERINDX': 999,
                        'DET': det, 'OBJTYPE': 'unknown', 'PYPELINE': 'Echelle'}

    allmask = slitmask > -1
    if inmask is None:
        inmask = allmask

    # Loop over orders and find objects
    sobjs = specobjs.SpecObjs()
    for iord, iorder in enumerate(order_vec):
        if not orders_gpm[iord]:
            continue
        #
        qa_title = 'Finding objects on order # {:d}'.format(iorder)
        msgs.info(qa_title)
        thisslit_gpm = slitmask == slit_spats[iord]
        inmask_iord = inmask & thisslit_gpm
        specobj_dict['SLITID'] = slit_spats[iord]
        specobj_dict['ECH_ORDERINDX'] = iord
        specobj_dict['ECH_ORDER'] = iorder
        std_in = None if std_trace is None else std_trace[:, iord]

        # Get SLTIORD_ID for the objfind QA
        ech_objfindQA_filename = objfindQA_filename.replace('S0999', 'S{:04d}'.format(order_vec[iord])) \
            if objfindQA_filename is not None else None
        # Run
        sobjs_slit = \
            objs_in_slit(
                image, ivar, thisslit_gpm, 
                slit_left[:,iord], slit_righ[:,iord], 
                spec_min_max=spec_min_max[:,iord],
                inmask=inmask_iord,std_trace=std_in, 
                ncoeff=ncoeff, fwhm=fwhm, use_user_fwhm=use_user_fwhm, maxdev=maxdev,
                hand_extract_dict=hand_extract_dict,  
                nperslit=nperorder, extract_maskwidth=extract_maskwidth,
                snr_thresh=snr_thresh, trim_edg=trim_edg, 
                boxcar_rad=box_radius/plate_scale_ord[iord],
                show_peaks=show_peaks, show_fits=show_single_fits,
                show_trace=show_single_trace, qa_title=qa_title, specobj_dict=specobj_dict,
                objfindQA_filename=ech_objfindQA_filename)
        sobjs.add_sobj(sobjs_slit)

    # Return
    return sobjs


def ech_fof_sobjs(sobjs:specobjs.SpecObjs, 
                  slit_left:np.ndarray, 
                  slit_righ:np.ndarray, 
                  plate_scale_ord:np.ndarray, 
                  fof_link:float=1.5):
    """
    Links together objects previously found using a 
    friends-of-friends algorithm on fractional order position.

    Each source from each order is then assigned an obj_id value

    Args:
        sobjs (:class:`~pypeit.specobj.SpecObj`):
            Previously found objects.  
            There needs to be at least 1.
        slit_left (`numpy.ndarray`_):
            Left boundary of orders to be extracted (given as floating point
            pixels).  Shape is (nspec, norders), where norders is the total
            number of traced echelle orders.
        slit_righ (`numpy.ndarray`_):
            Right boundary of orders to be extracted (given as floating point
            pixels).  Shape is (nspec, norders), where norders is the total
            number of traced echelle orders.
        plate_scale_ord (`numpy.ndarray`_):
            An array with shape (norders,) providing the plate 
            scale of each order in arcsec/pix, 
        fof_link (:obj:`float`, optional):
            Friends-of-friends linking length in arcseconds used to link
            together traces across orders. The routine links together at
            the same fractional slit position and links them together
            with a friends-of-friends algorithm using this linking
            length.

    Returns:
        `numpy.ndarray`_: An array containing the 
        object IDs of the objects linked together.
        This array is aligned with sobjs
    """
    # Prep
    norders = slit_left.shape[1]
    slit_width = slit_righ - slit_left
    nfound = len(sobjs)

    #
    FOF_frac = fof_link/(np.median(np.median(slit_width,axis=0)*plate_scale_ord))
    # Run the FOF. We use fake coordinates
    fracpos = sobjs.SPAT_FRACPOS
    ra_fake = fracpos/1000.0  # Divide all angles by 1000 to make geometry euclidian
    dec_fake = np.zeros_like(fracpos)
    if nfound>1:
        inobj_id, multobj_id, firstobj_id, nextobj_id \
                = pydl.spheregroup(ra_fake, dec_fake, FOF_frac/1000.0)
        # Modify to 1-based indexing
        obj_id_init = inobj_id + 1
    elif nfound==1:
        obj_id_init = np.ones(1,dtype='int')
    else:
        msgs.error('No objects found in ech_fof_sobjs. Should not have called this routine')

    uni_obj_id_init, uni_ind_init = np.unique(obj_id_init, return_index=True)

    # Now loop over the unique objects and check that there is only one object per order. If FOF
    # grouped > 1 objects on the same order, then this will be popped out as its own unique object
    obj_id = obj_id_init.copy()
    nobj_init = len(uni_obj_id_init)
    for iobj in range(nobj_init):
        for iord in range(norders):
            on_order = (obj_id_init == uni_obj_id_init[iobj]) & (
                sobjs.ECH_ORDERINDX == iord)
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

    # Finish
    uni_obj_id, uni_ind = np.unique(obj_id, return_index=True)
    nobj = len(uni_obj_id)
    msgs.info('FOF matching found {:d}'.format(nobj) + ' unique objects')

    return obj_id

def ech_fill_in_orders(sobjs:specobjs.SpecObjs, 
                  slit_left:np.ndarray, 
                  slit_righ:np.ndarray, 
                  order_vec:np.ndarray,
                  order_gpm:np.ndarray,
                  obj_id:np.ndarray,
                  slit_spat_id:np.ndarray,
                  std_trace:specobjs.SpecObjs=None,
                  show:bool=False):
    """
    For objects which were only found on some orders, the standard (or
        the slit boundaries) are placed at the appropriate fractional
        position along the order.

    This routine:

        - Assigns each specobj a fractional order position and an obj_id number.

        - Fills in missing objects. Fit the fraction slit position of the good
          orders where an object was found and use that fit to predict the
          fractional slit position on the bad orders where no object was found.

        - Now loop over the orders and add objects on the orders for 
          which the current object was not found.


    Args:
        sobjs (:class:`~pypeit.specobj.SpecObj`):
            Objects found on some orders thus far
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
        order_gpm (`numpy.ndarray`_):
            Boolean array indicating which orders are good (True),
            i.e. have good calibrations (wavelengths, etc.).  Shape is
        obj_id (`numpy.ndarray`_):
            Object IDs of the objects linked together.
        slit_spat_id (`numpy.ndarray`_):
            slit_spat values (spatial position 1/2 way up the detector)
            for the orders 
        std_trace (:class:`~pypeit.specobjs.SpecObjs`, optional): 
            Standard star objects (including the traces)
            Defaults to None.
        show (bool, optional): 
            Plot diagnostics related to filling the
            missing orders

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`:  A new SpecObjs object with the
        filled in orders
    """
    # Prep
    nfound = len(sobjs)
    uni_obj_id, uni_ind = np.unique(obj_id, return_index=True)
    nobj = len(uni_obj_id)
    fracpos = sobjs.SPAT_FRACPOS

    # Prep
    ngd_orders = np.sum(order_gpm)
    gd_orders = order_vec[order_gpm]
    slit_width = slit_righ - slit_left

    # Check standard star
    if std_trace is not None and std_trace.shape[1] != ngd_orders:
        msgs.error('Standard star trace does not match the number of orders in the echelle data.')

    # For traces
    nspec = slit_left.shape[0]
    spec_vec = np.arange(nspec)
    slit_spec_pos = nspec/2.0
    specmid = nspec // 2

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
        #for iord in range(norders):
        for iord in order_vec[order_gpm]:
            on_order = (obj_id == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDER == iord)
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
        if (nthisobj_id > 3) and (nthisobj_id<ngd_orders):
            thisorderindx = sobjs_align[indx_obj_id].ECH_ORDERINDX
            thisorder = sobjs_align[indx_obj_id].ECH_ORDER
            # Allow for masked orders
            xcen_good = (sobjs_align[indx_obj_id].TRACE_SPAT).T
            slit_frac_good = (xcen_good-slit_left[:,thisorderindx])/slit_width[:,thisorderindx]
            # Fractional slit position averaged across the spectral direction for each order
            frac_mean_good = np.mean(slit_frac_good, 0)
            # Perform a  linear fit to fractional slit position
            #TODO Do this as a S/N weighted fit similar to what is now in the pca_trace algorithm?
            #msk_frac, poly_coeff_frac = fitting.robust_fit(order_vec[goodorder], frac_mean_good, 1,
            pypeitFit = fitting.robust_fit(
                #order_vec[goodorder], frac_mean_good, 1,
                thisorder, frac_mean_good, 1,
                function='polynomial', maxiter=20, lower=2, upper=2,
                use_mad= True, sticky=False,
                minx = order_vec.min(), maxx=order_vec.max())
            # Fill
            goodorder = np.in1d(gd_orders, thisorder)
            badorder = np.invert(goodorder)
            frac_mean_new = np.zeros(gd_orders.size)
            frac_mean_new[badorder] = pypeitFit.eval(gd_orders[badorder])
            frac_mean_new[goodorder] = frac_mean_good
            # TODO This QA needs some work
            if show:
                frac_mean_fit = pypeitFit.eval(gd_orders)
                plt.plot(gd_orders[goodorder][pypeitFit.bool_gpm], frac_mean_new[goodorder][pypeitFit.bool_gpm], 'ko', mfc='k', markersize=8.0, label='Good Orders Kept')
                plt.plot(gd_orders[goodorder][np.invert(pypeitFit.bool_gpm)], frac_mean_new[goodorder][np.invert(pypeitFit.bool_gpm)], 'ro', mfc='k', markersize=8.0, label='Good Orders Rejected')
                plt.plot(gd_orders[badorder], frac_mean_new[badorder], 'ko', mfc='None', markersize=8.0, label='Predicted Bad Orders')
                plt.plot(gd_orders,frac_mean_new,'+',color='cyan',markersize=12.0,label='Final Order Fraction')
                plt.plot(gd_orders, frac_mean_fit, 'r-', label='Fractional Order Position Fit')
                plt.xlabel('Order Index', fontsize=14)
                plt.ylabel('Fractional Slit Position', fontsize=14)
                plt.title('Fractional Slit Position Fit')
                plt.legend()
                plt.show()
        else:
            frac_mean_new = np.full(gd_orders.size, uni_frac[iobj])


        # Now loop over the orders and add objects on the orders for 
        #  which the current object was not found
        for iord in range(order_vec.size):
            iorder = order_vec[iord]
            if iorder not in gd_orders:
                continue
            # Is the current object detected on this order?
            on_order = (sobjs_align.ECH_OBJID == uni_obj_id[iobj]) & (
                sobjs_align.ECH_ORDER == iorder)
            num_on_order = np.sum(on_order)
            if num_on_order == 0:
                msgs.info(f"Adding object={uni_obj_id[iobj]} to order={iorder}")
                # If it is not, create a new sobjs and add to sobjs_align and assign required tags
                thisobj = specobj.SpecObj('Echelle', sobjs_align[0].DET,
                                             OBJTYPE=sobjs_align[0].OBJTYPE,
                                             ECH_ORDERINDX=iord,
                                             ECH_ORDER=iorder)
                #thisobj.ECH_ORDERINDX = iord
                #thisobj.ech_order = order_vec[iord]
                thisobj.SPAT_FRACPOS = uni_frac[iobj]
                # Assign traces using the fractional position fit above
                if std_trace is not None:
                    x_trace = np.interp(slit_spec_pos, spec_vec, std_trace[:,iord])
                    shift = np.interp(
                        slit_spec_pos, spec_vec, slit_left[:,iord] + 
                        slit_width[:,iord]*frac_mean_new[iord]) - x_trace
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
                thisobj.hand_extract_flag = sobjs_align[imin].hand_extract_flag
                thisobj.maskwidth = sobjs_align[imin].maskwidth
                thisobj.smash_peakflux = sobjs_align[imin].smash_peakflux
                thisobj.smash_snr = sobjs_align[imin].smash_snr
                thisobj.BOX_RADIUS = sobjs_align[imin].BOX_RADIUS
                thisobj.ECH_FRACPOS = uni_frac[iobj]
                thisobj.ECH_OBJID = uni_obj_id[iobj]
                thisobj.OBJID = uni_obj_id[iobj]
                thisobj.SLITID = slit_spat_id[iord]
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
    # Return
    return sobjs_align

def ech_cutobj_on_snr(
    sobjs_align, image:np.ndarray, ivar:np.ndarray, 
    slitmask:np.ndarray, gd_orders:np.ndarray, 
    plate_scale_ord:np.ndarray, max_snr:float=2.0,
    nperorder:int=2, min_snr:float=1.0, 
    nabove_min_snr:int=2,
    box_radius:float=2.0, inmask:np.ndarray=None):
    """Cut down objects based on S/N

    This routine:

        - Loops over the objects and perform a quick and dirty extraction to
          assess S/N.

        - Purge objects with low SNR that don't show up in enough orders, sort
          the list of objects with respect to obj_id and orderindx.
        
        - Loop over objects from highest SNR to lowest SNR. Apply the S/N
          constraints.  Once we hit the maximum number objects requested exit,
          except keep any hand apertures that were requested.

    Args:
        sobjs_align (:class:`~pypeit.specobj.SpecObj`):
            Previously found objects
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
        gd_orders (`numpy.ndarray`_):
            `int` array of good orders 
        plate_scale_ord (`numpy.ndarray`_):
            An array with shape (norders,) providing the plate 
            scale of each order in arcsec/pix, 
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
        box_radius (:obj:`float`, optional):
            Box_car extraction radius in arcseconds to assign to each detected
            object and to be used later for boxcar extraction. In this method
            ``box_radius`` is converted into pixels using ``plate_scale``.
            ``box_radius`` is also used for SNR calculation and trimming.
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Must have the same shape as
            ``image``.  If None, all pixels in ``slitmask`` with non-negative
            values are considered good.

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`: The final set of objects
    """

    allmask = slitmask > -1
    if inmask is None:
        inmask = allmask
    # Prep
    nspec = image.shape[0]
    norders = gd_orders.size                     
    uni_obj_id = np.unique(sobjs_align.ECH_OBJID)
    nobj = uni_obj_id.size

    # Loop over the objects and perform a quick and dirty extraction to assess S/N.
    varimg = utils.calc_ivar(ivar)
    flux_box = np.zeros((nspec, norders, nobj))
    ivar_box = np.zeros((nspec, norders, nobj))
    mask_box = np.zeros((nspec, norders, nobj))
    SNR_arr = np.zeros((norders, nobj))
    slitfracpos_arr = np.zeros((norders, nobj))
    hand_flag = np.zeros(nobj, dtype=bool)

    for iobj in range(nobj):
        # Hand extraction?
        mt_obj = sobjs_align.ECH_OBJID == uni_obj_id[iobj]
        if np.any(sobjs_align[mt_obj].hand_extract_flag):
            hand_flag[iobj] = True
        # SNR
        for iord in range(norders):
            iorder_vec = gd_orders[iord]
            indx = sobjs_align.slitorder_objid_indices(
                iorder_vec, uni_obj_id[iobj])
            #indx = (sobjs_align.ECH_OBJID == uni_obj_id[iobj]) & (sobjs_align.ECH_ORDERINDX == iord)
            #spec = sobjs_align[indx][0]
            inmask_iord = inmask & (slitmask == sobjs_align[indx].SLITID)# gdslit_spat[iord])
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
            mean, med_sn, stddev = astropy.stats.sigma_clipped_stats(
                flux_box[mask_tmp,iord,iobj]*np.sqrt(ivar_box[mask_tmp,iord,iobj]),
                sigma_lower=5.0,sigma_upper=5.0
            )
            # ToDO assign this to sobjs_align for use in the extraction
            SNR_arr[iord,iobj] = med_sn
            sobjs_align[indx][0].ech_snr = med_sn
            # For hand extractions
            slitfracpos_arr[iord,iobj] = sobjs_align[indx][0].SPAT_FRACPOS

    # Purge objects with low SNR that don't show up in enough orders, sort the list of objects with respect to obj_id
    # and orderindx
    keep_obj = np.zeros(nobj,dtype=bool)
    # Empty specobjs object to hold the final objects
    sobjs_trim = specobjs.SpecObjs()
    # objids are 1 based so that we can easily asign the negative to negative objects
    iobj_keep = 1
    iobj_keep_not_hand = 1

    ## Loop over objects from highest SNR to lowest SNR. Apply the S/N constraints. Once we hit the maximum number
    # objects requested exit, except keep the hand apertures that were requested.
    isort_SNR_max = np.argsort(np.median(SNR_arr,axis=0))[::-1]
    for iobj in isort_SNR_max:
        hand_ap_flag = hand_flag[iobj]
        SNR_constraint = (SNR_arr[:,iobj].max() > max_snr) or (
            np.sum(SNR_arr[:,iobj] > min_snr) >= nabove_min_snr)
        nperorder_constraint = (iobj_keep-1) < nperorder
        if (SNR_constraint and nperorder_constraint) or hand_ap_flag:
            keep_obj[iobj] = True
            ikeep = sobjs_align.ECH_OBJID == uni_obj_id[iobj]
            sobjs_keep = sobjs_align[ikeep].copy()
            sobjs_keep.ECH_OBJID = iobj_keep
            sobjs_keep.OBJID = iobj_keep
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

    return sobjs_trim


def ech_pca_traces(
    sobjs_final:specobjs.SpecObjs, image:np.ndarray, 
    slitmask:np.ndarray, inmask:np.ndarray, 
    gd_orders:np.ndarray, spec_min_max,
    npca:int=None, coeff_npoly:int=None,
    pca_explained_var:float=99.0, 
    ncoeff:int=5, maxdev:float=2.0, fwhm:float=3.0,
    show_trace:bool=False, show_fits:bool=False, 
    show_pca:bool=False):
    """
    A PCA fit to the traces is performed using the routine pca_fit
    It then applies iterative flux-weighted centroiding to refine the traces

    Args:
        sobjs_final (:class:`~pypeit.specobj.SpecObj`):
            Final set of objects ready for tracing
        image (`numpy.ndarray`_):
            (Floating-point) Image to use for object search with shape (nspec,
            nspat).  The first dimension (nspec) is spectral, and second
            dimension (nspat) is spatial. Note this image can either have the
            sky background in it, or have already been sky subtracted.  Object
            finding works best on sky-subtracted images. Ideally, object finding
            is run in another routine, global sky-subtraction performed, and
            then this code should be run. However, it is also possible to run
            this code on non-sky-subtracted images.
        slitmask (`numpy.ndarray`_):
            Integer image indicating the pixels that belong to each order.
            Pixels that are not on an order have value -1, and those that are on
            an order have a value equal to the slit number (i.e. 0 to nslits-1
            from left to right on the image).  Shape must match ``image``,
            (nspec, nspat).
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Must have the same shape as
            ``image``.  If None, all pixels in ``slitmask`` with non-negative
            values are considered good.
        gd_orders (`numpy.ndarray`_):
            `int` array of good orders 
        spec_min_max (`numpy.ndarray`_): _description_
            `float` array of shape (2, norders) with the minimum and maximum
            spectral value for each order
        npca (:obj:`int`, optional):
            Number of PCA components to keep during PCA decomposition of the
            object traces.  If None, the number of components set by requiring
            the PCA accounts for approximately 99% of the variance.
        coeff_npoly (:obj:`int`, optional):
            Order of polynomial used for PCA coefficients fitting.  If None,
            value set automatically, see
            :func:`~pypeit.tracepca.pca_trace_object`.
        pca_explained_var (:obj:`float`, optional):
            The percentage (i.e., not the fraction) of the variance in the data
            accounted for by the PCA used to truncate the number of PCA
            coefficients to keep (see ``npca``). Ignored if ``npca`` is provided
            directly; see :func:`~pypeit.tracepca.pca_trace_object`.
        ncoeff (:obj:`int`, optional):
            Order of polynomial fit to traces.
        maxdev (:obj:`float`, optional):
            Maximum deviation of pixels from polynomial fit to trace
            used to reject bad pixels in trace fitting.
        fwhm (:obj:`float`, optional):
            Estimated fwhm of the objects in pixels
        show_trace (:obj:`bool`, optional):
            Display the object traces on top of the image.
        show_fits (:obj:`bool`, optional):
            Plot trace fitting for final fits using PCA as crutch.
        show_pca (:obj:`bool`, optional):
            Display debugging plots for the PCA decomposition.

    Returns:
        :class:`~pypeit.specobj.SpecObj`: Final set of objects, traced 
    """

    # Prep
    nspec, nspat = image.shape
    norders = gd_orders.size                     
    uni_obj = np.unique(sobjs_final.ECH_OBJID)
    nobj_trim = uni_obj.size
    spec_vec = np.arange(nspec)
    specmid = nspec // 2

    allmask = slitmask > -1
    if inmask is None:
        inmask = allmask

    # Checks
    if gd_orders.size != spec_min_max.shape[1]:
        msgs.error("Number of good orders does not match the number of orders in spec_min_max")

    # Loop over the objects one by one and adjust/predict the traces
    pca_fits = np.zeros((nspec, norders, nobj_trim))

    # Create the trc_inmask for iterative fitting below
    trc_inmask = np.zeros((nspec, norders), dtype=bool)
    for iord in range(norders):
        trc_inmask[:,iord] = (
            spec_vec >= spec_min_max[0,iord]) & (
                spec_vec <= spec_min_max[1,iord])

    for iobj in range(nobj_trim):
        indx_obj_id = sobjs_final.ECH_OBJID == (iobj + 1)
        # PCA predict all the orders now (where we have used the standard or slit boundary for the bad orders above)
        msgs.info('Fitting echelle object finding PCA for object {:d}/{:d} with median SNR = {:5.3f}'.format(
            iobj + 1,nobj_trim,np.median(sobjs_final[indx_obj_id].ech_snr)))
        pca_fits[:,:,iobj] \
                = tracepca.pca_trace_object(
                    sobjs_final[indx_obj_id].TRACE_SPAT.T,
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
            if spec.ech_frac_was_fit and (spec.ech_snr > 1.0):
                    spec.TRACE_SPAT = xfit_gweight[:,iord]
                    spec.SPAT_PIXPOS = spec.TRACE_SPAT[specmid]

    #TODO Put in some criterion here that does not let the fractional position change too much during the iterative
    # tracefitting. The problem is spurious apertures identified on one slit can be pulled over to the center of flux
    # resulting in a bunch of objects landing on top of each other.

    # Set the IDs
    #sobjs_final[:].ECH_ORDER = order_vec[sobjs_final[:].ECH_ORDERINDX]
    #for spec in sobjs_final:
    #    spec.ech_order = order_vec[spec.ECH_ORDERINDX]
    sobjs_final.set_names()

    if show_trace:
        viewer, ch = display.show_image(image*allmask)

        #for spec in sobjs_trim:
        for spec in sobjs_final:
            color = 'red' if spec.ech_frac_was_fit else 'magenta'
            ## Showing the final flux weighted centroiding from PCA predictions
            display.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        for iobj in range(nobj_trim):
            obj_idx = sobjs_final.ECH_OBJID == (iobj + 1) 
            frac = sobjs_final[obj_idx].SPAT_FRACPOS[0]
            for iord in range(norders):
                ## Showing PCA predicted locations before recomputing flux/gaussian weighted centroiding
                display.show_trace(
                    viewer, ch, pca_fits[:,iord, iobj], 
                    str(frac), color='yellow')
                    #str(uni_frac[iobj]), color='yellow')
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


def ech_objfind(image, ivar, slitmask, slit_left, slit_righ, order_vec, slits_bpm, 
                slit_spat_id, spec_min_max,
                det='DET01', inmask=None, 
                fof_link=1.5, plate_scale=0.2,
                std_trace=None, ncoeff=5, npca=None, 
                coeff_npoly=None, max_snr=2.0, min_snr=1.0,
                nabove_min_snr=2, pca_explained_var=99.0, 
                box_radius=2.0, fwhm=3.0,
                use_user_fwhm=False, maxdev=2.0, 
                nperorder=2,
                extract_maskwidth=3.0, snr_thresh=10.0,
                specobj_dict=None, trim_edg=(5,5),
                show_peaks=False, show_fits=False, 
                show_single_fits=False,
                show_trace=False, show_single_trace=False, 
                show_pca=False,
                debug_all=False, objfindQA_filename=None,
                manual_extract_dict=None):
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
        slits_bpm (`numpy.ndarray`_):
            Boolean array selecting orders that should be ignored (i.e., good
            orders are False, bad orders are True).  Shape must be (norders,).
        slit_spat_id (`numpy.ndarray`_):
            slit_spat values (spatial position 1/2 way up the detector)
            for the orders 
        spec_min_max (`numpy.ndarray`_):
            2D array defining the minimum and maximum pixel in the spectral
            direction with useable data for each order.  Shape must be (2,
            norders).  This should only be used for echelle spectrographs for
            which the orders do not entirely cover the detector. PCA tracing
            will re-map the traces such that they all have the same length,
            compute the PCA, and then re-map the orders back.  This improves
            performance for echelle spectrographs by removing the nonlinear
            shrinking of the orders so that the linear pca operation can better
            predict the traces.  Otherwise set the values to -1 and nspec
            where nspec is the number of pixels in the spectral direction.
        det (:obj:`str`, optional):
            The name of the detector containing the object.  Only used if
            ``specobj_dict`` is None.
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Must have the same shape as
            ``image``.  If None, all pixels in ``slitmask`` with non-negative
            values are considered good.
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
        manual_extract_dict : :obj:`dict`, optional
            Dict guiding the manual extraction

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`: Object containing the objects
        detected.
    """
        #debug_all=True
    if debug_all:
        show_peaks = True
        #show_fits = True
        #show_single_fits = True
        show_trace = True
        show_pca = True
        #show_single_trace = True

    if specobj_dict is None:
        specobj_dict = {'SLITID': 999, 
                        'ECH_ORDERINDX': 999,
                        'DET': det, 'OBJTYPE': 'unknown', 
                        'PYPELINE': 'Echelle'}

    # Loop over the orders and find the objects within them (if any)
    order_gpm = np.invert(slits_bpm)
    sobjs_in_orders = ech_findobj_ineach_order(
        image, ivar, slitmask, slit_left, 
        slit_righ, slit_spat_id,
        order_vec, order_gpm,
        spec_min_max, plate_scale,
        det=det,
        inmask=inmask, 
        std_trace=std_trace,
        specobj_dict=specobj_dict,
        snr_thresh=snr_thresh,
        show_peaks=show_peaks, 
        show_single_fits=show_single_fits,
        show_single_trace=show_single_trace,
        extract_maskwidth=extract_maskwidth,
        trim_edg=trim_edg,
        fwhm=fwhm,
        use_user_fwhm=use_user_fwhm,
        nperorder=nperorder,
        maxdev=maxdev,
        box_radius=box_radius,
        objfindQA_filename=objfindQA_filename,
        hand_extract_dict=manual_extract_dict)

    # No sources and no manual?
    if len(sobjs_in_orders) == 0: 
        return sobjs_in_orders

    # Additional work for slits with sources (found or input manually)

    # Friend of friend algorithm to group objects
    obj_id = ech_fof_sobjs(
        sobjs_in_orders, slit_left,
        slit_righ, plate_scale,
        fof_link=fof_link)

    # Fill in Orders
    sobjs_filled = ech_fill_in_orders(
        sobjs_in_orders, 
        slit_left, slit_righ,
        order_vec, order_gpm,
        obj_id, #obj_id[tmp], 
        slit_spat_id,
        std_trace=std_trace)

    # Cut on SNR and number of objects
    sobjs_pre_final = ech_cutobj_on_snr(
        sobjs_filled, image, ivar, slitmask,
        order_vec[order_gpm],
        plate_scale, 
        inmask=inmask,
        nperorder=nperorder,
        max_snr=max_snr,
        min_snr=min_snr,
        nabove_min_snr=nabove_min_snr,
        box_radius=box_radius)

    # PCA
    sobjs_ech = ech_pca_traces(
        sobjs_pre_final, 
        image, slitmask, inmask, 
        order_vec[order_gpm],
        spec_min_max[:, order_gpm],
        coeff_npoly=coeff_npoly,
        ncoeff=ncoeff, npca=npca,
        pca_explained_var=pca_explained_var,
        maxdev=maxdev,
        fwhm=fwhm,
        show_trace=show_trace, show_fits=show_fits, 
        show_pca=show_pca)

    return sobjs_ech


def orig_ech_objfind(image, ivar, slitmask, slit_left, slit_righ, order_vec, maskslits, det='DET01',
                inmask=None, spec_min_max=None, fof_link=1.5, plate_scale=0.2,
                std_trace=None, ncoeff=5, npca=None, coeff_npoly=None, max_snr=2.0, min_snr=1.0,
                nabove_min_snr=2, pca_explained_var=99.0, box_radius=2.0, fwhm=3.0,
                use_user_fwhm=False, maxdev=2.0, hand_extract_dict=None, nperorder=2,
                extract_maskwidth=3.0, snr_thresh=10.0,
                specobj_dict=None, trim_edg=(5,5),
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

    #debug_all=True
    if debug_all:
        show_peaks = True
        #show_fits = True
        #show_single_fits = True
        show_trace = True
        show_pca = True
        #show_single_trace = True
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
        #msgs.error('There is a mismatch between the number of valid orders found by PypeIt and '
        #           'the number expected for this spectrograph.  Unable to continue.  Please '
        #           'submit an issue on Github: https://github.com/pypeit/PypeIt/issues .')

    if spec_min_max is None:
        spec_min_max = np.zeros((2,norders), dtype=int)
        for iord in range(norders):
            ispec, ispat = np.where(slitmask == gdslit_spat[iord])
            spec_min_max[:,iord] = ispec.min(), ispec.max()

    # Setup the plate scale
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
            ispat = int(np.clip(np.round(spat),0,nspat-1))
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
            objs_in_slit(image, ivar, thisslit_gpm, slit_left[:,iord], slit_righ[:,iord], spec_min_max=spec_min_max[:,iord],
                    inmask=inmask_iord,std_trace=std_in, ncoeff=ncoeff, fwhm=fwhm, use_user_fwhm=use_user_fwhm, maxdev=maxdev,
                    hand_extract_dict=new_hand_extract_dict,  nperslit=nperorder, extract_maskwidth=extract_maskwidth,
                    snr_thresh=snr_thresh, trim_edg=trim_edg, boxcar_rad=box_radius/plate_scale_ord[iord],
                    show_peaks=show_peaks, show_fits=show_single_fits,
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
            mean, med_sn, stddev = astropy.stats.sigma_clipped_stats(
                flux_box[mask_tmp,iord,iobj]*np.sqrt(ivar_box[mask_tmp,iord,iobj]),
                sigma_lower=5.0,sigma_upper=5.0
            )
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
        msgs.info('Fitting echelle object finding PCA for object {:d}/{:d} with median SNR = {:5.3f}'.format(
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



def objfind_QA(spat_peaks, snr_peaks, spat_vector, snr_vector, snr_thresh, qa_title, peak_gpm,
               near_edge_bpm, nperslit_bpm, objfindQA_filename=None, show=False):
    """
    Utility routine for making object finding QA plots.

    Args:
        spat_peaks (`numpy.ndarray`_):
            Array of locations in the spatial direction at which objects were
            identified. Shape is ``(npeaks,)``,  where ``npeaks`` is the number
            of peaks identified.
        snr_peaks (`numpy.ndarray`_):
            S/N ratio in the spectral direction after collapsing along the
            spectral direction, evaluated at the location of each spatial peak.
            Shape must match ``spat_peaks``.
        spat_vector (`numpy.ndarray`_):
            A 1D array of spatial locations along the slit. Shape is
            ``(nsamp,)``, where ``nsamp`` is the number of spatial pixels
            defined by the slit edges.
        snr_vector (`numpy.ndarray`_):
            A 1D array with the S/N ratio sampled along the slit at each spatial
            location (i.e., spectral direction has been smashed out) defined by
            ``spat_vector``.  Shape must match ``spat_vector``.
        snr_thresh (:obj:`float`):
            The threshold S/N ratio adopted by the object finding.
        qa_title (:obj:`str`):
            Title for the QA file plot.
        peak_gpm (`numpy.ndarray`_):
            Boolean array containing a good pixel mask for each peak indicating
            whether it will be used as an object (True) or not (False).  Shape
            must match ``spat_peaks``.
        near_edge_bpm (`numpy.ndarray`_):
            A bad pixel mask (True is masked, False is unmasked) indicating
            which objects are masked because they are near the slit edges.
            Shape must match ``spat_peaks``.
        nperslit_bpm (`numpy.ndarray`_):
            A bad pixel mask (True is masked, False is unmasked) indicating
            which objects are masked because they exceed the maximum number of
            objects (see :func:`objs_in_slit` parameter ``nperslit``) that were
            specified as being on this slit.
        objfindQA_filename (:obj:`str`, optional):
            Output filename for the QA plot.  If None, plot is not saved.
        show (:obj:`bool`, optional):
            If True, show the plot as a matplotlib interactive plot.

    """

    plt.plot(spat_vector, snr_vector, drawstyle='steps-mid', color='black', label = 'Collapsed SNR (FWHM convol)')
    plt.hlines(snr_thresh,spat_vector.min(),spat_vector.max(), color='red',linestyle='--',
               label='SNR_THRESH={:5.3f}'.format(snr_thresh))
    if np.any(peak_gpm):
        plt.plot(spat_peaks[peak_gpm], snr_peaks[peak_gpm], color='red', marker='o', markersize=10.0,
                 mfc='lawngreen', fillstyle='full',linestyle='None', zorder = 10,label='{:d} Good Objects'.format(np.sum(peak_gpm)))
    if np.any(near_edge_bpm):
        plt.plot(spat_peaks[near_edge_bpm], snr_peaks[near_edge_bpm], color='red', marker='o', markersize=10.0,
                 mfc='cyan', fillstyle='full', linestyle='None', zorder = 10,label='{:d} Rejected: Near Edge'.format(np.sum(near_edge_bpm)))
    if np.any(nperslit_bpm):
        plt.plot(spat_peaks[nperslit_bpm], snr_peaks[nperslit_bpm], color='red', marker='o', markersize=10.0,
                 mfc='yellow', fillstyle='full', linestyle='None', zorder = 10,label='{:d} Rejected: Nperslit'.format(np.sum(nperslit_bpm)))
    plt.legend()
    plt.xlabel('Approximate Spatial Position (pixels)')
    plt.ylabel('SNR')
    plt.title(qa_title)
    #plt.ylim(np.fmax(snr_vector.min(), -20.0), 1.3*snr_vector.max())
    fig = plt.gcf()
    if show:
        plt.show()
    if objfindQA_filename is not None:
        fig.savefig(objfindQA_filename, dpi=400)
    plt.close('all')

def get_fwhm(fwhm_in, nsamp, smash_peakflux, spat_fracpos, flux_smash_smth):
    """
    Utility routine to measure the FWHM of an object trace from the spectrally
    collapsed flux profile by determining the locations along the spatial
    direction where this profile reaches have its peak value.

    Args:
        fwhm_in (:obj:`float`):
           Best guess for the FWHM of this object.
        nsamp (:obj:`int`):
           Number of pixels along the spatial direction.
        smash_peakflux (:obj:`float`):
            The peak flux in the 1d (spectrally collapsed) flux profile at the
            object location.
        spat_fracpos (:obj:`float`):
            Fractional spatial position along the slit where the object is
            located and at which the ``flux_smash_smth`` array has values
            provided by ``smash_peakflux`` (see above and below).
        flux_smash_smth (`numpy.ndarray`_):
            A 1D array with the flux averaged along the spectral direction at
            each location along the slit in the spatial direction location.
            Shape is ``(nsamp,)``.

    Returns:
        :obj:`float`: The FWHM determined from the object flux profile, unless
        the FWHM could not be found from the profile, in which case the input
        guess (``fwhm_in``) is simply returned.
    """

    # Determine the fwhm max
    yhalf = 0.5*smash_peakflux
    xpk = spat_fracpos*nsamp
    x0 = int(np.rint(xpk))
    # TODO It seems we have two codes that do similar things, i.e. findfwhm in arextract.py. Could imagine having one
    # Find right location where smash profile croses yhalf
    if x0 < (int(nsamp) - 1):
        ind_righ, = np.where(flux_smash_smth[x0:] < yhalf)
        if len(ind_righ) > 0:
            i2 = ind_righ[0]
            if i2 == 0:
                xrigh = None
            else:
                xarr_righ = x0 + np.array([i2 - 1, i2], dtype=float)
                xrigh_int = scipy.interpolate.interp1d(flux_smash_smth[x0 + i2 - 1:x0 + i2 + 1], xarr_righ,
                                                       assume_sorted=False, bounds_error=False,
                                                       fill_value=(xarr_righ[0], xarr_righ[1]))
                xrigh = xrigh_int([yhalf])[0]
        else:
            xrigh = None
    else:
        xrigh = None
    # Find left location where smash profile crosses yhalf
    if x0 > 0:
        ind_left, = np.where(flux_smash_smth[0:np.fmin(x0 + 1, int(nsamp) - 1)] < yhalf)
        if len(ind_left) > 0:
            i1 = (ind_left[::-1])[0]
            if i1 == (int(nsamp) - 1):
                xleft = None
            else:
                xarr_left = np.array([i1, i1 + 1], dtype=float)
                xleft_int = scipy.interpolate.interp1d(flux_smash_smth[i1:i1 + 2], xarr_left,
                                                       assume_sorted=False, bounds_error=False,
                                                       fill_value=(xarr_left[0], xarr_left[1]))
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
        fwhm_out = np.sqrt(np.fmax(fwhm_measure ** 2 - fwhm_in ** 2, (fwhm_in / 2.0) ** 2))  # Set a floor of fwhm/2 on fwhm
    else:
        fwhm_out = fwhm_in

    return fwhm_out


def objs_in_slit(image, ivar, thismask, slit_left, slit_righ, 
                 inmask=None, fwhm=3.0,
                 sigclip_smash=5.0, use_user_fwhm=False, boxcar_rad=7.,
                 maxdev=2.0, spec_min_max=None, hand_extract_dict=None, std_trace=None,
                 ncoeff=5, nperslit=None, snr_thresh=10.0, trim_edg=(5,5),
                 extract_maskwidth=4.0, specobj_dict=None, find_min_max=None,
                 show_peaks=False, show_fits=False, show_trace=False,
                 debug_all=False, qa_title='objfind', objfindQA_filename=None):
    """
    Find the location of objects in a slitmask slit or a echelle order.

    The algorithm for this function is:

        - Rectify the image by extracting along the edge traces.

        - Compute the sigma-clipped mean spatial profile and its variance by
          collapsing the image along the spectral direction to constuct a S/N
          spatial profile of the slit/order.

        - Smooth the S/N profile by a Gaussian with the provided FWHM (see
          ``fwhm``) and detect peaks in the smoothed profile using
          :func:`~pypeit.core.arc.detect_lines`.  Ignore peaks found near the
          slit edges, and limit the number of peaks to the number requested (see
          ``nperslit``).

        - Instantiate a :class:`~pypeit.specobj.SpecObj` object for each valid
          detection, construct preliminary spectral traces for them, and estimate
          the object spatial FWHM if one is not provided (see
          ``use_user_fwhm``).

        - For automatically identified objects (i.e., not manual extractions),
          improve the object trace by fitting the spatial position of the peak
          as a function of wavelength.  For manual apertures, use either the
          form of the brightest object on the slit, the trace of a standard star
          (see ``std_trace``), or the left edge trace to set the trace for the
          object.  Finally, remove automatically identified objects that overlap
          with manually defined extraction apertures.

    At the end of this function, the list of objects is ready for extraction.

    **Revision History:**
        
        - 10-Mar-2005 -- First version written by D. Schlegel, LBL
        - 2005-2018 -- Improved by J. F. Hennawi and J. X. Prochaska
        - 23-June-2018 -- Ported to python by J. F. Hennawi and significantly
          improved
        - 01-Feb-2022 -- Skymask stripped out by JXP

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
        thismask (`numpy.ndarray`_):
            Boolean mask image selecting pixels associated with the slit/order
            to search for objects on (True means on the slit/order).  Shape must
            match ``image``.
        slit_left (`numpy.ndarray`_):
            Left boundary of a single slit/orders to be extracted (given as
            floating point pixels).  Shape is (nspec,).
        slit_righ (`numpy.ndarray`_):
            Right boundary of a single slit/orders to be extracted (given as
            floating point pixels).  Shape is (nspec,).
        inmask (`numpy.ndarray`_, optional):
            Good-pixel mask for input image.  Shape must match ``image``.  If
            None, set to be the same as ``thismask``.
        fwhm (:obj:`float`, optional):
            Estimated FWHM of the objects in pixels
        sigclip_smash (:obj:`float`, optional):
            Sigma clipping threshold passed to
            `astropy.stats.sigma_clipped_stats`_ to compute average slit
            emission profile by averaging the (rectified) image along the
            spatial direction.
        use_user_fwhm (:obj:`bool`, optional):
            If True, ``PypeIt`` will use the spatial profile FWHM input by the
            user (see ``fwhm``) rather than determine the spatial FWHM from the
            smashed spatial profile via the automated algorithm.
        boxcar_rad (:obj:`float`, optional):
            Box_car extraction radius *in pixels* to assign to each detected
            object and to be used later for boxcar extraction. 
        maxdev (:obj:`float`, optional):
            Maximum deviation of pixels from polynomial fit to trace
            used to reject bad pixels in trace fitting.
        spec_min_max (:obj:`tuple`, optional):
            2-tuple defining the minimum and maximum pixel in the spectral
            direction with useable data for this slit/order.  If None, the
            values will be determined automatically from ``thismask``. Either
            element of the tuple can also None, which will then default to using
            the full min or max over which the slit is defined from
            ``thismask``.
        hand_extract_dict (:obj:`dict`, optional):
            Dictionary with info on manual extraction; see
            :class:`~pypeit.manual_extract.ManualExtractionObj`.
        std_trace (`numpy.ndarray`_, optional):
            Vector with the standard star trace, which is used as a crutch for
            tracing.  Shape must be (nspec,).  If None, the slit boundaries are
            used as the crutch.
        ncoeff (:obj:`int`, optional):
            Order of polynomial fit to traces.
        nperslit (:obj:`int`, optional):
            Maximum number of objects to find.  If there are more detections
            than this number, the code will select the ``nperslit`` most
            significant detections. However, hand apertures will always be
            returned and do not count against this budget.
        snr_thresh (:obj:`float`, optional):
            SNR threshold for finding objects
        trim_edg (:obj:`tuple`, optional):
            A two-tuple of integers or floats used to ignore objects within this
            many pixels of the left and right slit boundaries, respectively.
        extract_maskwidth (:obj:`float`, optional):
            Determines the initial size of the region in units of FWHM that will
            be used for local sky subtraction; See :func:`objs_in_slit` and
            :func:`~pypeit.core.skysub.local_skysub_extract`.
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
                                'OBJTYPE': 'unknown', 'PYPELINE': 'Multislit'}

        find_min_max (:obj:`tuple`, optional):
            2-tuple of integers that defines the minimum and maximum pixel
            location of your *object* in the spectral direction on the detector.
            It is only used for object finding.  This parameter is helpful if
            your object only has emission lines or at high redshift and the
            trace only shows in part of the detector.  Either element of the
            tuple can be None, which will then default to using the full range
            over which the slit is defined. This is distinct from
            ``spec_min_max`` in that ``spec_min_max`` indicates the range of the
            slit/order on the detector, whereas ``find_min_max`` indicates the
            range to be used for object finding. If None or if either member of
            the tuple is None, it will default to the values of
            ``spec_min_max``.
        show_peaks (:obj:`bool`, optional):
            Plot the QA of the object peak finding in each order.
        show_fits (:obj:`bool`, optional):
            Plot trace fitting for final fits using PCA as crutch.
        show_trace (:obj:`bool`, optional):
            Display the object traces on top of the image.
        debug_all (:obj:`bool`, optional):
            Show all the debugging plots.  If True, this also overrides any
            provided values for ``show_peaks``, ``show_fits``, and
            ``show_trace``, setting them to True.
        qa_title (:obj:`str`, optional):
            Title to be printed in the QA plots
        objfindQA_filename (:obj:`str`, optional):
            Full path (directory and filename) for the object profile QA plot.
            If None, not plot is produced and saved.

    Returns:
        :class:`~pypeit.specobjs.SpecObjs`: Object containing the objects
        detected.
    """

    #debug_all=True
    if debug_all:
        show_peaks=True
        show_fits = True
        show_trace = True

    if specobj_dict is None:
        specobj_dict = dict(SLITID=999, DET='DET01', OBJTYPE='unknown', PYPELINE='MultiSlit')

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
    ispec, ispat = np.where(thismask)
    spec_min = ispec.min()
    spec_max = ispec.max()
    if spec_min_max is None or np.any([s is None or np.isinf(s) for s in spec_min_max]):
        if spec_min_max is None:
            spec_min_max_out = np.array([spec_min, spec_max])
        else:
            spec_min_max_out = np.array(spec_min_max).copy()
            if spec_min_max_out[0] is None or np.isinf(spec_min_max_out[0]):
                spec_min_max_out[0] = spec_min
            if spec_min_max_out[1] is None or np.isinf(spec_min_max_out[1]):
                spec_min_max_out[1] = spec_max
        spec_min_max_out = np.array(spec_min_max_out).astype(int)
    else:
        spec_min_max_out = np.array(spec_min_max).astype(int)


    # If find_min_max was not passed in, set it to the values for spec_min_max
    if find_min_max is None or np.any([f is None or np.isinf(f) for f in find_min_max]):
        if find_min_max is None:
            find_min_max_out = spec_min_max_out
        else:
            find_min_max_out = np.array(find_min_max).copy()
            if find_min_max_out[0] is None or np.isinf(find_min_max_out[0]):
                find_min_max_out[0] = spec_min_max_out[0]
            if find_min_max_out[1] is None or np.isinf(find_min_max_out[1]):
                find_min_max_out[1] = spec_min_max_out[1]
        find_min_max_out = np.array(find_min_max_out).astype(int)
    else:
        find_min_max_out = np.array(find_min_max).astype(int)

    #totmask = thismask & inmask & np.invert(edgmask)
    #  Smash the image (for this slit) into a single flux vector.  How many pixels wide is the slit at each Y?
    xsize = slit_righ - slit_left
    #nsamp = np.ceil(np.median(xsize)) # JFH Changed 07-07-19
    nsamp = np.ceil(xsize.max())
    # Mask skypixels with 2 fwhm of edge
    left_asym = slit_left[:,None] + np.outer(xsize/nsamp, np.arange(nsamp))
    righ_asym = left_asym + np.outer(xsize/nsamp, np.ones(int(nsamp)))
    # This extract_asymbox_boxcar call rectifies the image along the curved object traces
    gpm_tot = thismask & inmask & (ivar > 0.0)

    image_rect, gpm_rect, npix_rect, ivar_rect = extract.extract_asym_boxcar(image, left_asym, righ_asym, gpm=gpm_tot, ivar=ivar)

    # This smashes out the spatial direction to construct an aggregate sky model
    #sky_mean, sky_median, sky_sig = stats.sigma_clipped_stats(image_rect, mask=np.logical_not(gpm_rect), axis=1, sigma=3.0,
    #                                                          cenfunc='median', stdfunc=utils.nan_mad_std)
    #gpm_sky = np.sum(gpm_rect,axis=1) != 0 & np.isfinite(sky_median)
    #sky_fill_value = np.median(sky_median[gpm_sky]) if np.any(gpm_sky) else 0.0
    #sky_median[np.logical_not(gpm_sky)] = sky_fill_value

    #sky_rect = np.repeat(sky_median[:, np.newaxis], nsamp, axis=1)
    #sky_rect = 0.0*image_rect


    # Apply find_min_max_out
    find_min_max_gpm = np.zeros_like(image_rect, dtype=bool)
    find_min_max_gpm[find_min_max_out[0]: find_min_max_out[1], :] = True
    data = np.ma.MaskedArray(
        image_rect, mask=np.logical_not(gpm_rect & find_min_max_gpm)) # the total gpm = gpm_rect & find_min_max_gpm
    sigclip = astropy.stats.SigmaClip(
        sigma=sigclip_smash, maxiters=25, cenfunc='median', stdfunc=utils.nan_mad_std
    )
    data_clipped, lower, upper = sigclip(data, axis=0, masked=True, return_bounds=True)
    gpm_sigclip = np.logical_not(data_clipped.mask)


    # Compute the average flux over the set of pixels that are not masked by gpm_sigclip
    nsmash = find_min_max_out[1] - find_min_max_out[0] + 1
    npix_smash = np.sum(gpm_sigclip[find_min_max_out[0]:find_min_max_out[1]], axis=0)
    gpm_smash = npix_smash > 0.3*nsmash
    flux_sum_smash = np.sum((image_rect*gpm_sigclip)[find_min_max_out[0]:find_min_max_out[1]], axis=0)
    flux_smash = flux_sum_smash*gpm_smash/(npix_smash + (npix_smash == 0.0))
    flux_smash_mean, flux_smash_med, flux_smash_std = astropy.stats.sigma_clipped_stats(
        flux_smash, mask=np.logical_not(gpm_smash), sigma_lower=3.0, sigma_upper=3.0
    )
    flux_smash_recen = flux_smash - flux_smash_med

    # Return if none found and no hand extraction
    if not np.any(gpm_smash): 
        sobjs = specobjs.SpecObjs()
        if hand_extract_dict is None:
            # Instantiate a null specobj and return
            msgs.info('No objects found automatically.  Consider manual extraction.')
            return sobjs
        else:
            msgs.info('No objects found automatically.')
            
    else:
        # Compute the formal corresponding variance over the set of pixels that are not masked by gpm_sigclip
        var_rect = utils.inverse(ivar_rect)
        var_sum_smash = np.sum((var_rect*gpm_sigclip)[find_min_max_out[0]:find_min_max_out[1]], axis=0)
        var_smash = var_sum_smash/(npix_smash**2 + (npix_smash == 0.0))
        ivar_smash = utils.inverse(var_smash)*gpm_smash
        snr_smash = flux_smash_recen*np.sqrt(ivar_smash)

        # Smooth this SNR image with a Gaussian set by the input fwhm
        gauss_smth_sigma = (fwhm/2.3548)
        snr_smash_smth = scipy.ndimage.gaussian_filter1d(snr_smash, gauss_smth_sigma, mode='nearest')
        flux_smash_smth = scipy.ndimage.gaussian_filter1d(flux_smash_recen, gauss_smth_sigma, mode='nearest')
        # Search for spatial direction peaks in the smoothed snr image
        _, _, x_peaks_out, x_width, x_err, igood, _, _ = arc.detect_lines(
            snr_smash_smth, input_thresh=snr_thresh, fit_frac_fwhm=1.5, fwhm=fwhm, min_pkdist_frac_fwhm=0.75,
            max_frac_fwhm=10.0, cont_subtract=False, debug_peak_find=False)

        x_peaks_all = x_peaks_out[igood]
        #x_peaks_all = arc.detect_peaks(snr_smash_smth, mph=snr_thresh, mpd=fwhm*0.75, show=False)
        snr_peaks_all = np.interp(x_peaks_all, np.arange(nsamp), snr_smash_smth)
        flux_peaks_all = np.interp(x_peaks_all, np.arange(nsamp), flux_smash_smth)
        npeaks_all = len(x_peaks_all)

        near_edge_bpm = (x_peaks_all < trim_edg[0]) | (x_peaks_all > (nsamp - trim_edg[1]))
        npeak_not_near_edge = np.sum(np.logical_not(near_edge_bpm))

        if np.any(near_edge_bpm):
            msgs.warn('Discarding {:d}'.format(np.sum(near_edge_bpm)) +
                    ' at spatial pixels spat = {:}'.format(x_peaks_all[near_edge_bpm]) +
                    ' which land within trim_edg = (left, right) = {:}'.format(trim_edg) +
                    ' pixels from the slit boundary for this nsamp = {:5.2f}'.format(nsamp) + ' wide slit')
            msgs.warn('You must decrease from the current value of trim_edg in order to keep them')
            msgs.warn('Such edge objects are often spurious')


        # If the user requested the nperslit most significant peaks have been requested, then only return these
        if nperslit is not None:
            # If the requested number is less than (the non-edge) number found, mask them out
            if nperslit < npeak_not_near_edge:
                snr_peaks_not_edge = np.sort(snr_peaks_all[np.logical_not(near_edge_bpm)])[::-1]
                snr_thresh_perslit = snr_peaks_not_edge[nperslit-1]
                nperslit_bpm = np.logical_not(near_edge_bpm) & (snr_peaks_all < snr_thresh_perslit)
            else:
                nperslit_bpm = np.zeros(npeaks_all, dtype=bool)
        else:
            nperslit_bpm = np.zeros(npeaks_all, dtype=bool)

        if np.any(nperslit_bpm):
            msgs.warn('Discarding {:d}'.format(np.sum(nperslit_bpm)) +
                    ' at spatial pixels spat = {:} and SNR = {:}'.format(
                        x_peaks_all[nperslit_bpm], snr_peaks_all[nperslit_bpm]) +
                    ' which are below SNR_thresh={:5.3f} set because the maximum number of objects '.format(snr_thresh_perslit) +
                    'requested nperslit={:d} was exceeded'.format(nperslit))

        peaks_gpm = np.logical_not(near_edge_bpm) & np.logical_not(nperslit_bpm)

        spat_vector = slit_left[specmid] + xsize[specmid] * np.arange(nsamp) / nsamp
        spat_peaks = slit_left[specmid] + xsize[specmid] * x_peaks_all/ nsamp

        # TODO: Change this to show_image or something
        if show_peaks:
            # Show rectified image here? Add this to QA
            viewer, ch = display.show_image(image_rect*gpm_rect*np.sqrt(ivar_rect), chname='objs_in_slit_show', cuts=(-5.0,5.0))
        # QA
        objfind_QA(spat_peaks, snr_peaks_all, spat_vector, snr_smash_smth, snr_thresh, qa_title, peaks_gpm,
                near_edge_bpm, nperslit_bpm, objfindQA_filename=objfindQA_filename, show=show_peaks) #show_peaks)

        nobj_reg = np.sum(peaks_gpm)

        # Instantiate a null specobj
        sobjs = specobjs.SpecObjs()
        # Trim to the good peaks
        x_peaks = x_peaks_all[peaks_gpm]
        snr_peaks = snr_peaks_all[peaks_gpm]
        flux_peaks = flux_peaks_all[peaks_gpm]

        # Now create SpecObj objects for all of these and assign preliminary traces to them.
        for iobj in range(nobj_reg):
            thisobj = specobj.SpecObj(**specobj_dict)
            thisobj.SPAT_FRACPOS = x_peaks[iobj]/nsamp
            thisobj.smash_peakflux = flux_peaks[iobj]
            thisobj.smash_snr = snr_peaks[iobj]
            sobjs.add_sobj(thisobj)

        for iobj in range(nobj_reg):
            # Was a standard trace provided? If so, use that as a crutch.
            if std_trace is not None:
                # Print a status message for the first object
                if iobj == 0:
                    msgs.info('Using input STANDARD star trace as crutch for object tracing')

                x_trace = np.interp(specmid, spec_vec, std_trace)
                shift = np.interp(specmid, spec_vec,
                slit_left + xsize * sobjs[iobj].SPAT_FRACPOS) - x_trace
                sobjs[iobj].TRACE_SPAT = std_trace + shift
                # If no standard trace is provided shift left slit boundary over to be initial trace
            else:
                # ToDO make this the average left and right boundary instead. That would be more robust.
                sobjs[iobj].TRACE_SPAT = slit_left + xsize*sobjs[iobj].SPAT_FRACPOS

            sobjs[iobj].trace_spec = spec_vec
            sobjs[iobj].SPAT_PIXPOS = sobjs[iobj].TRACE_SPAT[specmid]
            # Set the idx for any prelminary outputs we print out. These will be updated shortly
            sobjs[iobj].set_name()

            # assign FWHM to all these objects
            if use_user_fwhm:
                sobjs[iobj].FWHM = fwhm
            else:
                sobjs[iobj].FWHM = get_fwhm(fwhm, nsamp, sobjs[iobj].smash_peakflux, sobjs[iobj].SPAT_FRACPOS, flux_smash_smth)

            # assign BOX_RADIUS
            sobjs[iobj].BOX_RADIUS = boxcar_rad

        if len(sobjs) == 0 and hand_extract_dict is None:
            # TODO: Why is this not done way above?
            #  It appears possible to have an initial object detection, but then
            #  have it go away..
            msgs.info('No objects found automatically.  Consider manual extraction.')
            return specobjs.SpecObjs()

    msgs.info("Automatic finding routine found {0:d} objects".format(len(sobjs)))

    # Fit the object traces
    if len(sobjs) > 0:
        msgs.info('Fitting the object traces')
        # Note the transpose is here to pass in the TRACE_SPAT correctly.
        xinit_fweight = np.copy(sobjs.TRACE_SPAT.T)
        spec_mask = (spec_vec >= spec_min_max_out[0]) & (spec_vec <= spec_min_max_out[1])
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
        hand_extract_spec, hand_extract_spat, hand_extract_det, hand_extract_fwhm, \
            hand_extract_boxcar = [hand_extract_dict[key] for key in [
                'spec', 'spat', 'detname', 'fwhm', 'boxcar_rad']]

        # Determine if these hand apertures land on the slit in question
        hand_on_slit = np.where(np.array(thismask[np.rint(hand_extract_spec).astype(int),
                                                  np.rint(hand_extract_spat).astype(int)]))
        hand_extract_spec = hand_extract_spec[hand_on_slit]
        hand_extract_spat = hand_extract_spat[hand_on_slit]
        hand_extract_det  = hand_extract_det[hand_on_slit]
        hand_extract_fwhm = hand_extract_fwhm[hand_on_slit]
        hand_extract_boxcar = hand_extract_boxcar[hand_on_slit]
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
            thisobj.smash_peakflux = np.interp(thisobj.SPAT_FRACPOS*nsamp,np.arange(nsamp), flux_smash_smth) # interpolate from fluxconv
            thisobj.smash_snr = np.interp(thisobj.SPAT_FRACPOS*nsamp,np.arange(nsamp), snr_smash_smth) # interpolate from fluxconv
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
            # assign BOX_RADIUS (pixels!)
            if hand_extract_boxcar[iobj] > 0.:
                thisobj.BOX_RADIUS = hand_extract_boxcar[iobj]
            else:
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

    # Assign the maskwidth
    for iobj in range(nobj):
        sobjs[iobj].maskwidth = extract_maskwidth*sobjs[iobj].FWHM*(1.0 + 0.5*np.log10(np.fmax(sobjs[iobj].smash_snr,1.0)))

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
        # Vet
        if not sobj.ready_for_extraction():
            # embed(header=utils.embed_header())
            msgs.error("Bad SpecObj.  Can't proceed")

    # Return
    return sobjs

