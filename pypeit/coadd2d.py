"""
Module for performing two-dimensional coaddition of spectra.

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst
"""
import os
import copy

from IPython import embed

import numpy as np
from scipy import ndimage
from matplotlib import pyplot as plt
from astropy.table import Table, vstack

from pypeit import msgs
from pypeit import specobjs
from pypeit import slittrace
from pypeit import extraction
from pypeit import find_objects
from pypeit.images import pypeitimage
from pypeit.core import findobj_skymask
from pypeit.core.wavecal import wvutils
from pypeit.core import coadd
#from pypeit.core import parse
from pypeit import calibrations
from pypeit import spec2dobj
from pypeit.core.moment import moment1d
from pypeit.manual_extract import ManualExtractionObj



class CoAdd2D:

    """
    Main routine to run the extraction for 2d coadds.

    Algorithm steps are as follows:
        - Fill this in.

    This performs 2d coadd specific tasks, and then also performs some
    of the tasks analogous to the pypeit.extract_one method. Docs coming
    soon....
    """
    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, spec2dfiles, spectrograph, par, det=1, offsets=None, weights='auto',
                     spec_samp_fact=1.0, spat_samp_fact=1.0,
                     sn_smooth_npix=None, bkg_redux=False, find_negative=False, show=False,
                     show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """
        Instantiate the subclass appropriate for the provided spectrograph.

        The class to instantiate must match the ``pypeline``
        attribute of the provided ``spectrograph``, and must be a
        subclass of :class:`CoAdd2D`; see the parent class
        instantiation for parameter descriptions.

        Returns:
            :class:`CoAdd2D`: One of the subclasses with
            :class:`CoAdd2D` as its base.
        """

        return next(c for c in cls.__subclasses__() 
                    if c.__name__ == (spectrograph.pypeline + 'CoAdd2D'))(
                        spec2dfiles, spectrograph, par, det=det, offsets=offsets, weights=weights,
                        spec_samp_fact=spec_samp_fact, spat_samp_fact=spat_samp_fact,
                        sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative,
                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug,
                        **kwargs_wave)

    def __init__(self, spec2d, spectrograph, par, det=1, offsets=None, weights='auto',
                 spec_samp_fact=1.0, spat_samp_fact=1.0,
                 sn_smooth_npix=None, bkg_redux=False, find_negative=False, show=False,
                 show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """

        Args:
            spec2d_files (:obj:`list`):
                List of spec2d files or a list of
                :class:`~pypeit.spec2dobj.Spec2dObj` objects.
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
                The instrument used to collect the data to be reduced.
            par (:class:`~pypeit.par.parset.ParSet`):
                Processing parameters.
            det (:obj:`int`, :obj:`tuple`, optional):
                The 1-indexed detector number(s) to process.  If a tuple, it
                must include detectors viable as a mosaic for the provided
                spectrograph; see
                :func:`~pypeit.spectrographs.spectrograph.Spectrograph.allowed_mosaics`.
            offsets (`numpy.ndarray`_, optional):
                Spatial offsets to be applied to each image before coadding. For
                the default mode of None, images are registered automatically
                using the trace of the brightest object. Input offsets are not
                yet supported.
            weights (:obj:`str`, :obj:`list`, `numpy.ndarray`_):
                Mode for the weights used to coadd images. Options are
                ``'auto'`` (default), ``'uniform'``, or a list/array of weights
                with ``shape = (nexp,)`` to apply to the image. Note ``'auto'``
                is not allowed if offsets are input.
            spec_samp_fact (:obj:`float`, optional):
                Make the wavelength grid sampling finer (``spec_samp_fact`` less
                than 1.0) or coarser (``spec_samp_fact`` greather than 1.0) by
                this sampling factor. This basically multiples the 'native'
                spectral pixel size by ``spec_samp_fact``, i.e. the units of
                ``spec_samp_fact`` are pixels.
            spat_samp_fact (:obj:`float`, optional):
                Make the spatial sampling finer (``spat_samp_fact`` less
                than 1.0) or coarser (``spat_samp_fact`` greather than 1.0) by
                this sampling factor. This basically multiples the 'native'
                spatial pixel size by ``spat_samp_fact``, i.e. the units of
                ``spat_samp_fact`` are pixels.
            sn_smooth_npix (:obj:`int`, optional):
                Number of pixels to median filter by when computing S/N used to
                decide how to scale and weight spectra. If set to None, the code
                will simply take 10% of the image size in the spectral
                direction.  TODO: for truncated echelle orders we should be
                doing something more intelligent.
            bkg_redux (:obj:`bool`, optional):
                If True, the sciImg has been subtracted by
                a background image (e.g. standard treatment in the IR)
                This parameter is passed to pypeit.reduce for determining the reduction steps.
            find_negative (:obj:`bool`, optional):
                Do the images have negative trace as would be generated by
                differencing two science frames? This parameter is passed to
                pypeit.reduce for determining the reduction steps. If True, then
                find and mask negative object traces.
            show (:obj:`bool`, optional):
                Show results in ginga
            show_peaks (:obj:`bool`, optional):
                Plot the QA for the object finding algorithm peak finding to the
                screen.
            debug_offset (:obj:`bool`, optional):
                Plot QA for debugging the automatic determination of offsets to
                the screen.
            debug (:obj:`bool`, optional):
                Show QA for debugging.
            **kwargs_wave
                Keyword arguments passed to
                :func:`pypeit.core.coadd.get_wave_grid`, which determine how the
                wavelength grid is created for the 2d coadding.
        """

        # Use Cases:
        # offsets
        #    1) offsets is None -- auto compute offsets from brightest object (if exists)
        #    2) offsets not None (i.e. a list) - use them
        #    -------------- only for Multislit --------------
        #    3) offsets = 'maskdef_offsets' - use `maskdef_offset` saved in SlitTraceSet
        # ===============================================================================
        # weights
        #    1) weights = 'auto' -- if brightest object exists auto compute weights,
        #                           otherwise use uniform weights
        #    2) weights = 'uniform' -- use uniform weights
        #    3) weights is a list - use them

        self.spec2d = spec2d
        self.spectrograph = spectrograph
        self.par = par

        # This can be a single integer for a single detector or a tuple for
        # multiple detectors placed in a mosaic.
        self.det = det

        # This is the string name of the detector or mosaic used when saving the
        # processed data to PypeIt's main output files
        self.detname = self.spectrograph.get_det_name(self.det)

        self.weights = weights
        self.spec_samp_fact = spec_samp_fact
        self.spat_samp_fact = spat_samp_fact
        self.bkg_redux = bkg_redux
        self.find_negative = find_negative
        self.show = show
        self.show_peaks = show_peaks
        self.debug_offsets = debug_offsets
        self.debug = debug
        self.offsets = None
        self.stack_dict = None
        self.pseudo_dict = None

        self.objid_bri = None
        self.slitid_bri  = None
        self.snr_bar_bri = None
        self.use_weights = None
        self.wave_grid = None
        self.good_slits = None
        self.maskdef_offset = None


        # Load the stack_dict
        self.stack_dict = self.load_coadd2d_stacks(self.spec2d)
        self.pypeline = self.spectrograph.pypeline

        self.nexp = len(self.spec2d)

        # Check that there are the same number of slits on every exposure
        nslits_list = [slits.nslits for slits in self.stack_dict['slits_list']]
        if not len(set(nslits_list)) == 1:
            msgs.error('Not all of your exposures have the same number of slits. Check your inputs')
        # This is the number of slits of the single (un-coadded) frames
        self.nslits_single = nslits_list[0]

        # Check that nspec is the same for all the exposures
        nspec_list = [slits.nspec for slits in self.stack_dict['slits_list']]
        if not len(set(nspec_list)) == 1:
            msgs.error('Not all of your exposures have the same spectral dimension. Check your inputs')
        self.nspec = self.stack_dict['slits_list'][0].nspec

        # Check that binning is the same for all the exposures
        binspec_list = [slits.binspec for slits in self.stack_dict['slits_list']]
        binspat_list = [slits.binspat for slits in self.stack_dict['slits_list']]
        if not len(set(binspec_list)) == 1:
            msgs.error('Not all of your exposures have the same spectral binning. Check your inputs')
        if not len(set(binspat_list)) == 1:
            msgs.error('Not all of your exposures have the same spatial binning. Check your inputs')
        self.binning = np.array([self.stack_dict['slits_list'][0].binspec,
                                 self.stack_dict['slits_list'][0].binspat])

        self.spat_ids = self.stack_dict['slits_list'][0].spat_id

        # If smoothing is not input, smooth by 10% of the spectral dimension
        self.sn_smooth_npix = sn_smooth_npix if sn_smooth_npix is not None else 0.1*self.nspec

    def good_slitindx(self, only_slits=None):
        """
        This provides an array of index of slits in the un-coadded frames that are considered good for 2d coadding.
        A bitmask common to all the un-coadded frames is used to determine which slits are good. Also,
        If the `only_slits` parameter is provided only those slits are considered good for 2d coadding.
        
        Args:
            only_slits (:obj:`list`, optional):
                List of slits to combine. It must be `slitord_id`

        Returns:
            `numpy.ndarray`_: array of index of good slits in the un-coadded frames
        """

        only_slits = [only_slits] if (only_slits is not None and
                                      isinstance(only_slits, (int, np.integer))) else only_slits

        # This creates a unified bpm common to all frames
        slits0 = self.stack_dict['slits_list'][0]
        # bpm for the first frame
        reduce_bpm = (slits0.mask > 0) & (np.invert(slits0.bitmask.flagged(slits0.mask,
                                                                           flag=slits0.bitmask.exclude_for_reducing)))
        for i in range(1, self.nexp):
            # update bpm with the info from the other frames
            slits = self.stack_dict['slits_list'][i]
            reduce_bpm &= (slits.mask > 0) & (np.invert(slits.bitmask.flagged(slits.mask,
                                                                             flag=slits.bitmask.exclude_for_reducing)))
        # this are the good slit index according to the bpm mask
        good_slitindx = np.where(np.logical_not(reduce_bpm))[0]

        # If we want to coadd all the good slits
        if only_slits is None:
            return good_slitindx

        # If instead we want to coadd only a selected (by the user) number of slits
        # this are the `slitord_id` of the slits that we want to coadd
        only_slits = np.atleast_1d(only_slits)
        # create an array of slit index that are selected by the user and are also good slits
        good_onlyslits = np.array([], dtype=int)
        for islit in only_slits:
            if islit not in slits0.slitord_id[good_slitindx]:
                # Warnings for the slits that are selected by the user but NOT good slits
                msgs.warn('Slit {} cannot be coadd because masked'.format(islit))
            else:
                indx = np.where(slits0.slitord_id[good_slitindx] == islit)[0]
                good_onlyslits = np.append(good_onlyslits, good_slitindx[indx])
        return good_onlyslits

    def optimal_weights(self, slitorderid, objid, const_weights=False):
        """
        Determine optimal weights for 2d coadds. This script grabs the information from SpecObjs list for the
        object with specified slitid and objid and passes to coadd.sn_weights to determine the optimal weights for
        each exposure.

        Parameters
        ----------
        slitorderid : :obj:`int`
           The slit or order id that has the brightest object whose
           S/N will be used to determine the weight for each frame.
        objid : `numpy.ndarray`_
           Array of object indices with shape = (nexp,) of the
           brightest object whose S/N will be used to determine the
           weight for each frame.
        const_weights : :obj:`bool`
           Use constant weights for coadding the exposures.
           Default=False

        Returns
        -------
        rms_sn : ndarray, shape = (len(specobjs_list),)
            Root mean square S/N value for each input spectra
        weights : ndarray, shape (len(specobjs_list),)
            Weights to be applied to the spectra. These are
            signal-to-noise squared weights.
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

    def coadd(self, only_slits=None, interp_dspat=True):
        """
        Construct a 2d co-add of a stack of PypeIt spec2d reduction outputs.
        This method calls loops over slits/orders and performs the 2d-coadd by
        calling coadd.compute.coadd2d, which 'rectifies' images by coadding them
        about the reference_trace_stack.

        Parameters
        ----------
        only_slits : list, optional
           List of slits to operate on. Not currently supported, i.e. the code
           can currently only stack everything because the slit/reduction
           bitmask checking is not yet implemented. Default = None
        interp_dspat : bool, optional
           Interpolate in the spatial coordinate image to faciliate running
           through core.extract.local_skysub_extract.  Default=True


        Returns
        -------
        coadd_list : list
            List of dictionaries, one for each slit, containing the 2d stack.
            # TODO Make this a PypeIt object, with data model yada-yada.

        """
        # get slit index that indicates which slits are good for coadding
        self.good_slits = self.good_slitindx(only_slits=only_slits)
        # get the number of slits that are going to be coadded
        self.nslits_coadded = self.good_slits.size

        coadd_list = []
        for slit_idx in self.good_slits:
            msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(slit_idx, self.nslits_single - 1))
            ref_trace_stack = self.reference_trace_stack(slit_idx, offsets=self.offsets,
                                                         objid=self.objid_bri)
            thismask_stack = np.abs(self.stack_dict['slitmask_stack'] - self.stack_dict['slits_list'][0].spat_id[slit_idx]) <= self.par['coadd2d']['spat_toler']

            # maskdef info
            maskdef_dict = self.get_maskdef_dict(slit_idx, ref_trace_stack)

            # weights
            if not isinstance(self.use_weights, str) and self.use_weights.ndim > 2:
                weights = self.use_weights[slit_idx]
            else:
                weights = self.use_weights
            # Perform the 2d coadd
            coadd_dict = coadd.compute_coadd2d(ref_trace_stack, self.stack_dict['sciimg_stack'],
                                               self.stack_dict['sciivar_stack'],
                                               self.stack_dict['skymodel_stack'],
                                               self.stack_dict['mask_stack'] == 0,
                                               self.stack_dict['tilts_stack'],
                                               thismask_stack,
                                               self.stack_dict['waveimg_stack'],
                                               self.wave_grid, self.spat_samp_fact,
                                               maskdef_dict=maskdef_dict,
                                               weights=weights, interp_dspat=interp_dspat)
            coadd_list.append(coadd_dict)

        return coadd_list

    def create_pseudo_image(self, coadd_list):
        """
        ..todo.. see below

        THIS UNDOCUMENTED CODE PROBABLY SHOULD GENERATE AND RETURN
        STANDARD PYPEIT OBJCTS INSTEAD OF SOME UNDEFINED DICT"""

        # Check that self.nslit is equal to len(coadd_list)
        if self.nslits_coadded != len(coadd_list):
            msgs.error('Wrong number of slits for the 2d coadded frame')

        nspec_vec = np.zeros(self.nslits_coadded,dtype=int)
        nspat_vec = np.zeros(self.nslits_coadded,dtype=int)
        for islit, cdict in enumerate(coadd_list):
            nspec_vec[islit]=cdict['nspec']
            nspat_vec[islit]=cdict['nspat']

        # Determine the size of the pseudo image
        nspat_pad = 10
        nspec_pseudo = nspec_vec.max()
        nspat_pseudo = int(np.sum(nspat_vec) + (self.nslits_coadded + 1)*nspat_pad)  # Cast for SlitTraceSet
        spec_vec_pseudo = np.arange(nspec_pseudo)
        shape_pseudo = (nspec_pseudo, nspat_pseudo)
        imgminsky_pseudo = np.zeros(shape_pseudo)
        sciivar_pseudo = np.zeros(shape_pseudo)
        waveimg_pseudo = np.zeros(shape_pseudo)
        tilts_pseudo = np.zeros(shape_pseudo)
        spat_img_pseudo = np.zeros(shape_pseudo)
        nused_pseudo = np.zeros(shape_pseudo, dtype=int)
        inmask_pseudo = np.zeros(shape_pseudo, dtype=bool)
        wave_mid = np.zeros((nspec_pseudo, self.nslits_coadded))
        wave_mask = np.zeros((nspec_pseudo, self.nslits_coadded),dtype=bool)
        wave_min = np.zeros((nspec_pseudo, self.nslits_coadded))
        wave_max = np.zeros((nspec_pseudo, self.nslits_coadded))
        dspat_mid = np.zeros((nspat_pseudo, self.nslits_coadded))

        spat_left = nspat_pad
        slit_left = np.zeros((nspec_pseudo, self.nslits_coadded))
        slit_righ = np.zeros((nspec_pseudo, self.nslits_coadded))
        spec_min1 = np.zeros(self.nslits_coadded)
        spec_max1 = np.zeros(self.nslits_coadded)

        # maskdef info
        all_maskdef_ids = np.array([cc['maskdef_id'] for cc in coadd_list])
        if None not in all_maskdef_ids:
            maskdef_id = np.zeros(self.nslits_coadded, dtype=int)
            maskdef_objpos = np.zeros(self.nslits_coadded)
            maskdef_slitcen = np.zeros((nspec_pseudo, self.nslits_coadded))
            maskdef_designtab = Table()
        else:
            maskdef_id = None
            maskdef_objpos = None
            maskdef_slitcen = None
            maskdef_designtab = None

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
            image_temp = (coadd_dict['dspat'] - coadd_dict['dspat_mid'][0] + spat_left)*coadd_dict['outmask']
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

            # maskdef info
            if None not in all_maskdef_ids:
                maskdef_id[islit] = coadd_dict['maskdef_id']
                maskdef_objpos[islit] = coadd_dict['maskdef_objpos']
                maskdef_slitcen[:, islit] = np.full(nspec_pseudo, coadd_dict['maskdef_slitcen'])
                if coadd_dict['maskdef_designtab'] is not None:
                    maskdef_designtab = vstack([maskdef_designtab, coadd_dict['maskdef_designtab']])

        slits_pseudo \
                = slittrace.SlitTraceSet(slit_left, slit_righ, self.pypeline, detname=self.detname, nspat=nspat_pseudo,
                                         PYP_SPEC=self.spectrograph.name, specmin=spec_min1, specmax=spec_max1,
                                         maskdef_id=maskdef_id, maskdef_objpos=maskdef_objpos, maskdef_offset=0.,
                                         maskdef_slitcen=maskdef_slitcen, maskdef_designtab=maskdef_designtab)
                                         #master_key=self.stack_dict['master_key_dict']['trace'],
                                         #master_dir=self.master_dir)

        # change value of spat_id in maskdef_designtab
        # needs to be done here because spat_id is computed in slittrace
        if maskdef_designtab is not None:
            slits_pseudo.maskdef_designtab['SPAT_ID'] = slits_pseudo.spat_id

        # assign ech_order if exist
        slits_pseudo.ech_order = self.stack_dict['slits_list'][0].ech_order[self.good_slits] \
            if self.stack_dict['slits_list'][0].ech_order is not None else None
        slitmask_pseudo = slits_pseudo.slit_img()
        # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
        # which wreaks havoc on the local sky-subtraction
        min_slit_frac = 0.70
        spec_min = np.zeros(self.nslits_coadded)
        spec_max = np.zeros(self.nslits_coadded)
        for islit in range(self.nslits_coadded):
            spat_id = slits_pseudo.spat_id[islit]
            slit_width = np.sum(inmask_pseudo*(slitmask_pseudo == spat_id),axis=1)
            slit_width_img = np.outer(slit_width, np.ones(nspat_pseudo))
            med_slit_width = np.median(slit_width_img[slitmask_pseudo == spat_id])
            # TODO -- need inline docs
            nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
            nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
            slit_width_sm = ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
            igood = (slit_width_sm > min_slit_frac*med_slit_width)
            # TODO -- need inline docs
            spec_min[islit] = spec_vec_pseudo[igood].min()
            spec_max[islit] = spec_vec_pseudo[igood].max()
            bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_pseudo == spat_id)
            inmask_pseudo[bad_pix] = False

        # Update slits_pseudo
        slits_pseudo.specmin = spec_min
        slits_pseudo.specmax = spec_max

        return dict(nspec=nspec_pseudo, nspat=nspat_pseudo, imgminsky=imgminsky_pseudo,
                    sciivar=sciivar_pseudo, inmask=inmask_pseudo, tilts=tilts_pseudo,
                    waveimg=waveimg_pseudo, spat_img=spat_img_pseudo, slits=slits_pseudo,
                    wave_mask=wave_mask, wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max)

    def reduce(self, pseudo_dict, show=None, show_peaks=None, basename=None):
        """
        ..todo.. Please document me

        Args:
            pseudo_dict:
            show:
            show_peaks:

        Returns:

        """

        show = self.show if show is None else show
        show_peaks = self.show_peaks if show_peaks is None else show_peaks
        sciImage = pypeitimage.PypeItImage(image=pseudo_dict['imgminsky'],
                                           ivar=pseudo_dict['sciivar'],
                                           bpm=np.zeros_like(pseudo_dict['inmask'].astype(int)),  # Dummy bpm
                                           rn2img=np.zeros_like(pseudo_dict['inmask']).astype(float),  # Dummy rn2img
                                           crmask=np.invert(pseudo_dict['inmask'].astype(bool)))
        sciImage.detector = self.stack_dict['detectors'][0]
        #
        slitmask_pseudo = pseudo_dict['slits'].slit_img()
        sciImage.build_mask(slitmask=slitmask_pseudo)

        # Make changes to parset specific to 2d coadds
        parcopy = copy.deepcopy(self.par)
        parcopy['reduce']['findobj']['trace_npoly'] = 3        # Low order traces since we are rectified

        # Build the Calibrate object
        caliBrate = calibrations.Calibrations(None, self.par['calibrations'], self.spectrograph, None)
        caliBrate.slits = pseudo_dict['slits']
        caliBrate.det = self.det
        caliBrate.binning = self.binning

        # Manual extraction. Select only the manual extractions input by the user that match the current detector
        manual_dict = None
        if self.par['coadd2d']['manual'] is not None:
            manual_string = ','.join(np.atleast_1d(self.par['coadd2d']['manual']))
            if len(manual_string.strip()) > 0:
                manual_obj = ManualExtractionObj.by_fitstbl_input('None', manual_string)
                # manual_obj indx for this det
                thisdet = np.where(manual_obj.det == sciImage.detector.det)[0]
                manual_obj.det = manual_obj.det[thisdet]
                manual_obj.spec = manual_obj.spec[thisdet]
                manual_obj.spat = manual_obj.spat[thisdet]
                manual_obj.fwhm = manual_obj.fwhm[thisdet]
                if manual_obj.det.size > 0:
                    # TODO: Leaving `neg=False`, the default, consider changing to neg=self.find_negative.
                    manual_dict = manual_obj.dict_for_objfind(manual_obj.det[0])

        # Get bpm mask. There should not be any masked slits because we excluded those already
        # before the coadd, but we need to pass a bpm to FindObjects and Extract
        slits = pseudo_dict['slits']
        pseudo_reduce_bpm = (slits.mask > 0) & (np.invert(slits.bitmask.flagged(slits.mask,
                                                                         flag=slits.bitmask.exclude_for_reducing)))

        # Initiate FindObjects object
        objFind = find_objects.FindObjects.get_instance(sciImage, self.spectrograph, parcopy, caliBrate,
                                           'science_coadd2d', bkg_redux=self.bkg_redux,
                                           find_negative=self.find_negative, show=show)

        # Set the tilts and waveimg attributes from the psuedo_dict here, since we generate these dynamically from fits
        # normally, but this is not possible for coadds
        objFind.tilts = pseudo_dict['tilts']
        objFind.waveimg = pseudo_dict['waveimg']
        objFind.binning = self.binning
        objFind.basename = basename
        objFind.reduce_bpm = pseudo_reduce_bpm

        if show:
            gpm = sciImage.select_flag(invert=True)
            objFind.show('image', image=pseudo_dict['imgminsky']*gpm.astype(float),
                       chname='imgminsky', slits=True, clear=True)

        sobjs_obj, nobj = objFind.find_objects(sciImage.image, show_peaks=show_peaks,
                                               save_objfindQA=True, manual_extract_dict=manual_dict)

        # maskdef stuff
        if parcopy['reduce']['slitmask']['assign_obj'] and slits.maskdef_designtab is not None:
            # Get plate scale
            platescale = sciImage.detector.platescale * self.spat_samp_fact

            # Assign slitmask design information to detected objects
            slits.assign_maskinfo(sobjs_obj, platescale, None, TOLER=parcopy['reduce']['slitmask']['obj_toler'])

            if parcopy['reduce']['slitmask']['extract_missing_objs'] is True:
                # Set the FWHM for the extraction of missing objects
                fwhm = slits.get_maskdef_extract_fwhm(sobjs_obj, platescale,
                                                      parcopy['reduce']['slitmask']['missing_objs_fwhm'],
                                                      parcopy['reduce']['findobj']['find_fwhm'])
                # Assign undetected objects
                sobjs_obj = slits.mask_add_missing_obj(sobjs_obj, None, fwhm,
                                                       parcopy['reduce']['slitmask']['missing_objs_boxcar_rad']/platescale)

        # Initiate Extract object
        exTract = extraction.Extract.get_instance(sciImage, sobjs_obj, self.spectrograph, parcopy, caliBrate,
                                                  'science_coadd2d', bkg_redux=self.bkg_redux, show=show)

        # Set the tilts and waveimg attributes from the psuedo_dict here, since we generate these dynamically from fits
        # normally, but this is not possible for coadds
        exTract.tilts = pseudo_dict['tilts']
        exTract.waveimg = pseudo_dict['waveimg']
        exTract.binning = self.binning
        exTract.basename = basename
        exTract.reduce_bpm = pseudo_reduce_bpm

        # Local sky-subtraction
        global_sky_pseudo = np.zeros_like(pseudo_dict['imgminsky']) # No global sky for co-adds since we go straight to local
        skymodel_pseudo, objmodel_pseudo, \
            ivarmodel_pseudo, outmask_pseudo, sobjs = exTract.local_skysub_extract(global_sky_pseudo, sobjs_obj,
                                                                                   spat_pix=pseudo_dict['spat_img'],
                                                                                   model_noise=False, show_profile=show,
                                                                                   show=show)

        if self.find_negative:
            sobjs.purge_neg()

        # Add the rest to the pseudo_dict
        pseudo_dict['skymodel'] = skymodel_pseudo
        pseudo_dict['objmodel'] = objmodel_pseudo
        pseudo_dict['ivarmodel'] = ivarmodel_pseudo
        pseudo_dict['outmask'] = outmask_pseudo
        pseudo_dict['sobjs'] = sobjs
        self.pseudo_dict=pseudo_dict

        return pseudo_dict['imgminsky'], pseudo_dict['sciivar'], skymodel_pseudo, \
               objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs, sciImage.detector, slits, \
               pseudo_dict['tilts'], pseudo_dict['waveimg']



#    def save_masters(self):
#
#        # Write out the pseudo master files to disk
#        master_key_dict = self.stack_dict['master_key_dict']
#
#        # TODO: These saving operations are a temporary kludge
#        # spectrograph is needed for header
#        waveImage = WaveImage(self.pseudo_dict['waveimg'], PYP_SPEC=self.spectrograph.spectrograph)
#        wave_filename = masterframe.construct_file_name(WaveImage, master_key_dict['arc'], self.master_dir)
#        waveImage.to_master_file(wave_filename)
#
#        # TODO: Assumes overwrite=True
#        slit_filename = masterframe.construct_file_name(self.pseudo_dict['slits'], master_key_dict['trace'], self.master_dir)
#        self.pseudo_dict['slits'].to_master_file(slit_filename) #self.master_dir, master_key_dict['trace'], self.spectrograph.spectrograph)
#    '''

    def snr_report(self, snr_bar, slitid=None):
        """
        ..todo.. I need a doc string

        Args:
            snr_bar:
            slitid:

        Returns:

        """

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

    def offsets_report(self, offsets, offsets_method):
        """
        Print out a report on the offsets

        Args:
            offsets:
            offsets_method:

        Returns:

        """

        if offsets_method is not None and offsets is not None:
            msg_string = msgs.newline() + '---------------------------------------------'
            msg_string += msgs.newline() + ' Summary of offsets from {}     '.format(offsets_method)
            msg_string += msgs.newline() + '---------------------------------------------'
            msg_string += msgs.newline() + '           exp#      offset                  '
            for iexp, off in enumerate(offsets):
                msg_string += msgs.newline() + '            {:d}        {:5.2f}'.format(iexp, off)
            msg_string += msgs.newline() + '-----------------------------------------------'
            msgs.info(msg_string)

    def offset_slit_cen(self, slitid, offsets):
        """
        ..todo.. I need a doc string

        Args:
            slitid:
            offsets:

        Returns:

        """
        ref_trace_stack = np.zeros((self.stack_dict['slits_list'][0].nspec, len(offsets)),
                                   dtype=float)
        for iexp, slits in enumerate(self.stack_dict['slits_list']):
            ref_trace_stack[:,iexp] = slits.center[:,slitid] - offsets[iexp]
        return ref_trace_stack

    def get_wave_grid(self, **kwargs_wave):
        """
        Routine to create a wavelength grid for 2d coadds using all of the
        wavelengths of the extracted objects. Calls
        :func:`~pypeit.core.wavecal.wvutils.get_wave_grid`.

        Args:
            **kwargs_wave (dict):
                Optional argumments for
                :func:`~pypeit.core.wavecal.wvutils.get_wave_grid`.

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
        nobjs_tot = int(np.array([len(spec) for spec in self.stack_dict['specobjs_list']]).sum())
        # TODO: Do we need this flag since we can determine whether or not we have specobjs from nobjs_tot?
        #  This all seems a bit hacky
        if self.par['coadd2d']['use_slits4wvgrid'] or nobjs_tot==0:
            nslits_tot = np.sum([slits.nslits for slits in self.stack_dict['slits_list']])
            waves = np.zeros((self.nspec, nslits_tot*3))
            gpm = np.zeros_like(waves, dtype=bool)
            box_radius = 3.
            indx = 0
            # Loop on the exposures
            for waveimg, slitmask, slits in zip(self.stack_dict['waveimg_stack'],
                                                self.stack_dict['slitmask_stack'],
                                                self.stack_dict['slits_list']):
                slits_left, slits_righ, _ = slits.select_edges()
                row = np.arange(slits_left.shape[0])
                # Loop on the slits
                for kk, spat_id in enumerate(slits.spat_id):
                    mask = slitmask == spat_id
                    # Create apertures at 5%, 50%, and 95% of the slit width to cover full range of wavelengths
                    # on this slit
                    trace_spat = slits_left[:, kk][:,np.newaxis] +  np.outer((slits_righ[:,kk] - slits_left[:,kk]),[0.05,0.5,0.95])
                    box_denom = moment1d(waveimg * mask > 0.0, trace_spat, 2 * box_radius, row=row)[0]
                    wave_box = moment1d(waveimg * mask, trace_spat, 2 * box_radius,
                                    row=row)[0] / (box_denom + (box_denom == 0.0))
                    waves[:, indx:indx+3] = wave_box
                    # TODO -- This looks a bit risky
                    gpm[:, indx: indx+3] = wave_box > 0.
                    indx += 3
        else:
            waves = np.zeros((self.nspec, nobjs_tot))
            gpm = np.zeros_like(waves, dtype=bool)
            indx = 0
            for spec_this in self.stack_dict['specobjs_list']:
                for spec in spec_this:
                    waves[:, indx] = spec.OPT_WAVE
                    # TODO -- OPT_MASK is likely to become a bpm with int values
                    gpm[:, indx] = spec.OPT_MASK
                    indx += 1

        wave_grid, wave_grid_mid, dsamp = wvutils.get_wave_grid(waves, masks=gpm,
                                                                spec_samp_fact=self.spec_samp_fact,
                                                                **kwargs_wave)
        return wave_grid, wave_grid_mid, dsamp

    def load_coadd2d_stacks(self, spec2d, chk_version=False):
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

        # Get the master dir
        redux_path = os.getcwd()

        # Grab the files
        #head2d_list = []
        specobjs_list = []
        slits_list = []
        nfiles =len(spec2d)
        detectors_list = []
        maskdef_designtab_list = []
        spat_flexure_list = []
        for ifile, f in enumerate(spec2d):
            if isinstance(f, spec2dobj.Spec2DObj):
                # If spec2d is a list of objects
                s2dobj = f
            else:
                # If spec2d is a list of files, option to also use spec1ds
                s2dobj = spec2dobj.Spec2DObj.from_file(f, self.detname, chk_version=chk_version)
                spec1d_file = f.replace('spec2d', 'spec1d')
                if os.path.isfile(spec1d_file):
                    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file, chk_version=chk_version)
                    this_det = sobjs.DET == self.detname
                    specobjs_list.append(sobjs[this_det])
            # TODO the code should run without a spec1d file, but we need to implement that
            slits_list.append(s2dobj.slits)
            detectors_list.append(s2dobj.detector)
            maskdef_designtab_list.append(s2dobj.maskdef_designtab)
            spat_flexure_list.append(s2dobj.sci_spat_flexure)

            if ifile == 0:
                sciimg_stack = np.zeros((nfiles,) + s2dobj.sciimg.shape, dtype=float)
                waveimg_stack = np.zeros_like(sciimg_stack, dtype=float)
                tilts_stack = np.zeros_like(sciimg_stack, dtype=float)
                skymodel_stack = np.zeros_like(sciimg_stack, dtype=float)
                sciivar_stack = np.zeros_like(sciimg_stack, dtype=float)
                mask_stack = np.zeros_like(sciimg_stack, dtype=float)
                slitmask_stack = np.zeros_like(sciimg_stack, dtype=int)

            sciimg_stack[ifile, :, :] = s2dobj.sciimg
            waveimg_stack[ifile, :, :] = s2dobj.waveimg
            skymodel_stack[ifile, :, :] = s2dobj.skymodel
            sciivar_stack[ifile, :, :] = s2dobj.ivarmodel
            mask_stack[ifile, :, :] = s2dobj.bpmmask
            # TODO -- Set back after done testing

            slitmask_stack[ifile, :, :] = s2dobj.slits.slit_img(flexure=s2dobj.sci_spat_flexure)
            #slitmask_stack[ifile, :, :] = spec2DObj.slits.slit_img(flexure=0.)
            _spat_flexure = 0. if s2dobj.sci_spat_flexure is None else s2dobj.sci_spat_flexure
            #_tilt_flexure_shift = _spat_flexure - spec2DObj.tilts.spat_flexure if spec2DObj.tilts.spat_flexure is not None else _spat_flexure
            tilts_stack[ifile,:,:] = s2dobj.tilts #.fit2tiltimg(slitmask_stack[ifile, :, :], flexure=_tilt_flexure_shift)

        return dict(specobjs_list=specobjs_list, slits_list=slits_list,
                    slitmask_stack=slitmask_stack,
                    sciimg_stack=sciimg_stack, sciivar_stack=sciivar_stack,
                    skymodel_stack=skymodel_stack, mask_stack=mask_stack,
                    tilts_stack=tilts_stack, waveimg_stack=waveimg_stack,
                    redux_path=redux_path,
                    detectors=detectors_list,
                    spectrograph=self.spectrograph.name,
                    pypeline=self.spectrograph.pypeline,
                    maskdef_designtab_list=maskdef_designtab_list,
                    spat_flexure_list=spat_flexure_list)

    def check_input(self, input, type='weights'):
        """
        Check that the number of input values (weights or offsets) is the same as the number of exposures
        Args:
            input (:obj:`list` or `numpy.ndarray`_): User input values (e.g., weights or offsets)
            type (:obj:`str`): String defining what the quantities are

        Returns:
            :obj:`list` or `numpy.ndarray`_: User input values
        """
        if isinstance(input, (list, np.ndarray)):
            if len(input) != self.nexp:
                msgs.error(f'If {type} are input it must be a list/array with same number of elements as exposures')
            return np.atleast_1d(input)
        msgs.error(f'Unrecognized format for {type}')

    def compute_offsets(self, offsets):
        """
         Dummy method to compute offsets and weights. Overloaded by child methods.

        Determine self.offsets, the offset of the frames to be coadded with respect to the first frame

        Args:
            offsets (:obj:`list` or :obj:`str`):
                Value that guides the determination of the offsets.
                It could be a list of offsets, or a string, or None.

        """

        pass

    def compute_weights(self, weights):
        """
         Dummy method to compute weights. Overloaded by child methods.

        Determine self.use_weights, the weights to be used in the coadd2d

        Args:
            weights (:obj:`list` or :obj:`str`):
                Value that guides the determination of the weights.
                It could be a list of weights or a string. If 'auto' the weight will be computed using
                the brightest trace, if 'uniform' uniform weights will be used.

        """

        pass

    def get_brightest_object(self, specobjs_list, spat_ids):
        """
        Dummy method to identify the brightest object. Overloaded by child methods.

        Parameters
        ----------
        specobjs_list
        spat_ids

        Returns
        -------

        """
        pass

    def reference_trace_stack(self, slitid, offsets=None, objid=None):
        """
        Dummy method to obtain the stack of reference traces. Overloaded by child methods.

        Args:
            slitid:
            offsets:
            objid:

        Returns:

        """
        pass

    def get_maskdef_dict(self, slit_idx, ref_trace_stack):
        """
        Dummy method to get maskdef info. Overloaded by child methods.

        Args:
            slit_idx:
            ref_trace_stack:

        Returns:

        """

        return dict(maskdef_id=None, maskdef_objpos=None, maskdef_designtab=None)


# Multislit can coadd with:
# 1) input offsets or if offsets is None, it will find the brightest trace and compute them
# 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object

# Echelle can either stack with:
# 1) input offsets or if offsets is None, it will find the objid of brightest trace and stack all orders relative to the trace of this object.
# 2) specified weights, or if weights is None and auto_weights=True,
#    it will use wavelength dependent weights determined from the spectrum of the brightest objects objid on each order

class MultiSlitCoAdd2D(CoAdd2D):
    """
    Child of Coadd2d for Multislit and Longslit reductions. For documentation see CoAdd2d parent class above.

        # Multislit can coadd with:
        # 1) input offsets or if offsets is None, it will find the brightest trace and compute them
        # 2) specified weights, or if weights is None and auto_weights=True, it will compute weights using the brightest object


    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto',
                 spec_samp_fact=1.0, spat_samp_fact=1.0, sn_smooth_npix=None,
                 bkg_redux=False, find_negative=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(MultiSlitCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                               spec_samp_fact=spec_samp_fact, spat_samp_fact=spat_samp_fact,
                                               sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative, par=par,
                                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                                        debug=debug, **kwargs_wave)

        # maskdef offset
        self.maskdef_offset = np.array([slits.maskdef_offset for slits in self.stack_dict['slits_list']])

        # Default wave_method for Multislit is linear
        kwargs_wave['wave_method'] = 'linear' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        # find if there is a bright object we could use
        if len(self.stack_dict['specobjs_list']) > 0:
            self.objid_bri, self.slitidx_bri, self.spatid_bri, self.snr_bar_bri = \
                self.get_brightest_obj(self.stack_dict['specobjs_list'], self.spat_ids)
        else:
            self.objid_bri, self.slitidx_bri, self.spatid_bri, self.snr_bar_bri = None, None, None, None

        # get self.offsets
        self.compute_offsets(offsets)

        # get self.use_weights
        self.compute_weights(weights)

    # TODO When we run multislit, we actually compute the rebinned images twice. Once here to compute the offsets
    # and another time to weighted_combine the images in compute2d. This could be sped up
    # TODO The reason we rebin the images for the purposes of computing the offsets is to deal with combining
    # data that are dithered in the spectral direction. In these situations you cannot register the two dithered
    # reference objects into the same frame without first rebinning them onto the same grid.

    def compute_offsets(self, offsets):
        """
        Determine self.offsets, the offset of the frames to be coadded with respect to the first frame

        Args:
            offsets (:obj:`list` or :obj:`str`):
                Value that guides the determination of the offsets.
                It could be a list of offsets, or a string, or None. If equal to 'maskdef_offsets' the
                offsets computed during the slitmask design matching will be used.
        """

        # If offsets are not provided by the user and there is no a bright object
        # to be used to compute the offsets, we cannot continue
        if (self.objid_bri is None) and (offsets is None):
            msgs.error('Offsets cannot be computed because no unique reference object '
                       'with the highest S/N was found. To continue, provide offsets in `Coadd2DPar`')

        offsets_method = None
        # 1) parset `offsets` is = 'maskdef_offsets' (no matter if we have a bright object or not)
        if offsets == 'maskdef_offsets':
            if self.maskdef_offset is not None:
                # the offsets computed during the main reduction (`run_pypeit`) are used
                msgs.info('Determining offsets using maskdef_offset recoded in SlitTraceSet')
                self.offsets = self.maskdef_offset[0] - self.maskdef_offset
                offsets_method = 'maskdef_offset'
            else:
                # if maskdef_offsets were not computed during the main reduction, we cannot continue
                msgs.error('No maskdef_offset recoded in SlitTraceSet')
        # 2) a list of offsets is provided by the user (no matter if we have a bright object or not)
        elif offsets is not None:
            msgs.info('Using user input offsets')
            # use them
            self.offsets = self.check_input(offsets, type='offsets')
            offsets_method = 'user input'
        # 3) parset `offsets` is None but we have a bright object
        else:
            # Compute offsets using the bright object
            offsets_method = 'highest S/N object found on slitid = {:d}'.format(self.spatid_bri)
            msgs.info('Determining offsets using brightest object on slit: {:d} with avg SNR={:5.2f}'.format(self.spatid_bri,np.mean(self.snr_bar_bri)))
            thismask_stack = np.abs(self.stack_dict['slitmask_stack'] - self.spatid_bri) <= self.par['coadd2d']['spat_toler']
            trace_stack_bri = np.zeros((self.nspec, self.nexp))
            # TODO Need to think abbout whether we have multiple tslits_dict for each exposure or a single one
            for iexp in range(self.nexp):
                trace_stack_bri[:,iexp] = self.stack_dict['slits_list'][iexp].center[:,self.slitidx_bri]
            #            trace_stack_bri[:,iexp] = (self.stack_dict['tslits_dict_list'][iexp]['slit_left'][:,slitid_bri] +
            #                                       self.stack_dict['tslits_dict_list'][iexp]['slit_righ'][:,slitid_bri])/2.0
            # Determine the wavelength grid that we will use for the current slit/order

            ## TODO: Should the spatial and spectral samp_facts here match those of the final coadded data, or she would
            ## compute offsets at full resolution??
            wave_bins = coadd.get_wave_bins(thismask_stack, self.stack_dict['waveimg_stack'], self.wave_grid)
            dspat_bins, dspat_stack = coadd.get_spat_bins(thismask_stack, trace_stack_bri)

            sci_list = [self.stack_dict['sciimg_stack'] - self.stack_dict['skymodel_stack']]
            var_list = []

            msgs.info('Rebinning Images')
            sci_list_rebin, var_list_rebin, norm_rebin_stack, nsmp_rebin_stack = coadd.rebin2d(
                wave_bins, dspat_bins, self.stack_dict['waveimg_stack'], dspat_stack, thismask_stack,
                (self.stack_dict['mask_stack'] == 0), sci_list, var_list)
            thismask = np.ones_like(sci_list_rebin[0][0,:,:],dtype=bool)
            nspec_pseudo, nspat_pseudo = thismask.shape
            slit_left = np.full(nspec_pseudo, 0.0)
            slit_righ = np.full(nspec_pseudo, nspat_pseudo)
            inmask = norm_rebin_stack > 0
            traces_rect = np.zeros((nspec_pseudo, self.nexp))
            sobjs = specobjs.SpecObjs()
            for iexp in range(self.nexp):
                sobjs_exp = findobj_skymask.objs_in_slit(sci_list_rebin[0][iexp,:,:], thismask, slit_left, slit_righ,
                                               inmask=inmask[iexp,:,:], has_negative=self.find_negative,
                                               fwhm=self.par['reduce']['findobj']['find_fwhm'],
                                               trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                                               npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
                                               maxdev=self.par['reduce']['findobj']['find_maxdev'],
                                               ncoeff=3, sig_thresh=self.par['reduce']['findobj']['sig_thresh'], nperslit=1,
                                               find_min_max=self.par['reduce']['findobj']['find_min_max'],
                                               show_trace=self.debug_offsets, show_peaks=self.debug_offsets)
                sobjs.add_sobj(sobjs_exp)
                traces_rect[:, iexp] = sobjs_exp.TRACE_SPAT
            # Now deterimine the offsets. Arbitrarily set the zeroth trace to the reference
            med_traces_rect = np.median(traces_rect,axis=0)
            offsets = med_traces_rect[0] - med_traces_rect
            # TODO create a QA with this
            if self.debug_offsets:
                for iexp in range(self.nexp):
                    plt.plot(traces_rect[:, iexp], linestyle='--', label='original trace')
                    plt.plot(traces_rect[:, iexp] + offsets[iexp], label='shifted traces')
                    plt.legend()
                plt.show()

            self.offsets = offsets

        self.offsets_report(self.offsets, offsets_method)

    def compute_weights(self, weights):
        """
        Determine self.use_weights, the weights to be used in the coadd2d

        Args:
            weights (:obj:`list` or :obj:`str`):
                Value that guides the determination of the weights.
                It could be a list of weights or a string. If equal to 'auto', the weight will be computed
                using the brightest trace, if 'uniform' uniform weights will be used.

        """

        # 1) No bright object and parset `weights` is 'auto' or 'uniform',
        # or Yes bright object but the user wants still to use uniform weights
        if ((self.objid_bri is None) and (weights in ['auto', 'uniform'])) or \
                ((self.objid_bri is not None) and (weights == 'uniform')):
            # use uniform weights
            self.use_weights = 'uniform'
            if weights == 'auto':
                # warn if the user had put `auto` in the parset
                msgs.warn('Weights cannot be computed because no unique reference object '
                          'with the highest S/N was found. Using uniform weights instead.')
            elif weights == 'uniform':
                msgs.info('Using uniform weights')
        # 2) Yes bright object and parset `weights` is equal to 'auto'
        elif (self.objid_bri is not None) and (weights == 'auto'):
            # compute weights using bright object
            _, self.use_weights = self.optimal_weights(self.spatid_bri, self.objid_bri, const_weights=True)
            msgs.info('Computing weights using a unique reference object with the highest S/N')
            self.snr_report(self.snr_bar_bri, slitid=self.spatid_bri)
        # 3) User inputs weights
        else:
            # use those inputs
            self.use_weights = self.check_input(weights, type='weights')
            msgs.info('Using user input weights')

    def get_brightest_obj(self, specobjs_list, spat_ids):

        """
        Utility routine to find the brightest object in each exposure given a specobjs_list for MultiSlit reductions.

        Args:
            specobjs_list: list
               List of SpecObjs objects.
            spat_ids (`numpy.ndarray`_):

        Returns:
            tuple: Returns the following:
                - objid: ndarray, int, shape (len(specobjs_list),):
                  Array of object ids representing the brightest object
                  in each exposure
                - slit_idx (int): 0-based index
                - spat_id (int): SPAT_ID for slit that highest S/N ratio object is on
                  (only for pypeline=MultiSlit)
                - snr_bar: ndarray, float, shape (len(list),): Average
                  S/N over all the orders for this object
        """
        nexp = len(specobjs_list)
        nspec = specobjs_list[0][0].TRACE_SPAT.shape[0]
        nslits = spat_ids.size

        slit_snr_max = np.full((nslits, nexp), -np.inf)
        objid_max = np.zeros((nslits, nexp), dtype=int)
        # Loop over each exposure, slit, find the brighest object on that slit for every exposure
        for iexp, sobjs in enumerate(specobjs_list):
            msgs.info("Working on exposure {}".format(iexp))
            for islit, spat_id in enumerate(spat_ids):
                ithis = np.abs(sobjs.SLITID - spat_id) <= self.par['coadd2d']['spat_toler']
                nobj_slit = np.sum(ithis)
                # ignore objects for which optimal extraction could not be performed (i.e., OPT_COUNTS is None)
                if np.any(ithis) and np.any(sobjs[ithis].OPT_COUNTS):
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
            msgs.warn('You do not appear to have a unique reference object that was traced as the highest S/N '
                       'ratio on the same slit of every exposure')
            return None, None, None, None
        return objid, slitid, spat_ids[slitid], snr_bar

    # TODO add an option here to actually use the reference trace for cases where they are on the same slit and it is
    # single slit???
    def reference_trace_stack(self, slitid, offsets=None, objid=None):
        """
        ..todo..  I need a doc string

        Args:
            slitid:
            offsets:
            objid:

        Returns:

        """

        return self.offset_slit_cen(slitid, offsets)

    def get_maskdef_dict(self, slit_idx, ref_trace_stack):
        """

        Args:
            slit_idx (:obj:`int`): index of a slit in the uncoadded frames
            ref_trace_stack(`numpy.ndarray`_): Stack of reference traces about
             which the images are rectified and coadded.  It is the slitcen appropriately
            shifted according the frames offsets. Shape is (nspec, nimgs).

        Returns:
            :obj:`dict`: Dictionary containing all the maskdef info. The quantities saved
            are: maskdef_id, maskdef_objpos, maskdef_slitcen, maskdef_designtab. To learn what
            they are see :class:`~pypeit.slittrace.SlitTraceSet` datamodel.

        """
        # maskdef info
        if self.par['calibrations']['slitedges']['use_maskdesign'] and \
                self.stack_dict['slits_list'][0].maskdef_id is not None and \
                self.stack_dict['slits_list'][0].maskdef_objpos is not None and \
                self.stack_dict['maskdef_designtab_list'][0] is not None:
            # maskdef_designtab info for only this slit
            this_idx = self.stack_dict['maskdef_designtab_list'][0]['SPAT_ID'] == self.stack_dict['slits_list'][0].spat_id[slit_idx]
            this_maskdef_designtab = self.stack_dict['maskdef_designtab_list'][0][this_idx]
            # remove columns that that are irrelevant in the coadd2d frames
            this_maskdef_designtab.remove_columns(['TRACEID', 'TRACESROW', 'TRACELPIX', 'TRACERPIX',
                                                   'SLITLMASKDEF', 'SLITRMASKDEF'])
            this_maskdef_designtab.meta['MASKRMSL'] = 0.
            this_maskdef_designtab.meta['MASKRMSR'] = 0.

            # maskdef_id for this slit
            imaskdef_id = self.stack_dict['slits_list'][0].maskdef_id[slit_idx]

            # maskdef_slitcenters. This trace the slit center along the spectral direction.
            # But here we take only the value at the mid point
            maskdef_slitcen_pixpos = self.stack_dict['slits_list'][0].maskdef_slitcen[self.nspec//2, slit_idx] + self.maskdef_offset
            # binned maskdef_slitcenters position with respect to the center of the slit in ref_trace_stack
            # this value should be the same for each exposure, but in case there are differences we take the mean value
            imaskdef_slitcen_dspat = np.mean((maskdef_slitcen_pixpos - ref_trace_stack[self.nspec//2, :])/self.spat_samp_fact)

            # expected position of the targeted object from slitmask design (as distance from left slit edge)
            imaskdef_objpos = self.stack_dict['slits_list'][0].maskdef_objpos[slit_idx]

            # find left edge
            slits_left, _, _ = self.stack_dict['slits_list'][0].select_edges(flexure=self.stack_dict['spat_flexure_list'][0])
            # targeted object spat pix
            maskdef_obj_pixpos = imaskdef_objpos + self.maskdef_offset + slits_left[self.nspec//2, slit_idx]
            # binned expected object position with respect to the center of the slit in ref_trace_stack
            # this value should be the same for each exposure, but in case there are differences we take the mean value
            imaskdef_objpos_dspat = np.mean((maskdef_obj_pixpos - ref_trace_stack[self.nspec//2, :])/self.spat_samp_fact)

        else:
            this_maskdef_designtab = None
            imaskdef_id = None
            imaskdef_slitcen_dspat = None
            imaskdef_objpos_dspat = None

        return dict(maskdef_id=imaskdef_id, maskdef_objpos=imaskdef_objpos_dspat,
                    maskdef_slitcen=imaskdef_slitcen_dspat, maskdef_designtab=this_maskdef_designtab)


class EchelleCoAdd2D(CoAdd2D):
    """
    Coadd Echelle reductions.
    
    For documentation see :class:`CoAdd2D`.

    Echelle can either stack with:

        - input ``offsets`` or if ``offsets`` is None, it will find
          the ``objid`` of brightest trace and stack all orders
          relative to the trace of this object.

        - specified ``weights``, or if ``weights`` is None and
          ``auto_weights`` is True, it will use wavelength dependent
          weights determined from the spectrum of the brightest
          objects ``objid`` on each order

    """
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto',
                 spec_samp_fact=1.0, spat_samp_fact=1.0, sn_smooth_npix=None,
                 bkg_redux=False, find_negative=False, show=False, show_peaks=False, debug_offsets=False, debug=False,
                 **kwargs_wave):
        super(EchelleCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                             spec_samp_fact=spec_samp_fact, spat_samp_fact=spat_samp_fact,
                                             sn_smooth_npix=sn_smooth_npix, bkg_redux=bkg_redux, find_negative=find_negative,
                                             par=par, show=show, show_peaks=show_peaks, debug_offsets=debug_offsets,
                                             debug=debug, **kwargs_wave)

        # Default wave_method for Echelle is log10
        kwargs_wave['wave_method'] = 'log10' if 'wave_method' not in kwargs_wave else kwargs_wave['wave_method']
        self.wave_grid, self.wave_grid_mid, self.dsamp = self.get_wave_grid(**kwargs_wave)

        # find if there is a bright object we could use
        if len(self.stack_dict['specobjs_list']) > 0:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri = \
                self.get_brightest_obj(self.stack_dict['specobjs_list'], self.nslits_single)
        else:
            self.objid_bri, self.slitid_bri, self.snr_bar_bri = None, None, None

        # get self.offsets
        self.compute_offsets(offsets)

        # get self.use_weights
        self.compute_weights(weights)

    # TODO move the part of the compute_offsets compute_weights common to MultiSlit and Echelle to the parent class
    def compute_offsets(self, offsets):
        """
        Determine self.offsets, the offset of the frames to be coadded with respect to the first frame

        Args:
            offsets (:obj:`list` or :obj:`str`):
                Value that guides the determination of the offsets.
                It could be a list of offsets, or a string, or None.

        """

        # If offsets are not provided by the user and there is no a bright object
        # to be used to compute the offsets, we cannot continue
        if (self.objid_bri is None) and (offsets is None):
            msgs.error('Offsets cannot be computed because no unique reference object '
                       'with the highest S/N was found. To continue, provide offsets in `Coadd2DPar`')

        offsets_method = None
        # 1) a list of offsets is provided by the user (no matter if we have a bright object or not)
        if offsets is not None:
            msgs.info('Using user input offsets')
            # use them
            self.offsets = self.check_input(offsets, type='offsets')
            # if the user inputs the offsets, `objid_bri` cannot be used to compute the reference stack,
            # so it needs to be None
            self.objid_bri = None
            offsets_method = 'user input'
        # 2) parset `offsets` is None but we have a bright object
        else:
            # offsets are not determined, but the bright object is used to construct
            # a reference trace (this is done in coadd using method `reference_trace_stack`)
            msgs.info('Reference trace about which 2d coadd is performed is computed using the brightest object')
            self.offsets = None

        self.offsets_report(self.offsets, offsets_method)

    def compute_weights(self, weights):
        """
        Determine self.use_weights, the weights to be used in the coadd2d

        Args:
            weights (:obj:`list` or :obj:`str`):
                Value that guides the determination of the weights.
                It could be a list of weights or a string. If 'auto' the weight will be computed using
                the brightest trace, if 'uniform' uniform weights will be used.

        """

        # 1) No bright object and parset `weights` is 'auto' or 'uniform',
        # or Yes bright object but the user wants still to use uniform weights
        if ((self.objid_bri is None) and (weights in ['auto', 'uniform'])) or \
                ((self.objid_bri is not None) and (weights == 'uniform')):
            # use uniform weights
            self.use_weights = 'uniform'
            if weights == 'auto':
                # warn if the user had put `auto` in the parset
                msgs.warn('Weights cannot be computed because no unique reference object '
                          'with the highest S/N was found. Using uniform weights instead.')
            elif weights == 'uniform':
                msgs.info('Using uniform weights')
        # 2) Yes bright object and parset `weights` is equal to 'auto'
        elif (self.objid_bri is not None) and (weights == 'auto'):
            # computing a list of weights for all the slitord_ids that we than parse in coadd
            slitord_ids = self.stack_dict['slits_list'][0].slitord_id
            use_weights = []
            for id in slitord_ids:
                _, iweights = self.optimal_weights(id, self.objid_bri)
                use_weights.append(iweights)
            self.use_weights = np.array(use_weights)
            msgs.info('Computing weights using a unique reference object with the highest S/N')
            self.snr_report(self.snr_bar_bri)
        # 3) User inputs weights
        else:
            # use those inputs
            self.use_weights = self.check_input(weights, type='weights')
            msgs.info('Using user input weights')

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
                    # ignore objects for which optimal extraction could not be performed (i.e., OPT_COUNTS is None)
                    if np.any(sobjs[ind][0].OPT_COUNTS):
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
        if 0 in snr_bar:
            msgs.warn('You do not appear to have a unique reference object that was traced as the highest S/N '
                      'ratio for every exposure')
            return None, None, None
        return objid, None, snr_bar

    def reference_trace_stack(self, slitid, offsets=None, objid=None):
        """
        Utility function for determining the reference trace about
        which 2d coadds are performed.

        There are two modes of operation to determine the reference
        trace for the 2d coadd of a given slit/order:

            #. ``offsets``: We stack about the center of the slit for
               the slit in question with the input offsets added

            #. ``ojbid``: We stack about the trace of a reference
               object for this slit given for each exposure by the
               input objid

        Either offsets or objid must be provided, but the code will
        raise an exception if both are provided.

        Args:
            slitid (int):
                The slit or order that we are currently considering
            stack_dict (dict):
                Dictionary containing all the images and keys
                required for performing 2d coadds.
            offsets (list, `numpy.ndarray`_):
                An array of offsets with the same dimensionality as
                the nexp, the numer of images being coadded.
            objid (list, `numpy.ndarray`_):
                An array of objids with the same dimensionality as
                the nexp, the number of images being coadded.

        Returns:
            `numpy.ndarray`: An array with shape (nspec, nexp)
            containing the reference trace for each of the nexp
            exposures.

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

