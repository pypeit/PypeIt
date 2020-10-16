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

from pypeit import msgs
from pypeit import specobjs
from pypeit import slittrace
from pypeit import reduce
from pypeit.images import pypeitimage
from pypeit.core import extract
from pypeit.core import coadd
from pypeit.core import parse
from pypeit import calibrations
from pypeit import spec2dobj
from pypeit.core.moment import moment1d


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
                     sn_smooth_npix=None, ir_redux=False, show=False,
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
                        sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux,
                        show=show, show_peaks=show_peaks, debug_offsets=debug_offsets, debug=debug,
                        **kwargs_wave)

    def __init__(self, spec2d, spectrograph, par, det=1, offsets=None, weights='auto',
                 sn_smooth_npix=None, ir_redux=False, show=False,
                 show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        """

        Parameters
        ----------
            spec2d_files (list):
               List of spec2d files or a list of Spec2dObj objects
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
                The instrument used to collect the data to be reduced.
            par (:class:`pypeit.par.parset.ParSet`):
                Parset object
            det (int): optional
                Detector to reduce
            offsets (ndarray): default=None
                Spatial offsets to be applied to each image before coadding. For the default mode of None, images
                are registered automatically using the trace of the brightest object. Input offsets are not yet supported.
            weights (str, list or ndarray):
                Mode for the weights used to coadd images. Options are 'auto' (default), 'uniform', or list/array of
                weights with shape = (nexp,) can be input and will be applied to the image. Note 'auto' is not allowed
                if offsets are input, and if set this will cause an exception.
            sn_smooth_npix (int): optional, default=None
                Number of pixels to median filter by when computing S/N used to decide how to scale and weight spectra. If
                set to None, the code will simply take 10% of the image size in the spectral direction.
                TODO: for truncated echelle orders we should be doing something more intelligent.
            ir_redux (bool): optional, default=False
                Is this an near-IR reduction, True=yes. This parameter is passed to pypeit.reduce for determining the
                reduction steps.
            show (bool): optional, default=False
                Show results to ginga
            show_peaks (bool): optional, default=False
                Show the QA for object finding algorithm peak finding to the screen.
            debug_offset (bool): optional, default=False
                Show QA for debugging the automatic determination of offsets to the screen.
            debug (bool): optional, default=False
                Show QA for debugging.
            **kwargs_wave
                Keyword arguments pass to `pypeit.core.coadd.get_wvae_grid` which determine how the wavelength grid
                is created for the 2d coadding.
        """

        ## Use Cases:
        #  1) offsets is None -- auto compute offsets from brightest object, so then default to auto_weights=True
        #  2) offsets not None, weights = None (uniform weighting) or weights is not None (input weights)
        #  3) offsets not None, auto_weights=True (Do not support)
        if offsets is not None and 'auto' in weights:
            msgs.error("Automatic weights cannot be computed for input offsets. "
                       "Set weights='uniform' or input an array of weights with shape (nexp,)")
        self.spec2d = spec2d
        self.spectrograph = spectrograph
        self.par = par
        self.det = det
        self.offsets = offsets
        self.weights = weights
        self.ir_redux = ir_redux
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
        self.stack_dict = self.load_coadd2d_stacks(self.spec2d)
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
        self.spat_ids = self.stack_dict['slits_list'][0].spat_id

        # If smoothing is not input, smooth by 10% of the spectral dimension
        self.sn_smooth_npix = sn_smooth_npix if sn_smooth_npix is not None else 0.1*self.nspec

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
        ..todo.. We need a proper doc string

        Args:
            only_slits:

        Returns:

        """

        only_slits = [only_slits] if (only_slits is not None and
                                      isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        # TODO We should be checking the bitmask for the reductions or something here??
        #  Yes, definitely
        if only_slits is None:
            slits = self.stack_dict['slits_list'][0]
            reduce_bpm = (slits.mask > 0) & (np.invert(slits.bitmask.flagged(
                slits.mask, flag=slits.bitmask.exclude_for_reducing)))
            good_slits = np.where(np.invert(reduce_bpm))[0]
        else:
            embed(header='DEAL WITH bitmask')

        coadd_list = []
        for slit_idx in good_slits:
            slitord_id = self.stack_dict['slits_list'][0].slitord_id[slit_idx]
            msgs.info('Performing 2d coadd for slit: {:d}/{:d}'.format(slit_idx, self.nslits - 1))
            ref_trace_stack = self.reference_trace_stack(slit_idx, offsets=self.offsets,
                                                         objid=self.objid_bri)
            thismask_stack = self.stack_dict['slitmask_stack'] == self.stack_dict['slits_list'][0].spat_id[slit_idx]
            # TODO Can we get rid of this one line simply making the weights returned by parse_weights an
            # (nslit, nexp) array?
            # This one line deals with the different weighting strategies between MultiSlit echelle. Otherwise, we
            # would need to copy this method twice in the subclasses
            if 'auto_echelle' in self.use_weights:
                rms_sn, weights = self.optimal_weights(slitord_id, self.objid_bri)
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
                                           self.wave_grid, weights=weights, interp_dspat=interp_dspat)
            coadd_list.append(coadd_dict)

        return coadd_list


    def create_pseudo_image(self, coadd_list):
        """
        ..todo.. see below

        THIS UNDOCUMENTED CODE PROBABLY SHOULD GENERATE AND RETURN
        STANDARD PYPEIT OBJCTS INSTEAD OF SOME UNDEFINED DICT"""


        # Masking
        # TODO -- Make this a method or something
        slits = self.stack_dict['slits_list'][0]
        reduce_bpm = (slits.mask > 0) & (np.invert(slits.bitmask.flagged(
            slits.mask, flag=slits.bitmask.exclude_for_reducing)))
        good_slits = np.where(np.invert(reduce_bpm))[0]

        nspec_vec = np.zeros(self.nslits,dtype=int)
        nspat_vec = np.zeros(self.nslits,dtype=int)
        for kk, cdict in enumerate(coadd_list):
            islit = good_slits[kk]
            nspec_vec[islit]=cdict['nspec']
            nspat_vec[islit]=cdict['nspat']

        # Determine the size of the pseudo image
        nspat_pad = 10
        nspec_pseudo = nspec_vec.max()
        nspat_pseudo = int(np.sum(nspat_vec) + (self.nslits + 1)*nspat_pad)  # Cast for SlitTraceSet
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
        for kk, coadd_dict in enumerate(coadd_list):
            islit = good_slits[kk]
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
                = slittrace.SlitTraceSet(slit_left, slit_righ, self.pypeline, nspat=nspat_pseudo,
                                         PYP_SPEC=self.spectrograph.spectrograph,
                                         specmin=spec_min1, specmax=spec_max1, ech_order=slits.ech_order)
                                         #master_key=self.stack_dict['master_key_dict']['trace'],
                                         #master_dir=self.master_dir)
        slitmask_pseudo = slits_pseudo.slit_img()
        # This is a kludge to deal with cases where bad wavelengths result in large regions where the slit is poorly sampled,
        # which wreaks havoc on the local sky-subtraction
        min_slit_frac = 0.70
        spec_min = np.zeros(self.nslits)
        spec_max = np.zeros(self.nslits)
        for slit_idx in good_slits:
            spat_id = slits_pseudo.spat_id[slit_idx]
            slit_width = np.sum(inmask_pseudo*(slitmask_pseudo == spat_id),axis=1)
            slit_width_img = np.outer(slit_width, np.ones(nspat_pseudo))
            med_slit_width = np.median(slit_width_img[slitmask_pseudo == spat_id])
            # TODO -- need inline docs
            nspec_eff = np.sum(slit_width > min_slit_frac*med_slit_width)
            nsmooth = int(np.fmax(np.ceil(nspec_eff*0.02),10))
            slit_width_sm = ndimage.filters.median_filter(slit_width, size=nsmooth, mode='reflect')
            igood = (slit_width_sm > min_slit_frac*med_slit_width)
            # TODO -- need inline docs
            spec_min[slit_idx] = spec_vec_pseudo[igood].min()
            spec_max[slit_idx] = spec_vec_pseudo[igood].max()
            bad_pix = (slit_width_img < min_slit_frac*med_slit_width) & (slitmask_pseudo == spat_id)
            inmask_pseudo[bad_pix] = False

        # Update slits_pseudo
        slits_pseudo.specmin = spec_min
        slits_pseudo.specmax = spec_max

        return dict(nspec=nspec_pseudo, nspat=nspat_pseudo, imgminsky=imgminsky_pseudo,
                    sciivar=sciivar_pseudo, inmask=inmask_pseudo, tilts=tilts_pseudo,
                    waveimg=waveimg_pseudo, spat_img=spat_img_pseudo, slits=slits_pseudo,
                    wave_mask=wave_mask, wave_mid=wave_mid, wave_min=wave_min, wave_max=wave_max)

    def reduce(self, pseudo_dict, show=None, show_peaks=None):
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
        #parcopy['calibrations']['save_masters'] = False
        #parcopy['scienceimage']['find_extrap_npoly'] = 1  # Use low order for trace extrapolation

        # Build the Calibrate object
        caliBrate = calibrations.Calibrations(None, self.par['calibrations'], self.spectrograph, None)
        caliBrate.slits = pseudo_dict['slits']


        redux=reduce.Reduce.get_instance(sciImage, self.spectrograph, parcopy, caliBrate,
                                         'science_coadd2d', ir_redux=self.ir_redux, det=self.det, show=show)
        #redux=reduce.Reduce.get_instance(sciImage, self.spectrograph, parcopy, pseudo_dict['slits'],
        #                                 None, None, 'science_coadd2d', ir_redux=self.ir_redux, det=self.det, show=show)
        # Set the tilts and waveimg attributes from the psuedo_dict here, since we generate these dynamically from fits
        # normally, but this is not possible for coadds
        redux.tilts = pseudo_dict['tilts']
        redux.waveimg = pseudo_dict['waveimg']
        redux.binning = self.binning

        # Masking
        #  TODO: Treat the masking of the slits objects
        #   from every exposure, come up with an aggregate mask (if it is masked on one slit,
        #   mask the slit for all) and that should be propagated into the slits object in the psuedo_dict
        slits = self.stack_dict['slits_list'][0]
        reduce_bpm = (slits.mask > 0) & (np.invert(slits.bitmask.flagged(
            slits.mask, flag=slits.bitmask.exclude_for_reducing)))
        redux.reduce_bpm = reduce_bpm

        if show:
            redux.show('image', image=pseudo_dict['imgminsky']*(sciImage.fullmask == 0), chname = 'imgminsky', slits=True, clear=True)

        # TODO:
        #  Object finding, this appears inevitable for the moment, since we need to be able to call find_objects
        #  outside of reduce. I think the solution here is to create a method in reduce for that performs the modified
        #  2d coadd reduce
        sobjs_obj, nobj, skymask_init = redux.find_objects(
            sciImage.image, show_peaks=show_peaks,
            manual_extract_dict=self.par['reduce']['extraction']['manual'].dict_for_objfind())

        # Local sky-subtraction
        global_sky_pseudo = np.zeros_like(pseudo_dict['imgminsky']) # No global sky for co-adds since we go straight to local
        skymodel_pseudo, objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs = redux.local_skysub_extract(
            global_sky_pseudo, sobjs_obj, spat_pix=pseudo_dict['spat_img'], model_noise=False,
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

        return pseudo_dict['imgminsky'], pseudo_dict['sciivar'], skymodel_pseudo, \
               objmodel_pseudo, ivarmodel_pseudo, outmask_pseudo, sobjs, sciImage.detector, pseudo_dict['slits'], \
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

    def get_good_slits(self, only_slits):
        """
        ..todo.. I need a doc string

        Args:
            only_slits:

        Returns:

        """

        only_slits = [only_slits] if (only_slits is not None and
                                        isinstance(only_slits, (int, np.int, np.int64, np.int32))) else only_slits
        good_slits = np.arange(self.nslits) if only_slits is None else only_slits
        return good_slits

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

        wave_grid, wave_grid_mid, dsamp = coadd.get_wave_grid(waves, masks=gpm, **kwargs_wave)
        return wave_grid, wave_grid_mid, dsamp

    def load_coadd2d_stacks(self, spec2d):
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
        #head2d_list = []
        specobjs_list = []
        slits_list = []
        nfiles =len(spec2d)
        detectors_list = []
        for ifile, f in enumerate(spec2d):
            if isinstance(f, spec2dobj.Spec2DObj):
                # If spec2d is a list of objects
                s2dobj = f
            else:
                # If spec2d is a list of files, option to also use spec1ds
                s2dobj = spec2dobj.Spec2DObj.from_file(f, self.det)
                spec1d_file = f.replace('spec2d', 'spec1d')
                if os.path.isfile(spec1d_file):
                    sobjs = specobjs.SpecObjs.from_fitsfile(spec1d_file)
                    this_det = sobjs.DET == self.det
                    specobjs_list.append(sobjs[this_det])
            # TODO the code should run without a spec1d file, but we need to implement that
            slits_list.append(s2dobj.slits)
            detectors_list.append(s2dobj.detector)
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
                    spectrograph=self.spectrograph.spectrograph,
                    pypeline=self.spectrograph.pypeline)

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
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(MultiSlitCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                        sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, par=par,
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
            self.objid_bri, self.spatid_bri, self.snr_bar_bri, self.offsets = self.compute_offsets()

        self.use_weights = self.parse_weights(weights)

    def parse_weights(self, weights):

        if 'auto' in weights:
            rms_sn, use_weights = self.optimal_weights(self.spatid_bri, self.objid_bri, const_weights=True)
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

        objid_bri, slitidx_bri, spatid_bri, snr_bar_bri = self.get_brightest_obj(self.stack_dict['specobjs_list'],
                                                                    self.spat_ids)
        msgs.info('Determining offsets using brightest object on slit: {:d} with avg SNR={:5.2f}'.format(spatid_bri,np.mean(snr_bar_bri)))
        thismask_stack = self.stack_dict['slitmask_stack'] == spatid_bri
        trace_stack_bri = np.zeros((self.nspec, self.nexp))
        # TODO Need to think abbout whether we have multiple tslits_dict for each exposure or a single one
        for iexp in range(self.nexp):
            trace_stack_bri[:,iexp] = self.stack_dict['slits_list'][iexp].center[:,slitidx_bri]
#            trace_stack_bri[:,iexp] = (self.stack_dict['tslits_dict_list'][iexp]['slit_left'][:,slitid_bri] +
#                                       self.stack_dict['tslits_dict_list'][iexp]['slit_righ'][:,slitid_bri])/2.0
        # Determine the wavelength grid that we will use for the current slit/order
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
        #specobj_dict = {'setup': 'unknown', 'slitid': 999, 'orderindx': 999, 'det': self.det, 'objtype': 'unknown',
        #                'pypeline': 'MultiSLit' + '_coadd_2d'}
        for iexp in range(self.nexp):
            sobjs_exp, _ = extract.objfind(sci_list_rebin[0][iexp,:,:], thismask, slit_left, slit_righ,
                                           inmask=inmask[iexp,:,:], ir_redux=self.ir_redux,
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
        # Print out a report on the offsets
        msg_string = msgs.newline()  + '---------------------------------------------'
        msg_string += msgs.newline() + ' Summary of offsets for highest S/N object   '
        msg_string += msgs.newline() + '         found on slitid = {:d}              '.format(spatid_bri)
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

        return objid_bri, spatid_bri, snr_bar_bri, offsets

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
                ithis = sobjs.SLITID == spat_id
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
    def __init__(self, spec2d_files, spectrograph, par, det=1, offsets=None, weights='auto', sn_smooth_npix=None,
                 ir_redux=False, show=False, show_peaks=False, debug_offsets=False, debug=False, **kwargs_wave):
        super(EchelleCoAdd2D, self).__init__(spec2d_files, spectrograph, det=det, offsets=offsets, weights=weights,
                                      sn_smooth_npix=sn_smooth_npix, ir_redux=ir_redux, par=par,
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
        else:
            # TODO -- Check the input offsets list matches the science images
            pass

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
        Utility function for determining the reference trace about
        which 2d coadds are performed.

        There are two modes of operation to determine the reference
        trace for the 2d coadd of a given slit/order:

            #. ``offsets``: We stack about the center of the slit for
               the slit in question with the input offsets added

            #. ``ojbid``: We stack about the trace ofa reference
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

