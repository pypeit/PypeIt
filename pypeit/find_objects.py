"""
Main driver class for object finding, global skysubtraction and skymask construction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import inspect
import numpy as np
import os

from astropy import stats
from abc import ABCMeta


from pypeit import specobjs
from pypeit import msgs, utils
from pypeit import flatfield
from pypeit.display import display
from pypeit.core import skysub, qa, parse, flat, flexure
from pypeit.core import procimg
from pypeit.images import buildimage
from pypeit.core import findobj_skymask

from IPython import embed


class FindObjects:
    """
    Base class used to find objects and perform global sky subtraction for
    science or standard-star exposures.

    Args:
        sciImg (:class:`~pypeit.images.scienceimage.ScienceImage`):
            Image to reduce.
        slits (:class:`~pypeit.slittrace.SlitTracSet`):
            Object providing slit traces for the image to reduce.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            PypeIt Sspectrograph class
        par (:class:`~pypeit.par.pyepeitpar.PypeItPar`):
            Reduction parameters class
        objtype (:obj:`str`):
            Specifies object being reduced.  Should be 'science',
            'standard', or 'science_coadd2d'.
        wv_calib (:class:`~pypeit.wavetilts.WaveCalib`, optional):
            This is only used for the IFU child when a joint sky subtraction
            is requested.
        waveTilts (:class:`~pypeit.wavetilts.WaveTilts`, optional):
            Calibration frame with arc/sky line tracing of the wavelength
            tilt.  Only waveTilts or tilts is needed (not both).
        tilts (`numpy.ndarray`_, optional):
            Tilts frame produced by
            :func:`~pypeit.wavetilts.WaveTilts.fit2tiltimg` for given a
            spatial flexure.  Only waveTilts or tilts is needed (not both).
        initial_skymask (`numpy.ndarray`_, optional):
            Boolean array that selects (array elements are True) image
            pixels in sky regions.  If provided, the 2nd pass on the global
            sky subtraction is omitted.
        bkg_redux (:obj:`bool`, optional):
            If True, the sciImg has been subtracted by
            a background image (e.g. standard treatment in the IR)
        find_negative (:obj:`bool`, optional):
            If True, the negative objects are found
        std_redux (:obj:`bool`, optional):
            If True, the object being extracted is a standard star,
            so that the reduction parameters can be adjusted accordingly.
        show (:obj:`bool`, optional):
            Show plots along the way?
        clear_ginga (:obj:`bool`, optional):
            Clear the ginga window before showing the object finding results.
        basename (:obj:`str`, optional):
            Base name for output files
        manual (:class:`~pypeit.manual_extract.ManualExtractObj`, optional):
            Object containing manual extraction instructions/parameters.

    Attributes:
        ivarmodel (`numpy.ndarray`_):
            Model of inverse variance
        objimage (`numpy.ndarray`_):
            Model of object
        skyimage (`numpy.ndarray`_):
            Final model of sky
        initial_sky (`numpy.ndarray`_):
            Initial sky model after first pass with global_skysub()
        global_sky (`numpy.ndarray`_):
            Fit to global sky
        skymask (`numpy.ndarray`_):
            Mask of the sky fit
        outmask (`numpy.ndarray`_):
            Final output mask
        extractmask (`numpy.ndarray`_):
            Extraction mask
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
        sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
            Objects found
        spat_flexure_shift (:obj:`float`):
        tilts (`numpy.ndarray`_):
            WaveTilts images generated on-the-spot
        waveimg (`numpy.ndarray`_):
            WaveImage image generated on-the-spot
        slitshift (`numpy.ndarray`_):
            Global spectral flexure correction for each slit (in pixels)
        vel_corr (:obj:`float`):
            Relativistic reference frame velocity correction (e.g. heliocentyric/barycentric/topocentric)

    """

    __metaclass__ = ABCMeta

    # TODO Consider removing objtype argument and simply have an optional parameter which regulates the flexure
    # behavior which is all objtype seems to do. But we should consider consistency with Extract.

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, sciImg, slits, spectrograph, par, objtype, wv_calib=None, waveTilts=None,
                     tilts=None, initial_skymask=None, bkg_redux=False, find_negative=False,
                     std_redux=False, show=False, clear_ginga=True, basename=None, manual=None):
        """
        Instantiate and return the :class:`FindObjects` subclass appropriate for
        the provided spectrograph.

        For argument descriptions, see :class:`FindObjects`.
        """
        return next(c for c in utils.all_subclasses(FindObjects)
                    if c.__name__ == (spectrograph.pypeline + 'FindObjects'))(
            sciImg, slits, spectrograph, par, objtype, wv_calib=wv_calib, waveTilts=waveTilts,
            tilts=tilts, initial_skymask=initial_skymask, bkg_redux=bkg_redux,
            find_negative=find_negative, std_redux=std_redux, show=show, clear_ginga=clear_ginga,
            basename=basename, manual=manual)

    def __init__(self, sciImg, slits, spectrograph, par, objtype, wv_calib=None, waveTilts=None,
                 tilts=None, initial_skymask=None, bkg_redux=False, find_negative=False,
                 std_redux=False, show=False, clear_ginga=True, basename=None, manual=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!
        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        self.scaleimg = np.array([1.0], dtype=float)  # np.array([1]) applies no scale
        self.basename = basename
        self.manual = manual
        self.initial_skymask = initial_skymask
        self.wv_calib = wv_calib  # TODO :: Ideally, we want to avoid this if possible. Find a better way to do joint_skysub fitting outside of the find_objects class.
        # Parse
        # Slit pieces
        #   WARNING -- It is best to unpack here then pass around self.slits
        #      Otherwise you have to keep in mind flexure, tweaking, etc.

        # TODO: The spatial flexure is not copied to the PypeItImage object if
        # the image (science or otherwise) is from a combination of multiple
        # frames.  Is that okay for this usage?
        # Flexure
        self.spat_flexure_shift = None
        if (objtype == 'science' and self.par['scienceframe']['process']['spat_flexure_correct']) or \
           (objtype == 'standard' and self.par['calibrations']['standardframe']['process']['spat_flexure_correct']):
            self.spat_flexure_shift = self.sciImg.spat_flexure
        elif objtype == 'science_coadd2d':
            self.spat_flexure_shift = None

        # Initialise the slits
        msgs.info("Initialising slits")
        self.initialise_slits(slits)

        # Internal bpm mask
        # We want to keep the 'BOXSLIT', which has bpm=2. But we don't want to keep 'BOXSLIT'
        # with other bad flag (for which bpm>2)
        self.reduce_bpm = (self.slits.mask > 2) & (np.invert(self.slits.bitmask.flagged(
                        self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing)))
        self.reduce_bpm_init = self.reduce_bpm.copy()

        # Load up other input items
        self.bkg_redux = bkg_redux
        self.find_negative = find_negative

        self.std_redux = std_redux
        # This can be a single integer for a single detector or a tuple for
        # multiple detectors placed in a mosaic.
        self.det = self.sciImg.detector.det
        # This is the string name of the detector or mosaic used when saving the
        # processed data to PypeIt's main output files
        self.detname = self.spectrograph.get_det_name(self.det)

        self.binning = self.sciImg.detector.binning
        self.pypeline = spectrograph.pypeline
        self.findobj_show = show

        self.steps = []

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.initial_sky = None
        self.skymask = None
        # TODO: Is this ever used?
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None
        self.slitshift = np.zeros(self.slits.nslits)  # Global spectral flexure slit shifts (in pixels) that are applied to all slits.
        self.vel_corr = None

        # Deal with dynamically generated calibrations, i.e. the tilts.
        if waveTilts is None and tilts is None:
            msgs.error("Must provide either waveTilts or tilts to FindObjects")
        elif waveTilts is not None and tilts is not None:
            msgs.error("Cannot provide both waveTilts and tilts to FindObjects")
        elif waveTilts is not None and tilts is None:
            self.waveTilts = waveTilts
            self.waveTilts.is_synced(self.slits)
            #   Deal with Flexure
            if self.par['calibrations']['tiltframe']['process']['spat_flexure_correct']:
                _spat_flexure = 0. if self.spat_flexure_shift is None else self.spat_flexure_shift
                # If they both shifted the same, there will be no reason to shift the tilts
                tilt_flexure_shift = _spat_flexure - self.waveTilts.spat_flexure
            else:
                tilt_flexure_shift = self.spat_flexure_shift
            msgs.info("Generating tilts image from fit in waveTilts")
            self.tilts = self.waveTilts.fit2tiltimg(self.slitmask, flexure=tilt_flexure_shift)
        elif waveTilts is None and tilts is not None:
            msgs.info("Using user input tilts image")
            self.tilts = tilts

        # Show?
        if self.findobj_show:
            self.show('image', image=sciImg.image, chname='processed', slits=True, clear=clear_ginga)


    def create_skymask(self, sobjs_obj):
        r"""
        Creates a skymask from a SpecObjs object

        Args:
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                Objects for which you would like to create the mask

        Returns:
            `numpy.ndarray`_: Boolean image with shape :math:`(N_{\rm spec},
            N_{\rm spat})` indicating which pixels are usable for global sky
            subtraction.  True = usable for sky subtraction, False = should be
            masked when sky subtracting.
        """
        # Masking options
        boxcar_rad_pix = None

        skymask = np.ones_like(self.sciImg.image, dtype=bool)
        gdslits = np.where(np.invert(self.reduce_bpm))[0]
        if sobjs_obj.nobj > 0:
            for slit_idx in gdslits:
                slit_spat = self.slits.spat_id[slit_idx]
                qa_title ="Generating skymask for slit # {:d}".format(slit_spat)
                msgs.info(qa_title)
                thismask = self.slitmask == slit_spat
                this_sobjs = sobjs_obj.SLITID == slit_spat
                # Boxcar mask?
                if self.par['reduce']['skysub']['mask_by_boxcar']:
                    boxcar_rad_pix = self.par['reduce']['extraction']['boxcar_radius'] / \
                                     self.get_platescale(slitord_id=self.slits.slitord_id[slit_idx])
                # Do it
                skymask[thismask] = findobj_skymask.create_skymask(sobjs_obj[this_sobjs], thismask,
                                                                   self.slits_left[:,slit_idx],
                                                                   self.slits_right[:,slit_idx],
                                                                   box_rad_pix=boxcar_rad_pix,
                                                                   trim_edg=self.par['reduce']['findobj']['find_trim_edge'])
        # Return
        return skymask

    def initialise_slits(self, slits, initial=False):
        """
        Gather all the :class:`SlitTraceSet` attributes
        that we'll use here in :class:`FindObjects`

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                SlitTraceSet object containing the slit boundaries that will be initialized.
            initial (:obj:`bool`, optional):
                Use the initial definition of the slits. If False,
                tweaked slits are used.
        """
        # Slits
        self.slits = slits
        # Select the edges to use
        self.slits_left, self.slits_right, _ \
            = self.slits.select_edges(initial=initial, flexure=self.spat_flexure_shift)

        # Slitmask
        self.slitmask = self.slits.slit_img(initial=initial, flexure=self.spat_flexure_shift,
                                            exclude_flag=self.slits.bitmask.exclude_for_reducing+['BOXSLIT'])
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)
#        # For echelle
#        self.spatial_coo = self.slits.spatial_coordinates(initial=initial, flexure=self.spat_flexure_shift)

    def run(self, std_trace=None, show_peaks=False, show_skysub_fit=False):
        """
        Primary code flow for object finding in PypeIt reductions

        Parameters
        ----------
        std_trace : `numpy.ndarray`_, optional
            Trace of the standard star
        show_peaks : :obj:`bool`, optional
            Show peaks in find_objects methods
        show_skysub_fit : :obj:`bool`, optional
            Show the fits for the global sky subtraction

        Returns
        -------
        initial_sky : `numpy.ndarray`_
            Initial global sky model
        sobjs_obj : :class:`~pypeit.specobjs.SpecObjs`
            List of objects found
        """

        # If the skip_skysub is set (i.e. image is already sky-subtracted), simply find objects
        if self.par['reduce']['findobj']['skip_skysub']:
            msgs.info("Skipping global sky sub as per user request")
            sobjs_obj, self.nobj = self.find_objects(self.sciImg.image, self.sciImg.ivar,
                                                     std_trace=std_trace, show=self.findobj_show,
                                                     show_peaks=show_peaks)
            return np.zeros_like(self.sciImg.image), sobjs_obj

        # Perform a first pass sky-subtraction.  The mask is either empty or
        # uses the mask specified by the user.
        
        # TODO: Should we make this no_poly=True to have fewer degrees of freedom in
        # the with with-object global sky fits??
        initial_sky0 = self.global_skysub(skymask=self.initial_skymask, update_crmask=False,
                                          objs_not_masked=True, show_fit=show_skysub_fit)
        # First pass object finding
        save_objfindQA = self.par['reduce']['findobj']['skip_second_find'] or self.std_redux \
                            or self.initial_skymask is not None 
        sobjs_obj, self.nobj = \
            self.find_objects(self.sciImg.image-initial_sky0, self.sciImg.ivar, std_trace=std_trace,
                              show_peaks=show_peaks, show=self.findobj_show and not self.std_redux,
                              save_objfindQA=save_objfindQA)

        if self.nobj == 0 or self.initial_skymask is not None:
            # Either no objects were found, or the initial sky mask was provided by the user.
            # Either way, don't don't redo global sky subtraction
            msgs.info('Either no objects were found or a user-provided sky mask was used.  '
                      'Skipping second pass of sky-subtraction and object finding.')
            return initial_sky0, sobjs_obj

        # If objects were found, create skymask using first pass objects that
        # were identified, sobjs_obj
        skymask_init = self.create_skymask(sobjs_obj)
        # Global sky subtract now using the skymask defined by object positions
        initial_sky = self.global_skysub(skymask=skymask_init, show_fit=show_skysub_fit)
        # Second pass object finding on sky-subtracted image with updated sky
        # created after masking objects
        if not self.std_redux and not self.par['reduce']['findobj']['skip_second_find']:
            sobjs_obj, self.nobj = self.find_objects(self.sciImg.image - initial_sky,
                                                     self.sciImg.ivar, std_trace=std_trace,
                                                     show=self.findobj_show, show_peaks=show_peaks)
        else:
            msgs.info("Skipping 2nd run of finding objects")
        # TODO I think the final global should go here as well from the pypeit.py class lines 837
        return initial_sky, sobjs_obj

    def find_objects(self, image, ivar, std_trace=None,
                     show_peaks=False, show_fits=False,
                     show_trace=False, show=False, save_objfindQA=True,
                     manual_extract_dict=None, debug=False):
        """
        Single pass at finding objects in the input image

        If self.find_negative is True, do a search for negative objects too

        Parameters
        ----------
        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC
            modules window
        show : :obj:`bool`, optional
            Show all the QA
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        manual_extract_dict : :obj:`dict`, optional
            This is only used by 2D coadd
        debug : :obj:`bool`, optional
            Show debugging plots?

        Returns
        -------
        sobjs_obj_single : :class:`~pypeit.specobjs.SpecObjs`
            Objects found
        nobj_single : :obj:`int`
            Number of objects found
        """
        # Positive image
        if manual_extract_dict is None:
            manual_extract_dict= self.manual.dict_for_objfind(self.detname, neg=False) if self.manual is not None else None

        sobjs_obj_single, nobj_single = \
            self.find_objects_pypeline(image, ivar,
                                       std_trace=std_trace,
                                       show_peaks=show_peaks, show_fits=show_fits,
                                       show_trace=show_trace, save_objfindQA=save_objfindQA,
                                       manual_extract_dict=manual_extract_dict, 
                                       neg=False, debug=debug)

        # Find negative objects
        if self.find_negative:
            msgs.info("Finding objects in the negative image")
            # Parses
            manual_extract_dict = self.manual.dict_for_objfind(self.detname, neg=True) if self.manual is not None else None
            sobjs_obj_single_neg, nobj_single_neg = \
                self.find_objects_pypeline(-image, ivar, std_trace=std_trace,
                                           show_peaks=show_peaks, show_fits=show_fits,
                                           show_trace=show_trace, save_objfindQA=save_objfindQA,
                                           manual_extract_dict=manual_extract_dict, neg=True,
                                           debug=debug)
            # Add (if there are any)
            sobjs_obj_single.append_neg(sobjs_obj_single_neg)

        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='objfind',
                      sobjs=sobjs_obj_single, slits=True)

        # For nobj we take only the positive objects
        return sobjs_obj_single, nobj_single

    # TODO maybe we don't need parent and children for this method. But IFU has a bunch of extra methods.
    def find_objects_pypeline(self, image, ivar, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):

        """
         Dummy method for object finding. Overloaded by class specific object finding.

         Returns:

         """
        return None, None

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Over-loaded by the children

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        pass


    def global_skysub(self, skymask=None, update_crmask=True, trim_edg = (0, 0),
                      previous_sky=None, show_fit=False, show=False, show_objs=False, objs_not_masked=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Args:
            skymask (`numpy.ndarray`_, None):
                A 2D image indicating sky regions (1=sky)
            update_crmask (bool, optional):
            trim_edg (tuple, optional):
                 A two tuple of ints that specify the number of pixels to trim from the slit edges
            show_fit (bool, optional):
            show (bool, optional):
            show_objs (bool, optional):
            previous_sky (`numpy.ndarray`_, optional):
                Sky model estimate from a previous run of global_sky
                Used to generate an improved estimated of the variance
            objs_not_masked (bool, optional):
                Set this to be True if there are objects on the slit/order that are not being masked
                by the skymask. This is typically the case for the first pass sky-subtraction
                before object finding, since a skymask has not yet been created.

        Returns:
            `numpy.ndarray`_: image of the the global sky model

        """
        # reset bpm since global sky is run several times and reduce_bpm is here updated.
        self.reduce_bpm = self.reduce_bpm_init.copy()
        # Prep
        global_sky = np.zeros_like(self.sciImg.image)
        # Parameters for a standard star
        if self.std_redux:
            sigrej = 7.0
            update_crmask = False
            if not self.par['reduce']['skysub']['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return global_sky
        else:
            sigrej = 3.0

        # We use this tmp bpm so that we exclude the BOXSLITS during the global_skysub
        tmp_bpm = (self.slits.mask > 0) & \
                          (np.invert(self.slits.bitmask.flagged(self.slits.mask,
                                                                flag=self.slits.bitmask.exclude_for_reducing)))
        gdslits = np.where(np.invert(tmp_bpm))[0]

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)

        # Allow for previous sky to better estimate ivar
        #  Unless we used a background image (i.e. bkg_redux=True)
        if (previous_sky is not None) and (not self.bkg_redux):
            # Estimate the variance using the input sky model
            var = procimg.variance_model(self.sciImg.base_var,
                                          counts=previous_sky,
                                          count_scale=self.sciImg.img_scale,
                                          noise_floor=self.sciImg.noise_floor)
            skysub_ivar = utils.inverse(var)
        else:
            skysub_ivar = self.sciImg.ivar


        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            msgs.info("Global sky subtraction for slit: {:d}".format(slit_spat))
            thismask = self.slitmask == slit_spat
            inmask = self.sciImg.select_flag(invert=True) & thismask & skymask_now
            # All masked?
            if not np.any(inmask):
                msgs.warn("No pixels for fitting sky.  If you are using mask_by_boxcar=True, your radius may be too large.")
                self.reduce_bpm[slit_idx] = True
                continue

            # Find sky
            global_sky[thismask] = skysub.global_skysub(
                self.sciImg.image, skysub_ivar, self.tilts, thismask,
                self.slits_left[:,slit_idx], self.slits_right[:,slit_idx],
                inmask=inmask, sigrej=sigrej,
                bsp=self.par['reduce']['skysub']['bspline_spacing'],
                trim_edg=tuple(self.par['reduce']['trim_edge']),
                no_poly=self.par['reduce']['skysub']['no_poly'],
                pos_mask=not self.bkg_redux and not objs_not_masked,
                max_mask_frac=self.par['reduce']['skysub']['max_mask_frac'],
                show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(global_sky[thismask]) == 0.:
                msgs.warn("Bad fit to sky.  Rejecting slit: {:d}".format(slit_spat))
                self.reduce_bpm[slit_idx] = True

        if update_crmask and self.par['scienceframe']['process']['mask_cr']:
            # Find CRs with sky subtraction
            # TODO: Shouldn't the saturation flagging account for the
            # subtraction of the sky?
            self.sciImg.build_crmask(self.par['scienceframe']['process'],
                                     subtract_img=global_sky)
            # TODO: This mask update is done *inside* build_crmask.
#            # Update the fullmask
#            self.sciImg.update_mask_cr(self.sciImg.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', global_sky=global_sky, slits=True, sobjs=sobjs_show, clear=False)

        # Return
        return global_sky

    def show(self, attr, image=None, global_sky=None, showmask=False, sobjs=None,
             chname=None, slits=False,clear=False):
        """
        Show one of the internal images

        .. todo::
            Should probably put some of these in ProcessImages

        Parameters
        ----------
        attr : str
          global -- Sky model (global)
          sci -- Processed science image
          rawvar -- Raw variance image
          modelvar -- Model variance image
          crmasked -- Science image with CRs set to 0
          skysub -- Science image with global sky subtracted
          image -- Input image
        display : str, optional
        image : ndarray, optional
          User supplied image to display

        Returns
        -------

        """
        mask_in = self.sciImg.fullmask if showmask else None

        img_gpm = self.sciImg.select_flag(invert=True)

        if attr == 'global' and all([a is not None for a in [self.sciImg.image, global_sky,
                                                             self.sciImg.fullmask]]):
            # global sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - global_sky) * img_gpm.astype(float)
            mean, med, sigma = stats.sigma_clipped_stats(image[img_gpm], sigma_lower=5.0,
                                                         sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            ch_name = chname if chname is not None else f'global_sky_{self.detname}'
            viewer, ch = display.show_image(image, chname=ch_name, mask=mask_in, clear=clear,
                                            wcs_match=True)
        elif attr == 'image':
            ch_name = chname if chname is not None else 'image'
            viewer, ch = display.show_image(image, chname=ch_name, clear=clear, wcs_match=True)
        else:
            msgs.warn("Not an option for show")

        if sobjs is not None:
            for spec in sobjs:
                color = 'magenta' if spec.hand_extract_flag else 'orange'
                display.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        if slits and self.slits_left is not None:
            display.show_slits(viewer, ch, self.slits_left, self.slits_right)

    def __repr__(self):
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                        self.nsci)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt


class MultiSlitFindObjects(FindObjects):
    """
    Child of Reduce for Multislit and Longslit reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, spectrograph, par, objtype, **kwargs)

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        bin_spec, bin_spat = parse.parse_binning(self.binning)
        return self.sciImg.detector.platescale * bin_spec

    def find_objects_pypeline(self, image, ivar, std_trace=None,
                              manual_extract_dict=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False):
        """
        Pipeline specific find objects routine

        Parameters
        ----------

        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        manual_extract_dict : :obj:`dict`, optional
            Dict guiding the manual extraction
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC
            modules window
        show : :obj:`bool`, optional
            Show all the QA
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        neg : :obj:`bool`, optional
            Is this a negative image?
        debug : :obj:`bool`, optional
            Show debugging plots?

        Returns
        -------
        specobjs : :class:`~pypeot.specobjs.Specobjs`
            Container holding Specobj objects
        nobj : :obj:`int`
            Number of objects identified
        """
        gdslits = np.where(np.invert(self.reduce_bpm))[0]

        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            qa_title ="Finding objects on slit # {:d}".format(slit_spat)
            msgs.info(qa_title)
            thismask = self.slitmask == slit_spat
            inmask = self.sciImg.select_flag(invert=True) & thismask
            specobj_dict = {'SLITID': slit_spat,
                            'DET': self.sciImg.detector.name,
                            'OBJTYPE': self.objtype,
                            'PYPELINE': self.pypeline}

            # This condition allows to not use a threshold to find objects in alignment boxes
            # because these boxes are smaller than normal slits and the stars are very bright,
            # the detection threshold would be too high and the star not detected.
            if self.slits.bitmask.flagged(self.slits.mask[slit_idx], flag='BOXSLIT'):
                snr_thresh = 0.
            else:
                snr_thresh = self.par['reduce']['findobj']['snr_thresh']

            # Set objfind QA filename
            objfindQA_filename = None
            if save_objfindQA and (self.basename is not None):
                out_dir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])
                if self.find_negative:
                    basename = 'neg_' + self.basename if neg else 'pos_' + self.basename
                else:
                    basename = self.basename
                objfindQA_filename = qa.set_qa_filename(basename, 'obj_profile_qa', slit=slit_spat,
                                                        det=self.detname, out_dir=out_dir)

            maxnumber =  self.par['reduce']['findobj']['maxnumber_std'] if self.std_redux \
                else self.par['reduce']['findobj']['maxnumber_sci']
            sobjs_slit = \
                    findobj_skymask.objs_in_slit(image, ivar, thismask,
                                self.slits_left[:,slit_idx],
                                self.slits_right[:,slit_idx],
                                inmask=inmask,
                                ncoeff=self.par['reduce']['findobj']['trace_npoly'],
                                std_trace=std_trace,
                                snr_thresh=snr_thresh,
                                hand_extract_dict=manual_extract_dict,
                                specobj_dict=specobj_dict, show_peaks=show_peaks,
                                show_fits=show_fits, show_trace=show_trace,
                                trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                                fwhm=self.par['reduce']['findobj']['find_fwhm'],
                                use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
                                boxcar_rad=self.par['reduce']['extraction']['boxcar_radius'] / self.get_platescale(),  #pixels
                                maxdev=self.par['reduce']['findobj']['find_maxdev'],
                                find_min_max=self.par['reduce']['findobj']['find_min_max'],
                                extract_maskwidth=self.par['reduce']['skysub']['local_maskwidth'],
                                qa_title=qa_title, nperslit=maxnumber,
                                objfindQA_filename=objfindQA_filename,
                                debug_all=debug)
            # Record
            sobjs.add_sobj(sobjs_slit)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='objfind', sobjs=sobjs,
                      slits=True)

        # Return
        return sobjs, len(sobjs)


class EchelleFindObjects(FindObjects):
    """
    Child of Reduce for Echelle reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, spectrograph, par, objtype, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
                            else self.slits.ech_order
#                            else self.spectrograph.order_vec(self.spatial_coo)
        if self.order_vec is None:
            msgs.error('Unable to set Echelle orders, likely because they were incorrectly '
                       'assigned in the relevant SlitTraceSet.')

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, IFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        if slitord_id is None:
            msgs.error('slitord_id is missing. Plate scale for current echelle order cannot be determined.')
        return self.spectrograph.order_platescale(slitord_id, binning=self.binning)[0]


    def find_objects_pypeline(self, image, ivar, std_trace=None,
                              show=False, show_peaks=False, show_fits=False,
                              show_trace=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):
        """
        Pipeline specific find objects routine

        Parameters
        ----------
        image : `numpy.ndarray`_
            Image to search for objects from. This floating-point image has shape
            (nspec, nspat) where the first dimension (nspec) is
            spectral, and second dimension (nspat) is spatial.
        std_trace : `numpy.ndarray`_, optional
            This is a one dimensional float array with shape = (nspec,) containing the standard star
            trace which is used as a crutch for tracing. If the no
            standard star is provided the code uses the the slit
            boundaries as the crutch.
        manual_extract_dict : :obj:`dict`, optional
            Dict guiding the manual extraction
        show_peaks : :obj:`bool`, optional
            Generate QA showing peaks identified by object finding
        show_fits : :obj:`bool`, optional
            Generate QA  showing fits to traces
        show_trace : :obj:`bool`, optional
            Generate QA  showing traces identified. Requires an open ginga RC modules window
        save_objfindQA : :obj:`bool`, optional
            Save to disk (png file) QA showing the object profile
        neg : :obj:`bool`, optional
            Is this a negative image?
        show : :obj:`bool`, optional
        debug : :obj:`bool`, optional

        Returns
        -------
        specobjs : :class:`~pypeit.specobjs.Specobjs`
            Container holding Specobj objects
        nobj : :obj:`int`
            Number of objects identified
        """

        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        inmask = self.sciImg.select_flag(invert=True)
        # Find objects
        # TODO: Not sure how this fairs if self.det is a tuple...
        specobj_dict = {'SLITID': 999, 'DET': self.sciImg.detector.name, 
                        'ECH_ORDERINDX': 999,
                        'OBJTYPE': self.objtype,
                        'PYPELINE': self.pypeline}

        # Set objfind QA filename
        objfindQA_filename = None
        if save_objfindQA and (self.basename is not None):
            out_dir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])
            if self.find_negative:
                basename = 'neg_' + self.basename if neg else 'pos_' + self.basename
            else:
                basename = self.basename
            objfindQA_filename = qa.set_qa_filename(basename, 'obj_profile_qa', slit=999,
                                                    det=self.detname, out_dir=out_dir)

        #This could cause problems if there are more than one object on the echelle slit, i,e, this tacitly
        #assumes that the standards for echelle have a single object. If this causes problems, we could make an
        #nperorder_std as a parameter in the parset that the user can adjust.
        nperorder =  self.par['reduce']['findobj']['maxnumber_std'] if self.std_redux \
            else self.par['reduce']['findobj']['maxnumber_sci']

        sobjs_ech = findobj_skymask.ech_objfind(
            image, ivar, self.slitmask, self.slits_left, self.slits_right,
            self.order_vec, self.reduce_bpm, 
            self.slits.spat_id,
            np.vstack((self.slits.specmin, self.slits.specmax)),
            det=self.det,
            inmask=inmask, 
            ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            manual_extract_dict=manual_extract_dict, 
            plate_scale=plate_scale,
            std_trace=std_trace,
            specobj_dict=specobj_dict,
            snr_thresh=self.par['reduce']['findobj']['snr_thresh'],
            show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
            fwhm=self.par['reduce']['findobj']['find_fwhm'],
            use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
            maxdev=self.par['reduce']['findobj']['find_maxdev'],
            nperorder=nperorder,
            max_snr=self.par['reduce']['findobj']['ech_find_max_snr'],
            min_snr=self.par['reduce']['findobj']['ech_find_min_snr'],
            nabove_min_snr=self.par['reduce']['findobj']['ech_find_nabove_min_snr'],
            box_radius=self.par['reduce']['extraction']['boxcar_radius'],  # arcsec
            show_trace=show_trace, objfindQA_filename=objfindQA_filename)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float), chname='ech_objfind',
                      sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech)


class IFUFindObjects(MultiSlitFindObjects):
    """
    Child of Reduce for IFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, spectrograph, par, objtype, **kwargs)
        self.initialise_slits(slits, initial=True)

    def find_objects_pypeline(self, image, ivar, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False,
                              manual_extract_dict=None):
        """
        See MultiSlitReduce for slit-based IFU reductions
        """
        if self.par['reduce']['cube']['slit_spec']:
            return super().find_objects_pypeline(image, ivar, std_trace=std_trace,
                                                 show_peaks=show_peaks, show_fits=show_fits, show_trace=show_trace,
                                                 show=show, save_objfindQA=save_objfindQA, neg=neg,
                                                 debug=debug, manual_extract_dict=manual_extract_dict)
        return None, None, None

    def apply_relative_scale(self, scaleImg):
        """Apply a relative scale to the science frame (and correct the varframe, too)

         Args:
             scaleImg (`numpy.ndarray`_):
                scale image to divide the science frame by
        """
        # Check that scaleimg is set to the correct shape
        if self.scaleimg.size == 1:
            self.scaleimg = np.ones_like(self.sciImg.image)
        # Correct the relative illumination of the science frame
        msgs.info("Correcting science frame for relative illumination")
        self.scaleimg *= scaleImg.copy()
        self.sciImg.image, _bpm, varImg = flat.flatfield(self.sciImg.image, scaleImg,
                                                         varframe=utils.inverse(self.sciImg.ivar))
        if np.any(_bpm):
            self.sciImg.update_mask('BADSCALE', indx=_bpm)
        self.sciImg.ivar = utils.inverse(varImg)

    # TODO :: THIS FUNCTION IS NOT CURRENTLY USED, BUT RJC REQUESTS TO KEEP THIS CODE HERE FOR THE TIME BEING.
    # def illum_profile_spatial(self, skymask=None, trim_edg=(0, 0), debug=False):
    #     """
    #     Calculate the residual spatial illumination profile using the sky regions.
    #
    #     The redisual is calculated using the differential:
    #
    #     .. code-block:: python
    #
    #         correction = amplitude * (1 + spatial_shift * (dy/dx)/y)
    #
    #     where ``y`` is the spatial profile determined from illumflat, and
    #     spatial_shift is the residual spatial flexure shift in units of pixels.
    #
    #      Args:
    #         skymask (`numpy.ndarray`_):
    #             Mask of sky regions where the spatial illumination will be determined
    #         trim_edg (:obj:`tuple`):
    #             A tuple of two ints indicated how much of the slit edges should be
    #             trimmed when fitting to the spatial profile.
    #         debug (:obj:`bool`):
    #             Show debugging plots?
    #     """
    #
    #     msgs.info("Performing spatial sensitivity correction")
    #     # Setup some helpful parameters
    #     skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
    #     hist_trim = 0  # Trim the edges of the histogram to take into account edge effects
    #     gpm = self.sciImg.select_flag(invert=True)
    #     slitid_img_init = self.slits.slit_img(pad=0, initial=True, flexure=self.spat_flexure_shift)
    #     spatScaleImg = np.ones_like(self.sciImg.image)
    #     # For each slit, grab the spatial coordinates and a spline
    #     # representation of the spatial profile from the illumflat
    #     rawimg = self.sciImg.image.copy()
    #     numbins = int(np.max(self.slits.get_slitlengths(initial=True, median=True)))
    #     spatbins = np.linspace(0.0, 1.0, numbins + 1)
    #     spat_slit = 0.5 * (spatbins[1:] + spatbins[:-1])
    #     slitlength = np.median(self.slits.get_slitlengths(median=True))
    #     coeff_fit = np.zeros((self.slits.nslits, 2))
    #     for sl, slitnum in enumerate(self.slits.spat_id):
    #         msgs.info("Deriving spatial correction for slit {0:d}/{1:d}".format(sl + 1, self.slits.spat_id.size))
    #         # Get the initial slit locations
    #         onslit_b_init = (slitid_img_init == slitnum)
    #
    #         # Synthesize ximg, and edgmask from slit boundaries. Doing this outside this
    #         # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.
    #         spatcoord, edgmask = pixels.ximg_and_edgemask(self.slits_left[:, sl], self.slits_right[:, sl],
    #                                                       onslit_b_init, trim_edg=trim_edg)
    #
    #         # Make the model histogram
    #         xspl = np.linspace(0.0, 1.0, 10 * int(slitlength))  # Sub sample each pixel with 10 subpixels
    #         # TODO: caliBrate is no longer a dependency. If you need these b-splines pass them in.
    #         modspl = self.caliBrate.flatimages.illumflat_spat_bsplines[sl].value(xspl)[0]
    #         gradspl = interpolate.interp1d(xspl, np.gradient(modspl) / modspl, kind='linear', bounds_error=False,
    #                                        fill_value='extrapolate')
    #
    #         # Ignore skymask
    #         coord_msk = onslit_b_init & gpm
    #         hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
    #         cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
    #         hist_slit_all = hist / (cntr + (cntr == 0))
    #         histmod, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=gradspl(spatcoord[coord_msk]))
    #         hist_model = histmod / (cntr + (cntr == 0))
    #
    #         # Repeat with skymask
    #         coord_msk = onslit_b_init & gpm & skymask_now
    #         hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
    #         cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
    #         hist_slit = hist / (cntr + (cntr == 0))
    #
    #         # Prepare for fit - take the non-zero elements and trim slit edges
    #         if hist_trim == 0:
    #             ww = (hist_slit != 0)
    #             xfit = spat_slit[ww]
    #             yfit = hist_slit_all[ww]
    #             mfit = hist_model[ww]
    #         else:
    #             ww = (hist_slit[hist_trim:-hist_trim] != 0)
    #             xfit = spat_slit[hist_trim:-hist_trim][ww]
    #             yfit = hist_slit_all[hist_trim:-hist_trim][ww]
    #             mfit = hist_model[hist_trim:-hist_trim][ww]
    #
    #         # Fit the function
    #         spat_func = lambda par, ydata, model: par[0]*(1 + par[1] * model) - ydata
    #         res_lsq = least_squares(spat_func, [np.median(yfit), 0.0], args=(yfit, mfit))
    #         spatnorm = spat_func(res_lsq.x, 0.0, gradspl(spatcoord[onslit_b_init]))
    #         spatnorm /= spat_func(res_lsq.x, 0.0, gradspl(0.5))
    #         # Set the scaling factor
    #         spatScaleImg[onslit_b_init] = spatnorm
    #         coeff_fit[sl, :] = res_lsq.x
    #
    #     if debug:
    #         from matplotlib import pyplot as plt
    #         xplt = np.arange(24)
    #         plt.subplot(121)
    #         plt.plot(xplt[0::2], coeff_fit[::2, 0], 'rx')
    #         plt.plot(xplt[1::2], coeff_fit[1::2, 0], 'bx')
    #         plt.subplot(122)
    #         plt.plot(xplt[0::2], coeff_fit[::2, 1]/10, 'rx')
    #         plt.plot(xplt[1::2], coeff_fit[1::2, 1]/10, 'bx')
    #         plt.show()
    #         plt.imshow(spatScaleImg, vmin=0.99, vmax=1.01)
    #         plt.show()
    #         plt.subplot(133)
    #         plt.plot(xplt[0::2], coeff_fit[::2, 2], 'rx')
    #         plt.plot(xplt[1::2], coeff_fit[1::2, 2], 'bx')
    #         plt.show()
    #     # Apply the relative scale correction
    #     self.apply_relative_scale(spatScaleImg)

    def illum_profile_spectral(self, global_sky, skymask=None):
        """Calculate the residual spectral illumination profile using the sky regions.
        This uses the same routine as the flatfield spectral illumination profile.

         Args:
             global_sky (`numpy.ndarray`_):
                Model of the sky
             skymask (`numpy.ndarray`_, optional):
                Mask of sky regions where the spectral illumination will be determined
        """
        trim = self.par['calibrations']['flatfield']['slit_trim']
        ref_idx = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
        smooth_npix = self.par['calibrations']['flatfield']['slit_illum_smooth_npix']
        gpm = self.sciImg.select_flag(invert=True)
        # Note :: Need to provide wavelength to illum_profile_spectral (not the tilts) so that the
        # relative spectral sensitivity is calculated at a given wavelength for all slits simultaneously.
        scaleImg = flatfield.illum_profile_spectral(self.sciImg.image.copy(), self.waveimg, self.slits,
                                                    slit_illum_ref_idx=ref_idx, model=global_sky, gpmask=gpm,
                                                    skymask=skymask, trim=trim, flexure=self.spat_flexure_shift,
                                                    smooth_npix=smooth_npix)
        # Now apply the correction to the science frame
        self.apply_relative_scale(scaleImg)

    def joint_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                     show_fit=False, show=False, show_objs=False, adderr=0.01, objs_not_masked=False):
        """ Perform a joint sky model fit to the data. See Reduce.global_skysub()
        for parameter definitions.
        """
        msgs.info("Performing joint global sky subtraction")
        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        global_sky = np.zeros_like(self.sciImg.image)
        thismask = (self.slitmask > 0)
        inmask = (self.sciImg.select_flag(invert=True) & thismask & skymask_now).astype(bool)
        # Convert the wavelength image to A/pixel, registered at pixel 0 (this gives something like
        # the tilts frame, but conserves wavelength position in each slit)
        wavemin = self.waveimg[self.waveimg != 0.0].min()
        tilt_wave = (self.waveimg - wavemin) / (self.waveimg.max() - wavemin)

        # Parameters for a standard star
        sigrej = 3.0
        if self.std_redux:
            sigrej = 7.0
            update_crmask = False
            if not self.par['reduce']['skysub']['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return global_sky

        # Iterate to use a model variance image
        numiter = 4
        model_ivar = self.sciImg.ivar
        for nn in range(numiter):
            msgs.info("Performing iterative joint sky subtraction - ITERATION {0:d}/{1:d}".format(nn+1, numiter))
            # TODO trim_edg is in the parset so it should be passed in here via trim_edg=tuple(self.par['reduce']['trim_edge']),
            global_sky[thismask] = skysub.global_skysub(self.sciImg.image, model_ivar, tilt_wave,
                                                        thismask, self.slits_left, self.slits_right, inmask=inmask,
                                                        sigrej=sigrej, trim_edg=trim_edg,
                                                        bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                                        no_poly=self.par['reduce']['skysub']['no_poly'],
                                                        pos_mask=not self.bkg_redux and not objs_not_masked,
                                                        max_mask_frac=self.par['reduce']['skysub']['max_mask_frac'],
                                                        show_fit=show_fit)
            # Update the ivar image used in the sky fit
            msgs.info("Updating sky noise model")
            # Choose the highest counts out of sky and object
            counts = global_sky
            _scale = None if self.sciImg.img_scale is None else self.sciImg.img_scale[thismask]
            # NOTE: darkcurr must be a float for the call below to work.
            var = procimg.variance_model(self.sciImg.base_var[thismask], counts=counts[thismask],
                                         count_scale=_scale, noise_floor=adderr)
            model_ivar[thismask] = utils.inverse(var)
            # Redo the relative spectral illumination correction with the improved sky model
            self.illum_profile_spectral(global_sky, skymask=thismask)

        if update_crmask:
            # Find CRs with sky subtraction
            # NOTE: There's no need to run `sciImg.update_mask_cr` after this.
            # This operation updates the mask directly!
            self.sciImg.build_crmask(self.par['scienceframe']['process'],
                                     subtract_img=global_sky)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', global_sky=global_sky, slits=True, sobjs=sobjs_show, clear=False)
        return global_sky

    def global_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                      previous_sky=None, show_fit=False, show=False, show_objs=False, objs_not_masked=False):
        """
        Perform global sky subtraction. This IFU-specific routine ensures that the
        edges of the slits are not trimmed, and performs a spatial and spectral
        correction using the sky spectrum, if requested. See Reduce.global_skysub()
        for parameter definitions.
        """
        # Generate a global sky sub for all slits separately
        global_sky_sep = super().global_skysub(skymask=skymask, update_crmask=update_crmask,
                                               trim_edg=trim_edg, show_fit=show_fit, show=show,
                                               show_objs=show_objs)

        # Check if flexure or a joint fit is requested
        if not self.par['reduce']['skysub']['joint_fit'] and self.par['flexure']['spec_method'] == 'skip':
            return global_sky_sep
        if self.wv_calib is None:
            msgs.error("A wavelength calibration is needed (wv_calib) if a joint sky fit is requested.")
        msgs.info("Generating wavelength image")
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)
        # Calculate spectral flexure
        method = self.par['flexure']['spec_method']
        if method in ['slitcen']:
            trace_spat = 0.5 * (self.slits_left + self.slits_right)
            gd_slits = np.ones(self.slits.nslits, dtype=bool)
            flex_list = flexure.spec_flexure_slit_global(self.sciImg, self.waveimg, global_sky_sep, self.par,
                                                         self.slits, self.slitmask, trace_spat, gd_slits,
                                                         self.wv_calib, self.pypeline, self.det)
            for sl in range(self.slits.nslits):
                self.slitshift[sl] = flex_list[sl]['shift'][0]
                msgs.info("Flexure correction of slit {0:d}: {1:.3f} pixels".format(1 + sl, self.slitshift[sl]))

        # If the joint fit or spec/spat sensitivity corrections are not being performed, return the separate slits sky
        if not self.par['reduce']['skysub']['joint_fit']:
            return global_sky_sep

        # Do the spatial scaling first
        # if self.par['scienceframe']['process']['use_illumflat']:
        #     # Perform the correction
        #     self.illum_profile_spatial(skymask=skymask)
        #     # Re-generate a global sky sub for all slits separately
        #     global_sky_sep = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
        #                                           show_fit=show_fit, show=show, show_objs=show_objs)

        # Recalculate the wavelength image, and the global sky taking into account the spectral flexure
        msgs.info("Generating wavelength image, accounting for spectral flexure")
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spec_flexure=self.slitshift,
                                                   spat_flexure=self.spat_flexure_shift)

        self.illum_profile_spectral(global_sky_sep, skymask=skymask)

        # Use sky information in all slits to perform a joint sky fit
        global_sky = self.joint_skysub(skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
                                       show_fit=show_fit, show=show, show_objs=show_objs,
                                       objs_not_masked=objs_not_masked)

        return global_sky

