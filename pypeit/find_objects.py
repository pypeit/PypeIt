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
from pypeit.core import findobj_skymask

from IPython import embed


class FindObjects:
    """
    Base class used to find objects and perform global sky subtraction for
    science or standard-star exposures.

    Args:
        sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
            Image to reduce.
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
            Object providing slit traces for the image to reduce.
        spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            PypeIt Spectrograph class
        par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
            Reduction parameters class
        objtype (:obj:`str`):
            Specifies object being reduced.  Should be 'science',
            'standard', or 'science_coadd2d'.
        wv_calib (:class:`~pypeit.wavecalib.WaveCalib`, optional):
            This is only used for the :class:`SlicerIFUFindObjects` child when a joint sky subtraction
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
        manual (:class:`~pypeit.manual_extract.ManualExtractionObj`, optional):
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
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
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
        self.waveimg = None
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
        msgs.info("Initializing slits")
        self.initialize_slits(slits)

        # Internal bpm mask
        # We want to keep the 'BOXSLIT', which has bpm=2. But we don't want to keep 'BOXSLIT'
        # with other bad flag (for which bpm>2)
        # TODO: To my mind, we should never be using the value of the bit to
        # check for flags.  We should be using the BitMask functions.  I *think*
        # what you want is this:
        self.reduce_bpm = self.slits.bitmask.flagged(
                                self.slits.mask,
                                exclude='BOXSLIT',
                                and_not=self.slits.bitmask.exclude_for_reducing)
        # I.e., mask anything *except* slits flagged by only 'BOXSLIT', and also
        # make sure any of the `exclude_for_reducing` flags are not on.  This
        # previous code may also have included slits that were flagged as
        # SHORTSLIT.  Was that on purpose?
#        self.reduce_bpm = (self.slits.mask > 2) & (np.invert(self.slits.bitmask.flagged(
#                        self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing)))
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

    # TODO Make this a method possibly in slittrace.py. Almost identical code is in extraction.py
    def initialize_slits(self, slits, initial=False):
        """
        Gather all the :class:`~pypeit.slittrace.SlitTraceSet` attributes
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
        # TODO JFH: his is an ugly hack for the present moment until we get the slits object sorted out
        self.slits_left, self.slits_right, _ \
            = self.slits.select_edges(initial=initial, flexure=self.spat_flexure_shift)
        # This matches the logic below that is being applied to the slitmask. Better would be to clean up slits to
        # to return a new slits object with the desired selection criteria which would remove the ambiguity
        # about whether the slits and the slitmask are in sync.
        #bpm = self.slits.mask.astype(bool)
        #bpm &= np.invert(self.slits.bitmask.flagged(self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing + ['BOXSLIT']))
        #gpm = np.logical_not(bpm)
        #self.slits_left = slits_left[:, gpm]
        #self.slits_right = slits_right[:, gpm]


        # Slitmask
        self.slitmask = self.slits.slit_img(initial=initial, flexure=self.spat_flexure_shift,
                                            exclude_flag=self.slits.bitmask.exclude_for_reducing+['BOXSLIT'])
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)
#        # For echelle
#        self.spatial_coo = self.slits.spatial_coordinates(initial=initial, flexure=self.spat_flexure_shift)

    # TODO There are going to be problems with std_trace not being aligned with whatever orders are getting masked in
    # this routine.
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
        sobjs_obj, self.nobj = \
            self.find_objects(self.sciImg.image-initial_sky0, self.sciImg.ivar, std_trace=std_trace,
                              show_peaks=show_peaks, show=self.findobj_show and not self.std_redux)

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

    # TODO maybe we don't need parent and children for this method. But SlicerIFU has a bunch of extra methods.
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
                slit spat_id (:class:`MultiSlitFindObjects`, :class:`SlicerIFUFindObjects`)
                or ech_order (:class:`EchelleFindObjects`) value.

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        pass

    def global_skysub(self, skymask=None, bkg_redux_sciimg=None,
                      update_crmask=True,
                      previous_sky=None, show_fit=False, show=False, 
                      show_objs=False, objs_not_masked=False,
                      reinit_bpm:bool=True):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Args:
            skymask (`numpy.ndarray`_, optional):
                A 2D image indicating sky regions (1=sky)
            bkg_redux_sciimg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                PypeIt image of the science image before background subtraction
                if self.bkg_redux is True, otherwise None.
                It's used to generate a global sky  model without bkg subtraction.
            update_crmask (bool, optional):
                Update the crmask in the science image
            show_fit (bool, optional):
                Show the sky fits?
            show (bool, optional):
                Show the sky image generated?
            show_objs (bool, optional):
                If show=True, show the objects on the sky image?
            previous_sky (`numpy.ndarray`_, optional):
                Sky model estimate from a previous run of global_sky
                Used to generate an improved estimated of the variance
            objs_not_masked (bool, optional):
                Set this to be True if there are objects on the slit/order that are not being masked
                by the skymask. This is typically the case for the first pass sky-subtraction
                before object finding, since a skymask has not yet been created.
            reinit_bpm (:obj:`bool`, optional):
                If True (default), the bpm is reinitialized to the initial bpm 
                Should be False on the final run in case there was a failure
                upstream and no sources were found in the slit/order

        Returns:
            `numpy.ndarray`_: image of the the global sky model

        """
        # reset bpm since global sky is run several times and reduce_bpm is here updated.
        if reinit_bpm:
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
        tmp_bpm = self.slits.bitmask.flagged(self.slits.mask,
                                             and_not=self.slits.bitmask.exclude_for_reducing)
        gdslits = np.where(np.logical_not(tmp_bpm))[0]

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
            skysub_ivar = self.sciImg.ivar if bkg_redux_sciimg is None else bkg_redux_sciimg.ivar

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
            _image = self.sciImg.image if bkg_redux_sciimg is None else bkg_redux_sciimg.image
            global_sky[thismask] = skysub.global_skysub(
                _image, skysub_ivar, self.tilts, thismask,
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

            - This docstring is incomplete!

        Parameters
        ----------
        attr : str
            String specifying the image to show.  Options are:
                - global -- Sky model (global)
                - sci -- Processed science image
                - rawvar -- Raw variance image
                - modelvar -- Model variance image
                - crmasked -- Science image with CRs set to 0
                - skysub -- Science image with global sky subtracted
                - image -- Input image
        image : ndarray, optional
            User supplied image to display
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
                slit spat_id (MultiSlit, SlicerIFU) or ech_order (Echelle) value

        Returns:
            :obj:`float`: plate scale in binned pixels

        """
        _, bin_spat = parse.parse_binning(self.binning)
        return self.sciImg.detector.platescale * bin_spat

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
        specobjs : :class:`~pypeit.specobjs.SpecObjs`
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
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype and spectrograph.orders is not None \
                            else self.slits.ech_order
        if self.order_vec is None:
            msgs.error('Unable to set Echelle orders, likely because they were incorrectly '
                       'assigned in the relevant SlitTraceSet.')

    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector/echelle order

        Args:
            slitord_id (:obj:`int`, optional):
                slit spat_id (MultiSlit, SlicerIFU) or ech_order (Echelle) value

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
        specobjs : :class:`~pypeit.specobjs.SpecObjs`
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

        reduce_gpm = np.logical_not(self.reduce_bpm)
        sobjs_ech = findobj_skymask.ech_objfind(
            image, ivar, self.slitmask, self.slits_left[:, reduce_gpm], self.slits_right[:, reduce_gpm],
            self.slits.spat_id[reduce_gpm], self.order_vec[reduce_gpm],
            np.vstack((self.slits.specmin, self.slits.specmax))[:, reduce_gpm],
            det=self.det,
            inmask=inmask, 
            ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            manual_extract_dict=manual_extract_dict, 
            plate_scale=plate_scale[reduce_gpm],
            std_trace=std_trace,
            specobj_dict=specobj_dict,
            snr_thresh=self.par['reduce']['findobj']['snr_thresh'],
            show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
            fwhm=self.par['reduce']['findobj']['find_fwhm'],
            use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
            fof_link = self.par['reduce']['findobj']['fof_link'],
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


class SlicerIFUFindObjects(MultiSlitFindObjects):
    """
    Child of Reduce for SlicerIFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, spectrograph, par, objtype, **kwargs)

    def initialize_slits(self, slits, initial=True):
        """
        Gather all the :class:`~pypeit.slittrace.SlitTraceSet` attributes that
        we'll use here in :class:`FindObjects`. Identical to the parent but the
        slits are not trimmed.

        Args:
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                SlitTraceSet object containing the slit boundaries that will be
                initialized.
            initial (:obj:`bool`, optional):
                Use the initial definition of the slits. If False,
                tweaked slits are used.
        """
        super().initialize_slits(slits, initial=True)

    def global_skysub(self, skymask=None, bkg_redux_sciimg=None, update_crmask=True,
                      previous_sky=None, show_fit=False, show=False, show_objs=False, objs_not_masked=False,
                      reinit_bpm: bool = True):
        """
        Perform global sky subtraction. This SlicerIFU-specific routine ensures that the
        edges of the slits are not trimmed, and performs a spatial and spectral
        correction using the sky spectrum, if requested. See Reduce.global_skysub()
        for parameter definitions.

        See base class method for description of parameters.

        Args:
            reinit_bpm (:obj:`bool`, optional):
                If True (default), the bpm is reinitialized to the initial bpm
                Should be False on the final run in case there was a failure
                upstream and no sources were found in the slit/order
        """

        global_sky_sep = super().global_skysub(skymask=skymask, bkg_redux_sciimg=bkg_redux_sciimg, update_crmask=update_crmask,
                                               previous_sky=previous_sky, show_fit=show_fit, show=show,
                                               show_objs=show_objs,
                                               objs_not_masked=objs_not_masked, reinit_bpm=reinit_bpm)

        # Check if flexure or a joint fit is requested. If not return this
        if not self.par['reduce']['skysub']['joint_fit'] and self.par['flexure']['spec_method'] == 'skip':
            return global_sky_sep

        if self.wv_calib is None:
            msgs.error("A wavelength calibration is needed (wv_calib) if a joint sky fit is requested.")
        msgs.info("Generating wavelength image")

        # Generate the waveimg which is needed if flexure is being computed
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)

        # Calculate spectral flexure, so that we can align all slices of the IFU. We need to do this on the model
        # sky spectrum. It must be performed in this class if the joint sky fit is requested, because all of
        # the wavelengths need to be aligned for different slits before the sky is fit.
        method = self.par['flexure']['spec_method']
        # TODO :: Perhaps include a new label for IFU flexure correction - e.g. 'slitcen_relative' or 'slitcenIFU' or 'IFU'
        #      :: If a new label is introduced, change the other instances of 'method' (see below), and in flexure.spec_flexure_qa()
        if method in ['slitcen']:
            self.slitshift = self.calculate_flexure(global_sky_sep)
            # Recalculate the wavelength image, and the global sky taking into account the spectral flexure
            msgs.info("Generating wavelength image, accounting for spectral flexure")
            self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spec_flexure=self.slitshift,
                                                       spat_flexure=self.spat_flexure_shift)

        # If the joint fit or spec/spat sensitivity corrections are not being performed, return the separate slits sky
        if not self.par['reduce']['skysub']['joint_fit'] or bkg_redux_sciimg is not None:
            return global_sky_sep

        # If we reach this point in the code, a joint skysub has been requested.
        # Use sky information in all slits to perform a joint sky fit
        global_sky = self.joint_skysub(skymask=skymask, update_crmask=update_crmask,
                                       show_fit=show_fit, show=show, show_objs=show_objs,
                                       objs_not_masked=objs_not_masked)

        return global_sky

    def joint_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                     show_fit=False, show=False, show_objs=False, adderr=0.01, objs_not_masked=False):
        """ Perform a joint sky model fit to the data. See Reduce.global_skysub()
        for parameter definitions.
        """
        msgs.info("Performing joint global sky subtraction")
        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        _global_sky = np.zeros_like(self.sciImg.image)
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
                return _global_sky

        # Use the FWHM map determined from the arc lines to convert the science frame
        # to have the same effective spectral resolution.
        fwhm_map = self.wv_calib.build_fwhmimg(self.tilts, self.slits, initial=True, spat_flexure=self.spat_flexure_shift)
        thismask = thismask & (fwhm_map != 0.0)
        # Need to include S/N for deconvolution
        sciimg = skysub.convolve_skymodel(self.sciImg.image, fwhm_map, thismask)
        # Iterate to use a model variance image
        numiter = 4  # This is more than enough, and will probably break earlier than this
        model_ivar = self.sciImg.ivar
        sl_ref = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
        # Prepare the slitmasks for the relative spectral illumination
        slitmask = self.slits.slit_img(pad=0, initial=True, flexure=self.spat_flexure_shift)
        slitmask_trim = self.slits.slit_img(pad=-3, initial=True, flexure=self.spat_flexure_shift)
        for nn in range(numiter):
            msgs.info("Performing iterative joint sky subtraction - ITERATION {0:d}/{1:d}".format(nn+1, numiter))
            # TODO trim_edg is in the parset so it should be passed in here via trim_edg=tuple(self.par['reduce']['trim_edge']),
            _global_sky[thismask] = skysub.global_skysub(sciimg, model_ivar, tilt_wave,
                                                         thismask, self.slits_left, self.slits_right, inmask=inmask,
                                                         sigrej=sigrej, trim_edg=trim_edg,
                                                         bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                                         no_poly=self.par['reduce']['skysub']['no_poly'],
                                                         pos_mask=not self.bkg_redux and not objs_not_masked,
                                                         max_mask_frac=self.par['reduce']['skysub']['max_mask_frac'],
                                                         show_fit=show_fit)

            # Calculate the relative spectral illumination
            scaleImg = flat.illum_profile_spectral_poly(sciimg, self.waveimg, slitmask, slitmask_trim, _global_sky,
                                                        slit_illum_ref_idx=sl_ref, gpmask=inmask, thismask=thismask)
            # Apply this scale image to the temporary science frame
            sciimg /= scaleImg

            # Update the ivar image used in the sky fit
            msgs.info("Updating sky noise model")
            # Choose the highest counts out of sky and object
            counts = _global_sky
            _scale = None if self.sciImg.img_scale is None else self.sciImg.img_scale[thismask]
            # NOTE: darkcurr must be a float for the call below to work.
            if not self.bkg_redux:
                var = procimg.variance_model(self.sciImg.base_var[thismask], counts=counts[thismask],
                                             count_scale=_scale, noise_floor=adderr)
                model_ivar[thismask] = utils.inverse(var)
            else:
                model_ivar[thismask] = self.sciImg.ivar[thismask]
            # RJC :: Recalculating the global sky and flexure is probably overkill... but please keep this code in for now
            # Recalculate the sky on each individual slit and redetermine the spectral flexure
            # global_sky_sep = super().global_skysub(skymask=skymask, update_crmask=update_crmask,
            #                                        trim_edg=trim_edg, show_fit=show_fit, show=show,
            #                                        show_objs=show_objs)
            # self.calculate_flexure(global_sky_sep)

            # Check if the relative scaling isn't changing much after at least 4 iterations
            minv, maxv = np.min(scaleImg[thismask]), np.max(scaleImg[thismask])
            if nn >= 3 and max(abs(1/minv), abs(maxv)) < 1.005:  # Relative accuracy of 0.5% is sufficient
                break

        if update_crmask:
            # Find CRs with sky subtraction
            # NOTE: There's no need to run `sciImg.update_mask_cr` after this.
            # This operation updates the mask directly!
            self.sciImg.build_crmask(self.par['scienceframe']['process'], subtract_img=_global_sky)

        # Now we have a correct scale, apply it to the original science image
        self.apply_relative_scale(scaleImg)

        # Recalculate the joint sky using the original image
        _global_sky[thismask] = skysub.global_skysub(self.sciImg.image, model_ivar, tilt_wave,
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
        counts = _global_sky
        _scale = None if self.sciImg.img_scale is None else self.sciImg.img_scale[thismask]
        # NOTE: darkcurr must be a float for the call below to work.
        var = procimg.variance_model(self.sciImg.base_var[thismask], counts=counts[thismask],
                                     count_scale=_scale, noise_floor=adderr)
        model_ivar[thismask] = utils.inverse(var)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', global_sky=_global_sky, slits=True, sobjs=sobjs_show, clear=False)
        return _global_sky

    # TODO :: This function should be removed from the find_objects() class, once the flexure code has been tidied up.
    def calculate_flexure(self, global_sky):
        """
        Convenience function to calculate the flexure of an IFU. The flexure is calculated by cross-correlating the
        sky model of a reference slit with an archival sky spectrum. This gives an "absolute" flexure correction for
        the reference slit in pixels. Then, the flexure for all other slits is calculated by cross-correlating the
        sky model of each slit with the sky model of the reference slit. This gives a "relative" flexure correction
        for each slit in pixels. The relative flexure is then added to the absolute flexure to give the total flexure
        correction for each slit in pixels.

        Parameters
        ----------
        global_sky : ndarray
            Sky model

        Returns
        -------
        new_slitshift: ndarray
            The flexure in pixels
        """
        sl_ref = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
        box_rad = self.par['reduce']['extraction']['boxcar_radius']
        trace_spat = 0.5 * (self.slits_left + self.slits_right)
        iwv = np.where(self.wv_calib.spat_ids == self.slits.spat_id[sl_ref])[0][0]
        ref_fwhm_pix = self.wv_calib.wv_fits[iwv].fwhm
        # Extract a spectrum of the sky
        thismask = (self.slitmask == self.slits.spat_id[sl_ref])
        ref_skyspec = flexure.get_sky_spectrum(self.sciImg.image, self.sciImg.ivar, self.waveimg, thismask,
                                               global_sky, box_rad, self.slits, trace_spat[:, sl_ref],
                                               self.pypeline, self.det)
        # Calculate the flexure
        flex_dict_ref = flexure.spec_flex_shift(ref_skyspec, sky_file=self.par['flexure']['spectrum'], spec_fwhm_pix=ref_fwhm_pix,
                                            mxshft=self.par['flexure']['spec_maxshift'],
                                            excess_shft=self.par['flexure']['excessive_shift'],
                                            method="slitcen",
                                            minwave=self.par['flexure']['minwave'],
                                            maxwave=self.par['flexure']['maxwave'])
        this_slitshift = np.zeros(self.slits.nslits)
        if flex_dict_ref is not None:
            msgs.warn("Only a relative spectral flexure correction will be performed")
            this_slitshift = np.ones(self.slits.nslits) * flex_dict_ref['shift']
        # Now loop through all slits to calculate the additional shift relative to the reference slit
        flex_list = []
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            thismask = (self.slitmask == slit_spat)
            # Extract sky spectrum for this slit
            this_skyspec = flexure.get_sky_spectrum(self.sciImg.image, self.sciImg.ivar, self.waveimg, thismask,
                                                    global_sky, box_rad, self.slits, trace_spat[:, slit_idx],
                                                    self.pypeline, self.det)
            # Calculate the flexure
            flex_dict = flexure.spec_flex_shift(this_skyspec, arx_skyspec=ref_skyspec, arx_fwhm_pix=ref_fwhm_pix * 1.01,
                                                spec_fwhm_pix=ref_fwhm_pix,
                                                mxshft=self.par['flexure']['spec_maxshift'],
                                                excess_shft=self.par['flexure']['excessive_shift'],
                                                method="slitcen",
                                                minwave=self.par['flexure']['minwave'],
                                                maxwave=self.par['flexure']['maxwave'])
            this_slitshift[slit_idx] += flex_dict['shift']
            flex_list.append(flex_dict.copy())
        # Replace the reference slit with the absolute shift
        flex_list[sl_ref] = flex_dict_ref.copy()
        # Add this flexure to the previous flexure correction
        new_slitshift = self.slitshift + this_slitshift
        # Now report the flexure values
        for slit_idx, slit_spat in enumerate(self.slits.spat_id):
            msgs.info("Flexure correction, slit {0:d} (spat id={1:d}): {2:.3f} pixels".format(1+slit_idx, slit_spat,
                                                                                              self.slitshift[slit_idx]))
        # Save QA
        # TODO :: Need to implement QA once the flexure code has been tidied up, and this routine has been moved
        #         out of the find_objects() class.
        msgs.work("QA is not currently implemented for the flexure correction")
        if False:#flex_list is not None:
            basename = f'{self.basename}_global_{self.spectrograph.get_det_name(self.det)}'
            out_dir = os.path.join(self.par['rdx']['redux_path'], 'QA')
            slit_bpm = np.zeros(self.slits.nslits, dtype=bool)
            flexure.spec_flexure_qa(self.slits.slitord_id, slit_bpm, basename, flex_list, out_dir=out_dir)
        return new_slitshift

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
