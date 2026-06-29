"""
Main driver class for object finding, global skysubtraction and skymask construction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import inspect
import numpy as np
import os
from pathlib import Path

from abc import ABCMeta

from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d

from pypeit import specobj
from pypeit import specobjs
from pypeit import log, utils
from pypeit import PypeItError
from pypeit.display import display
from pypeit.core import skysub, qa, parse, flat, flexure
from pypeit.core import procimg
from pypeit.core import findobj_skymask
from pypeit.core import extract
from pypeit.core.moment import moment1d
from pypeit.core.fitting import iterfit
from pypeit.flatfield import FiberFlatImages, FlatImages
from pypeit.extraction import FiberExtract


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
            Currently only used with the IFU
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
        log.info("Initializing slits")
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
#        self.reduce_bpm = (self.slits.mask > 2) & (np.logical_not(self.slits.bitmask.flagged(
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
            raise PypeItError("Must provide either waveTilts or tilts to FindObjects")
        elif waveTilts is not None and tilts is not None:
            raise PypeItError("Cannot provide both waveTilts and tilts to FindObjects")
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
            log.info("Generating tilts image from fit in waveTilts")
            self.tilts = self.waveTilts.fit2tiltimg(self.slitmask, flexure=tilt_flexure_shift)
        elif waveTilts is None and tilts is not None:
            log.info("Using user input tilts image")
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
        # Instantiate the mask
        skymask = np.ones_like(self.sciImg.image, dtype=bool)
        if sobjs_obj.nobj == 0:
            # No objects found, so entire image contains sky
            return skymask

        # Build the mask for each slit
        boxcar_rad_pix = None
        gdslits = np.where(np.logical_not(self.reduce_bpm))[0]
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            log.info(f'Generating skymask for slit # {slit_spat}')
            thismask = self.slitmask == slit_spat
            this_sobjs = sobjs_obj.SLITID == slit_spat
            # Boxcar mask?
            if self.par['reduce']['skysub']['mask_by_boxcar']:
                boxcar_rad_pix = self.par['reduce']['extraction']['boxcar_radius'] / \
                                    self.get_platescale(slitord_id=self.slits.slitord_id[slit_idx])
            # Do it
            skymask[thismask] = findobj_skymask.create_skymask(
                                    sobjs_obj[this_sobjs], thismask, self.slits_left[:,slit_idx],
                                    self.slits_right[:,slit_idx], box_rad_pix=boxcar_rad_pix,
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
        #bpm &= np.logical_not(self.slits.bitmask.flagged(self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing + ['BOXSLIT']))
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
        std_trace : `astropy.table.Table`_, optional
            Table with the trace of the standard star on the input detector,
            which is used as a crutch for tracing. For MultiSlit reduction,
            the table has a single column: `TRACE_SPAT`.
            For Echelle reduction, the table has two columns: `ECH_ORDER` and `TRACE_SPAT`.
            The shape of each row must be (nspec,). For SlicerIFU reduction, std_trace is None.
            If None, the slit boundaries are used as the crutch.
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
            log.info("Skipping global sky sub as per user request")
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
            log.info('Either no objects were found or a user-provided sky mask was used.  '
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
            log.info("Skipping 2nd run of finding objects")
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
        std_trace : `astropy.table.Table`_, optional
            Table with the trace of the standard star on the input detector,
            which is used as a crutch for tracing. For MultiSlit reduction,
            the table has a single column: `TRACE_SPAT`.
            For Echelle reduction, the table has two columns: `ECH_ORDER` and `TRACE_SPAT`.
            The shape of each row must be (nspec,). For SlicerIFU reduction, std_trace is None.
            If None, the slit boundaries are used as the crutch.
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
            log.info("Finding objects in the negative image")
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
                log.info('Skipping global sky-subtraction for standard star.')
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
            log.info("Global sky subtraction for slit: {:d}".format(slit_spat))
            thismask = self.slitmask == slit_spat
            inmask = self.sciImg.select_flag(invert=True) & thismask & skymask_now
            # All masked?
            if not np.any(inmask):
                log.warning("No pixels for fitting sky.  If you are using mask_by_boxcar=True, your radius may be too large.")
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
                log.warning("Bad fit to sky.  Rejecting slit: {:d}".format(slit_spat))
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
            ch_name = chname if chname is not None else f'global_sky_{self.detname}'
            viewer, ch = display.show_image(image, chname=ch_name, mask=mask_in, clear=clear,
                                            wcs_match=True)
        elif attr == 'image':
            ch_name = chname if chname is not None else 'image'
            viewer, ch = display.show_image(image, chname=ch_name, clear=clear, wcs_match=True)
        else:
            log.warning("Not an option for show")

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

    def get_findobj_qa_filename(self, slit:int, neg:bool, save_objfindQA:bool) -> str:
        # Set objfind QA filename
        objfindQA_filename = None
        if save_objfindQA and (self.basename is not None):
            out_dir = os.path.join(self.par['rdx']['redux_path'], self.par['rdx']['qadir'])
            if self.find_negative:
                basename = 'neg_' + self.basename if neg else 'pos_' + self.basename
            else:
                basename = self.basename
            objfindQA_filename = qa.set_qa_filename(basename, 'obj_profile_qa', slit=slit,
                                                    det=self.detname, out_dir=out_dir)

        return objfindQA_filename

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
        std_trace : `astropy.table.Table`_, optional
            Table with the trace of the standard star on the input detector,
            which is used as a crutch for tracing. For MultiSlit reduction,
            the table has a single column: `TRACE_SPAT`.
            The shape of each row must be (nspec,). If None,
            the slit boundaries are used as the crutch.
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
        gdslits = np.where(np.logical_not(self.reduce_bpm))[0]

        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            qa_title ="Finding objects on slit # {:d}".format(slit_spat)
            log.info(qa_title)
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

            objfindQA_filename = self.get_findobj_qa_filename(slit_spat, neg, save_objfindQA)

            maxnumber =  self.par['reduce']['findobj']['maxnumber_std'] if self.std_redux \
                else self.par['reduce']['findobj']['maxnumber_sci']
            # standard star
            std_in = std_trace[0]['TRACE_SPAT'] if std_trace is not None else None

            # set the find_min_max and trace_min_max parameters
            find_minmax = [self.slits.specmin[slit_idx], self.slits.specmax[slit_idx]] if \
                self.par['reduce']['findobj']['find_min_max'] is None else \
                self.par['reduce']['findobj']['find_min_max']
            trace_minmax = [self.slits.specmin[slit_idx], self.slits.specmax[slit_idx]] if \
                self.par['reduce']['findobj']['trace_min_max'] is None else \
                self.par['reduce']['findobj']['trace_min_max']

            # Find objects
            sobjs_slit = \
                    findobj_skymask.objs_in_slit(image, ivar, thismask,
                                self.slits_left[:,slit_idx],
                                self.slits_right[:,slit_idx],
                                inmask=inmask,
                                ncoeff=self.par['reduce']['findobj']['trace_npoly'],
                                std_trace=std_in,
                                snr_thresh=snr_thresh,
                                hand_extract_dict=manual_extract_dict,
                                specobj_dict=specobj_dict, show_peaks=show_peaks,
                                show_fits=show_fits, show_trace=show_trace,
                                trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                                fwhm=self.par['reduce']['findobj']['find_fwhm'],
                                use_user_fwhm=self.par['reduce']['extraction']['use_user_fwhm'],
                                boxcar_rad=self.par['reduce']['extraction']['boxcar_radius'] / self.get_platescale(),  #pixels
                                maxshift=self.par['reduce']['findobj']['trace_maxshift'],
                                maxdev=self.par['reduce']['findobj']['trace_maxdev'],
                                numiterfit=self.par['reduce']['findobj']['find_numiterfit'],
                                find_min_max=find_minmax,
                                spec_min_max=trace_minmax,
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
            raise PypeItError('Unable to set Echelle orders, likely because they were incorrectly '
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
            raise PypeItError('slitord_id is missing. Plate scale for current echelle order cannot be determined.')
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
        std_trace : `astropy.table.Table`_, optional
            Table with the trace of the standard star on the input detector,
            which is used as a crutch for tracing. For Echelle reduction,
            the table has two columns: `ECH_ORDER` and `TRACE_SPAT`.
            The shape of each row must be (nspec,). If None,
            the slit boundaries are used as the crutch.
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

        objfindQA_filename = self.get_findobj_qa_filename(999, neg, save_objfindQA)

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
            maxshift=self.par['reduce']['findobj']['trace_maxshift'],
            maxdev=self.par['reduce']['findobj']['trace_maxdev'],
            numiterfit=self.par['reduce']['findobj']['find_numiterfit'],
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
            raise PypeItError("A wavelength calibration is needed (wv_calib) if a joint sky fit is requested.")
        log.info("Generating wavelength image")
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
            log.info("Generating wavelength image, accounting for spectral flexure")
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
        log.info("Performing joint global sky subtraction")
        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        _global_sky = np.zeros_like(self.sciImg.image)
        thismask = (self.slitmask != -1)
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
                log.info('Skipping global sky-subtraction for standard star.')
                return _global_sky

        # Use the FWHM map determined from the arc lines to convert the science frame
        # to have the same effective spectral resolution.
        fwhm_map = self.wv_calib.build_fwhmimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)
        thismask = thismask & (fwhm_map != 0.0)
        # Need to include S/N for deconvolution
        sciimg = skysub.convolve_skymodel(self.sciImg.image, fwhm_map, thismask)
        # Iterate to use a model variance image
        numiter = 4  # This is more than enough, and will probably break earlier than this
        model_ivar = self.sciImg.ivar
        sl_ref = self.par['calibrations']['flatfield']['slit_illum_ref_idx']
        # Prepare the slitmasks for the relative spectral illumination
        slitmask = self.slits.slit_img(pad=0, flexure=self.spat_flexure_shift)
        slitmask_trim = self.slits.slit_img(pad=-3, flexure=self.spat_flexure_shift)
        for nn in range(numiter):
            log.info("Performing iterative joint sky subtraction - ITERATION {0:d}/{1:d}".format(nn+1, numiter))
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
            log.info("Updating sky noise model")
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
        log.info("Updating sky noise model")
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
        box_rad = self.par['reduce']['extraction']['boxcar_radius']/ self.get_platescale()
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
            log.warning("Only a relative spectral flexure correction will be performed")
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
            log.info("Flexure correction, slit {0:d} (spat id={1:d}): {2:.3f} pixels".format(1+slit_idx, slit_spat,
                                                                                              self.slitshift[slit_idx]))
        # Save QA
        # TODO :: Need to implement QA once the flexure code has been tidied up, and this routine has been moved
        #         out of the find_objects() class.
        log.debug("QA is not currently implemented for the flexure correction")
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
        log.info("Correcting science frame for relative illumination")
        self.scaleimg *= scaleImg.copy()
        self.sciImg.image, _bpm, varImg = flat.flatfield(self.sciImg.image, scaleImg,
                                                         varframe=utils.inverse(self.sciImg.ivar))
        if np.any(_bpm):
            self.sciImg.update_mask('BADSCALE', indx=_bpm)
        self.sciImg.ivar = utils.inverse(varImg)


class FiberFindObjects(FindObjects):
    """
    Child of FindObjects for fiber-fed spectrographs.

    For fiber spectrographs, each fiber IS the object — there is no need
    for peak-detection object finding. This class creates one SpecObj per
    fiber within each block-slit, using reference fiber positions from
    the spectrograph.

    Follows the IDL Binospec pipeline approach (Fabricant et al. 2025,
    Section 6) for sky subtraction: extract all fibers, divide by the
    globally-normalized fiber flat correction, then fit
    a 2D B-spline sky model (wavelength + spatial Legendre polynomial)
    to the dedicated sky fibers and subtract from all fibers.

    See parent doc string for Args and Attributes.
    """
    def get_platescale(self, slitord_id=None):
        """
        Return the platescale in binned pixels for the current detector.

        Define get_platescale here (rather than inherit it) so the class is
        self-contained: the base :class:`FindObjects.get_platescale` is an
        unimplemented stub, and the inherited
        :meth:`FindObjects.create_skymask` calls this method.

        Args:
            slitord_id (:obj:`int`, optional):
                Unused; present for interface compatibility with the parent.

        Returns:
            :obj:`float`: Plate scale in binned pixels.
        """
        _, bin_spat = parse.parse_binning(self.binning)
        return self.sciImg.detector.platescale * bin_spat

    def run(self, std_trace=None, show_peaks=False, show_skysub_fit=False):
        """
        Primary code flow for fiber object finding and 2D sky model build.

        For fiber spectrographs each fiber IS the object -- no peak detection
        is needed.  We boxcar-pre-extract every fiber, fit one cross-detector
        2D bspline to the throughput-normalized sky-fiber sums, predict a 1D
        sky spectrum per fiber, and stamp it into a 2D ``SKYMODEL`` weighted
        by the empirical flat profile.  Standard ``FiberExtract`` then runs
        on ``sciimg - skymodel``.

        1. Build wavelength image (if not already available).
        2. Load ``FiberFlatImages`` for the per-fiber throughput scalars.
        3. Create one ``SpecObj`` per fiber via ``find_objects_pypeline``.
        4. Attach a per-fiber throughput scalar to each ``SpecObj``.
        5. Boxcar pre-extract every fiber from the un-subtracted image so
           the sky-model fit has aperture-summed counts and a wavelength
           solution to evaluate against.
        6. Fit a 2D bspline (wave + Legendre in normalized spatial position)
           to the throughput-normalized sky-fiber spectra, project per-fiber
           predictions, and distribute via the empirical flat profile to
           build the per-pixel 2D ``SKYMODEL``.

        Parameters
        ----------
        std_trace : `astropy.table.Table`_, optional
            Ignored for fiber reductions.
        show_peaks : :obj:`bool`, optional
            Ignored for fiber reductions.
        show_skysub_fit : :obj:`bool`, optional
            Ignored for fiber reductions.

        Returns
        -------
        initial_sky : `numpy.ndarray`_
            Per-pixel sky model in detector counts.
        sobjs_obj : :class:`~pypeit.specobjs.SpecObjs`
            List of objects found (one per fiber).
        """
        if self.waveimg is None:
            if self.wv_calib is None:
                raise PypeItError("Wavelength calibration required for "
                                  "fiber sky subtraction")
            self.waveimg = self.wv_calib.build_waveimg(
                self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)

        fiber_flatimages = self._load_fiber_flatimages()

        self.reduce_bpm = self.reduce_bpm_init.copy()
        # Fiber data never needs manual extraction or an A-B negative-image
        # pass, so skip the base find_objects wrapper (which only adds those)
        # and create the per-fiber SpecObjs directly.  find_objects_pypeline
        # ignores std_trace/show_peaks for fiber reductions.
        sobjs_obj, self.nobj = self.find_objects_pypeline(
            self.sciImg.image, self.sciImg.ivar, show=self.findobj_show)

        if len(sobjs_obj) == 0:
            return np.zeros_like(self.sciImg.image), sobjs_obj

        # Snap each trace to its actual fiber peak in the flat.  For Binospec
        # IFU, the bulk-shifted reference positions can still be off by 2-3 px
        # in individual blocks; without this refinement the end-fiber aperture
        # can miss its real peak and leave a large per-pixel residual at the
        # block edge.
        self._refine_fiber_traces(sobjs_obj)

        # Attach throughput scalar to each SpecObj so the sky-model build
        # and FiberExtract see one consistent per-fiber correction.
        self._attach_fiber_throughput(
            sobjs_obj, fiber_flatimages, self.spectrograph,
            self.sciImg.detector.det)

        # Boxcar pre-extract every fiber so we have BOX_WAVE + BOX_COUNTS
        # for the 2D sky bspline.  FiberExtract overwrites these from
        # ``sciimg - skymodel`` later.
        self._preextract_fibers(sobjs_obj)

        initial_sky = self._fiber_global_skysub(sobjs_obj)
        return initial_sky, sobjs_obj

    def _preextract_fibers(self, sobjs):
        """
        Boxcar-extract every fiber from the un-subtracted image.

        Populates ``BOX_WAVE``, ``BOX_COUNTS``, ``BOX_COUNTS_IVAR`` and
        ``BOX_MASK`` on each fiber SpecObj using a zero sky input, so the
        sky-model fit operates on aperture-summed sky counts in detector
        units.  These BOX_* fields are subsequently overwritten by
        :class:`~pypeit.extraction.FiberExtract` on ``sciimg - skymodel``.

        Fast path: single :func:`~pypeit.core.moment.moment1d` call per
        boxcar quantity with the per-fiber traces stacked into a
        ``(nspec, nfib)`` column array, replacing one round trip through
        :func:`~pypeit.core.extract.extract_boxcar` per fiber.
        Bit-identical to the legacy per-fiber loop, ~30x faster.

        Two preconditions guard the fast path, with a fallback to the
        per-fiber loop if either fails:

        * Uniform ``BOX_R_PIX`` across fibers, so a single scalar boxcar
          width can be used.  ``moment1d`` accepts a per-element ``width``
          array, but it sizes its integration column grid to the *minimum*
          aperture width across all elements (``moment.py``, the
          ``np.amin(i2-i1)`` term), so a varying ``width`` would silently
          truncate the wider apertures rather than reproduce the per-fiber
          boxcars.  The loop fallback uses correct per-fiber scalar widths.
        * No aperture crosses into a different block-slit, which would
          break the fast path's single global-mask assumption (never
          happens for Binospec's ~70 px inter-block gaps; see
          :meth:`_fiber_apertures_within_slit`).
        """
        if len(sobjs) == 0:
            return

        radii = np.array([float(s.BOX_R_PIX) for s in sobjs])
        if not np.allclose(radii, radii[0]):
            log.warning("Non-uniform BOX_R_PIX detected; "
                        "falling back to per-fiber boxcar pre-extract")
            self._preextract_fibers_loop(sobjs)
            return

        if not self._fiber_apertures_within_slit(sobjs):
            log.warning("A fiber aperture crosses a block-slit boundary; "
                        "falling back to per-fiber boxcar pre-extract")
            self._preextract_fibers_loop(sobjs)
            return

        gpm = self.sciImg.select_flag(invert=True)
        trace_2d = np.column_stack(
            [s.TRACE_SPAT.astype(float) for s in sobjs])
        box_width = 2.0 * float(radii[0])

        # Pre-mask images so moment1d zeroes masked pixels uniformly.
        # Pseudo-slit edges aren't hard apertures (they sit on the outermost
        # fiber, so the aperture extends into unassigned pixels in
        # inter-block gaps), and the runtime check above confirmed that
        # no aperture extends into pixels of a *different* slit -- so the
        # global gpm alone is equivalent to the per-fiber
        # ``gpm & ~other_slit`` mask the loop fallback uses.
        image_m = self.sciImg.image * gpm
        var_img = 1.0 / (self.sciImg.ivar + (self.sciImg.ivar == 0.0))
        var_m = var_img * gpm
        wave_m = self.waveimg * gpm
        wave_pos = (wave_m > 0.0).astype(float)

        flux = moment1d(image_m, trace_2d, box_width)[0]
        box_denom = moment1d(wave_pos, trace_2d, box_width)[0]
        wave_num = moment1d(wave_m, trace_2d, box_width)[0]
        var_box = moment1d(var_m, trace_2d, box_width)[0]
        pixtot = moment1d(np.ones_like(self.sciImg.image),
                          trace_2d, box_width)[0]
        pixmsk = moment1d(((self.sciImg.ivar == 0.0)
                           | ~gpm).astype(float),
                          trace_2d, box_width)[0]

        fwhmimg = None
        if self.wv_calib is not None:
            fwhmimg = self.wv_calib.build_fwhmimg(
                self.tilts, self.slits, initial=True,
                spat_flexure=self.spat_flexure_shift)
        fwhm_box = None

        wave_box = wave_num / (box_denom + (box_denom == 0.0))
        ivar_box = 1.0 / (var_box + (var_box == 0.0))
        # Compute the box mask on the *pre-replacement* wave_box so bad
        # rows (where wave_box would be zero) are correctly flagged.
        box_gpm = ((pixmsk != pixtot)
                   & np.isfinite(wave_box)
                   & (wave_box > 0))

        bad = ~np.isfinite(wave_box) | (wave_box <= 0)
        if bad.any():
            wave_pos_um = (self.waveimg > 0.0).astype(float)
            denom_um = moment1d(wave_pos_um, trace_2d, box_width)[0]
            wave_um = moment1d(self.waveimg, trace_2d, box_width)[0]
            wave_um_fixed = wave_um / (denom_um + (denom_um == 0.0))
            wave_box = np.where(bad, wave_um_fixed, wave_box)

        if fwhmimg is not None:
            fwhm_m = fwhmimg * gpm
            fwhm_sum = moment1d(fwhm_m, trace_2d, box_width)[0]
            n_good = pixtot - pixmsk
            fwhm_box = fwhm_sum * utils.inverse(n_good)
            for j in range(fwhm_box.shape[1]):
                fwhm_box[:, j] *= np.gradient(wave_box[:, j])

        # Match legacy extract_boxcar return semantics: zero out flux
        # and ivar at masked rows so low-throughput rows don't carry
        # spurious finite values into downstream fits.
        box_gpm_f = box_gpm.astype(float)
        flux = flux * box_gpm_f
        ivar_box = ivar_box * box_gpm_f

        for i, sobj in enumerate(sobjs):
            sobj.BOX_WAVE = wave_box[:, i].copy()
            sobj.BOX_COUNTS = flux[:, i].copy()
            sobj.BOX_COUNTS_IVAR = ivar_box[:, i].copy()
            sobj.BOX_MASK = box_gpm[:, i].copy()
            if fwhm_box is not None:
                sobj.BOX_FWHM = fwhm_box[:, i].copy()

    def _preextract_fibers_loop(self, sobjs):
        """
        Per-fiber fallback for :meth:`_preextract_fibers`.

        Preserves the original per-fiber ``gpm & ~other_slit`` mask and
        handles non-uniform ``BOX_R_PIX``.  Used when the vectorized fast
        path's preconditions fail.
        """

        gpm = self.sciImg.select_flag(invert=True)
        slitmask = self.slits.slit_img(initial=True)
        zero_sky = np.zeros_like(self.sciImg.image)

        fwhmimg = None
        if self.wv_calib is not None:
            fwhmimg = self.wv_calib.build_fwhmimg(
                self.tilts, self.slits, initial=True,
                spat_flexure=self.spat_flexure_shift)

        for sobj in sobjs:
            other_slit = (slitmask != sobj.SLITID) & (slitmask != -1)
            inmask = gpm & ~other_slit
            try:
                wave, flux, flux_ivar, _, _, box_gpm, fwhm, _, _, _, _ = \
                    extract.extract_boxcar(
                        sobj.BOX_R_PIX, sobj.TRACE_SPAT,
                        self.sciImg.image, self.sciImg.ivar,
                        inmask, self.waveimg, zero_sky,
                        fwhmimg=fwhmimg)
            except Exception as e:
                log.warning(f"Pre-extract failed for fiber "
                            f"{getattr(sobj, 'MASKDEF_ID', '?')}: {e}")
                continue
            sobj.BOX_WAVE = wave
            sobj.BOX_COUNTS = flux
            sobj.BOX_COUNTS_IVAR = flux_ivar
            sobj.BOX_MASK = box_gpm
            if fwhm is not None:
                sobj.BOX_FWHM = fwhm

    def _fiber_apertures_within_slit(self, sobjs):
        """
        Check that every fiber's boxcar aperture stays in own slit / gap.

        Returns ``True`` when, for every fiber and every spectral row, the
        integer aperture columns are either inside the fiber's own slit
        or in an unassigned pixel (``slitmask == -1``, i.e. inter-block
        gap).  If even one fiber crosses into a different slit, the
        vectorized boxcar's "global mask only" assumption breaks and the
        caller should fall back to the per-fiber loop.  For Binospec IFU
        with ~70 px inter-block gaps and ~4 px apertures this is always
        true; the check is the cheap guard that keeps the fast path safe
        for future fiber spectrographs.
        """
        nspec, nspat = self.sciImg.image.shape
        slitmask = self.slits.slit_img(initial=True)
        rows_full = np.arange(nspec)[:, None]
        for sobj in sobjs:
            trace = np.asarray(sobj.TRACE_SPAT, dtype=float)
            box_r = float(sobj.BOX_R_PIX)
            half = int(np.ceil(box_r))
            offsets = np.arange(-half, half + 1)
            cols_2d = (np.rint(trace).astype(int)[:, None]
                       + offsets[None, :])
            cols_clipped = np.clip(cols_2d, 0, nspat - 1)
            rows_2d = np.broadcast_to(rows_full, cols_2d.shape)
            in_aper = np.abs(cols_2d - trace[:, None]) <= box_r
            slit_at = slitmask[rows_2d, cols_clipped]
            if (in_aper
                    & (slit_at != sobj.SLITID)
                    & (slit_at != -1)).any():
                return False
        return True

    def _refine_fiber_traces(self, sobjs, search_radius=4.0):
        """
        Snap each fiber's trace to the nearest peak in the flat field.

        The static reference-profile fiber positions, even after the bulk
        cross-correlation shift, can be off by 2-3 px
        due to optical distortion or mechanical effects that a single global shift
        can't capture.  Left uncorrected, apertures can miss their real peaks and leave large
        per-pixel sky-subtraction residuals.

        For each fiber:

        * detect peaks in the flat-field spatial profile within its block-slit
        * pick the peak nearest the current ``SPAT_PIXPOS``
        * sub-pixel-refine that peak with a parabolic fit on the three
          flat samples around it
        * update ``TRACE_SPAT``, ``SPAT_PIXPOS``, and ``SPAT_PIXPOS_ID``

        Fibers whose nearest detected peak lies more than
        ``search_radius`` px away (dead fibers, masked regions, or
        anything pathological) are left at their bulk-shifted positions.

        Parameters
        ----------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            One SpecObj per fiber.
        search_radius : :obj:`float`, optional
            Maximum allowed shift in spatial pixels.  Defaults to 4.0 --
            comfortably larger than the worst per-block offsets observed
            (~3 px on DET02), but small enough to avoid snapping onto a
            neighbouring fiber's peak (inter-fiber spacing ~6.6 px).
        """

        if len(sobjs) == 0:
            return

        flatimg = self._load_flatimg()
        if flatimg is None:
            log.warning("No flat image for trace refinement; "
                        "using bulk-shifted reference positions")
            return

        nspec, nspat = flatimg.shape
        # Build a stable spatial profile by collapsing the central half of
        # the spectral axis; sky lines and CRs in the science image would
        # be a problem here but the flat is smooth in spectral.
        central = flatimg[nspec // 4:3 * nspec // 4, :]
        flat_strip = np.nanmedian(central, axis=0)

        # Group fibers by slit so peak detection only sees one block at a
        # time -- distance / height thresholds calibrate to the local
        # block's fibers without being confused by neighbours.
        by_slit = {}
        for sobj in sobjs:
            by_slit.setdefault(sobj.SLITID, []).append(sobj)

        n_refined = n_skipped = 0
        for slitid, fibers in by_slit.items():
            positions = np.array([s.SPAT_PIXPOS for s in fibers])
            box_r = float(fibers[0].BOX_R_PIX)
            lo = int(max(0, positions.min() - 2 * box_r - search_radius))
            hi = int(min(nspat, positions.max() + 2 * box_r + search_radius + 1))
            if hi - lo < 5:
                continue

            local_strip = flat_strip[lo:hi]
            finite_local = local_strip[np.isfinite(local_strip)]
            if finite_local.size == 0 or np.nanmax(local_strip) <= 0:
                continue

            # Minimum peak separation: 1.5 * box_r ~= 1.5 * half-spacing
            # is comfortably less than the inter-fiber spacing, so no two
            # real fiber peaks merge but spurious noise bumps are
            # rejected.  Height threshold 10% of the brightest peak.
            min_dist = max(2, int(box_r * 1.5))
            height = 0.1 * float(np.nanmax(local_strip))
            peaks_local, _ = find_peaks(local_strip, distance=min_dist,
                                        height=height)
            if peaks_local.size == 0:
                continue
            peaks = peaks_local + lo

            for sobj in fibers:
                ref_pos = float(sobj.SPAT_PIXPOS)
                dists = np.abs(peaks - ref_pos)
                k = int(np.argmin(dists))
                if dists[k] > search_radius:
                    n_skipped += 1
                    continue

                # Parabolic sub-pixel fit on the three flat samples around
                # the peak: y(x) = y1 - 0.5*(y2-y0)*dx + 0.5*(y0-2y1+y2)*dx^2
                # peak at dx = (y0 - y2) / (2 * (y0 - 2y1 + y2)).
                p_int = int(peaks[k])
                p_loc = p_int - lo
                if 0 < p_loc < len(local_strip) - 1:
                    y0 = float(local_strip[p_loc - 1])
                    y1 = float(local_strip[p_loc])
                    y2 = float(local_strip[p_loc + 1])
                    denom = y0 - 2.0 * y1 + y2
                    if denom < 0 and np.isfinite(denom):
                        sub_offset = 0.5 * (y0 - y2) / denom
                        # Reject runaway parabolic fits (e.g. when one of
                        # the flanking pixels is masked-fill).
                        if abs(sub_offset) <= 1.0:
                            new_center = float(p_int) + sub_offset
                        else:
                            new_center = float(p_int)
                    else:
                        new_center = float(p_int)
                else:
                    new_center = float(p_int)

                if abs(new_center - ref_pos) > search_radius:
                    n_skipped += 1
                    continue

                sobj.TRACE_SPAT = np.full(nspec, new_center, dtype=float)
                sobj.SPAT_PIXPOS = new_center
                sobj.SPAT_PIXPOS_ID = int(np.rint(new_center))
                n_refined += 1

        log.info(f"Trace refinement: snapped {n_refined}/{len(sobjs)} "
                 f"fibers to flat peaks ({n_skipped} skipped, "
                 f"search radius {search_radius} px)")

    @staticmethod
    def _attach_fiber_throughput(sobjs, fiber_flatimages, spectrograph, det):
        """
        Combine flat-derived and static per-fiber throughput.

        The flat field contributes the geometric/transmission ratio; the static
        ``fiber_illumination.fits`` vector (when the spectrograph supplies
        one) is a refinement that captures per-fiber variations the lamp
        cannot due to illumination effects.  Both are multiplied together and stashed on each SpecObj as
        a single scalar so the sky-model build and the post-extraction
        division apply the same value.
        """
        flat_lookup = {}
        if fiber_flatimages is not None \
                and getattr(fiber_flatimages, 'fiber_throughput', None) is not None:
            flat_lookup = {int(fid): float(t) for fid, t in
                           zip(fiber_flatimages.fiber_ids,
                               fiber_flatimages.fiber_throughput)}

        illum_lookup = {}
        if hasattr(spectrograph, 'load_fiber_illumination'):
            try:
                f_illum = spectrograph.load_fiber_illumination(det)
                ref = spectrograph.load_fiber_ref_profile(det)
                ref_ids = ref['FIB_ID']
                for i, fid in enumerate(ref_ids):
                    if i >= len(f_illum):
                        continue
                    val = float(f_illum[i])
                    if np.isfinite(val) and val > 0.1:
                        illum_lookup[int(fid)] = val
            except Exception as e:
                log.warning(f"Could not load fiber_illumination: {e}; "
                            f"using flat-only throughput")

        for sobj in sobjs:
            fid = sobj.MASKDEF_ID
            if fid is None or fid < 0:
                sobj.fiber_throughput = 1.0
                continue
            fid_int = int(fid)
            sobj.fiber_throughput = (flat_lookup.get(fid_int, 1.0)
                                     * illum_lookup.get(fid_int, 1.0))

    def _load_fiber_flatimages(self):
        """
        Load :class:`~pypeit.flatfield.FiberFlatImages` from the calibrations
        directory.

        Uses the ``calib_dir`` stored on the slits object (inherited from
        :class:`~pypeit.calibframe.CalibFrame`) to search for files matching
        the ``FiberFlat_*.fits`` pattern.

        Returns
        -------
        fiber_flatimages : :class:`~pypeit.flatfield.FiberFlatImages` or None
            Loaded fiber flat calibration, or ``None`` if no file is found.
        """

        calib_dir = getattr(self.slits, 'calib_dir', None)
        if calib_dir is None:
            log.warning("No calib_dir available on slits; "
                        "proceeding without fiber flat equalization")
            return None

        # Search for FiberFlat files matching the detector
        files = sorted(Path(calib_dir).glob(f'FiberFlat_*{self.detname}*.fits'))

        if files:
            log.info(f"Loading FiberFlatImages from {files[0].name}")
            return FiberFlatImages.from_file(files[0])

        log.warning(f"No FiberFlatImages found for {self.detname}; "
                     "proceeding without equalization")
        return None

    def _fiber_global_skysub(self, sobjs):
        """
        Fit one 2D bspline to throughput-normalized fiber spectra, predict
        per-fiber sky, and stamp it into a per-pixel 2D ``SKYMODEL``.

        We use the boxcar-extracted ``BOX_COUNTS`` (aperture sums per row),
        divide by each fiber's scalar throughput so all fibers share a
        common surface, and fit a single 2D bspline in
        ``(wavelength, normalized spat_pos)`` with a Legendre polynomial
        along the spatial axis (``npoly=2``).  Each fiber's predicted 1D
        sky spectrum is then evaluated, re-scaled by its throughput to
        recover detector counts, and the per-row sum is distributed across
        the fiber's aperture using the empirical spatial profile derived
        from the flat field.

        The fit runs as **two passes** with per-fiber wavelength refinement
        in between.  The per-block wavelength solution carries ~0.2 px
        (1 sigma) of per-fiber offset that the joint bspline can't resolve
        and that shows up as dipole residuals on bright OH lines.  After
        the first-pass fit, each fiber's ``BOX_WAVE`` is shifted to the
        Delta-lambda that minimises chi2 against the bspline model on the
        line wings (see :meth:`_refine_fiber_wavelengths`), then the joint
        fit is re-run.

        By default the fit is joint over dedicated sky fibers and all
        science fibers, relying on the iterative sigma rejection in
        :func:`~pypeit.core.fitting.iterfit` (asymmetric, tighter on the
        upper side) to down-weight pixels that contain real source flux.
        Two ``[reduce][skysub]`` parameters control which fibers
        contribute:

        * ``joint_fit_use_sci`` (default ``True``) -- when ``False``,
          only fibers whose ``MASKDEF_OBJNAME`` starts with ``SKY`` are
          used (legacy behaviour).
        * ``sci_exclude_radius`` (default ``None``) -- when set and the
          joint fit is enabled, any science fiber whose on-sky position
          lies within this radius (in arcsec) of the IFU geometric centre
          is dropped from the fit.  Use this to exclude inner science
          fibers most likely to contain source flux while retaining
          fibers away from the source for additional sky background
          information.

        Parameters
        ----------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            One SpecObj per fiber with ``BOX_WAVE``, ``BOX_COUNTS``,
            ``BOX_COUNTS_IVAR``, ``fiber_throughput`` and ``SPAT_PIXPOS``
            populated.  ``BOX_COUNTS_SKY`` is set on each SpecObj as a
            side effect; ``BOX_WAVE`` is updated in place with the
            per-fiber wavelength refinement.

        Returns
        -------
        sky_model : `numpy.ndarray`_
            Sky model in detector counts, shape ``self.sciImg.image.shape``.
        """

        nspec, nspat = self.sciImg.image.shape
        sky_model = np.zeros_like(self.sciImg.image)

        skysub_par = self.par['reduce']['skysub']
        # SkySubPar gained joint_fit_use_sci / sci_exclude_radius in this
        # change.  Guard against being driven by an older ParSet that
        # doesn't carry these keys.
        try:
            use_sci = bool(skysub_par['joint_fit_use_sci'])
        except KeyError:
            use_sci = True
        try:
            excl_r = float(skysub_par['sci_exclude_radius'] or 0.0)
        except KeyError:
            excl_r = 0.0

        sci_xy = self._get_sci_fiber_xy(sobjs) \
            if (use_sci and excl_r > 0) else {}

        bsp = skysub_par['bspline_spacing']
        if bsp is None:
            bsp = 1.2

        # First-pass fit on the as-extracted BOX_WAVE so we have a reference
        # bspline model to cross-correlate each fiber against.
        fit1 = self._fit_fiber_sky_bspline(
            sobjs, use_sci, excl_r, sci_xy, bsp, label='first pass')
        if fit1 is None:
            return sky_model
        sset1, wave_min, wave_max, _ = fit1

        # Per-fiber wavelength refinement: ~0.2 px (1 sigma) of per-block
        # solution scatter shows up as dipole residuals on sharp OH lines.
        # Snap each fiber's BOX_WAVE to the Delta-lambda that minimises
        # chi2 against sset1 on the line wings.
        self._refine_fiber_wavelengths(sobjs, sset1, wave_min, wave_max)

        # Second-pass fit with refined wavelengths.  Falls back to sset1
        # if the re-fit fails for any reason.
        fit2 = self._fit_fiber_sky_bspline(
            sobjs, use_sci, excl_r, sci_xy, bsp, label='second pass')
        if fit2 is not None:
            sset, wave_min, wave_max, _ = fit2
        else:
            log.warning("Second-pass fiber sky bspline failed; "
                        "predicting with first-pass model")
            sset = sset1

        # Predict per-fiber 1D sky spectrum in detector counts (aperture sums
        # per row) and stash it on each SpecObj.
        for sobj in sobjs:
            if sobj.BOX_WAVE is None:
                continue
            thru = float(getattr(sobj, 'fiber_throughput', 1.0))
            if not np.isfinite(thru) or thru <= 0:
                thru = 1.0
            in_range = ((sobj.BOX_WAVE >= wave_min)
                        & (sobj.BOX_WAVE <= wave_max))
            sky_spec = np.zeros_like(sobj.BOX_WAVE, dtype=float)
            if np.any(in_range):
                x2_arr = np.full(int(np.sum(in_range)),
                                 float(sobj.SPAT_PIXPOS) / nspat)
                vals, _ = sset.value(sobj.BOX_WAVE[in_range], x2=x2_arr)
                sky_spec[in_range] = vals * thru
            sobj.BOX_COUNTS_SKY = sky_spec

        # Per-fiber LSF broadening (IDL-style).  Binospec IFU sky fibers have a larger
        # output aperture and a measurably broader spectral LSF than
        # science fibers (~0.10 px sigma gap from the arc-line FWHM on
        # Binospec).  The joint bspline pools all fibers and dominantly
        # sees the sci LSF (sci fibers outnumber sky 8:1), so subtracting
        # the un-broadened model leaves negative cores at the broader
        # (sky) fibers.  Mitigate by Gaussian-broadening each fiber's
        # predicted 1D sky to match its own arc-line LSF before stamping.
        # Reference sigma is the sci-fiber median; broadening sigma is
        # ``sqrt(sig_fiber^2 - sig_ref^2)`` and only applied above a
        # 0.2 px threshold.  Models the IDL pipeline
        # ``bino_sub_sky_ms.pro`` per-row broadening step.
        if self._measure_fiber_lsf(sobjs, sset, wave_min, wave_max):
            self._broaden_sky_per_fiber(sobjs)

        # Distribute each fiber's per-row sky total across its aperture using
        # the empirical flat profile (unit sum per row).  Profiles are built
        # per-fiber inside the loop and discarded immediately afterwards;
        # caching a dense (nspec, nspat) profile per fiber would cost ~134 MB
        # per fiber on a 4k detector.  Falls back to a uniform aperture stamp
        # if the flat image isn't available or no fiber profile succeeds.
        flatimg = self._load_flatimg()
        if flatimg is None:
            log.warning("No flat image for empirical profiles; using uniform "
                        "per-aperture sky stamps")
            return self._stamp_skymodel_uniform(sobjs)

        slitmask = self.slits.slit_img(initial=True)
        n_stamped = 0
        for sobj in sobjs:
            sky_spec = sobj.BOX_COUNTS_SKY
            if sky_spec is None:
                continue
            prof = FiberExtract._build_empirical_profile(
                flatimg, slitmask, sobj, nspec, nspat)
            if prof is None:
                continue
            good = np.isfinite(sky_spec)
            if not np.any(good):
                del prof
                continue
            sky_model[good, :] += prof[good, :] * sky_spec[good, None]
            n_stamped += 1
            del prof

        if n_stamped == 0:
            log.warning("No empirical flat profiles built; using uniform "
                        "per-aperture sky stamps")
            return self._stamp_skymodel_uniform(sobjs)

        log.info(f"Stamped 2D sky model from empirical profiles for "
                 f"{n_stamped}/{len(sobjs)} fibers")
        return sky_model

    def _fit_fiber_sky_bspline(self, sobjs, use_sci, excl_r, sci_xy, bsp,
                                label=''):
        """
        Collect throughput-normalised fiber spectra and fit one 2D bspline.

        Used by :meth:`_fiber_global_skysub` for both the first-pass fit on
        as-extracted ``BOX_WAVE`` and the second-pass fit on per-fiber
        refined wavelengths.

        Returns
        -------
        result : tuple or None
            ``(sset, wave_min, wave_max, outmask)`` on success, or ``None``
            if there are no good fiber inputs or ``iterfit`` raises.
        """

        nspec, nspat = self.sciImg.image.shape
        all_wave, all_flux, all_ivar, all_x2 = [], [], [], []
        n_sky_fibers = n_sci_fibers = n_excluded = 0
        for sobj in sobjs:
            name = sobj.MASKDEF_OBJNAME
            is_sky = name is not None and str(name).upper().startswith('SKY')
            if not is_sky:
                if not use_sci:
                    continue
                if excl_r > 0 and sci_xy:
                    xy = sci_xy.get(int(sobj.MASKDEF_ID)) \
                        if sobj.MASKDEF_ID is not None else None
                    if xy is not None and float(np.hypot(*xy)) < excl_r:
                        n_excluded += 1
                        continue
            if sobj.BOX_COUNTS is None or sobj.BOX_WAVE is None:
                continue
            thru = float(getattr(sobj, 'fiber_throughput', 1.0))
            if not np.isfinite(thru) or thru <= 0:
                thru = 1.0
            good = ((sobj.BOX_WAVE > 0)
                    & np.isfinite(sobj.BOX_COUNTS)
                    & np.isfinite(sobj.BOX_COUNTS_IVAR)
                    & (sobj.BOX_COUNTS_IVAR > 0))
            if sobj.BOX_MASK is not None:
                good &= sobj.BOX_MASK
            if not np.any(good):
                continue
            n = int(np.sum(good))
            all_wave.append(sobj.BOX_WAVE[good])
            all_flux.append(sobj.BOX_COUNTS[good] / thru)
            all_ivar.append(sobj.BOX_COUNTS_IVAR[good] * thru**2)
            all_x2.append(np.full(n, float(sobj.SPAT_PIXPOS) / nspat))
            if is_sky:
                n_sky_fibers += 1
            else:
                n_sci_fibers += 1

        if not all_wave:
            log.warning("No fiber data for sky fit; returning zero sky model")
            return None

        all_wave = np.concatenate(all_wave)
        all_flux = np.concatenate(all_flux)
        all_ivar = np.concatenate(all_ivar)
        all_x2 = np.concatenate(all_x2)
        srt = np.argsort(all_wave)
        all_wave = all_wave[srt]
        all_flux = all_flux[srt]
        all_ivar = all_ivar[srt]
        all_x2 = all_x2[srt]

        # When science fibers participate the spatial axis is densely
        # sampled, so use a higher-order Legendre polynomial (cubic) to
        # let the bspline capture the smooth spatial throughput variation
        # the static fiber_illumination scalar misses.  Sky-only mode keeps
        # npoly=2 (matches the IDL Binospec IFU pipeline).
        # Symmetric sigma rejection on the joint fit: the previously
        # asymmetric (upper=3, lower=5) rejection was found to bias the
        # bspline model upward at bright OH-line cores, leaving negative
        # residuals (model > data) for both sky and sci fibers.  Source
        # flux still gets rejected at upper=3.
        if n_sci_fibers > 0:
            npoly = 4
            upper, lower = 3.0, 3.0
        else:
            npoly = 2
            upper, lower = 3.0, 3.0

        tag = f" ({label})" if label else ''
        log.info(
            f"Fiber sky bspline{tag}: {len(all_wave)} pts from "
            f"{n_sky_fibers} sky + {n_sci_fibers} sci fibers "
            f"(excluded {n_excluded} sci within r<{excl_r}\"), "
            f"bsp={bsp}, npoly={npoly}, sigrej=({lower},{upper})")
        try:
            sset, outmask = iterfit(
                all_wave, all_flux, invvar=all_ivar, x2=all_x2,
                upper=upper, lower=lower, maxiter=10, nord=4,
                kwargs_bspline={'bkspace': bsp, 'npoly': npoly})
        except Exception as e:
            log.warning(f"Fiber sky bspline fit{tag} failed: {e}")
            return None

        n_rej = int(np.sum(~outmask))
        log.info(f"Fiber sky bspline{tag}: {n_rej}/{len(outmask)} pixels "
                 f"rejected ({100 * n_rej / max(len(outmask), 1):.1f}%)")
        return sset, float(all_wave[0]), float(all_wave[-1]), outmask

    def _measure_fiber_lsf(self, sobjs, sset, wave_min, wave_max, n_lines=5):
        """
        Per-fiber spectral LSF sigma (in pixels) from the arc-line FWHM.

        Uses ``sobj.BOX_FWHM`` (wavecal-derived spectral FWHM in Angstroms,
        already populated along each fiber's trace).  For each fiber,
        samples ``BOX_FWHM`` at the brightest sky lines in the joint
        bspline model and converts FWHM_AA -> sigma_pix via the local
        dispersion ``median(diff(BOX_WAVE))``.

        Per-fiber empirical estimators (2nd-moment / Gaussian fit on
        ``BOX_COUNTS``) were tried first and discarded -- they are
        biased high by OH-line wings/blending and have 40-50% within-type
        scatter, swamping the real ~10% per-type LSF gap.  ``BOX_FWHM``
        is per-block-slit (wavecal granularity), which collapses this to
        a per-type measurement -- adequate because the sky/sci LSF gap
        comes from the fiber aperture size , which is uniform within a fiber type.

        Returns
        -------
        n_measured : int
            Number of fibers for which a usable LSF sigma was measured and
            stored on ``sobj.fiber_lsf_sigma_pix`` (pixels).  Fibers without
            a usable value are left with ``fiber_lsf_sigma_pix = None``.
        """

        n_wave_search = 4000
        grid_w = np.linspace(wave_min + 5.0, wave_max - 5.0, n_wave_search)
        try:
            mdl, _ = sset.value(grid_w, x2=np.full_like(grid_w, 0.5))
        except Exception as e:
            log.warning(f"Could not evaluate bspline for LSF line search: {e}")
            return 0
        finite = np.isfinite(mdl)
        if not finite.any():
            return 0
        base = float(np.nanmedian(mdl[finite]))
        noise = 1.4826 * float(np.nanmedian(np.abs(mdl[finite] - base)))
        if not np.isfinite(noise) or noise <= 0:
            return 0
        peaks, _ = find_peaks(mdl, height=base + 25.0 * noise, distance=20)
        if peaks.size == 0:
            log.warning("No bright lines found in bspline model for LSF "
                        "measurement; skipping per-fiber sky broadening")
            return 0
        order = np.argsort(mdl[peaks])[::-1]
        line_waves = grid_w[peaks[order][:n_lines]]

        n_measured = 0
        for sobj in sobjs:
            fwhm = getattr(sobj, 'BOX_FWHM', None)
            wave = sobj.BOX_WAVE
            if fwhm is None or wave is None:
                continue
            mask = sobj.BOX_MASK if sobj.BOX_MASK is not None \
                else np.ones_like(wave, dtype=bool)
            good = (mask & np.isfinite(wave) & np.isfinite(fwhm)
                    & (wave > 0) & (fwhm > 0))
            if good.sum() < 100:
                continue
            wave_g, fwhm_g = wave[good], fwhm[good]
            srt = np.argsort(wave_g)
            wave_g, fwhm_g = wave_g[srt], fwhm_g[srt]
            dlam = float(np.median(np.diff(wave_g)))
            if not np.isfinite(dlam) or dlam <= 0:
                continue
            w_lo, w_hi = float(wave_g[0]), float(wave_g[-1])
            sigs = []
            for lw in line_waves:
                if lw < w_lo or lw > w_hi:
                    continue
                fwhm_aa = float(np.interp(lw, wave_g, fwhm_g))
                sig_px = (fwhm_aa / 2.355) / dlam
                if 0.3 <= sig_px <= 3.0:
                    sigs.append(sig_px)
            if sigs:
                sobj.fiber_lsf_sigma_pix = float(np.median(sigs))
                n_measured += 1
        return n_measured

    def _broaden_sky_per_fiber(self, sobjs, threshold_pix=0.2):
        """
        Gaussian-convolve each fiber's predicted 1D sky to its observed LSF.

        Reads the per-fiber LSF sigma stored on ``sobj.fiber_lsf_sigma_pix``
        by :meth:`_measure_fiber_lsf`.  Reference sigma is the median over
        **science** fibers (the dominant population in the joint bspline fit,
        so the model sits closest to their LSF).  Each fiber whose own sigma
        exceeds the reference is broadened by ``sqrt(sig_fiber^2 - sig_ref^2)``
        in row-index space (one row = one spectral pixel), only applied above
        ``threshold_pix`` (default 0.2 px, matches the IDL pipeline).
        ``scipy.ndimage.gaussian_filter1d`` is flux-preserving so per-fiber
        total sky counts are conserved.
        """

        all_sigs = []
        sci_sigs = []
        for sobj in sobjs:
            sig = getattr(sobj, 'fiber_lsf_sigma_pix', None)
            if sig is None:
                continue
            all_sigs.append(sig)
            name = sobj.MASKDEF_OBJNAME
            if name is None or not str(name).upper().startswith('SKY'):
                sci_sigs.append(sig)
        if sci_sigs:
            sig_ref = float(np.median(sci_sigs))
            ref_src = f"median of {len(sci_sigs)} sci fibers"
        else:
            sig_ref = float(np.min(all_sigs))
            ref_src = f"min of {len(all_sigs)} fibers (no sci available)"

        n_broad = 0
        deltas = []
        for sobj in sobjs:
            sig = getattr(sobj, 'fiber_lsf_sigma_pix', None)
            if sig is None or sobj.BOX_COUNTS_SKY is None:
                continue
            diff_sq = sig ** 2 - sig_ref ** 2
            if diff_sq <= 0:
                deltas.append(0.0)
                continue
            sig_broad = float(np.sqrt(diff_sq))
            deltas.append(sig_broad)
            if sig_broad < threshold_pix:
                continue
            sobj.BOX_COUNTS_SKY = gaussian_filter1d(
                sobj.BOX_COUNTS_SKY, sigma=sig_broad, mode='nearest')
            n_broad += 1
        arr = np.asarray(deltas)
        if arr.size:
            log.info(
                f"Per-fiber LSF broadening: ref_sig={sig_ref:.3f} px "
                f"({ref_src}); broadened {n_broad}/{arr.size} fibers "
                f"(threshold {threshold_pix} px); delta_sig "
                f"median={float(np.median(arr)):.3f} px "
                f"max={float(arr.max()):.3f} px")

    def _refine_fiber_wavelengths(self, sobjs, sset, wave_min, wave_max,
                                  search_pix=1.5, probe_pix=0.05,
                                  relinearize_pix=0.3,
                                  min_npts=100):
        """
        Per-fiber wavelength refinement against the first-pass bspline.

        For spectrographs with many fibers, performing a full wavelength
        calibration for each fiber is very time-consuming. As an optimization,
        fibers that are closely spaced are grouped into "blocks" and share
        a common wavelength solution derived from the arc lines. Tilt-fitting
        is used to capture the intra-block variation, but the per-block solution still
        leaves ~0.2 px (1 sigma) of
        per-fiber offset within a block.  Pooling fibers with mismatched
        wavelengths smears the bspline at sharp OH lines and shows up as
        a dipole residual on each fiber's bright lines.

        For each fiber, locally linearize the bspline in pixel-shift
        ``Delta-p`` via a one-sided finite difference of width
        ``probe_pix``: ``dbase/dp = (base(lambda + probe*dlam) - base) /
        probe``.  Minimizing the weighted chi2

            chi2(Delta-p) = sum_i w_i^2 (r_i - Delta-p * dbase/dp_i)^2

        in closed form gives

            Delta-p = sum w_i^2 dbase/dp_i r_i / sum w_i^2 (dbase/dp_i)^2,

        where ``w_i^2 = (dbase/dp_i)^2 * ivar_i`` -- the same
        gradient-weighted norm the legacy grid scan used (continuum rows
        contribute nothing; line wings dominate; source flux on sci fibers
        averages out because dbase/dp swaps sign across each line peak).
        If the first-pass shift exceeds ``relinearize_pix`` the
        linearization is repeated centered on the new estimate; one Newton
        step is enough for the largest shifts seen in practice (<= 0.7 px
        with knot spacing ~0.8 px).  The whole solve costs 2--4 bspline
        evaluations per fiber instead of the 60+ required by a grid scan.

        Parameters
        ----------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            Fiber SpecObjs from the first pass; ``BOX_WAVE`` is updated in
            place.
        sset : :class:`~pypeit.bspline.bspline.bspline`
            First-pass bspline.
        wave_min, wave_max : float
            Valid wavelength range for the bspline.
        search_pix : float, optional
            Maximum allowed ``|Delta-p|``; the solve is clipped to this.
        probe_pix : float, optional
            Pixel offset used for the finite-difference derivative.
        relinearize_pix : float, optional
            If ``|Delta-p|`` from the first solve exceeds this, re-linearize
            around the new estimate and re-solve.
        min_npts : int, optional
            Skip fibers with fewer good rows than this.
        """
        shifts = []
        n_skip_box = n_skip_npts = n_skip_dlam = n_skip_signal = 0
        for sobj in sobjs:
            if sobj.BOX_WAVE is None or sobj.BOX_COUNTS is None:
                n_skip_box += 1
                continue
            thru = float(getattr(sobj, 'fiber_throughput', 1.0))
            if not np.isfinite(thru) or thru <= 0:
                thru = 1.0
            good = ((sobj.BOX_WAVE > wave_min)
                    & (sobj.BOX_WAVE < wave_max)
                    & np.isfinite(sobj.BOX_COUNTS)
                    & np.isfinite(sobj.BOX_COUNTS_IVAR)
                    & (sobj.BOX_COUNTS_IVAR > 0))
            if sobj.BOX_MASK is not None:
                good &= sobj.BOX_MASK
            if int(np.sum(good)) < min_npts:
                n_skip_npts += 1
                continue

            wave = sobj.BOX_WAVE[good].astype(float)
            data = sobj.BOX_COUNTS[good].astype(float) / thru
            ivar = sobj.BOX_COUNTS_IVAR[good].astype(float) * thru**2
            x2 = np.full(int(good.sum()),
                         float(sobj.SPAT_PIXPOS) /
                         self.sciImg.image.shape[1])

            dlam = float(np.median(np.diff(np.sort(wave))))
            if not np.isfinite(dlam) or dlam <= 0:
                n_skip_dlam += 1
                continue

            delta = self._solve_fiber_shift(
                sset, wave, data, ivar, x2, dlam, probe_pix, center=0.0)
            if delta is None:
                n_skip_signal += 1
                continue

            # The linearization is around delta=0; if the first-pass
            # estimate moves the center far enough that local curvature
            # matters, re-linearize there and add the Newton correction.
            if abs(delta) > relinearize_pix:
                delta2 = self._solve_fiber_shift(
                    sset, wave, data, ivar, x2, dlam, probe_pix,
                    center=delta)
                if delta2 is not None:
                    delta = delta + delta2

            delta = float(np.clip(delta, -search_pix, search_pix))
            shift_AA = delta * dlam
            sobj.BOX_WAVE = sobj.BOX_WAVE + shift_AA
            # Cache the shift so FiberExtract can re-apply it after
            # extract_boxcar/extract_optimal overwrite BOX_WAVE/OPT_WAVE
            # from the global wavelength image -- the extracted spectrum
            # then shares the same per-fiber wavelength solution that the
            # joint sky model was built against.
            sobj.wave_refine_shift_AA = shift_AA
            shifts.append(delta)

        skipped = (f"skipped: box={n_skip_box} npts={n_skip_npts} "
                   f"dlam={n_skip_dlam} signal={n_skip_signal}")
        if shifts:
            arr = np.array(shifts)
            log.info(
                f"Fiber wavelength refinement: N={arr.size} refined / "
                f"{len(sobjs)} total ({skipped}); "
                f"median={np.median(arr):+.3f}px "
                f"mean={np.mean(arr):+.3f}px "
                f"std={np.std(arr):.3f}px "
                f"range=({arr.min():+.3f}, {arr.max():+.3f})px")
        else:
            log.warning(f"Fiber wavelength refinement: no fibers refined "
                        f"({skipped})")

    @staticmethod
    def _solve_fiber_shift(sset, wave, data, ivar, x2, dlam, probe_pix,
                           center):
        """
        Closed-form WLS solve for the per-fiber Delta-pixel shift.

        Linearizes the bspline around ``center`` (in pixel units) via a
        one-sided finite-difference derivative of width ``probe_pix``, then
        returns the Delta-p that minimizes
        ``sum w_i^2 (r_i - Delta-p * dbase/dp_i)^2`` with weights
        ``w_i^2 = (dbase/dp_i)^2 * ivar_i``.  Returns ``None`` if the model
        has no gradient on this fiber (e.g. continuum-only after masking)
        or the system is otherwise ill-conditioned.

        Parameters
        ----------
        sset : :class:`~pypeit.bspline.bspline.bspline`
            The reference (joint) bspline being refined against.
        wave : :class:`numpy.ndarray`
            Per-row wavelength of the fiber's boxcar spectrum (good rows
            only).
        data : :class:`numpy.ndarray`
            Throughput-corrected per-row counts.
        ivar : :class:`numpy.ndarray`
            Throughput-corrected per-row inverse variance.
        x2 : :class:`numpy.ndarray`
            Normalized spatial coordinate fed to ``sset.value`` for the
            Legendre dependence.
        dlam : float
            Local wavelength dispersion (A / pixel).
        probe_pix : float
            Finite-difference step (in pixels) used to estimate
            ``dbase/dp``.
        center : float
            Pixel offset (relative to the input ``wave``) at which to
            linearize.

        Returns
        -------
        float or None
            The closed-form Delta-p that adds to ``center`` to give the
            estimated shift, or ``None`` if the solve was degenerate.
        """
        base, _ = sset.value(wave + center * dlam, x2=x2)
        base_plus, _ = sset.value(
            wave + (center + probe_pix) * dlam, x2=x2)
        dbase = (base_plus - base) / probe_pix
        # weight^2 = (dbase/dp)^2 * ivar; matches the legacy
        # grad*sqrt(ivar) weighting (the 1/dlam^2 conversion factor cancels
        # in numer/denom so we drop it).
        w2 = (dbase ** 2) * np.where(ivar > 0, ivar, 0.0)
        denom = float(np.sum(w2 * dbase ** 2))
        if not np.isfinite(denom) or denom <= 0:
            return None
        numer = float(np.sum(w2 * dbase * (data - base)))
        delta = numer / denom
        return float(delta) if np.isfinite(delta) else None

    def _stamp_skymodel_uniform(self, sobjs):
        """
        Build a 2D SKYMODEL by stamping per-row sky totals uniformly.

        Fallback when an empirical flat profile is unavailable.  For each
        fiber, every pixel in the row's aperture gets ``sky_spec[row] /
        n_aperture_pixels``; the per-row sum still equals ``sky_spec[row]``
        so boxcar extraction recovers the correct sky.
        """
        nspec, nspat = self.sciImg.image.shape
        sky_model = np.zeros_like(self.sciImg.image)
        slitmask = self.slits.slit_img(
            pad=0, flexure=self.spat_flexure_shift)
        spat_pix = np.arange(nspat)
        for sobj in sobjs:
            sky_spec = sobj.BOX_COUNTS_SKY
            if sky_spec is None or sobj.TRACE_SPAT is None \
                    or sobj.BOX_R_PIX is None:
                continue
            rows = np.where(np.isfinite(sky_spec)
                            & np.isfinite(sobj.TRACE_SPAT))[0]
            for row in rows:
                # Allow aperture into unassigned pixels (inter-block gap),
                # never into a different slit.  Matches the relaxed
                # aperture used in extraction.
                other_slit = ((slitmask[row, :] != sobj.SLITID)
                              & (slitmask[row, :] != -1))
                in_aperture = ((np.abs(spat_pix - sobj.TRACE_SPAT[row])
                                <= sobj.BOX_R_PIX)
                               & ~other_slit)
                cols = np.where(in_aperture)[0]
                if cols.size == 0:
                    continue
                pix_sky = float(sky_spec[row]) / cols.size
                if not np.isfinite(pix_sky):
                    continue
                sky_model[row, cols] = pix_sky
        return sky_model

    def _get_sci_fiber_xy(self, sobjs):
        """
        Look up on-sky (x, y) in arcsec for each science fiber.

        Used by :meth:`_fiber_global_skysub` to optionally exclude science
        fibers within a configured radius of an IFU geometric center from
        the joint sky fit.

        Relies on the spectrograph exposing ``load_sky_layout()`` (returning
        ``targetx, targety`` arrays in arcsec) and
        ``get_science_fiber_layout_indices(det, fiber_ids, fiber_types)``.
        Spectrographs that do not implement these get an empty mapping and
        no exclusion is applied.

        Parameters
        ----------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            Fiber SpecObjs with ``MASKDEF_ID`` and ``MASKDEF_OBJNAME`` set.

        Returns
        -------
        :obj:`dict`
            Mapping ``MASKDEF_ID -> (x_arcsec, y_arcsec)`` for science
            fibers that resolve to a layout entry.  Empty if the lookup is
            not supported.
        """
        spec = self.spectrograph
        if not (hasattr(spec, 'load_sky_layout')
                and hasattr(spec, 'get_science_fiber_layout_indices')):
            return {}

        try:
            targetx, targety = spec.load_sky_layout()
        except Exception as e:
            log.warning(f"Could not load IFU sky layout for "
                        f"sci_exclude_radius: {e}")
            return {}

        det = self.sciImg.detector.det
        fiber_ids = np.array(
            [int(s.MASKDEF_ID) if s.MASKDEF_ID is not None else -1
             for s in sobjs], dtype=int)
        fiber_types = np.array(
            ['SKY' if (s.MASKDEF_OBJNAME is not None
                       and str(s.MASKDEF_OBJNAME).upper().startswith('SKY'))
             else 'SCI' for s in sobjs])

        try:
            layout_idx = spec.get_science_fiber_layout_indices(
                det, fiber_ids, fiber_types)
        except Exception as e:
            log.warning(f"Could not map fibers to layout indices for "
                        f"sci_exclude_radius: {e}")
            return {}

        xy = {}
        for sobj, lidx in zip(sobjs, layout_idx):
            if lidx is None or lidx < 0:
                continue
            if sobj.MASKDEF_ID is None:
                continue
            xy[int(sobj.MASKDEF_ID)] = (float(targetx[lidx]),
                                        float(targety[lidx]))
        return xy

    def _load_flatimg(self):
        """
        Load the processed flat image used to build empirical fiber profiles.

        Returns
        -------
        flatimg : `numpy.ndarray`_ or None
            Processed flat image in detector-pixel space, or ``None`` if no
            matching flat calibration is available.
        """

        calib_dir = getattr(self.slits, 'calib_dir', None)
        if calib_dir is None:
            return None

        files = sorted(Path(calib_dir).glob(f'Flat_*{self.detname}*.fits'))
        if not files:
            log.warning(f"No FlatImages found for {self.detname}; "
                        "using uniform per-aperture sky stamps")
            return None

        log.info(f"Loading FlatImages from {files[0].name} for fiber profiles")
        flat_images = FlatImages.from_file(files[0])
        flat_raw = flat_images.pixelflat_raw
        if flat_raw is None:
            return None

        if flat_images.pixelflat_norm is None:
            return flat_raw.copy()

        flatimg, _ = flat.flatfield(flat_raw, flat_images.pixelflat_norm)
        return flatimg

    def find_objects_pypeline(self, image, ivar, std_trace=None,
                              manual_extract_dict=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, save_objfindQA=False, neg=False, debug=False):
        """
        Create one SpecObj per fiber within each block-slit.

        For each pseudo-slit, uses fiber positions from the spectrograph's
        reference profile to create SpecObjs at known fiber locations.
        FWHM and BOX_R_PIX are set from inter-fiber spacing so adjacent
        boxcar apertures touch but don't overlap.

        Parameters
        ----------
        image : `numpy.ndarray`_
            Image to search for objects from. Shape ``(nspec, nspat)``.
        ivar : `numpy.ndarray`_
            Inverse variance of ``image``.
        std_trace : `astropy.table.Table`_, optional
            Ignored for fiber reductions.
        manual_extract_dict : :obj:`dict`, optional
            Ignored for fiber reductions.
        show_peaks : :obj:`bool`, optional
            Ignored for fiber reductions.
        show_fits : :obj:`bool`, optional
            Ignored for fiber reductions.
        show_trace : :obj:`bool`, optional
            Ignored for fiber reductions.
        show : :obj:`bool`, optional
            Show QA plots.
        save_objfindQA : :obj:`bool`, optional
            Ignored for fiber reductions.
        neg : :obj:`bool`, optional
            Ignored for fiber reductions.
        debug : :obj:`bool`, optional
            Ignored for fiber reductions.

        Returns
        -------
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            Container holding one SpecObj per fiber.
        nobj : :obj:`int`
            Number of objects identified.
        """
        gdslits = np.where(np.logical_not(self.reduce_bpm))[0]
        sobjs = specobjs.SpecObjs()
        nspec = image.shape[0]
        obj_counter = 0

        # Get block structure from spectrograph
        blocks = self.spectrograph.get_fiber_blocks(self.det)
        fiber_shift = self.spectrograph.get_fiber_position_shift(
            self.slits, self.det) if hasattr(self.spectrograph, 'get_fiber_position_shift') else 0.0

        for slit_idx in gdslits:
            slit_spat_id = self.slits.spat_id[slit_idx]
            left = self.slits_left[:, slit_idx]
            right = self.slits_right[:, slit_idx]

            # Get reference fiber positions for this block
            if slit_idx >= len(blocks):
                continue
            block = blocks[slit_idx]
            fiber_centers = block['fiber_positions'] + fiber_shift

            if len(fiber_centers) == 0:
                continue

            # Identify fibers using spectrograph reference
            fiber_meta = self.spectrograph.identify_fibers_in_block(
                self.det, slit_idx, fiber_centers)

            # Compute inter-fiber spacings for BOX_R_PIX
            spacings = np.diff(fiber_centers)
            half_spacings = np.zeros(len(fiber_centers))
            if len(spacings) > 0:
                half_spacings[0] = spacings[0] / 2.0
                half_spacings[-1] = spacings[-1] / 2.0
                half_spacings[1:-1] = np.minimum(spacings[:-1], spacings[1:]) / 2.0
            else:
                half_spacings[0] = np.median(right - left) / 2.0

            for j, center_pix in enumerate(fiber_centers):
                obj_counter += 1
                trace_center = np.full(nspec, center_pix)

                thisobj = specobj.SpecObj(
                    PYPELINE='Fiber',
                    DET=self.sciImg.detector.name,
                    OBJTYPE=self.objtype,
                    SLITID=slit_spat_id,
                )
                thisobj.TRACE_SPAT = trace_center.astype(float)
                thisobj.trace_spec = np.arange(nspec)
                thisobj.SPAT_PIXPOS = float(center_pix)
                thisobj.SPAT_PIXPOS_ID = int(np.rint(center_pix))
                # Median slit/block width, floored to 1 px to avoid division
                # by zero on a degenerate or mis-traced slit (left == right).
                slit_width = max(float(np.median(right - left)), 1.0)
                thisobj.SPAT_FRACPOS = (center_pix - np.median(left)) / \
                    slit_width
                thisobj.FWHM = 2.0 * half_spacings[j]
                thisobj.maskwidth = half_spacings[j] / (slit_width / 2.0)
                thisobj.BOX_R_PIX = half_spacings[j]
                thisobj.smash_peakflux = 1.0
                thisobj.smash_snr = 100.0
                thisobj.OBJID = obj_counter

                # Assign fiber metadata
                if fiber_meta is not None:
                    thisobj.MASKDEF_ID = int(fiber_meta['fiber_id'][j])
                    thisobj.MASKDEF_OBJNAME = fiber_meta['fiber_name'][j]

                thisobj.set_name()
                sobjs.add_sobj(thisobj)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            gpm = self.sciImg.select_flag(invert=True)
            self.show('image', image=image*gpm.astype(float),
                      chname='objfind', sobjs=sobjs, slits=True)

        return sobjs, len(sobjs)
