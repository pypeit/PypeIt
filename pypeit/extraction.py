"""
Main driver class for skysubtraction and extraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import inspect
import numpy as np
import os

from abc import ABCMeta

from pypeit import log, utils
from pypeit import PypeItError
from pypeit.display import display
from pypeit.core import skysub, flexure, flat
from pypeit.core.moment import moment1d


class Extract:
    """
    This class will organize and run actions relatedt to sky subtraction, and extraction for
    a Science or Standard star exposure


    Attributes:
        ivarmodel (`numpy.ndarray`_):
            Model of inverse variance
        objimage (`numpy.ndarray`_):
            Model of object
        skyimage (`numpy.ndarray`_):
            Final model of sky
        global_sky (`numpy.ndarray`_):
            Fit to global sky
        outmask (:class:`~pypeit.images.imagebitmask.ImageBitMaskArray`):
            Final output mask
        extractmask (`numpy.ndarray`_):
            Extraction mask
        slits (:class:`~pypeit.slittrace.SlitTraceSet`):
        sobjs_obj (:class:`~pypeit.specobjs.SpecObjs`):
            Only object finding but no extraction
        sobjs (:class:`~pypeit.specobjs.SpecObjs`):
            Final extracted object list with trace corrections applied
        spat_flexure_shift (float):
        tilts (`numpy.ndarray`_):
            WaveTilts images generated on-the-spot
        waveimg (`numpy.ndarray`_):
            WaveImage image generated on-the-spot
        slitshift (`numpy.ndarray`_):
            Global spectral flexure correction for each slit (in pixels)
        vel_corr (float):
            Relativistic reference frame velocity correction (e.g. heliocentyric/barycentric/topocentric)
        extract_bpm (`numpy.ndarray`_):
            Bad pixel mask for extraction

    """

    __metaclass__ = ABCMeta

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, sciImg, slits, sobjs_obj, spectrograph, par, objtype, global_sky=None,
                     bkg_redux_global_sky=None, waveTilts=None, tilts=None, wv_calib=None, waveimg=None,
                     flatimages=None, bkg_redux=False, return_negative=False, std_redux=False, show=False, basename=None):
        """
        Instantiate the Extract subclass appropriate for the provided
        spectrograph.

        The class must be subclassed from Extract.  See :class:`Extract` for
        the description of the valid keyword arguments.

        Args:
            sciImg (:class:`~pypeit.images.pypeitimage.PypeItImage`):
                Image to reduce.
            slits (:class:`~pypeit.slittrace.SlitTraceSet`):
                Slit trace set object
            sobjs_obj (:class:`~pypeit.specobjs.SpecObjs`):
                Objects found but not yet extracted
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            par (:class:`~pypeit.par.pypeitpar.PypeItPar`):
                Parameter set for Extract
            objtype (:obj:`str`):
                Specifies object being reduced 'science', 'standard', or
                'science_coadd2d'.  This is used only to determine the
                spat_flexure_shift and ech_order for coadd2d.
            global_sky (`numpy.ndarray`_, optional):
                Fit to global sky. If None, an array of zeroes is generated the
                same size as ``sciImg``.
            bkg_redux_global_sky (`numpy.ndarray`_, optional):
                Sky estimate without background subtraction. This is used for 1d sky spectrum extraction
                in the case bkg_redux=True. Default is None.
            waveTilts (:class:`~pypeit.wavetilts.WaveTilts`, optional):
                This is waveTilts object which is optional, but either waveTilts or tilts must
                be provided.
            tilts (`numpy.ndarray`_, optional):
                Tilts image. Either a tilts image or waveTilts object (above) must be provided.
            wv_calib (:class:`~pypeit.wavecalib.WaveCalib`, optional):
                This is the waveCalib object which is optional, but either wv_calib or waveimg must be provided.
            waveimg (`numpy.ndarray`_, optional):
                Wave image. Either a wave image or wv_calib object (above) must be provided
            flatimages (:class:`~pypeit.flatfield.FlatImages`, optional):
                FlatImages class. This is optional, but if provided, it is used to extract the
                normalized blaze profile.
            bkg_redux (:obj:`bool`, optional):
                If True, the sciImg has been subtracted by
                a background image (e.g. standard treatment in the IR)
            return_negative (:obj:`bool`, optional):
                If True, negative objects from difference imaging will also be
                extracted and returned. Default=False.  This option only applies
                to the case where bkg_redux=True, i.e. typically a near-IR
                reduction where difference imaging has been employed to perform
                a first-pass at sky-subtraction. The default behavior is to not
                extract these objects, although they are masked in global
                sky-subtraction (performed in the find_objects class), and
                modeled in local sky-subtraction (performed by this class).
            std_redux (:obj:`bool`, optional):
                If True the object being extracted is a standards star
                so that the reduction parameters can be adjusted accordingly.
            basename (str, optional):
                Output filename used for spectral flexure QA
            show (:obj:`bool`, optional):
                Show plots along the way?

        Returns:
            :class:`~pypeit.extraction.Extract`:  Extraction object.
        """
        return next(c for c in utils.all_subclasses(Extract)
                    if c.__name__ == (spectrograph.pypeline + 'Extract'))(
            sciImg, slits, sobjs_obj, spectrograph, par, objtype, global_sky=global_sky,
            bkg_redux_global_sky=bkg_redux_global_sky, waveTilts=waveTilts, tilts=tilts,
            wv_calib=wv_calib, waveimg=waveimg, flatimages=flatimages, bkg_redux=bkg_redux,
            return_negative=return_negative, std_redux=std_redux, show=show, basename=basename)

    def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, global_sky=None,
                 bkg_redux_global_sky=None, waveTilts=None, tilts=None, wv_calib=None, waveimg=None,
                 flatimages=None, bkg_redux=False, return_negative=False, std_redux=False, show=False,
                 basename=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.sobjs_obj = sobjs_obj.copy() # This guarantees that stuff below does not mute this input variable
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        self.global_sky = global_sky if global_sky is not None else np.zeros_like(sciImg.image)
        self.bkg_redux_global_sky = bkg_redux_global_sky

        self.basename = basename
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

        # Initialize the slits
        log.info("Initializing slits")
        self.initialize_slits(slits)

        # Internal bpm mask
        self.extract_bpm = self.slits.bitmask.flagged(
                                self.slits.mask,
                                and_not=self.slits.bitmask.exclude_for_reducing)
        self.extract_bpm_init = self.extract_bpm.copy()

        # These may be None (i.e. COADD2D)
        #self.waveTilts = caliBrate.wavetilts
        #self.wv_calib = caliBrate.wv_calib

        # Load up other input items
        self.bkg_redux = bkg_redux
        self.return_negative=return_negative

        self.std_redux = std_redux
        self.det = sciImg.detector.det
        self.binning = sciImg.detector.binning
        self.pypeline = spectrograph.pypeline
        self.extract_show = show

        self.steps = []

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objmodel = None
        self.skymodel = None
        self.bkg_redux_skymodel = None
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs = None  # Final extracted object list with trace corrections applied
        self.slitshift = np.zeros(self.slits.nslits)  # Global spectral flexure slit shifts (in pixels) that are applied to all slits.
        self.vel_corr = None

        # Deal with dynamically generated calibrations, i.e. the tilts. If the tilts are not input generate
        # them from the  fits in caliBrate, otherwise use the input tilts
        if waveTilts is None and tilts is None:
            raise PypeItError("Must provide either waveTilts or tilts to Extract")
        elif waveTilts is not None and tilts is not None:
            raise PypeItError("Cannot provide both waveTilts and tilts to Extract")
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

        # Now generate the wavelength image
        log.info("Generating wavelength image")
        if wv_calib is None and waveimg is None:
            raise PypeItError("Must provide either wv_calib or waveimg to Extract")
        if wv_calib is not None and waveimg is not None:
            raise PypeItError("Cannot provide both wv_calib and waveimg to Extract")
        if wv_calib is not None and waveimg is None:
            self.wv_calib = wv_calib
            self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)
        elif wv_calib is None and waveimg is not None:
            self.waveimg = waveimg

        log.info("Generating spectral FWHM image")
        self.fwhmimg = None
        if wv_calib is not None:
            self.fwhmimg = wv_calib.build_fwhmimg(self.tilts, self.slits, initial=True, spat_flexure=self.spat_flexure_shift)
        else:
            log.warning("Spectral FWHM image could not be generated")

        # get flatfield image for blaze function
        self.flatimg = None
        if flatimages is not None:
            # This way of getting the flat image (instead of just reading flatimages.pixelflat_model) ensures that the
            # flat image is available also when an archival pixel flat is used.
            flat_raw = flatimages.pixelflat_raw if flatimages.pixelflat_raw is not None else flatimages.illumflat_raw
            if flat_raw is not None and flatimages.pixelflat_norm is not None:
                # TODO: Can we just use flat_raw if flatimages.pixelflat_norm is None?
                self.flatimg, _ = flat.flatfield(flat_raw, flatimages.pixelflat_norm)
        if self.flatimg is None:
            log.warning("No flat image was found. A spectrum of the flatfield will not be extracted!")

        # Now apply a global flexure correction to each slit provided it's not a standard star
        if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
            # Update slitshift values
            self.spec_flexure_correct(mode='global')
            # Apply?
            for iobj in range(self.sobjs_obj.nobj):
                # Ignore negative sources
                if self.sobjs_obj[iobj].sign < 0 and not self.return_negative:
                    continue
                islit = self.slits.spatid_to_zero(self.sobjs_obj[iobj].SLITID)
                self.sobjs_obj[iobj].update_flex_shift(self.slitshift[islit], flex_type='global')

        # Remove objects rejected for one or another reasons (we don't want to extract those)
        remove_idx = []
        for i, sobj in enumerate(self.sobjs_obj):
            if sobj.SLITID in list(self.slits.spat_id[self.extract_bpm]):
                remove_idx.append(i)
        # remove
        self.sobjs_obj.remove_sobj(remove_idx)


    @property
    def nsobj_to_extract(self):
        """
        Number of sobj objects in sobjs_obj taking into account whether or not we are returning negative traces

        Returns:

        """

        if len(self.sobjs_obj) > 0:
            return len(self.sobjs_obj) if self.return_negative else np.sum(self.sobjs_obj.sign > 0)
        else:
            return 0

    # TODO Make this a method possibly in slittrace.py. Almost identical code is in find_objects.py
    def initialize_slits(self, slits, initial=False):
        """
        Gather all the :class:`~pypeit.slittrace.SlitTraceSet` attributes
        that we'll use here in :class:`Extract`

        Args
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
        #bpm &= np.logical_not(self.slits.bitmask.flagged(self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing))
        #gpm = np.logical_not(bpm)
        #self.slits_left = slits_left[:, gpm]
        #self.slits_right = slits_right[:, gpm]

        # Slitmask
        self.slitmask = self.slits.slit_img(initial=initial, flexure=self.spat_flexure_shift,
                                            exclude_flag=self.slits.bitmask.exclude_for_reducing)
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)

    def extract(self, global_sky, bkg_redux_global_sky=None, model_noise=None, spat_pix=None):
        """
        Main method to extract spectra from the ScienceImage

        Args:
            global_sky (`numpy.ndarray`_):
                Sky estimate
            sobjs_obj (:class:`~pypeit.specobjs.SpecObjs`):
                List of SpecObj that have been found and traced
            bkg_redux_global_sky (`numpy.ndarray`_):
                Sky estimate without background subtraction. This is used for 1d sky spectrum extraction
                in the case bkg_redux=True. Default is None.
            model_noise (bool):
                If True, construct and iteratively update a model inverse variance image
                using :func:`~pypeit.core.procimg.variance_model`. If False, a
                variance model will not be created and instead the input sciivar will
                always be taken to be the inverse variance. See :func:`~pypeit.core.skysub.local_skysub_extract`
                for more info. Default is None, which is to say pypeit will use the bkg_redux attribute to
                decide whether or not to model the noise.
            spat_pix (`numpy.ndarray`_):
                 Image containing the spatial coordinates. This option is used for 2d coadds
                 where the spat_pix image is generated as a coadd of images. For normal reductions
                 spat_pix is not required as it is trivially created from the image itself. Default is None.

        """
        # This holds the objects, pre-extraction
        # JFH Commenting this out. Not sure why we need this. It overwrites the previous stuff from the init
        #self.sobjs_obj = sobjs_obj

        if self.par['reduce']['extraction']['skip_optimal']:  # Boxcar only with global sky subtraction
            log.info("Skipping optimal extraction")

            # This will hold the extracted objects
            self.sobjs = self.sobjs_obj.copy()

            # Quick loop over the objects
            for sobj in self.sobjs:
                # True  = Good, False = Bad for inmask
                thismask = self.slitmask == sobj.SLITID  # pixels for this slit
                inmask = self.sciImg.select_flag(invert=True) & thismask
                extract_sky = global_sky if bkg_redux_global_sky is None else bkg_redux_global_sky
                # Do it
                sobj.extract_boxcar(
                    self.sciImg.image-global_sky, self.sciImg.ivar, inmask, self.waveimg,
                    extract_sky, fwhmimg=self.fwhmimg, base_var=self.sciImg.base_var,
                    count_scale=self.sciImg.img_scale, noise_floor=self.sciImg.noise_floor
                )

            # Fill up extra bits and pieces
            self.objmodel = np.zeros_like(self.sciImg.image)
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # NOTE: fullmask is a bit mask, make sure it's treated as such, not
            # a boolean (e.g., bad pixel) mask.
            self.outmask = self.sciImg.fullmask.copy()
            self.skymodel = global_sky.copy()
            self.bkg_redux_skymodel = bkg_redux_global_sky.copy() if bkg_redux_global_sky is not None else None

        else:  # Local sky subtraction and optimal extraction.
            model_noise_1 = not self.bkg_redux if model_noise is None else model_noise
            self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
                self.local_skysub_extract(global_sky, self.sobjs_obj,
                                          bkg_redux_global_sky=bkg_redux_global_sky,
                                          model_noise=model_noise_1,
                                          spat_pix = spat_pix,
                                          show_profile=self.extract_show,
                                          show=self.extract_show)

        # Remove sobjs that don't have either OPT_COUNTS or BOX_COUNTS
        remove_idx = []
        for idx, sobj in enumerate(self.sobjs):
            # Find them
            if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
                remove_idx.append(idx)
                log.warning(f'Removing object at pixel {sobj.SPAT_PIXPOS} because '
                          f'both optimal and boxcar extraction could not be performed')
            elif sobj.OPT_COUNTS is None:
                log.warning(f'Optimal extraction could not be performed for object at pixel {sobj.SPAT_PIXPOS}')

        # Remove them
        if len(remove_idx) > 0:
            self.sobjs.remove_sobj(remove_idx)

        # Add the DETECTOR container, S/N ratio, and FWHM in ARCSEC and BOX_R_ASEC for each extracted object
        for sobj in self.sobjs:
            # this is an internal attribute
            sobj.spectrograph = self.spectrograph
            # these are datamodel attributes
            sobj.PYP_SPEC = self.spectrograph.name
            sobj.DETECTOR = self.sciImg.detector
            sobj.S2N = sobj.med_s2n()
            sobj.SPAT_FWHM = sobj.med_fwhm()
            sobj.BOX_R_ASEC = sobj.boxcar_arcsec()
            if self.sciImg.det_img is not None:
                _det_spec = moment1d(self.sciImg.det_img, sobj.TRACE_SPAT, 1, row=sobj.trace_spec)[0].astype(int)
                sobj.SPEC_DET = np.clip(_det_spec, 0, len(sobj.DETECTOR.detectors))

        # Return
        return self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    def run(self, model_noise=None, spat_pix=None):
        """
        Primary code flow for PypeIt reductions

        Args:
            model_noise (bool):
                If True, construct and iteratively update a model inverse variance image
                using :func:`~pypeit.core.procimg.variance_model`. If False, a
                variance model will not be created and instead the input sciivar will
                always be taken to be the inverse variance. See :func:`~pypeit.core.skysub.local_skysub_extract`
                for more info. Default is None, which is to say pypeit will use the bkg_redux attribute to
                decide whether or not to model the noise.
            spat_pix (`numpy.ndarray`_):
                 Image containing the spatial coordinates. This option is used for 2d coadds
                 where the spat_pix image is generated as a coadd of images. For normal reductions
                 spat_pix is not required as it is trivially created from the image itself. Default is None.

        Returns:
            tuple: skymodel (ndarray), bkg_redux_skymodel (ndarray), objmodel (ndarray), ivarmodel (ndarray),
               outmask (ndarray), sobjs (SpecObjs), waveimg (`numpy.ndarray`_),
               tilts (`numpy.ndarray`_), slits (:class:`~pypeit.slittrace.SlitTraceSet`).
               See main doc string for description

        """
        # Do we have any detected objects to extract?
        if self.nsobj_to_extract > 0:
            # Extract + Return
            self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = self.extract(self.global_sky, bkg_redux_global_sky=self.bkg_redux_global_sky,
                               model_noise=model_noise, spat_pix=spat_pix)
            if self.bkg_redux:
                # purge negative objects if not return_negative otherwise keep them
                self.sobjs.make_neg_pos() if self.return_negative else self.sobjs.purge_neg()

            # Correct for local spectral flexure
            if self.par['flexure']['spec_method'] not in ['skip', 'slitcen'] and not self.std_redux:
                # Apply a refined estimate of the flexure to objects
                self.spec_flexure_correct(mode='local', sobjs=self.sobjs)

        else:  # No objects, pass back what we have
            # Could have negative objects but no positive objects so purge them if not return_negative
            if self.bkg_redux:
                self.sobjs_obj.make_neg_pos() if self.return_negative else self.sobjs_obj.purge_neg()
            self.skymodel = self.global_sky
            self.bkg_redux_skymodel = self.bkg_redux_global_sky
            self.objmodel = np.zeros_like(self.sciImg.image)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # Set to the initial mask in case no objects were found
            # NOTE: fullmask is a bit mask, make sure it's treated as such, not
            # a boolean (e.g., bad pixel) mask.
            self.outmask = self.sciImg.fullmask.copy()
            # empty specobjs object from object finding
            self.sobjs = self.sobjs_obj

        # Update the mask
        # TODO: change slits.mask > 2 to use named flags.
        reduce_masked = np.where(np.logical_not(self.extract_bpm_init) & self.extract_bpm & (self.slits.mask > 2))[0]
        if len(reduce_masked) > 0:
            self.slits.mask[reduce_masked] = self.slits.bitmask.turn_on(
                self.slits.mask[reduce_masked], 'BADEXTRACT')

        # Return
        return self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, \
            self.outmask, self.sobjs, self.waveimg, self.tilts, self.slits

    def local_skysub_extract(self, global_sky, sobjs, bkg_redux_global_sky=None,
                             model_noise=True, spat_pix=None,
                             show_profile=False, show_resids=False, show=False):
        """
        Dummy method for local sky-subtraction and extraction.

        Overloaded by class specific skysub and extraction.
        """
        return None, None, None, None, None, None

    def spec_flexure_correct(self, mode="local", sobjs=None):
        """ Correct for spectral flexure

        Spectra are modified in place (wavelengths are shifted)

        Args:
            mode (str):
                "local" - Use sobjs to determine flexure correction
                "global" - Use waveimg and global_sky to determine flexure correction at the centre of the slit
            sobjs (:class:`~pypeit.specobjs.SpecObjs`, None):
                Spectrally extracted objects
        """
        if self.par['flexure']['spec_method'] == 'skip':
            log.info('Skipping flexure correction.')
            return

        # Perform some checks
        if mode == "local" and sobjs is None:
            raise PypeItError("No spectral extractions provided for flexure, using slit center instead")
        elif mode not in ["local", "global"]:
            raise PypeItError("Flexure mode must be 'global' or 'local'.")

        # initialize flex_list
        flex_list = None

        # Prepare a list of slit spectra, if required.
        if mode == "global":
            log.info('Performing global spectral flexure correction')
            gd_slits = np.logical_not(self.extract_bpm)
            trace_spat = 0.5 * (self.slits_left + self.slits_right)
            _global_sky = self.global_sky if self.bkg_redux_global_sky is None else self.bkg_redux_global_sky
            flex_list = flexure.spec_flexure_slit_global(self.sciImg, self.waveimg, _global_sky, self.par,
                                                         self.slits, self.slitmask, trace_spat, gd_slits,
                                                         self.wv_calib, self.pypeline, self.det)
            # Store the slit shifts that were applied to each slit
            # These corrections are later needed so the specobjs metadata contains the total spectral flexure
            self.slitshift = np.zeros(self.slits.nslits)
            for islit in range(self.slits.nslits):
                if gd_slits[islit] and len(flex_list[islit]['shift']) > 0:
                    self.slitshift[islit] = flex_list[islit]['shift'][0]
            # Apply flexure to the new wavelength solution
            log.info("Regenerating wavelength image")
            self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits,
                                                       spat_flexure=self.spat_flexure_shift,
                                                       spec_flexure=self.slitshift)
        elif mode == "local":
            log.info('Performing local spectral flexure correction')
            # Measure flexure:
            # If mode == local: specobjs != None and slitspecs = None
            flex_list = flexure.spec_flexure_slit(self.slits, self.slits.slitord_id, self.extract_bpm,
                                                  self.par['flexure']['spectrum'],
                                                  method=self.par['flexure']['spec_method'],
                                                  mxshft=self.par['flexure']['spec_maxshift'],
                                                  excess_shft=self.par['flexure']['excessive_shift'],
                                                  specobjs=sobjs, slit_specs=None, wv_calib=self.wv_calib,
                                                  minwave=self.par['flexure']['minwave'],
                                                  maxwave=self.par['flexure']['maxwave'])
            # Apply flexure to objects
            for islit in range(self.slits.nslits):
                i_slitord = self.slits.slitord_id[islit]
                indx = sobjs.slitorder_indices(i_slitord)
                this_specobjs = sobjs[indx]
                this_flex_dict = flex_list[islit]
                # Loop through objects
                for ss, sobj in enumerate(this_specobjs):
                    if sobj is None or sobj['BOX_WAVE'] is None:  # Nothing extracted; only the trace exists
                        continue
                    # Interpolate
                    if len(this_flex_dict['shift']) > 0 and this_flex_dict['shift'][ss] is not None:
                        new_sky = sobj.apply_spectral_flexure(this_flex_dict['shift'][ss],
                                                              this_flex_dict['sky_spec'][ss])
                        flex_list[islit]['sky_spec'][ss] = new_sky#.copy()

        # Save QA
        if flex_list is not None:
            basename = f'{self.basename}_{mode}_{self.spectrograph.get_det_name(self.det)}'
            out_dir = os.path.join(self.par['rdx']['redux_path'], 'QA')
            flexure.spec_flexure_qa(self.slits.slitord_id, self.extract_bpm, basename, flex_list,
                                    specobjs=sobjs, out_dir=out_dir)

    def show(self, attr, image=None, showmask=False, sobjs=None,
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

        """

        mask_in = self.sciImg.fullmask if showmask else None

        img_gpm = self.sciImg.select_flag(invert=True)
        detname = self.spectrograph.get_det_name(self.det)

        # TODO Do we still need this here?
        if attr == 'global' and all([a is not None for a in [self.sciImg.image, self.global_sky,
                                                             self.sciImg.fullmask]]):
            # global sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - self.global_sky) * img_gpm.astype(float)
            ch_name = chname if chname is not None else f'global_sky_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, mask=mask_in, clear=clear,
                                            wcs_match=True)
        elif attr == 'local' and all([a is not None for a in [self.sciImg.image, self.skymodel,
                                                              self.sciImg.fullmask]]):
            # local sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - self.skymodel) * img_gpm.astype(float)
            ch_name = chname if chname is not None else f'local_sky_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, mask=mask_in, clear=clear,
                                            wcs_match=True)
        elif attr == 'sky_resid' and all([a is not None for a in [self.sciImg.image, self.skymodel,
                                                                  self.objmodel, self.ivarmodel,
                                                                  self.sciImg.fullmask]]):
            # sky residual map with object included
            image = (self.sciImg.image - self.skymodel) * np.sqrt(self.ivarmodel)
            image *= img_gpm.astype(float)
            ch_name = chname if chname is not None else f'sky_resid_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                            mask=mask_in, clear=clear, wcs_match=True)
        elif attr == 'resid' and all([a is not None for a in [self.sciImg.image, self.skymodel,
                                                              self.objmodel, self.ivarmodel,
                                                              self.sciImg.fullmask]]):
            # full residual map with object model subtractede
            # full model residual map
            image = (self.sciImg.image - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel)
            image *= img_gpm.astype(float)
            ch_name = chname if chname is not None else f'resid_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, cuts=(-5.0, 5.0), mask=mask_in,
                                            clear=clear, wcs_match=True)
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


class MultiSlitExtract(Extract):
    """
    Child of Extract for Multislit and Longslit reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs)


    # TODO: JFH Should we reduce the number of iterations for standards or
    # near-IR redux where the noise model is not being updated?
    def local_skysub_extract(self, global_sky, sobjs, bkg_redux_global_sky=None,
                             spat_pix=None, model_noise=True,
                             show_resids=False, show_profile=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction
        slit by slit.

        Wrapper to :func:`~pypeit.core.skysub.local_skysub_extract`.

        Args:
            global_sky (`numpy.ndarray`_):
                Global sky model
            sobjs (:class:`~pypeit.specobjs.SpecObjs`):
                Class containing the information about the objects found
            bkg_redux_global_sky (`numpy.ndarray`_, optional):
                Sky estimate without background subtraction. This is used for 1d sky spectrum extraction
                in the case bkg_redux=True. Default is None.
            spat_pix (`numpy.ndarray`_, optional):
                Image containing the spatial location of pixels. If not
                input, it will be computed from ``spat_img =
                np.outer(np.ones(nspec), np.arange(nspat))``.
            model_noise (:obj:`bool`, optional):
                If True, construct and iteratively update a model inverse variance image
                using :func:`~pypeit.core.procimg.variance_model`. If False, a
                variance model will not be created and instead the input sciivar will
                always be taken to be the inverse variance. See
                :func:`~pypeit.core.skysub.local_skysub_extract` for more info.
            show_resids (:obj:`bool`, optional):
                Show the model fits and residuals.
            show_profile (:obj:`bool`, optional):
                Show QA for the object profile fitting to the screen. Note
                that this will show interactive matplotlib plots which will
                block the execution of the code until the window is closed.
            show (:obj:`bool`, optional):
                Show debugging plots

        Returns:
            :obj:`tuple`: Return the model sky flux, object flux, inverse
            variance, and mask as `numpy.ndarray`_ objects, and returns a
            :class:`~pypeit.specobjs.SpecObjs`: instance c containing the
            information about the objects found.
        """
        self.global_sky = global_sky

        # get the good slits
        gdslits = np.where(np.logical_not(self.extract_bpm))[0]

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        # NOTE: fullmask is a bit mask, make sure it's treated as such, not a
        # boolean (e.g., bad pixel) mask.
        self.outmask = self.sciImg.fullmask.copy()
        # Initialize to input mask in case no objects were found
        self.extractmask = self.sciImg.select_flag(invert=True)
        # Initialize to zero in case no objects were found
        self.objmodel = np.zeros_like(self.sciImg.image)
        # Set initially to global sky in case no objects were found
        self.skymodel = np.copy(self.global_sky)
        self.bkg_redux_skymodel = bkg_redux_global_sky.copy() if bkg_redux_global_sky is not None else None
        # Set initially to sciivar in case no obects were found.
        self.ivarmodel = np.copy(self.sciImg.ivar)

        # Could actually create a model anyway here, but probably
        # overkill since nothing is extracted
        self.sobjs = sobjs.copy()  # WHY DO WE CREATE A COPY HERE?

        base_gpm = self.sciImg.select_flag(invert=True)

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            log.info("Local sky subtraction and extraction for slit: {:d}".format(slit_spat))
            thisobj = self.sobjs.SLITID == slit_spat    # indices of objects for this slit
            if not np.any(thisobj):
                continue
            # Setup to run local skysub
            thismask = self.slitmask == slit_spat   # pixels for this slit
            # True  = Good, False = Bad for inmask
            ingpm = base_gpm & thismask

            # ... Just for readability
            model_full_slit = self.par['reduce']['extraction']['model_full_slit']
            sigrej = self.par['reduce']['skysub']['sky_sigrej']
            bsp = self.par['reduce']['skysub']['bspline_spacing']
            force_gauss = self.par['reduce']['extraction']['use_user_fwhm']
            sn_gauss = self.par['reduce']['extraction']['sn_gauss']
            use_2dmodel_mask = self.par['reduce']['extraction']['use_2dmodel_mask']
            no_local_sky = self.par['reduce']['skysub']['no_local_sky']
            min_frac_prof = self.par['reduce']['extraction']['min_frac_prof']

            # TODO: skysub.local_skysub_extract() accepts a `prof_nsigma` parameter, but none
            #       is provided here.  Additionally, the ExtractionPar keyword std_prof_nsigma
            #       is not used anywhere in the code.  Should it be be used here, in conjunction
            #       with whether this object IS_STANDARD?
            # prof_nsigma = self.par['reduce']['extraction']['std_prof_nsigma'] if IS_STANDARD else None

            # Local sky subtraction and extraction
            self.skymodel[thismask], _this_bkg_redux_skymodel, self.objmodel[thismask], self.ivarmodel[thismask], self.extractmask[thismask] \
                = skysub.local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                              self.tilts, self.waveimg, self.global_sky,
                                              thismask, self.slits_left[:,slit_idx],
                                              self.slits_right[:, slit_idx],
                                              self.sobjs[thisobj], min_frac_use=min_frac_prof, ingpm=ingpm,
                                              bkg_redux_global_sky=bkg_redux_global_sky,
                                              fwhmimg=self.fwhmimg, flatimg=self.flatimg, spat_pix=spat_pix,
                                              model_full_slit=model_full_slit,
                                              sigrej=sigrej, model_noise=model_noise,
                                              std=self.std_redux, bsp=bsp,
                                              force_gauss=force_gauss, sn_gauss=sn_gauss,
                                              show_profile=show_profile,
                                              # prof_nsigma=prof_nsigma,
                                              use_2dmodel_mask=use_2dmodel_mask,
                                              no_local_sky=no_local_sky,
                                              base_var=self.sciImg.base_var,
                                              count_scale=self.sciImg.img_scale,
                                              adderr=self.sciImg.noise_floor)
            if self.bkg_redux_skymodel is not None:
                self.bkg_redux_skymodel[thismask] = _this_bkg_redux_skymodel

        # Set the bit for pixels which were masked by the extraction.
        # For extractmask, True = Good, False = Bad
        self.outmask.turn_on('EXTRACT', select=base_gpm & np.logical_not(self.extractmask))

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True)
            self.show('resid', sobjs = self.sobjs, slits= True)

        # Return
        return self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs


class EchelleExtract(Extract):
    """
    Child of Extract for Echelle reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs):
        super(EchelleExtract, self).__init__(sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        # TODO Need to test this for 2d coadds
        #self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
        #                    else self.slits.ech_order
        #if self.order_vec is None:
        #    raise PypeItError('Unable to set Echelle orders, likely because they were incorrectly '
        #               'assigned in the relevant SlitTraceSet.')

    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, global_sky, sobjs, bkg_redux_global_sky=None,
                             spat_pix=None, model_noise=True, min_snr=2.0, fit_fwhm=False,
                             show_profile=False, show_resids=False, show_fwhm=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to :func:`~pypeit.core.skysub.local_skysub_extract`.

        Args:
            global_sky (`numpy.ndarray`_):
                Global sky model
            sobjs (:class:`~pypeit.specobjs.SpecObjs`):
                Class containing the information about the objects found
            bkg_redux_global_sky (`numpy.ndarray`_, optional):
                Sky estimate without background subtraction. This is used for 1d sky spectrum extraction
                in the case bkg_redux=True. Default is None.
            spat_pix (`numpy.ndarray`_, optional):
                Image containing the spatial location of pixels. If not
                input, it will be computed from ``spat_img =
                np.outer(np.ones(nspec), np.arange(nspat))``.
            model_noise (:obj:`bool`, optional):
                If True, construct and iteratively update a model inverse variance image
                using :func:`~pypeit.core.procimg.variance_model`. If False, a
                variance model will not be created and instead the input sciivar will
                always be taken to be the inverse variance. See
                `~pypeit.core.skysub.local_skysub_extract` for more info.
            show_resids (:obj:`bool`, optional):
                Show the model fits and residuals.
            show_profile (:obj:`bool`, optional):
                Show QA for the object profile fitting to the screen. Note
                that this will show interactive matplotlib plots which will
                block the execution of the code until the window is closed.
            show (:obj:`bool`, optional):
                Show debugging plots

        Returns:
            :obj:`tuple`: Return the model sky flux, object flux, inverse
            variance, and mask as `numpy.ndarray`_ objects, and returns a
            :class:`~pypeit.specobjs.SpecObjs`: instance c containing the
            information about the objects found.
        """
        self.global_sky = global_sky

        # get the good slits
        gdorders = np.logical_not(self.extract_bpm)

        # Pulled out some parameters to make the method all easier to read
        bsp = self.par['reduce']['skysub']['bspline_spacing']
        sigrej = self.par['reduce']['skysub']['sky_sigrej']
        sn_gauss = self.par['reduce']['extraction']['sn_gauss']
        model_full_slit = self.par['reduce']['extraction']['model_full_slit']
        force_gauss = self.par['reduce']['extraction']['use_user_fwhm']
        use_2dmodel_mask = self.par['reduce']['extraction']['use_2dmodel_mask']
        no_local_sky = self.par['reduce']['skysub']['no_local_sky']
        min_frac_prof = self.par['reduce']['extraction']['min_frac_prof']

        self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
            = skysub.ech_local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                              self.sciImg.fullmask, self.tilts, self.waveimg,
                                              self.global_sky, self.slits_left[:, gdorders],
                                              self.slits_right[:, gdorders], self.slitmask, sobjs,
                                              min_frac_use=min_frac_prof,
                                              bkg_redux_global_sky=bkg_redux_global_sky,
                                              spat_pix=spat_pix,
                                              fwhmimg=self.fwhmimg, flatimg=self.flatimg,
                                              std=self.std_redux, fit_fwhm=fit_fwhm,
                                              min_snr=min_snr, bsp=bsp, sigrej=sigrej,
                                              force_gauss=force_gauss, sn_gauss=sn_gauss,
                                              use_2dmodel_mask=use_2dmodel_mask,
                                              model_full_slit=model_full_slit,
                                              model_noise=model_noise,
                                              show_profile=show_profile,
                                              show_resids=show_resids, show_fwhm=show_fwhm,
                                              no_local_sky=no_local_sky,
                                              base_var=self.sciImg.base_var,
                                              count_scale=self.sciImg.img_scale,
                                              adderr=self.sciImg.noise_floor)
        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs=self.sobjs, slits=True, chname='ech_local')
            self.show('resid', sobjs=self.sobjs, slits=True, chname='ech_resid')

        return self.skymodel, self.bkg_redux_skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs


class SlicerIFUExtract(MultiSlitExtract):
    """
    Child of Extract for IFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs)

    #def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs):
    #    # IFU doesn't extract, and there's no need for a super call here.
    #    return


class FiberExtract(Extract):
    """
    Child of Extract for fiber-fed spectrographs.

    Extracts 1D spectra from each fiber using boxcar and optimal (Horne 1986)
    extraction. Skips local sky subtraction — the global sky model is used
    directly.

    See parent doc string for Args and Attributes.
    """
    def __init__(self, sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs):
        super().__init__(sciImg, slits, sobjs_obj, spectrograph, par, objtype, **kwargs)

    def local_skysub_extract(self, global_sky, sobjs, bkg_redux_global_sky=None,
                             spat_pix=None, model_noise=True,
                             show_resids=False, show_profile=False, show=False):
        """
        Extract fiber spectra without local sky subtraction.

        For each pseudo-slit group of fibers, performs Horne (1986) optimal extraction
        for every fiber SpecObj using flat-derived empirical profiles. The global
        sky model is used directly (no local sky subtraction). Boxcar extraction is
        skipped for fibers whose ``BOX_COUNTS`` attribute is already populated (e.g., by
        :class:`~pypeit.find_objects.FiberFindObjects`), which performs boxcar
        extraction and 1D sky subtraction as part of object finding.

        Parameters
        ----------
        global_sky : `numpy.ndarray`_
            Global sky model, shape ``(nspec, nspat)``.
        sobjs : :class:`~pypeit.specobjs.SpecObjs`
            Objects to extract (multiple per block-slit from FiberFindObjects).
        bkg_redux_global_sky : `numpy.ndarray`_, optional
            Sky estimate without background subtraction.
        spat_pix : `numpy.ndarray`_, optional
            Image containing the spatial location of pixels.
        model_noise : :obj:`bool`, optional
            Not used for fiber extraction.
        show_resids : :obj:`bool`, optional
            Not used for fiber extraction.
        show_profile : :obj:`bool`, optional
            Not used for fiber extraction.
        show : :obj:`bool`, optional
            Show debugging plots.

        Returns
        -------
        :obj:`tuple`
            skymodel, bkg_redux_skymodel, objmodel, ivarmodel, outmask, sobjs
        """
        self.global_sky = global_sky
        nspec, nspat = self.sciImg.image.shape

        # Initialize output arrays
        self.outmask = self.sciImg.fullmask.copy()
        self.extractmask = self.sciImg.select_flag(invert=True)
        self.objmodel = np.zeros_like(self.sciImg.image)
        self.skymodel = np.copy(self.global_sky)
        self.bkg_redux_skymodel = bkg_redux_global_sky.copy() \
            if bkg_redux_global_sky is not None else None
        self.ivarmodel = np.copy(self.sciImg.ivar)
        self.sobjs = sobjs.copy()

        # Build slit ID image
        slitid_img = self.slits.slit_img(pad=0, flexure=self.spat_flexure_shift)

        # Science image minus sky
        inmask = self.sciImg.select_flag(invert=True)
        imgminsky = self.sciImg.image - global_sky
        extract_sky = global_sky if bkg_redux_global_sky is None \
            else bkg_redux_global_sky

        # Extract each fiber.  For the Fiber pypeline the pseudo-slit edges
        # sit on the outermost fibers, not on a physical edge, so the
        # extraction aperture is allowed to spill into unassigned pixels
        # (slitmask == -1, i.e. inter-block gaps and detector edges) to
        # capture the fiber wings.  Pixels belonging to a different slit
        # are always excluded.  For non-fiber pypelines the original
        # in-slit-only mask is preserved.
        #
        # The flat-derived empirical profile for each fiber is built
        # lazily per-iteration and discarded at the end of the iteration:
        # For Binospec IFU, caching all of them at once costs ~134 MB per
        # fiber on a 4k detector (~48 GB for 360 fibers).
        is_fiber = self.spectrograph.pypeline == 'Fiber'
        n_opt = 0
        for sobj in self.sobjs:
            if is_fiber:
                other_slit = (slitid_img != sobj.SLITID) \
                    & (slitid_img != -1)
                thismask = ~other_slit
            else:
                thismask = slitid_img == sobj.SLITID
            sobj_inmask = inmask & thismask

            sobj.extract_boxcar(
                imgminsky, self.sciImg.ivar, sobj_inmask,
                self.waveimg, extract_sky,
                fwhmimg=self.fwhmimg, flatimg=self.flatimg,
                base_var=self.sciImg.base_var,
                count_scale=self.sciImg.img_scale,
                noise_floor=self.sciImg.noise_floor)

            oprof = None
            if self.flatimg is not None:
                oprof = self._build_empirical_profile(
                    self.flatimg, slitid_img, sobj, nspec, nspat)
            if oprof is not None:
                sobj.extract_optimal(
                    imgminsky, self.sciImg.ivar, sobj_inmask,
                    self.waveimg, extract_sky, thismask, oprof,
                    fwhmimg=self.fwhmimg, flatimg=self.flatimg,
                    base_var=self.sciImg.base_var,
                    count_scale=self.sciImg.img_scale,
                    noise_floor=self.sciImg.noise_floor)
                n_opt += 1
            elif not is_fiber:
                # Fallback: per-row Horne extraction with Gaussian profile
                # (not used for Fiber pypeline — empirical profiles required)
                self._optimal_extract_fiber(
                    sobj, slitid_img, inmask, global_sky, None)
            # Explicitly drop the per-fiber profile so peak memory stays
            # constant rather than scaling with fiber count.
            del oprof

        if is_fiber and self.flatimg is not None:
            log.info(f"Empirical-profile optimal extraction succeeded for "
                     f"{n_opt}/{len(self.sobjs)} fibers")

        # Divide out the per-fiber throughput scalar attached by
        # FiberFindObjects.  No wavelength-dependent flat division is
        # applied; spectral flattening is performed downstream by flux
        # calibration.
        if self.spectrograph.pypeline == 'Fiber':
            n_thru = 0
            for sobj in self.sobjs:
                thru = float(getattr(sobj, 'fiber_throughput', 1.0) or 1.0)
                if not np.isfinite(thru) or thru <= 0 or thru == 1.0:
                    continue
                if sobj.BOX_COUNTS is not None:
                    sobj.BOX_COUNTS = sobj.BOX_COUNTS / thru
                    sobj.BOX_COUNTS_IVAR = sobj.BOX_COUNTS_IVAR * thru**2
                    if sobj.BOX_COUNTS_SKY is not None:
                        sobj.BOX_COUNTS_SKY = sobj.BOX_COUNTS_SKY / thru
                if sobj.OPT_COUNTS is not None:
                    sobj.OPT_COUNTS = sobj.OPT_COUNTS / thru
                    sobj.OPT_COUNTS_IVAR = sobj.OPT_COUNTS_IVAR * thru**2
                    if sobj.OPT_COUNTS_SKY is not None:
                        sobj.OPT_COUNTS_SKY = sobj.OPT_COUNTS_SKY / thru
                n_thru += 1
            log.info(f"Applied per-fiber throughput scalar to {n_thru} fibers")

        if self.spectrograph.pypeline != 'Fiber' and \
                hasattr(self.spectrograph, 'apply_throughput_corrections'):
            self.spectrograph.apply_throughput_corrections(self.sobjs, self.det)

        # Re-apply the per-fiber wavelength refinement cached by
        # FiberFindObjects._refine_fiber_wavelengths.  extract_boxcar /
        # extract_optimal both reassign BOX_WAVE / OPT_WAVE from the global
        # waveimg, so without this step the extracted spec1d wavelengths
        # would not match the per-fiber wavelengths that the joint sky
        # model was fit against.  Shift is the same scalar dlam offset for
        # both BOX_WAVE and OPT_WAVE because they share the same fiber
        # trace.
        if self.spectrograph.pypeline == 'Fiber':
            n_shift = 0
            for sobj in self.sobjs:
                shift = float(getattr(sobj, 'wave_refine_shift_AA', 0.0)
                              or 0.0)
                if not np.isfinite(shift) or shift == 0.0:
                    continue
                if sobj.BOX_WAVE is not None:
                    sobj.BOX_WAVE = sobj.BOX_WAVE + shift
                if sobj.OPT_WAVE is not None:
                    sobj.OPT_WAVE = sobj.OPT_WAVE + shift
                n_shift += 1
            if n_shift:
                log.info(f"Applied cached per-fiber wavelength refinement "
                         f"to {n_shift} fibers")

        # Set the bit for pixels masked by extraction
        base_gpm = self.sciImg.select_flag(invert=True)
        self.outmask.turn_on('EXTRACT',
                             select=base_gpm & np.logical_not(self.extractmask))

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs=self.sobjs, slits=True)
            self.show('resid', sobjs=self.sobjs, slits=True)

        return self.skymodel, self.bkg_redux_skymodel, self.objmodel, \
            self.ivarmodel, self.outmask, self.sobjs

    @staticmethod
    def _build_empirical_profile(flatimg, slitmask, sobj, nspec, nspat):
        """
        Build one fiber's empirical spatial profile from the flat field.

        Returns a fresh ``(nspec, nspat)`` profile each call so callers can
        process and discard one fiber at a time -- caching the dense profile
        for every fiber simultaneously would cost ~134 MB per fiber on a
        4k detector (~48 GB for 360 fibers).

        Fiber pseudo-slit edges are not hard apertures (they sit on the
        outermost fiber, not a physical edge), so the per-row aperture
        extends into unassigned pixels (``slitmask == -1`` -- inter-block
        gaps and detector edges) to capture the fiber's PSF wings.  Pixels
        belonging to a different slit are always excluded.

        Per row, if too few pixels in the aperture have a usable flat value
        (NaN / non-positive), the profile falls back to a uniform aperture
        weight (``1 / n_aperture_pixels``).  This guards against BPM rows
        and low-throughput wavelength regions, where a profile built from
        a few positive flat outliers would concentrate the sky on 1-2
        pixels and produce huge spurious spikes.

        Vectorized across rows: assembles the per-row aperture in a
        compact ``(nspec, max_aper)`` form, evaluates the slit and
        flat-validity predicates with numpy, and scatters the normalized
        weights into a fresh dense array.  Bit-identical to the legacy
        row-loop implementation.

        Parameters
        ----------
        flatimg : `numpy.ndarray`_
            Raw (unnormalized) flat field image, shape ``(nspec, nspat)``.
        slitmask : `numpy.ndarray`_
            Slit ID image from ``slits.slit_img(pad=0)``, shape
            ``(nspec, nspat)``.
        sobj : :class:`~pypeit.specobj.SpecObj`
            One fiber with ``SLITID``, ``TRACE_SPAT`` and ``BOX_R_PIX``
            set.
        nspec, nspat : int
            Detector image shape.

        Returns
        -------
        `numpy.ndarray`_ or None
            Per-pixel profile normalized to unit sum per row in the
            aperture, shape ``(nspec, nspat)``.  ``None`` when the fiber's
            slit footprint is empty or no row has any valid data.
        """
        trace = np.asarray(sobj.TRACE_SPAT, dtype=float)
        box_r = float(sobj.BOX_R_PIX)

        # Compact per-row aperture (nspec, max_aper): center on the nearest
        # integer column, then re-filter by the exact |col - trace|
        # predicate so fractional traces match the legacy aperture
        # selection bit-for-bit.
        half = int(np.ceil(box_r))
        offsets = np.arange(-half, half + 1)
        trace_int = np.rint(trace).astype(int)
        cols_2d = trace_int[:, None] + offsets[None, :]
        in_image = (cols_2d >= 0) & (cols_2d < nspat)
        cols_clipped = np.clip(cols_2d, 0, nspat - 1)

        in_aper = (np.abs(cols_2d - trace[:, None]) <= box_r) & in_image
        rows_idx = np.broadcast_to(np.arange(nspec)[:, None],
                                   cols_2d.shape)
        slit_at = slitmask[rows_idx, cols_clipped]
        in_aper &= (slit_at == sobj.SLITID) | (slit_at == -1)

        n_aper = in_aper.sum(axis=1).astype(int)
        flat_vals = flatimg[rows_idx, cols_clipped]
        good_flat = in_aper & np.isfinite(flat_vals) & (flat_vals > 0)
        n_good = good_flat.sum(axis=1).astype(int)

        # Per-row branch: flat-weighted if at least max(3, n_aper // 2)
        # pixels carry a usable flat value; uniform aperture fallback
        # otherwise.  Matches the legacy `good.sum() < max(3, ...)` check.
        flat_min = np.maximum(3, n_aper // 2)
        use_flat = (n_aper > 0) & (n_good >= flat_min)
        use_uniform = (n_aper > 0) & ~use_flat

        weights = np.zeros_like(flat_vals)

        if np.any(use_flat):
            flat_for_norm = np.where(good_flat, flat_vals, 0.0)
            row_sum = flat_for_norm.sum(axis=1)
            valid = use_flat & (row_sum > 0)
            if np.any(valid):
                weights[valid] = (flat_for_norm[valid]
                                  / row_sum[valid][:, None])
            # Rows where the flat-weighted aperture sums to zero (e.g. all
            # good vals == 0) are skipped, matching the legacy
            # `if row_sum > 0` guard.
            use_flat = valid

        if np.any(use_uniform):
            uniform = 1.0 / np.maximum(n_aper, 1)[:, None]
            weights[use_uniform] = np.where(
                in_aper[use_uniform], uniform[use_uniform], 0.0)

        if not np.any(use_flat | use_uniform):
            return None

        # Scatter compact (nspec, max_aper) weights into a fresh dense
        # (nspec, nspat) image.  Only in-image entries are written; cols
        # clipped to nspat-1 from off-image entries would otherwise
        # overwrite a real aperture pixel with zero.
        prof = np.zeros((nspec, nspat), dtype=float)
        mask_flat = in_image.ravel()
        rs = rows_idx.ravel()[mask_flat]
        cs = cols_2d.ravel()[mask_flat]
        ws = weights.ravel()[mask_flat]
        prof[rs, cs] = ws
        return prof

    def _optimal_extract_fiber(self, sobj, slitid_img, inmask, skymodel,
                               empirical_profiles):
        """
        Perform per-row Horne (1986) optimal extraction for a single fiber.

        Fallback method when flat-derived empirical profiles are not available.
        Uses a Gaussian profile centered on the fiber trace with sigma derived
        from the slit width.

        Sets OPT_WAVE, OPT_COUNTS, OPT_COUNTS_IVAR, OPT_MASK,
        OPT_COUNTS_SKY, OPT_COUNTS_SIG on the SpecObj.

        Parameters
        ----------
        sobj : :class:`~pypeit.specobj.SpecObj`
            The SpecObj for this fiber. Modified in place.
        slitid_img : `numpy.ndarray`_
            Slit ID image.
        inmask : `numpy.ndarray`_
            Good pixel mask (True = good).
        skymodel : `numpy.ndarray`_
            Sky model image.
        empirical_profiles : :obj:`dict` or None
            Not used in fallback mode (kept for interface compatibility).
        """
        nspec, nspat = self.sciImg.image.shape
        onslit = slitid_img == sobj.SLITID
        good = onslit & inmask

        # Use the fiber trace and BOX_R_PIX for Gaussian profile
        trace_center = sobj.TRACE_SPAT
        # Sigma ~ half the fiber aperture / 2.3548 (FWHM -> sigma)
        trace_sigma = np.full(nspec, sobj.BOX_R_PIX / 1.1774)  # half-FWHM -> sigma

        opt_flux = np.zeros(nspec)
        opt_ivar = np.zeros(nspec)
        opt_wave = np.zeros(nspec)
        opt_sky = np.zeros(nspec)
        opt_mask = np.zeros(nspec, dtype=bool)

        for row in range(nspec):
            pix = good[row, :]
            if not np.any(pix):
                continue

            cols = np.where(pix)[0]
            # Restrict to pixels within BOX_R_PIX of the fiber trace
            in_aperture = np.abs(cols - trace_center[row]) <= sobj.BOX_R_PIX
            cols = cols[in_aperture]
            if len(cols) == 0:
                continue

            opt_wave[row] = np.median(self.waveimg[row, cols])

            sig = trace_sigma[row]
            if sig < 0.1:
                sig = 1.0
            cen = trace_center[row]
            profile = np.exp(-0.5 * ((cols - cen) / sig) ** 2)

            psum = np.sum(profile)
            if psum <= 0:
                continue
            profile = profile / psum

            iv = self.sciImg.ivar[row, cols]
            good_iv = iv > 0
            if not np.any(good_iv):
                continue

            # Horne Eq. 8: flux = sum(P * ivar * data) / sum(P^2 * ivar)
            denom = np.sum(profile[good_iv] ** 2 * iv[good_iv])
            if denom <= 0:
                continue

            imgminsky_row = self.sciImg.image[row, cols] - skymodel[row, cols]
            opt_flux[row] = np.sum(
                profile[good_iv] * iv[good_iv]
                * imgminsky_row[good_iv]) / denom
            opt_sky[row] = np.sum(
                profile[good_iv] * iv[good_iv]
                * skymodel[row, cols[good_iv]]) / denom
            # Horne Eq. 9: var = sum(P) / sum(P^2 * ivar)
            opt_ivar[row] = denom / np.sum(profile[good_iv])
            opt_mask[row] = True

        sobj.OPT_WAVE = opt_wave
        sobj.OPT_COUNTS = opt_flux
        sobj.OPT_COUNTS_IVAR = opt_ivar
        sobj.OPT_COUNTS_SIG = np.sqrt(utils.inverse(opt_ivar))
        sobj.OPT_MASK = opt_mask
        sobj.OPT_COUNTS_SKY = opt_sky
