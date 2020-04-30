"""
Main driver class for skysubtraction and extraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""

import inspect
import numpy as np
import os

from astropy import stats
from abc import ABCMeta

from linetools import utils as ltu

from pypeit import specobjs
from pypeit import ginga, msgs
from pypeit.core import skysub, extract, wave, flexure
from pypeit import wavecalib

from IPython import embed

class Reduce(object):
    """
    This class will organize and run actions related to
    finding objects, sky subtraction, and extraction for
    a Science or Standard star exposure

    Args:
        sciImg (pypeit.images.scienceimage.ScienceImage):
        spectrograph (pypeit.spectrograph.Spectrograph):
        par (:class:`pypeit.par.pyepeitpar.PypeItPar`):
        caliBrate (:class:`pypeit.calibrations.Calibrations`):
        objtype (str):
           Specifies object being reduced 'science' 'standard' 'science_coadd2d'
        det (int, optional):
           Detector indice
        setup (str, optional):
           Used for naming
        maskslits (ndarray, optional):
          Specifies masked out slits
          True = Masked
        show (bool, optional):
           Show plots along the way?

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
        sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
            Only object finding but no extraction
        sobjs (SpecObsj):
            Final extracted object list with trace corrections applied
        spat_flexure_shift (float):
        tilts (`numpy.ndarray`_):
            WaveTilts images generated on-the-spot
        waveimg (`numpy.ndarray`_):
            WaveImage image generated on-the-spot

    """

    __metaclass__ = ABCMeta



    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, sciImg, spectrograph, par, caliBrate,
                 objtype, ir_redux=False, det=1, std_redux=False, show=False,
                 binning=None, setup=None):
        """
        Instantiate the Reduce subclass appropriate for the provided
        spectrograph.

        The class must be subclassed from Reduce.  See :class:`Reduce` for
        the description of the valid keyword arguments.

        Args:
            sciImg (pypeit.images.scienceimage.ScienceImage):
            spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            par (pypeit.par.pyepeitpar.PypeItPar):
            caliBrate (:class:`pypeit.calibrations.Calibrations`):

            **kwargs
                Passed to Parent init

        Returns:
            :class:`pypeit.reduce.Reduce`:
        """
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (spectrograph.pypeline + 'Reduce'))(
            sciImg, spectrograph, par, caliBrate, objtype, ir_redux=ir_redux, det=det,
            std_redux=std_redux, show=show,binning=binning, setup=setup)


    def __init__(self, sciImg, spectrograph, par, caliBrate,
                 objtype, ir_redux=False, det=1, std_redux=False, show=False,
                 binning=None, setup=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        # Parse
        # Slit pieces
        #   WARNING -- It is best to unpack here then pass around self.slits
        #      Otherwise you have to keep in mind flexure, tweaking, etc.

        # Flexure
        self.spat_flexure_shift = None
        if objtype == 'science':
            if self.par['scienceframe']['process']['spat_flexure_correct']:
                self.spat_flexure_shift = self.sciImg.spat_flexure
        elif objtype == 'standard':
            if self.par['calibrations']['standardframe']['process']['spat_flexure_correct']:
                self.spat_flexure_shift = self.sciImg.spat_flexure
        elif objtype == 'science_coadd2d':
            self.spat_flexure_shift = None
        else:
            msgs.error("Not ready for this objtype in Reduce")

        # Slits
        self.slits = caliBrate.slits
        # Select the edges to use
        self.slits_left, self.slits_right, _ \
                = self.slits.select_edges(flexure=self.spat_flexure_shift)

        # Slitmask
        self.slitmask = self.slits.slit_img(flexure=self.spat_flexure_shift,
                                           exclude_flag=self.slits.bitmask.exclude_for_reducing)
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)
        # For echelle
        self.spatial_coo = self.slits.spatial_coordinates(flexure=self.spat_flexure_shift)

        # Internal bpm mask
        self.reduce_bpm = (self.slits.mask > 0) & (np.invert(self.slits.bitmask.flagged(
                        self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing)))
        self.reduce_bpm_init = self.reduce_bpm.copy()

        # These may be None (i.e. COADD2D)
        self.waveTilts = caliBrate.wavetilts
        self.wv_calib = caliBrate.wv_calib


        # Load up other input items
        self.ir_redux = ir_redux
        self.std_redux = std_redux
        self.det = det
        self.binning = binning
        self.setup = setup
        self.pypeline = spectrograph.pypeline
        self.reduce_show = show

        self.steps = []

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.initial_sky = None
        self.global_sky = None
        self.skymask = None
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None  # Only object finding but no extraction
        self.sobjs = None  # Final extracted object list with trace corrections applied



    def parse_manual_dict(self, manual_dict, neg=False):
        """
        Parse the manual dict
        This method is here mainly to deal with negative images

        Args:
            manual_dict (dict or None):
            neg (bool, optional):

        Returns:
            None or dict:  None if no matches; dict if there are for manual extraction

        """
        if manual_dict is None:
            return None
        #
        dets = manual_dict['hand_extract_det']
        # Grab the ones we want
        gd_det = dets > 0
        if neg:
            gd_det = np.invert(gd_det)
        # Any?
        if not np.any(gd_det):
            return None
        # Fill
        manual_extract_dict = {}
        for key in manual_dict.keys():
            sgn = 1
            if key == 'hand_extract_det':
                sgn = -1
            manual_extract_dict[key] = sgn*manual_dict[key][gd_det]
        # Return
        return manual_extract_dict

    def extract(self, global_sky, sobjs_obj):
        """
        Main method to extract spectra from the ScienceImage

        Args:
            global_sky (ndarray):
                Sky estimate
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                List of SpecObj that have been found and traced
        """
        # This holds the objects, pre-extraction
        self.sobjs_obj = sobjs_obj

        if self.par['reduce']['extraction']['skip_optimal']:  # Boxcar only with global sky subtraction
            msgs.info("Skipping optimal extraction")

            # This will hold the extracted objects
            self.sobjs = self.sobjs_obj.copy()
            # Only extract positive objects
            self.sobjs.purge_neg()

            # Quick loop over the objects
            for iobj in range(self.sobjs.nobj):
                sobj = self.sobjs[iobj]
                plate_scale = self.get_platescale(sobj)
                # True  = Good, False = Bad for inmask
                thismask = self.slitmask == sobj.SLITID  # pixels for this slit
                inmask = (self.sciImg.fullmask == 0) & thismask
                # Do it
                extract.extract_boxcar(self.sciImg.image, self.sciImg.ivar,
                                               inmask, self.waveimg,
                                               global_sky, self.sciImg.rn2img,
                                               self.par['reduce']['extraction']['boxcar_radius']/plate_scale,
                                               sobj)
            # Fill up extra bits and pieces
            self.objmodel = np.zeros_like(self.sciImg.image)
            self.ivarmodel = np.copy(self.sciImg.ivar)
            self.outmask = self.sciImg.fullmask
            self.skymodel = global_sky.copy()
        else:  # Local sky subtraction and optimal extraction.
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
                self.local_skysub_extract(global_sky, self.sobjs_obj,
                                          model_noise=(not self.ir_redux),
                                          show_profile=self.reduce_show,
                                          show=self.reduce_show)
        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    def get_platescale(self, sobj):
        """
        Return the platescale for the current detector/echelle order

        Over-loaded by the children

        Args:
            sobj (:class:`pypeit.specobj.SpecObj`):

        Returns:
            float:

        """
        pass

    def run(self, basename=None, ra=None, dec=None, obstime=None,
            std_trace=None, manual_extract_dict=None, show_peaks=False):
        """
        Primary code flow for PypeIt reductions

        *NOT* used by COADD2D

        Args:
            basename (str, optional):
                Required if flexure correction is to be applied
            ra (str, optional):
                Required if helio-centric correction is to be applied
            dec (str, optional):
                Required if helio-centric correction is to be applied
            obstime (:obj:`astropy.time.Time`, optional):
                Required if helio-centric correction is to be applied
            std_trace (np.ndarray, optional):
                Trace of the standard star
            manual_extract_dict (dict, optional):
            show_peaks (bool, optional):
                Show peaks in find_objects methods

        Returns:
            tuple: skymodel (ndarray), objmodel (ndarray), ivarmodel (ndarray),
               outmask (ndarray), sobjs (SpecObjs), waveimg (`numpy.ndarray`_),
               tilts (`numpy.ndarray`_).
               See main doc string for description

        """

        # Deal with dynamic calibrations
        # Tilts
        self.waveTilts.is_synced(self.slits)
        #   Deal with Flexure
        if self.par['calibrations']['tiltframe']['process']['spat_flexure_correct']:
            _spat_flexure = 0. if self.spat_flexure_shift is None else self.spat_flexure_shift
            # If they both shifted the same, there will be no reason to shift the tilts
            tilt_flexure_shift = _spat_flexure - self.waveTilts.spat_flexure
        else:
            tilt_flexure_shift = self.spat_flexure_shift
        self.tilts = self.waveTilts.fit2tiltimg(self.slitmask, flexure=tilt_flexure_shift)

        # Wavelengths (on unmasked slits)
        self.waveimg = wavecalib.build_waveimg(self.spectrograph, self.tilts, self.slits,
                                               self.wv_calib, spat_flexure=self.spat_flexure_shift)

        # First pass object finding
        self.sobjs_obj, self.nobj, skymask_init = \
            self.find_objects(self.sciImg.image, std_trace=std_trace,
                              show_peaks=show_peaks,
                              show=self.reduce_show & (not self.std_redux),
                              manual_extract_dict=manual_extract_dict)

        # Global sky subtract
        self.initial_sky = \
            self.global_skysub(skymask=skymask_init).copy()

        # Second pass object finding on sky-subtracted image
        if (not self.std_redux) and (not self.par['reduce']['findobj']['skip_second_find']):
            self.sobjs_obj, self.nobj, self.skymask = \
                self.find_objects(self.sciImg.image - self.initial_sky,
                                  std_trace=std_trace,
                                  show=self.reduce_show,
                                  show_peaks=show_peaks,
                                  manual_extract_dict=manual_extract_dict)
        else:
            msgs.info("Skipping 2nd run of finding objects")

        # Do we have any positive objects to proceed with?
        if self.nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            if (self.std_redux or self.par['reduce']['extraction']['skip_optimal'] or
                    self.par['reduce']['findobj']['skip_second_find']):
                self.global_sky = self.initial_sky.copy()
            else:
                self.global_sky = self.global_skysub(skymask=self.skymask,
                                                     show=self.reduce_show)
            # Extract + Return
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = self.extract(self.global_sky, self.sobjs_obj)
        else:  # No objects, pass back what we have
            self.skymodel = self.initial_sky
            self.objmodel = np.zeros_like(self.sciImg.image)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # Set to the initial mask in case no objects were found
            self.outmask = self.sciImg.fullmask
            # empty specobjs object from object finding
            self.sobjs = self.sobjs_obj

        # Purge out the negative objects if this was a near-IR reduction.
        if self.ir_redux:
            self.sobjs.purge_neg()

        # Finish up
        if self.sobjs.nobj == 0:
            msgs.warn('No objects to extract!')
        else:
            # TODO -- Should we move these to redux.run()?
            # Flexure correction if this is not a standard star
            if not self.std_redux:
                self.spec_flexure_correct(self.sobjs, basename)
            # Heliocentric
            radec = ltu.radec_to_coord((ra, dec))
            self.helio_correct(self.sobjs, radec, obstime)

        # Update the mask
        reduce_masked = np.where(np.invert(self.reduce_bpm_init) & self.reduce_bpm)[0]
        if len(reduce_masked) > 0:
            self.slits.mask[reduce_masked] = self.slits.bitmask.turn_on(
                self.slits.mask[reduce_masked], 'BADREDUCE')

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs, \
               self.waveimg, self.tilts

    def find_objects(self, image, std_trace=None,
                     show_peaks=False, show_fits=False,
                     show_trace=False, show=False, manual_extract_dict=None,
                     debug=False):
        """
        Single pass at finding objects in the input image

        If self.ir_redux is True, do a search for negative objects too

        Args:
            image (np.ndarray):
                Input image
            std_trace (ndarray, optional):
            show_peaks (bool, optional):
            show_fits (bool, optional):
            show_trace (bool, optional):
            show (bool, optional):
            manual_extract_dict (dict, optional):
            debug (bool, optional):

        Returns:
            specobjs (:class:`pypeit.specobjs.SpecObjs`), int, np.ndarray:
               Objects found,  number of objects found, skymask

        """

        # Positive image
        parse_manual = self.parse_manual_dict(manual_extract_dict, neg=False)
        sobjs_obj_single, nobj_single, skymask_pos = \
            self.find_objects_pypeline(image,
                                       std_trace=std_trace,
                                       show_peaks=show_peaks, show_fits=show_fits,
                                       show_trace=show_trace,
                                       manual_extract_dict=parse_manual, debug=debug)

        # For nobj we take only the positive objects
        if self.ir_redux:
            msgs.info("Finding objects in the negative image")
            # Parses
            parse_manual = self.parse_manual_dict(manual_extract_dict, neg=True)
            sobjs_obj_single_neg, nobj_single_neg, skymask_neg = \
                self.find_objects_pypeline(-image, std_trace=std_trace,
                                           show_peaks=show_peaks, show_fits=show_fits,
                                           show_trace=show_trace,
                                           manual_extract_dict=parse_manual,
                                           debug=debug)
            # Mask
            skymask = skymask_pos & skymask_neg
            # Add (if there are any)
            sobjs_obj_single.append_neg(sobjs_obj_single_neg)
        else:
            skymask = skymask_pos

        if show:
            self.show('image', image=image*(self.sciImg.fullmask == 0), chname='objfind',sobjs=sobjs_obj_single, slits=True)

        # For nobj we take only the positive objects
        return sobjs_obj_single, nobj_single, skymask

    def find_objects_pypeline(self, image, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, debug=False,
                              manual_extract_dict=None):

        """
         Dummy method for object finding. Overloaded by class specific object finding.

         Returns:

         """
        return None, None, None

    def global_skysub(self, skymask=None, update_crmask=True,
                      show_fit=False, show=False, show_objs=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Args:
            skymask (np.ndarray):
            update_crmask (bool, optional):
            show_fit (bool, optional):
            show (bool, optional):
            show_objs (bool, optional):

        Returns:
            numpy.ndarray: image of the the global sky model

        """
        # Prep
        self.global_sky = np.zeros_like(self.sciImg.image)
        # Parameters for a standard star
        if self.std_redux:
            sigrej = 7.0
            update_crmask = False
            if not self.par['reduce']['skysub']['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return self.global_sky
        else:
            sigrej = 3.0

        gdslits = np.where(np.invert(self.reduce_bpm))[0]

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            msgs.info("Global sky subtraction for slit: {:d}".format(slit_spat))
            thismask = self.slitmask == slit_spat
            inmask = (self.sciImg.fullmask == 0) & thismask & skymask_now
            # Find sky
            self.global_sky[thismask] \
                    = skysub.global_skysub(self.sciImg.image, self.sciImg.ivar, self.tilts,
                                           thismask, self.slits_left[:,slit_idx],
                                           self.slits_right[:,slit_idx],
                                           inmask=inmask, sigrej=sigrej,
                                           bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                           no_poly=self.par['reduce']['skysub']['no_poly'],
                                           pos_mask=(not self.ir_redux), show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.reduce_bpm[slit_idx] = True

        if update_crmask:
            # Find CRs with sky subtraction
            self.sciImg.build_crmask(self.par['scienceframe']['process'],
                                   subtract_img=self.global_sky)
            # Update the fullmask
            self.sciImg.update_mask_cr(self.sciImg.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', slits=True, sobjs=sobjs_show, clear=False)

        # Return
        return self.global_sky

    def local_skysub_extract(self, global_sky, sobjs,
                             model_noise=True, spat_pix=None,
                             show_profile=False, show_resids=False, show=False):
        """
        Dummy method for locak skysubtraction and extraction.

        Overloaded by class specific skysub and extraction.

        Returns:

        """

        return None, None, None, None, None

    def spec_flexure_correct(self, sobjs, basename):
        """ Correct for spectral flexure

        Spectra are modified in place (wavelengths are shifted)

        Wrapper to flexure.flexure_obj()

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
            basename (str):

        """

        if self.par['flexure']['spec_method'] != 'skip':
            # Measure
            flex_list = flexure.spec_flexure_obj(sobjs, self.slits.slitord_id, self.reduce_bpm,
                                                 self.par['flexure']['spec_method'],
                                                 self.par['flexure']['spectrum'],
                                                 mxshft=self.par['flexure']['spec_maxshift'])
            # QA
            flexure.spec_flexure_qa(sobjs, self.slits.slitord_id, self.reduce_bpm, basename, self.det, flex_list,
                                    out_dir=os.path.join(self.par['rdx']['redux_path'], 'QA'))
        else:
            msgs.info('Skipping flexure correction.')

    def helio_correct(self, sobjs, radec, obstime):
        """ Perform a heliocentric correction

        Wrapper to wave.geomotion_correct()

        Input objects are modified in place

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
            radec (astropy.coordiantes.SkyCoord):
            obstime (:obj:`astropy.time.Time`):

        """
        # Helio, correct Earth's motion
        if (self.par['calibrations']['wavelengths']['frame'] in ['heliocentric', 'barycentric']) \
                and (self.par['calibrations']['wavelengths']['reference'] != 'pixel'):
            # TODO change this keyword to refframe instead of frame
            msgs.info("Performing a {0} correction".format(self.par['calibrations']['wavelengths']['frame']))
            # Good slitord
            gd_slitord = self.slits.slitord_id[np.invert(self.reduce_bpm)]
            vel, vel_corr = wave.geomotion_correct(sobjs, radec, obstime, gd_slitord,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   self.par['calibrations']['wavelengths']['frame'])
        else:
            msgs.info('A wavelength reference-frame correction will not be performed.')
            vel_corr = None

        return

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

        Returns
        -------

        """

        if showmask:
            mask_in = self.sciImg.fullmask
            bitmask_in = self.sciImg.bitmask
        else:
            mask_in = None
            bitmask_in = None

        if attr == 'global':
            # global sky subtraction
            if self.sciImg.image is not None and self.global_sky is not None and self.sciImg.fullmask is not None:
                # sky subtracted image
                image = (self.sciImg.image - self.global_sky)*(self.sciImg.fullmask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.sciImg.fullmask == 0], sigma_lower=5.0,
                                                       sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                ch_name = chname if chname is not None else 'global_sky_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                              mask=mask_in, clear=clear, wcs_match=True)
                                              #, cuts=(cut_min, cut_max))
        elif attr == 'local':
            # local sky subtraction
            if self.sciImg.image is not None and self.skymodel is not None and self.sciImg.fullmask is not None:
                # sky subtracted image
                image = (self.sciImg.image - self.skymodel)*(self.sciImg.fullmask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.sciImg.fullmask == 0], sigma_lower=5.0,
                                                       sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                ch_name = chname if chname is not None else 'local_sky_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                              mask=mask_in, clear=clear, wcs_match=True)
                                              #, cuts=(cut_min, cut_max))
        elif attr == 'sky_resid':
            # sky residual map with object included
            if self.sciImg.image is not None and self.skymodel is not None \
                    and self.objmodel is not None and self.ivarmodel is not None \
                    and self.sciImg.fullmask is not None:
                image = (self.sciImg.image - self.skymodel) * np.sqrt(self.ivarmodel)
                image *= (self.sciImg.fullmask == 0)
                ch_name = chname if chname is not None else 'sky_resid_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                              bitmask=bitmask_in, mask=mask_in, clear=clear,
                                              wcs_match=True)
        elif attr == 'resid':
            # full residual map with object model subtractede
            if self.sciImg.image is not None and self.skymodel is not None \
                    and self.objmodel is not None and self.ivarmodel is not None \
                    and self.sciImg.fullmask is not None:
                # full model residual map
                image = (self.sciImg.image - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel)
                image *= (self.sciImg.fullmask == 0)
                ch_name = chname if chname is not None else 'resid_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                              bitmask=bitmask_in, mask=mask_in, clear=clear,
                                              wcs_match=True)
        elif attr == 'image':
            ch_name = chname if chname is not None else 'image'
            viewer, ch = ginga.show_image(image, chname=ch_name, clear=clear, wcs_match=True)
        else:
            msgs.warn("Not an option for show")

        if sobjs is not None:
            for spec in sobjs:
                color = 'magenta' if spec.hand_extract_flag else 'orange'
                ginga.show_trace(viewer, ch, spec.TRACE_SPAT, spec.NAME, color=color)

        if slits and self.slits_left is not None:
            ginga.show_slits(viewer, ch, self.slits_left, self.slits_right)

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


class MultiSlitReduce(Reduce):
    """
    Child of Reduce for Multislit and Longslit reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super(MultiSlitReduce, self).__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)

    def get_platescale(self, dummy):
        """
        Return the platescale for multislit.
        The input argument is ignored

        Args:
            dummy:
                ignored
                Keeps argument lists the same amongst the children

        Returns:
            float:

        """
        plate_scale = self.sciImg.detector.platescale
        return plate_scale

    def find_objects_pypeline(self, image, std_trace=None,
                              manual_extract_dict=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, debug=False):
        """
        Pipeline specific find objects routine

        Args:
            image (np.ndarray):
            std_trace (np.ndarray, optional):
            manual_extract_dict (dict, optional):
            show_peaks (bool, optional):
              Generate QA showing peaks identified by object finding
            show_fits (bool, optional):
              Generate QA  showing fits to traces
            show_trace (bool, optional):
              Generate QA  showing traces identified. Requires an open ginga RC modules window
            show (bool, optional):
            debug (bool, optional):

        Returns:
            tuple:
                specobjs : Specobjs object
                    Container holding Specobj objects
                nobj (int):
                    Number of objects identified
                skymask : ndarray
                    Boolean image indicating which pixels are useful for global sky subtraction

        """
        gdslits = np.where(np.invert(self.reduce_bpm))[0]

        # create the ouptut image for skymask
        skymask = np.zeros_like(image, dtype=bool)
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            qa_title ="Finding objects on slit # {:d}".format(slit_spat)
            msgs.info(qa_title)
            thismask = self.slitmask == slit_spat
            inmask = (self.sciImg.fullmask == 0) & thismask
            # Find objects
            specobj_dict = {'SLITID': slit_spat,
                            'DET': self.det, 'OBJTYPE': self.objtype,
                            'PYPELINE': self.pypeline}

            # TODO we need to add QA paths and QA hooks. QA should be
            # done through objfind where all the relevant information
            # is. This will be a png file(s) per slit.

            sobjs_slit, skymask[thismask] = \
                    extract.objfind(image, thismask,
                                self.slits_left[:,slit_idx],
                                self.slits_right[:,slit_idx],
                                inmask=inmask, ir_redux=self.ir_redux,
                                ncoeff=self.par['reduce']['findobj']['trace_npoly'],
                                std_trace=std_trace,
                                sig_thresh=self.par['reduce']['findobj']['sig_thresh'],
                                hand_extract_dict=manual_extract_dict,
                                specobj_dict=specobj_dict, show_peaks=show_peaks,
                                show_fits=show_fits, show_trace=show_trace,
                                trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
                                cont_fit=self.par['reduce']['findobj']['find_cont_fit'],
                                npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
                                fwhm=self.par['reduce']['findobj']['find_fwhm'],
                                maxdev=self.par['reduce']['findobj']['find_maxdev'],
                                qa_title=qa_title, nperslit=self.par['reduce']['findobj']['maxnumber'],
                                debug_all=debug)
            sobjs.add_sobj(sobjs_slit)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.sciImg.fullmask == 0), chname = 'objfind',
                      sobjs=sobjs, slits=True)

        # Return
        return sobjs, len(sobjs), skymask


    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, global_sky, sobjs, spat_pix=None, model_noise=True, show_resids=False,
                             show_profile=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Args:
            global_sky (np.ndarray):
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
            spat_pix (np.ndarray, optional):
            model_noise (bool, optional):
            show_resids (bool, optional):
            show_profile (bool, optional):
            show (bool, optional):

        Returns:
            tuple: skymodel (np.ndarray), objmodel (np.ndarray), ivarmodel (np.ndarray), outmask (np.ndarray), sobjs

        """
        self.global_sky = global_sky

        # get the good slits
        gdslits = np.where(np.invert(self.reduce_bpm))[0]

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        self.outmask = np.copy(self.sciImg.fullmask)
        # Initialize to input mask in case no objects were found
        self.extractmask = (self.sciImg.fullmask == 0)
        # Initialize to zero in case no objects were found
        self.objmodel = np.zeros_like(self.sciImg.image)
        # Set initially to global sky in case no objects were found
        self.skymodel  = np.copy(self.global_sky)
        # Set initially to sciivar in case no obects were found.
        self.ivarmodel = np.copy(self.sciImg.ivar)

        # Could actually create a model anyway here, but probably
        # overkill since nothing is extracted
        self.sobjs = sobjs.copy()  # WHY DO WE CREATE A COPY HERE?
        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            msgs.info("Local sky subtraction and extraction for slit: {:d}".format(slit_spat))
            thisobj = self.sobjs.SLITID == slit_spat    # indices of objects for this slit
            if np.any(thisobj):
                thismask = self.slitmask == slit_spat   # pixels for this slit
                # True  = Good, False = Bad for inmask
                ingpm = (self.sciImg.fullmask == 0) & thismask
                # Local sky subtraction and extraction
                self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], \
                    self.extractmask[thismask] = skysub.local_skysub_extract(
                    self.sciImg.image, self.sciImg.ivar, self.tilts, self.waveimg,
                    self.global_sky, self.sciImg.rn2img,
                    thismask, self.slits_left[:,slit_idx], self.slits_right[:, slit_idx],
                    self.sobjs[thisobj], ingpm,
                    spat_pix=spat_pix,
                    model_full_slit=self.par['reduce']['extraction']['model_full_slit'],
                    box_rad=self.par['reduce']['extraction']['boxcar_radius']/self.get_platescale(None),
                    sigrej=self.par['reduce']['skysub']['sky_sigrej'],
                    model_noise=model_noise, std=self.std_redux,
                    bsp=self.par['reduce']['skysub']['bspline_spacing'],
                    sn_gauss=self.par['reduce']['extraction']['sn_gauss'],
                    show_profile=show_profile,
                    use_2dmodel_mask=self.par['reduce']['extraction']['use_2dmodel_mask'])

        # Set the bit for pixels which were masked by the extraction.
        # For extractmask, True = Good, False = Bad
        iextract = (self.sciImg.fullmask == 0) & (self.extractmask == False)
        self.outmask[iextract] = self.sciImg.bitmask.turn_on(self.outmask[iextract], 'EXTRACT')

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True)
            self.show('resid', sobjs = self.sobjs, slits= True)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs


class EchelleReduce(Reduce):
    """
    Child of Reduce for Echelle reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super(EchelleReduce, self).__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
                            else self.spectrograph.order_vec(self.spatial_coo)

    def get_platescale(self, sobj):
        """
        Return the plate scale for the given current echelle order
        based on the order index

        Args:
            sobj (:class:`pypeit.specobj.SpecObj`):

        Returns:
            float:

        """
        return self.spectrograph.order_platescale(sobj.ECH_ORDER, binning=self.binning)[0]

    def get_positive_sobj(self, specobjs, iord):
        """
        Return the current object from self.sobjs_obj

        Args:
            iord (int):
                Echelle order index

        Returns:
            :class:`pypeit.specobj.SpecObj`:

        """
        # pos indices of objects for this slit
        thisobj = (self.sobjs_obj.ech_orderindx == iord) & (self.sobjs_obj.ech_objid > 0)
        return self.sobjs_obj[np.where(thisobj)[0][0]]

    def find_objects_pypeline(self, image, std_trace=None,
                              show=False, show_peaks=False, show_fits=False,
                              show_trace=False, debug=False,
                              manual_extract_dict=None):
        """
         Pipeline specific find objects routine

         Args:
             image (np.ndarray):
             std_trace (np.ndarray, optional):
             manual_extract_dict (dict, optional):
             show_peaks (bool, optional):
               Generate QA showing peaks identified by object finding
             show_fits (bool, optional):
               Generate QA  showing fits to traces
             show_trace (bool, optional):
               Generate QA  showing traces identified. Requires an open ginga RC modules window
             show (bool, optional):
             debug (bool, optional):

         Returns:
             tuple:
                 specobjs : Specobjs object
                     Container holding Specobj objects
                 nobj (int):
                     Number of objects identified
                 skymask : ndarray
                     Boolean image indicating which pixels are useful for global sky subtraction

        """
        # create the ouptut image for skymask
        skymask = np.zeros_like(image, dtype=bool)

        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        inmask = self.sciImg.fullmask == 0
        # Find objects
        specobj_dict = {'SLITID': 999, #'orderindx': 999,
                        'DET': self.det, 'OBJTYPE': self.objtype, 'PYPELINE': self.pypeline}

        sobjs_ech, skymask[self.slitmask > -1] = extract.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, self.slits_left, self.slits_right,
            self.order_vec, self.reduce_bpm,
            spec_min_max=np.vstack((self.slits.specmin, self.slits.specmax)),
            inmask=inmask, ir_redux=self.ir_redux, ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            hand_extract_dict=manual_extract_dict, plate_scale=plate_scale,
            std_trace=std_trace,
            specobj_dict=specobj_dict,
            sig_thresh=self.par['reduce']['findobj']['sig_thresh'],
            show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.par['reduce']['findobj']['find_trim_edge'],
            cont_fit=self.par['reduce']['findobj']['find_cont_fit'],
            npoly_cont=self.par['reduce']['findobj']['find_npoly_cont'],
            fwhm=self.par['reduce']['findobj']['find_fwhm'],
            maxdev=self.par['reduce']['findobj']['find_maxdev'],
            max_snr=self.par['reduce']['findobj']['ech_find_max_snr'],
            min_snr=self.par['reduce']['findobj']['ech_find_min_snr'],
            nabove_min_snr=self.par['reduce']['findobj']['ech_find_nabove_min_snr'],
            show_trace=show_trace, debug=debug)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.sciImg.fullmask == 0), chname='ech_objfind',sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech), skymask


    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, global_sky, sobjs,
                             spat_pix=None, model_noise=True, min_snr=2.0, fit_fwhm=False,
                             show_profile=False, show_resids=False, show_fwhm=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Args:
            global_sky (np.ndarray):
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
            spat_pix (np.ndarray, optional):
            model_noise (bool, optional):
            show_resids (bool, optional):
            show_profile (bool, optional):
            show (bool, optional):

        Returns:
            tuple: skymodel (np.ndarray), objmodel (np.ndarray), ivarmodel (np.ndarray), outmask (np.ndarray), sobjs

        """
        self.global_sky = global_sky

        # Pulled out some parameters to make the method all easier to read
        bsp = self.par['reduce']['skysub']['bspline_spacing']
        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        box_rad_order = self.par['reduce']['extraction']['boxcar_radius']/plate_scale
        sigrej = self.par['reduce']['skysub']['sky_sigrej']
        sn_gauss = self.par['reduce']['extraction']['sn_gauss']
        model_full_slit = self.par['reduce']['extraction']['model_full_slit']

        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = skysub.ech_local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                                  self.sciImg.fullmask, self.tilts, self.waveimg,
                                                  self.global_sky, self.sciImg.rn2img,
                                                  self.slits_left, self.slits_right,
                                                  self.slitmask,
                                                  sobjs, self.order_vec, spat_pix=spat_pix,
                                                  std=self.std_redux, fit_fwhm=fit_fwhm,
                                                  min_snr=min_snr, bsp=bsp,
                                                  box_rad_order=box_rad_order, sigrej=sigrej,
                                                  sn_gauss=sn_gauss,
                                                  model_full_slit=model_full_slit,
                                                  model_noise=model_noise,
                                                  show_profile=show_profile,
                                                  show_resids=show_resids, show_fwhm=show_fwhm)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True, chname='ech_local')
            self.show('resid', sobjs = self.sobjs, slits= True, chname='ech_resid')

        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

