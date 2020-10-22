"""
Main driver class for skysubtraction and extraction

.. include common links, assuming primary doc root is up one directory
.. include:: ../include/links.rst

"""

import inspect
import numpy as np
import os

from astropy import stats
from abc import ABCMeta

from linetools import utils as ltu
from scipy import interpolate
from scipy.optimize import least_squares

from pypeit import specobjs
from pypeit import msgs, utils
from pypeit import masterframe, flatfield
from pypeit.display import display
from pypeit.core import skysub, extract, pixels, wave, flexure, flat
from pypeit.images import buildimage
from pypeit.core.moment import moment1d

from linetools.spectra import xspectrum1d

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
        std_outfile (str):
           Filename of the standard star output

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
                 binning=None, setup=None, std_outfile=None, basename=None):
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
            basename (str, optional):
                Output filename used for spectral flexure QA
            **kwargs
                Passed to Parent init

        Returns:
            :class:`pypeit.reduce.Reduce`:
        """
        return next(c for c in cls.__subclasses__()
                    if c.__name__ == (spectrograph.pypeline + 'Reduce'))(
            sciImg, spectrograph, par, caliBrate, objtype, ir_redux=ir_redux, det=det,
            std_redux=std_redux, show=show,binning=binning, setup=setup,
            std_outfile=std_outfile, basename=basename)

    def __init__(self, sciImg, spectrograph, par, caliBrate,
                 objtype, ir_redux=False, det=1, std_redux=False, show=False,
                 binning=None, setup=None, std_outfile=None, basename=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        self.caliBrate = caliBrate
        self.std_outfile = std_outfile
        self.scaleimg = np.array([1.0], dtype=np.float)  # np.array([1]) applies no scale
        self.basename = basename
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

        # Initialise the slits
        msgs.info("Initialising slits")
        self.initialise_slits()

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
        self.slitshift = np.zeros(self.slits.nslits)  # Global spectral flexure slit shifts (in pixels) that are applied to all slits.

    def initialise_slits(self, initial=False):
        """
        Initialise the slits

        Args:
            initial (:obj:`bool`, optional):
                Use the initial definition of the slits. If False,
                tweaked slits are used.
        """
        # Slits
        self.slits = self.caliBrate.slits
        # Select the edges to use
        self.slits_left, self.slits_right, _ \
                = self.slits.select_edges(initial=initial, flexure=self.spat_flexure_shift)

        # Slitmask
        self.slitmask = self.slits.slit_img(initial=initial, flexure=self.spat_flexure_shift,
                                            exclude_flag=self.slits.bitmask.exclude_for_reducing)
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.sciImg.update_mask_slitmask(self.slitmask)
#        # For echelle
#        self.spatial_coo = self.slits.spatial_coordinates(initial=initial, flexure=self.spat_flexure_shift)

    def parse_manual_dict(self, manual_dict, neg=False):
        """
        Parse the manual dict
        This method is here mainly to deal with negative images

        Args:
            manual_dict (dict or None):
            neg (bool, optional):
                Negative image

        Returns:
            None or dict:  None if no matches; dict if there are for manual extraction

        """
        if manual_dict is None:
            return None
        #
        dets = np.atleast_1d(manual_dict['hand_extract_det'])
        # Grab the ones we want
        gd_det = dets > 0
        if not neg:
            gd_det = np.invert(gd_det)
        # Any?
        if not np.any(gd_det):
            return manual_dict
        # Fill
        manual_extract_dict = {}
        for key in manual_dict.keys():
            sgn = 1
            if key == 'hand_extract_det':
                sgn = -1
            manual_extract_dict[key] = sgn*np.atleast_1d(manual_dict[key])[gd_det]
        # Return
        return manual_extract_dict

    def extract(self, global_sky, sobjs_obj):
        """
        Main method to extract spectra from the ScienceImage

        Args:
            global_sky (`numpy.ndarray`_):
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
            # Purge out the negative objects if this was a near-IR reduction unless negative objects are requested

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

    def run(self, ra=None, dec=None, obstime=None, std_trace=None, show_peaks=False, return_negative=False):
        """
        Primary code flow for PypeIt reductions

        *NOT* used by COADD2D

        Args:
            ra (str, optional):
                Required if helio-centric correction is to be applied
            dec (str, optional):
                Required if helio-centric correction is to be applied
            obstime (:obj:`astropy.time.Time`, optional):
                Required if helio-centric correction is to be applied
            std_trace (np.ndarray, optional):
                Trace of the standard star
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
        msgs.info("Generating tilts image")
        self.tilts = self.waveTilts.fit2tiltimg(self.slitmask, flexure=tilt_flexure_shift)

        # Wavelengths (on unmasked slits)
        msgs.info("Generating wavelength image")
        self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits, spat_flexure=self.spat_flexure_shift)

        # First pass object finding
        self.sobjs_obj, self.nobj, skymask_init = \
            self.find_objects(self.sciImg.image, std_trace=std_trace,
                              show_peaks=show_peaks,
                              show=self.reduce_show & (not self.std_redux),
                              manual_extract_dict=self.par['reduce']['extraction']['manual'].dict_for_objfind())

        # Check if the user wants to overwrite the skymask with a pre-defined sky regions file
        skymask_init, usersky = self.load_skyregions(skymask_init)

        # Global sky subtract
        self.initial_sky = self.global_skysub(skymask=skymask_init).copy()

        # Second pass object finding on sky-subtracted image
        if (not self.std_redux) and (not self.par['reduce']['findobj']['skip_second_find']):
            self.sobjs_obj, self.nobj, self.skymask = \
                self.find_objects(self.sciImg.image - self.initial_sky,
                                  std_trace=std_trace,
                                  show=self.reduce_show,
                                  show_peaks=show_peaks,
                                  manual_extract_dict=self.par['reduce']['extraction']['manual'].dict_for_objfind())
        else:
            msgs.info("Skipping 2nd run of finding objects")


        # Do we have any positive objects to proceed with?
        if self.nobj > 0:
            # Global sky subtraction second pass. Uses skymask from object finding
            if (self.std_redux or self.par['reduce']['extraction']['skip_optimal'] or
                    self.par['reduce']['findobj']['skip_second_find'] or usersky):
                self.global_sky = self.initial_sky.copy()
            else:
                self.global_sky = self.global_skysub(skymask=self.skymask, show=self.reduce_show)

            # Apply a global flexure correction to each slit
            # provided it's not a standard star
            if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
                self.spec_flexure_correct(mode='global')

            # Extract + Return
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = self.extract(self.global_sky, self.sobjs_obj)
            if self.ir_redux:
                self.sobjs.make_neg_pos() if return_negative else self.sobjs.purge_neg()
        else:  # No objects, pass back what we have
            # Apply a global flexure correction to each slit
            # provided it's not a standard star
            if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
                self.spec_flexure_correct(mode='global')
            #Could have negative objects but no positive objects so purge them
            if self.ir_redux:
                self.sobjs_obj.make_neg_pos() if return_negative else self.sobjs_obj.purge_neg()
            self.skymodel = self.initial_sky
            self.objmodel = np.zeros_like(self.sciImg.image)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # Set to the initial mask in case no objects were found
            self.outmask = self.sciImg.fullmask
            # empty specobjs object from object finding
            self.sobjs = self.sobjs_obj

        # If a global spectral flexure has been applied to all slits, store this correction as metadata in each specobj
        if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
            for iobj in range(self.sobjs.nobj):
                islit = self.slits.spatid_to_zero(self.sobjs[iobj].SLITID)
                self.sobjs[iobj].update_flex_shift(self.slitshift[islit], flex_type='global')

        # Correct for local spectral flexure
        if self.sobjs.nobj == 0:
            msgs.warn('No objects to extract!')
        elif self.par['flexure']['spec_method'] not in ['skip', 'slitcen'] and not self.std_redux:
            # Apply a refined estimate of the flexure to objects, and then apply reference frame correction to objects
            self.spec_flexure_correct(mode='local', sobjs=self.sobjs)

        # Apply a reference frame correction to each object and the waveimg
        self.refframe_correct(ra, dec, obstime, sobjs=self.sobjs)

        # Update the mask
        reduce_masked = np.where(np.invert(self.reduce_bpm_init) & self.reduce_bpm)[0]
        if len(reduce_masked) > 0:
            self.slits.mask[reduce_masked] = self.slits.bitmask.turn_on(
                self.slits.mask[reduce_masked], 'BADREDUCE')

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs, \
               self.scaleimg, self.waveimg, self.tilts

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

    def global_skysub(self, skymask=None, update_crmask=True, trim_edg=(3,3),
                      show_fit=False, show=False, show_objs=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Args:
            skymask (np.ndarray, None):
                A 2D image indicating sky regions (1=sky)
            update_crmask (bool, optional):
            show_fit (bool, optional):
            show (bool, optional):
            show_objs (bool, optional):

        Returns:
            `numpy.ndarray`_: image of the the global sky model

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
            msgs.info("Global sky subtraction for slit: {:d}".format(slit_idx))
            thismask = self.slitmask == slit_spat
            inmask = (self.sciImg.fullmask == 0) & thismask & skymask_now
            # All masked?
            if not np.any(inmask):
                msgs.warn("No pixels for fitting sky.  If you are using mask_by_boxcar=True, your radius may be too large.")
                self.reduce_bpm[slit_idx] = True
                continue

            # Find sky
            self.global_sky[thismask] = skysub.global_skysub(self.sciImg.image, self.sciImg.ivar, self.tilts,
                                                             thismask, self.slits_left[:,slit_idx],
                                                             self.slits_right[:,slit_idx],
                                                             inmask=inmask, sigrej=sigrej,
                                                             bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                                             no_poly=self.par['reduce']['skysub']['no_poly'],
                                                             pos_mask=(not self.ir_redux), show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.reduce_bpm[slit_idx] = True

        if update_crmask and self.par['scienceframe']['process']['mask_cr']:
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

    def load_skyregions(self, skymask_init):
        """
        Load or generate the sky regions

        Parameters
        ----------
        skymask_init :  `numpy.ndarray`_
            A boolean array of sky pixels (True is pixel is a sky region)

        Returns
        -------
        skymask_init :  `numpy.ndarray`_
            A boolean array of sky pixels (True is pixel is a sky region)
        usersky : bool
            If the user has defined the sky, set this variable to True (otherwise False).
        """
        usersky = False
        if self.par['reduce']['skysub']['load_mask']:
            # Check if a master Sky Regions file exists for this science frame
            file_base = os.path.basename(self.sciImg.files[0])
            prefix = os.path.splitext(file_base)
            if prefix[1] == ".gz":
                sciName = os.path.splitext(prefix[0])[0]
            else:
                sciName = prefix[0]

            # Setup the master frame name
            master_dir = self.caliBrate.master_dir
            master_key = self.caliBrate.fitstbl.master_key(0, det=self.det) + "_" + sciName

            regfile = masterframe.construct_file_name(buildimage.SkyRegions,
                                                      master_key=master_key,
                                                      master_dir=master_dir)
            # Check if a file exists
            if os.path.exists(regfile):
                msgs.info("Loading SkyRegions file for: {0:s} --".format(sciName) + msgs.newline() + regfile)
                skyreg = buildimage.SkyRegions.from_file(regfile)
                skymask_init = skyreg.image.astype(np.bool)
                usersky = True
            else:
                msgs.warn("SkyRegions file not found:" + msgs.newline() + regfile)
        elif self.par['reduce']['skysub']['user_regions'] is not None:
            if len(self.par['reduce']['skysub']['user_regions']) != 0:
                skyregtxt = self.par['reduce']['skysub']['user_regions']
                if type(skyregtxt) is list:
                    skyregtxt = ",".join(skyregtxt)
                msgs.info("Generating skysub mask based on the user defined regions   {0:s}".format(skyregtxt))
                # The resolution probably doesn't need to be a user parameter
                maxslitlength = np.max(self.slits_right-self.slits_left)
                # Get the regions
                status, regions = skysub.read_userregions(skyregtxt, self.slits.nslits, maxslitlength)
                # Generate image
                skymask_init = skysub.generate_mask(self.pypeline, regions, self.slits, self.slits_left, self.slits_right, spat_flexure=self.spat_flexure_shift)
                usersky = True
        return skymask_init, usersky

    def spec_flexure_correct(self, mode="local", sobjs=None):
        """ Correct for spectral flexure

        Spectra are modified in place (wavelengths are shifted)

        Args:
            mode (str):
                "local" - Use sobjs to determine flexure correction
                "global" - Use waveimg and global_sky to determine flexure correction at the centre of the slit
            sobjs (:class:`pypeit.specobjs.SpecObjs`, None):
                Spectrally extracted objects
        """
        if self.par['flexure']['spec_method'] == 'skip':
            msgs.info('Skipping flexure correction.')
        else:
            # Perform some checks
            if mode == "local" and sobjs is None:
                msgs.error("No spectral extractions provided for flexure, using slit center instead")
            elif mode not in ["local", "global"]:
                msgs.error("mode must be 'global' or 'local'. Assuming 'global'.")

            # Prepare a list of slit spectra, if required.
            if mode == "global":
                gd_slits = np.logical_not(self.reduce_bpm)
                # TODO :: Need to think about spatial flexure - is the appropriate spatial flexure already included in trace_spat via left/right slits?
                trace_spat = 0.5 * (self.slits_left + self.slits_right)
                trace_spec = np.arange(self.slits.nspec)
                slit_specs = []
                for ss in range(self.slits.nslits):
                    if not gd_slits[ss]:
                        slit_specs.append(None)
                        continue
                    slit_spat = self.slits.spat_id[ss]
                    thismask = (self.slitmask == slit_spat)
                    box_denom = moment1d(self.waveimg * thismask > 0.0, trace_spat[:, ss], 2, row=trace_spec)[0]
                    wghts = (box_denom + (box_denom == 0.0))
                    slit_sky = moment1d(self.global_sky * thismask, trace_spat[:, ss], 2, row=trace_spec)[0] / wghts
                    # Denom is computed in case the trace goes off the edge of the image
                    slit_wave = moment1d(self.waveimg * thismask, trace_spat[:, ss], 2, row=trace_spec)[0] / wghts
                    # TODO :: Need to remove this XSpectrum1D dependency - it is required in:  flexure.spec_flex_shift
                    slit_specs.append(xspectrum1d.XSpectrum1D.from_tuple((slit_wave, slit_sky)))

                # Measure flexure
                # If mode == global: specobjs = None and slitspecs != None
                flex_list = flexure.spec_flexure_slit(self.slits, self.slits.slitord_id, self.reduce_bpm,
                                                      self.par['flexure']['spectrum'],
                                                      method=self.par['flexure']['spec_method'],
                                                      mxshft=self.par['flexure']['spec_maxshift'],
                                                      specobjs=sobjs, slit_specs=slit_specs)

                # Store the slit shifts that were applied to each slit
                # These corrections are later needed so the specobjs metadata contains the total spectral flexure
                self.slitshift = np.zeros(self.slits.nslits)
                for islit in range(self.slits.nslits):
                    if (not gd_slits[islit]) or len(flex_list[islit]['shift']) == 0:
                        continue
                    self.slitshift[islit] = flex_list[islit]['shift'][0]
                # Apply flexure to the new wavelength solution
                msgs.info("Regenerating wavelength image")
                self.waveimg = self.wv_calib.build_waveimg(self.tilts, self.slits,
                                                           spat_flexure=self.spat_flexure_shift,
                                                           spec_flexure=self.slitshift)
            elif mode == "local":
                # Measure flexure:
                # If mode == local: specobjs != None and slitspecs = None
                flex_list = flexure.spec_flexure_slit(self.slits, self.slits.slitord_id, self.reduce_bpm,
                                                      self.par['flexure']['spectrum'],
                                                      method=self.par['flexure']['spec_method'],
                                                      mxshft=self.par['flexure']['spec_maxshift'],
                                                      specobjs=sobjs, slit_specs=None)
                # Apply flexure to objects
                for islit in range(self.slits.nslits):
                    i_slitord = self.slits.slitord_id[islit]
                    indx = sobjs.slitorder_indices(i_slitord)
                    this_specobjs = sobjs[indx]
                    this_flex_dict = flex_list[islit]
                    # Loop through objects
                    cntr = 0
                    for ss, sobj in enumerate(this_specobjs):
                        if sobj is None:
                            continue
                        if sobj['BOX_WAVE'] is None:  # Nothing extracted; only the trace exists
                            continue
                        # Interpolate
                        new_sky = sobj.apply_spectral_flexure(this_flex_dict['shift'][cntr],
                                                              this_flex_dict['sky_spec'][cntr])
                        flex_list[islit]['sky_spec'][cntr] = new_sky.copy()
                        cntr += 1

            # Save QA
            basename = "{0:s}_{1:s}".format(self.basename, mode)
            flexure.spec_flexure_qa(self.slits.slitord_id, self.reduce_bpm, basename, self.det, flex_list,
                                    specobjs=sobjs, out_dir=os.path.join(self.par['rdx']['redux_path'], 'QA'))

    def refframe_correct(self, ra, dec, obstime, sobjs=None):
        """ Correct the calibrated wavelength to the user-supplied reference frame

        Args:
            radec (astropy.coordiantes.SkyCoord):
                Sky Coordinate of the observation
            obstime (:obj:`astropy.time.Time`):
                Observation time
            sobjs (:class:`pypeit.specobjs.Specobjs`, None):
                Spectrally extracted objects

        """
        radec = ltu.radec_to_coord((ra, dec))
        # Correct Telescope's motion
        refframe = self.par['calibrations']['wavelengths']['refframe']
        if (refframe in ['heliocentric', 'barycentric']) \
                and (self.par['calibrations']['wavelengths']['reference'] != 'pixel'):
            msgs.info("Performing a {0} correction".format(self.par['calibrations']['wavelengths']['refframe']))
            # Calculate correction
            vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   refframe)
            # Apply correction to objects
            msgs.info('Applying {0} correction = {1:0.5f} km/s'.format(refframe, vel))
            if (sobjs is not None) and (sobjs.nobj != 0):
                # Loop on slits to apply
                gd_slitord = self.slits.slitord_id[np.logical_not(self.reduce_bpm)]
                for slitord in gd_slitord:
                    indx = sobjs.slitorder_indices(slitord)
                    this_specobjs = sobjs[indx]
                    # Loop on objects
                    for specobj in this_specobjs:
                        if specobj is None:
                            continue
                        specobj.apply_helio(vel_corr, refframe)

            # Apply correction to wavelength image
            self.waveimg *= vel_corr

        else:
            msgs.info('A wavelength reference frame correction will not be performed.')
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
                viewer, ch = display.show_image(image, chname=ch_name, bitmask=bitmask_in,
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
                viewer, ch = display.show_image(image, chname=ch_name, bitmask=bitmask_in,
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
                viewer, ch = display.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
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
                viewer, ch = display.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                                bitmask=bitmask_in, mask=mask_in, clear=clear,
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

        # Masking options
        if self.par['reduce']['skysub']['mask_by_boxcar']:
            boxcar_rad_skymask = self.par['reduce']['extraction']['boxcar_radius'] / self.get_platescale(None),
        else:
            boxcar_rad_skymask = None

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
                                boxcar_rad_skymask=boxcar_rad_skymask,
                                maxdev=self.par['reduce']['findobj']['find_maxdev'],
                                find_min_max=self.par['reduce']['findobj']['find_min_max'],
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
                    use_2dmodel_mask=self.par['reduce']['extraction']['use_2dmodel_mask'],
                    no_local_sky=self.par['reduce']['skysub']['no_local_sky'])

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
                            else self.slits.ech_order
#                            else self.spectrograph.order_vec(self.spatial_coo)
        if self.order_vec is None:
            msgs.error('Unable to set Echelle orders, likely because they were incorrectly '
                       'assigned in the relevant SlitTraceSet.')


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
        # TODO -- Eliminate this specobj_dict thing
        specobj_dict = {'SLITID': 999, #'orderindx': 999,
                        'DET': self.det, 'OBJTYPE': self.objtype, 'PYPELINE': self.pypeline}

        sobjs_ech, skymask[self.slitmask > -1] = extract.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, self.slits_left, self.slits_right,
            self.order_vec, self.reduce_bpm, det=self.det,
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
            skymask_by_boxcar=self.par['reduce']['skysub']['mask_by_boxcar'],
            boxcar_rad=self.par['reduce']['extraction']['boxcar_radius'],  # arcsec
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


class IFUReduce(MultiSlitReduce, Reduce):
    """
    Child of Reduce for IFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, spectrograph, par, caliBrate, objtype, **kwargs):
        super(IFUReduce, self).__init__(sciImg, spectrograph, par, caliBrate, objtype, **kwargs)
        self.initialise_slits(initial=True)

    def find_objects_pypeline(self, image, std_trace=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, debug=False,
                              manual_extract_dict=None):
        """
        See MultiSlitReduce for slit-based IFU reductions
        """
        if self.par['reduce']['cube']['slit_spec']:
            return super().find_objects_pypeline(image, std_trace=std_trace,
                                                 show_peaks=show_peaks, show_fits=show_fits, show_trace=show_trace,
                                                 show=show, debug=debug, manual_extract_dict=manual_extract_dict)
        else:
            return None, None, None

    def apply_relative_scale(self, scaleImg):
        """Apply a relative scale to the science frame (and correct the varframe, too)

         Args:
             scaleImg (np.ndarray):
                scale image to divide the science frame by
        """
        # Check that scaleimg is set to the correct shape
        if self.scaleimg.size == 1:
            self.scaleimg = np.ones_like(self.sciImg.image)
        # Correct the relative illumination of the science frame
        msgs.info("Correcting science frame for relative illumination")
        scaleFact = scaleImg + (scaleImg == 0)
        self.scaleimg *= scaleFact
        sciImg, varImg = flat.flatfield(self.sciImg.image.copy(), scaleFact, self.sciImg.fullmask,
                                        varframe=utils.inverse(self.sciImg.ivar.copy()))
        self.sciImg.image = sciImg.copy()
        self.sciImg.ivar = utils.inverse(varImg)
        return

    def illum_profile_spatial(self, skymask=None, trim_edg=(0, 0)):
        """ Calculate the residual spatial illumination profile using the sky regions.
        The redisual is calculated using the differential:
        correction = amplitude * (1 + spatial_shift * (dy/dx)/y)
        where y is the spatial profile determined from illumflat, and spatial_shift
        is the residual spatial flexure shift in units of pixels.

         Args:
             skymask (np.ndarray):
                Mask of sky regions where the spatial illumination will be determined
             trim_edg (tuple):
                A tuple of two ints indicated how much of the slit edges should be
                trimmed when fitting to the spatial profile.
        """

        msgs.info("Performing spatial sensitivity correction")
        # Setup some helpful parameters
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        hist_trim = 0  # Trim the edges of the histogram to take into account edge effects
        gpm = (self.sciImg.fullmask == 0)
        slitid_img_init = self.slits.slit_img(pad=0, initial=True, flexure=self.spat_flexure_shift)
        spatScaleImg = np.ones_like(self.sciImg.image)
        # For each slit, grab the spatial coordinates and a spline
        # representation of the spatial profile from the illumflat
        rawimg = self.sciImg.image.copy()
        numbins = int(np.max(self.slits.get_slitlengths(initial=True, median=True)))
        spatbins = np.linspace(0.0, 1.0, numbins + 1)
        spat_slit = 0.5 * (spatbins[1:] + spatbins[:-1])
        slitlength = np.median(self.slits.get_slitlengths(median=True))
        coeff_fit = np.zeros((self.slits.nslits, 2))
        for sl, slitnum in enumerate(self.slits.spat_id):
            msgs.info("Deriving spatial correction for slit {0:d}/{1:d}".format(sl + 1, self.slits.spat_id.size))
            # Get the initial slit locations
            onslit_b_init = (slitid_img_init == slitnum)

            # Synthesize ximg, and edgmask from slit boundaries. Doing this outside this
            # routine would save time. But this is pretty fast, so we just do it here to make the interface simpler.
            spatcoord, edgmask = pixels.ximg_and_edgemask(self.slits_left[:, sl], self.slits_right[:, sl],
                                                          onslit_b_init, trim_edg=trim_edg)

            # Make the model histogram
            xspl = np.linspace(0.0, 1.0, 10 * int(slitlength))  # Sub sample each pixel with 10 subpixels
            modspl = self.caliBrate.flatimages.illumflat_spat_bsplines[sl].value(xspl)[0]
            gradspl = interpolate.interp1d(xspl, np.gradient(modspl) / modspl, kind='linear', bounds_error=False,
                                           fill_value='extrapolate')

            # Ignore skymask
            coord_msk = onslit_b_init & gpm
            hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
            cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
            hist_slit_all = hist / (cntr + (cntr == 0))
            histmod, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=gradspl(spatcoord[coord_msk]))
            hist_model = histmod / (cntr + (cntr == 0))

            # Repeat with skymask
            coord_msk = onslit_b_init & gpm & skymask_now
            hist, _ = np.histogram(spatcoord[coord_msk], bins=spatbins, weights=rawimg[coord_msk])
            cntr, _ = np.histogram(spatcoord[coord_msk], bins=spatbins)
            hist_slit = hist / (cntr + (cntr == 0))

            # Prepare for fit - take the non-zero elements and trim slit edges
            if hist_trim == 0:
                ww = (hist_slit != 0)
                xfit = spat_slit[ww]
                yfit = hist_slit_all[ww]
                mfit = hist_model[ww]
            else:
                ww = (hist_slit[hist_trim:-hist_trim] != 0)
                xfit = spat_slit[hist_trim:-hist_trim][ww]
                yfit = hist_slit_all[hist_trim:-hist_trim][ww]
                mfit = hist_model[hist_trim:-hist_trim][ww]

            # Fit the function
            spat_func = lambda par, ydata, model: par[0]*(1 + par[1] * model) - ydata
            res_lsq = least_squares(spat_func, [np.median(yfit), 0.0], args=(yfit, mfit))
            spatnorm = spat_func(res_lsq.x, 0.0, gradspl(spatcoord[onslit_b_init]))
            spatnorm /= spat_func(res_lsq.x, 0.0, gradspl(0.5))
            # Set the scaling factor
            spatScaleImg[onslit_b_init] = spatnorm
            coeff_fit[sl, :] = res_lsq.x

        debug=False
        if debug:
            from matplotlib import pyplot as plt
            xplt = np.arange(24)
            plt.subplot(121)
            plt.plot(xplt[0::2], coeff_fit[::2, 0], 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 0], 'bx')
            plt.subplot(122)
            plt.plot(xplt[0::2], coeff_fit[::2, 1]/10, 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 1]/10, 'bx')
            plt.show()
            plt.imshow(spatScaleImg, vmin=0.99, vmax=1.01)
            plt.show()
            plt.subplot(133)
            plt.plot(xplt[0::2], coeff_fit[::2, 2], 'rx')
            plt.plot(xplt[1::2], coeff_fit[1::2, 2], 'bx')
            plt.show()
        # Apply the relative scale correction
        self.apply_relative_scale(spatScaleImg)

    def illum_profile_spectral(self, global_sky, skymask=None):
        """Calculate the residual spectral illumination profile using the sky regions.
        This uses the same routine as the flatfield spectral illumination profile.

         Args:
             global_sky (np.ndarray):
                Model of the sky
             skymask (np.ndarray, None):
                Mask of sky regions where the spatial illumination will be determined
        """
        trim = self.par['calibrations']['flatfield']['slit_trim']
        gpm = (self.sciImg.fullmask == 0)
        scaleImg = flatfield.illum_profile_spectral(self.sciImg.image.copy(), self.waveimg, self.slits,
                                                    model=global_sky, gpmask=gpm, skymask=skymask, trim=trim,
                                                    flexure=self.spat_flexure_shift)
        # Now apply the correction to the science frame
        self.apply_relative_scale(scaleImg)

    def joint_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                     show_fit=False, show=False, show_objs=False, adderr=0.01):
        """ Perform a joint sky model fit to the data. See Reduce.global_skysub()
        for parameter definitions.
        """
        msgs.info("Performing joint global sky subtraction")
        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        nslits = self.slits.spat_id.size
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        self.global_sky = np.zeros_like(self.sciImg.image)
        thismask = (self.slitmask > 0)
        inmask = ((self.sciImg.fullmask == 0) & thismask & skymask_now).astype(np.bool)
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
                return self.global_sky

        # Iterate to use a model variance image
        numiter = 4
        model_ivar = self.sciImg.ivar.copy()
        for nn in range(numiter):
            msgs.info("Performing iterative joint sky subtraction - ITERATION {0:d}/{1:d}".format(nn+1, numiter))
            self.global_sky[thismask] = skysub.global_skysub(self.sciImg.image, model_ivar, tilt_wave,
                                                             thismask, self.slits_left, self.slits_right, inmask=inmask,
                                                             sigrej=sigrej, trim_edg=trim_edg,
                                                             bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                                             no_poly=self.par['reduce']['skysub']['no_poly'],
                                                             pos_mask=(not self.ir_redux), show_fit=show_fit)
            # Update the ivar image used in the sky fit
            msgs.info("Updating sky noise model")
            var = np.abs(self.global_sky - np.sqrt(2.0) * np.sqrt(self.sciImg.rn2img)) + self.sciImg.rn2img
            var = var + adderr ** 2 * (np.abs(self.global_sky)) ** 2
            model_ivar = utils.inverse(var)

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
        return self.global_sky

    def global_skysub(self, skymask=None, update_crmask=True, trim_edg=(0,0),
                      show_fit=False, show=False, show_objs=False):
        """
        Perform global sky subtraction. This IFU-specific routine ensures that the
        edges of the slits are not trimmed, and performs a spatial and spectral
        correction using the sky spectrum, if requested. See Reduce.global_skysub()
        for parameter definitions.
        """
        # Generate a global sky sub for all slits separately
        global_sky_sep = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
                                              show_fit=show_fit, show=show, show_objs=show_objs)
        # If the joint fit or spec/spat sensitivity corrections are not being performed, return the separate slits sky
        if not self.par['reduce']['skysub']['joint_fit'] and \
                not self.par['scienceframe']['process']['use_specillum'] and \
                not self.par['scienceframe']['process']['use_illumflat']:
            return global_sky_sep

        # Do the spatial scaling first
        # if self.par['scienceframe']['process']['use_illumflat']:
        #     # Perform the correction
        #     self.illum_profile_spatial(skymask=skymask)
        #     # Re-generate a global sky sub for all slits separately
        #     global_sky_sep = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
        #                                           show_fit=show_fit, show=show, show_objs=show_objs)

        if self.par['scienceframe']['process']['use_specillum']:
            self.illum_profile_spectral(global_sky_sep, skymask=skymask)

        # Fit to the sky
        if self.par['reduce']['skysub']['joint_fit']:
            # Use sky information in all slits to perform a joint sky fit
            self.global_sky = self.joint_skysub(skymask=skymask, update_crmask=update_crmask, trim_edg=trim_edg,
                                                show_fit=show_fit, show=show, show_objs=show_objs)
        else:
            # Re-run global skysub on individual slits, with the science frame now scaled
            self.global_sky = Reduce.global_skysub(self, skymask=skymask, update_crmask=update_crmask,
                                                   trim_edg=trim_edg, show_fit=show_fit, show=show, show_objs=show_objs)

        debug = False
        if debug:
            embed()
            wavefull = np.linspace(3950, 4450, 10000)
            import matplotlib.pylab as pl
            from matplotlib import pyplot as plt
            colors = pl.cm.jet(np.linspace(0, 1, gdslits.size))
            plt.subplot(121)
            for sl, slit_idx in enumerate(gdslits):
                slit_spat = self.slits.spat_id[slit_idx]
                thismask = self.slitmask == slit_spat
                wav = self.waveimg[thismask]
                flx = global_sky_sep[thismask]
                argsrt = np.argsort(wav)
                spl = interpolate.interp1d(wav[argsrt], flx[argsrt], bounds_error=False)
                if sl == 0:
                    ref = spl(wavefull)
                    plt.plot(wavefull, ref / np.nanmedian(ref), color=colors[sl], linestyle=':')
                plt.plot(wavefull, spl(wavefull) / ref, color=colors[sl])
            plt.subplot(122)
            for sl, slit_idx in enumerate(gdslits):
                slit_spat = self.slits.spat_id[slit_idx]
                thismask = self.slitmask == slit_spat
                wav = self.waveimg[thismask]
                flx = self.global_sky[thismask]
                argsrt = np.argsort(wav)
                spl = interpolate.interp1d(wav[argsrt], flx[argsrt], bounds_error=False)
                if sl == 0:
                    ref = spl(wavefull)
                    plt.plot(wavefull, ref / np.nanmedian(ref), color=colors[sl], linestyle=':')
                plt.plot(wavefull, spl(wavefull) / ref, color=colors[sl])
                print(sl, np.median(spl(wavefull) / ref))
                # plt.plot(wavefull, spl(wavefull), color=colors[sl])

            plt.show()

        return self.global_sky
