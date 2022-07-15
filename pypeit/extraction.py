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

from pypeit import msgs, utils
from pypeit.display import display
from pypeit.core import skysub, extract, wave, flexure
from pypeit.core.moment import moment1d

from linetools.spectra import xspectrum1d

from IPython import embed


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
        outmask (`numpy.ndarray`_):
            Final output mask
        extractmask (`numpy.ndarray`_):
            Extraction mask
        slits (:class:`pypeit.slittrace.SlitTraceSet`):
        sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
            Only object finding but no extraction
        sobjs (:class:`pypeit.specobjs.SpecObjs`):
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

    """

    __metaclass__ = ABCMeta

    # Superclass factory method generates the subclass instance
    @classmethod
    def get_instance(cls, sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, bkg_redux=False,
                     return_negative=False, std_redux=False, show=False, basename=None):
        """
        Instantiate the Reduce subclass appropriate for the provided
        spectrograph.

        The class must be subclassed from Reduce.  See :class:`Reduce` for
        the description of the valid keyword arguments.

        Args:
            sciImg (:class:`~pypeit.images.scienceimage.ScienceImage`):
                Image to reduce.
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                Objects found but not yet extracted
            spectrograph (:class:`~pypeit.spectrographs.spectrograph.Spectrograph`):
            par (pypeit.par.pyepeitpar.PypeItPar):
            caliBrate (:class:`pypeit.calibrations.Calibrations`):
            objtype (:obj:`str`):
                Specifies object being reduced 'science' 'standard'
                'science_coadd2d'.  This is used only to determine the
                spat_flexure_shift and ech_order for coadd2d.
            bkg_redux (:obj:`bool`, optional):
                If True, the sciImg has been subtracted by
                a background image (e.g. standard treatment in the IR)
            return_negative (:obj:`bool`, optional):
                If True, negative objects from difference imaging will also be extracted and returned. Default=False.
                This option only applies to the case where bkg_redux=True, i.e. typically a near-IR reduction
                where difference imaging has been employed to perform a first-pass at sky-subtraction. The
                default behavior is to not extract these objects, although they are masked in global sky-subtraction
                (performed in the find_objects class), and modeled in local sky-subtraction (performed by this class).

            std_redux (:obj:`bool`, optional):
                If True the object being extracted is a standards star
                so that the reduction parameters can be adjusted accordingly.
            basename (str, optional):
                Output filename used for spectral flexure QA
            show (:obj:`bool`, optional):
                Show plots along the way?

        Returns:
            :class:`~pypeit.extraction.Extract`:
        """
        return next(c for c in utils.all_subclasses(Extract)
                    if c.__name__ == (spectrograph.pypeline + 'Extract'))(
                            sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, 
                            bkg_redux=bkg_redux, return_negative=return_negative, std_redux=std_redux, show=show,
                            basename=basename)

    def __init__(self, sciImg, sobjs_obj, spectrograph, par, caliBrate,
                 objtype, bkg_redux=False, return_negative=False, std_redux=False, show=False,
                 basename=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.sobjs_obj = sobjs_obj.copy() # This guarantees that stuff below does not mute this input variable
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        self.caliBrate = caliBrate
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

        # Initialise the slits
        msgs.info("Initialising slits")
        self.initialise_slits()

        # Internal bpm mask
        self.extract_bpm = (self.slits.mask > 0) & (np.invert(self.slits.bitmask.flagged(
                        self.slits.mask, flag=self.slits.bitmask.exclude_for_reducing)))
        self.extract_bpm_init = self.extract_bpm.copy()

        # These may be None (i.e. COADD2D)
        self.waveTilts = caliBrate.wavetilts
        self.wv_calib = caliBrate.wv_calib

        # Load up other input items
        self.bkg_redux = bkg_redux
        self.return_negative=return_negative

        self.std_redux = std_redux
        self.det = caliBrate.det
        self.binning = caliBrate.binning
        self.pypeline = spectrograph.pypeline
        self.extract_show = show

        self.steps = []

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.global_sky = None
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs = None  # Final extracted object list with trace corrections applied
        self.slitshift = np.zeros(self.slits.nslits)  # Global spectral flexure slit shifts (in pixels) that are applied to all slits.
        self.vel_corr = None

        # remove objects found in `BOXSLIT` (we don't want to extract those)
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

    @property
    def nobj_to_extract(self):
        """
        Number of objects to extract. Defined in children
        Returns:

        """
        return None

    def initialise_slits(self, initial=False):
        """
        Gather all the :class:`SlitTraceSet` attributes
        that we'll use here in :class:`Extract`

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

    def extract(self, global_sky, model_noise=None, spat_pix=None):
        """
        Main method to extract spectra from the ScienceImage

        Args:
            global_sky (`numpy.ndarray`_):
                Sky estimate
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                List of SpecObj that have been found and traced
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
            msgs.info("Skipping optimal extraction")

            # This will hold the extracted objects
            self.sobjs = self.sobjs_obj.copy()

            # Quick loop over the objects
            for iobj in range(self.sobjs.nobj):
                sobj = self.sobjs[iobj]
                # True  = Good, False = Bad for inmask
                thismask = self.slitmask == sobj.SLITID  # pixels for this slit
                inmask = self.sciImg.select_flag(invert=True) & thismask
                # Do it
                extract.extract_boxcar(self.sciImg.image, self.sciImg.ivar, inmask, self.waveimg,
                                       global_sky, sobj, base_var=self.sciImg.base_var,
                                       count_scale=self.sciImg.img_scale,
                                       noise_floor=self.sciImg.noise_floor)

            # Fill up extra bits and pieces
            self.objmodel = np.zeros_like(self.sciImg.image)
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # NOTE: fullmask is a bit mask, make sure it's treated as such, not
            # a boolean (e.g., bad pixel) mask.
            self.outmask = self.sciImg.fullmask
            self.skymodel = global_sky.copy()
        else:  # Local sky subtraction and optimal extraction.
            model_noise_1 = not self.bkg_redux if model_noise is None else model_noise
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
                self.local_skysub_extract(global_sky, self.sobjs_obj,
                                          model_noise=model_noise_1,
                                          spat_pix = spat_pix,
                                          show_profile=self.extract_show,
                                          show=self.extract_show)

        # Remove sobjs that don't have both OPT_COUNTS and BOX_COUNTS
        remove_idx = []
        for idx, sobj in enumerate(self.sobjs):
            # Find them
            if sobj.OPT_COUNTS is None and sobj.BOX_COUNTS is None:
                remove_idx.append(idx)
                msgs.warn(f'Removing object at pixel {sobj.SPAT_PIXPOS} because '
                          f'both optimal and boxcar extraction could not be performed')
            elif sobj.OPT_COUNTS is None:
                msgs.warn(f'Optimal extraction could not be performed for object at pixel {sobj.SPAT_PIXPOS}')

        # Remove them
        if len(remove_idx) > 0:
            self.sobjs.remove_sobj(remove_idx)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    def prepare_extraction(self):
        """ Prepare the masks and wavelength image for extraction.
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




    def run(self, global_sky, prepare_extraction=True, model_noise=None, spat_pix=None, ra=None, dec=None, obstime=None):
        """
        Primary code flow for PypeIt reductions

        *NOT* used by COADD2D

        Args:
            global_sky (`numpy.ndarray`_):
                Global sky model
            sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
                List of objects found during `run_objfind`
            prepare_extraction (bool):
                If True, generate the tilts image and the waveimg using fits that are stored in their
                respective objects. This allows spatial flexure compensation to optionally be taken into account.
                This option is set to False for coadd2d, since
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
            ra (:obj:`float`, optional):
                Required if helio-centric correction is to be applied
            dec (:obj:`float`, optional):
                Required if helio-centric correction is to be applied
            obstime (:obj:`astropy.time.Time`, optional):
                Required if helio-centric correction is to be applied


            return_negative (:obj:`bool`, optional):
                Do you want to extract the negative objects?

        Returns:
            tuple: skymodel (ndarray), objmodel (ndarray), ivarmodel (ndarray),
               outmask (ndarray), sobjs (SpecObjs), waveimg (`numpy.ndarray`_),
               tilts (`numpy.ndarray`_).
               See main doc string for description

        """
        # Start by preparing some masks and the wavelength image, ready for extraction
        # TODO this should return things to make the control flow less opqaque.
        if prepare_extraction:
            self.prepare_extraction()

        self.global_sky = global_sky

        # Do we have any positive objects to proceed with?
        if self.nobj_to_extract > 0:
            # Apply a global flexure correction to each slit
            # provided it's not a standard star
            if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
                self.spec_flexure_correct(mode='global')

            # Extract + Return
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = self.extract(self.global_sky, model_noise=model_noise, spat_pix=spat_pix)

            if self.bkg_redux:
                self.sobjs.make_neg_pos() if self.return_negative else self.sobjs.purge_neg()
        else:  # No objects, pass back what we have
            # Apply a global flexure correction to each slit
            # provided it's not a standard star
            if self.par['flexure']['spec_method'] != 'skip' and not self.std_redux:
                self.spec_flexure_correct(mode='global')
            #Could have negative objects but no positive objects so purge them
            if self.bkg_redux:
                self.sobjs_obj.make_neg_pos() if self.return_negative else self.sobjs_obj.purge_neg()
            self.skymodel = global_sky 
            self.objmodel = np.zeros_like(self.sciImg.image)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # Set to the initial mask in case no objects were found
            # NOTE: fullmask is a bit mask, make sure it's treated as such, not
            # a boolean (e.g., bad pixel) mask.
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
        # TODO avoid modifying arguments to a class or function in place. If slits is mutable, it should be a return
        # value for the run function
        reduce_masked = np.where(np.invert(self.extract_bpm_init) & self.extract_bpm & (self.slits.mask > 2))[0]
        if len(reduce_masked) > 0:
            # TODO Change BADREDUCE to BADEXTRACT
            self.slits.mask[reduce_masked] = self.slits.bitmask.turn_on(
                self.slits.mask[reduce_masked], 'BADREDUCE')

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs, \
               self.waveimg, self.tilts

    def local_skysub_extract(self, global_sky, sobjs, model_noise=True, spat_pix=None,
                             show_profile=False, show_resids=False, show=False):
        """
        Dummy method for local sky-subtraction and extraction.

        Overloaded by class specific skysub and extraction.
        """
        return None, None, None, None, None

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
            return

        # Perform some checks
        if mode == "local" and sobjs is None:
            msgs.error("No spectral extractions provided for flexure, using slit center instead")
        elif mode not in ["local", "global"]:
            msgs.error("mode must be 'global' or 'local'. Assuming 'global'.")

        # Prepare a list of slit spectra, if required.
        if mode == "global":
            gd_slits = np.logical_not(self.extract_bpm)
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
            flex_list = flexure.spec_flexure_slit(self.slits, self.slits.slitord_id, self.extract_bpm,
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
            flex_list = flexure.spec_flexure_slit(self.slits, self.slits.slitord_id, self.extract_bpm,
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
                    if sobj is None or sobj['BOX_WAVE'] is None:  # Nothing extracted; only the trace exists
                        continue
                    # Interpolate
                    new_sky = sobj.apply_spectral_flexure(this_flex_dict['shift'][cntr],
                                                          this_flex_dict['sky_spec'][cntr])
                    flex_list[islit]['sky_spec'][cntr] = new_sky.copy()
                    cntr += 1

        # Save QA
        basename = f'{self.basename}_{mode}_{self.spectrograph.get_det_name(self.det)}'
        out_dir = os.path.join(self.par['rdx']['redux_path'], 'QA')
        flexure.spec_flexure_qa(self.slits.slitord_id, self.extract_bpm, basename, flex_list,
                                specobjs=sobjs, out_dir=out_dir)

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
        # Correct Telescope's motion
        refframe = self.par['calibrations']['wavelengths']['refframe']
        if refframe in ['heliocentric', 'barycentric'] \
                and self.par['calibrations']['wavelengths']['reference'] != 'pixel':
            msgs.info("Performing a {0} correction".format(self.par['calibrations']['wavelengths']['refframe']))
            # Calculate correction
            radec = ltu.radec_to_coord((ra, dec))
            vel, vel_corr = wave.geomotion_correct(radec, obstime,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   refframe)
            # Apply correction to objects
            msgs.info('Applying {0} correction = {1:0.5f} km/s'.format(refframe, vel))
            if (sobjs is not None) and (sobjs.nobj != 0):
                # Loop on slits to apply
                gd_slitord = self.slits.slitord_id[np.logical_not(self.extract_bpm)]
                for slitord in gd_slitord:
                    indx = sobjs.slitorder_indices(slitord)
                    this_specobjs = sobjs[indx]
                    # Loop on objects
                    for specobj in this_specobjs:
                        if specobj is None:
                            continue
                        specobj.apply_helio(vel_corr, refframe)

            # Apply correction to wavelength image
            self.vel_corr = vel_corr
            self.waveimg *= vel_corr

        else:
            msgs.info('A wavelength reference frame correction will not be performed.')

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

        if showmask:
            mask_in = self.sciImg.fullmask
            bitmask_in = self.sciImg.bitmask
        else:
            mask_in = None
            bitmask_in = None

        img_gpm = self.sciImg.select_flag(invert=True)
        detname = self.spectrograph.get_det_name(self.det)

        # TODO Do we still need this here?
        if attr == 'global' and all([a is not None for a in [self.sciImg.image, self.global_sky, self.sciImg.fullmask]]):
            # global sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - self.global_sky) * img_gpm.astype(float)
            mean, med, sigma = stats.sigma_clipped_stats(image[img_gpm], sigma_lower=5.0,
                                                         sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            ch_name = chname if chname is not None else f'global_sky_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                            mask=mask_in, clear=clear, wcs_match=True)
                                          #, cuts=(cut_min, cut_max))
        elif attr == 'local' and all([a is not None for a in [self.sciImg.image, self.skymodel, self.sciImg.fullmask]]):
            # local sky subtraction
            # sky subtracted image
            image = (self.sciImg.image - self.skymodel) * img_gpm.astype(float)
            mean, med, sigma = stats.sigma_clipped_stats(image[img_gpm], sigma_lower=5.0,
                                                         sigma_upper=5.0)
            cut_min = mean - 1.0 * sigma
            cut_max = mean + 4.0 * sigma
            ch_name = chname if chname is not None else f'local_sky_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                            mask=mask_in, clear=clear, wcs_match=True)
                                          #, cuts=(cut_min, cut_max))
        elif attr == 'sky_resid' and all([a is not None for a in [self.sciImg.image, self.skymodel, self.objmodel,
                                                                  self.ivarmodel, self.sciImg.fullmask]]):
            # sky residual map with object included
            image = (self.sciImg.image - self.skymodel) * np.sqrt(self.ivarmodel)
            image *= img_gpm.astype(float)
            ch_name = chname if chname is not None else f'sky_resid_{detname}'
            viewer, ch = display.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                            bitmask=bitmask_in, mask=mask_in, clear=clear,
                                            wcs_match=True)
        elif attr == 'resid' and all([a is not None for a in [self.sciImg.image, self.skymodel, self.objmodel,
                                                                      self.ivarmodel, self.sciImg.fullmask]]):
            # full residual map with object model subtractede
            # full model residual map
            image = (self.sciImg.image - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel)
            image *= img_gpm.astype(float)
            ch_name = chname if chname is not None else f'resid_{detname}'
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


class MultiSlitExtract(Extract):
    """
    Child of Reduce for Multislit and Longslit reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs):
        super().__init__(sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs)


    @property
    def nobj_to_extract(self):
        """
        See parent method docs.

        Returns:

        """
        return self.nsobj_to_extract

    # TODO: JFH Should we reduce the number of iterations for standards or
    # near-IR redux where the noise model is not being updated?
    def local_skysub_extract(self, global_sky, sobjs, spat_pix=None, model_noise=True,
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
        gdslits = np.where(np.invert(self.extract_bpm))[0]

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        # NOTE: fullmask is a bit mask, make sure it's treated as such, not a
        # boolean (e.g., bad pixel) mask.
        self.outmask = np.copy(self.sciImg.fullmask)
        # Initialize to input mask in case no objects were found
        self.extractmask = self.sciImg.select_flag(invert=True)
        # Initialize to zero in case no objects were found
        self.objmodel = np.zeros_like(self.sciImg.image)
        # Set initially to global sky in case no objects were found
        self.skymodel  = np.copy(self.global_sky)
        # Set initially to sciivar in case no obects were found.
        self.ivarmodel = np.copy(self.sciImg.ivar)

        # Could actually create a model anyway here, but probably
        # overkill since nothing is extracted
        self.sobjs = sobjs.copy()  # WHY DO WE CREATE A COPY HERE?

        base_gpm = self.sciImg.select_flag(invert=True)

        # Loop on slits
        for slit_idx in gdslits:
            slit_spat = self.slits.spat_id[slit_idx]
            msgs.info("Local sky subtraction and extraction for slit: {:d}".format(slit_spat))
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

            # Local sky subtraction and extraction
            self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], self.extractmask[thismask] \
                = skysub.local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                              self.tilts, self.waveimg, self.global_sky,
                                              thismask, self.slits_left[:,slit_idx],
                                              self.slits_right[:, slit_idx],
                                              self.sobjs[thisobj], ingpm=ingpm,
                                              spat_pix=spat_pix,
                                              model_full_slit=model_full_slit,
                                              sigrej=sigrej, model_noise=model_noise,
                                              std=self.std_redux, bsp=bsp,
                                              force_gauss=force_gauss, sn_gauss=sn_gauss,
                                              show_profile=show_profile,
                                              use_2dmodel_mask=use_2dmodel_mask,
                                              no_local_sky=no_local_sky,
                                              base_var=self.sciImg.base_var,
                                              count_scale=self.sciImg.img_scale,
                                              adderr=self.sciImg.noise_floor)

        # Set the bit for pixels which were masked by the extraction.
        # For extractmask, True = Good, False = Bad
        iextract = base_gpm & np.logical_not(self.extractmask)
        # TODO: Change this to use the update_mask method?
        self.outmask[iextract] = self.sciImg.bitmask.turn_on(self.outmask[iextract], 'EXTRACT')

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True)
            self.show('resid', sobjs = self.sobjs, slits= True)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs


class EchelleExtract(Extract):
    """
    Child of Reduce for Echelle reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs):
        super(EchelleExtract, self).__init__(sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
                            else self.slits.ech_order
#                            else self.spectrograph.order_vec(self.spatial_coo)
        if self.order_vec is None:
            msgs.error('Unable to set Echelle orders, likely because they were incorrectly '
                       'assigned in the relevant SlitTraceSet.')


    @property
    def nobj_to_extract(self):
        """
        See parent method docs.

        Returns:

        """

        norders = self.order_vec.size
        if (self.nsobj_to_extract % norders) == 0:
            return int(self.nsobj_to_extract/norders)
        else:
            msgs.error('Number of specobjs in sobjs is not an integer multiple of the number or ordres!')

    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, global_sky, sobjs,
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

        # Pulled out some parameters to make the method all easier to read
        bsp = self.par['reduce']['skysub']['bspline_spacing']
        sigrej = self.par['reduce']['skysub']['sky_sigrej']
        sn_gauss = self.par['reduce']['extraction']['sn_gauss']
        model_full_slit = self.par['reduce']['extraction']['model_full_slit']
        force_gauss = self.par['reduce']['extraction']['use_user_fwhm']


        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = skysub.ech_local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                                  self.sciImg.fullmask, self.tilts, self.waveimg,
                                                  self.global_sky, self.slits_left,
                                                  self.slits_right, self.slitmask, sobjs,
                                                  self.order_vec, spat_pix=spat_pix,
                                                  std=self.std_redux, fit_fwhm=fit_fwhm,
                                                  min_snr=min_snr, bsp=bsp, sigrej=sigrej,
                                                  force_gauss=force_gauss, sn_gauss=sn_gauss,
                                                  model_full_slit=model_full_slit,
                                                  model_noise=model_noise,
                                                  show_profile=show_profile,
                                                  show_resids=show_resids, show_fwhm=show_fwhm,
                                                  base_var=self.sciImg.base_var,
                                                  count_scale=self.sciImg.img_scale,
                                                  adderr=self.sciImg.noise_floor)
        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True, chname='ech_local')
            self.show('resid', sobjs = self.sobjs, slits= True, chname='ech_resid')

        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

# TODO Should this be removed? I think so.
class IFUExtract(MultiSlitExtract):
    """
    Child of Reduce for IFU reductions

    See parent doc string for Args and Attributes

    """
    def __init__(self, sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs):
        super(IFUExtract, self).__init__(sciImg, sobjs_obj, spectrograph, par, caliBrate, objtype, **kwargs)
        self.initialise_slits(initial=True)


