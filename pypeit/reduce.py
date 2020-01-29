
import inspect
import numpy as np

from astropy import stats
from abc import ABCMeta

from linetools import utils as ltu

from pypeit import specobjs
from pypeit import ginga, msgs, edgetrace
from pypeit.core import skysub, extract, pixels, wave

from IPython import embed

class Reduce(object):
    """
    This class will organize and run actions related to
    finding objects, sky subtraction, and extraction for
    a Science or Standard star exposure

    Args:
        sciImg (pypeit.images.scienceimage.ScienceImage):
        spectrograph (pypeit.spectrograph.Spectrograph):
        par (pypeit.par.pyepeitpar.PypeItPar):
        caliBrate (pypeit.calibrations.Calibrations):
           This is only used as a container and it must contain the main products
           of WaveTilts, WaveImage, and EdgeTrace
        det (int, optional):
           Detector indice
        setup (str, optional):
           Used for naming
        maskslits (ndarray, optional):
          Specifies masked out slits
          True = Masked
        objtype (str, optional):
           Specifies object being reduced 'science' 'standard' 'science_coadd2d'
        show (bool, optional):
           Show plots along the way?

    Attributes:
        ivarmodel (np.ndarray):
            Model of inverse variance
        objimage (np.ndarray):
            Model of object
        skyimage (np.ndarray):
            Final model of sky
        initial_sky (np.ndarray):
            Initial sky model after first pass with global_skysub()
        global_sky (np.ndarray):
            Fit to global sky
        skymask (np.ndarray):
            Mask of the sky fit
        outmask (np.ndarray):
            Final output mask
        extractmask (np.ndarray):
            Extraction mask
        sobjs_obj (:class:`pypeit.specobjs.SpecObjs`):
            Only object finding but no extraction
        sobjs (SpecObsj):
            Final extracted object list with trace corrections applied

    """

    __metaclass__ = ABCMeta

    def __init__(self, sciImg, spectrograph, par, caliBrate,
                 ir_redux=False, det=1, std_redux=False, show=False,
                 objtype='science', binning=None, setup=None, maskslits=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.objtype = objtype
        self.par = par
        #self.extraction_par = self.par['reduce']['extraction']
        #self.wave_par = self.par['calibrations']['wavelengths']
        #self.flex_par = self.par['flexure']
        # Parse
        self.caliBrate = caliBrate
        self.slits = self.caliBrate.slits
        self.tilts = self.caliBrate.tilts_dict['tilts']
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        # TODO: We keep creating this image...
        # NOTE: this uses the par defined by EdgeTraceSet; this will
        # use the tweaked traces if they exist
        self.slitmask = self.slits.slit_img()
        self.sciImg.update_mask_slitmask(self.slitmask)
        self.maskslits = self._get_goodslits(maskslits)
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
                '''
                if self.pypeline == 'Echelle':
                    # Grab the positive object only
                    thisobj = (self.sobjs.ech_orderindx == iord) & (
                            self.sobjs.ech_objid > 0)  # pos indices of objects for this slit
                    sobj = self.sobjs[np.where(thisobj)[0][0]]
                    # Plate scale
                    plate_scale = self.spectrograph.order_platescale(sobj.ECH_ORDER,
                                                                 binning=self.binning)[0]
                else:
                    sobj = self.sobjs[iord]
                    plate_scale = self.spectrograph.detector[self.det - 1]['platescale']
                '''
                # True  = Good, False = Bad for inmask
                thismask = (self.slitmask == iobj)  # pixels for this slit
                inmask = (self.sciImg.mask == 0) & thismask
                # Do it
                extract.extract_boxcar(self.sciImg.image, self.sciImg.ivar,
                                               inmask, self.caliBrate.mswave,
                                               global_sky, self.sciImg.rn2img,
                                               self.par['reduce']['extraction']['boxcar_radius']/plate_scale,
                                               sobj)
            # Fill up extra bits and pieces
            self.objmodel = np.zeros_like(self.sciImg.image)
            self.ivarmodel = np.copy(self.sciImg.ivar)
            self.outmask = self.sciImg.mask
            self.skymodel = global_sky.copy()
        else:  # Local sky subtraction and optimal extraction.
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = \
                self.local_skysub_extract(self.caliBrate.mswave, global_sky, self.sobjs_obj,
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
               outmask (ndarray), sobjs (SpecObjs).  See main doc string for description

        """
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
                                  std_trace=std_trace, show=self.reduce_show,
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
            self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = self.extract(self.global_sky,
                                                                                                  self.sobjs_obj)
        else:  # No objects, pass back what we have
            self.skymodel = self.initial_sky
            self.objmodel = np.zeros_like(self.sciImg.image)
            # Set to sciivar. Could create a model but what is the point?
            self.ivarmodel = np.copy(self.sciImg.ivar)
            # Set to the initial mask in case no objects were found
            self.outmask = self.sciImg.mask
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
                self.flexure_correct(self.sobjs, basename)
            # Heliocentric
            radec = ltu.radec_to_coord((ra, dec))
            self.helio_correct(self.sobjs, radec, obstime)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

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
            self.show('image', image=image*(self.sciImg.mask == 0), chname='objfind',sobjs=sobjs_obj_single, slits=True)

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

        #self.maskslits = self.maskslits if maskslits is None else maskslits
        gdslits = np.where(np.invert(self.maskslits))[0]

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)

        # Select the edges to use: Selects the edges tweaked by the
        # illumination profile if they're present; otherwise, it
        # selects the original edges from EdgeTraceSet. To always
        # select the latter, use the method with `original=True`.
        left, right = self.slits.select_edges()

        # Loop on slits
        for slit in gdslits:
            msgs.info("Global sky subtraction for slit: {:d}".format(slit))
            thismask = (self.slitmask == slit)
            inmask = (self.sciImg.mask == 0) & thismask & skymask_now
            # Find sky
            self.global_sky[thismask] \
                    = skysub.global_skysub(self.sciImg.image, self.sciImg.ivar, self.tilts,
                                           thismask, left[:,slit], right[:,slit], inmask=inmask,
                                           sigrej=sigrej,
                                           bsp=self.par['reduce']['skysub']['bspline_spacing'],
                                           no_poly=self.par['reduce']['skysub']['no_poly'],
                                           pos_mask=(not self.ir_redux), show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.maskslits[slit] = True

        if update_crmask:
            self.sciImg.update_mask_cr(subtract_img=self.global_sky)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', slits=True, sobjs=sobjs_show, clear=False)

        # Return
        return self.global_sky

    def local_skysub_extract(self, waveimg, global_sky, sobjs,
                             model_noise=True, spat_pix=None,
                             show_profile=False, show_resids=False, show=False):
        """
        Dummy method for locak skysubtraction and extraction.

        Overloaded by class specific skysub and extraction.

        Returns:

        """

        return None, None, None, None, None

    def flexure_correct(self, sobjs, basename):
        """ Correct for flexure

        Spectra are modified in place (wavelengths are shifted)

        Wrapper to wave.flexure_obj()

        Args:
            sobjs (:class:`pypeit.specobjs.SpecObjs`):
            basename (str):

        """

        if self.par['flexure']['method'] != 'skip':
            flex_list = wave.flexure_obj(sobjs, self.maskslits, self.par['flexure']['method'],
                                         self.par['flexure']['spectrum'],
                                         mxshft=self.par['flexure']['maxshift'])
            # QA
            wave.flexure_qa(sobjs, self.maskslits, basename, self.det, flex_list,out_dir=self.par['rdx']['redux_path'])
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
            vel, vel_corr = wave.geomotion_correct(sobjs, radec, obstime, self.maskslits,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   self.par['calibrations']['wavelengths']['frame'])
        else:
            msgs.info('A wavelength reference-frame correction will not be performed.')
            vel_corr = None

        return

    def _get_goodslits(self, maskslits):
        """
        Return the slits to be reduce by going through the maskslits
        logic below. If the input maskslits is None it uses previously
        assigned maskslits

        Args:
            maskslits (np.ndarray or None):

        Returns:
            np.ndarray: slit numbers to be reduced
        """

        # Identify the slits that we want to consider.
        if maskslits is not None:
            # If maskslits was passed in use it, and update self
            self.maskslits = maskslits
            return self.maskslits
        else:
            try:
                return self.maskslits
            except AttributeError:
                # If maskslits was not passed, and it does not exist in self, reduce all slits
                self.maskslits = np.zeros(self.slits.nslits, dtype=bool)
                return self.maskslits

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
            mask_in = self.sciImg.mask
            bitmask_in = self.sciImg.bitmask
        else:
            mask_in = None
            bitmask_in = None

        if attr == 'global':
            # global sky subtraction
            if self.sciImg.image is not None and self.global_sky is not None and self.sciImg.mask is not None:
                # sky subtracted image
                image = (self.sciImg.image - self.global_sky)*(self.sciImg.mask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.sciImg.mask == 0], sigma_lower=5.0,
                                                       sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                ch_name = chname if chname is not None else 'global_sky_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                              mask=mask_in, clear=clear, wcs_match=True)
                                              #, cuts=(cut_min, cut_max))
        elif attr == 'local':
            # local sky subtraction
            if self.sciImg.image is not None and self.skymodel is not None and self.sciImg.mask is not None:
                # sky subtracted image
                image = (self.sciImg.image - self.skymodel)*(self.sciImg.mask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.sciImg.mask == 0], sigma_lower=5.0,
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
                    and self.sciImg.mask is not None:
                image = (self.sciImg.image - self.skymodel) * np.sqrt(self.ivarmodel)
                image *= (self.sciImg.mask == 0)
                ch_name = chname if chname is not None else 'sky_resid_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                              bitmask=bitmask_in, mask=mask_in, clear=clear,
                                              wcs_match=True)
        elif attr == 'resid':
            # full residual map with object model subtractede
            if self.sciImg.image is not None and self.skymodel is not None \
                    and self.objmodel is not None and self.ivarmodel is not None \
                    and self.sciImg.mask is not None:
                # full model residual map
                image = (self.sciImg.image - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel)
                image *= (self.sciImg.mask == 0)
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
                ginga.show_trace(viewer, ch, spec.TRACE_SPAT, spec.name, color=color)

        if slits:
            if self.slits is not None:
                # TODO: IDs are always set by the original edge traces
                # produced by EdgeTraceSet, not the tweaked ones
                # produced by FlatField. Is that the desired behavior?
                left, right = self.slits.select_edges()
                ginga.show_slits(viewer, ch, left, right, self.slits.id)

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
    def __init__(self, sciImg, spectrograph, par, caliBrate, **kwargs):
        super(MultiSlitReduce, self).__init__(sciImg, spectrograph, par, caliBrate, **kwargs)

    def get_platescale(self, dummy):
        """
        Return the platescale for multislit.
        The input argument is ignored

        Args:
            dummy (:class:`pypeit.specobj.SpecObj`):
                ignored

        Returns:
            float:

        """
        plate_scale = self.spectrograph.detector[self.det - 1]['platescale']
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
        gdslits = np.where(np.invert(self.maskslits))[0]

        # create the ouptut image for skymask
        skymask = np.zeros_like(image, dtype=bool)
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Select the edges to use: Selects the edges tweaked by the
        # illumination profile if they're present; otherwise, it
        # selects the original edges from EdgeTraceSet. To always
        # select the latter, use the method with `original=True`.
        left, right = self.slits.select_edges()

        # Loop on slits
        for slit in gdslits:
            qa_title ="Finding objects on slit # {:d}".format(slit)
            msgs.info(qa_title)
            thismask = (self.slitmask == slit)
            inmask = (self.sciImg.mask == 0) & thismask
            # Find objects
            specobj_dict = {'setup': self.setup, 'slitid': slit, #'orderindx': 999,
                            'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}

            # TODO we need to add QA paths and QA hooks. QA should be
            # done through objfind where all the relevant information
            # is. This will be a png file(s) per slit.

            sobjs_slit, skymask[thismask] = \
                extract.objfind(image, thismask, left[:,slit], right[:,slit], inmask=inmask,
                                ir_redux=self.ir_redux,
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
            self.show('image', image=image*(self.sciImg.mask == 0), chname = 'objfind',
                      sobjs=sobjs, slits=True)

        # Return
        return sobjs, len(sobjs), skymask


    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, waveimg, global_sky, sobjs,
                             spat_pix=None, model_noise=True,
                             show_resids=False,
                             show_profile=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Args:
            waveimg (np.ndarray):
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
        self.waveimg = waveimg
        self.global_sky = global_sky

        # get the good slits and assign self.maskslits
        gdslits = np.where(np.invert(self.maskslits))[0]

        # Select the edges to use: Selects the edges tweaked by the
        # illumination profile if they're present; otherwise, it
        # selects the original edges from EdgeTraceSet. To always
        # select the latter, use the method with `original=True`.
        left, right = self.slits.select_edges()

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        self.outmask = np.copy(self.sciImg.mask)
        # Initialize to input mask in case no objects were found
        self.extractmask = (self.sciImg.mask == 0)
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
        for slit in gdslits:
            msgs.info("Local sky subtraction and extraction for slit: {:d}".format(slit))
            thisobj = (self.sobjs.SLITID == slit) # indices of objects for this slit
            if np.any(thisobj):
                thismask = (self.slitmask == slit) # pixels for this slit
                # True  = Good, False = Bad for inmask
                inmask = (self.sciImg.mask == 0) & thismask
                # Local sky subtraction and extraction
                self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], \
                    self.extractmask[thismask] = skysub.local_skysub_extract(
                    self.sciImg.image, self.sciImg.ivar, self.tilts, self.waveimg, self.global_sky, self.sciImg.rn2img,
                    thismask, left[:,slit], right[:, slit],
                    self.sobjs[thisobj], spat_pix=spat_pix,
                    model_full_slit=self.par['reduce']['extraction']['model_full_slit'],
                    box_rad=self.par['reduce']['extraction']['boxcar_radius']/self.get_platescale(0), #self.spectrograph.detector[self.det-1]['platescale'],
                    sigrej=self.par['reduce']['skysub']['sky_sigrej'],
                    model_noise=model_noise, std=self.std_redux,
                    bsp=self.par['reduce']['skysub']['bspline_spacing'],
                    sn_gauss=self.par['reduce']['extraction']['sn_gauss'],
                    inmask=inmask, show_profile=show_profile)

        # Set the bit for pixels which were masked by the extraction.
        # For extractmask, True = Good, False = Bad
        iextract = (self.sciImg.mask == 0) & (self.extractmask == False)
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
    def __init__(self, sciImg, spectrograph, par, caliBrate, **kwargs):
        super(EchelleReduce, self).__init__(sciImg, spectrograph, par, caliBrate, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        self.order_vec = spectrograph.orders if 'coadd2d' in self.objtype \
                            else self.spectrograph.order_vec(self.slits.spatial_coordinates())

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
        inmask = self.sciImg.mask == 0
        # Find objects
        specobj_dict = {'setup': self.setup, 'slitid': 999, #'orderindx': 999,
                        'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}

        # Select the edges to use: Selects the edges tweaked by the
        # illumination profile if they're present; otherwise, it
        # selects the original edges from EdgeTraceSet. To always
        # select the latter, use the method with `original=True`.
        left, right = self.slits.select_edges()

        # TODO This is a bad idea -- we want to find everything for standards
        #sig_thresh = 30.0 if std else self.redux_par['sig_thresh']
        sobjs_ech, skymask[self.slitmask > -1] = extract.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, left, right, self.order_vec, self.maskslits,
            spec_min_max=np.vstack((self.slits.specmin, self.slits.specmax)),
            inmask=inmask, ir_redux=self.ir_redux, ncoeff=self.par['reduce']['findobj']['trace_npoly'],
            hand_extract_dict=manual_extract_dict, plate_scale=plate_scale,
            std_trace=std_trace,
            specobj_dict=specobj_dict,sig_thresh=self.par['reduce']['findobj']['sig_thresh'],
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
            self.show('image', image=image*(self.sciImg.mask == 0), chname='ech_objfind',sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech), skymask


    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, waveimg, global_sky, sobjs,
                             spat_pix=None, model_noise=True, min_snr=2.0, fit_fwhm=False,
                             show_profile=False, show_resids=False, show_fwhm=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Args:
            waveimg (np.ndarray):
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
        self.waveimg = waveimg
        self.global_sky = global_sky

        # TODO: Is this already available from the __init__ or could it have been overwritten?
        self.slitmask = self.slits.slit_img()

        # Select the edges to use: Selects the edges tweaked by the
        # illumination profile if they're present; otherwise, it
        # selects the original edges from EdgeTraceSet. To always
        # select the latter, use the method with `original=True`.
        left, right = self.slits.select_edges()

        # Pulled out some parameters to make the method all easier to read
        bsp = self.par['reduce']['skysub']['bspline_spacing']
        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        box_rad_order = self.par['reduce']['extraction']['boxcar_radius']/plate_scale
        sigrej = self.par['reduce']['skysub']['sky_sigrej']
        sn_gauss = self.par['reduce']['extraction']['sn_gauss']
        model_full_slit = self.par['reduce']['extraction']['model_full_slit']

        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs \
                = skysub.ech_local_skysub_extract(self.sciImg.image, self.sciImg.ivar,
                                                  self.sciImg.mask, self.tilts, self.waveimg,
                                                  self.global_sky, self.sciImg.rn2img,
                                                  self.slits.nslits, left, right, self.slitmask,
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


def instantiate_me(sciImg, spectrograph, par, caliBrate, **kwargs):
    """
    Instantiate the Reduce subclass appropriate for the provided
    spectrograph.

    The class must be subclassed from Reduce.  See :class:`Reduce` for
    the description of the valid keyword arguments.

    Args:
        sciImg (pypeit.images.scienceimage.ScienceImage):
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
        par (pypeit.par.pyepeitpar.PypeItPar):
        caliBrate (pypeit.calibrations.Calibrations):
        **kwargs
            Passed to Parent init

    Returns:
        :class:`pypeit.reduce.Reduce`:
    """
    indx = [c.__name__ == spectrograph.pypeline+'Reduce' for c in Reduce.__subclasses__()]
    if not np.any(indx):
        msgs.error('Pipeline {0} is not defined!'.format(spectrograph.pypeline))
    return Reduce.__subclasses__()[np.where(indx)[0][0]](sciImg, spectrograph,
                                                         par, caliBrate, **kwargs)



