
import inspect
import numpy as np

from astropy import stats
from abc import ABCMeta

from pypeit import specobjs
from pypeit import ginga, msgs, edgetrace
from pypeit.core import skysub, extract, pixels, wave

from IPython import embed

class Reduce(object):
    """
     This class will organize and run actions related to
     a Science or Standard star exposure

     Args:
         file_list : list
           List of raw files to produce the flat field
         spectrograph : str
         tslits_dict : dict
           dict from TraceSlits class
         par :
         tilts : ndarray
           tilts from WaveTilts class
           used for sky subtraction and object finding
         det : int
         setup : str
         datasec_img : ndarray
           Identifies pixels to amplifiers
         bpm : ndarray
           Bad pixel mask
         maskslits : ndarray (bool)
           Specifies masked out slits
         pixlocn : ndarray
         objtype : str
           'science'
           'standard'
           'science_coadd2d'
         scidx : int
           Row in the fitstbl corresponding to the exposure

     Attributes:
         frametype : str
           Set to 'science'
         sciframe : ndarray
           Processed 2D frame
         rawvarframe : ndarray
           Variance generated without a sky (or object) model
         modelvarframe : ndarray
           Variance generated with a sky model
         finalvar : ndarray
           Final variance frame
         global_sky : ndarray
           Sky model across the slit/order
         skycorr_box : ndarray
           Local corrections to the sky model
         final_sky : ndarray
           Final sky model; may include 'local' corrections
         obj_model : ndarray
           Model of the object flux
         trcmask : ndarray
           Masks of objects for sky subtraction
         tracelist : list
           List of traces for objects in slits
         inst_name : str
           Short name of the spectrograph, e.g. KASTb
         target_name : str
           Parsed from the Header
         basename : str
           Combination of camera, target, and time
           e.g. J1217p3905_KASTb_2015May20T045733.56
         time : Time
           time object
         specobjs : list
           List of specobjs
         bm: ScienceImageBitMask
           Object used to select bits of a given type
     """

    __metaclass__ = ABCMeta

    def __init__(self, sciImg, spectrograph, tslits_dict, par, tilts,
                 ir_redux=False, det=1,
                 objtype='science', binning=None, setup=None, maskslits=None):

        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!
        self.objtype = objtype
        self.par = par
        # TODO -- I find it very confusing to break apart the main parset
        self.proc_par = self.par['scienceframe']['process']
        # TODO Rename the scienceimage arset to reduce.
        self.redux_par = self.par['scienceimage']
        self.wave_par = self.par['calibrations']['wavelengths']
        self.flex_par = self.par['flexure']

        # Instantiation attributes for this object
        self.sciImg = sciImg
        self.spectrograph = spectrograph
        self.tslits_dict = tslits_dict
        self.slitmask = pixels.tslits2mask(self.tslits_dict)
        # Now add the slitmask to the mask (i.e. post CR rejection in proc)
        self.sciImg.update_mask_slitmask(self.slitmask)
        self.maskslits = self._get_goodslits(maskslits)
        self.ir_redux = ir_redux
        self.det = det
        self.tilts = tilts
        self.binning = binning
        self.setup = setup
        self.pypeline = spectrograph.pypeline

        self.steps = []

        # Other attributes that will be set later during object finding,
        # sky-subtraction, and extraction
        self.waveimage = None  # used by extract

        # Key outputs images for extraction
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.global_sky = None
        self.skymask = None
        self.outmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None  # Only object finding but no extraction
        self.sobjs = None  # Final extracted object list with trace corrections applied

    def _chk_objs(self, items):
        """

        Args:
            items:

        Returns:

        """
        for obj in items:
            if getattr(self, obj) is None:
                msgs.warn('You need to generate {:s} prior to this step..'.format(obj))
                if obj in ['sciimg', 'sciivar', 'rn2_img']:
                    msgs.warn('Run the process() method')
                elif obj in ['sobjs_obj']:
                    msgs.warn('Run the find_objects() method')
                elif obj in['global_sky']:
                    msgs.warn('Run the global_skysub() method')
                elif obj in ['tilts', 'tslits_dict'] :
                    msgs.warn('Calibrations missing: these were required to run find_objects() '
                              'and global_skysub()')
                elif obj in ['waveimg']:
                    msgs.warn('Calibrations missing: waveimg must be input as a parameter. Try '
                              'running calibrations')
                return False
        return True

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

    def find_objects(self, image, std=False, ir_redux=False, std_trace=None, maskslits=None,
                     show_peaks=False, show_fits=False, show_trace=False, show=False,
                     manual_extract_dict=None, debug=False):
        """

        Args:
            image:
            std:
            ir_redux:
            std_trace:
            maskslits:
            show_peaks:
            show_fits:
            show_trace:
            show:
            manual_extract_dict:
            debug:

        Returns:

        """

        # Positive image
        parse_manual = self.parse_manual_dict(manual_extract_dict, neg=False)
        sobjs_obj_init, nobj_init, skymask_pos = \
            self.find_objects_pypeline(image, std=std, ir_redux=ir_redux,
                                       std_trace=std_trace, maskslits=maskslits,
                                       show_peaks = show_peaks, show_fits = show_fits,
                                       show_trace = show_trace,
                                       manual_extract_dict=parse_manual, debug=debug)

        # For nobj we take only the positive objects
        if ir_redux:
            msgs.info("Finding objects in the negative image")
            # Parses
            parse_manual = self.parse_manual_dict(manual_extract_dict, neg=True)
            sobjs_obj_init_neg, nobj_init_neg, skymask_neg = \
                self.find_objects_pypeline(-image, std=std, ir_redux=ir_redux,
                                           std_trace=std_trace, maskslits=maskslits,
                                           show_peaks=show_peaks, show_fits=show_fits,
                                           show_trace=show_trace,
                                           manual_extract_dict=parse_manual,
                                           debug=debug)
            # Mask
            skymask = skymask_pos & skymask_neg
            # Add (if there are any)
            sobjs_obj_init.append_neg(sobjs_obj_init_neg)
        else:
            skymask = skymask_pos

        if show:
            self.show('image', image=image*(self.sciImg.mask == 0), chname='objfind',sobjs=sobjs_obj_init, slits=True)

        # For nobj we take only the positive objects
        return sobjs_obj_init, nobj_init, skymask

    def find_objects_pypeline(self, image, std=False, ir_redux=False, std_trace=None, maskslits=None,
                              show_peaks=False, show_fits=False, show_trace=False, show=False, debug=False,
                              manual_extract_dict=None):

        """
         Dummy method for object finding. Overloaded by class specific object finding.

         Returns:

         """
        return None, None, None

    def global_skysub(self, std=False, skymask=None, update_crmask=True, maskslits=None, show_fit=False,
                      show=False, show_objs=False):
        """
        Perform global sky subtraction, slit by slit

        Wrapper to skysub.global_skysub

        Parameters
        ----------
        tslits_dict: dict
           Dictionary containing information on the slits traced for this image

        Optional Parameters
        -------------------
        bspline_spaceing: (float):
           Break-point spacing for bspline

        use_skymask: (bool, optional):
           Mask objects using self.skymask if object finding has been run
           (This requires they were found previously, i.e. that find_objects was already run)

        Returns:
            global_sky: (numpy.ndarray) image of the the global sky model
        """


        # Prep
        self.global_sky = np.zeros_like(self.sciImg.image)
        if std:
            sigrej = 7.0
            update_crmask = False
            if not self.redux_par['global_sky_std']:
                msgs.info('Skipping global sky-subtraction for standard star.')
                return self.global_sky
        else:
            sigrej = 3.0

        self.maskslits = self.maskslits if maskslits is None else maskslits
        gdslits = np.where(np.invert(self.maskslits))[0]


        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciImg.image, dtype=bool)
        # Loop on slits
        for slit in gdslits:
            msgs.info("Global sky subtraction for slit: {:d}".format(slit))
            thismask = (self.slitmask == slit)
            inmask = (self.sciImg.mask == 0) & thismask & skymask_now
            # Find sky
            self.global_sky[thismask] = skysub.global_skysub(self.sciImg.image,
                                                             self.sciImg.ivar,
                                                             self.tilts, thismask,
                                                             self.tslits_dict['slit_left'][:,slit],
                                                             self.tslits_dict['slit_righ'][:,slit],
                                                             inmask=inmask,
                                                             sigrej=sigrej,
                                                             bsp=self.redux_par['bspline_spacing'],
                                                             no_poly=self.redux_par['no_poly'],
                                                             pos_mask = (not self.ir_redux),
                                                             show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.maskslits[slit] = True

        if update_crmask:
            self.sciImg.update_mask_cr(subtract_img=self.global_sky)
            #self.crmask = procimg.lacosmic(self.det, self.sciimg-self.global_sky,
            #                               self.spectrograph.detector[self.det-1]['saturation'],
            #                               self.spectrograph.detector[self.det-1]['nonlinear'],
            #                               varframe=utils.calc_ivar(self.sciivar),
            #                               maxiter=self.proc_par['lamaxiter'],
            #                               grow=self.proc_par['grow'],
            #                               remove_compact_obj=self.proc_par['rmcompact'],
            #                               sigclip=self.proc_par['sigclip'],
            #                               sigfrac=self.proc_par['sigfrac'],
            #                               objlim=self.proc_par['objlim'])
            # Rebuild the mask with this new crmask
            #self.mask = procimg.update_mask_cr(self.sciImg.bitmask, self.mask, self.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', slits=True, sobjs =sobjs_show, clear=False)

        # Return
        return self.global_sky

    def local_skysub_extract(self, waveimg, global_sky, sobjs,
                             maskslits=None, model_noise=True, std=False,
                             show_profile=False, show_resids=False, show=False):

        """
         Dummy method for locak skysubtraction and extraction. Overloaded by class specific skysub and extraction.

         Returns:

         """

        return None, None, None, None, None


    def flexure_correct(self, sobjs, basename):
        """ Correct for flexure

        Spectra are modified in place (wavelengths are shifted)

        Args:
            sobjs: SpecObjs object
            maskslits: ndarray

        """

        if self.flex_par['method'] != 'skip':
            flex_list = wave.flexure_obj(sobjs, self.maskslits, self.flex_par['method'],
                                         self.flex_par['spectrum'],
                                         mxshft=self.flex_par['maxshift'])
            # QA
            wave.flexure_qa(sobjs, self.maskslits, basename, self.det, flex_list,out_dir=self.par['rdx']['redux_path'])
        else:
            msgs.info('Skipping flexure correction.')


    def helio_correct(self, sobjs, radec, obstime):
        """ Perform a heliocentric correction """
        # Helio, correct Earth's motion
        if (self.wave_par['frame'] in ['heliocentric', 'barycentric']) \
                and (self.wave_par['reference'] != 'pixel'):
            # TODO change this keyword to refframe instead of frame
            msgs.info("Performing a {0} correction".format(self.wave_par['frame']))
            vel, vel_corr = wave.geomotion_correct(sobjs, radec, obstime, self.maskslits,
                                                   self.spectrograph.telescope['longitude'],
                                                   self.spectrograph.telescope['latitude'],
                                                   self.spectrograph.telescope['elevation'],
                                                   self.wave_par['frame'])
        else:
            msgs.info('A wavelength reference-frame correction will not be performed.')
            vel_corr = None

        return



    def _get_goodslits(self, maskslits):
        """
        Return the slits to be reduce by going through the maskslits
        logic below. If the input maskslits is None it uses previously
        assigned maskslits

        Returns
        -------
        gdslits
            numpy array of slit numbers to be reduced
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
                self.maskslits = np.zeros(self.tslits_dict['slit_left'].shape[1], dtype=bool)
                return self.maskslits


    def show(self, attr, image=None, showmask=False, sobjs=None, chname=None, slits=False,clear=False):
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
            if self.tslits_dict is not None:
                slit_ids = [edgetrace.get_slitid(self.sciImg.mask.shape,
                                                 self.tslits_dict['slit_left'],
                                                 self.tslits_dict['slit_righ'], ii)[0]
                                    for ii in range(self.tslits_dict['slit_left'].shape[1])]
                ginga.show_slits(viewer, ch, self.tslits_dict['slit_left'],
                                 self.tslits_dict['slit_righ'], slit_ids)

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




class MultiSlit(Reduce):
    """
    Child of Reduce for Multislit and Longslit reductions

    """
    def __init__(self, sciImg, spectrograph, tslits_dict, par, tilts, **kwargs):
        super(MultiSlit, self).__init__(sciImg, spectrograph, tslits_dict, par, tilts, **kwargs)

    def find_objects_pypeline(self, image, std=False, ir_redux=False, std_trace=None, maskslits=None,
                              manual_extract_dict=None,
                              show_peaks=False, show_fits=False, show_trace=False,
                              show=False, debug=False):

        """
        Find objects in the slits. This is currently setup only for ARMS

        Wrapper to extract.objfind

        Parameters
        ----------
        tslits_dict: dict
           Dictionary containing information on the slits traced for this image

        Optional Parameters
        -------------------
        SHOW_PEAKS:  bool
          Generate QA showing peaks identified by object finding

        SHOW_FITS:  bool
          Generate QA  showing fits to traces

        SHOW_TRACE:  bool
          Generate QA  showing traces identified. Requires an open ginga RC modules window

        Returns
        -------
        specobjs : Specobjs object
            Container holding Specobj objects
        nobj:
            Number of objects identified
        self.skymask : ndarray
                Boolean image indicating which pixels are useful for global sky subtraction

        """
        self.maskslits = self.maskslits if maskslits is None else maskslits
        gdslits = np.where(np.invert(self.maskslits))[0]

        # create the ouptut image for skymask
        skymask = np.zeros_like(image, dtype=bool)
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

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
                extract.objfind(image, thismask, self.tslits_dict['slit_left'][:,slit],
                                self.tslits_dict['slit_righ'][:,slit], inmask=inmask, ir_redux=ir_redux,
                                ncoeff=self.redux_par['trace_npoly'], std_trace=std_trace,
                                sig_thresh=self.redux_par['sig_thresh'], hand_extract_dict=manual_extract_dict,
                                specobj_dict=specobj_dict, show_peaks=show_peaks,
                                show_fits=show_fits, show_trace=show_trace,
                                trim_edg=self.redux_par['find_trim_edge'],
                                cont_fit=self.redux_par['find_cont_fit'],
                                npoly_cont=self.redux_par['find_npoly_cont'],
                                fwhm=self.redux_par['find_fwhm'],
                                maxdev=self.redux_par['find_maxdev'],
                                qa_title=qa_title, nperslit=self.redux_par['maxnumber'], debug_all=debug)
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
                             spat_pix=None, maskslits=None, model_noise=True, std = False,
                             show_profile=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Parameters
        ----------
        sobjs: object
           Specobjs object containing Specobj objects containing information about objects found.
        waveimg: ndarray, shape (nspec, nspat)
           Wavelength map

        Optional Parameters
        -------------------


        Returns:
            global_sky: (numpy.ndarray) image of the the global sky model
        """
        self.waveimg = waveimg
        self.global_sky = global_sky

        # get the good slits and assign self.maskslits
        self.maskslits = self.maskslits if maskslits is None else maskslits
        gdslits = np.where(np.invert(self.maskslits))[0]

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
                    thismask, self.tslits_dict['slit_left'][:,slit], self.tslits_dict['slit_righ'][:, slit],
                    self.sobjs[thisobj], spat_pix=spat_pix, model_full_slit=self.redux_par['model_full_slit'],
                    box_rad=self.redux_par['boxcar_radius']/self.spectrograph.detector[self.det-1]['platescale'],
                    sigrej=self.redux_par['sky_sigrej'],
                    model_noise=model_noise, std=std, bsp=self.redux_par['bspline_spacing'],
                    sn_gauss=self.redux_par['sn_gauss'], inmask=inmask, show_profile=show_profile)

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


class Echelle(Reduce):
    """
    Child of Reduce for Echelle reductions

    """
    def __init__(self, sciImg, spectrograph, tslits_dict, par, tilts, **kwargs):
        super(Echelle, self).__init__(sciImg, spectrograph, tslits_dict, par, tilts, **kwargs)

        # JFH For 2d coadds the orders are no longer located at the standard locations
        if 'coadd2d' in self.objtype:
            self.order_vec = spectrograph.orders
        else:
            slit_spat_pos = edgetrace.slit_spat_pos(self.tslits_dict)
            self.order_vec = self.spectrograph.order_vec(slit_spat_pos)


    def find_objects_pypeline(self, image, std=False, ir_redux=False, std_trace = None, maskslits=None,
                              show=False, show_peaks=False, show_fits=False, show_trace = False, debug=False,
                              manual_extract_dict=None):

        # create the ouptut image for skymask
        skymask = np.zeros_like(image, dtype=bool)

        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        inmask = self.sciImg.mask == 0
        # Find objects
        specobj_dict = {'setup': self.setup, 'slitid': 999, #'orderindx': 999,
                        'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}
        # TODO This is a bad idea -- we want to find everything for standards
        #sig_thresh = 30.0 if std else self.redux_par['sig_thresh']
        sobjs_ech, skymask[self.slitmask > -1] = extract.ech_objfind(
            image, self.sciImg.ivar, self.slitmask, self.tslits_dict['slit_left'], self.tslits_dict['slit_righ'],
            spec_min_max=np.vstack((self.tslits_dict['spec_min'],self.tslits_dict['spec_max'])),
            inmask=inmask, ir_redux=ir_redux, ncoeff=self.redux_par['trace_npoly'], order_vec=self.order_vec,
            hand_extract_dict=manual_extract_dict, plate_scale=plate_scale, std_trace=std_trace,
            specobj_dict=specobj_dict,sig_thresh=self.redux_par['sig_thresh'], show_peaks=show_peaks, show_fits=show_fits,
            trim_edg=self.redux_par['find_trim_edge'], cont_fit=self.redux_par['find_cont_fit'],
            npoly_cont=self.redux_par['find_npoly_cont'], fwhm=self.redux_par['find_fwhm'],
            maxdev=self.redux_par['find_maxdev'], max_snr=self.redux_par['ech_find_max_snr'],
            min_snr=self.redux_par['ech_find_min_snr'], nabove_min_snr=self.redux_par['ech_find_nabove_min_snr'],
            show_trace=show_trace, debug=debug)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.sciImg.mask == 0), chname = 'ech_objfind',sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech), skymask


    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, waveimg, global_sky, sobjs,
                             spat_pix=None, model_noise=True, min_snr=2.0, std = False, fit_fwhm=False,
                             maskslits=None, show_profile=False, show_resids=False, show_fwhm=False, show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Parameters
        ----------
        sobjs: object
           Specobjs object containing Specobj objects containing information about objects found.
        waveimg: ndarray, shape (nspec, nspat)
           Wavelength map

        Optional Parameters
        -------------------

        Returns:
            global_sky: (numpy.ndarray) image of the the global sky model
        """

        self.waveimg = waveimg
        self.global_sky = global_sky

        plate_scale = self.spectrograph.order_platescale(self.order_vec, binning=self.binning)
        self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs = skysub.ech_local_skysub_extract(
            self.sciImg.image, self.sciImg.ivar, self.sciImg.mask, self.tilts, self.waveimg, self.global_sky,
            self.sciImg.rn2img, self.tslits_dict, sobjs, self.order_vec, spat_pix=spat_pix,
            std=std, fit_fwhm=fit_fwhm, min_snr=min_snr, bsp=self.redux_par['bspline_spacing'],
            box_rad_order=self.redux_par['boxcar_radius']/plate_scale,
            sigrej=self.redux_par['sky_sigrej'],
            sn_gauss=self.redux_par['sn_gauss'], model_full_slit=self.redux_par['model_full_slit'],
            model_noise=model_noise, show_profile=show_profile, show_resids=show_resids, show_fwhm=show_fwhm)


        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True, chname='ech_local')
            self.show('resid', sobjs = self.sobjs, slits= True, chname='ech_resid')

        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs



def instantiate_me(sciImg, spectrograph, tslits_dict, par, tilts, **kwargs):
    """
    Instantiate the Reduce subclass appropriate for the provided
    spectrograph.

    The class must be subclassed from Reduce.  See :class:`Reduce` for
    the description of the valid keyword arguments.

    Args:
        spectrograph
            (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The instrument used to collect the data to be reduced.

        tslits_dict: dict
            dictionary containing slit/order boundary information
        tilts (np.ndarray):

    Returns:
        :class:`PypeIt`: One of the classes with :class:`PypeIt` as its
        base.
    """
    indx = [ c.__name__ == spectrograph.pypeline for c in Reduce.__subclasses__() ]
    if not np.any(indx):
        msgs.error('Pipeline {0} is not defined!'.format(spectrograph.pypeline))
    return Reduce.__subclasses__()[np.where(indx)[0][0]](sciImg, spectrograph, tslits_dict, par, tilts, **kwargs)



