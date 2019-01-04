# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

import time
import datetime

from astropy import stats

from pypeit import msgs
from pypeit import processimages
from pypeit import specobjs
from pypeit import utils
from pypeit import ginga
from pypeit.core import skysub
from pypeit.core import extract
from pypeit.core import trace_slits
from pypeit.par import pypeitpar
from pypeit.core import procimg

#
# from pypeit.bitmask import BitMask
#
#
# class ScienceImageBitMask(BitMask):
#     """
#     Define a bitmask used to set the reasons why each pixel in a science
#     image was masked.
#     """
#     def __init__(self):
#         # TODO:
#         #   - Can IVAR0 and IVAR_NAN be consolidated into a single bit?
#         #   - Is EXTRACT ever set?
#         mask = {       'BPM': 'Component of the instrument-specific bad pixel mask',
#                         'CR': 'Cosmic ray detected',
#                 'SATURATION': 'Saturated pixel',
#                  'MINCOUNTS': 'Pixel below the instrument-specific minimum counts',
#                   'OFFSLITS': 'Pixel does not belong to any slit',
#                     'IS_NAN': 'Pixel value is undefined',
#                      'IVAR0': 'Inverse variance is undefined',
#                   'IVAR_NAN': 'Inverse variance is NaN',
#                    'EXTRACT': 'Pixel masked during local skysub and extraction'
#                }
#         super(ScienceImageBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))

class ScienceImage(processimages.ProcessImages):
    """
    This class will organize and run actions related to
    a Science or Standard star exposure

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    tslits_dict : dict
      dict from TraceSlits class
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
    scidx : int
      Row in the fitstbl corresponding to the exposure

    Attributes
    ----------
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

    # Frametype is a class attribute
    frametype = 'science'

    # TODO: Merge into a single parset, one for procing, and one for scienceimage
    def __init__(self, tslits_dict, spectrograph, file_list, bg_file_list = [], ir_redux=False,
                 det=1, objtype='science', binning = None, setup=None,
                 par=None, frame_par=None):


        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!
        self.objtype = objtype
        self.par = pypeitpar.ScienceImagePar() if par is None else par
        self.frame_par = pypeitpar.FrameGroupPar(objtype) if frame_par is None else frame_par
        self.proc_par = self.frame_par['process']

        # Start up by instantiating the process images class for reading in the relevant science files
        processimages.ProcessImages.__init__(self, spectrograph, [], det=det,
                                             par=self.frame_par['process'])

        # Instantiation attributes for this object
        self.tslits_dict = tslits_dict
        self.spectrograph = spectrograph
        self.file_list = file_list
        self.nsci = len(file_list)
        self.bg_file_list = bg_file_list
        # Are we subtracing the sky using background frames? If yes, set ir_redux=True
        self.nbg = len(self.bg_file_list)
        self.ir_redux = ir_redux
        if self.ir_redux and self.nbg == 0:
            msgs.error('IR reductions require that bg files are specified')
        self.det = det
        self.binning = binning
        self.setup = setup
        self.pypeline = spectrograph.pypeline

        # Set some detector parameters that we will need
        self.saturation = self.spectrograph.detector[self.det - 1]['saturation']
        self.mincounts = self.spectrograph.detector[self.det - 1]['mincounts']

        # These attributes will be sert when the image(s) are processed
        self.bpm = None
        self.bias = None
        self.pixel_flat = None
        self.illum_flat = None

        self.steps = []
        # Other attributes that will be set later during object finding,
        # sky-subtraction, and extraction
        self.slitmask = self.spectrograph.slitmask(self.tslits_dict)
        self.tilts = None # used by extract
        self.mswave = None # used by extract
        self.maskslits = None # used in find_object and extract

        # Key outputs images for extraction
        self.sciimg = None
        self.sciivar = None
        self.ivarmodel = None
        self.objimage = None
        self.skyimage = None
        self.global_sky = None
        self.skymask = None
        self.outmask = None
        self.mask = None                        # The composite bit value array
#        self.bitmask = ScienceImageBitMask()    # The bit mask interpreter
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None # Only object finding but no extraction
        self.sobjs = None  # Final extracted object list with trace corrections applied

        # Other bookeeping internals
        self.crmask = None
        self.mask = None


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

    def find_objects(self, image, std=False, std_trace = None, maskslits=None, show_peaks=False,
                     show_fits=False, show_trace=False, show=False):
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

        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        # create the ouptut image for skymask
        skymask = np.zeros_like(self.sciimg,dtype=bool)
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit in gdslits:
            qa_title ="Finding objects on slit # {:d}".format(slit)
            msgs.info(qa_title)
            thismask = (self.slitmask == slit)
            inmask = (self.mask == 0) & (self.crmask == False) & thismask
            # Find objects
            specobj_dict = {'setup': self.setup, 'slitid': slit,
                            'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}

            # TODO we need to add QA paths and QA hooks. QA should be
            # done through objfind where all the relevant information
            # is. This will be a png file(s) per slit.
            sig_thresh = 30.0 if std else self.par['sig_thresh']
            sobjs_slit, skymask[thismask] = \
                extract.objfind(image, thismask, self.tslits_dict['lcen'][:,slit],self.tslits_dict['rcen'][:,slit],
                inmask=inmask, std_trace=std_trace, sig_thresh=sig_thresh, hand_extract_dict=self.par['manual'],specobj_dict=specobj_dict,
                show_peaks=show_peaks,show_fits=show_fits, show_trace=show_trace,qa_title=qa_title,
                nperslit=self.par['maxnumber'])
            sobjs.add_sobj(sobjs_slit)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.mask == 0), chname = 'objfind',
                      sobjs=sobjs, slits=True)

        # Return
        return sobjs, len(sobjs), skymask


    def find_objects_ech(self, image, std=False, snr_trim=False, std_trace = None, show=False, show_peaks=False, show_fits=False, show_trace = False, debug=False):

        # Did they run process?
        if not self._chk_objs(['sciivar']):
            msgs.error('All quantities necessary to run ech_objfind() have not been set.')

        # Check for global sky if it does not exist print out a warning

        # Somehow implmenent masking below? Not sure it is worth it
        #self.maskslits = self._get_goodslits(maskslits)
        #gdslits = np.where(~self.maskslits)[0]

        # create the ouptut image for skymask
        skymask = np.zeros_like(self.sciimg,dtype=bool)

        plate_scale = self.spectrograph.order_platescale(binning=self.binning)
        inmask = (self.mask == 0) & (self.crmask == False)
        # Find objects
        specobj_dict = {'setup': self.setup, 'slitid': 999,
                        'det': self.det, 'objtype': self.objtype, 'pypeline': self.pypeline}
        # ToDO implement parsets here!
        sig_thresh = 30.0 if std else self.par['sig_thresh']
        sobjs_ech, skymask[self.slitmask > -1] = \
            extract.ech_objfind(image, self.sciivar, self.slitmask, self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                snr_trim=snr_trim, inmask=inmask, plate_scale=plate_scale, std_trace=std_trace,
                                specobj_dict=specobj_dict,sig_thresh=sig_thresh,
                                show_peaks=show_peaks, show_fits=show_fits, show_trace=show_trace, debug=debug)



        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.mask == 0), chname = 'ech_objfind',sobjs=sobjs_ech, slits=False)

        return sobjs_ech, len(sobjs_ech), skymask


    def global_skysub(self, tilts, std = False, skymask=None, update_crmask=True, maskslits=None, show_fit=False,
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

        if std:
            sigrej = 7.0
            update_crmask = False
        else:
            sigrej = 3.0

        self.tilts = tilts
        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        # Prep
        self.global_sky = np.zeros_like(self.sciimg)

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask_now = skymask if (skymask is not None) else np.ones_like(self.sciimg, dtype=bool)
        # Loop on slits
        for slit in gdslits:
            msgs.info("Global sky subtraction for slit: {:d}".format(slit))
            thismask = (self.slitmask == slit)
            inmask = (self.mask == 0) & thismask & skymask_now
            # Find sky
            self.global_sky[thismask] = skysub.global_skysub(self.sciimg, self.sciivar,
                                                             self.tilts, thismask,
                                                             self.tslits_dict['lcen'][:,slit],
                                                             self.tslits_dict['rcen'][:,slit],
                                                             inmask=inmask,
                                                             sigrej=sigrej,
                                                             bsp=self.par['bspline_spacing'],
                                                             pos_mask = (not self.ir_redux),
                                                             show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.maskslits[slit] = True

        if update_crmask:
            self.crmask = self.build_crmask(self.sciimg - self.global_sky, self.proc_par,
                                                                   self.det, self.spectrograph, ivar = self.sciivar,
                                                                   binning=self.binning)
            # Rebuild the mask with this new crmask
            self.mask = self.update_mask_cr(self.mask, self.crmask)

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            sobjs_show = None if show_objs else self.sobjs_obj
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', slits=True, sobjs =sobjs_show, clear=False)

        # Return
        return self.global_sky
    #
    # def get_init_sky(self, tslits_dict, tilts, maskslits=None, update_crmask = True, show_fit = False, show = False):
    #
    #     self.sobjs_obj_init, self.nobj_init, self.skymask = \
    #         self.find_objects(self.sciimg, tslits_dict,maskslits=maskslits)
    #
    #     if self.ir_redux:
    #         self.sobjs_obj_init_neg, self.nobj_init_neg, self.skymask_neg = \
    #             self.find_objects(-self.sciimg, tslits_dict, maskslits=maskslits)
    #         skymask = self.skymask & self.skymask_neg
    #     else:
    #         skymask = self.skymask
    #
    #     # Global sky subtraction, first pass. Uses skymask from object finding step above
    #     self.global_sky = self.global_skysub(tslits_dict, tilts, skymask=skymask, update_crmask = update_crmask,
    #                                          maskslits=maskslits, show_fit = show_fit, show=show)
    #
    #     return self.global_sky

    # JFH TODO Should we reduce the number of iterations for standards or near-IR redux where the noise model is not
    # being updated?
    def local_skysub_extract(self, sobjs, waveimg, maskslits=None, model_noise=True, std = False,
                             show_profile=False, show_resids=False, show=False):
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

        if not self._chk_objs([ # Did they run process?
                                'sciimg', 'sciivar', 'rn2img',
                                # Did they run global sky subtraction, self.global_skysub()?
                                'global_sky',
                                # Did the input the right calibrations in prev steps?
                                'tilts', 'tslits_dict']):
            msgs.error('All quantities necessary to run local_skysub_extract() have not been set.')


        self.waveimg = waveimg
        # get the good slits and assign self.maskslits
        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        self.outmask = np.copy(self.mask)
        # Initialize to input mask in case no objects were found
        #self.extractmask = (self.mask == 0) & self.negmask
        self.extractmask = (self.mask == 0)
        # Initialize to zero in case no objects were found
        self.objmodel = np.zeros_like(self.sciimg)
        # Set initially to global sky in case no objects were found
        self.skymodel  = np.copy(self.global_sky)
        # Set initially to sciivar in case no obects were found.
        self.ivarmodel = np.copy(self.sciivar)

        # Could actually create a model anyway here, but probably
        # overkill since nothing is extracted

        self.sobjs = sobjs.copy()
        # Loop on slits
        for slit in gdslits:
            msgs.info("Local sky subtraction and extraction for slit: {:d}".format(slit))
            thisobj = (self.sobjs.slitid == slit) # indices of objects for this slit
            if np.any(thisobj):
                thismask = (self.slitmask == slit) # pixels for this slit
                # True  = Good, False = Bad for inmask
                inmask = (self.mask == 0) & (self.crmask == False) & thismask
                # Local sky subtraction and extraction
                self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], \
                    self.extractmask[thismask] \
                        = skysub.local_skysub_extract(self.sciimg, self.sciivar, self.tilts,
                                                      self.waveimg, self.global_sky, self.rn2img,
                                                      thismask, self.tslits_dict['lcen'][:,slit],
                                                      self.tslits_dict['rcen'][:, slit],
                                                      self.sobjs[thisobj], model_noise=model_noise,
                                                      std = std, bsp=self.par['bspline_spacing'],
                                                      sn_gauss=self.par['sn_gauss'],
                                                      inmask=inmask, show_profile=show_profile,
                                                      show_resids=show_resids)

        # Set the bit for pixels which were masked by the extraction.
        # For extractmask, True = Good, False = Bad
        iextract = (self.mask == 0) & (self.extractmask == False)
        self.outmask[iextract] += np.uint64(2**8)
        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            self.show('local', sobjs = self.sobjs, slits= True)
            self.show('resid', sobjs = self.sobjs, slits= True)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs


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
        elif (self.maskslits is None):
            # If maskslits was not passed, and it does not exist in self, reduce all slits
            self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
        else: # Otherwise, if self.maskslits exists, use the previously set maskslits
            pass
        return self.maskslits
    #
    #
    # # JFH TODO This stuff should be eventually moved to proc
    # def proc_old(self, bias, pixel_flat, bpm, illum_flat=None, sigma_clip=False, sigrej=None, maxiters=5, show=False):
    #     """ Process the image
    #
    #     Wrapper to ProcessImages.process()
    #
    #     Needed in part to set self.sciframe, although I could kludge it another way..
    #
    #     Returns
    #     -------
    #     self.sciframe
    #     self.rawvarframe
    #     self.crmask
    #
    #     """
    #     # Process
    #     self.bpm = bpm
    #     self.bias = bias
    #     self.pixel_flat = pixel_flat
    #     self.illum_flat = illum_flat
    #
    #     if self.ir_redux:
    #         if sigma_clip is True:
    #             msgs.error('You cannot sigma clip with difference imaging as this will reject objects')
    #         all_files = self.file_list + self.bg_file_list
    #         cosmics = False # If we are differencing CR reject after we difference for better performance
    #         # weights account for possibility of differing number of sci and bg images, i.e.
    #         #  stack = 1/n_sci \Sum sci  - 1/n_bg \Sum bg
    #         weights = np.hstack((np.ones(self.nsci)/float(self.nsci),-1.0*np.ones(self.nbg)/float(self.nbg)))
    #     else:
    #         all_files = self.file_list
    #         cosmics = True
    #         weights = np.ones(self.nsci)/float(self.nsci)
    #
    #     sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack = \
    #         self.read_stack(all_files, bias, pixel_flat, bpm, illum_flat, cosmics=cosmics)
    #     nfiles = len(all_files)
    #
    #     # ToDO The bitmask is not being properly propagated here!
    #     if self.nsci > 1 or self.ir_redux:
    #         if sigma_clip:
    #             msgs.error('Sigma clipping is not yet supported')
    #             if sigrej is None:
    #             if self.nsci <= 2:
    #                 sigrej = 100.0  # Irrelevant for only 1 or 2 files, we don't sigma clip below
    #             elif self.nsci == 3:
    #                 sigrej = 1.1
    #             elif self.nsci == 4:
    #                 sigrej = 1.3
    #             elif self.nsci == 5:
    #                 sigrej = 1.6
    #             elif self.nsci == 6:
    #                 sigrej = 1.9
    #             else:
    #                 sigrej = 2.0
    #             # sigma clip if we have enough images
    #             if self.nsci > 2: # cannot sigma clipo for <= 2 images
    #                 ## TODO THis is not tested!!
    #                 # JFH ToDO Should we be sigma clipping here at all? What if the two background frames are not
    #                 # at the same location, this then causes problems?
    #                 # mask_stack > 0 is a masked value. numpy masked arrays are True for masked (bad) values
    #                 data = np.ma.MaskedArray(sciimg_stack, (mask_stack > 0))
    #                 sigclip = stats.SigmaClip(sigma=sigrej, maxiters=maxiters,cenfunc='median')
    #                 data_clipped = sigclip(data, axis=0, masked=True)
    #                 outmask_stack = np.invert(data_clipped.mask) # outmask = True are good values
    #         else:
    #             outmask_stack = (mask_stack == 0)  # outmask = True are good values
    #
    #             var_stack = utils.calc_ivar(sciivar_stack)
    #             weights_stack = np.einsum('i,ijk->ijk',weights,outmask_stack)
    #             weights_sum = np.sum(weights_stack, axis=0)
    #             # Masked everwhere nused == 0
    #             self.crmask = np.sum(crmask_stack,axis=0) == nfiles # Was everywhere a CR
    #             self.sciimg = np.sum(sciimg_stack*weights_stack,axis=0)/(weights_sum + (weights_sum == 0.0))
    #             varfinal = np.sum(var_stack*weights_stack**2,axis=0)/(weights_sum + (weights_sum == 0.0))**2
    #             self.sciivar = utils.calc_ivar(varfinal)
    #             self.rn2img = np.sum(rn2img_stack*weights_stack**2,axis=0)/(weights_sum + (weights_sum == 0.0))**2
    #             # ToDO If I new how to add the bits, this is what I would do do create the mask. For now
    #             # we simply create the mask again using the stacked images and the stacked mask
    #             #nused = np.sum(outmask_stack,axis=0)
    #             #self.mask = (nused == 0) * np.sum(mask_stack, axis=0)
    #             self.mask = self._build_mask(self.sciimg, self.sciivar, self.crmask, saturation=self.saturation,
    #                                          mincounts=(not self.ir_redux))
    #     else:
    #         self.mask  = mask_stack[0,:,:]
    #         self.crmask = crmask_stack[0,:,:]
    #         self.sciimg = sciimg_stack[0,:,:]
    #         self.sciivar = sciivar_stack[0,:,:]
    #         self.rn2img = rn2img_stack[0,:,:]
    #
    #     # For an IR reduction we build the CR mask after differencing
    #     if self.ir_redux:
    #         self.crmask = self.build_crmask(self.sciimg, ivar=self.sciivar)
    #         self.mask = self._update_mask_cr(self.mask, self.crmask)
    #
    #     # Toggle the OFFSLIT bit at the very end
    #     self.mask = self._update_mask_slitmask(self.mask, self.slitmask)
    #
    #     # Show the science image if an interactive run, only show the crmask
    #     if show:
    #         # Only mask the CRs in this image
    #         self.show('image', image=self.sciimg*(self.crmask == 0), chname='sciimg')
    #
    #     return self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask

    def proc(self, bias, pixel_flat, bpm, illum_flat=None, reject_cr=True, sigma_clip=False, sigrej=None,
             maxiters=5,show=False):

        # Process
        self.bpm = bpm
        self.bias = bias
        self.pixel_flat = pixel_flat
        self.illum_flat = illum_flat

        if self.ir_redux:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_diff(
                self.file_list, self.bg_file_list, reject_cr=True, sigma_clip=False)
        else:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_sci(
                self.file_list, reject_cr=True, sigma_clip=False)

        # Now add the slitmask to the mask (i.e. post CR reject)
        self.mask = self.update_mask_slitmask(self.mask, self.slitmask)

        # Show the science image if an interactive run, only show the crmask
        if show:
            # Only mask the CRs in this image
            self.show('image', image=self.sciimg * (self.crmask == 0), chname='sciimg')

        return self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask

    def proc_sci(self, file_list, reject_cr=True, sigma_clip=False, sigrej=None, maxiters=5):
        """

        Args:
            bias:
            pixel_flat:
            bpm:
            illum_flat:
            sigma_clip:
            sigrej:
            maxiters:
            show:

        Returns:

        """

        nsci = len(file_list)
        weights = np.ones(nsci)/float(nsci)
        sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack = \
        self.read_stack(file_list, self.bias, self.pixel_flat, self.bpm, self.det, self.proc_par, self.spectrograph,
                            illum_flat=self.illum_flat, reject_cr=reject_cr, binning=self.binning)

        # ToDO The bitmask is not being properly propagated here!
        if nsci > 1:
            sciimg, sciivar, rn2img, outmask = procimg.weighted_combine(
                weights, sciimg_stack, sciivar_stack, rn2img_stack, (mask_stack == 0),
                sigma_clip=sigma_clip, sigrej=sigrej, maxiters=maxiters)
            # assumes everything masked in the outmask is a CR in the individual images
            crmask = np.invert(outmask)
            # Create a mask for this image now
            mask = self.build_mask(sciimg, sciivar, crmask, self.bpm, saturation=self.saturation, mincounts=self.mincounts)
        else:
            mask = mask_stack[0, :, :]
            crmask = crmask_stack[0, :, :]
            sciimg = sciimg_stack[0, :, :]
            sciivar = sciivar_stack[0, :, :]
            rn2img = rn2img_stack[0, :, :]

        return sciimg, sciivar, rn2img, mask, crmask


    # JFH TODO This stuff should be eventually moved to proc
    def proc_diff(self, file_list, bg_file_list, reject_cr = True,sigma_clip=False, sigrej=None, maxiters=5):
        """ Process the image

        Wrapper to ProcessImages.process()

        Needed in part to set self.sciframe, although I could kludge it another way..

        Returns
        -------
        self.sciframe
        self.rawvarframe
        self.crmask

        """

        sciimg_sci, sciivar_sci, rn2img_sci, mask_sci, crmask_sci = self.proc_sci(
            file_list, reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)
        sciimg_bg, sciivar_bg, rn2img_bg, mask_bg, crmask_bg = self.proc_sci(
            bg_file_list, reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)

        # Combine the images
        outmask_comb = (mask_sci == 0) & (mask_bg == 0)
        sciimg = sciimg_sci - sciimg_bg
        varcomb = utils.calc_ivar(sciivar_sci) + utils.calc_ivar(sciivar_bg)
        sciivar = utils.calc_ivar(varcomb)*outmask_comb
        rn2img = rn2img_sci + rn2img_bg
        # Now reject CRs again on the differenced image
        crmask_diff = self.build_crmask(sciimg, self.proc_par, self.det, self.spectrograph, ivar=sciivar, binning=self.binning)
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        crmask = crmask_diff | np.invert(outmask_comb)
        # Create a mask for this image now
        mask = self.build_mask(sciimg, sciivar, crmask, self.bpm, saturation=self.saturation)

        return sciimg, sciivar, rn2img, mask, crmask


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
            mask_in = self.mask
            bitmask_in = self.bitmask
        else:
            mask_in = None
            bitmask_in = None

        if attr == 'global':
            # global sky subtraction
            if self.sciimg is not None and self.global_sky is not None and self.mask is not None:
                # sky subtracted image
                image = (self.sciimg - self.global_sky)*(self.mask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.mask == 0], sigma_lower=5.0,
                                                       sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                ch_name = chname if chname is not None else 'global_sky_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                              mask=mask_in, clear=clear, wcs_match=True)
                                              #, cuts=(cut_min, cut_max))
        elif attr == 'local':
            # local sky subtraction
            if self.sciimg is not None and self.skymodel is not None and self.mask is not None:
                # sky subtracted image
                image = (self.sciimg - self.skymodel)*(self.mask == 0)
                mean, med, sigma = stats.sigma_clipped_stats(image[self.mask == 0], sigma_lower=5.0,
                                                       sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                ch_name = chname if chname is not None else 'local_sky_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in,
                                              mask=mask_in, clear=clear, wcs_match=True)
                                              #, cuts=(cut_min, cut_max))
        elif attr == 'sky_resid':
            # sky residual map with object included
            if self.sciimg is not None and self.skymodel is not None \
                    and self.objmodel is not None and self.ivarmodel is not None \
                    and self.mask is not None:
                image = (self.sciimg - self.skymodel) * np.sqrt(self.ivarmodel)
                image *= (self.mask == 0)
                ch_name = chname if chname is not None else 'sky_resid_{}'.format(self.det)
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0),
                                              bitmask=bitmask_in, mask=mask_in, clear=clear,
                                              wcs_match=True)
        elif attr == 'resid':
            # full residual map with object model subtractede
            if self.sciimg is not None and self.skymodel is not None \
                    and self.objmodel is not None and self.ivarmodel is not None \
                    and self.mask is not None:
                # full model residual map
                image = (self.sciimg - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel)
                image *= (self.mask == 0)
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
                ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)

        if slits:
            if self.tslits_dict is not None:
                slit_ids = [trace_slits.get_slitid(self.sciimg.shape, self.tslits_dict['lcen'],
                                                   self.tslits_dict['rcen'], ii)[0]
                                for ii in range(self.tslits_dict['lcen'].shape[1])]

                ginga.show_slits(viewer, ch, self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                 slit_ids)  # , args.det)

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



