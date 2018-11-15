# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

import time
import datetime

from multiprocessing import Process

from astropy.stats import sigma_clipped_stats

from pypeit import msgs
from pypeit import processimages
from pypeit import specobjs
from pypeit import utils
from pypeit import ginga
from pypeit.core import procimg
from pypeit.core import skysub
from pypeit.core import extract
from pypeit.core import trace_slits
from pypeit.par import pypeitpar

from pypeit.bitmask import BitMask

from pypeit import debugger

class ScienceImageBitMask(BitMask):
    """
    Define a bitmask used to set the reasons why each pixel in a science
    image was masked.
    """
    def __init__(self):
        # TODO:
        #   - Can IVAR0 and IVAR_NAN be consolidated into a single bit?
        #   - Is EXTRACT ever set?
        mask = {       'BPM': 'Component of the instrument-specific bad pixel mask',
                        'CR': 'Cosmic ray detected',
                'SATURATION': 'Saturated pixel',
                 'MINCOUNTS': 'Pixel below the instrument-specific minimum counts',
                  'OFFSLITS': 'Pixel does not belong to any slit',
                    'IS_NAN': 'Pixel value is undefined',
                     'IVAR0': 'Inverse variance is undefined',
                  'IVAR_NAN': 'Inverse variance is NaN',
                   'EXTRACT': 'Pixel included in spectral extraction'
               }
        super(ScienceImageBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))

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
      dict from TraceSlits class (e.g. slitpix)
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
    def __init__(self, spectrograph, file_list, det=None, objtype='science', scidx=0, setup=None,
                 par=None, frame_par=None):

        # Instantiation attributes for this object
        self.spectrograph = spectrograph
        self.file_list = file_list
        self.det = det
        self.objtype = objtype
        self.scidx = scidx
        self.setup = setup

        # Setup the parameters sets for this object
        # NOTE: This uses objtype, not frametype!
        self.par = pypeitpar.ScienceImagePar() if par is None else par
        self.frame_par = pypeitpar.FrameGroupPar(objtype) if frame_par is None else frame_par

        # These attributes will be sert when the image(s) are processed
        self.bpm = None
        self.bias = None
        self.pixflat = None

        # Start up by instantiating the process images class for reading
        # in the relevant science files
        processimages.ProcessImages.__init__(self, spectrograph, file_list, det=det,
                                             par=self.frame_par['process'])

        # Set atrributes for this file and detector using spectrograph class
        self.datasec_img = spectrograph.get_datasec_img(file_list[0], det = det)



        # Other attributes that will be set later during object finding,
        # sky-subtraction, and extraction
        self.tslits_dict = None # used by find_object
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
        self.objmask = None
        self.outmask = None
        self.mask = None                        # The composite bit value array
        self.bitmask = ScienceImageBitMask()    # The bit mask interpreter
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None # Only object finding but no extraction
        self.sobjs = None  # Final extracted object list with trace corrections applied
        self.qa_proc_list = []

        # Other bookeeping internals
        self.exptime = None
        self.binning = None
        self.time = None
        self.inst_name = None
        self.target_name = None
        self.basename = None

        # Child-specific Internals
        #    See ProcessImages
        self.crmask = None

    def _chk_objs(self, items):
        """

        Args:
            items:

        Returns:

        """
        for obj in items:
            if getattr(self, obj) is None:
                msgs.warn('You need to generate {:s} prior to this calibration..'.format(obj))
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

    def find_objects(self, tslits_dict, maskslits=None, skysub=True, show_peaks=False,
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
        self.specobjs : Specobjs object
                Container holding Specobj objects
        self.skymask : ndarray
                Boolean image indicating which pixels are useful for global sky subtraction
        self.objmask : ndarray
                Boolean image indicating which pixels have object flux on them

        """

        self.tslits_dict = tslits_dict
        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        # create the ouptut images skymask and objmask
        self.skymask = np.zeros_like(self.sciimg,dtype=bool)
        self.objmask = np.zeros_like(self.sciimg,dtype=bool)

        # If we are object finding on the sky subtracted image, then
        # check that the global sky exists
        if skysub is True:
            if self.global_sky is None:
                msgs.error('Object finding on sky subtracted image requested, but global_sky '
                           'is not set. Run global_skysub() first')
            image = self.sciimg - self.global_sky
        else:
            image = self.sciimg

        # Build and assign the input mask
        self.mask = self._build_mask()
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit in gdslits:
            qa_title ="Finding objects on slit # {:d}".format(slit +1)
            msgs.info(qa_title)
            thismask = (self.tslits_dict['slitpix'] == slit + 1)
            inmask = (self.mask == 0) & thismask
            # Find objects
            specobj_dict = {'setup': self.setup, 'slitid': slit+1, 'scidx': self.scidx,
                            'det': self.det, 'objtype': self.objtype}

            # TODO we need to add QA paths and QA hooks. QA should be
            # done through objfind where all the relevant information
            # is. This will be a png file(s) per slit.
            sobjs_slit, self.skymask[thismask], self.objmask[thismask], proc_list \
                    = extract.objfind(image, thismask, self.tslits_dict['lcen'][:,slit],
                                      self.tslits_dict['rcen'][:,slit], inmask=inmask,
                                      hand_extract_dict=self.par['manual'],
                                      specobj_dict=specobj_dict, show_peaks=show_peaks,
                                      show_fits=show_fits, show_trace=show_trace,
                                      qa_title=qa_title)
            sobjs.add_sobj(sobjs_slit)
            self.qa_proc_list += proc_list

        self.sobjs_obj = sobjs
        # Finish
        self.nobj = len(sobjs)

        # Steps
        self.steps.append(inspect.stack()[0][3])
        if show:
            self.show('image', image=image*(self.mask == 0), chname = 'objfind',
                      sobjs=self.sobjs_obj, slits=True)

        # Return
        return self.sobjs_obj, self.nobj

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

    def global_skysub(self, tslits_dict, tilts, use_skymask=True, maskslits=None, show_fit=False,
                      show=False):
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
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        # Prep
        self.global_sky = np.zeros_like(self.sciimg)

        # Mask objects using the skymask? If skymask has been set by
        # objfinding, and masking is requested, then do so
        skymask = self.skymask if ((self.skymask is not None) & use_skymask) \
                        else np.ones_like(self.sciimg, dtype=bool)

        # Build and assign the input mask
        self.mask = self._build_mask()
        # Loop on slits
        for slit in gdslits:
            msgs.info("Global sky subtraction for slit: {:d}".format(slit +1))
            thismask = (self.tslits_dict['slitpix'] == slit + 1)
            inmask = (self.mask == 0) & thismask & skymask
            # Find sky
            self.global_sky[thismask] =  skysub.global_skysub(self.sciimg, self.sciivar,
                                                              self.tilts, thismask,
                                                              self.tslits_dict['lcen'][:,slit],
                                                              self.tslits_dict['rcen'][:,slit],
                                                              inmask=inmask,
                                                              bsp=self.par['bspline_spacing'],
                                                              show_fit=show_fit)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.maskslits[slit] = True

        # Step
        self.steps.append(inspect.stack()[0][3])

        if show:
            # Global skysub is the first step in a new extraction so clear the channels here
            self.show('global', slits=True, sobjs =self.sobjs_obj, clear=False)


        # Return
        return self.global_sky

    def local_skysub_extract(self, waveimg, maskslits=None, show_profile=False, show_resids=False,
                             show=False):
        """
        Perform local sky subtraction, profile fitting, and optimal extraction slit by slit

        Wrapper to skysub.local_skysub_extract

        Parameters
        ----------

        Optional Parameters
        -------------------
        bspline_spaceing: (float):
           Break-point spacing for bspline

        Returns:
            global_sky: (numpy.ndarray) image of the the global sky model
        """

        self.waveimg = waveimg
        # get the good slits and assign self.maskslits
        self.maskslits = self._get_goodslits(maskslits)
        gdslits = np.where(~self.maskslits)[0]

        if not self._chk_objs([ # Did they run process?
                                'sciimg', 'sciivar', 'rn2img',
                                # Did they run object finding, self.find_objects() ?
                                'sobjs_obj',
                                # Did they run global sky subtraction, self.global_skysub()?
                                'global_sky',
                                # Did the input the right calibrations in prev steps?
                                'tilts', 'waveimg', 'tslits_dict']):
            msgs.error('All quantities necessary to run local_skysub_extract() have not been set.')

        # Build and assign the input mask
        self.mask = self._build_mask()

        # Allocate the images that are needed
        # Initialize to mask in case no objects were found
        self.outmask = np.copy(self.mask)
        # Initialize to input mask in case no objects were found
        self.extractmask = (self.mask == 0)
        # Initialize to zero in case no objects were found
        self.objmodel = np.zeros_like(self.sciimg)
        # Set initially to global sky in case no objects were found
        self.skymodel  = np.copy(self.global_sky)
        # Set initially to sciivar in case no obects were found.
        self.ivarmodel = np.copy(self.sciivar)

        # Could actually create a model anyway here, but probably
        # overkill since nothing is extracted

        self.sobjs = self.sobjs_obj.copy()
        # Loop on slits
        for slit in gdslits:
            msgs.info("Local sky subtraction and extraction for slit: {:d}".format(slit+1))
            thisobj = (self.sobjs.slitid == slit + 1) # indices of objects for this slit
            if np.any(thisobj):
                thismask = (self.tslits_dict['slitpix'] == slit + 1) # pixels for this slit
                # True  = Good, False = Bad for inmask
                inmask = (self.mask == 0) & thismask
                # Local sky subtraction and extraction
                self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], \
                    self.extractmask[thismask] \
                        = skysub.local_skysub_extract(self.sciimg, self.sciivar, self.tilts,
                                                      self.waveimg, self.global_sky, self.rn2img,
                                                      thismask, self.tslits_dict['lcen'][:,slit],
                                                      self.tslits_dict['rcen'][:, slit],
                                                      self.sobjs[thisobj],
                                                      bsp=self.par['bspline_spacing'],
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

        # Clean up any interactive windows that are still up
        for proc in self.qa_proc_list:
            proc.terminate()
            proc.join()

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    def process(self, bias_subtract, pixel_flat, bpm, illum_flat=None, apply_gain=True, trim=True,show=False):
        """ Process the image

        Wrapper to ProcessImages.process()

        Needed in part to set self.sciframe, although I could kludge it another way..

        Returns
        -------
        self.sciframe
        self.rawvarframe
        self.crmask

        """
        # Process
        self.bpm = bpm
        self.bias = bias_subtract
        self.pixflat = pixel_flat

        self.sciimg = super(ScienceImage, self).process(bias_subtract=bias_subtract,
                                                        apply_gain=apply_gain,
                                                        pixel_flat=pixel_flat,
                                                        illum_flat=illum_flat, bpm=self.bpm,
                                                        trim=trim)

        # Construct raw variance image
        rawvarframe = self.build_rawvarframe(trim=trim)
        self.sciivar = utils.calc_ivar(rawvarframe)
        # Build read noise squared image
        self.rn2img = self.build_rn2img()
        # Build CR mask
        self.crmask = self.build_crmask()

        # Show the science image if an interactive run, only show the crmask
        if show:
            # Only mask the CRs in this image
            self.show('image', image=self.sciimg*(self.crmask == 0), chname='sciimg', clear=True)
        return self.sciimg, self.sciivar, self.rn2img, self.crmask

    def _build_mask(self):
        """
        Return the bit value mask used during extraction.
        
        The mask keys are defined by :class:`ScienceImageBitMask`.  Any
        pixel with mask == 0 is valid, otherwise the pixel has been
        masked.  To determine why a given pixel has been masked::

            bitmask = ScienceImageBitMask()
            reasons = bm.flagged_bits(mask[i,j])

        To get all the pixel masked for a specific set of reasons::

            indx = bm.flagged(mask, flag=['CR', 'SATURATION'])

        Returns:
            numpy.ndarray: The bit value mask for the science image.
        """
        # Instatiate the mask
        mask = np.zeros_like(self.sciimg, dtype=self.bitmask.minimum_dtype(asuint=True))

        # Bad pixel mask
        indx = self.bpm.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'BPM')

        # Cosmic rays
        indx = self.crmask.astype(bool)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'CR')

        # Saturated pixels
        indx = self.sciimg >= self.spectrograph.detector[self.det - 1]['saturation']
        mask[indx] = self.bitmask.turn_on(mask[indx], 'SATURATION')

        # Minimum counts
        indx = self.sciimg <= self.spectrograph.detector[self.det - 1]['mincounts']
        mask[indx] = self.bitmask.turn_on(mask[indx], 'MINCOUNTS')

        # Pixels excluded from any slit.  Use a try/except block so that
        # the mask can still be created even if tslits_dict has not
        # been instantiated yet
        # TODO: Is this still necessary?
        try:
            indx = self.tslits_dict['slitpix'] == 0
            mask[indx] = self.bitmask.turn_on(mask[indx], 'OFFSLITS')
        except:
            pass

        # Undefined counts
        indx = np.invert(np.isfinite(self.sciimg))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IS_NAN')

        # Bad inverse variance values
        indx = np.invert(self.sciivar > 0.0)
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR0')

        # Undefined inverse variances
        indx = np.invert(np.isfinite(self.sciivar))
        mask[indx] = self.bitmask.turn_on(mask[indx], 'IVAR_NAN')

        return mask

    def run_the_steps(self):
        """
        Run full the full recipe of calibration steps

        Returns:

        """
        for step in self.steps:
            getattr(self, 'get_{:s}'.format(step))()

    def show(self, attr, image=None, showmask=False, sobjs=None, chname=None, slits=False,
             clear=False):
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
                mean, med, sigma = sigma_clipped_stats(image[self.mask == 0], sigma_lower=5.0,
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
                mean, med, sigma = sigma_clipped_stats(image[self.mask == 0], sigma_lower=5.0,
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
                                        self.nfiles)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt



