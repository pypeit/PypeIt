# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

#from importlib import reload
import datetime

from astropy.time import Time

from pypeit import msgs
from pypeit import processimages
from pypeit import specobjs
from pypeit.core import procimg
from pypeit.core import skysub
from pypeit.core import extract
from pypeit.core import trace_slits
from pypeit import utils
from pypeit import artrace
from pypeit import ginga
from astropy.stats import sigma_clipped_stats
import time


from pypeit.par import pypeitpar

from pypeit import debugger

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
    fitstbl : Table
      Header info
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


    """

    # Frametype is a class attribute
    frametype = 'science'

    # TODO: Merge into a single parset, one for procing, and one for scienceimage
    def __init__(self, spectrograph, file_list, det = None, objtype='science', scidx =0, setup = None,par = None, frame_par = None):

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

        # Start up by instantiating the process images class for reading in the relevant science files
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.frame_par['process'])

        # Set atrributes for this file and detector using spectrograph class
        self.datasec_img = spectrograph.get_datasec_img(file_list[0], det = det)



        # Other attributes that will be set later during object finding, sky-subtraction, and extraction
        self.fitstbl = None  # Set and used in init_time_names below
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
        self.bitmask = None
        self.extractmask = None
        # SpecObjs object
        self.sobjs_obj = None # Only object finding but no extraction
        self.sobjs = None  # Final extracted object list with trace corrections applied


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

    def init_time_names(self, fitstbl):
        """
        Setup the basename (for file output mainly)
        and time objects (for heliocentric)

        Parameters
        ----------
        camera : str
          Taken from settings['mosaic']['camera']
        timeunit : str
          mjd

        Returns
        -------
        self.time : Time
        self.basename : str

        """

        timeunit = self.spectrograph.timeunit
        camera = self.spectrograph.camera

        self.fitstbl =fitstbl
        # ToDo: Given that we just read the file header to get the datasec_img in the init function above, I don't see why
        # I need to access the fits table for exptime and binning. This information is also in the headers. By simply
        # pulling the stuff from the header, we would remove the fitstbl entirely. Another option would be to put the
        # datasec_img stuff in the fitstbl for each detector
        self.exptime = self.fitstbl['exptime'][self.scidx]
        self.binning = self.fitstbl['binning'][self.scidx]

        tbname = None
        try:
            if 'T' in self.fitstbl['date'][self.scidx]:
                tbname = self.fitstbl['date'][self.scidx]
        except IndexError:
            debugger.set_trace()
        else:
            if tbname is None:
                if timeunit == 'mjd':
                    # Not ideal, but convert MJD into a date+time
                    timval = Time(self.fitstbl['time'][self.scidx] / 24.0, scale='tt',
                                  format='mjd')
                    tbname = timval.isot
                else:
                    # Really not ideal... just append date and time
                    tbname = self.fitstbl['date'][self.scidx] + 'T' \
                                    + str(self.fitstbl['time'][self.scidx])
        # Time
        tval = Time(tbname, format='isot')#'%Y-%m-%dT%H:%M:%S.%f')
        dtime = datetime.datetime.strptime(tval.value, '%Y-%m-%dT%H:%M:%S.%f')
        self.time = tval
        # Basename
        self.inst_name = camera
        self.target_name = self.fitstbl['target'][self.scidx].replace(" ", "")
        self.basename = self.target_name+'_'+self.inst_name+'_'+ \
                         datetime.datetime.strftime(dtime, '%Y%b%dT') + \
                         tbname.split("T")[1].replace(':','')
        # Return
        return self.time, self.basename


    def find_objects(self, tslits_dict, maskslits = None, SKYSUB = True, SHOW_PEAKS= False, SHOW_FITS = False,
                     SHOW_TRACE=False, SHOW = False):
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

        # If we are object finding on the sky subtracted image, then check that the global sky exists
        if SKYSUB is True:
            if self.global_sky is None:
                msgs.error('Object finding on sky subtracted image requested, but global_sky is not set. Run global_skysub() first')
            image = self.sciimg - self.global_sky
        else:
            image = self.sciimg

        # Build and assign the input mask
        self.bitmask = self._build_bitmask()
        # Instantiate the specobjs container
        sobjs = specobjs.SpecObjs()

        # Loop on slits
        for slit in gdslits:
            msgs.info("Finding objects on slit: {:d}".format(slit))
            thismask = (self.tslits_dict['slitpix'] == slit + 1)
            inmask = (self.bitmask == 0) & thismask
            # Find objects
            specobj_dict = {'setup': self.setup, 'slitid': slit+1, 'scidx': self.scidx, 'det': self.det, 'objtype': self.objtype}
            # TODO we need to add QA paths and QA hooks. QA should be done through objfind where all the relevant information is. This will
            # be a png file(s) per slit.
            sobjs_slit, self.skymask[thismask], self.objmask[thismask] = extract.objfind(image, thismask,
                                                                                    self.tslits_dict['lcen'][:,slit],
                                                                                    self.tslits_dict['rcen'][:,slit],
                                                                                    inmask = inmask,
                                                                                    HAND_EXTRACT_DICT=self.par['manual'],
                                                                                    specobj_dict=specobj_dict,SHOW_PEAKS=SHOW_PEAKS,
                                                                                    SHOW_FITS = SHOW_FITS,SHOW_TRACE=SHOW_TRACE)
            sobjs.add_sobj(sobjs_slit)

        self.sobjs_obj = sobjs
        # Finish
        self.nobj = len(sobjs)

        # Steps
        self.steps.append(inspect.stack()[0][3])

        if SHOW:
            self.show('image',image = image*(self.bitmask == 0), chname = 'objfind', sobjs =self.sobjs_obj, slits = True)

        # Return
        return self.sobjs_obj, self.nobj

    def global_skysub(self, tslits_dict, tilts, USE_SKYMASK=True, maskslits = None, SHOW_FIT = False,
                      SHOW = False):
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

        # Mask objects using the skymask? If skymask has been set by objfinding, and masking is requested, then do so
        skymask = self.skymask if ((self.skymask is not None) & USE_SKYMASK) else np.ones_like(self.sciimg,dtype=bool)

        # Build and assign the input mask
        self.bitmask = self._build_bitmask()
        # Loop on slits
        for slit in gdslits:
            msgs.info("Working on slit: {:d}".format(slit +1))
            thismask = (self.tslits_dict['slitpix'] == slit + 1)
            inmask = (self.bitmask == 0) & thismask & skymask
            # Find sky
            self.global_sky[thismask] =  skysub.global_skysub(self.sciimg, self.sciivar, self.tilts, thismask,
                                                              self.tslits_dict['lcen'][:, slit], self.tslits_dict['rcen'][:,slit],
                                                              inmask=inmask, bsp=self.par['bspline_spacing'], SHOW_FIT=SHOW_FIT)
            # Mask if something went wrong
            if np.sum(self.global_sky[thismask]) == 0.:
                self.maskslits[slit] = True

        # Step
        self.steps.append(inspect.stack()[0][3])

        if SHOW:
            self.show('global', slits=True)

        # Return
        return self.global_sky

    def local_skysub_extract(self, waveimg, maskslits=None, SHOW_PROFILE=False, SHOW_RESIDS = False, SHOW = False):
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


        if not self._chk_objs(['sciimg','sciivar','rn2img',     # Did they run process?
                               'sobjs_obj',  #Did they run object finding, self.find_objects() ?
                               'global_sky', # Did they run global sky subtraction, self.global_skysub()?
                               'tilts', 'waveimg', 'tslits_dict']): # Did the input the right calibrations in prev steps?
            msgs.error('You do not have all the quantities set necessary to run local_skysub_extract()')

        # Build and assign the input mask
        self.bitmask = self._build_bitmask()

        # Allocate the images that are needed
        self.outmask = np.copy(self.bitmask) # Initialize to bitmask in case no objects were found
        self.extractmask = (self.bitmask == 0)   # Initialize to input mask in case no objects were found
        self.objmodel = np.zeros_like(self.sciimg)      # Initialize to zero in case no objects were found
        self.skymodel  = np.copy(self.global_sky)  # Set initially to global sky in case no objects were found
        self.ivarmodel = np.copy(self.sciivar)          # Set initially to sciivar in case no obects were found.
                                                   # Could actually create a model anyway here,
                                                   # but probalby overkill since nothing is extracted

        self.sobjs = self.sobjs_obj.copy()
        # Loop on slits
        for slit in gdslits:
            msgs.info("Working on slit: {:d}".format(slit))
            thisobj = (self.sobjs.slitid == slit + 1) # indices of objects for this slit
            if np.any(thisobj):
                thismask = (self.tslits_dict['slitpix'] == slit + 1) # pixels for this slit
                # True  = Good, False = Bad for inmask
                inmask = (self.bitmask == 0) & thismask
                # Local sky subtraction and extraction
                self.skymodel[thismask], self.objmodel[thismask], self.ivarmodel[thismask], self.extractmask[thismask] = \
                    skysub.local_skysub_extract(self.sciimg, self.sciivar, self.tilts,self.waveimg, self.global_sky, self.rn2img,
                                                thismask, self.tslits_dict['lcen'][:, slit], self.tslits_dict['rcen'][:, slit],
                                                self.sobjs[thisobj],bsp=self.par['bspline_spacing'], inmask = inmask,
                                                SHOW_PROFILE=SHOW_PROFILE, SHOW_RESIDS = SHOW_RESIDS)

        # Set the bit for pixels which were masked by the extraction
        iextract = (self.bitmask == 0) & (self.extractmask == False) # For extractmask, True = Good, False = Bad
        self.outmask[iextract] += np.uint64(2**8)
        # Step
        self.steps.append(inspect.stack()[0][3])

        if SHOW:
            self.show('local', sobjs = self.sobjs, slits= True)
            self.show('resid', sobjs = self.sobjs, slits= True)

        # Return
        return self.skymodel, self.objmodel, self.ivarmodel, self.outmask, self.sobjs

    # Chk attributest for running local_skysub_extract
    def _chk_objs(self, items):
        for obj in items:
            if getattr(self, obj) is None:
                msgs.warn("You need to generate {:s} prior to this calibration..".format(obj))
                if obj in ['sciimg', 'sciivar', 'rn2_img']:
                    msgs.warn("Run the process() method")
                elif obj in ['sobjs_obj']:
                    msgs.warn("Run the find_objects() method")
                elif obj in['global_sky']:
                    msgs.warn("Run the global_skysub() method")
                elif obj in ['tilts', 'tslits_dict'] :
                    msgs.warn("Calibraitons missing: these were required to run find_objects() and global_skysub()")
                elif obj in ['waveimg']:
                    msgs.warn("Calibraitons missing: waveimg must be input as a parameter. Try running calibrations")
                return False
        return True

    def process(self, bias_subtract, pixel_flat, bpm, apply_gain=True, trim=True):
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
                                                          pixel_flat=pixel_flat, bpm=self.bpm,
                                                          trim=trim)

        # Construct raw variance image
        rawvarframe = self.build_rawvarframe(trim=trim)
        self.sciivar = utils.calc_ivar(rawvarframe)
        # Build read noise squared image
        self.rn2img = self.build_rn2img()
        # Build CR mask
        self.crmask = self.build_crmask()

        return self.sciimg, self.sciivar, self.rn2img, self.crmask


    def _build_bitmask(self):
        """
        Create the input mask for various extraction routines. Here is the bitwise key for interpreting masks

        Bit Key
        ---
        BPM            0
        CR             1
        SATURATION     2
        MINCOUNTS      3
        OFFSLITS        4
        IS_NAN         5
        IVAR0          6
        IVAR_NAN       7
        EXTRACT        8

        bitmask = 0 is good, inmask > 0 has been masked.

        To figure out why it has been masked for example you can type

        bpm = (bitmask & np.uint64(2**0)) > 0
        crmask = (bitmask & np.uint64(2**1)) > 0

        etc.

        Returns
        -------
        bitmask

        """

        # Create and assign the inmask
        SATURATION = self.spectrograph.detector[self.det - 1]['saturation']
        MINCOUNTS    = self.spectrograph.detector[self.det - 1]['mincounts']

        bitmask = np.zeros_like(self.sciimg,dtype=np.uint64)
        bitmask[self.bpm == True] += np.uint64(2**0)
        bitmask[self.crmask == True] += np.uint64(2**1)
        bitmask[self.sciimg >= SATURATION] += np.uint64(2**2)
        bitmask[self.sciimg <= MINCOUNTS] += np.uint64(2**3)
        bitmask[self.tslits_dict['slitpix'] == 0] += np.uint64(2**4)
        bitmask[np.isfinite(self.sciimg) == False] += np.uint64(2**5)
        bitmask[self.sciivar <= 0.0] += np.uint64(2**6)
        bitmask[np.isfinite(self.sciivar) == False] += np.uint64(2**7)

        #inmask = (self.bpm == False) & (self.crmask == False) & \
        #         (self.sciivar > 0.0) & np.isfinite(self.sciivar) & \
        #         np.isfinite(self.sciimg) & (self.sciimg < SATURATION) & (self.sciimg > MINCOUNTS)

        return bitmask

    def _get_goodslits(self, maskslits):
        """
        Return the slits to be reduce by going through the maskslits logic below. If the input maskslits is None it
        uses previously assigned maskslits

        Returns
        -------
        gdslits
            numpy array of slit numbers to be reduced
        """

        # Identify the slits that we want to consider.
        if maskslits is not None:  # If maskslits was passed in use it, and update self
            self.maskslits = maskslits
        elif (self.maskslits is None):  # If maskslits was not passed, and it does not exist in self, reduce all slits
            self.maskslits = np.zeros(self.tslits_dict['lcen'].shape[1], dtype=bool)
        else: # Otherwise, if self.maskslits exists, use the previously set maskslits
            pass
        return self.maskslits


    def show(self, attr, image=None, showmask = False, sobjs = None, chname = None, slits = False):
        """
        Show one of the internal images
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

        if attr == 'global':
            # global sky subtraction
            if self.sciimg is not None and self.global_sky is not None and self.bitmask is not None:
                image = (self.sciimg - self.global_sky)*(self.bitmask == 0)  # sky subtracted image
                (mean, med, sigma) = sigma_clipped_stats(image[self.bitmask == 0], sigma_lower=5.0, sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                if showmask:
                    bitmask_in = self.bitmask
                else:
                    bitmask_in = None
                ch_name = chname if chname is not None else 'global_sky'
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask = bitmask_in)
                #cuts=(cut_min, cut_max)
        elif attr == 'local':
            # local sky subtraction
            if self.sciimg is not None and self.skymodel is not None and self.bitmask is not None:
                image = (self.sciimg - self.skymodel)*(self.bitmask == 0)  # sky subtracted image
                (mean, med, sigma) = sigma_clipped_stats(image[self.bitmask == 0], sigma_lower=5.0, sigma_upper=5.0)
                cut_min = mean - 1.0 * sigma
                cut_max = mean + 4.0 * sigma
                if showmask:
                    bitmask_in = self.bitmask
                else:
                    bitmask_in = None
                ch_name = chname if chname is not None else 'local_sky'
                viewer, ch = ginga.show_image(image, chname=ch_name, bitmask=bitmask_in)
                #cuts=(cut_min, cut_max),
        elif attr == 'sky_resid':
            # sky residual map with object included
            if self.sciimg is not None and self.skymodel is not None and \
                    self.objmodel is not None and self.ivarmodel is not None and self.bitmask is not None:
                image = (self.sciimg - self.skymodel) * np.sqrt(self.ivarmodel) * (self.bitmask == 0)
                if showmask:
                    bitmask_in = self.bitmask
                else:
                    bitmask_in = None
                ch_name = chname if chname is not None else 'sky_resid'
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0), bitmask=bitmask_in)
        elif attr == 'resid':
            # full residual map with object model subtractede
            if self.sciimg is not None and self.skymodel is not None and \
                    self.objmodel is not None and self.ivarmodel is not None and self.bitmask is not None:
                image = (self.sciimg - self.skymodel - self.objmodel) * np.sqrt(self.ivarmodel) * (self.bitmask == 0)  # full model residual map
                if showmask:
                    bitmask_in = self.bitmask
                else:
                    bitmask_in = None
                ch_name = chname if chname is not None else 'resid'
                viewer, ch = ginga.show_image(image, chname=ch_name, cuts=(-5.0, 5.0), bitmask=bitmask_in)
        elif attr == 'image':
            ch_name = chname if chname is not None else 'image'
            viewer, ch = ginga.show_image(image, chname = ch_name)
        else:
            msgs.warn("Not an option for show")

        if sobjs is not None:
            for spec in sobjs:
                if spec.HAND_EXTRACT_FLAG is True:
                    color = 'magenta'
                else:
                    color = 'orange'
                ginga.show_trace(viewer, ch, spec.trace_spat, spec.idx, color=color)

        if slits:
            if self.tslits_dict is not None:
                slit_ids = [trace_slits.get_slitid(self.sciimg.shape, self.tslits_dict['lcen'], self.tslits_dict['rcen'], ii)[0] for ii in
                            range(self.tslits_dict['lcen'].shape[1])]

                ginga.show_slits(viewer, ch,self.tslits_dict['lcen'], self.tslits_dict['rcen'], slit_ids)  # , args.det)

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


    # TODO Everything below here is deprecated
    def _build_modelvar(self, skyframe=None, objframe=None):
        """
        Generate a model variance image using the sky model
        and (optional) object model

        Wrapper to arprocimg.variance_frame

        Parameters
        ----------
        skyframe : ndarray
          Sky model
        objframe : ndarray
          Object model

        Returns
        -------
        self.modelvarframe : ndarray
          Model variance image

        """
        if skyframe is None:
            skyframe = self.global_sky
        self.modelvarframe = procimg.variance_frame(self.datasec_img, self.sciframe,
                                                    self.spectrograph.detector[self.det - 1]['gain'],
                                                    self.spectrograph.detector[self.det - 1]['ronoise'],
                                                    numamplifiers=self.spectrograph.detector[self.det - 1]['numamplifiers'],
                                                    darkcurr=self.spectrograph.detector[self.det - 1]['darkcurr'],
                                                    exptime=self.exptime, skyframe=skyframe, objframe=objframe)
        return self.modelvarframe


    def boxcar(self, mswave):
        """
        Perform boxcar extraction

        Wrapper to arextract.boxcar

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.skycorr_box : ndarray
          Local corrections to the sky model

        Extractions are ingested within self.specobjs

        """
        msgs.info("Performing boxcar extraction")
        self.skycorr_box = extract.boxcar(self.specobjs, self.sciframe,
                                          self.modelvarframe, self.bpm,
                                          self.global_sky, self.crmask,
                                          self.tracelist, mswave,
                                          self.maskslits, self.tslits_dict['slitpix'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.skycorr_box


    def original_optimal(self, mswave):
        """
        Perform optimal extraction using the 'original' PYPIT algorithm

        Wrapper to arextract.obj_profiles and arextract.optimal_extract

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.obj_model : ndarray
          Model of the object flux; used for improving the variance estimate

        """
        msgs.info("Attempting optimal extraction with model profile")
        # Profile
        extract.obj_profiles(self.det, self.specobjs,
                             self.sciframe - self.global_sky - self.skycorr_box,
                             self.modelvarframe, self.crmask, self.tracelist, self.tilts,
                             self.maskslits, self.tslits_dict['slitpix'], doqa=False)
        # Extract
        self.obj_model = extract.optimal_extract(self.specobjs,
                                                 self.sciframe - self.global_sky - self.skycorr_box,
                                                 self.modelvarframe, self.crmask,
                                                 self.tracelist, self.tilts, mswave,
                                                 self.maskslits, self.tslits_dict['slitpix'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.obj_model


    def extraction(self, mswave):
        """
        Perform the extraction
          Boxcar and optimal (as desired)

        Code flow:
          1. Instantiate the specobjs list -- This should be done when finding objects
          2. Boxcar extraction
            i. Update sky model
          3. Optimal extraction
            i. Iterative
            ii. With an update to the variance image
          4. One last update to the variance image

        Parameters
        ----------
        mswave : ndarray
          Wavelength image

        Returns
        -------
        self.specobjs
        self.finalvar
        self.finalsky

        """
        # Init specobjs
        #  Nested -- self.specobjs[slit][object]
        self.specobjs = self._build_specobj()

        # Boxcar -- Fills specobj.boxcar in place
        self.skycorr_box = self.boxcar(mswave)
        self.finalsky = self.global_sky + self.skycorr_box

        if self.objtype != 'standard':
            # Optimal (original recipe)
            self.obj_model = self.original_optimal(mswave)
            #
            msgs.info("Update model variance image (and trace?) and repeat")
            _ = self._build_modelvar(skyframe=self.finalsky, objframe=self.obj_model)
            self.obj_model = self.original_optimal(mswave)

            # Final variance image
            self.finalvar = self._build_modelvar(skyframe=self.finalsky, objframe=self.obj_model)
        else:
            self.finalvar = self.rawvarframe

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.specobjs, self.finalvar, self.finalsky

    def _build_specobj(self):
        """
        Initialize the specobjs for all slits
        Key input is self.tracelist

        Wrapper to arspecobj.init_exp

        Returns
        -------
        self.specobjs : list

        """
        self.specobjs = specobjs.init_exp(self.tslits_dict['lcen'], self.tslits_dict['rcen'],
                                          self.sciframe.shape, self.maskslits, self.det,self.scidx, self.fitstbl, self.tracelist,
                                          binning=self.binning, objtype=self.objtype)
            # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.specobjs



def unpack_bitmask(bitmask):
    """
    Utility function to unpack the bitmask into its respective masks.

    Parameters
    ----------
    bitmask  - ndarray dtype = uint64 created following _build_bitmask() method above

    Returns
    -------
        mask_tuple  = (bpm, crmask, satmask, minmask, offslitmask, nanmask, ivar0mask,ivarnanmask, extractmask)

        tuple of 8 masks corresponding to each bit that can be set in the bitmask


        Bit Key
        ---
        BPM            0
        CR             1
        SATURATION     2
        MINCOUNTS      3
        OFFSLITS        4
        IS_NAN         5
        IVAR0          6
        IVAR_NAN       7
        EXTRACT        8

        bitmask = 0 is good, bitmask > 0 has been masked.

        To figure out why it has been masked for example you can type


        crmask = (inmask & np.uint64(2**1)) > 0

        etc.
    
    """
    bpm         = (bitmask & np.uint64(2 ** 0)) > 0
    crmask      = (bitmask & np.uint64(2 ** 1)) > 0
    satmask     = (bitmask & np.uint64(2 ** 2)) > 0
    minmask     = (bitmask & np.uint64(2 ** 3)) > 0
    offslitmask    = (bitmask & np.uint64(2 ** 4)) > 0
    nanmask     = (bitmask & np.uint64(2 ** 5)) > 0
    ivar0mask   = (bitmask & np.uint64(2 ** 6)) > 0
    ivarnanmask = (bitmask & np.uint64(2 ** 7)) > 0
    extractmask = (bitmask & np.uint64(2 ** 8)) > 0

    return (bpm, crmask, satmask, minmask, offslitmask, nanmask, ivar0mask,ivarnanmask, extractmask)
