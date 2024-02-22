""" Module for the SciImgStack class"""
import numpy as np
from pypeit import msgs
from pypeit import utils
from pypeit import ginga
from pypeit.deprecated import coadd2d
from pypeit.par import pypeitpar

from pypeit.images import scienceimage
from pypeit.images import maskimage


class SciImgStack(object):
    """
    This class will organize and run actions related to
    a Science or Standard star exposure

    ..todo.. Clean this up JFH

    Parameters
    ----------
    spectrograph : pypeit.spectrograph.Spectrograph
    file_list : list
      List of raw files to produce the flat field
    par (PypeItPar):
    settings : dict-like
    tslits_dict : dict
      dict from TraceSlits class
    tilts : ndarray
      tilts from WaveTilts class
      used for sky subtraction and object finding
    det : int
    sci_bpm : ndarray
      Bad pixel mask for this science image;  can and often does differ
      from the default BPM
    objtype : str
      'science'
      'standard'

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
    frametype = 'scienceimages'

    def __init__(self, spectrograph, file_list, par, bg_file_list=[], ir_redux=False,
                 det=1, binning=None):


        # Instantiation attributes for this object
        self.spectrograph = spectrograph
        self.par = par
        if not isinstance(self.par, pypeitpar.FrameGroupPar):
            msgs.error("Bad par type")
        if not isinstance(file_list, list):
            msgs.error("Bad file_list type")
        self.file_list = file_list
        if not isinstance(bg_file_list, list):
            msgs.error("Bad file_list type")
        self.bg_file_list = bg_file_list

        # Are we subtracting the sky using background frames? If yes, set ir_redux=True
        self.ir_redux = ir_redux
        if self.ir_redux and (len(self.bg_file_list) == 0):
            msgs.error('IR reductions require that bg files are specified')
        self.det = det
        self.binning = binning

        # Set some detector parameters that we will need
        self.saturation = self.spectrograph.detector[self.det - 1]['saturation']
        self.mincounts = self.spectrograph.detector[self.det - 1]['mincounts']

        # These attributes will be set when the image(s) are processed
        self.bias = None
        self.sci_bpm = None
        self.pixel_flat = None
        self.illum_flat = None

        self.steps = []

        # Other bookeeping internals
        self.crmask = None
        self.mask = None

    @property
    def nfiles(self):
        """
        Returns:
            int: Number of science frames as per file_list

        """
        return len(self.file_list)

    def build_stack(self, files, bpm, reject_cr=False):
        """
        Generate a set of stack arrays useful for image processing

        These are generated from the internal ProcessImage objects
        held in self.pimages which need to have been previously
        generated/loaded

        Note:  To save on memory, the image attributes of each
        ProcessImage in self.pimages is reset after it contributes
        to the stack arrays.

        Args:
            bpm (np.ndarray):
                Bad pixel mask image
            reject_cr (bool, optional):
                If true, generate a CR image

        Returns:
            tuple: np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray
                sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack

        """
        # Get it ready
        nimages = len(files)
        shape = (nimages, bpm.shape[0], bpm.shape[1])
        sciimg_stack = np.zeros(shape)
        sciivar_stack= np.zeros(shape)
        rn2img_stack = np.zeros(shape)
        crmask_stack = np.zeros(shape, dtype=bool)

        # Mask
        bitmask = maskimage.ImageBitMask()
        mask_stack = np.zeros(shape, bitmask.minimum_dtype(asuint=True))

        # Loop on the files
        for kk, ifile in enumerate(files):
            # Instantiate
            sciImage = scienceimage.ScienceImage(self.spectrograph, self.det,
                                                 self.par['process'], bpm)
            # Process
            sciimg_stack[kk,:,:] = sciImage.process_raw(ifile, self.bias, self.pixel_flat, illum_flat=self.illum_flat)
            # Construct raw variance image and turn into inverse variance
            sciivar_stack[kk, :, :] = sciImage.build_ivar()
            # Mask cosmic rays
            if reject_cr:
                crmask_stack[kk, :, :] = sciImage.build_crmask()
            # Build read noise squared image
            rn2img_stack[kk, :, :] = sciImage.build_rn2img()
            # Final mask for this image
            mask_stack[kk, :, :] = sciImage.build_mask(
                saturation=self.spectrograph.detector[self.det - 1]['saturation'],
                mincounts=self.spectrograph.detector[self.det - 1]['mincounts'])

        # Return
        return sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack

    def proc(self, bias, pixel_flat, bpm, illum_flat=None, sigrej=None, maxiters=5, show=False):
        """
        Primary wrapper for processing one or more science frames or science frames with bgframes

        Args:
            bias (ndarray, None or str):  Specifies bias subtraction approach and/or provides bias image
            pixel_flat (ndarray):  Pixel flat image
            bpm (ndarray):  Bad pixel mask
            illum_flat (ndarray, optional): Illumination flat
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the numberr of images provided.
            maxiters (int, optional):
            show (bool, optional):

        Returns:
            ndarray, ndarray, ndarray, ndarray, ndarray:
              sciimg
              sciivar
              rn2img
              mask
              crmask

        """

        # Process
        self.sci_bpm = bpm
        self.bias = bias
        self.pixel_flat = pixel_flat
        self.illum_flat = illum_flat

        if self.ir_redux:
            sciImg = self.proc_diff(reject_cr=True, sigma_clip=False, sigrej=sigrej, maxiters=maxiters)
        else:
            sciImg = self.proc_list('sci', reject_cr=True, sigma_clip=False, sigrej=sigrej, maxiters=maxiters)

        # Show the science image if an interactive run, only show the crmask
        if show:
            # Only mask the CRs in this image
            self.show(sciImg.image * (sciImg.crmask == 0), chname='sciimg')

        #return self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask
        return sciImg

    def proc_list(self, ltype, reject_cr=True, sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images

        This includes stacking the images if there is more than 1

        Args:
            ltype (str): Type of images to process ('sci', 'bkg')
            reject_cr (bool, optional):
            sigrej (int or float, optional): Rejection threshold for sigma clipping.
                 Code defaults to determining this automatically based on the numberr of images provided.
            maxiters (int, optional):
            show (bool, optional):

        Returns:
            ndarray, ndarray, ndarray, ndarray, ndarray:
              sciimg
              sciivar
              rn2img
              mask
              crmask

        """
        # Init
        if ltype == 'sci':
            files = self.file_list
        elif ltype == 'bkg':
            files = self.bg_file_list
        else:
            msgs.error("Bad ltype for proc_list")
        nimg = len(files)
        weights = np.ones(nimg)/float(nimg)

        # Load
        img_stack, ivar_stack, rn2img_stack, crmask_stack, mask_stack = self.build_stack(
            files, self.sci_bpm, reject_cr=reject_cr)

        # ToDO The bitmask is not being properly propagated here!
        if nimg > 1:
            img_list = [img_stack]
            var_stack = utils.inverse(ivar_stack, positive=True)
            var_list = [var_stack, rn2img_stack]
            img_list_out, var_list_out, outmask, nused = coadd2d.weighted_combine(
                weights, img_list, var_list, (mask_stack == 0),
                sigma_clip=sigma_clip, sigma_clip_stack = img_stack, sigrej=sigrej, maxiters=maxiters)
            '''
            img = img_list_out[0]
            ivar = utils.calc_ivar(var_list_out[0])
            rn2img = var_list_out[1]
            '''
            sciImage = scienceimage.ScienceImage.from_images(self.spectrograph, self.det,
                                                             self.par['process'], self.sci_bpm,
                                                             img_list_out[0],
                                                             utils.inverse(var_list_out[0], positive=True),
                                                             var_list_out[1], np.invert(outmask),
                                                             files=files)
            sciImage.build_mask(saturation=self.saturation, mincounts=self.mincounts)
            '''
            # assumes everything masked in the outmask is a CR in the individual images
            crmask = np.invert(outmask)
            # Create a mask for this combined image
            #processImage = processimage.ProcessImage(None, self.spectrograph, self.det, self.proc_par)
            #processImage.image = img
            #processImage.rawvarframe = var_list_out[0]
            #processImage.crmask = crmask
            #mask = procimg.build_mask(bpm=self.sci_bpm, saturation=self.saturation, mincounts=self.mincounts)
            mask = procimg.build_mask(self.bitmask, img, ivar, self.sci_bpm, crmask,
                                           saturation=self.saturation, mincounts=self.mincounts)
            '''
        else:
            mask = mask_stack[0, :, :]
            crmask = crmask_stack[0, :, :]
            img = img_stack[0, :, :]
            ivar = ivar_stack[0, :, :]
            rn2img = rn2img_stack[0, :, :]
            sciImage = scienceimage.ScienceImage.from_images(self.spectrograph, self.det,
                                                             self.par['process'], self.sci_bpm,
                                                             img, ivar, rn2img, crmask=crmask,
                                                             mask=mask, files=files)

        return sciImage


    def proc_diff(self, reject_cr=True, sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images and their background frames
        Primarily for near-IR reductions

        Wrapper to proc_sci for

        Needed in part to set self.sciframe, although I could kludge it another way..

        Args:
            file_list:
            bg_file_list:
            reject_cr:
            sigma_clip:
            sigrej:
            maxiters:

        Returns:
            tuple: sciimg, sciivar, rn2img, mask, crmask

        """

        sciImg = self.proc_list('sci', reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)
        bgImg = self.proc_list('bkg', reject_cr=reject_cr, sigma_clip=sigma_clip, sigrej=sigrej,maxiters=maxiters)

        new_sciImg = sciImg - bgImg

        '''
        # Combine the images
        outmask_comb = (mask_sci == 0) & (mask_bg == 0)
        sciimg = sciimg_sci - sciimg_bg
        varcomb = utils.calc_ivar(sciivar_sci) + utils.calc_ivar(sciivar_bg)
        sciivar = utils.calc_ivar(varcomb)*outmask_comb
        rn2img = rn2img_sci + rn2img_bg
        # Let's do some more processing
        processImage = processimage.ProcessImage(None, self.spectrograph, self.det,
                                                 self.par['process'])
        processImage.image = sciimg
        processImage.rawvarframe = varcomb
        # Now reject CRs again on the differenced image
        crmask_diff = processImage.build_crmask()
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        crmask = crmask_diff | np.invert(outmask_comb)
        # Create a mask for this image now
        mask = processImage.build_mask(bpm=self.sci_bpm, saturation=self.saturation)#, mincounts=self.mincounts)
        #mask = self.build_mask(sciimg, sciivar, crmask, self.sci_bpm, saturation=self.saturation)
        '''
        return new_sciImg  # sciimg, sciivar, rn2img, mask, crmask


    def show(self, image, chname=None):
        """
        Show one of the internal images

        Args:
            image : ndarray, optional
              User supplied image to display

        """

        ch_name = chname if chname is not None else 'image'
        viewer, ch = ginga.show_image(image, chname=ch_name)


    def __repr__(self):
        txt = '<{:s}: nfiles={:d}'.format(self.__class__.__name__,
                                        self.nfiles)
        txt += '>'
        return txt



