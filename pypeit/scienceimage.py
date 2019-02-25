""" Module for the ScienceImage class"""
from __future__ import absolute_import, division, print_function

import numpy as np
from pypeit import msgs
from pypeit import processimages
from pypeit import utils
from pypeit import ginga
from pypeit.core import coadd2d

from pypeit import debugger


class ScienceImage(processimages.ProcessImages):
    """
    This class will organize and run actions related to
    a Science or Standard star exposure

    ..todo.. Clean this up JFH

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
    bpm : ndarray
      Bad pixel mask
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
    frametype = 'science'

    # TODO: Merge into a single parset, one for procing, and one for scienceimage
    def __init__(self, spectrograph, file_list, bg_file_list = [], ir_redux=False, det=1, binning=None, par=None):


        # Setup the parameters sets for this object. NOTE: This uses objtype, not frametype!
        #self.par = pypeitpar.FrameGroupPar(objtype) if par is None else par
        self.par = spectrograph.default_pypeit_par()['scienceframe'] if par is None else par

        # Start up by instantiating the process images class for reading in the relevant science files
        processimages.ProcessImages.__init__(self, spectrograph, self.par['process'],
                                             files=[], det=det)

        # Instantiation attributes for this object
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

        # Set some detector parameters that we will need
        self.saturation = self.spectrograph.detector[self.det - 1]['saturation']
        self.mincounts = self.spectrograph.detector[self.det - 1]['mincounts']

        # These attributes will be sert when the image(s) are processed
        self.bpm = None
        self.bias = None
        self.pixel_flat = None
        self.illum_flat = None

        self.steps = []

        # Other bookeeping internals
        self.crmask = None
        self.mask = None


    # JFH TODO This stuff should be eventually moved to processimages?
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
        self.bpm = bpm
        self.bias = bias
        self.pixel_flat = pixel_flat
        self.illum_flat = illum_flat

        if self.ir_redux:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_diff(
                self.file_list, self.bg_file_list, reject_cr=True, sigma_clip=False,
                sigrej=sigrej, maxiters=maxiters)
        else:
            self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask = self.proc_sci(
                self.file_list, reject_cr=True, sigma_clip=False, sigrej=sigrej, maxiters=maxiters)

        # Show the science image if an interactive run, only show the crmask
        if show:
            # Only mask the CRs in this image
            self.show(self.sciimg * (self.crmask == 0), chname='sciimg')

        return self.sciimg, self.sciivar, self.rn2img, self.mask, self.crmask

    def proc_sci(self, file_list, reject_cr=True, sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images

        This includes stacking the images if there is more than 1

        Args:
            file_list (list): List of filenames for science frames
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
        nsci = len(file_list)
        weights = np.ones(nsci)/float(nsci)
        sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack = \
        self.read_stack(file_list, self.bias, self.pixel_flat, self.bpm, self.det, self.par['process'], self.spectrograph,
                            illum_flat=self.illum_flat, reject_cr=reject_cr, binning=self.binning)

        # ToDO The bitmask is not being properly propagated here!

        if nsci > 1:
            sci_list = [sciimg_stack]
            var_stack = utils.calc_ivar(sciivar_stack)
            var_list = [var_stack, rn2img_stack]
            sci_list_out, var_list_out, outmask, nused = coadd2d.weighted_combine(
                weights, sci_list, var_list, (mask_stack == 0),
                sigma_clip=sigma_clip, sigma_clip_stack = sciimg_stack, sigrej=sigrej, maxiters=maxiters)
            sciimg = sci_list_out[0]
            sciivar = utils.calc_ivar(var_list_out[0])
            rn2img = var_list_out[1]
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


    def proc_diff(self, file_list, bg_file_list, reject_cr = True,sigma_clip=False, sigrej=None, maxiters=5):
        """
        Process a list of science images and their background frames
        Primarily for near-IR reductions

        Wrapper to proc_sci for

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
        crmask_diff = self.build_crmask(sciimg, self.par['process'], self.det, self.spectrograph, ivar=sciivar, binning=self.binning)
        # crmask_eff assumes evertything masked in the outmask_comb is a CR in the individual images
        crmask = crmask_diff | np.invert(outmask_comb)
        # Create a mask for this image now
        mask = self.build_mask(sciimg, sciivar, crmask, self.bpm, saturation=self.saturation)

        return sciimg, sciivar, rn2img, mask, crmask


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
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                        self.nsci)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt



