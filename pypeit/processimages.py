# Module for Processing Images, e.g. bias frames, arc frames

from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

#from importlib import reload

from astropy.io import fits

from pypeit import msgs
from pypeit import ginga
from pypeit.core import combine
from pypeit.core import procimg
from pypeit.core import flat
from pypeit.core import parse
from pypeit import utils

from pypeit.par import pypeitpar

from pypeit.spectrographs.util import load_spectrograph
from pypeit.bitmask import BitMask


from pypeit import debugger


class ProcessImagesBitMask(BitMask):
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
                   'EXTRACT': 'Pixel masked during local skysub and extraction'
               }
        super(ProcessImagesBitMask, self).__init__(list(mask.keys()), descr=list(mask.values()))


class ProcessImages(object):
    """
    Base class to guide image loading and processing.

    Args:
        file_list (list):
            List of files to read and process.
        spectrograph (:obj:`str`, :class:`Spectrograph`):
            The spectrograph from which the data was taken.  Must be
            provided as a string that can be interpreted by
            :func:`pypeit.spectrographs.util.load_spectrograph` or a
            preconstructed instance of :class:`Spectrograph`.
        det (:obj:`int`, optional):
            The 1-indexed number of the detector.  Default is 1.
        par (:obj:`pypeitpar.ProcessImagesPar`, optional):
            Parameters that dictate the processing of the images.  See
            `pypeit.par.pypeitpar.ProcessImagesPar` for the defaults. 

    Attributes:
        file_list (list):
        det (int):
        spectrograph (:obj:`str`, :class:`Spectrograph`):
        frametype (str):
        stack (:obj:`numpy.ndarray`):
        steps (list):
        raw_images (list):
        headers (list):
        proc_images (:obj:`numpy.ndarray`):
            3D array of processed, individual images
        datasec (list):
            List of **slice** objects that select the data section from
            the images.
        oscansec (list):
            List of **slice** objects that select the overscan section
            from the images.

    Raises:
        IOError: Raised if the provided file_list is not a list.
        TypeError: Raised if the spectrograph is not a :obj:`str` or
            :class:`Spectrograph`.
    """

    # ToDO add a utility function which acts as a wrapper for this class, and simply takes a filename, and returns the
    # image and the classor something, i.e.
    #
    # class = pypeit_proc(filename, biasfile = biasfile, pixflatfile=pixflatfile, illumflatfile=illumflatfile)
    # or maybe it should return the image and the class?

    # Class attribute is unknown.  Will be overwritten by children
    frametype='Unknown'
    bitmask = ProcessImagesBitMask()  # The bit mask interpreter

    def __init__(self, spectrograph, file_list, det=1, par=None):

        # Required parameters
        if not isinstance(file_list, list):
            raise IOError("file_list input to ProcessImages must be list. Empty is fine")
        self.file_list = file_list
        self.spectrograph = load_spectrograph(spectrograph)

        # Optional
        self.det = det

        if par is not None and not isinstance(par, pypeitpar.ProcessImagesPar):
            raise TypeError('Provided ParSet for processing images must be type ProcessImagesPar.')

        # TODO: This can't be called self.par because it may overwrite
        # the self.par of the derived classes (e.g. BiasFrame).  The may
        # be a better way to do this, but I'm punting for now.
        self.proc_par = pypeitpar.ProcessImagesPar() if par is None else par

        # Main (possible) outputs
        self.stack = None
        self.steps = []

        self.raw_images = []
        self.headers = []

        self.proc_images = None  # Will be an ndarray
        self.datasec = []
        self.oscansec = []

        # WARNING: Exposure time None by default here in the base class.
        # The exposure time is currently only defined by ScienceImage
        # and only used by build_rawvarframe
        self.exptime = None

        # Constructed by process:
        self.crmask = None          # build_crmask
        self.rawvarframe = None     # build_rawvarframe
        self.rn2img = None          # build_rn2img
        self.bpm = None             # passed as an argument to process(), flat_field()
        self.pixel_flat = None      # passed as an argument to process(), flat_field()
        self.illum_flat = None        # passed as an argument to process(), flat_field()

    # TODO: This is currently only use by BiasFrame as a test.  Now that
    # ProcessImages takes in these other parameter sets, we'll need to
    # be more clever as to how to do this, or create methods specific to
    # each child class.
#    @classmethod
#    def from_fits(cls, fits_file, **kwargs):
#        """
#        Instantiate from a FITS file
#
#        Parameters
#        ----------
#        fits_file : str
#        kwargs : passed to the __init__
#
#        Returns
#        -------
#        slf
#
#        """
#        if not os.path.isfile(fits_file):
#            msgs.error("FITS file not found: {:s}".format(fits_file))
#        # Load
#        hdul = fits.open(fits_file)
#        head0 = hdul[0].header
#        file_list = []
#        for key in head0:
#            if 'FRAME' in key:
#                file_list.append(head0[key])
#
#        # Initialize
#        slf = cls(head0['INSTRUME'], file_list=file_list,
#                  overscan_par=OverscanPar.from_header(head0),
#                  combine_par=CombineFramesPar.from_header(head0),
#                  lacosmic_par=LACosmicPar.from_header(head0), **kwargs)
#
#        # Fill
#        slf.stack = hdul[0].data
#        return slf

    # JFH I'm not following why these are properties and not simply attributes.
    @property
    def nfiles(self):
        """
        Number of files in the file_list

        Returns
        -------
        nfiles : int

        """
        return len(self.file_list)

    @property
    def nloaded(self):
        """
        Number of raw images loaded

        Returns
        -------
        nloaded : int

        """
        return len(self.raw_images)

    def load_images(self):
        """ Load raw images from the disk
        Also loads the datasec info

        Returns
        -------
        self.raw_images : list
        self.headers : list
        """
        # Load the image data and headers
        self.raw_images = []  # Zeros out any previous load
        self.headers = []
        for ifile in self.file_list:
            img, head = self.spectrograph.load_raw_frame(ifile, det=self.det)
            self.raw_images.append(img)
            self.headers.append(head)
        # Get the data sections
        self.datasec, one_indexed, include_end, transpose \
            = self.spectrograph.get_image_section(self.file_list[0], self.det,section='datasec')
        self.datasec = [ parse.sec2slice(sec, one_indexed=one_indexed,
                                           include_end=include_end, require_dim=2,
                                           transpose=transpose) for sec in self.datasec ]
        # Get the overscan sections
        self.oscansec, one_indexed, include_end, transpose = \
            self.spectrograph.get_image_section(self.file_list[0], self.det,section='oscansec')
        self.oscansec = [ parse.sec2slice(sec, one_indexed=one_indexed,
                                            include_end=include_end, require_dim=2,
                                            transpose=transpose) for sec in self.oscansec ]
        # Step
        self.steps.append(inspect.stack()[0][3])

    def apply_gain(self, trim=True):
        """
        Apply gain (instead of ampsec scale)

        Parameters
        ----------

        Returns
        -------
        self.mspixelflat -- Modified internally

        """
        # TODO: This is overkill when self.datasec is loaded, and this
        # call is made for a few of the steps.  Can we be more
        # efficient?
        datasec_img = self.spectrograph.get_datasec_img(self.file_list[0], det=self.det)
        if trim:
            datasec_img = procimg.trim_frame(datasec_img, datasec_img < 1)
        if self.stack.shape != datasec_img.shape:
            raise ValueError('Shape mismatch: {0} {1}'.format(self.stack.shape, datasec_img.shape))
        
        gain = self.spectrograph.detector[self.det-1]['gain']
        if self.spectrograph.detector[self.det-1]['numamplifiers'] == 1 \
                and not isinstance(gain, list):
            gain = [gain]
        self.stack *= procimg.gain_frame(datasec_img,
                                           self.spectrograph.detector[self.det-1]['numamplifiers'],
                                           gain)
        # Step
        self.steps.append(inspect.stack()[0][3])

        # Return
        return self.stack

    def bias_subtract(self, msbias, trim=True, force=False, par=None):
        """
        Subtract the bias.

        Parameters
        ----------
        msbias : ndarray, str (optional)
          If ndarray, the input is a Bias image
          If str, the input is guides the Bias subtraction method e.g.  'overscan'
        trim

        Returns
        -------

        """
        # Check if the bias has already been subtracted
        if (inspect.stack()[0][3] in self.steps) & (not force):
            msgs.warn("Images already bias subtracted.  Use force=True to reset proc_images "
                      "and do it again. Returning...")
            return

        # Set the parameters
        if par is not None and not isinstance(par, pypeitpar.ProcessImagesPar):
            raise TypeError('Provided ParSet for must be type ProcessImagesPar.')
        if par is not None:
            self.proc_par = par

        # If trimming, get the image identifying amplifier used for the
        # data section
        datasec_img = self.spectrograph.get_datasec_img(self.file_list[0], det=self.det)

        msgs.info("Bias subtracting your image(s)")
        # Reset proc_images -- Is there any reason we wouldn't??
        numamplifiers = self.spectrograph.detector[self.det-1]['numamplifiers']
        for kk,image in enumerate(self.raw_images):
            # Bias subtract (move here from procimg)
            if isinstance(msbias, np.ndarray):
                msgs.info("Subtracting bias image from raw frame")
                # Trim?
                if trim:
                    image = procimg.trim_frame(image, datasec_img < 1)
                temp = image-msbias
            elif isinstance(msbias, str) and msbias == 'overscan':
                msgs.info("Using overscan to subtact")
                #debugger.set_trace()
                temp = procimg.subtract_overscan(image, numamplifiers, self.datasec,
                                                   self.oscansec,
                                                   method=self.proc_par['overscan'],
                                                   params=self.proc_par['overscan_par'])
                # Trim?
                if trim:
                    temp = procimg.trim_frame(temp, datasec_img < 1)
            else:
                msgs.error('Could not subtract bias level with the input bias approach.')

            # Save
            if kk==0:
                # Instantiate proc_images
                self.proc_images = np.zeros((temp.shape[0], temp.shape[1], self.nloaded))
            self.proc_images[:,:,kk] = temp.copy()

        # Step
        self.steps.append(inspect.stack()[0][3])

    def combine(self, par=None):
        """
        Combine the processed images

        Returns
        -------
        self.stack : ndarray

        """
        # Set the parameters
        if par is not None and not isinstance(par, pypeitpar.ProcessImagesPar):
            raise TypeError('Provided ParSet for must be type ProcessImagesPar.')
        if par is not None:
            self.proc_par = par

        # Now we can combine
        saturation = self.spectrograph.detector[self.det-1]['saturation']
        self.stack = combine.comb_frames(self.proc_images, frametype=self.frametype,
                                             saturation=saturation,
                                             method=self.proc_par['combine'],
                                             satpix=self.proc_par['satpix'],
                                             cosmics=self.proc_par['sigrej'],
                                             n_lohi=self.proc_par['n_lohi'],
                                             sig_lohi=self.proc_par['sig_lohi'],
                                             replace=self.proc_par['replace'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.stack



    def flat_field(self, pixel_flat, bpm, illum_flat=None):
        """
        Flat field the stack image

        pixel_flat and illum_flat are passed here to force users to
        consider that they're needed when calling flat_field().

        Wrapper to arflat.flatfield()

        Returns
        -------
        self.stack : ndarray
          Flat fielded

        """
        # Assign the relevant data to self
        self.pixel_flat = pixel_flat
        self.bpm = bpm
        self.illum_flat = illum_flat

        # Check that the bad-pixel mask is available
        if self.bpm is None:
            msgs.error('No bpm for {0}'.format(self.spectrograph.spectrograph))

        msgs.info("Flat fielding your image")
        # Flat-field the data and return the result

        self.stack = flat.flatfield(self.stack, self.pixel_flat, self.bpm, illum_flat=self.illum_flat)
        return self.stack

    def process(self, bias_subtract=None, apply_gain=False, trim=True, overwrite=False,
                pixel_flat=None, bpm=None, illum_flat=None):
        """
        Process the images from loading to combining

        Parameters
        ----------
        bias_subtract : str or ndarray or None
          Guides bias subtraction
        apply_gain : bool, optional
          Apply gain to the various amplifier regions
        trim : bool, optional
        overwrite :
        pixel_flat : ndarray or None
          This is the normalized pixel flat (i.e. no blaze, no slit illumination profile).
          The values of this array should have a scatter about 1.0

        Returns
        -------
        self.stack : ndarray

        """

        # Over-write?
        if (inspect.stack()[0][3] in self.steps) & (not overwrite):
            msgs.warn("Images already combined.  Use overwrite=True to do it again.")
            return

        # Load images
        if 'load_images' not in self.steps:
            self.load_images()

        # Bias subtract
        if bias_subtract is not None:
            self.bias_subtract(bias_subtract, trim=trim)
        elif 'bias_subtract' not in self.steps:
            msgs.warn("Your images have not been bias subtracted!")

        # Create proc_images from raw_images if need be
        #   Mainly if no bias subtraction was performed
        if self.proc_images is None:
            # Trim even if not bias subtracting
            temp = self.raw_images[0]
            if trim:
                datasec_img = self.spectrograph.get_datasec_img(self.file_list[0], det=self.det)
                temp = procimg.trim_frame(temp, datasec_img < 1)
            self.proc_images = np.zeros((temp.shape[0], temp.shape[1], self.nloaded))
            for kk,image in enumerate(self.raw_images):
                self.proc_images[:,:,kk] = procimg.trim_frame(image, datasec_img < 1) \
                                                if trim else image

        # Combine
        self.stack = self.proc_images[:,:,0] if self.proc_images.shape[2] == 1 else self.combine()
        self.raw_stack = self.stack

        # Apply gain?
        if apply_gain:
            self.apply_gain(trim=trim)

        # Flat field?
        if pixel_flat is not None:
            self.stack = self.flat_field(pixel_flat, bpm, illum_flat=illum_flat)

        # Done
        return self.stack.copy()

    # TODO sort out dark current here. Need to pass exposure time for that.
    def build_rn2img(self, trim=True):
        """
        Generate the model read noise squared image

        Currently only used by ScienceImage.

        Wrapper to procimg.rn_frame

        Returns
        -------
        self.rn2img : ndarray

        """
        msgs.info("Generating read noise image from detector properties and amplifier layout)")
        datasec_img = self.spectrograph.get_datasec_img(self.file_list[0], det=self.det)
        if trim:
            datasec_img = procimg.trim_frame(datasec_img, datasec_img < 1)
        detector = self.spectrograph.detector[self.det-1]
        self.rn2img = procimg.rn_frame(datasec_img, detector['gain'], detector['ronoise'], numamplifiers=detector['numamplifiers'])

        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.rn2img

    def build_rawvarframe(self, trim=True):
        """
        Generate the Raw Variance frame

        Currently only used by ScienceImage.

        Wrapper to procimg.variance_frame

        Returns
        -------
        self.rawvarframe : ndarray

        """
        msgs.info("Generating raw variance frame (from detected counts [flat fielded])")
        datasec_img = self.spectrograph.get_datasec_img(self.file_list[0], det=self.det)
        if trim:
            datasec_img = procimg.trim_frame(datasec_img, datasec_img < 1)
        detector = self.spectrograph.detector[self.det-1]
        self.rawvarframe = procimg.variance_frame(datasec_img, self.stack,
                                                    detector['gain'], detector['ronoise'],
                                                    numamplifiers=detector['numamplifiers'],
                                                    darkcurr=detector['darkcurr'],
                                                    exptime=self.exptime)

        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.rawvarframe


    # This is a static method because I need to be able to run it from outside the class and would prefer
    # to not have to create an instance of the class everytime I want to do that.
    @staticmethod
    def build_crmask(stack, proc_par, det, spectrograph, ivar=None, binning=None):
        """
        Generate the CR mask frame

        Wrapper to procimg.lacosmic

        Parameters
        ----------
        varframe : ndarray, optional

        Returns
        -------
        self.crmask : ndarray
          1. = Masked CR

        """
        # Run LA Cosmic to get the cosmic ray mask
        varframe = utils.calc_ivar(ivar)
        saturation = spectrograph.detector[det-1]['saturation']
        nonlinear = spectrograph.detector[det-1]['nonlinear']
        sigclip, objlim = spectrograph.get_lacosmics_par(proc_par,binning=binning)
        crmask = procimg.lacosmic(det, stack, saturation, nonlinear,
                                  varframe=varframe, maxiter=proc_par['lamaxiter'],
                                  grow=proc_par['grow'],
                                  remove_compact_obj=proc_par['rmcompact'],
                                  sigclip=sigclip,
                                  sigfrac=proc_par['sigfrac'],
                                  objlim=objlim)

        # Return
        return crmask

    @classmethod
    def build_mask(cls, sciimg, sciivar, crmask, bpm, saturation=1e10, mincounts=-1e10, slitmask=None):
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
        mask = np.zeros_like(sciimg, dtype=cls.bitmask.minimum_dtype(asuint=True))

        # Bad pixel mask
        indx = bpm.astype(bool)
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'BPM')

        # Cosmic rays
        indx = crmask.astype(bool)
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'CR')

        # Saturated pixels
        indx = sciimg >= saturation
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'SATURATION')

        # Minimum counts
        indx = sciimg <= mincounts
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'MINCOUNTS')

        # Undefined counts
        indx = np.invert(np.isfinite(sciimg))
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'IS_NAN')

        # Bad inverse variance values
        indx = np.invert(sciivar > 0.0)
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'IVAR0')

        # Undefined inverse variances
        indx = np.invert(np.isfinite(sciivar))
        mask[indx] = cls.bitmask.turn_on(mask[indx], 'IVAR_NAN')

        if slitmask is not None:
            indx = slitmask == -1
            mask[indx] = cls.bitmask.turn_on(mask[indx], 'OFFSLITS')

        return mask

    @classmethod
    def update_mask_cr(cls, mask_old, crmask_new):

        # Unset the CR bit from all places where it was set
        CR_old = (cls.bitmask.unpack(mask_old, flag='CR'))[0]
        mask_new = np.copy(mask_old)
        mask_new[CR_old] = cls.bitmask.turn_off(mask_new[CR_old], 'CR')
        # Now set the CR bit using the new crmask
        indx = crmask_new.astype(bool)
        mask_new[indx] = cls.bitmask.turn_on(mask_new[indx], 'CR')
        return mask_new


    # Do we still need this function?
    @classmethod
    def update_mask_slitmask(cls, mask_old, slitmask):

        # Pixels excluded from any slit.
        mask_new = np.copy(mask_old)
        indx = slitmask == -1
        mask_new[indx] = cls.bitmask.turn_on(mask_new[indx], 'OFFSLITS')
        return mask_new

    @classmethod
    def read_stack(cls, files, bias, pixel_flat, bpm, det, proc_par, spectrograph, illum_flat=None, reject_cr=False,
                   binning=None):
        """  Utility function for reading in image stacks using ProcessImages
        Parameters
            file_list:
            bias:
            pixel_flat:
            bpm:
            illum_flat:
        Returns:
        """
        nfiles = len(files)
        for ifile in range(nfiles):
            this_proc = ProcessImages(spectrograph, [files[ifile]], det=det, par=proc_par)
            # TODO I think trim should be hard wired, and am not letting it be a free parameter
            sciimg = this_proc.process(bias_subtract=bias,pixel_flat=pixel_flat, illum_flat=illum_flat, bpm=bpm,
                                       apply_gain=True, trim=True)
            # Allocate the images
            if ifile == 0:
                # numpy is row major so stacking will be fastest with nfiles as the first dimensions
                shape = (nfiles, sciimg.shape[0],sciimg.shape[1])
                sciimg_stack  = np.zeros(shape)
                sciivar_stack = np.zeros(shape)
                rn2img_stack  = np.zeros(shape)
                crmask_stack  = np.zeros(shape,dtype=bool)
                mask_stack  = np.zeros(shape,this_proc.bitmask.minimum_dtype(asuint=True))

            # Construct raw variance image
            rawvarframe = this_proc.build_rawvarframe(trim=True)
            # Mask cosmic rays
            sciivar_stack[ifile,:,:] =  utils.calc_ivar(rawvarframe)
            if reject_cr:
                crmask_stack[ifile,:,:] = this_proc.build_crmask(sciimg, proc_par, det, spectrograph,
                                                                 ivar=sciivar_stack[ifile,:,:], binning=binning)
            sciimg_stack[ifile,:,:] = sciimg
            # Build read noise squared image
            rn2img_stack[ifile,:,:] = this_proc.build_rn2img()
            # Final mask for this image
            mask_stack[ifile,:,:] = this_proc.build_mask(sciimg, sciivar_stack[ifile,:,:], crmask_stack[ifile,:,:],
                                                         bpm, saturation = spectrograph.detector[det - 1]['saturation'],
                                                         mincounts = spectrograph.detector[det - 1]['mincounts'])

        return sciimg_stack, sciivar_stack, rn2img_stack, crmask_stack, mask_stack

    def show(self, attr='stack', idx=None, display='ginga'):
        """
        Show an image

        Parameters
        ----------
        attr : str, optional
          Internal name of the image to show
            proc_image, raw_image, stack
        idx : int, optional
          Specifies the index of the raw or processed image
          Required if proc_image or raw_image is called
        display : str

        Returns
        -------

        """
        if 'proc_image' in attr:
            img = self.proc_images[:,:,idx]
        elif 'raw_image' in attr:
            img = self.raw_images[idx]
        elif 'stack' in attr:
            img = self.stack
        else:
            msgs.warn("Options:  proc_image, raw_image, stack")
            return
        # Show
        viewer, ch = ginga.show_image(img)

    def write_stack_to_fits(self, outfile, overwrite=True):
        """
        Write the combined image to disk as a FITS file

        Parameters
        ----------
        outfile : str
        overwrite

        Returns
        -------

        """
        if self.stack is None:
            msgs.warn("You need to generate the stack before you can write it!")
            return
        #
        hdu = fits.PrimaryHDU(self.stack)
        # Add raw_files to header
        for i in range(self.nfiles):
            hdrname = "FRAME{0:03d}".format(i+1)
            hdu.header[hdrname] = self.file_list[i]
        # Spectrograph
        hdu.header['INSTRUME'] = self.spectrograph.spectrograph
        # Parameters
        self.proc_par.to_header(hdu.header)
        # Steps
        steps = ','
        hdu.header['STEPS'] = steps.join(self.steps)
        # Finish
        hlist = [hdu]
        hdulist = fits.HDUList(hlist)
        hdulist.writeto(outfile, overwrite=overwrite)
        msgs.info("Wrote stacked image to {:s}".format(outfile))

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





# TODO Add a function here that just reads in a fits file given a filename. Guess the instrument from the headers.