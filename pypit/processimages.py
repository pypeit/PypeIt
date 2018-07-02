# Module for Processing Images, e.g. bias frames, arc frames
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arcomb
from pypit.spectrographs import io
from pypit.core import arprocimg
from pypit.core import arflat
from pypit import ginga


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
def default_settings():
    # Hiding in a def because of the nested dicts
    #  I prefer this to deepcopy
    default_settings = dict(detector={'numamplifiers': 1,   # This dict is not complete; consider readnoise, binning
                                  'saturation': 60000.,  # Spectra aligned with columns
                                  'dispaxis': 0,  # Spectra aligned with columns
                                  'dataext': None,
                                  'gain': [1.],
                                      },
                        combine={'match': -1.,
                                 'satpix': 'reject',
                                 'method': 'weightmean',
                                 'reject': {'cosmics': 20.,
                                            'lowhigh': [0,0],
                                            'level': [3.,3.],
                                            'replace': 'maxnonsat'}}
                        )
    return default_settings


class ProcessImages(object):
    """Base class to guide image loading+processing

    Parameters
    ----------
    file_list : list
      List of filenames
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    det : int, optional
      Detector index, starts at 1
    settings : dict, optional
      Settings for trace slits
    user_settings : dict, optional
      Allow for user to over-ride individual internal/default settings
      without providing a full settings dict
    datasec_img : ndarray, optional
      Specifies which pixels go with which amplifier

    Attributes
    ----------
    stack : ndarray
    steps : list
    raw_images : list
    headers : list
    proc_images : ndarray
      3D array of processed, individual images
    datasec : list
    oscansec : list
    """
    def __init__(self, file_list, spectrograph=None, settings=None, det=1, user_settings=None,
                 datasec_img=None, bpm=None):

        # Required parameters
        if not isinstance(file_list, list):
            raise IOError("file_list input to ProcessImages must be list. Empty is fine")
        self.file_list = file_list

        # Optional
        self.det = det
        self.datasec_img = datasec_img
        self.bpm = bpm

        if settings is None:
            self.settings = default_settings()
        else:
            # The copy allows up to update settings with user settings without changing the original
            self.settings = settings.copy()
        if user_settings is not None:
            # This only works to replace entire dicts
            #    Hopefully parsets will work more cleanly..
            self.settings.update(user_settings)
        self.frametype='Unknown'

        # Main (possible) outputs
        self.stack = None
        self.steps = []

        # Key Internals
        if spectrograph is None:
            try:
                spectrograph = self.settings['run']['spectrograph']
            except:
                msgs.warn("No information on the spectrograph was given.  Do not attempt to (re)process the images")
        self.spectrograph = spectrograph
        self.raw_images = []
        self.proc_images = None  # Will be an ndarray
        self.headers = []
        self.datasec = []
        self.oscansec = []

    @classmethod
    def from_fits(cls, fits_file, **kwargs):
        """
        Instantiate from a FITS file

        Parameters
        ----------
        fits_file : str
        kwargs : passed to the __init__

        Returns
        -------
        slf

        """
        if not os.path.isfile(fits_file):
            msgs.error("FITS file not found: {:s}".format(fits_file))
        # Load
        hdul = fits.open(fits_file)
        head0 = hdul[0].header
        file_list = []
        for key in head0:
            if 'FRAME' in key:
                file_list.append(head0[key])
        slf = cls(file_list, **kwargs)
        # Fill
        slf.stack = hdul[0].data
        return slf

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
        self.raw_images = []  # Zeros out any previous load
        self.headers = []
        for ifile in self.file_list:
            img, head = io.load_raw_frame(self.spectrograph, ifile, self.det,
                                        dataext=self.settings['detector']['dataext'],
                                        disp_dir=self.settings['detector']['dispaxis'])
            # Save
            self.raw_images.append(img)
            self.headers.append(head)
        # Grab datasec (almost always desired)
        self._grab_datasec(redo=True)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _grab_datasec(self, redo=False):
        """
        Load the datasec parameters

        Parameters
        ----------
        redo : bool, optional

        Returns
        -------
        self.datasec : list
        self.oscansec : list
        """
        if (self.datasec is not None) and (not redo):
            return
        # Spectrograph specific
        self.datasec, self.oscansec, _, _ = io.get_datasec(self.spectrograph,
                                                     self.file_list[0], self.det,
                                                     self.settings['detector'])

    def apply_gain(self, datasec_img):
        """
        # Apply gain (instead of ampsec scale)

        Parameters
        ----------
        datasec_img : ndarray
          Defines which pixels belong to which amplifier

        Returns
        -------
        self.mspixelflat -- Modified internally

        """
        if datasec_img is None:
            msgs.error("Need to input datasec_img!")
        self.stack = self.stack*arprocimg.gain_frame(datasec_img,
                                                 self.settings['detector']['numamplifiers'],
                                                 self.settings['detector']['gain'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.stack


    def bias_subtract(self, msbias, trim=True, force=False):
        """

        Parameters
        ----------
        msbias : ndarray, str (optional)
          If ndarray, the input is a Bias image
          If str, the input is guides the Bias subtraction method e.g.  'overscan'
        trim

        Returns
        -------

        """
        if (inspect.stack()[0][3] in self.steps) & (not force):
            msgs.warn("Images already bias subtracted.  Use force=True to reset proc_images and do it again.")
            msgs.warn("Returning..")
            return
        msgs.info("Bias subtracting your image(s)")
        # Reset proc_images -- Is there any reason we wouldn't??
        for kk,image in enumerate(self.raw_images):
            # Bias subtract
            temp = arprocimg.bias_subtract(image, msbias,
                                        numamplifiers=self.settings['detector']['numamplifiers'],
                                        datasec=self.datasec, oscansec=self.oscansec)
            # Trim?
            if trim:
                temp = arprocimg.trim(temp, self.settings['detector']['numamplifiers'], self.datasec)

            # Save
            if kk==0:
                # Instantiate proc_images
                self.proc_images = np.zeros((temp.shape[0], temp.shape[1], self.nloaded))
            self.proc_images[:,:,kk] = temp.copy()

        # Step
        self.steps.append(inspect.stack()[0][3])

    def combine(self):
        """
        Combine the processed images

        Returns
        -------
        self.stack : ndarray

        """
        # Now we can combine
        self.stack = arcomb.core_comb_frames(self.proc_images, frametype=self.frametype,
                                             method=self.settings['combine']['method'],
                                             reject=self.settings['combine']['reject'],
                                             satpix=self.settings['combine']['satpix'],
                                             saturation=self.settings['detector']['saturation'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.stack

    def build_crmask(self, varframe=None):
        """
        Generate the CR mask frame

        Wrapper to arprocimg.lacosmic

        Parameters
        ----------
        varframe : ndarray, optional

        Returns
        -------
        self.crmask : ndarray
          1. = Masked CR

        """
        self.crmask = arprocimg.lacosmic(self.det, self.stack, self.settings['detector'],
                                    grow=1.5, varframe=varframe)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.crmask

    def flat_field(self):
        """
        Flat field the stack image

        Wrapper to arflat.flatfield()

        Returns
        -------
        self.stack : ndarray
          Flat fielded

        """
        if self.bpm is None:
            msgs.error("Need to set the BPM image, even if all zeros")
        self.stack = arflat.flatfield(self.stack, self.pixel_flat, self.bpm, slitprofile=self.slitprof)
        return self.stack


    def process(self, bias_subtract=None, apply_gain=False,
                trim=True, overwrite=False,
                pixel_flat=None, slitprof=None):
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
        if (bias_subtract is not None):
            self.bias_subtract(bias_subtract, trim=trim)
        else:
            if 'bias_subtract' not in self.steps:
                msgs.warn("Your images have not been bias subtracted!")

        # Create proc_images from raw_images if need be
        #   Mainly if no bias subtraction was performed
        if self.proc_images is None:
            self.proc_images = np.zeros((self.raw_images[0].shape[0],
                                         self.raw_images[0].shape[1],
                                         self.nloaded))
            for kk,image in enumerate(self.raw_images):
                self.proc_images[:,:,kk] = image
        # Combine
        if self.proc_images.shape[2] == 1:  # Nothing to combine
            self.stack = self.proc_images[:,:,0]
        else:
            self.stack = self.combine()

        # Apply gain?
        if apply_gain:
            self.apply_gain(self.datasec_img)

        # Flat field?
        if pixel_flat is not None:
            self.pixel_flat = pixel_flat
            self.slitprof = slitprof
            self.stack = self.flat_field()

        # Done
        return self.stack.copy()

    def build_rawvarframe(self, dnoise=None):
        """
        Generate the Raw Variance frame

        Wrapper to arprocimg.variance_frame

        Parameters
        ----------
        dnoise : float
          Noise related to dark current (generally 0.)

        Returns
        -------
        self.rawvarframe : ndarray

        """
        msgs.info("Generate raw variance frame (from detected counts [flat fielded])")
        self.rawvarframe = arprocimg.variance_frame(self.datasec_img, self.det, self.stack,
                                               self.settings['detector'], dnoise=dnoise)
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.rawvarframe

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


