# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from collections import OrderedDict

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arcomb
from pypit import arload
from pypit import arproc
from pypit import arparse
from pypit import ginga


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
default_settings = dict(run={'spectrograph': 'UNKNOWN'},
                        detector={'numamplifiers': 1,
                                  'saturation': 60000.,  # Spectra aligned with columns
                                  'dispaxis': 0,  # Spectra aligned with columns
                                  'dataext': None},
                        combine={'match': -1.,
                                 'satpix': 'reject',
                                 'method': 'weightmean',
                                 'reject': {'cosmics': 20.,
                                            'lowhigh': [0,0],
                                            'level': [3.,3.],
                                            'replace': 'maxnonsat'}}
                        )

# datasec kludge (until settings is Refactored)
# TODO -- Remove this eventually
do_sec_dict = dict(
    shane_kast_blue=OrderedDict([  # The ordering of these is important, ie. 1, 2,  hence the OrderedDict
        ('datasec01','[1:1024,:]'),        # Either the data sections (IRAF format) or the header keyword where the valid data sections can be obtained
        ('oscansec01', '[2050:2080,:]'),    # Either the overscan sections (IRAF format) or the header keyword where the valid overscan sections can be obtained
        ('datasec02', '[1025:2048,:]'),     # Either the data sections (IRAF format) or the header keyword where the valid data sections can be obtained
        ('oscansec02', '[2081:2111,:]')])    # Either the overscan sections (IRAF format) or the header keyword where the valid overscan sections can be obtained
)


class ProcessImages(object):
    """Base class to guide image loading+processing

    Parameters
    ----------
    file_list : list
      List of filenames
    det : int, optional
      Detector index, starts at 1
    settings : dict, optional
      Settings for trace slits
    user_settings : dict, optional
      Allow for user to over-ride individual internal/default settings
      without providing a full settings dict
    spectrograph : str
       Used to specify properties of the detector (for processing)
       Is set with settings['run']['spectrograph'] if not input

    Attributes
    ----------
    images : list
    stack : ndarray
    raw_images : list
    headers : list
    proc_images : ndarray
      3D array of processed, individual images
    datasec : list
    oscansec : list (optional)
    """
    def __init__(self, file_list, settings=None, det=1, user_settings=None, spectrograph=None):

        # Parameters
        self.file_list = file_list
        self.det = det

        # Settings
        if settings is None:
            self.settings = default_settings.copy()
        else:
            # The copy allows up to update settings with user settings without changing the original
            self.settings = settings.copy()
        if user_settings is not None:
            # This only works to replace entire dicts
            #    Hopefully parsets will work more cleanly..
            self.settings.update(user_settings)

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
        return len(self.file_list)

    @property
    def nloaded(self):
        return len(self.raw_images)

    def load_images(self):
        self.raw_images = []  # Zeros out any previous load
        self.headers = []
        for ifile in self.file_list:
            img, head = arload.load_raw_frame(self.spectrograph, ifile, self.det,
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
        if (self.datasec is not None) and (not redo):
            return
        # TODO -- Eliminate this instrument specific bit here. Probably by generating a Detector object
        if self.spectrograph in ['keck_lris_blue', 'keck_lris_red', 'keck_deimos']:
            self.datasec, self.oscansec, _, _ = arproc.get_datasec(
                self.spectrograph, self.file_list[0],
                numamplifiers=self.settings['detector']['numamplifiers'],
                det=self.det)
        else:
            self.datasec, self.oscansec = [], []
            if self.spectrograph not in do_sec_dict.keys():
                debugger.set_trace()
                msgs.error("NOT READY FOR THIS SPECTROGRAPH: {:s}".format(self.spectrograph))
            msgs.info("Parsing the datasec and oscansec values")
            for key in do_sec_dict[self.spectrograph]:
                if 'datasec' in key:
                    self.datasec.append(arparse.load_sections(do_sec_dict[self.spectrograph][key], fmt_iraf=False))
                elif 'oscansec' in key:
                    self.oscansec.append(arparse.load_sections(do_sec_dict[self.spectrograph][key], fmt_iraf=False))
                else:
                    pass

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
            temp = arproc.bias_subtract(image, msbias,
                                        numamplifiers=self.settings['detector']['numamplifiers'],
                                        datasec=self.datasec, oscansec=self.oscansec)
            # Trim?
            if trim:
                temp = arproc.trim(temp, self.settings['detector']['numamplifiers'], self.datasec)

            # Save
            if kk==0:
                # Instantiate proc_images
                self.proc_images = np.zeros((temp.shape[0], temp.shape[1], self.nloaded))
            self.proc_images[:,:,kk] = temp.copy()

        # Step
        self.steps.append(inspect.stack()[0][3])

    def combine(self, frametype='Unknown'):
        # Now we can combine
        self.stack = arcomb.core_comb_frames(self.proc_images, frametype=frametype,
                                             method=self.settings['combine']['method'],
                                             reject=self.settings['combine']['reject'],
                                             satpix=self.settings['combine']['satpix'],
                                             saturation=self.settings['detector']['saturation'])
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.stack

    def process(self, bias_subtract=None, trim=True, overwrite=False):
        # Over-write?
        if (inspect.stack()[0][3] in self.steps) & (not overwrite):
            msgs.warn("Images already combined.  Use overwrite=True to do it again.")
            return

        # Allow for one-stop-shopping
        if 'load_images' not in self.steps:
            self.load_images()
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
        self.stack = self.combine()
        return self.stack.copy()

    def flat_field(self):
        pass

    def show(self, attr, idx=None, display='ginga'):
        if 'proc_image' in attr:
            img = self.proc_images[:,:,idx]
        elif 'raw_image' in attr:
            img = self.raw_images[idx]
        elif 'stack' in attr:
            img = self.stack
        # Show
        viewer, ch = ginga.show_image(img)

    def write_stack_to_fits(self, outfile, overwrite=True):
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


