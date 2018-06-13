# Module for generating the BPM image
from __future__ import absolute_import, division, print_function

import numpy as np
import os

from pypit import msgs
from pypit.core import arprocimg
from pypit.core import arlris
from pypit.core import ardeimos

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'bpm'


class BPMImage(object):
    """
    This class is primarily designed to generate an Bad Pixel Image
      The master() method will return the image

    There are several ways to build a BPM:
       1. keck_lris_red
          Provide binning and detector number to generate this instrument specific BPM
       2. keck_deimos
          Provide detector number to generate this instrument specific BPM
       3. From a bias image (not well tested)
          Set reduce_badpix='bias'
          Provide msbias image
       4. Dummy image from shape

    Parameters
    ----------
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
      Settings for trace slits
    det : int, optional
      Detector index
      Required for LRISr and DEIMOS
    binning : Table (optional)
        Could remove if we put binning in settings['detector']
    shape : tuple
      Image shape;  used to construct a dummy BPM if all else fails
    msbias : ndarray, optional
      Used to construct the BPM if reduce_badpix=='bias'
    reduce_badpix : str, optional
      'bias' -- Build from bias images

    Attributes
    ----------
    bpm : ndarray

    """
    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph=None, settings=None, det=None,
                 binning=None, reduce_badpix=None, msbias=None, shape=None):

        # Parameters unique to this Object
        self.spectrograph = spectrograph
        self.binning = binning
        self.det = det

        self.reduce_badpix = reduce_badpix
        self.msbias = msbias

        self.shape = shape
        self.settings = settings

        # Checks
        if (self.reduce_badpix == 'bias') and (self.msbias is None):
            msgs.error("Need to supply msbias image with this option")
        if (self.spectrograph == 'keck_deimos') and (self.det is None):
            msgs.error("Need to supply det with this option")
        if (self.spectrograph == 'keck_lris_red') and (self.det is None):
            msgs.error("Need to supply det with this option")
        if (self.spectrograph == 'keck_lris_red') and (self.binning is None):
            msgs.error("Need to supply binning with this option")

        # Attributes (set after init)
        self.frametype = frametype

        # Output
        self.bpm = None


    def build(self):
        """
        Generate the BPM Image

        Returns
        -------
        self.bpm

        """
        if self.reduce_badpix == 'bias':
            # Get all of the bias frames for this science frame
            if self.msbias is None:
                msgs.warn("No bias frame provided!")
                msgs.info("Not preparing a bad pixel mask")
                return False
            # TODO -- Deal better with this datasec kludge
            datasec = []
            for i in range(self.settings['detector']['numamplifiers']):
                sdatasec = "datasec{0:02d}".format(i+1)
                datasec.append(self.settings['detector'][sdatasec])
            # Construct
            self.bpm = arprocimg.badpix(self.msbias, self.settings['detector']['numamplifiers'], datasec)
        else:
            # Instrument dependent
            if self.spectrograph in ['keck_lris_red']:
                # Index in fitstbl for binning
                xbin, ybin = [int(ii) for ii in self.binning.split(',')]
                self.bpm = arlris.bpm(xbin, ybin, 'red', self.det)
            elif self.spectrograph in ['keck_deimos']:
                self.bpm = ardeimos.bpm(self.det)
            elif self.spectrograph in ['keck_nirspec']:
                # Edges of the detector are junk
                msgs.info("Custom bad pixel mask for NIRSPEC")
                self.bpm = np.zeros((self.shape[0], self.shape[1]))
                self.bpm[:, :20] = 1.
                self.bpm[:, 1000:] = 1.
            else:
                ###############
                # Set the number of spectral and spatial pixels, and the bad pixel mask is it does not exist
                self.bpm = np.zeros((self.shape[0], self.shape[1]))
        # Return
        return self.bpm


