""" Module for guiding Bias subtraction including generating a Bias image as desired
"""
from __future__ import absolute_import, division, print_function

import inspect
import os


from pypeit import msgs
from pypeit import processimages
from pypeit import masterframe
from pypeit.par import pypeitpar

from pypeit import debugger


class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate a Bias frame for bias subtraction
      It also contains I/O methods for the Master frames of PypeIt
      The build_master() method will return a simple command (str) if that is the specified setting
      in settings['bias']['useframe']

    Instead have this comment and more description here:
        # Child-specific Internals
        #    See ProcessImages for the rest

    Parameters
    ----------
    file_list : list (optional)
      List of filenames
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Attempt to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
      Settings for trace slits
    setup : str (optional)
      Setup tag
    det : int, optional
      Detector index, starts at 1
    fitstbl : PypeItMetaData (optional)
      FITS info (mainly for filenames)
    sci_ID : int (optional)
      Science ID value
      used to match bias frames to the current science exposure
    par : ParSet
      PypitPar['calibrations']['biasframe']
    redux_path : str (optional)
      Path for reduction


    Attributes
    ----------
    frametype : str
      Set to 'bias'

    Inherited Attributes
    --------------------
    stack : ndarray
    """

    # Frame type is a class attribute
    frametype = 'bias'

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, file_list = [], det=1, par=None, setup=None, master_dir=None,
                 mode=None):

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, setup, mode=mode,
                                         master_dir=master_dir)

    def build_image(self, overwrite=False, trim=True):
        """
        Grab the bias files (as needed) and then
         process the input bias frames with ProcessImages.process()
          Avoid bias subtraction
          Avoid trim

        Parameters
        ----------
        overwrite : bool, optional

        Returns
        -------
        stack : ndarray

        """
        # Combine
        self.stack = self.process(bias_subtract=None, trim=trim, overwrite=overwrite)
        #
        return self.stack

    def determine_bias_mode(self):
        # How are we treating biases?
        # 1) No bias subtraction
        if self.par['useframe'].lower() == 'none':
            msgs.info("Will not perform bias/dark subtraction")
            self.msbias = None
        # 2) Use overscan
        elif self.par['useframe'] == 'overscan':
            self.msbias = 'overscan'
        # 3) User wants bias subtractions, use a Master biasframe?
        elif self.par['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists and user requested one to load it
            self.msbias = self.master()

        return self.msbias
