# Module for generating the Arc image
from __future__ import absolute_import, division, print_function

import os
import inspect
import numpy as np

from pypeit import msgs
from pypeit import processimages
from pypeit import masterframe
from pypeit.core import sort
from pypeit.par import pypeitpar

from pypeit import debugger

class ArcImage(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate an Arc Image from one or more arc frames
      The master() method returns the image (loaded or built)

    Parameters
    ----------
    file_list : list (optional)
      List of filenames
    spectrograph : str (optional)
       Used to specify properties of the detector (for processing)
       Passed to ProcessImages
       Attempts to set with settings['run']['spectrograph'] if not input
    settings : dict (optional)
       Passed to ProcessImages
       Settings for image combining+detector
    setup : str (optional)
      Setup tag;  required for MasterFrame functionality
    det : int, optional
      Detector index, starts at 1
    sci_ID : int (optional)
      Science ID value
      used to match bias frames to the current science exposure
    msbias : ndarray or str
      Guides bias subtraction
    fitstbl : Table (optional)
      FITS info (mainly for filenames)

    Attributes
    ----------
    frametype : str
      Set to 'arc'

    Inherited Attributes
    --------------------
    stack : ndarray
      Final output image
    """

    # Frametype is a class attribute
    frametype = 'arc'

    def __init__(self, spectrograph, file_list=[], det=1, par=None, setup=None, root_path=None,
                 mode=None, fitstbl=None, sci_ID=None, msbias=None):
    
        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        directory_path = None if root_path is None \
                                else root_path+'_'+self.spectrograph.spectrograph
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         directory_path=directory_path, mode=mode)


    def build_image(self):
        """
        Build the arc image from one or more arc files

        Returns
        -------

        """
        # Get list of arc frames for this science frame
        #  unless one was input already
        if self.nfiles == 0:
            self.file_list = sort.list_of_files(self.fitstbl, self.frametype, self.sci_ID)
        # Combine
        self.stack = self.process(bias_subtract=self.msbias)
        #
        return self.stack

    # TODO: There is no master() method.  Does this mean useframe is
    # always 'arc'?


