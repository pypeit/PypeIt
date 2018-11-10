# Module for generating the Arc image
from __future__ import absolute_import, division, print_function, unicode_literals

import os
import inspect
import numpy as np

from pypeit import msgs
from pypeit import processimages
from pypeit import masterframe
from pypeit.par import pypeitpar

from pypeit import debugger

class ArcImage(processimages.ProcessImages, masterframe.MasterFrame):
    """
    Generate an Arc Image from one or more arc frames.

    Args:
        spectrograph (:obj:`str`,
            :class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The string or `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        file_list (:obj:`list`, optional):
            The list of files to process.  Can be an empty list.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.
        setup (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.

        root_path (:obj:`str`, optional):
            
    sci_ID : int (optional)
      Science ID value
      used to match bias frames to the current science exposure
    msbias : ndarray or str
      Guides bias subtraction
    fitstbl : PypeItMetaData (optional)
      FITS info (mainly for filenames)
    redux_path : str (optional)
      Path for reduction

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

    def __init__(self, spectrograph, file_list, det=1, par=None, setup=None,
                 master_dir=None, mode=None, msbias=None):
    
        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, setup,
                                         mode=mode, master_dir=master_dir)


    def build_image(self, overwrite=False, trim=True):
        """
        Build the arc image from one or more arc files

        Returns
        -------

        """
        # Get list of arc frames for this science frame
        #  unless one was input already
        # Combine
        self.stack = self.process(bias_subtract=self.msbias, overwrite=overwrite, trim=True)
        #
        return self.stack

    # TODO: There is no master() method.  Does this mean useframe is
    # always 'arc'?


