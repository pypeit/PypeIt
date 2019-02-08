""" Module for generating the Arc image"""
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
    Generate an Arc Image by processing and combining one or more arc frames.

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        files (:obj:`list`, optional):
            The list of files to process.  Can be an empty list.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (str, optional): Path to master frames
        reuse_masters (bool, optional): Load from disk if possible
        msbias (ndarray or str, optional): Guides bias subtraction

    Attributes:
        stack (ndarray): Final output image

    """

    # Frametype is a class attribute
    frametype = 'arc'

    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False, msbias=None):
    
        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, self.par['process'],
                                             files=files, det=det)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, master_key, master_dir,
                                         reuse_masters=reuse_masters)

    def build_image(self, overwrite=False):
        """ Build the arc image from one or more arc files

        Args:
            overwrite: (:obj: `bool`, optional):
                Regenerate the stack image

        Returns:
            ndarray: :attr:`stack` Combined, processed image
            
        """
        # Combine
        self.stack = self.process(bias_subtract=self.msbias, overwrite=overwrite, trim=True)
        #
        return self.stack


