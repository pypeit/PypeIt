"""
Module for generating the Arc image.

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
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
        msbias (ndarray):
            Bias image or bias-subtraction method; see
            :func:`pypeit.processimages.ProcessImages.process`.
    """

    # Frametype is a class attribute
    frametype = 'arc'
    master_type = 'Arc'

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
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

    # TODO: Allow trim to be a keyword argument?
    def build_image(self, overwrite=False):
        """
        Build the arc image from one or more arc files.

        Args:
            overwrite: (:obj: `bool`, optional):
                Regenerate the stack image

        Returns:
            `numpy.ndarray`_: Combined, processed image
            
        """
        return self.process(bias_subtract=self.msbias, overwrite=overwrite, trim=True)

    def save(self, outfile=None, overwrite=True):
        """
        Save the arc master data.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        super(ArcImage, self).save(self.stack, 'ARC', outfile=outfile, overwrite=overwrite,
                                   raw_files=self.files, steps=self.steps)

    # TODO: it would be better to have this instantiate the full class
    # as a classmethod.
    def load(self, ifile=None, return_header=False):
        """
        Load the arc frame data from a saved master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            Returns a `numpy.ndarray`_ with the arc master frame image.
            Also returns the primary header, if requested.
        """
        return super(ArcImage, self).load('ARC', ifile=ifile, return_header=return_header)

