""" Module for guiding Bias subtraction including generating a Bias image as desired
"""
from __future__ import absolute_import, division, print_function

from pypeit import msgs
from pypeit import processimages
from pypeit import masterframe
from pypeit.par import pypeitpar

class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    Class to generate/load the Bias image or instructions on how to deal with the bias

    This class is primarily designed to generate a Bias frame for bias subtraction
      It also contains I/O methods for the Master frames of PypeIt
      The build_master() method will return a simple command (str) if that is the specified setting
      in settings['bias']['useframe']

    Instead have this comment and more description here:
        # Child-specific Internals
        #    See ProcessImages for the rest

    Args:
        spectrograph (str):
           Used to specify properties of the detector (for processing)
        files (list,optional): List of filenames to process
            master_key (:obj:`str`, optional):
                The string identifier for the instrument configuration.  See
                :class:`pypeit.masterframe.MasterFrame`.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (str, optional): Path to master frames
        reuse_masters (bool, optional): Load from disk if possible

    Attributes:

    """

    # Frame type is a class attribute
    frametype = 'bias'

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False):

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph,
                                             self.par['process'],
                                             files=files, det=det)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, master_key, master_dir, reuse_masters=reuse_masters)

    def build_image(self, overwrite=False, trim=True):
        """
        Grab the bias files (as needed) and then
         process the input bias frames with ProcessImages.process()
          Avoid bias subtraction

        Args:
            overwrite: (:obj: `bool`, optional):
                Regenerate the stack image
            trim (bool, optional): If True, trim the image

        Returns:
            ndarray: :attr:`stack` Combined, processed image

        """
        # Combine
        self.stack = self.process(bias_subtract=None, trim=trim, overwrite=overwrite)
        #
        return self.stack

    def determine_bias_mode(self, prev_build=False):
        """
        Determine the bias mode to use in this reduction
          - None -- No bias subtraction
          - 'overscan' -- Overscan subtract
          - msbias -- Use a generated bias image

        Args:
            prev_build (bool, optional):  Load the master frame if it exists and was
            built on this run.

        Returns:
            ndarray, str or None: :attr:`msbias` str, np.ndarray or None

        """
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
            self.msbias = self.master(prev_build=prev_build)

        return self.msbias
