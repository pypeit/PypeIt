"""
Module for guiding Bias subtraction including generating a Bias image as desired

.. _numpy.ndarray: https://docs.scipy.org/doc/numpy/reference/generated/numpy.ndarray.html

"""
import numpy as np

from pypeit import msgs
from pypeit import processimages
from pypeit import masterframe
from pypeit.par import pypeitpar

class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    Class to generate/load the Bias image or instructions on how to deal
    with the bias.

    This class is primarily designed to generate a Bias frame for bias
    subtraction.  It also contains I/O methods for the Master frames of
    PypeIt.  The build_master() method will return a simple command
    (str) if that is the specified parameter (`par['useframe']`).

    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            Spectrograph used to take the data.

        files (:obj:`list`, optional):
            List of filenames to process.
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):
            The parameters used to process the frames.  If None, set
            to::
                
                pypeitpar.FrameGroupPar('bias')

        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        reuse_masters (:obj:`bool`, optional):
            Load from disk if possible
    """

    # Frame type is a class attribute
    frametype = 'bias'
    master_type = 'Bias'

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False):

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, self.par['process'],
                                             files=files, det=det)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

    def build_image(self, overwrite=False, trim=True):
        """
        Grab the bias files (as needed) and then process the input bias
        frames with :func:`pypeit.processimages.ProcessImages.process`.

        Args:
            overwrite: (:obj: `bool`, optional):
                Regenerate the stack image
            trim (:obj:`bool`, optional):
                If True, trim the image

        Returns:
            `numpy.ndarray`_: Combined, processed image.
        """
        if self.par['useframe'] == 'none':
            return None
        return self.process(bias_subtract=None, trim=trim, overwrite=overwrite)

    def save(self, outfile=None, overwrite=True):
        """
        Save the bias master data.

        Args:
            outfile (:obj:`str`, optional):
                Name for the output file.  Defaults to
                :attr:`file_path`.
            overwrite (:obj:`bool`, optional):
                Overwrite any existing file.
        """
        if self.stack is None:
            msgs.warn('No MasterBias to save!')
            return
        if not isinstance(self.stack, np.ndarray):
            msgs.warn('MasterBias is not an image.')
            return
        super(BiasFrame, self).save(self.stack, 'BIAS', outfile=outfile, overwrite=overwrite,
                                    raw_files=self.files, steps=self.steps)

    # TODO: it would be better to have this instantiate the full class
    # as a classmethod.
    def load(self, ifile=None, return_header=False):
        """
        Load the bias frame.
        
        This overwrites :func:`pypeit.masterframe.MasterFrame.load` so
        that the bias can be returned as a string as necessary.

        The bias mode to use in this reduction is either
          - None -- No bias subtraction
          - 'overscan' -- Overscan subtract
          - stack -- Use a generated bias image

        The result is *not* saved internally.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            Returns either the `numpy.ndarray`_ with the bias image, a
            string specifying that the bias should be subtracted using
            the overscan, or None if no bias is to be subtracted.
        """
        # Check input
        if self.par['useframe'].lower() in ['none', 'overscan'] and return_header:
            msgs.warn('No image data to read.  Header returned as None.')

        # How are we treating biases?
        # 1) No bias subtraction
        if self.par['useframe'].lower() == 'none':
            msgs.info("Will not perform bias/dark subtraction")
            return (None,None) if return_header else None

        # 2) Use overscan
        if self.par['useframe'] == 'overscan':
            return ('overscan', None) if return_header else 'overscan'

        # 3) User wants bias subtractions, use a Master biasframe?
        if self.par['useframe'] in ['bias', 'dark']:
            return super(BiasFrame, self).load('BIAS', ifile=ifile, return_header=return_header)

