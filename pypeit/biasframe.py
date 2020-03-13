"""
Module for guiding Bias subtraction including generating a Bias image as desired

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import numpy as np
import os
from IPython import embed

from pypeit import msgs
from pypeit import masterframe
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage


class BiasImage(calibrationimage.CalibrationImage):
    # Set the version of this class
    version = '1.0.0'

    # Output to disk
    output_to_disk = ('BIAS_IMAGE',)
    hdu_prefix = 'BIAS_'
    master_type = 'Bias'
    frametype = 'bias'


class BiasFrame(calibrationimage.BuildCalibrationImage):
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
    master_version = '1.0.0'
    image_type = BiasImage

    #@classmethod
    #def from_master_file(cls, master_file, par=None):
    #    """
    #    Instantiate from a master file
#
#        Args:
#            master_file (str):
#            par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):
#
#        Returns:
#            biasframe.BiasFrame:
#                The PypeItImage is loaded into self.pypeitImage
#
#        """
#        # Spectrograph
#        spectrograph, extras = masterframe.items_from_master_file(master_file)
#        head0 = extras[0]
#        # Master info
#        master_dir = head0['MSTRDIR']
#        master_key = head0['MSTRKEY']
#        # Instantiate
#        slf = cls(spectrograph, par=par, master_dir=master_dir, master_key=master_key,
#                  reuse_masters=True)
#        slf.pypeitImage = slf.load(ifile=master_file)
#        # Return
#        return slf

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False):

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.BuildCalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)

        # Processing steps
        self.process_steps = []
        if self.par['process']['overscan'].lower() != 'none':
            self.process_steps.append('subtract_overscan')
        self.process_steps += ['trim']
        self.process_steps += ['orient']

    def build_image(self, overwrite=False, trim=True):
        """
        Grab the bias files (as needed) and then process the input bias
        frames with :func:`pypeit.processimages.ProcessImages.process`.

        Args:
            overwrite: (:obj: `bool`, optional):
                Regenerate the combined image
            trim (:obj:`bool`, optional):
                If True, trim the image

        Returns:
            `numpy.ndarray`_: Combined, processed image.
        """
        # Nothing?
        if self.par['useframe'].lower() == 'none':
            msgs.info("Bias image subtraction not activated.")
            return None
        if self.nfiles == 0:
            msgs.info("No bias frames provided.  No bias image will be generated or used")
            return None
        # Build
        pypeitimage = super(BiasFrame, self).build_image(ignore_saturation=True)
        biasImage = BiasImage.from_pypeitimage(pypeitimage)
        # Return
        return biasImage

    #def save(self, outfile=None, overwrite=True):
    #    """
    #    Save the bias master data.
#
#        Args:
#            outfile (:obj:`str`, optional):
#                Name for the output file.  Defaults to
#                :attr:`file_path`.
#            overwrite (:obj:`bool`, optional):
#                Overwrite any existing file.
#        """
#        # Some checks
#        if self.pypeitImage is None:
#            msgs.warn('No MasterBias to save!')
#            return
#        # Proceed
#        _outfile = self.master_file_path if outfile is None else outfile
#        # Check if it exists
#        if os.path.exists(_outfile) and not overwrite:
#            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
#                      + 'Set overwrite=True to overwrite it.')
#            return
#        # Save
#        hdr = self.build_master_header(steps=self.process_steps, raw_files=self.file_list)
#        self.pypeitImage.to_file(_outfile, primary_hdr=hdr, hdu_prefix='BIAS_')#, iext='BIAS')
#        msgs.info('Master frame written to {0}'.format(_outfile))
#        #super(BiasFrame, self).save(self.pypeitImage, 'BIAS', outfile=outfile, overwrite=overwrite,
#        #                            raw_files=self.file_list, steps=self.process_steps)

    def load(self, ifile, reuse_masters=False):
        """
        Load the bias frame according to how par['useframe'] is set.
        
        Args:
            ifile (:obj:`str`):
                Name of the master frame file.

        Returns:
            Returns either the `numpy.ndarray`_ with the bias image
            or None if no bias is to be subtracted.
        """
        # How are we treating biases?
        # 1) No bias subtraction
        if self.par['useframe'].lower() == 'none':
            msgs.info("Will not perform bias/dark subtraction")
            return None

        # 2) Use overscan
        if self.par['useframe'] == 'overscan':
            msgs.error("useframe=overscan was Deprecated. Remove it from your pypeit file")

        # 3) User wants bias subtractions
        if self.par['useframe'] in ['bias', 'dark']:
            # Check on whether to reuse and whether the file exists
            if os.path.isfile(ifile) and reuse_masters:
                self.pypeitImage = BiasImage.from_file(ifile, hdu_prefix='BIAS_')
                return self.pypeitImage
            else:
                return None

