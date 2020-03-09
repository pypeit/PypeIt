"""
Module for generating the Arc image.

.. include common links, assuming primary doc root is up one directory
.. include:: ../links.rst
"""
import os

from pypeit import msgs
from pypeit.par import pypeitpar
from pypeit.images import calibrationimage
from pypeit.images import pypeitimage
from pypeit.core import procimg

from IPython import embed

class ArcImage(calibrationimage.CalibrationImage):

    # Fix the version of this class by PypeItImage
    version = pypeitimage.PypeItImage.version

    # Output to disk
    output_to_disk = ('ARC_IMAGE', 'ARC_FULLMASK')
    hdu_prefix = 'ARC_'
    master_type = 'Arc'
    frametype = 'arc'


class BuildArcImage(calibrationimage.BuildCalibrationImage):#, masterframe.MasterFrame):
    """
    Generate an ArcImage by processing and combining one or more arc frames.


    Args:
        spectrograph (:class:`pypeit.spectrographs.spectrograph.Spectrograph`):
            The `Spectrograph` instance that sets the
            instrument used to take the observations.  Used to set
            :attr:`spectrograph`.
        files (:obj:`list`, optional):
            The list of files to process.
            Can be an empty list or None
        det (:obj:`int`, optional):
            The 1-indexed detector number to process.
        par (:class:`pypeit.par.pypeitpar.FrameGroupPar`):
            The parameters used to type and process the arc frames.
        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (str, optional):
            Path to master frames
        reuse_masters (bool, optional):
            Load from disk if possible
        msbias (ndarray or str, optional):
            Guides bias subtraction

    Attributes:
        msbias (ndarray):
            Bias image or bias-subtraction method; see
            :func:`pypeit.processimages.ProcessImages.process`.
    """
    # Frametype is a class attribute
    frametype = 'arc'
    image_type = ArcImage

    def __init__(self, spectrograph, files=None, det=1, par=None, msbias=None):
    
        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.BuildCalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        #masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
        #                                 master_key=master_key, reuse_masters=reuse_masters)
        # Process steps
        self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['orient']
        self.process_steps += ['apply_gain']

    '''
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
        _outfile = self.master_file_path if outfile is None else outfile
        # Check if it exists
        if os.path.exists(_outfile) and not overwrite:
            msgs.warn('Master file exists: {0}'.format(_outfile) + msgs.newline()
                      + 'Set overwrite=True to overwrite it.')
            return
        #
        hdr = self.build_master_header(steps=self.process_steps, raw_files=self.file_list)
        self.pypeitImage.to_file(_outfile, primary_hdr=hdr, hdu_prefix='ARC_')#, iext='ARC')
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
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
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        # Load it up
        self.pypeitImage = pypeitimage.PypeItImage.from_file(master_file, hdu_prefix='ARC_')
        return self.pypeitImage
    '''

