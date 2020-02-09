"""
Module for generating a Bar image to map constant spatial locations

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
from pypeit.core import procimg


class BarFrame(calibrationimage.CalibrationImage, masterframe.MasterFrame):
    """
    Class to generate/load the bar image

    This class is primarily designed to generate a Bar frame to map constant
    spatial locations.  It also contains I/O methods for the Master frames
    of PypeIt.  The build_master() method will return a simple command
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

                pypeitpar.FrameGroupPar('bar')

        master_key (:obj:`str`, optional):
            The string identifier for the instrument configuration.  See
            :class:`pypeit.masterframe.MasterFrame`.
        master_dir (:obj:`str`, optional):
            Path to master frames
        reuse_masters (:obj:`bool`, optional):
            Load from disk if possible
    """

    # Frame type is a class attribute
    frametype = 'bar'
    master_type = 'Bar'

    @classmethod
    def from_master_file(cls, master_file, par=None):
        """
        Instantiate from a master file

        Args:
            master_file (str):
            par (:class:`pypeit.par.pypeitpar.FrameGroupPar`, optional):

        Returns:
            barframe.BarFrame:
                The PypeItImage is loaded into self.pypeitImage

        """
        # Spectrograph
        spectrograph, extras = masterframe.items_from_master_file(master_file)
        head0 = extras[0]
        # Master info
        master_dir = head0['MSTRDIR']
        master_key = head0['MSTRKEY']
        # Instantiate
        slf = cls(spectrograph, par=par, master_dir=master_dir, master_key=master_key,
                  reuse_masters=True)
        slf.pypeitImage = slf.load(ifile=master_file)
        # Return
        return slf

    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, files=None, det=1, par=None, master_key=None,
                 master_dir=None, reuse_masters=False, msbias=None):

        # Parameters unique to this Object
        self.msbias = msbias

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        calibrationimage.CalibrationImage.__init__(self, spectrograph, det, self.par['process'], files=files)

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.master_type, master_dir=master_dir,
                                         master_key=master_key, reuse_masters=reuse_masters)
        # Process steps
        self.process_steps = procimg.init_process_steps(self.msbias, self.par['process'])
        self.process_steps += ['trim']
        self.process_steps += ['orient']
        self.process_steps += ['apply_gain']

    def save(self, outfile=None, overwrite=True):
        """
        Save the master bar data.

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
        self.pypeitImage.write(_outfile, hdr=hdr, iext='BAR')
        msgs.info('Master frame written to {0}'.format(_outfile))

    def load(self, ifile=None):
        """
        Load the bar frame data from a saved master frame.

        Args:
            ifile (:obj:`str`, optional):
                Name of the master frame file.  Defaults to
                :attr:`file_path`.
            return_header (:obj:`bool`, optional):
                Return the header

        Returns:
            Returns a `numpy.ndarray`_ with the bar master frame image.
            Also returns the primary header, if requested.
        """
        # Check on whether to reuse and whether the file exists
        master_file = self.chk_load_master(ifile)
        if master_file is None:
            return
        # Load it up
        self.pypeitImage = pypeitimage.PypeItImage.from_file(master_file)
        return self.pypeitImage
