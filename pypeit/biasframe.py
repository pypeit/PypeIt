""" Module for guiding Bias subtraction including generating a Bias image as desired
"""
from __future__ import absolute_import, division, print_function

import inspect
import os


from pypeit import msgs
from pypeit.core import masters
from pypeit.core import fsort
from pypeit import processimages
from pypeit import masterframe
from pypeit.par import pypeitpar

from pypeit import debugger


class BiasFrame(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class is primarily designed to generate a Bias frame for bias subtraction
      It also contains I/O methods for the Master frames of PYPIT
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
    fitstbl : Table (optional)
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
    def __init__(self, spectrograph, file_list=[], det=1, par=None, setup=None, master_dir=None,
                 mode=None, fitstbl=None, sci_ID=None):

        # Parameters
        self.par = pypeitpar.FrameGroupPar(self.frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, spectrograph, file_list=file_list, det=det,
                                             par=self.par['process'])

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        masterframe.MasterFrame.__init__(self, self.frametype, setup, mode=mode, master_dir=master_dir)

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID

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
        # Get all of the bias frames for this science frame
        if self.nfiles == 0:
            self.file_list = fsort.list_of_files(self.fitstbl, 'bias', self.sci_ID)
        # Combine
        self.stack = self.process(bias_subtract=None, trim=trim, overwrite=overwrite)
        #
        return self.stack

    def master(self):
        """
        Load the master frame from disk, as the settings allow
        or return the command
        or return None

        Note that the user-preference currently holds court, e.g.
          'userframe' = 'overscan' will do an Overscan analysis instead
          of loading an existing MasterFrame bias image

        Returns
        -------
        msframe : ndarray or str or None

        """
        # (KBW) Not sure this is how it should be treated if loaded is
        # being deprecated

        # Generate a bias or dark image (or load a pre-made Master by PYPIT)?
        if self.par['useframe'] is None:
            msgs.info("Will not perform bias/dark subtraction")
            return None

        # Simple command?
        if self.par['useframe'] == 'overscan':
            return self.par['useframe']

        if self.par['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists and user requested one to load it
            msframe, header, raw_files = self.load_master_frame()
            if msframe is None:
                return None
        else:
            # It must be a user-specified file the user wishes to load
            msframe_name = os.path.join(self.directory_path, self.par['useframe'])
            msframe, head, _ = masters._core_load(msframe_name, frametype=self.frametype)

        # Put in
        self.stack = msframe
        return msframe.copy()

