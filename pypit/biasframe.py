""" Module for guiding Bias subtraction including generating a Bias image as desired
"""
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os


from pypit import msgs
from pypit.core import armasters
from pypit.core import arsort
from pypit import processimages
from pypit import masterframe

from pypit import ardebug as debugger

# Does not need to be global, but I prefer it
# ----------------------------------------------------------------------
# (KBW) I don't like doing this, but I realized it's not so bad after
# playing around a bit.  E.g., if you do:
#   from pypit import arcimage
#   from pypit import biasframe
# You get:
#   print(arcimage.frametype)   # prints arc
#   print(biasframe.frametype)  # prints bias
# So even though you're defining it here, it still respects the
# namespace.  The only time you'll get into trouble is if you do:
#   from pypit.arcimage import frametype
#   from pypit.biasframe import frametype
# But I don't think you'd ever do that.
# ----------------------------------------------------------------------
frametype = 'bias'

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
#  These are settings beyond those in the Parent class (ProcessImages)
# additional_default_settings = {frametype: {'useframe': 'none'}}


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

    Attributes
    ----------
    frametype : str
      Set to 'bias'

    Inherited Attributes
    --------------------
    stack : ndarray
    """
    # Keep order same as processimages (or else!)
    def __init__(self, spectrograph, file_list=[], det=1, par=None, setup=None, root_path=None,
                 mode=None, fitstbl=None, sci_ID=None):

        # Parameters unique to this Object
        self.fitstbl = fitstbl
        self.sci_ID = sci_ID

        # Parameters
        self.par = pypitpar.FrameGroupPar(frametype) if par is None else par

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph, det=det,
                                             combine_par=self.par['combine'],
                                             lacosmic_par=self.par['lacosmic'])

        # Attributes  (set after ProcessImages call)
        # (KBW) Even given above, I don't understand why this is
        # preferable to `self.frametype = 'bias'`
        # (KBW) This copying to self also happens in MasterFrame and
        # doesn't need to happen here.
        self.frametype = frametype

        # MasterFrames: Specifically pass the ProcessImages-constructed
        # spectrograph even though it really only needs the string name
        directory_path = None if root_path is None \
                                else root_path+'_'+self.spectrograph.spectrograph
        masterframe.MasterFrame.__init__(self, frametype, setup, directory_path=directory_path,
                                         mode=mode)

    def build_image(self, overwrite=False):
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
            self.file_list = arsort.list_of_files(self.fitstbl, 'bias', self.sci_ID)
        # Combine
        self.stack = self.process(bias_subtract=None, trim=False, overwrite=overwrite)
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
        elif self.par['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists and user requested one to load it
            msframe, header, raw_files = self.load_master_frame()
            if msframe is None:
                return None
        # Simple command?
        elif self.par['useframe'] == 'overscan':
            return self.par['useframe']
        # It must be a user-specified file the user wishes to load
        else:
            msframe_name = os.path.join(self.directory_path, self.par['useframe'])
            msframe, head, _ = armasters._core_load(msframe_name, frametype=self.frametype)
        # Put in
        self.stack = msframe
        return msframe.copy()

