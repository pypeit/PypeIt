""" Abstract class for Master Image frames
      This could go inside masters.py
"""
from __future__ import absolute_import, division, print_function

import os
import warnings

from pypeit import msgs
from pypeit.core import masters

from pypeit import debugger

from abc import ABCMeta


class MasterFrame(object):
    """
    This class is designed to gather a set of like methods
    for Master calibration frames

    Parameters
    ----------
    frametype : str
    setup : str
    settings : dict

    Attributes
    ----------
    """
    __metaclass__ = ABCMeta

    def __init__(self, frametype, setup, spectrograph=None, directory_path=None, mode=None):

        # Set the output directory
        spec_name = None
        if spectrograph is not None:
            # Get the name of the spectrograph used for the directory name.
            if isinstance(spectrograph, str):
                spec_name = spectrograph
            else:
                # Try treating the object as a Spectrograph object
                try:
                    spec_name = spectrograph.spectrograph
                except:
                    pass
        self.directory_path = directory_path
        if self.directory_path is None:
            self.directory_path = os.path.join(os.getcwd(), 'MF')
            if spec_name is not None:
                self.directory_path += '_'+spec_name
        
        # Other parameters
        self.frametype = frametype
        self.setup = setup
        self.mode = mode
        self.msframe = None

    @property
    def ms_name(self):
        return masters.master_name(self.frametype, self.setup, self.directory_path)

    @property
    def mdir(self):
        return self.directory_path

    def _masters_load_chk(self):
        # Logic on whether to load the masters frame
        return self.mode == 'reuse' or self.mode == 'force'

    def load_master_frame(self, force=False):
        """
        Load a MasterFrame

        Parameters
        ----------

        Returns
        -------
        master_frame : ndarray or dict or None
        head0 : Header or None
        file_list : list or None
        """
        if self._masters_load_chk() or force:
            return masters.load_master_frame(self.frametype, self.setup, self.mdir, force=force)
        else:
            return None, None, None

    def master(self):
        """
        Load the master frame from disk, as settings allows

        Returns
        -------
        msframe : ndarray or None
          arc image

        """
        # (KBW) Not sure this is how it should be treated if loaded is
        # being deprecated
        msframe, header, raw_files = self.load_master_frame()
        if msframe is None:
            return None
        self.msframe = msframe
        return msframe.copy()

    def save_master(self, data, outfile=None, raw_files=None, steps=None):
        """
        Save the input data as a MasterFrame file
          Primarily a wrapper to masters.core_save_master

        Intended for simple images only; more complex objects need their own method

        Parameters
        ----------
        data : ndarray or dict
        outfile : str (optional)
        raw_files : list (optional)
        steps : list (optional)
        """
        _outfile = self.ms_name if outfile is None else outfile
        keywds = None if steps is None else dict(steps=','.join(steps))
        masters.save_master(data, filename=_outfile, raw_files=raw_files, keywds=keywds,
                              frametype=self.frametype)

