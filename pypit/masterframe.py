""" Abstract class for Master Image frames
      This could go inside armasters.py
"""
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from astropy.io import fits

from pypit import msgs
from pypit import ardebug as debugger
from pypit import armasters

from abc import ABCMeta

default_settings = dict(masters={
    'directory': 'MF/',  # local to the run directory
    'force': False,
    'loaded': [],
    'setup': None,
})


class MasterFrame(object):
    """
    This class is designed to gather a set of like methods
    for Master calibration images, e.g. bias, arc

    Parameters
    ----------
    frametype : str

    Attributes
    ----------
    """
    __metaclass__ = ABCMeta

    def __init__(self, frametype, setup, settings):
        # Start us up

        # Parameters
        self.frametype = frametype
        self.setup = setup
        self.settings = settings
        # Kludge settings a bit for now
        if 'masters' not in self.settings.keys():
            self.settings['masters']['directory'] = self.settings['run']['directory']['master']+'_'+self.settings['run']['spectrograph']
            self.settings['masters']['loaded'] = settings['reduce']['masters']['loaded']

    @property
    def ms_name(self):
        ms_name = armasters.core_master_name(self.frametype, self.setup,
                                             self.settings['masters']['directory'])
        return ms_name

    @property
    def mdir(self):
        return self.settings['masters']['directory']

    def load_master(self, exten=0):
        """

        Parameters
        ----------
        exten

        Returns
        -------

        """
        return armasters.core_load_master_frame(self.ms_name, frametype=self.frametype, exten=exten)

    def save_master(self, image, outfile=None, raw_files=None, steps=None):
        """
        Save the stacked image as a MasterFrame FITS file
          Primarily a wrapper to armasters.save_master

        Parameters
        ----------
        image : ndarray
        outfile : str (optional)
        raw_files : list (optional)
        steps : list (optional)
        """
        if outfile is None:
            outfile = self.ms_name
        # Steps
        if steps is not None:
            jsteps = ','
            jsteps.join(steps)
            keywds=dict(steps=jsteps)
        else:
            keywds = None
        # Finish
        armasters.core_save_master(image, filename=outfile,
                                   raw_files=raw_files, keywds=keywds,
                                   frametype=self.frametype)

