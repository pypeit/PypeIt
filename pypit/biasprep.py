# Module for guiding Slit/Order tracing
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np

from pypit import msgs
from pypit import ardebug as debugger
from pypit import arcomb
from pypit import arload
from pypit import armasters
from pypit import processimages
from pypit import ginga


# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Place these here or elsewhere?
#  Wherever they be, they need to be defined, described, etc.
#    This assumes the descendents of settings['bias']
default_settings = dict(bias={'useframe': 'none'})


class BiasPrep(processimages.ProcessImages):
    """
    Prepare for Bias Subtraction

    Generate or load a bias image, if desired
    Or return a command for bias removal, e.g.
      'overscan'

    Parameters
    ----------
    settings : dict
      Settings for trace slits
    file_list : list
      List of filenames
    det : int, optional
      Detector index, starts at 1
    user_settings : dict, optional
      Allow for user to over-ride individual internal/default settings
      without providing a full settings dict
    ind : list (optional)
      Indices for bias frames (if a Bias image may be generated)
    fitsdict : dict
      FITS info (mainly for filenames)

    Attributes
    ----------
    images : list
    stack : ndarray
    """
    def __init__(self, setup, settings=None, file_list=[], det=1, user_settings=None,
                 ind=[], fitsdict=None):

        # Start us up
        processimages.ProcessImages.__init__(self, file_list)

        # Parameters
        self.det = det
        self.setup = setup
        self.ind = ind
        self.fitsdict = fitsdict

        # Settings
        if settings is None:
            self.settings = default_settings.copy()
        else:
            # The copy allows up to update settings with user settings without changing the original
            self.settings = settings.copy()
        if user_settings is not None:
            # This only works to replace entire dicts
            #    Hopefully parsets will work more cleanly..
            self.settings.update(user_settings)


    def combine(self, overwrite=False):
        # Over-write?
        if (inspect.stack()[0][3] in self.steps) & (not overwrite):
            msgs.warn("Images already combined.  Use overwrite=True to do it again.")
            return

        # Allow for one-stop-shopping
        if 'load_images' not in self.steps:
            self.load_images()

        # Create proc_images from raw_images if need be
        self.proc_images = np.zeros((self.raw_images[0].shape[0],
                                         self.raw_images[0].shape[1],
                                         self.nloaded))
        for kk,image in enumerate(self.raw_images):
                self.proc_images[:,:,kk] = image
        # Combine
        self.stack = self._combine()
        # Step
        self.steps.append(inspect.stack()[0][3])
        return self.stack.copy()

    def run(self):

        # Generate/load an image??
        if self.settings['bias']['useframe'] in ['bias', 'dark']:
            # Load the MasterFrame if it exists
            #ms_name = armasters.master_name('bias', self.setup)
            msbias = armasters.core_get_master_frame("bias")
            if msbias is None:
                msgs.info("Preparing a master {0:s} frame".format(self.settings['bias']['useframe']))
                # Get all of the bias frames for this science frame
                if self.nfiles == 0:
                    for i in range(len(self.ind)):
                        self.file_list.append(self.fitsdict['directory'][self.ind[i]]+
                                              self.fitsdict['filename'][self.ind[i]])
                # Combine
                msbias = self.combine()
                #frames = arload.load_frames(fitsdict, self.ind, det,
                #                            frametype=settings.argflag['bias']['useframe'], trim=False)
                #msbias = arcomb.comb_frames(frames, det, 'bias', printtype=settings.argflag['bias']['useframe'])
        elif self.settings['bias']['useframe'] in ['overscan', 'none']:
            if self.settings['bias']['useframe'] == 'none':
                msgs.info("Will not perform bias/dark subtraction")
            return self.settings['bias']['useframe']
        else:  # It must be the name of a file the user wishes to load
            msbias_name = self.settings['run']['directory']['master']+u'/'+self.settings['bias']['useframe']
            msbias, head = arload.load_master(msbias_name, frametype="bias")
            # TODO -- Might remove 'loaded' option
            self.settings['reduce']['masters']['loaded'].append('bias')

        self.stack = msbias
        # Save to Masters
        return msbias.copy()


    def show(self, attr, idx=None, display='ginga'):
        if 'proc_image' in attr:
            img = self.proc_images[:,:,idx]
        elif 'raw_image' in attr:
            img = self.raw_images[idx]
        elif 'stack' in attr:
            img = self.stack
        # Show
        viewer, ch = ginga.show_image(img)

    def __repr__(self):
        txt = '<{:s}: nimg={:d}'.format(self.__class__.__name__,
                                         self.nfiles)
        if len(self.steps) > 0:
            txt+= ' steps: ['
            for step in self.steps:
                txt += '{:s}, '.format(step)
            txt = txt[:-2]+']'  # Trim the trailing comma
        txt += '>'
        return txt


