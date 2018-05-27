# Module for generating the Flat Field
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from pypit import msgs
from pypit import processimages
from pypit import masterframe
from pypit.core import arprocimg
from pypit.core import arflat

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'normpixelflat'


class FlatField(processimages.ProcessImages, masterframe.MasterFrame):
    """
    This class will generat the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------

    Attributes
    ----------
    frametype : str
      Set to 'arc'

    Inherited Attributes
    --------------------
    stack : ndarray
      Final output image
    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None, msbias=None,
                 det=None):

        # Parameters unique to this Object
        self.msbias = msbias
        self.det = det

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph, settings=settings, det=det)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs
        self.mspixelflat = None

        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings.copy()
        else:
            self.settings = settings.copy()
            # The following is somewhat kludgy and the current way we do settings may
            #   not touch all the options (not sure .update() would help)
            if 'combine' not in settings.keys():
                if self.frametype in settings.keys():
                    self.settings['combine'] = settings[self.frametype]['combine']

        # Child-specific Internals
        #    See ProcessImages for the rest

        # MasterFrames
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    def apply_gain(self, datasec_img, settings_det):
        # Apply gain (instead of ampsec scale)
        self.mspixelflat *= arprocimg.gain_frame(datasec_img, settings_det['numamplifiers'], settings_det['gain'])

    def build_pixflat(self, trim=True):
        # Generate the image
        self.mspixelflat = self.process(bias_subtract=self.msbias, trim=trim)

    def run(self, datasec_img, lordloc, rordloc, pixwid, slitpix, tilts, settings_det, armed=False):
        #ntcky=settings.argflag['reduce']['flatfield']['params'][0])

        # Build the pixel flat (as needed)
        if self.mspixelflat is None:
            self.mspixelflat = self.build_pixflat()

        # Apply gain
        self.apply_gain(datasec_img, settings_det)

        # Normalize the flat field
        msgs.info("Normalizing the pixel flat")
        slit_profiles, mstracenrm, msblaze, flat_ext1d, extrap_slit = \
            arflat.slit_profile(self.mspixelflat, datasec_img, lordloc, rordloc, pixwid, slitpix,
                                self.det, tilts, ntcky=self.settings['flatfield']['params'][0])

        # If some slit profiles/blaze functions need to be extrapolated, do that now
        if armed:  #settings.spect['mosaic']['reduction'] == 'AMRED':
            if np.sum(extrap_slit) != 0.0:
                slit_profiles, mstracenrm, msblaze = arflat.slit_profile_pca(
                    self.mspixelflat, tilts, msblaze, extrap_slit, slit_profiles,
                    lordloc, rordloc, pixwid, slitpix)

        mspixelflatnrm = mstracenrm.copy()
        winpp = np.where(slit_profiles != 0.0)
        mspixelflatnrm[winpp] /= slit_profiles[winpp]

        # QA
        '''
        if np.array_equal(self._idx_flat, self._idx_trace):
            # The flat field frame is also being used to trace the slit edges and determine the slit
            # profile. Avoid recalculating the slit profile and blaze function and save them here.
            self.SetFrame(self._msblaze, msblaze, det)
            self.SetFrame(self._slitprof, slit_profiles, det)
            armasters.save_masters(self, det, mftype='slitprof')
            if settings.argflag["reduce"]["slitprofile"]["perform"]:
                msgs.info("Preparing QA of each slit profile")
                #                            arqa.slit_profile(self, mstracenrm, slit_profiles, self._lordloc[det - 1], self._rordloc[det - 1],
                #                                              self._slitpix[det - 1], desc="Slit profile")
                arproc.slit_profile_qa(self, mstracenrm, slit_profiles,
                                       self._lordloc[det - 1], self._rordloc[det - 1],
                                       self._slitpix[det - 1], desc="Slit profile")
            msgs.info("Saving blaze function QA")
            #                        arqa.plot_orderfits(self, msblaze, flat_ext1d, desc="Blaze function")
            artracewave.plot_orderfits(self.setup, msblaze, flat_ext1d, desc="Blaze function")
        #
        '''
