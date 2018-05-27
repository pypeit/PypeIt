# Module for generating the Flat Field
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload

from pypit import msgs
from pypit import processimages
from pypit import masterframe
from pypit.core import arprocimg
from pypit.core import arflat
from pypit import ginga

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'normpixelflat'

default_settings = dict(flatfield={'method': 'bspline',
                                   "params": [20],
                                   },
                        slitprofile={'perform': True,
                                     },
                        )

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
                 slits_dict=None, tilts=None, det=None, setup=None):

        # Parameters unique to this Object
        self.msbias = msbias
        self.det = det
        self.setup = setup
        self.slits_dict = slits_dict
        self.tilts = tilts

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph, settings=settings, det=det)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs
        self.mspixelflat = None
        self.mspixelflatnrm = None

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
        if 'flatfield' not in settings.keys():
            self.settings.update(default_settings)

        # Child-specific Internals
        #    See ProcessImages and MasterFrame for others
        self.extrap_slit = None
        self.msblaze = None
        self.blazeext = None
        self.slit_profiles = None

        # MasterFrames
        masterframe.MasterFrame.__init__(self, self.frametype, self.setup, self.settings)

    @property
    def nslits(self):
        if self.slits_dict is not None:
            return self.slits_dict['lcen'].shape[1]
        else:
            return 0

    def apply_gain(self, datasec_img):
        # Apply gain (instead of ampsec scale)
        self.mspixelflat *= arprocimg.gain_frame(datasec_img,
                                                 self.settings['detector']['numamplifiers'],
                                                 self.settings['detector']['gain'])
        # Step
        self.steps.append(inspect.stack()[0][3])

    def build_pixflat(self, trim=True):
        # Generate the image
        self.mspixelflat = self.process(bias_subtract=self.msbias, trim=trim)
        # Step
        self.steps.append(inspect.stack()[0][3])

    def _prep_tck(self):
        reload(arflat)
        # Step
        self.steps.append(inspect.stack()[0][3])
        self.ntckx, self.ntcky = arflat.prep_ntck(self.slits_dict['pixwid'], self.settings)

    def slit_profile(self, slit, ntckx=3, ntcky=20):
        reload(arflat)
        # Wrap me
        slordloc = self.slits_dict['lcen'][:,slit]
        srordloc = self.slits_dict['rcen'][:,slit]
        modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = arflat.slit_profile(
            slit, self.mspixelflat, self.tilts, slordloc, srordloc,
            self.slits_dict['slitpix'], self.slits_dict['pixwid'],
            ntckx=ntckx, ntcky=ntcky)
        # Step
        step = inspect.stack()[0][3]
        if step not in self.steps:  # Only do it once
            self.steps.append(step)
        return modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit

    def run(self, datasec_img, armed=False):
        #ntcky=settings.argflag['reduce']['flatfield']['params'][0]) : Default = 20

        # Build the pixel flat (as needed)
        if self.mspixelflat is None:
            self.mspixelflat = self.build_pixflat()

        # Apply gain
        self.apply_gain(datasec_img)

        # Prep tck
        self._prep_tck()

        # Setup
        self.extrap_slit = np.zeros(self.nslits, dtype=np.int)
        self.mspixelflatnrm = self.mspixelflat.copy()
        self.msblaze = np.ones_like(self.slits_dict['lcen'])
        self.blazeext = np.ones_like(self.slits_dict['lcen'])
        self.slit_profiles = np.ones_like(self.mspixelflat)

        # Loop on slits
        for slit in range(self.nslits):
            # Normalize a single slit
            modvals, nrmvals, msblaze_slit, blazeext_slit, iextrap_slit = self.slit_profile(
                slit, ntckx=self.ntckx, ntcky=self.ntcky)

            word = np.where(self.slits_dict['slitpix'] == slit+1)
            self.extrap_slit[slit] = iextrap_slit
            if modvals is None:
                continue
            # Fill
            if self.settings["slitprofile"]["perform"]:
                # Leave slit_profiles as ones if the slitprofile is not being determined, otherwise, set the model.
                self.slit_profiles[word] = modvals/nrmvals
            self.mspixelflatnrm[word] /= nrmvals
            # Fill
            self.msblaze[:,slit] = msblaze_slit
            self.blazeext[:,slit] = blazeext_slit


        # If some slit profiles/blaze functions need to be extrapolated, do that now
        if armed:
            if np.sum(self.extrap_slit) != 0.0:
                slit_profiles, mstracenrm, msblaze = arflat.slit_profile_pca(
                    self.mspixelflat, self.tilts, self.msblaze, self.extrap_slit, self.slit_profiles,
                    self.slits_dict['lcen'], self.slits_dict['rcen'], self.slits_dict['pixwid'],
                    self.slits_dict['slitpix'])

        # Apply slit profile
        winpp = np.where(self.slit_profiles != 0.0)
        self.mspixelflatnrm[winpp] /= self.slit_profiles[winpp]

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

        # Return
        return self.mspixelflatnrm

    def show(self, attr, slit=None, display='ginga', cname=None):
        if attr == 'mspixelflat':
            if self.mspixelflat is not None:
                ginga.show_image(self.mspixelflat)
        elif attr == 'norm':
            if self.mspixelflatnrm is not None:
                ginga.show_image(self.mspixelflatnrm)
