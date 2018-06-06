# Module for the ScienceImage class
from __future__ import absolute_import, division, print_function

import inspect
import numpy as np
import os

from importlib import reload

from pypit import msgs
from pypit import processimages
from pypit.core import arprocimg
from pypit.core import arskysub
from pypit import artrace
from pypit import ginga

from pypit import ardebug as debugger

# For out of PYPIT running
if msgs._debug is None:
    debug = debugger.init()
    debug['develop'] = True
    msgs.reset(debug=debug, verbosity=2)

# Does not need to be global, but I prefer it
frametype = 'science'


class ScienceImage(processimages.ProcessImages):
    """
    This class will generat the pixel-level FlatField
      The master() method returns the image

    Parameters
    ----------
    file_list : list
      List of raw files to produce the flat field
    spectrograph : str
    settings : dict-like
    msbias : ndarray or str or None
    tslits_dict : dict
      dict from TraceSlits class (e.g. slitpix)
    tilts : ndarray
      tilts from WaveTilts class
    det : int
    setup : str

    Attributes
    ----------
    frametype : str
      Set to 'pixelflat'
    mspixelflat : ndarray
      Stacked image
    mspixelflatnrm : ndarray
      Normalized flat
    extrap_slit
    msblaze : ndarray
      Blaze function fit to normalize
    blazeext :
    slit_profiles : ndarray
      Slit profile(s)
    self.ntckx : int
      Number of knots in the spatial dimension
    self.ntcky : int
      Number of knots in the spectral dimension

    """
    # Keep order same as processimages (or else!)
    def __init__(self, file_list=[], spectrograph=None, settings=None,
                 tslits_dict=None, tilts=None, det=None, setup=None, datasec_img=None,
                 bpm=None, maskslits=None, pixlocn=None, standard=False):

        # Parameters unique to this Object
        self.det = det
        self.setup = setup
        self.tslits_dict = tslits_dict
        self.tilts = tilts
        self.maskslits = maskslits
        self.pixlocn = pixlocn
        self.standard = standard

        # Start us up
        processimages.ProcessImages.__init__(self, file_list, spectrograph=spectrograph,
                                             settings=settings, det=det, datasec_img=datasec_img,
                                             bpm=bpm)

        # Attributes (set after init)
        self.frametype = frametype

        # Key outputs/internals
        self.sciframe = None
        self.varframe = None
        self.modelvarframe = None
        self.global_sky = None
        self.tracelist = []


        # Settings
        # The copy allows up to update settings with user settings without changing the original
        if settings is None:
            # Defaults
            self.settings = processimages.default_settings()
        else:
            self.settings = settings.copy()

        # Child-specific Internals
        #    See ProcessImages
        self.crmask = None

    def find_objects(self, doqa=False):
        reload(artrace)

        nslit = len(self.maskslits)
        gdslits = np.where(~self.maskslits)[0]
        self.tracelist = []

        # Loop on good slits
        for slit in range(nslit):
            if slit not in gdslits:
                self.tracelist.append({})
                continue
            tlist = artrace.trace_objects_in_slit(
                self.det, slit, self.sciframe-self.bgframe,
                self.modelvarframe, self.crmask)
            # Append
            self.tracelist += tlist

        # QA?
        if doqa: # and (not msgs._debug['no_qa']):
            obj_trace_qa(slf, sciframe, det, tracelist, root="object_trace", normalize=False)



    def global_skysub(self, settings_skysub):
        reload(arskysub)

        self.global_sky = np.zeros_like(self.sciframe)
        gdslits = np.where(~self.maskslits)[0]

        for slit in gdslits:
            msgs.info("Working on slit: {:d}".format(slit))
            # TODO -- Replace this try/except when a more stable b-spline is used..
            try:
                slit_bgframe = arskysub.bg_subtraction_slit(self.tslits_dict, self.pixlocn,
                                                            slit, self.tilts, self.sciframe,
                                                            self.varframe, self.bpm,
                                                            self.crmask, settings_skysub)
            except ValueError:  # Should have been bspline..
                msgs.warn("B-spline sky subtraction failed.  Slit {:d} will no longer be processed..".format(slit))
                #msgs.warn("Continue if you wish..")
                debugger.set_trace()
                self.maskslits[slit] = True
            else:
                self.global_sky += slit_bgframe
        # Step
        self.steps.append(inspect.stack()[0][3])
        # Return
        return self.global_sky

    def build_modelvar(self):
        self.modelvarframe = arprocimg.variance_frame(
                self.datasec_img, self.det, self.sciframe, skyframe=self.bgframe)
        return self.modelvarframe

    def _process(self, bias_subtract, pixel_flat, apply_gain=True):
        self.sciframe = self.process(bias_subtract=bias_subtract, apply_gain=apply_gain, pixel_flat=pixel_flat)
        return self.sciframe


    def show(self, attr, image=None, display='ginga'):
        """
        Show one of the internal images

        Parameters
        ----------
        attr : str
          bgframe -- Show the combined flat image, unnormalized
          norm -- Show the combined normalized flat image
        display : str, optional

        Returns
        -------

        """
        if attr == 'global':
            if self.global_sky is not None:
                ginga.show_image(self.global_sky, chname='Global')
        elif attr == 'sci':
            if self.sciframe is not None:
                ginga.show_image(self.sciframe)
        elif attr == 'rawvar':
            if self.varframe is not None:
                ginga.show_image(self.varframe)
        elif attr == 'modelvar':
            if self.modelvarframe is not None:
                ginga.show_image(self.modelvarframe)
        elif attr == 'crmasked':
            ginga.show_image(self.sciframe*(1-self.crmask), chname='CR')
        elif attr == 'image':
            ginga.show_image(image)
        else:
            msgs.warn("Not an option for show")
